"""
Traversing ELTIM scanner on a UHMWPE UD film: what the sensor really sees, and how to analyse it.

The ELTIM head moves across the film (CD) while the film moves (MD), so the samples lie on a
zigzag path, not on a grid:

    film speed  v_w = 20 m/min   head speed v_s = 5 m/min   sample rate 120 /s   spot 10 mm
    -> one traverse of 1.2 m takes 14.4 s, the film moves 4.8 m MD per traverse
    -> sample spacing 0.7 mm CD x 2.8 mm MD along the path

The film model (rollers along MD, spreading pitch along CD, noise) is the same as in
simulate_eltim.py, but it is evaluated at the sample positions of the scanner. The 10 mm
measuring spot averages the film: each sinusoidal component is multiplied by the spot transfer
function 2*J1(pi*d*f)/(pi*d*f).

Two analyses are compared:
  naive_2dfft()      interpolate the scan samples onto a grid and take a 2D FFT. This fails:
                     the scans are 4.8 m apart along MD while the rollers repeat every 0.3-0.6 m,
                     so the roller peaks are aliased.
  robust_recovery()  works on the raw sample stream (time, CD position, MD position, g/m2):
     1. time spectrum: rollers at f = v_w / circumference, spreading at f = v_s / pitch
        (CD peaks split into sidebands by the head reversals; they are merged)
     2. MD/CD classification: with a second run at another head speed the CD peaks move and
        the MD peaks stay; with a single run the CD pitch from the CD profile is used instead
     3. CD profile: average of all traverses per CD position -> spreading pitch
     4. period refinement + least-squares fit of mean + MD sinusoids + CD sinusoid at the known
        sample positions; amplitudes corrected for the spot size

robust_recovery() can be applied to real ELTIM exports directly: pass the sample arrays and the
line settings.

Simplifications: instant reversal at the edges (no dwell), formation noise generated as
correlated noise along the scan path.

Run:  python scanner_simulation.py
"""
import os
from dataclasses import dataclass

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks as find_peaks_1d
from scipy.special import j1

from simulate_eltim import Params, fft2_amplitude, find_peaks as find_peaks_2d, plot_spectrum, OUT


# ---------------------------------------------------------------- scanner parameters
@dataclass
class Scanner:
    web_speed: float = 20 / 60       # m/s  (20 m/min film speed, MD)
    head_speed: float = 5 / 60       # m/s  (5 m/min traverse speed, CD)
    fs: float = 120.0                # samples per second
    spot_d: float = 0.010            # m, measuring spot diameter
    duration_s: float = 1800.0       # 30 min of production (600 m of film)
    width_m: float = 1.2             # traversed width

    @property
    def traverse_s(self):
        return self.width_m / self.head_speed

    @property
    def md_per_traverse(self):
        return self.web_speed * self.traverse_s


def spot_transfer(f, spot_d=Scanner.spot_d):
    """Amplitude transfer of a circular averaging spot for spatial frequency f (cycles/m)."""
    a = np.pi * spot_d * np.asarray(f, float)
    return np.where(a < 1e-9, 1.0, 2 * j1(a) / np.where(a < 1e-9, 1.0, a))


# ---------------------------------------------------------------- simulated measurement
def scan_path(sc: Scanner):
    t = np.arange(0, sc.duration_s, 1 / sc.fs)
    period = 2 * sc.width_m / sc.head_speed
    phase = (t % period) / period
    x = np.where(phase < 0.5, phase * 2 * sc.width_m, (1 - phase) * 2 * sc.width_m)   # triangle wave
    y = sc.web_speed * t
    return t, x, y


def film_components(p: Params, rng):
    """The periodic components of the synthetic film (true amplitudes and phases)."""
    comps = []
    for r in p.rollers:
        ph = rng.uniform(0, 2 * np.pi, 2)
        comps.append(dict(kind="MD", period_m=r.period_m, amplitude=r.amplitude, phase=ph[0]))
        if r.harmonic2:
            comps.append(dict(kind="MD", period_m=r.period_m / 2, amplitude=r.amplitude * r.harmonic2,
                              phase=ph[1]))
    comps.append(dict(kind="CD", period_m=p.cd_period_m, amplitude=p.cd_amplitude,
                      phase=rng.uniform(0, 2 * np.pi)))
    return comps


def evaluate(comps, x, y, spot_d=None):
    """Sum of periodic components at (x, y); with spot_d the spot averaging is applied."""
    w = np.zeros(np.broadcast(x, y).shape)
    for c in comps:
        pos = y if c["kind"] == "MD" else x
        damp = float(spot_transfer(1 / c["period_m"], spot_d)) if spot_d else 1.0
        w = w + c["amplitude"] * damp * np.cos(2 * np.pi * pos / c["period_m"] + c["phase"])
    return w


def measure(sc: Scanner, p: Params, seed=1):
    """Simulated ELTIM sample stream: (t, x, y, w, true components)."""
    rng = np.random.default_rng(seed)
    t, x, y = scan_path(sc)
    comps = film_components(p, rng)
    w = p.mean_gsm + evaluate(comps, x, y, sc.spot_d)
    path_speed = np.hypot(sc.web_speed, sc.head_speed)
    form = gaussian_filter1d(rng.normal(0, 1, t.size), sigma=p.formation_md_m / path_speed * sc.fs)
    w += p.formation_sigma * form / form.std() + rng.normal(0, p.white_sigma, t.size)
    return t, x, y, w, comps


# ---------------------------------------------------------------- spectral helpers
def amplitude_spectrum_1d(sig, fs):
    sig = sig - sig.mean()
    win = np.hanning(sig.size)
    A = 2 * np.abs(np.fft.rfft(sig * win)) / win.sum()
    f = np.fft.rfftfreq(sig.size, 1 / fs)
    return f, A


def spectrum_peaks(f, A, fmin=0.2, fmax=20.0, factor=12.0, min_amp=0.08):
    band = (f > fmin) & (f < fmax)
    floor = np.median(A[band])
    idx, _ = find_peaks_1d(np.where(band, A, 0), height=max(factor * floor, min_amp),
                           distance=max(1, int(0.05 / (f[1] - f[0]))))
    return idx


def merge_sidebands(cd, head_speed):
    """The head reverses at every edge, which splits a CD peak into sidebands spaced
    1 / (2 x traverse time). Merge each cluster into one peak (amplitude-weighted centre)."""
    out = []
    while cd:
        f0, _ = max(cd, key=lambda c: c[1])
        cluster = [c for c in cd if abs(c[0] - f0) < 0.1 * f0]
        cd = [c for c in cd if c not in cluster]
        fc = sum(f * a for f, a in cluster) / sum(a for _, a in cluster)
        out.append(dict(kind="CD", freq=fc, period_m=head_speed / fc,
                        amp_raw=max(a for _, a in cluster), sidebands=len(cluster)))
    return out


def cd_profile(x, w, width, bin_m=0.002):
    edges = np.arange(0, width + bin_m, bin_m)
    idx = np.clip(np.digitize(x, edges) - 1, 0, len(edges) - 2)
    prof = np.bincount(idx, w, len(edges) - 1) / np.maximum(np.bincount(idx, None, len(edges) - 1), 1)
    return edges[:-1] + bin_m / 2, prof


def refine_period(pos, w, period, rel=0.004, n=161):
    """Fine search around `period` for the best-fitting sinusoid along `pos` (MD or CD).

    A tiny period error accumulates over hundreds of metres and lowers the fitted amplitude,
    so the period from the spectrum is refined before the final fit.
    """
    z = w - w.mean()
    best = (0.0, period)
    for per in period * (1 + np.linspace(-rel, rel, n)):
        ph = 2 * np.pi * pos / per
        a = np.hypot(z @ np.cos(ph), z @ np.sin(ph))
        if a > best[0]:
            best = (a, per)
    return best[1]


# ---------------------------------------------------------------- the two methods
def naive_2dfft(x, y, w, p: Params, first_m=120.0, dx=0.005, dy=0.05):
    """Interpolate the scan samples onto a regular grid and take a 2D FFT (what NOT to do)."""
    m = y < y[0] + first_m
    gx = np.arange(0, x.max(), dx)
    gy = np.arange(y[m][0], y[m][-1], dy)
    GX, GY = np.meshgrid(gx, gy)
    grid = griddata((x[m], y[m]), w[m], (GX, GY), method="linear")
    grid = np.where(np.isnan(grid), np.nanmean(grid), grid)
    fx, fy, A = fft2_amplitude(grid, dx, dy)
    peaks = find_peaks_2d(fx, fy, A, p)
    found = [dict(kind=pk["kind"], period_m=pk["period_mm"] / 1000, amplitude=pk["amplitude"]) for pk in peaks]
    return dict(grid=grid, gx=gx, gy=gy, fx=fx, fy=fy, A=A, peaks=peaks, found=found)


def robust_recovery(x, y, w, *, fs, web_speed, head_speed, width, spot_d, second=None):
    """Recover periodic MD and CD disturbances from a raw traversing-scanner sample stream.

    x, y, w  CD position (m), MD position (m) and areal weight (g/m2) of every sample, in time order
    second   optional (w2, head_speed2): a second run at another head speed, used to classify
             peaks as MD (stay) or CD (move with the head speed)
    Returns a dict with the spectra, the CD profile and the fitted components
    (kind, period_m, amp_measured, amplitude = spot-corrected, phase).
    """
    # 1. time spectrum
    f, A = amplitude_spectrum_1d(w, fs)
    idx = spectrum_peaks(f, A)

    # 3. CD profile (needed for the CD pitch, and for classification with a single run)
    cx, prof = cd_profile(x, w, width)
    fcd, Acd = amplitude_spectrum_1d(prof - prof.mean(), 1 / (cx[1] - cx[0]))
    icd = np.argmax(np.where(fcd > 5, Acd, 0))
    cd_pitch = 1 / fcd[icd]

    # 2. classify the time-spectrum peaks
    md, cd_raw = [], []
    if second is not None:
        w2, head_speed2 = second
        f2, A2 = amplitude_spectrum_1d(w2, fs)
        p2 = f2[spectrum_peaks(f2, A2)]
        ratio = head_speed2 / head_speed
        for i in idx:
            fr = f[i]
            if p2.size and np.min(np.abs(p2 - fr)) < 0.02 * fr:
                md.append(dict(kind="MD", freq=fr, period_m=web_speed / fr, amp_raw=A[i]))
            elif p2.size and np.min(np.abs(p2 - fr * ratio)) < 0.02 * fr * ratio:
                cd_raw.append((fr, A[i]))
    else:
        f2 = A2 = None
        f_cd = head_speed / cd_pitch                  # where the spreading pattern must appear
        for i in idx:
            fr = f[i]
            if abs(fr - f_cd) < 0.1 * f_cd:
                cd_raw.append((fr, A[i]))
            else:
                md.append(dict(kind="MD", freq=fr, period_m=web_speed / fr, amp_raw=A[i]))
    peaks = md + merge_sidebands(cd_raw, head_speed)

    # 4. refine the periods and fit all components together
    md_per = [refine_period(y, w, pk["period_m"]) for pk in md]
    cd_per = [refine_period(x, w, cd_pitch, rel=0.03, n=301)]
    cols, spec = [np.ones_like(w)], []
    for per in md_per:
        cols += [np.cos(2 * np.pi * y / per), np.sin(2 * np.pi * y / per)]
        spec.append(("MD", per))
    for per in cd_per:
        cols += [np.cos(2 * np.pi * x / per), np.sin(2 * np.pi * x / per)]
        spec.append(("CD", per))
    coef, *_ = np.linalg.lstsq(np.column_stack(cols), w, rcond=None)
    comps = []
    for k, (kind, per) in enumerate(spec):
        a, b = coef[1 + 2 * k], coef[2 + 2 * k]
        meas = float(np.hypot(a, b))
        damp = float(spot_transfer(1 / per, spot_d))
        # a cos(th) + b sin(th) = R cos(th + phase) with phase = atan2(-b, a)
        comps.append(dict(kind=kind, period_m=per, amp_measured=meas, amplitude=meas / damp,
                          phase=float(np.arctan2(-b, a))))
    return dict(f=f, A=A, f2=f2, A2=A2, peaks=peaks, cx=cx, prof=prof, fcd=fcd, Acd=Acd,
                cd_pitch=cd_pitch, mean=float(coef[0]), components=comps)


def match(truth, found, tol=0.05):
    """Amplitude of the found component matching each true component (None if not found)."""
    out = []
    for tr in truth:
        cands = [c for c in found if c["kind"] == tr["kind"] and abs(c["period_m"] - tr["period_m"]) < tol * tr["period_m"]]
        out.append(max(cands, key=lambda c: c["amplitude"]) if cands else None)
    return out


# ---------------------------------------------------------------- main
def main():
    os.makedirs(OUT, exist_ok=True)
    p = Params()
    sc = Scanner(width_m=p.width_m)
    sc_fast = Scanner(width_m=p.width_m, head_speed=7.5 / 60)      # second run at another head speed
    print(f"scanner: film {sc.web_speed * 60:.0f} m/min, head {sc.head_speed * 60:.1f} m/min, "
          f"{sc.fs:.0f} samples/s, spot {sc.spot_d * 1000:.0f} mm")
    print(f"  traverse {sc.traverse_s:.1f} s, MD advance per traverse {sc.md_per_traverse:.2f} m, "
          f"sample step {sc.head_speed / sc.fs * 1000:.2f} mm CD x {sc.web_speed / sc.fs * 1000:.2f} mm MD")
    print(f"  spot transfer at the 20 mm pitch: {float(spot_transfer(1 / p.cd_period_m)):.2f}")

    t, x, y, w, truth = measure(sc, p, seed=1)
    _, x2, y2, w2, _ = measure(sc_fast, p, seed=1)
    kw = dict(fs=sc.fs, web_speed=sc.web_speed, head_speed=sc.head_speed, width=sc.width_m, spot_d=sc.spot_d)

    naive = naive_2dfft(x, y, w, p)
    rob1 = robust_recovery(x, y, w, **kw)                                       # single run
    rob2 = robust_recovery(x, y, w, **kw, second=(w2, sc_fast.head_speed))      # two head speeds

    # report
    labels = [f"{c['kind']} {c['period_m'] * 1000:.0f} mm" for c in truth]
    m_naive, m_r1, m_r2 = match(truth, naive["found"]), match(truth, rob1["components"]), match(truth, rob2["components"])
    print("\nNaive 2D FFT peaks (interpolated grid, first 120 m):")
    for pk in naive["peaks"][:6]:
        print(f"   {pk['kind']:5s} period {pk['period_mm']:7.1f} mm  amplitude {pk['amplitude']:.2f} g/m2")
    print("\nRobust recovery (two head speeds), time-spectrum peaks:")
    for pk in rob2["peaks"]:
        extra = f"  ({pk['sidebands']} sidebands merged)" if pk.get("sidebands", 1) > 1 else ""
        print(f"   {pk['kind']}  {pk['freq']:6.3f} Hz -> period {pk['period_m'] * 1000:7.1f} mm{extra}")
    print(f"   CD profile ({int(sc.duration_s / sc.traverse_s)} traverses): pitch {rob2['cd_pitch'] * 1000:.1f} mm")

    def fmt(mm):
        return "not found" if mm is None else f"{mm['period_m'] * 1000:6.1f} mm {mm['amplitude']:5.2f}"
    print("\nComparison (amplitudes in g/m2; robust amplitudes are spot-corrected):")
    print(f"   {'injected':22s} | {'naive 2D FFT':17s} | {'robust, 1 run':17s} | {'robust, 2 head speeds':17s}")
    for tr, a, b, c in zip(truth, m_naive, m_r1, m_r2):
        print(f"   {tr['kind']} {tr['period_m'] * 1000:6.1f} mm {tr['amplitude']:5.2f}    | {fmt(a):17s} | "
              f"{fmt(b):17s} | {fmt(c):17s}")
    false_naive = [f for f in naive["found"] if not any(
        f["kind"] == tr["kind"] and abs(f["period_m"] - tr["period_m"]) < 0.05 * tr["period_m"] for tr in truth)]
    print(f"   naive 2D FFT: {len(false_naive)} false peaks (aliases), e.g. "
          + ", ".join(f"{f['kind']} {f['period_m'] * 1000:.0f} mm" for f in false_naive[:3]))

    # ---------------------------------------------------------------- figures 04-07 (analysis steps)
    vlim = (np.percentile(w, 1), np.percentile(w, 99))
    fig, ax = plt.subplots(figsize=(7, 7))
    mk = y < 30
    s_ = ax.scatter(x[mk], y[mk], c=w[mk], s=0.5, cmap="viridis", vmin=vlim[0], vmax=vlim[1])
    ax.set_xlabel("CD position (m)")
    ax.set_ylabel("MD position (m)")
    ax.set_title(f"Scanner path on the film (first 30 m)\n{sc.md_per_traverse:.1f} m MD per traverse, "
                 f"rollers repeat every 0.31-0.63 m")
    plt.colorbar(s_, ax=ax, label="measured areal weight (g/m²)")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "04_scan_path.png"), dpi=130)

    fig, axs = plt.subplots(1, 2, figsize=(15, 6))
    g = naive
    im = axs[0].imshow(g["grid"], origin="lower", aspect="auto", cmap="viridis",
                       extent=[g["gx"][0], g["gx"][-1], g["gy"][0], g["gy"][-1]], vmin=vlim[0], vmax=vlim[1])
    axs[0].set_ylim(0, 40)
    axs[0].set_xlabel("CD position (m)")
    axs[0].set_ylabel("MD position (m)")
    axs[0].set_title("Scan samples interpolated onto a grid (zoom 40 m)")
    plt.colorbar(im, ax=axs[0], label="areal weight (g/m²)")
    plot_spectrum(axs[1], g["fx"], g["fy"], g["A"], g["peaks"][:3],
                  "Naive 2D FFT of the interpolated grid: aliased peaks", fmax_cd=60, fmax_md=8)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "05_naive_2dfft_scanner.png"), dpi=130)

    fig, ax = plt.subplots(figsize=(12, 5))
    ax.semilogy(rob2["f"], rob2["A"], lw=0.7, label=f"head {sc.head_speed * 60:.1f} m/min")
    ax.semilogy(rob2["f2"], rob2["A2"], lw=0.7, alpha=0.8, label=f"head {sc_fast.head_speed * 60:.1f} m/min")
    for pk in rob2["peaks"]:
        ax.annotate(f"{pk['kind']} {pk['period_m'] * 1000:.0f} mm", (pk["freq"], pk["amp_raw"]),
                    xytext=(0, 10), textcoords="offset points", ha="center", fontsize=8)
    ax.set_xlim(0, 8)
    ax.set_ylim(1e-3, 2)
    ax.set_xlabel("frequency (Hz)")
    ax.set_ylabel("amplitude (g/m²)")
    ax.set_title("Time spectrum of the raw sample stream: MD (roller) peaks stay, CD (spreading) peak moves "
                 "with the head speed")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "06_time_spectrum.png"), dpi=130)

    fig, axs = plt.subplots(1, 2, figsize=(14, 5))
    axs[0].plot(rob2["cx"] * 1000, rob2["prof"], lw=0.8)
    axs[0].set_xlim(0, 300)
    axs[0].set_xlabel("CD position (mm)")
    axs[0].set_ylabel("average areal weight (g/m²)")
    axs[0].set_title(f"CD profile: average of {int(sc.duration_s / sc.traverse_s)} traverses (zoom 300 mm)")
    axs[1].plot(rob2["fcd"], rob2["Acd"], lw=0.8)
    axs[1].set_xlim(0, 120)
    axs[1].set_xlabel("spatial frequency (cycles/m)")
    axs[1].set_ylabel("amplitude (g/m²)")
    axs[1].set_title(f"Spectrum of the CD profile: pitch {rob2['cd_pitch'] * 1000:.1f} mm")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "07_cd_profile.png"), dpi=130)

    # ---------------------------------------------------------------- figure 08: amplitude comparison
    fig, ax = plt.subplots(figsize=(12, 5.5))
    xs = np.arange(len(truth))
    series = [("injected (true film)", [tr["amplitude"] for tr in truth], "#4d4d4d"),
              ("naive 2D FFT of interpolated grid", [0 if m is None else m["amplitude"] for m in m_naive], "#e63946"),
              ("robust, 1 run", [0 if m is None else m["amplitude"] for m in m_r1], "#457b9d"),
              ("robust, 2 head speeds", [0 if m is None else m["amplitude"] for m in m_r2], "#2a9d8f")]
    bw = 0.2
    for k, (name, vals, col) in enumerate(series):
        bars = ax.bar(xs + (k - 1.5) * bw, vals, bw, label=name, color=col)
        for bar, v in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width() / 2, v + 0.02, "not\nfound" if v == 0 else f"{v:.2f}",
                    ha="center", va="bottom", fontsize=7)
    ax.set_xticks(xs)
    ax.set_xticklabels(labels)
    ax.set_ylabel("amplitude (g/m²)")
    ax.set_ylim(0, max(tr["amplitude"] for tr in truth) * 1.35)
    ax.set_title(f"Recovered disturbances from traversing-scanner data ({sc.duration_s / 60:.0f} min, "
                 f"{sc.web_speed * sc.duration_s:.0f} m)\nnaive 2D FFT also reports {len(false_naive)} false "
                 f"(aliased) peaks; robust amplitudes are corrected for the {sc.spot_d * 1000:.0f} mm spot")
    ax.legend(loc="upper left", fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "08_comparison_amplitudes.png"), dpi=130)

    # ---------------------------------------------------------------- figure 09: reconstruction vs truth
    md_true = [c for c in truth if c["kind"] == "MD"]
    cd_true = [c for c in truth if c["kind"] == "CD"]
    md_fit = [c for c in rob2["components"] if c["kind"] == "MD"]
    cd_fit = [c for c in rob2["components"] if c["kind"] == "CD"]
    yy = np.linspace(100, 105, 2000)
    xx = np.linspace(0.3, 0.5, 2000)
    g = naive
    rows = (g["gy"] >= 100 - 1e-9) & (g["gy"] <= 105 + 1e-9)
    naive_md = g["grid"][rows].mean(axis=1) - g["grid"].mean()

    fig, axs = plt.subplots(1, 2, figsize=(15, 5))
    axs[0].plot(yy, evaluate(md_true, 0, yy), color="#4d4d4d", lw=2.2, label="true film (rollers)")
    axs[0].plot(yy, evaluate(md_fit, 0, yy), "--", color="#2a9d8f", lw=1.6, label="robust recovery")
    axs[0].plot(g["gy"][rows], naive_md, color="#e63946", lw=1, alpha=0.8,
                label="interpolated grid, CD-averaged")
    axs[0].set_xlabel("MD position (m)")
    axs[0].set_ylabel("deviation from mean (g/m²)")
    axs[0].set_title("MD variation (rollers), 5 m of film")
    axs[0].legend(fontsize=8)
    cmask = (rob2["cx"] >= 0.3) & (rob2["cx"] <= 0.5)
    axs[1].plot(xx * 1000, evaluate(cd_true, xx, 0), color="#4d4d4d", lw=2.2, label="true film (spreading)")
    axs[1].plot(xx * 1000, evaluate(cd_fit, xx, 0), "--", color="#2a9d8f", lw=1.6,
                label="robust recovery (spot-corrected)")
    axs[1].plot(rob2["cx"][cmask] * 1000, rob2["prof"][cmask] - rob2["prof"].mean(), color="#457b9d", lw=1,
                alpha=0.8, label="measured CD profile (spot-damped)")
    axs[1].set_xlabel("CD position (mm)")
    axs[1].set_ylabel("deviation from mean (g/m²)")
    axs[1].set_title("CD variation (spreading), 200 mm of width")
    axs[1].legend(fontsize=8, loc="lower right")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "09_reconstruction_vs_truth.png"), dpi=130)
    print(f"\nwrote figures 04-09 to {OUT}")


if __name__ == "__main__":
    main()
