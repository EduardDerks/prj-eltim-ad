"""
Synthetic ELTIM ultrasonic basis-weight map of a UHMWPE UD film, with 2D FFT analysis.

The ELTIM gauge measures the areal weight (g/m2) of the moving film. This script builds a
synthetic 2D map  W(MD, CD)  that mimics such a measurement:

  W = mean + random noise                       (white sensor noise + correlated "formation")
        + periodic MD disturbance               (rollers: period = roller circumference)
        + periodic CD disturbance               (fibre spreading: tow / spreader-bar pitch)

and then computes two 2D FFTs:
  1. FFT of the measured map  -> the periodic disturbances show up as sharp peaks
                                 (MD rollers on the f_MD axis, CD spreading on the f_CD axis)
  2. FFT of the residual map  -> after notch-filtering the detected peaks, only noise remains

Axes: MD = machine direction (along the film, y), CD = cross direction (across the width, x).
Spatial frequencies are in cycles/m; periods in mm.

Run:  python simulate_eltim.py            (writes figures and data to outputs/)
      python simulate_eltim.py --seed 7   (different random realisation)
      python simulate_eltim.py --show     (also open the figures)
"""
import argparse
import os
from dataclasses import dataclass, field

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter, maximum_filter

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "outputs")


# ---------------------------------------------------------------- parameters
@dataclass
class Roller:
    name: str
    diameter_m: float          # MD period = pi * diameter
    amplitude: float           # g/m2 (amplitude of the fundamental)
    harmonic2: float = 0.0     # relative amplitude of the 2nd harmonic

    @property
    def period_m(self):
        return np.pi * self.diameter_m


@dataclass
class Params:
    # film and scan geometry
    width_m: float = 1.2               # CD: film width
    length_m: float = 12.0             # MD: scanned length
    dx_m: float = 0.005                # CD resolution (5 mm)
    dy_m: float = 0.010                # MD resolution (10 mm)
    mean_gsm: float = 50.0             # nominal areal weight of the UD film (g/m2)

    # random noise
    white_sigma: float = 0.8           # sensor / measurement noise (g/m2)
    formation_sigma: float = 0.4       # spatially correlated material variation (g/m2)
    formation_md_m: float = 0.15       # correlation length along MD (fibre direction: long)
    formation_cd_m: float = 0.01       # correlation length across CD (short)

    # periodic MD disturbance (rollers): a roller with a flat spot / eccentricity repeats
    # once per revolution -> period = circumference
    rollers: list = field(default_factory=lambda: [
        Roller("roller A (D 200 mm)", 0.200, 0.70, harmonic2=0.25),
        Roller("roller B (D 120 mm)", 0.120, 0.40),
    ])

    # periodic CD disturbance (spreading): tow / spreader-bar pitch leaves streaks along MD
    cd_period_m: float = 0.020         # 20 mm spreading pitch
    cd_amplitude: float = 1.2          # g/m2
    cd_wander_m: float = 0.0           # optional slow lateral wander of the streaks (m, peak)

    # analysis
    peak_factor: float = 8.0           # peak must exceed this times the median spectral level
    min_amplitude: float = 0.10        # g/m2, ignore smaller peaks
    notch_radius: int = 2              # bins removed around each detected peak

    seed: int = 1


# ---------------------------------------------------------------- simulation
def simulate(p: Params):
    """Return (x_cd, y_md, W, parts) with W the areal-weight map in g/m2, shape (n_md, n_cd)."""
    rng = np.random.default_rng(p.seed)
    x = np.arange(0, p.width_m, p.dx_m)         # CD positions
    y = np.arange(0, p.length_m, p.dy_m)        # MD positions
    X, Y = np.meshgrid(x, y)

    white = rng.normal(0, p.white_sigma, X.shape)
    form = gaussian_filter(rng.normal(0, 1, X.shape),                  # elongated along the fibres (MD)
                           sigma=(p.formation_md_m / p.dy_m, p.formation_cd_m / p.dx_m))
    form *= p.formation_sigma / form.std()

    md = np.zeros_like(X)
    for r in p.rollers:
        ph = rng.uniform(0, 2 * np.pi, 2)
        md += r.amplitude * np.cos(2 * np.pi * Y / r.period_m + ph[0])
        md += r.amplitude * r.harmonic2 * np.cos(4 * np.pi * Y / r.period_m + ph[1])

    shift = p.cd_wander_m * np.sin(2 * np.pi * Y / p.length_m)       # optional lateral wander
    cd = p.cd_amplitude * np.cos(2 * np.pi * (X - shift) / p.cd_period_m + rng.uniform(0, 2 * np.pi))

    W = p.mean_gsm + white + form + md + cd
    return x, y, W, dict(white=white, formation=form, md=md, cd=cd)


# ---------------------------------------------------------------- 2D FFT
def fft2_amplitude(Z, dx, dy, window=True):
    """Centred 2D amplitude spectrum, scaled so a sinusoid of amplitude A gives a peak of ~A.

    Returns (fx, fy, A) with fx (CD) and fy (MD) in cycles/m.
    """
    Z = Z - Z.mean()
    ny, nx = Z.shape
    w = np.outer(np.hanning(ny), np.hanning(nx)) if window else np.ones_like(Z)
    F = np.fft.fftshift(np.fft.fft2(Z * w))
    A = 2 * np.abs(F) / w.sum()                 # a real sinusoid splits over +f and -f
    fx = np.fft.fftshift(np.fft.fftfreq(nx, dx))
    fy = np.fft.fftshift(np.fft.fftfreq(ny, dy))
    return fx, fy, A


def find_peaks(fx, fy, A, p: Params):
    """Local maxima well above the noise floor (one of each +/- pair), strongest first."""
    floor = np.median(A)
    local_max = (A == maximum_filter(A, size=5))
    cand = local_max & (A > p.peak_factor * floor) & (A > p.min_amplitude)
    iy, ix = np.nonzero(cand)
    peaks = []
    for j, i in zip(iy, ix):
        f_cd, f_md = fx[i], fy[j]
        if f_cd == 0 and f_md == 0:
            continue
        if f_md < 0 or (f_md == 0 and f_cd < 0):             # keep one half-plane
            continue
        kind = "MD" if abs(f_cd) < 1e-9 else ("CD" if abs(f_md) < 1e-9 else "mixed")
        f = np.hypot(f_cd, f_md)
        peaks.append(dict(kind=kind, f_cd=f_cd, f_md=f_md, period_mm=1000.0 / f, amplitude=A[j, i],
                          ix=i, iy=j))
    return sorted(peaks, key=lambda d: -d["amplitude"])


def notch_filter(W, peaks, radius):
    """Remove the detected periodic components (and their mirror bins) from the map."""
    Z = W - W.mean()
    ny, nx = Z.shape
    F = np.fft.fftshift(np.fft.fft2(Z))
    cy, cx = ny // 2, nx // 2
    yy, xx = np.ogrid[:ny, :nx]
    mask = np.ones_like(F, dtype=float)
    for pk in peaks:
        for sy, sx in ((pk["iy"], pk["ix"]), (2 * cy - pk["iy"], 2 * cx - pk["ix"])):
            mask[(yy - sy) ** 2 + (xx - sx) ** 2 <= radius ** 2] = 0.0
    periodic = np.real(np.fft.ifft2(np.fft.ifftshift(F * (1 - mask))))
    residual = Z - periodic
    return residual + W.mean(), periodic


# ---------------------------------------------------------------- plotting
def plot_map(ax, x, y, W, title, cmap="viridis", vlim=None):
    im = ax.imshow(W, origin="lower", aspect="auto", cmap=cmap,
                   extent=[x[0], x[-1], y[0], y[-1]], vmin=None if vlim is None else vlim[0],
                   vmax=None if vlim is None else vlim[1])
    ax.set_xlabel("CD position (m)")
    ax.set_ylabel("MD position (m)")
    ax.set_title(title)
    plt.colorbar(im, ax=ax, label="areal weight (g/m²)")


def plot_spectrum(ax, fx, fy, A, peaks, title, fmax_cd=60, fmax_md=8, vlim=None):
    Adb = 20 * np.log10(A + 1e-6)
    vmin, vmax = vlim if vlim else (np.percentile(Adb, 50), Adb.max())
    im = ax.imshow(Adb, origin="lower", aspect="auto", cmap="magma",
                   extent=[fx[0], fx[-1], fy[0], fy[-1]], vmin=vmin, vmax=vmax)
    ax.set_xlim(-fmax_cd, fmax_cd)
    ax.set_ylim(-fmax_md, fmax_md)
    ax.set_xlabel("f CD (cycles/m)")
    ax.set_ylabel("f MD (cycles/m)")
    ax.set_title(title)
    plt.colorbar(im, ax=ax, label="amplitude (dB re 1 g/m²)")
    for k, pk in enumerate(peaks):                    # alternate label sides to avoid overlap
        right = k % 2 == 0
        ax.annotate(f"{pk['kind']} {pk['period_mm']:.0f} mm\n{pk['amplitude']:.2f} g/m²",
                    (pk["f_cd"], pk["f_md"]), xytext=(10 if right else -10, 0), textcoords="offset points",
                    ha="left" if right else "right", va="center", color="cyan", fontsize=8)
        ax.plot(pk["f_cd"], pk["f_md"], "o", mfc="none", mec="cyan", ms=10)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--show", action="store_true")
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)

    p = Params(seed=args.seed)
    x, y, W, parts = simulate(p)
    print(f"map: {W.shape[1]} CD x {W.shape[0]} MD points "
          f"({p.width_m} m x {p.length_m} m, {p.dx_m * 1000:.0f} x {p.dy_m * 1000:.0f} mm)")
    print(f"areal weight: mean {W.mean():.2f} g/m2, std {W.std():.2f} g/m2, "
          f"min {W.min():.1f}, max {W.max():.1f}")

    # FFT 1: measured map
    fx, fy, A = fft2_amplitude(W, p.dx_m, p.dy_m)
    peaks = find_peaks(fx, fy, A, p)

    print("\nInjected disturbances:")
    for r in p.rollers:
        print(f"  MD {r.name:22s} period {r.period_m * 1000:6.1f} mm  amplitude {r.amplitude:.2f} g/m2"
              + (f"  (+2nd harmonic {r.amplitude * r.harmonic2:.2f} at {r.period_m * 500:.1f} mm)"
                 if r.harmonic2 else ""))
    print(f"  CD spreading               period {p.cd_period_m * 1000:6.1f} mm  amplitude {p.cd_amplitude:.2f} g/m2")
    print("\nDetected peaks in the 2D FFT (FFT 1):")
    for pk in peaks:
        print(f"  {pk['kind']:5s} period {pk['period_mm']:6.1f} mm  (f_CD {pk['f_cd']:6.2f}, "
              f"f_MD {pk['f_md']:5.2f} cycles/m)  amplitude {pk['amplitude']:.2f} g/m2")

    # FFT 2: residual after removing the detected periodic components
    R, periodic = notch_filter(W, peaks, p.notch_radius)
    fx2, fy2, A2 = fft2_amplitude(R, p.dx_m, p.dy_m)
    peaks2 = find_peaks(fx2, fy2, A2, p)
    noise_true = (parts["white"] + parts["formation"]).std()
    print(f"\nResidual after notch filtering: std {R.std():.2f} g/m2 "
          f"(true random noise std {noise_true:.2f} g/m2), remaining peaks: {len(peaks2)}")

    # figures
    vlim = (np.percentile(W, 1), np.percentile(W, 99))
    fig, axs = plt.subplots(1, 3, figsize=(17, 6), gridspec_kw=dict(width_ratios=[1.3, 1, 1]))
    plot_map(axs[0], x, y, W, "Synthetic ELTIM map, UHMWPE UD film", vlim=vlim)
    axs[1].plot(y, W.mean(axis=1), lw=0.8)
    axs[1].set_xlabel("MD position (m)")
    axs[1].set_ylabel("CD-averaged areal weight (g/m²)")
    axs[1].set_title("MD profile (rollers)")
    axs[1].set_xlim(0, 4)
    axs[2].plot(x * 1000, W.mean(axis=0), lw=0.8)
    axs[2].set_xlabel("CD position (mm)")
    axs[2].set_ylabel("MD-averaged areal weight (g/m²)")
    axs[2].set_title("CD profile (spreading)")
    axs[2].set_xlim(0, 300)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "01_map_and_profiles.png"), dpi=130)

    fig, axs = plt.subplots(1, 2, figsize=(15, 6))
    Adb = 20 * np.log10(A + 1e-6)
    spec_lim = (np.percentile(Adb, 50), Adb.max())               # same colour scale for both FFTs
    plot_spectrum(axs[0], fx, fy, A, peaks, "2D FFT 1: measured map", vlim=spec_lim)
    plot_spectrum(axs[1], fx2, fy2, A2, peaks2, "2D FFT 2: residual after notch filter", vlim=spec_lim)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "02_fft2d_measured_vs_residual.png"), dpi=130)

    fig, axs = plt.subplots(1, 3, figsize=(18, 6))
    plot_map(axs[0], x, y, W, "Measured", vlim=vlim)
    plot_map(axs[1], x, y, periodic + W.mean(), "Periodic part (detected peaks)", vlim=vlim)
    plot_map(axs[2], x, y, R, "Residual (random noise)", vlim=vlim)
    for ax in axs:
        ax.set_ylim(0, 4)
        ax.set_xlim(0, 0.3)
    fig.suptitle("Zoom 0.3 m CD x 4 m MD")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "03_decomposition_zoom.png"), dpi=130)

    # data export
    np.savez_compressed(os.path.join(OUT, "eltim_synthetic.npz"), x_cd_m=x, y_md_m=y, gsm=W,
                        periodic=periodic, residual=R, fx=fx, fy=fy, fft_amplitude=A)
    np.savetxt(os.path.join(OUT, "eltim_synthetic_gsm.csv"), W, delimiter=",", fmt="%.3f",
               header=f"areal weight g/m2; rows = MD ({p.dy_m} m step), columns = CD ({p.dx_m} m step)")
    print(f"\nwrote figures and data to {OUT}")
    if args.show:
        matplotlib.use("TkAgg", force=True)
        plt.show()


if __name__ == "__main__":
    main()
