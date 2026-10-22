# ELTIM areal-weight simulation: UHMWPE UD film
Synthetic 2D areal-weight data (g/m²) of a UHMWPE unidirectional (UD) film, as measured by an **ELTIM ultrasonic basis-weight detector**. The data is analysed with a **2D FFT** to find periodic weight variations caused by the production line.

![2D FFT](outputs/02_fft2d_measured_vs_residual.png)

## What is simulated
The map `W(MD, CD)` covers 1.2 m film width (CD, cross direction) × 12 m length (MD, machine direction), on a 5 mm (CD) × 10 mm (MD) grid with a nominal weight of 50 g/m². It contains:

- **Random noise:**
  - white sensor/measurement noise (σ 0.8 g/m²)
  - spatially correlated material variation, the "formation" (σ 0.4 g/m²). It is elongated along the fibres (150 mm along MD, 10 mm across CD), as typical for a UD film.
- **Periodic MD disturbance (rollers):** a roller with eccentricity or a flat spot repeats once per revolution, so the period is the roller circumference.
  - Roller A, Ø 200 mm: period 628 mm, amplitude 0.70 g/m², plus a 2nd harmonic of 0.17 g/m² at 314 mm
  - Roller B, Ø 120 mm: period 377 mm, amplitude 0.40 g/m²
- **Periodic CD disturbance (spreading):** the fibre spreading pitch leaves streaks along MD, with a period of 20 mm and an amplitude of 1.2 g/m².

All values are parameters in the `Params` class at the top of `simulate_eltim.py`. Set them to match your line: roller diameters, spreader pitch, scan resolution, film width and weight.

## Analysis: two 2D FFTs
1. **2D FFT of the measured map:**
   - The map is mean-removed and Hann-windowed, and the amplitude is scaled so that a sinusoid of amplitude A gives a peak of A g/m².
   - The periodic disturbances appear as sharp peaks: rollers on the f_MD axis, spreading on the f_CD axis.
   - The peaks are detected automatically and reported with their period and amplitude.
2. **2D FFT of the residual:**
   - The detected peaks are removed with a notch filter, and the result is transformed again.
   - No peaks remain, and the residual spread equals the true random noise. This confirms that the periodic components were found and removed correctly.

Result with the default settings (`--seed 1`):
- **Detected:** CD 20.0 mm (1.20 g/m²), MD 632 mm (0.69 g/m²), MD 375 mm (0.40 g/m²) and MD 316 mm (0.17 g/m²). These match the injected disturbances within the FFT bin resolution (1/12 m⁻¹ along MD).
- **Spread:** the total σ is 1.36 g/m². The residual σ after filtering is 0.89 g/m², equal to the true random noise.

## Outputs (`outputs/`)
- `01_map_and_profiles.png`: the map, the MD profile (rollers) and the CD profile (spreading)
- `02_fft2d_measured_vs_residual.png`: both 2D FFTs on the same colour scale, with the detected peaks labelled
- `03_decomposition_zoom.png`: measured map = periodic part + residual (zoom 0.3 m × 4 m)
- `eltim_synthetic.npz`:
  - positions (`x_cd_m`, `y_md_m`)
  - maps (`gsm`, `periodic`, `residual`)
  - FFT axes (`fx`, `fy`) and FFT amplitude (`fft_amplitude`)
- `eltim_synthetic_gsm.csv`: the areal-weight map (rows = MD, columns = CD)

## Usage
Requirements: Python 3.10+.
```powershell
pip install -r requirements.txt
python simulate_eltim.py            # default realisation (seed 1)
python simulate_eltim.py --seed 7   # different random noise and phases
python simulate_eltim.py --show     # also open the figures
```

## Applying it to real ELTIM data
Replace the output of `simulate()` with a measured map:
- shape `(n_md, n_cd)`, in g/m²
- equal spacing `dx_m` (CD) and `dy_m` (MD)

Then use `fft2_amplitude`, `find_peaks` and `notch_filter` as in `main()`. A detected MD period divided by π gives the diameter of the roller causing it, which helps trace the source.
