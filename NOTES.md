# slitdec_ox: slitdec with oversampling in x and y

`slitdec_ox.c` is CharSlit `slitdec.c` (c878d00, as vendored in PyReduce) with the
spectrum solved on `osample_x` fine bins per column. Each y-subpixel of the slit
function, shifted by the curvature, contributes to a detector pixel by the overlap
of each shifted fine bin with that pixel (box integration). The spectrum penalty is a
k-th difference in physical units, `mu_x * F * ∫ (d^k sP/dx^k)^2 dx`, with F the median
data curvature per pixel, so `mu_x` does not depend on `osample_x`, flux or weights.
`sP` comes back on the detector grid (mean of the fine bins = what a vertical slit
records); `sP_fine` optionally returns the fine solution.

Build: `sh build.sh` (two shared libraries; `libslitdec_orig.so` is the unchanged slitdec.c for
comparison). Python: `pyslitdec_ox.slitdec_ox(im, pix_unc, mask, ycen, slitcurve,
slitdeltas, osample=8, osample_x=8, mu_x=1.64e-4, sp_order=3, ...)`, otherwise as
`pyreduce.cwrappers.slitdec`.

## Choosing osample_x and mu_x (ox_geometry.py)

- **mu_x needs no geometry.** When the tilt spreads the row phases (|c_n| = |Σ q e^{2πi n s}|/Σq
  small for n = 1..3, true for shear ≳ 1 px), the aliases decouple and noise and
  transfer are independent of the tilt. Then var/mean = (1/2π)∫R(ν)² dν with
  R = B²/(B² + μν^{2k}); var/mean = 1 gives **μ = 1.64e-4 for k = 3** (1.6e-3 for k = 2,
  1.67e-5 for k = 4). Transfer at that μ: 0.997 at π/2, 0.956 at 3π/4, 0.72 at Nyquist
  (unavoidable, shared by ±π). For shear ≲ 0.5 px the output stays ~Poisson (0.94).
- **osample_x is set by the narrowest line, not by the tilt**: rms line error ≈
  0.025 (h/σ)², h = 1/osample_x. FWHM 4.5 px: 8 → 1e-4, 16 → 2.5e-5.
- **osample (y)**: no new requirement; the shift change over one subpixel, c1/osample,
  only adds c1²/(12 osample²) px² to the line variance (1e-6 relative for ANDES).

## Results (synthetic 400 × 136 swath, flat-topped slit, joint sL/sP fit, D3, μ = 1.64e-4)

Noise relative to a vertical slit; line bias = rms(sP − truth)/peak on a noise-free FP
frame (FWHM 4.5 px).

| shear 8 px, weights `none` | var/mean | lag-1 | line bias | peak | s/call |
|---|---|---|---|---|---|
| slitdec.c (linear) | 1.72 | −0.26 | 8.0e-3 | 1.0184 | 0.05 |
| osample_x = 4 | 1.07 | −0.04 | 4.8e-4 | 1.0011 | 0.14 |
| osample_x = 8 | 1.02 | −0.01 | 1.2e-4 | 1.0002 | 0.26 |
| osample_x = 16 | 1.01 | 0.00 | 2.7e-5 | 1.0000 | 0.52 |

Shear 2 and 16 px and `model` weights give the same picture (osample_x = 8: var/mean
1.02–1.04, bias 1e-4). The reported SIG (residual based) is ≈ √S in all cases, which
is now correct for osample_x ≥ 8 and 1.3× too small for slitdec.c.

## Regression

`slitdec_ox(osample_x=1, mu_x=0)` with `slitdec_ox_t_offset = 0` is bit-identical to
slitdec.c (spectrum, slit function, model, uncertainty, mask; tilted and sloped traces,
kappa on, `none`/`model` weights): `test_regression.py`.

## Finding in slitdec.c: tilt evaluated half a subpixel off

slitdec.c evaluates the curvature at t = (lower edge of the y-subpixel) − ycen instead
of the subpixel centre. Every line then shifts by c1/(2·osample): 3.7 mpx for c1 = 8/136
at osample 8, scaling exactly with 1/osample (`test_tilt_reference.py`). It is smooth in
x and calibrates out if science and calibration share osample and tilt, but it changes
with osample and with any tilt difference (~5 m/s-equivalent for c1 ≈ 0.12 at osample 8 on
ANDES). slitdec_ox evaluates at the centre (`slitdec_ox_t_offset = 0.5`, a global for
now). The same one-line fix applies to slitdec.c.

## IDL REDUCE interface (slitdec_ox_idl.c)

`libslitdec_ox_idl.so` (built by `build.sh`) exports `slit_func_2d` with the 25-argument
CALL_EXTERNAL list of mkslitf.pro, so it replaces `reduce.so` without IDL changes:
`mkslitf, ..., LIBRARY_NAME='~/Slit_func/slitdec_ox/libslitdec_ox_idl.so'` (or link
`slitdec_ox.c` + `slitdec_ox_idl.c` into reduce.so instead of `slitdec.c` + `slitdec_idl.c`).
Optional arguments 25-33 are maxiter, osample_x, mu1_x, mu2_x, mu3_x, mu1_y, mu2_y (all
scalars), sP_fine (out) and info (out, dblarr(11)); when absent, the env vars
`SLITDEC_OX_MAXITER/OSAMPLE_X/MU1_X/MU2_X/MU3_X/MU1_Y/MU2_Y` apply
(`setenv,'SLITDEC_OX_OSAMPLE_X=16'`), then the defaults 100 / 8 / 0 / 0 / 1.64e-4 / 0 / 0.
mu1_y = mu2_y = 0 means REDUCE's `lambda_sf` as before. `lambda_sP` is ignored (mu_x
replaces it); kappa = 10, weights `none`, as in ../slitdec_idl.c.

The trace row is an argument of the general entry point `slitdec_ox2(..., y_lower_lim)`;
mkslitf.pro passes asymmetric values (19 for 40 rows in test.sav, any split for
two-fibre slits) where slitdec_ox() assumes nrows/2. Asymmetric cuts of test.sav
(0-19 rows below, 0-20 above) match the REDUCE library to <= 6e-6 at osample_x = 1,
clean under AddressSanitizer.

`test_idl_interface.py [test.sav] [ref.so ...]` emulates CALL_EXTERNAL on the swath that
mkslitf.pro saved (../test.sav: 412 x 40, osamp 10, shear -0.15 px/row). Result: with
osample_x = 1, mu_x = 0, t_offset = 0 it matches the current REDUCE library
(../slitdec.c + ../slitdec_idl.c) to 9e-16 and the original slit_func_2d to 3e-4 rms. The
default (ox 8) differs from both by 2.9% rms in sP: ~1e-3 in the continuum, up to a few % on
steep line flanks (sP spans 2.4k-57k there), i.e. the linear-interpolation line bias above.

## Penalties: mu_k in x and y (penalty_design.py, sl_penalty_study.py)

`slitdec_ox2` takes `mu_x[k-1]`, `mu_y[k-1]` (k = 1..15): penalty
sum_k mu_k F ∫ (d^k s)^2 in physical units, i.e. the filter 1/(1 + sum_k mu_k nu^2k),
nu in rad/px, independent of osample, osample_x, flux, weights (F = median diagonal
x os^2; coefficient on k-th differences mu_k · diag_med · os^(2k+1)). No mu_y: the
legacy lambda_sL (mean diagonal, 1st differences of subpixels). `slitdec_ox()` keeps
its (mu_x, sp_order, lambda_sL) API. info[4] = -11 on input returns in info[7..10]
the damping ratio ζ = mu1/(2√mu2) and nominal resolution π/ν½ [px] for x and y.
ζ ≥ 1: the penalty kernel is a convolution of two exponentials, so it cannot ring.
For the legacy lambda_sL, info[9] = −(equivalent mu1_y): REDUCE's lambda_sf = 2 at
osamp 10 is mu1_y = 1.9e-3, resolution 0.14 px, i.e. hardly any smoothing.

**x (penalty_design.py):** at the strength that whitens the noise (var/mean = 1),
all shapes (k = 1, 2, 3, critically damped 1+2) give the extracted pixel spectrum
the same ~2-2.5 % negative kernel lobes: in x the ringing comes from the pixel
band limit, not the penalty. k = 3 has the flattest transfer (0.997 at π/2 vs
0.956 for k = 1) and stays the default: mu3_x = 1.64e-4.

**y (sl_penalty_study.py, C code, 40 x 200 synthetic swaths, osample 10):**
penalties compared at equal nominal resolution r.
- Edges blurred to ~1 px (realistic): no penalty rings (≤ 0.5 %).
- Edges sharper than a pixel (0.25 px): every penalty rings with ~1 px period,
  4 % with good phase coverage, 10 % with none (|c1| = 1); the cause is information
  the pixels never recorded, and only smoothing (larger r) reduces it.
- Pure k = 2 is worse in both respects: still rings 4-7 % at r = 1 px where
  k = 1 drops to 0.4-3 %, and amplifies noise 2.5x at small r. This matches
  REDUCE experience that the 1st-derivative term must be kept.
- Critically damped 1+2 ≈ k = 1 within a few %: no measurable gain from mu2_y.
- Realistic 1 px edges: r = 0.3 px (mu1_y = (0.3/π)^2 = 9e-3) keeps the edge width of
  legacy (r ≈ 0.1-0.14) with 1.6x lower noise; r = 0.6 widens edges by 0.2 px.
So use k = 1 in y; mu1_y ≈ 1e-2 is a reasonable start. Not default yet (IDL keeps
lambda_sf). Choosing r per swath from the data (GCV or discrepancy principle) is the
objective next step.

## Not done

- Real frames: the ANDES flat/FP in `~/ANDES/E2E/fp_slitAB` are not on this machine.
- CRIRES+ through IDL REDUCE: interface done and tested in Python only, not yet run inside IDL.
- Memory: the zeta list holds ≈ (osample+1)(osample_x+1) entries per pixel
  (≈ 90 MB for a 500 × 136 swath at 8/8; the list length is now counted exactly).
- Speed: no flat fast path for osample_x > 1; 5× (ox 8) to 10× (ox 16) slower than
  slitdec.c on a tilted swath.
- The uncertainty estimate is unchanged (residual based, per column).
