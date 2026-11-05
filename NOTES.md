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

## Not done

- Real frames: the ANDES flat/FP in `~/ANDES/E2E/fp_slitAB` are not on this machine.
- Memory: the zeta list holds ≈ (osample+1)(osample_x+1) entries per pixel
  (≈ 90 MB for a 500 × 136 swath at 8/8; the list length is now counted exactly).
- Speed: no flat fast path for osample_x > 1; 5× (ox 8) to 10× (ox 16) slower than
  slitdec.c on a tilted swath.
- The uncertainty estimate is unchanged (residual based, per column).
