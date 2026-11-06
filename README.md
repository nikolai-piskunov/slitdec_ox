# slitdec_ox

Optimal extraction of echelle orders by slit decomposition, with the spectrum
oversampled along the dispersion (x) as well as the slit function along the slit (y).

slitdec_ox is a C library derived from the `slitdec` slit-decomposition code used
in [PyReduce](https://github.com/ivh/PyReduce) and
[CharSlit](https://github.com/ivh/CharSlit), with interfaces for Python (ctypes) and
for IDL REDUCE (`CALL_EXTERNAL`).

## What is different from slitdec

- **Oversampled spectrum.** The spectrum is solved on `osample_x` fine bins per
  detector column. Each y-subpixel of the slit function, shifted by the slit
  curvature, contributes to a detector pixel through the exact overlap of each
  shifted fine bin with that pixel. For a tilted slit this removes the correlated
  noise and line-profile bias of the two-column linear weighting: on a synthetic
  swath with an 8 px shear, noise relative to a vertical slit drops from
  var/mean 1.72 to 1.02 and the line-profile error from 8e-3 to 1.2e-4
  (`osample_x = 8`).
- **Penalties in physical units.** Spectrum and slit function are each regularised by
  a sum of derivative penalties, `Σ_k mu_k · F · ∫ (d^k s)^2`, normalised so that
  the result acts as the filter `1 / (1 + Σ_k mu_k ν^2k)` (ν in rad per detector
  pixel), whatever the oversampling, flux or weights. One set of values therefore
  serves every order and instrument setting:
  - spectrum: `mu3_x = 1.64e-4` makes the extraction noise Poisson-like (white) for
    any tilt;
  - slit function: a first-derivative penalty `mu1_y` (around `1e-2`) is
    recommended. Higher orders alone ring at sharp slit edges and amplify noise.

  `mu1 ≥ 2 sqrt(mu2)` guarantees that a mixed penalty cannot ring. The library
  reports the damping ratio and nominal resolution of the penalties it used.
- **Tilt evaluated at the subpixel centre.** slitdec evaluates the slit curvature at
  the lower edge of each y-subpixel, which shifts every line by `c1 / (2 osample)`
  pixels; slitdec_ox removes this shift.
- **Asymmetric extraction windows.** The position of the trace in the swath
  (`y_lower_lim`) is an argument, as IDL REDUCE requires (for example two-fibre slits
  extracted with fractional order separations).

With `osample_x = 1`, no spectrum penalty and the old tilt convention, slitdec_ox
reproduces slitdec bit for bit; `test_regression.py` checks this.

[NOTES.md](NOTES.md) has the derivations, test results and open issues.

## Building

```sh
sh build.sh
```

This builds three shared libraries with the system C compiler:

| library | contents |
|---|---|
| `libslitdec_ox.so` | slitdec_ox, for Python |
| `libslitdec_ox_idl.so` | `slit_func_ox` entry point for IDL REDUCE |
| `libslitdec_orig.so` | unchanged slitdec, for comparisons |

## Python

```python
import pyslitdec_ox as P

r = P.slitdec_ox(im, pix_unc=unc, ycen=ycen, slitcurve=slitcurve,
                 osample=10, osample_x=8,
                 mu_x=(0, 0, 1.64e-4),   # penalties on 1st, 2nd, 3rd derivative of the spectrum
                 mu_y=(1e-2,),           # slit function; None = legacy lambda_sL
                 weights="model")
spectrum, slitfunction = r["spectrum"], r["slitfunction"]
```

Arguments and returned values otherwise follow `pyreduce.cwrappers.slitdec`.
`r["spectrum_fine"]` holds the oversampled spectrum. `r["info"][7:11]` gives the
damping ratio and nominal resolution [px] of the penalties in x and y.

## IDL REDUCE

`libslitdec_ox_idl.so` exports `slit_func_ox`. All 24 arguments are required and there
are no defaults, so every call states the parameters it uses (list in the header of
`slitdec_ox_idl.c`):

```idl
iret = CALL_EXTERNAL(lib, 'slit_func_ox', ncols, nrows, osample, osample_x, im, $
         ycen_frac, ycen_offset, y_lower_lim, PSF_curve, mu1_x, mu2_x, mu3_x, $
         mu1_y, mu2_y, maxiter, kappa, sL_preset, mask, model, sP, sP_fine, sL, unc, info)
```

Compared with REDUCE's `slit_func_2d`, `lambda_sP` and `lambda_sF` are replaced by the
`mu` penalties, and `delta_x` (returned in `info[4]`), `sP_old` and the scratch arrays
are gone. In `mkslitf.pro` the call is selected with a structure that must hold every
parameter:

```idl
mkslitf, ..., SLITDEC_OX={LIBRARY:'/path/to/libslitdec_ox_idl.so', OSAMPLE_X:8L, $
                          MU1_X:0d, MU2_X:0d, MU3_X:1.64d-4, MU1_Y:1d-2, MU2_Y:0d, $
                          MAXITER:100L, KAPPA:10d}
```

`idl_test_ox.pro` compares the library with an existing `reduce.so` on a swath
saved from `mkslitf.pro`.

## Files

| file | purpose |
|---|---|
| `slitdec_ox.c`, `slitdec_ox.h` | the algorithm |
| `slitdec_ox_idl.c` | IDL `CALL_EXTERNAL` wrapper |
| `pyslitdec_ox.py` | Python ctypes wrapper |
| `slitdec_orig.c`, `slitdec.h` | reference slitdec (CharSlit c878d00) |
| `synth.py` | synthetic swaths in slitdec's geometry |
| `ox_geometry.py` | choosing `osample_x` and `mu_x`; checking a given tilt |
| `penalty_design.py`, `sl_penalty_study.py` | comparison of penalty shapes in x and y |
| `test_*.py`, `idl_test_ox.pro` | regression, noise, tilt and interface tests |

## Status

Work in progress. Tested on synthetic swaths, on one CRIRES+ swath saved from
`mkslitf.pro`, and in full SALT-HRS red-arm reductions with IDL REDUCE (white extraction
noise where REDUCE's is blue; differences in the cores of deep narrow lines still being
investigated). Not yet done: speed (about 3× slower than REDUCE per SALT-HRS frame),
memory use at high oversampling.

## Reference

The slit-decomposition method:
N. Piskunov, A. Wehrhahn, T. Marquart, *Optimal extraction of echelle spectra:
Getting the most out of observations*, A&A 646, A32 (2021),
[doi:10.1051/0004-6361/202038293](https://doi.org/10.1051/0004-6361/202038293).

## License

GNU General Public License v3.0, as for PyReduce, from which `slitdec.c` (the basis of
`slitdec_ox.c`, and copied unchanged as `slitdec_orig.c`) is taken. See [LICENSE](LICENSE).
