"""ctypes wrappers for slitdec_ox (and the original slitdec, for comparison).
Arguments and returned dict follow pyreduce.cwrappers.slitdec."""
import ctypes as ct
import os
import numpy as np

_here = os.path.dirname(os.path.abspath(__file__))
WEIGHT_MODES = {"none": 0, "unc": 1, "model": 2, "flux": 3}
_D = ct.POINTER(ct.c_double)
_U = ct.POINTER(ct.c_ubyte)


def _lib(name):
    return ct.CDLL(os.path.join(_here, name))


def _p(a):
    return a.ctypes.data_as(_D)


def _prep(im, pix_unc, mask, ycen, slitcurve, slitdeltas, osample, preset_slitfunc):
    im = np.ascontiguousarray(im, dtype=np.float64)
    nrows, ncols = im.shape
    pix_unc = np.ascontiguousarray(pix_unc if pix_unc is not None else np.ones_like(im), dtype=np.float64)
    mask = np.ascontiguousarray(mask if mask is not None else np.ones(im.shape), dtype=np.uint8).copy()
    ycen = np.ascontiguousarray(ycen, dtype=np.float64).copy()
    sc = np.zeros((ncols, 6))
    slitcurve = np.atleast_2d(np.asarray(slitcurve, dtype=np.float64))
    sc[:, :slitcurve.shape[1]] = slitcurve
    ny = osample * (nrows + 1) + 1
    sd = np.zeros(ny) if slitdeltas is None else np.asarray(slitdeltas, dtype=np.float64).ravel()
    if sd.size == nrows:
        sd = np.interp(np.arange(ny) * (nrows - 1.0) / (ny - 1.0), np.arange(nrows), sd)
    sd = np.ascontiguousarray(sd)
    use_preset = 0
    if preset_slitfunc is not None:
        use_preset = 1
        ps = np.asarray(preset_slitfunc, dtype=np.float64).ravel()
        sL = ps.copy() if ps.size == ny else np.interp(np.arange(ny) * (nrows - 1.0) / (ny - 1.0), np.arange(nrows), ps)
    else:
        sL = np.full(ny, 1.0 / osample)
    return im, pix_unc, mask, ycen, np.ascontiguousarray(sc), sd, np.ascontiguousarray(sL), use_preset, nrows, ncols


def slitdec_ox(im, pix_unc=None, mask=None, ycen=None, slitcurve=None, slitdeltas=None,
               osample=6, osample_x=8, mu_x=1.64e-4, sp_order=3, mu_y=None, lambda_sL=1.0,
               maxiter=20, kappa=10.0, preset_slitfunc=None, weights="none", noise_floor=None,
               sP0=None, y_lower_lim=None, lib="libslitdec_ox.so"):
    """mu_x: scalar = penalty on the sp_order-th derivative only, or a sequence
    (mu1, mu2, ...) of penalties on the 1st, 2nd, ... derivative of the spectrum.
    mu_y: None = slitdec's lambda_sL, or a sequence (mu1, mu2, ...) for the slit
    function. Physical units: transfer 1 / (1 + sum mu_k nu^2k), nu in rad/px.
    y_lower_lim: row of the trace in the swath (None: nrows // 2).
    info[7:11]: damping ratio and nominal resolution [px] in x and in y."""
    L = _lib(lib)
    im, pix_unc, mask, ycen, sc, sd, sL, use_preset, nrows, ncols = _prep(
        im, pix_unc, mask, ycen, slitcurve, slitdeltas, osample, preset_slitfunc)
    if np.ndim(mu_x) == 0:
        mux = np.zeros(sp_order); mux[sp_order - 1] = mu_x
    else:
        mux = np.array(mu_x, dtype=np.float64)
    muy = np.zeros(1) if mu_y is None else np.array(mu_y, dtype=np.float64)
    sP = np.ones(ncols) if sP0 is None else np.array(sP0, dtype=np.float64)
    sPf = np.zeros(ncols * osample_x)
    model = np.zeros((nrows, ncols)); unc = np.zeros(ncols); info = np.zeros(11); info[4] = -11
    L.slitdec_ox2.restype = ct.c_int
    rc = L.slitdec_ox2(ct.c_int(ncols), ct.c_int(nrows), _p(im), _p(pix_unc), mask.ctypes.data_as(_U), _p(ycen),
                       _p(sc), _p(sd), ct.c_int(osample), ct.c_int(osample_x),
                       _p(mux), ct.c_int(mux.size), _p(muy), ct.c_int(0 if mu_y is None else muy.size),
                       ct.c_double(lambda_sL), ct.c_int(maxiter), ct.c_double(kappa),
                       ct.c_int(use_preset), ct.c_int(WEIGHT_MODES[weights]),
                       ct.c_double(-1.0 if noise_floor is None else noise_floor),
                       _p(sP), _p(sPf), _p(sL), _p(model), _p(unc), _p(info),
                       ct.c_int(-1 if y_lower_lim is None else y_lower_lim))
    return dict(spectrum=sP, spectrum_fine=sPf, slitfunction=sL, model=model, uncertainty=unc,
                info=info, mask=mask, return_code=rc)


def slitdec_orig(im, pix_unc=None, mask=None, ycen=None, slitcurve=None, slitdeltas=None,
                 osample=6, lambda_sP=0.0, lambda_sL=1.0, maxiter=20, kappa=10.0,
                 preset_slitfunc=None, weights="none", noise_floor=None, sP0=None,
                 lib="libslitdec_orig.so"):
    L = _lib(lib)
    im, pix_unc, mask, ycen, sc, sd, sL, use_preset, nrows, ncols = _prep(
        im, pix_unc, mask, ycen, slitcurve, slitdeltas, osample, preset_slitfunc)
    sP = np.ones(ncols) if sP0 is None else np.array(sP0, dtype=np.float64)
    model = np.zeros((nrows, ncols)); unc = np.zeros(ncols); info = np.zeros(5)
    rc = L.slitdec(ct.c_int(ncols), ct.c_int(nrows), _p(im), _p(pix_unc), mask.ctypes.data_as(_U), _p(ycen),
                   _p(sc), _p(sd), ct.c_int(osample), ct.c_double(lambda_sP), ct.c_double(lambda_sL),
                   ct.c_int(maxiter), ct.c_double(kappa), ct.c_int(use_preset),
                   ct.c_int(WEIGHT_MODES[weights]), ct.c_double(-1.0 if noise_floor is None else noise_floor),
                   _p(sP), _p(sL), _p(model), _p(unc), _p(info))
    return dict(spectrum=sP, slitfunction=sL, model=model, uncertainty=unc, info=info, mask=mask, return_code=rc)
