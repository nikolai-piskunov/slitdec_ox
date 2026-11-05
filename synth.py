"""Synthetic swaths in slitdec's geometry convention.

Packed swath: nrows rows, the trace centre of column x at row coordinate
y_lower_lim + frac(ycen[x]) (pixel y spans [y, y+1)); integer parts of ycen are
the per-column packing offsets, as in PyReduce. A point at distance t along the
slit from the centre of source column x0 lands at detector column x0 + delta(t),
delta(t) = c1 t + c2 t^2. Pixel values integrate the slit profile in y (nsub
samples per row) and the spectrum exactly in x."""
import numpy as np
from scipy.special import erf


def slit_profile(t, half=60.0, edge=4.0):
    return 0.5 * (erf((t + half) / edge) - erf((t - half) / edge)) / (2 * half)


def box_int(spec, a, b):
    if spec["kind"] == "flat":
        return (b - a)
    s = spec["fwhm"] / 2.3548
    out = spec.get("cont", 0.0) * (b - a)
    for c0, A in zip(spec["c"], spec["A"]):
        out = out + A * 0.5 * (erf((b - c0) / (np.sqrt(2) * s)) - erf((a - c0) / (np.sqrt(2) * s)))
    return out


def fp_spec(ncols, fwhm=4.5, spacing=13.7, cont=0.05):
    c = np.arange(-30, ncols + 30, spacing) + 0.37
    return dict(kind="lines", fwhm=fwhm, c=c, A=np.ones_like(c), cont=cont)


def make_swath(spec, ncols=400, nrows=136, c1=8 / 136, c2=0.0, ycen=None, flux=3000.0,
               profile=slit_profile, nsub=32):
    """Returns expected counts (nrows, ncols), ycen (ncols,), slitcurve (ncols, 3),
    truth spectrum (ncols,) = flux * pixel-integrated spectrum (vertical-slit equivalent)."""
    yll = nrows // 2
    if ycen is None:
        ycen = np.full(ncols, 0.3)
    off = np.floor(ycen).astype(int)
    X = np.arange(ncols)
    # detector-frame centre of column x (packing offset included)
    Yc = lambda x: yll + np.interp(x, X, ycen)
    img = np.zeros((nrows, ncols))
    ysub = (np.arange(nsub) + 0.5) / nsub
    for y in range(nrows):
        for k in ysub:
            Ypk = y + k                       # packed row coordinate
            Ydet = Ypk + off[None, :] * 1.0   # detector row for each column
            Ydet = Ydet.ravel() if np.ndim(Ydet) > 1 else Ydet
            t = Ydet - Yc(X)
            for _ in range(4):                # source column x0 = X - delta(t), t wrt its centre
                d = c1 * t + c2 * t * t
                t = Ydet - Yc(X - d)
            d = c1 * t + c2 * t * t
            img[y] += flux * profile(t) / nsub * box_int(spec, X - 0.5 - d, X + 0.5 - d)
    slitcurve = np.zeros((ncols, 3)); slitcurve[:, 1] = c1; slitcurve[:, 2] = c2
    truth = flux * box_int(spec, X - 0.5, X + 0.5)
    return img, ycen.astype(float), slitcurve, truth
