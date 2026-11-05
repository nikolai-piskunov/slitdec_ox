"""Slit-function penalty: ringing, resolution and noise of sL from the C code.

Synthetic swaths (40 rows, 200 columns, flat continuum, vertical slit) with a
flat-topped slit with sharp edges and a two-fibre slit. The y-subpixels are
sampled by the drift of the trace across the swath: phase coverage
|c1| = |mean exp(2 pi i ycen)| (1 = every column sees the same phase, the
sub-pixel structure is then set by the penalty alone).

Penalties are compared at equal nominal resolution r = pi / nu_half [px]
(transfer 1 / (1 + P(nu)) = 1/2 at nu_half):
    k=1   mu1 = (r/pi)^2
    k=2   mu2 = (r/pi)^4
    crit  mu1 = 2 l^2, mu2 = l^4, l^2 = (sqrt 2 - 1)(r/pi)^2   (critically damped)
usage: python sl_penalty_study.py [flux=3e4] [edge=0.25]

Metrics on sL, all relative to the plateau / peak and free of grid alignment:
    over   largest overshoot above the plateau (flat top) / of the peaks
    under  most negative value of sL
    edge   10-90 % rise of the lower slit edge [px] (flat top)
    dip    valley / peak between the fibres (two-fibre; truth from the noise-free
           osample 40 profile in the first line)
    noise  rms of sL(noisy) - sL(noise-free) over the slit, Poisson at `flux`.
"""
import sys
import numpy as np
from scipy.special import erf
from synth import make_swath
import pyslitdec_ox as P

NR, NC, OS = 40, 200, 10


EDGE = 0.25


def flat_top(t, half=13.0, edge=None):
    edge = EDGE if edge is None else edge
    return 0.5 * (erf((t + half) / edge) - erf((t - half) / edge)) / (2 * half)


def two_fibre(t, sep=6.0, fwhm=2.5):
    s = fwhm / 2.3548
    return (np.exp(-0.5 * ((t - sep / 2) / s) ** 2) + np.exp(-0.5 * ((t + sep / 2) / s) ** 2)) / (2 * s * np.sqrt(2 * np.pi))


def penalty(kind, r):
    a = r / np.pi
    if kind == "k=1":
        return (a ** 2, 0.0)
    if kind == "k=2":
        return (0.0, a ** 4)
    if kind == "crit":
        l2 = (np.sqrt(2) - 1) * a ** 2
        return (2 * l2, l2 ** 2)
    raise ValueError(kind)


def metrics(sL, prof, os_):
    y = (np.arange(sL.size) + 0.5) / os_
    if prof is flat_top:
        mid = slice(sL.size // 2 - 5 * os_, sL.size // 2 + 5 * os_)
        pl = np.median(sL[mid])
        s = sL / pl
        lo = y[np.argmax(s > 0.1)]; hi = y[np.argmax(s > 0.9)]
        return dict(over=s.max() - 1, under=-s.min(), edge=hi - lo)
    s = sL / sL.max()
    c = sL.size // 2
    w = slice(c - 3 * os_, c + 3 * os_ + 1)
    return dict(over=0.0, under=-s.min(), dip=s[w].min())


def run(img, ycen, sc, mu_y=None, lam=None, os_=OS):
    kw = dict(ycen=ycen, slitcurve=sc, osample=os_, osample_x=1, mu_x=0.0, kappa=0, maxiter=30)
    r = P.slitdec_ox(img, mu_y=mu_y, lambda_sL=lam if lam is not None else 1.0, **kw)
    return r["slitfunction"], r["info"]


if __name__ == "__main__":
    flux = float(sys.argv[1]) if len(sys.argv) > 1 else 3e4
    EDGE = float(sys.argv[2]) if len(sys.argv) > 2 else EDGE    # slit edge width [px]
    rng = np.random.default_rng(1)
    spec = dict(kind="flat")
    X = np.arange(NC)
    for prof in (flat_top, two_fibre):
        print(f"\n=== {prof.__name__}, osample {OS}, flux {flux:g}" + (f", edge {EDGE} px" if prof is flat_top else ""))
        for drift in (0.0, 0.3, 2.0):
            ycen = 0.31 + drift * X / NC
            c1 = abs(np.mean(np.exp(2j * np.pi * ycen)))
            img, ycen, sc, _ = make_swath(spec, NC, NR, c1=0.0, ycen=ycen, flux=flux, profile=prof)
            noisy = rng.poisson(img).astype(float)
            if prof is two_fibre and drift == 2.0:
                ref, _ = run(img, ycen, sc, mu_y=(1e-6, 0.0), os_=40)
                print(f"   truth (osample 40, noise free): dip {metrics(ref, prof, 40)['dip']:.3f}")
            print(f" trace drift {drift} px, |c1| = {c1:.3f}")
            cases = [("legacy lambda_sL=1", None, 1.0), ("legacy lambda_sL=2", None, 2.0)]
            for r in (0.15, 0.3, 0.6):
                for kind in ("k=1", "crit", "k=2"):
                    cases.append((f"{kind:4s} r={r}", penalty(kind, r), None))
            for name, mu, lam in cases:
                s0, info = run(img, ycen, sc, mu_y=mu, lam=lam)
                s1, _ = run(noisy, ycen, sc, mu_y=mu, lam=lam)
                m = metrics(s0, prof, OS)
                sl = s0 > 0.05 * s0.max()
                noise = np.std((s1 - s0)[sl]) / np.median(s0[sl])
                extra = f"edge {m['edge']:.2f} px" if "edge" in m else f"dip {m['dip']:.3f}"
                print(f"   {name:20s} res_y {info[10]:5.2f} px  over {m['over']:+.3f}  under {m['under']:+.3f}  "
                      f"{extra}  noise {noise:.3f}")
