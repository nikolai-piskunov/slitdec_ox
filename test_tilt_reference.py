"""The tilt is evaluated at the y-subpixel centre in slitdec_ox (t_offset 0.5).
slitdec.c (t_offset 0) shifts every line by c1/(2*osample)."""
import numpy as np, ctypes as ct
from synth import make_swath, fp_spec, box_int
import pyslitdec_ox as P
NC, NR = 400, 136; TR = slice(30, NC - 30)
fp = fp_spec(NC)

def shift_fit(e, t):
    X = np.arange(NC)[TR]; best = None
    for s in np.linspace(-0.012, 0.012, 481):
        m = 3000 * box_int(fp, X - 0.5 - s, X + 0.5 - s); a = (m @ e[TR]) / (m @ m); c = ((e[TR] - a * m) ** 2).sum()
        if best is None or c < best[0]: best = (c, s)
    return best[1]

if __name__ == "__main__":
    L = P._lib("libslitdec_ox.so")
    mfp, ycen, sc, tfp = make_swath(fp, NC, NR, c1=8 / NR)
    for os_ in (4, 8, 16):
        r0 = P.slitdec_orig(mfp, ycen=ycen, slitcurve=sc, osample=os_, kappa=0)
        r = P.slitdec_ox(mfp, ycen=ycen, slitcurve=sc, osample=os_, osample_x=16, kappa=0)
        print(f"osample {os_:2}: line shift slitdec.c {shift_fit(r0['spectrum'], tfp)*1e3:+.2f} mpx, "
              f"slitdec_ox {shift_fit(r['spectrum'], tfp)*1e3:+.2f} mpx  (c1/(2 osample) = {8/NR/os_/2*1e3:.2f})")
