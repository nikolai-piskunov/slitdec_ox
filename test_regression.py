import numpy as np, ctypes as ct
from synth import *
import pyslitdec_ox as P
L = P._lib("libslitdec_ox.so"); ct.c_double.in_dll(L, "slitdec_ox_t_offset").value = 0.0
fp = fp_spec(400)
for c1, yc in ((8/136, None), (0.0, None), (-12/136, 0.37 + 0.004*np.arange(400))):
    img, ycen, sc, truth = make_swath(fp, c1=c1, ycen=yc)
    for w in ("none", "model"):
        kw = dict(ycen=ycen, slitcurve=sc, osample=8, kappa=10, weights=w, pix_unc=np.sqrt(img+1))
        r0 = P.slitdec_orig(img, **kw); r1 = P.slitdec_ox(img, osample_x=1, mu_x=0, **kw)
        same = all(np.array_equal(r0[k], r1[k]) for k in ("spectrum", "slitfunction", "model", "uncertainty", "mask"))
        print(f"c1={c1:+.3f} sloped={yc is not None} weights={w}: identical to slitdec.c: {same}")
