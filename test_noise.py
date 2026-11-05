import numpy as np, time, sys
from synth import *
from pyslitdec_ox import slitdec_ox, slitdec_orig
NC, NR = 400, 136
c1 = float(sys.argv[1]) / NR if len(sys.argv) > 1 else 8 / NR
weights = sys.argv[2] if len(sys.argv) > 2 else "none"
nreal = int(sys.argv[3]) if len(sys.argv) > 3 else 60
TR = slice(30, NC - 30)
flat = dict(kind="flat")
mflat, ycen, sc, tflat = make_swath(flat, NC, NR, c1=c1)
m0, _, sc0, _ = make_swath(flat, NC, NR, c1=0.0)
fp = fp_spec(NC); mfp, _, _, tfp = make_swath(fp, NC, NR, c1=c1)
rng = np.random.default_rng(3)
noises = [rng.standard_normal(mflat.shape) for _ in range(nreal)]

def stats(fn, model, truth):
    R = []
    for z in noises:
        d = model + np.sqrt(model) * z
        R.append(fn(d)["spectrum"][TR] - truth[TR])
    R = np.array(R); R -= R.mean(0)
    v = R.var(); l1 = (R[:, 1:] * R[:, :-1]).mean() / v
    kb = []
    for k in (2, 4, 8):
        n = R.shape[1] // k; kb.append(R[:, :n*k].reshape(len(R), n, k).sum(2).var() / k)
    return v, l1, kb

ref = lambda d: slitdec_orig(d, pix_unc=np.sqrt(np.maximum(d, 1)), ycen=ycen, slitcurve=sc0, osample=8, kappa=0, weights=weights)
v0, _, kb0 = stats(ref, m0, tflat)
print(f"shear {c1*NR:.1f} px, weights {weights}: vertical-slit reference var {v0:.1f} (expected ~{tflat.mean():.0f})")
confs = [("linear (current)", lambda d: slitdec_orig(d, pix_unc=np.sqrt(np.maximum(d,1)), ycen=ycen, slitcurve=sc, osample=8, kappa=0, weights=weights)),
         ("ox=4  D3 mu=1.64e-4", lambda d: slitdec_ox(d, pix_unc=np.sqrt(np.maximum(d,1)), ycen=ycen, slitcurve=sc, osample=8, osample_x=4, mu_x=1.64e-4, kappa=0, weights=weights)),
         ("ox=8  D3 mu=1.64e-4", lambda d: slitdec_ox(d, pix_unc=np.sqrt(np.maximum(d,1)), ycen=ycen, slitcurve=sc, osample=8, osample_x=8, mu_x=1.64e-4, kappa=0, weights=weights)),
         ("ox=16 D3 mu=1.64e-4", lambda d: slitdec_ox(d, pix_unc=np.sqrt(np.maximum(d,1)), ycen=ycen, slitcurve=sc, osample=8, osample_x=16, mu_x=1.64e-4, kappa=0, weights=weights))]
for name, fn in confs:
    t = time.time(); v, l1, kb = stats(fn, mflat, tflat); dt = (time.time() - t) / nreal
    e = fn(mfp)["spectrum"]; r = (e - tfp)[TR]
    sig = fn(mflat + np.sqrt(mflat) * noises[0])
    print(f"{name:20} var={v/v0:5.3f} lag1={l1:+.2f} k2/4/8={'/'.join(f'{a/b:.3f}' for a,b in zip(kb,kb0))} "
          f"| FP bias rms/peak={np.sqrt((r**2).mean())/tfp.max():.1e} peak={e[TR].max()/tfp[TR].max():.4f} "
          f"| SIG/sqrt(S)={np.median(sig['uncertainty'][TR]/np.sqrt(sig['spectrum'][TR])):.3f} | {dt:.2f} s/call")
