"""Call slit_func_2d in libslitdec_ox_idl.so the way IDL CALL_EXTERNAL does, on the
swath saved by mkslitf.pro (test.sav), and compare with other REDUCE libraries.

usage: python test_idl_interface.py [test.sav] [reference_lib.so ...]
"""
import ctypes as ct
import os
import sys
import numpy as np
from scipy.io import readsav

here = os.path.dirname(os.path.abspath(__file__))
sav = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "..", "test.sav")
refs = sys.argv[2:]
d = readsav(sav)


def native(a, dt):
    return np.ascontiguousarray(np.asarray(a), dtype=dt).copy()


def call(lib, ox=None, mux=(0.0, 0.0, 1.64e-4), muy=(0.0, 0.0), maxiter=100, t_offset=None):
    L = ct.CDLL(lib)
    # the IDL wrapper reads the tilt reference point from the environment on every call
    if t_offset is None:
        os.environ.pop("SLITDEC_OX_T_OFFSET", None)
    else:
        os.environ["SLITDEC_OX_T_OFFSET"] = repr(t_offset)
    nc, nysf, osamp = int(d["nc"]), int(d["nysf"]), int(d["osamp"])
    ny = (nysf + 1) * osamp + 1
    a = [np.array(nc, np.int32), np.array(nysf, np.int32), np.array(osamp, np.int32),
         native(d["sf"], np.float64), native(d["ycen_fract"], np.float64),
         native(d["ycen_int"], np.int32), np.array(int(d["y_lower_lim"]), np.int32),
         native(d["psf_curve"], np.float64), np.array(int(d["idelta_x"]), np.int32),
         np.array(float(d["lambda_sp"]), np.float64), np.array(float(d["lambda_sf"]), np.float64),
         native(d["msk"], np.uint8), np.zeros((nysf, nc)), np.zeros(nc), np.zeros(nc),
         np.zeros(ny), np.zeros(nc)]
    # scratch arrays, sized as mkslitf.pro allocates them (the original code uses them)
    a += [native(d[k], d[k].dtype.newbyteorder("=")) * 0
          for k in ("xi", "zeta", "m_zeta", "l_aij", "l_bj", "p_aij", "p_bj")]
    a += [np.array(0, np.int32)]
    if ox is not None:
        a += [np.array(maxiter, np.int32), np.array(ox, np.int32)]
        a += [np.array(v, np.float64) for v in tuple(mux) + tuple(muy)]
        a += [np.zeros(nc * ox), np.zeros(11)]
    argv = (ct.c_void_p * len(a))(*[x.ctypes.data for x in a])
    L.slit_func_2d.restype = ct.c_int
    rc = L.slit_func_2d(ct.c_int(len(a)), argv)
    return dict(rc=rc, sp=a[13], sL=a[15], model=a[12], unc=a[16], mask=a[11],
                sp_fine=a[-2] if ox is not None else None, info=a[-1] if ox is not None else None)


def report(name, r, ref=None):
    dx = int(d["idelta_x"])
    sl = slice(dx, int(d["nc"]) - dx)
    resid = (r["model"] - d["sf"])[:, sl][r["mask"][:, sl] > 0]
    s = (f"{name:34s} rc={r['rc']:2d} median sP={np.median(r['sp'][sl]):10.1f} "
         f"rms(model-im)={resid.std():8.2f} masked={np.sum(r['mask'] == 0):4d}")
    if ref is not None:
        s += f"  rms(sP/ref-1)={np.std(r['sp'][sl] / ref['sp'][sl] - 1):.2e}"
    print(s)


lib = os.path.join(here, "libslitdec_ox_idl.so")
print(f"{sav}: nc={int(d['nc'])} nysf={int(d['nysf'])} osamp={int(d['osamp'])} "
      f"y_lower_lim={int(d['y_lower_lim'])} shear={np.median(d['psf_curve'][:, 1]):+.4f}")
r_def = call(lib)
r8 = call(lib, 8)
r16 = call(lib, 16)
r1 = call(lib, 1, (0.0, 0.0, 0.0), t_offset=0.0)
assert r_def["rc"] == 0 and np.array_equal(r_def["sp"], r8["sp"]), "defaults != explicit defaults"
assert np.allclose(r8["sp_fine"].reshape(-1, 8).mean(1)[5:-5], r8["sp"][5:-5]), "sP != mean(sP_fine)"
report("slitdec_ox ox=8 (default)", r8)
report("slitdec_ox ox=16", r16, r8)
i = r8["info"]
print(f"  info: damping x {i[7]:.3g}, resolution x {i[8]:.3f} px; lambda_sf={float(d['lambda_sf'])} "
      f"= mu1_y {-i[9]:.3g}, resolution y {i[10]:.3f} px")
ry = call(lib, 8, muy=(-i[9], 0.0))
report("slitdec_ox ox=8, mu1_y = that", ry, r8)
report("slitdec_ox ox=1 mu=0 t_offset=0", r1, r8)
for ref in refs:
    rr = call(ref)
    report(os.path.basename(ref), rr, r8)
    report("  slitdec_ox ox=1 mu=0 t_offset=0", r1, rr)
