"""Choosing osample_x and mu_x for slitdec_ox, and checking a given tilt geometry.

Model (stationary, continuum limit). Row r (weight q_r = w_r sL_r^2) records the
spectrum shifted by s_r and integrated over a pixel. For detector frequency w the
fine spectrum has aliases nu_m = w + 2 pi m. Their information matrix is

    G_mm' = F B(nu_m) B(nu_m') c_(m - m'),   c_n = sum_r q_r exp(2 pi i n s_r) / sum_r q_r,

B the pixel box transfer, and the penalty adds mu F nu_m^(2k). The output (the
pixel integral) has noise N(w) and transfer T(w) given in predict(). The tilt
enters only through the phase coefficients c_n:

  |c_n| << 1 for n = 1..3  -> aliases decouple; N and T do not depend on the tilt
                              at all, and var/mean = (1/2 pi) int R(nu)^2 dnu with
                              R = B^2 / (B^2 + mu nu^2k). var/mean = 1 fixes
                              mu = 1.64e-4 (k = 3), 1.6e-3 (k = 2), 1.67e-5 (k = 4).
  |c_1| -> 1 (shear << 1 px) -> rows see the same phase; the fine bins are only held
                              by the penalty and the output stays ~Poisson (0.94).

So mu_x needs no geometry: one value for every order, column and flux. osample_x
is set by how well piecewise-constant bins represent the narrowest line
(measured with the C code: rms line error ~ 0.025 (h / sigma)^2, h = 1/osample_x),
not by the tilt. The geometry check below tells whether a particular order is in
the decoupled regime and what to expect.
"""
import numpy as np

MU_WHITE = {2: 1.6e-3, 3: 1.64e-4, 4: 1.67e-5}


def choose_osample_x(fwhm_min, eps=1e-4, c=0.025):
    """Smallest power of two with rms line-profile error < eps for lines of
    FWHM >= fwhm_min [px]."""
    sigma = fwhm_min / 2.3548
    ox = np.sqrt(c / eps) / sigma
    return int(2 ** np.ceil(np.log2(max(ox, 1))))


def subpixel_shifts(slitcurve_col, ycen_frac, nrows, osample, t_offset=0.5):
    """Shift of every y-subpixel of one column, and its height weight, in slitdec's
    convention (pixel y covers [y, y+1), centre at nrows//2 + ycen_frac)."""
    step = 1.0 / osample
    yc = nrows // 2 + ycen_frac
    # subpixel boundaries aligned to the centre line, as in slitdec
    grid = yc + (np.arange(-int(np.ceil(yc / step)) - 1, int(np.ceil((nrows - yc) / step)) + 1)) * step
    lower = grid[:-1]; upper = grid[1:]
    lower_c = np.clip(lower, 0, nrows); upper_c = np.clip(upper, 0, nrows)
    h = upper_c - lower_c
    keep = h > 0
    t = (lower + t_offset * step - yc)[keep]
    c = np.zeros(6); c[:len(slitcurve_col)] = slitcurve_col
    s = sum(c[k] * t ** k for k in range(1, 6))
    return t, s, h[keep]


def phase_coefficients(shifts, q, nmax=3):
    q = np.asarray(q, float) / np.sum(q)
    return np.array([np.abs(np.sum(q * np.exp(2j * np.pi * n * shifts))) for n in range(1, nmax + 1)])


def predict(shifts, q, mu=MU_WHITE[3], k=3, M=6, nw=257):
    """Predicted var/mean, lag-1 correlation and transfer of the extracted spectrum."""
    q = np.asarray(q, float) / np.sum(q)
    c = np.array([np.sum(q * np.exp(2j * np.pi * n * shifts)) for n in range(-2 * M, 2 * M + 1)])
    ms = np.arange(-M, M + 1)
    C = c[(ms[:, None] - ms[None, :]) + 2 * M]
    ws = np.linspace(-np.pi, np.pi, nw, endpoint=False) + np.pi / nw
    N = np.zeros(nw); T = np.zeros(nw)
    for i, w in enumerate(ws):
        nu = w + 2 * np.pi * ms
        B = np.sinc(nu / 2 / np.pi)
        G = B[:, None] * C * B[None, :]
        Hi = np.linalg.pinv(G + np.diag(mu * nu ** (2 * k)), rcond=1e-13)
        N[i] = np.real(B @ Hi @ G @ Hi.conj().T @ B)
        T[i] = np.real(B @ Hi @ G[:, M]) / B[M]
    var = N.mean()
    return dict(var=var, lag1=(N * np.cos(ws)).mean() / var, w=ws, T=T, N=N)


def check_column(slitcurve_col, ycen_frac, sL_rows, osample=8, mu=MU_WHITE[3], k=3, weights="model"):
    """Diagnostics for one column: sL_rows is the slit illumination per detector row."""
    nrows = len(sL_rows)
    t, s, hgt = subpixel_shifts(slitcurve_col, ycen_frac, nrows, osample)
    yc = nrows // 2 + ycen_frac
    prof = np.interp(t + yc, np.arange(nrows) + 0.5, sL_rows, left=0, right=0) * hgt
    q = prof if weights == "model" else prof ** 2     # w sL^2: photon-limited vs unweighted
    p = predict(s, q, mu, k)
    return dict(shear=s.max() - s.min(), c=phase_coefficients(s, q), var=p["var"], lag1=p["lag1"],
                T_half=np.interp(np.pi / 2, p["w"], p["T"]))


if __name__ == "__main__":
    from scipy.special import erf
    nrows = 136
    y = np.arange(nrows) + 0.5 - nrows / 2
    sL = 0.5 * (erf((y + 60) / 4) - erf((y - 60) / 4))
    print("osample_x for rms line error 1e-4:", {f: choose_osample_x(f) for f in (3.0, 4.5, 6.5)})
    print(f"{'shear [px]':>10} {'|c1|':>6} {'|c2|':>6} {'|c3|':>6} {'var/mean':>9} {'lag1':>6} {'T(pi/2)':>8}")
    for shear in (0.0, 0.25, 0.5, 1, 2, 4, 8, 16):
        d = check_column([0, shear / nrows], 0.3, sL)
        print(f"{shear:10.2f} {d['c'][0]:6.3f} {d['c'][1]:6.3f} {d['c'][2]:6.3f} {d['var']:9.3f} {d['lag1']:+6.2f} {d['T_half']:8.3f}")
