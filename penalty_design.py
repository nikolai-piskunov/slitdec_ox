"""Choosing the mixed penalties mu1, mu2, mu3 (spectrum) objectively.

Penalty P(nu) = sum_k mu_k nu^2k in physical units (nu in rad per detector pixel),
so the data/penalty filter is R = B^2 / (B^2 + P), B the pixel box transfer.
For a penalty mu1 nu^2 + mu2 nu^4 the impulse response of 1 / (1 + P) is a
convolution of two exponentials, i.e. positive (no ringing), iff
    zeta = mu1 / (2 sqrt(mu2)) >= 1          (critical damping at zeta = 1),
the critically damped penalty being (1 + l^2 nu^2)^2 - 1, mu1 = 2 l^2, mu2 = l^4.

For every penalty shape this script finds the strength that whitens the extracted
spectrum (var/mean = 1 in the decoupled-alias regime, the criterion that gave
mu3 = 1.64e-4; see ox_geometry.py) and reports, for the extracted pixel spectrum:
transfer at pi/2, 3pi/4, pi, the lag-1 noise correlation, and the ringing of its
resolution kernel on the pixel grid (sum of negative / sum of positive weights,
and the largest negative weight relative to the central one).
"""
import numpy as np
from scipy.optimize import brentq


def predict_P(P, M=6, nw=256):
    """ox_geometry.predict for decoupled aliases (|c_n| = 0) and penalty P(nu)."""
    ms = np.arange(-M, M + 1)
    ws = np.linspace(-np.pi, np.pi, nw, endpoint=False) + np.pi / nw
    N = np.zeros(nw); T = np.zeros(nw)
    for i, w in enumerate(ws):
        nu = w + 2 * np.pi * ms
        B = np.sinc(nu / 2 / np.pi)
        h = 1.0 / (B ** 2 + P(nu))          # G diagonal = B^2 when the aliases decouple
        N[i] = np.sum(B ** 2 * h * B ** 2 * h)
        T[i] = B[M] ** 2 * h[M]
    var = N.mean()
    return dict(var=var, lag1=(N * np.cos(ws)).mean() / var, w=ws, T=T)


def kernel_ringing(T, ws, nlag=12):
    """Resolution kernel of the extracted pixel spectrum, from its transfer T(w)."""
    j = np.arange(-nlag, nlag + 1)
    k = np.array([np.mean(T * np.cos(ws * jj)) for jj in j])
    neg = -k[k < 0].sum() / k[k > 0].sum()
    return k, neg, -k.min() / k[nlag]


FAMILIES = {
    "k=1": lambda s: (lambda nu: s * nu ** 2),
    "k=2": lambda s: (lambda nu: s * nu ** 4),
    "k=3": lambda s: (lambda nu: s * nu ** 6),
    "crit (1,2)": lambda s: (lambda nu: 2 * s ** 2 * nu ** 2 + s ** 4 * nu ** 4),  # s = l
    "crit (1,2) zeta=2": lambda s: (lambda nu: 4 * s ** 2 * nu ** 2 + s ** 4 * nu ** 4),
    "k=3 + crit(1,2)": None,  # filled in below
}


def whiten(make, lo, hi):
    return brentq(lambda s: predict_P(make(s))["var"] - 1.0, lo, hi, xtol=1e-12, rtol=1e-6)


if __name__ == "__main__":
    print("check: k=3 whitening strength", f"{whiten(FAMILIES['k=3'], 1e-7, 1e-1):.3g}", "(ox_geometry: 1.64e-4)")
    rows = []
    for name, make in FAMILIES.items():
        if make is None:
            continue
        s = whiten(make, 1e-8, 10.0)
        p = predict_P(make(s))
        k, neg, nmax = kernel_ringing(p["T"], p["w"])
        T = lambda w0: np.interp(w0, p["w"], p["T"])
        rows.append((name, s, T(np.pi / 2), T(3 * np.pi / 4), T(np.pi - 1e-3), p["lag1"], neg, nmax, k))
    print(f"{'penalty':>20} {'strength':>10} {'T(pi/2)':>8} {'T(3pi/4)':>8} {'T(pi)':>6} {'lag1':>6} "
          f"{'neg/pos':>8} {'min/ctr':>8}")
    for r in rows:
        print(f"{r[0]:>20} {r[1]:10.3g} {r[2]:8.3f} {r[3]:8.3f} {r[4]:6.3f} {r[5]:+6.3f} {r[6]:8.4f} {r[7]:8.4f}")
    print("\nkernel (lags 0..4):")
    for r in rows:
        print(f"{r[0]:>20}", " ".join(f"{v:+.4f}" for v in r[8][12:17]))
