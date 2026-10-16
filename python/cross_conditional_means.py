#!/usr/bin/env python3
"""
cross_conditional_means.py — cross-conditional group means from the published
coefficient files: E[outcome | conditioning variable in group], e.g. average
wealth of the bottom 50% / next 40% / top 10% of the INCOME distribution.

Why a separate method. The coefficient files carry an order-12 Legendre
representation of the copula, i.e. the moments m_kl = E[Q_k(u) Q_l(v)] of the
pairwise copula density. Plugging the truncated series sum_kl m_kl Q_k Q_l into
the group integral is equivalent to replacing the group indicator 1{u in G} by
its own order-12 approximation, which rings (Gibbs) and lets every group pick
up a few percent of the other variable's top tail. For heavy-tailed outcomes
this biases the means (SCF 2019 microdata: wealth by middle-40% income -9%;
stocks of the bottom 50% of wealth $9,077 instead of $619).

Default method (--method maxent): fit the positive maximum-entropy density
    f(u, v) = exp( sum_kl lambda_kl Q_k(u) Q_l(v) ) / Z
whose moments match m_kl, and integrate the EXACT indicator over it. Same
information (the published coefficients), no ringing; on microdata the group
means come out within about 1% (wealth by income) to a few percent (very
concentrated items). The fit is a damped Newton method on the convex dual with
a small ridge (eps) so it converges even when the smoothed coefficients are not
exactly the moments of a positive density. --method series gives the truncated
series for comparison (see docs/methods.md).

    python python/cross_conditional_means.py                                 # PSID, 2019-Q4
    python python/cross_conditional_means.py --dataset SCF --date all
    python python/cross_conditional_means.py --groups 40,40,20 --method series

Output: python/cross_conditional_means_output.csv, long format
  date, outcome, conditioning, group, rel_mean, level, moment_error
where rel_mean is relative to the outcome's per-HH mean and level multiplies by
<outcome>_per_hh from aggregate_anchors.csv (when available).

Dependencies: numpy, pandas.
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, HERE)
from reconstruct import Reconstruction, Q, I  # noqa: E402

MEASURES = ["consum", "income", "wealth"]      # axis order of the copula tensor
K = 12
DENSE_TOP = [0.0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 0.999, 1.0]
GL_X, GL_W = np.polynomial.legendre.leggauss(12)


# ----------------------------------------------------------------------------
# Grid (Gauss-Legendre per segment, dense at the top, aligned with group edges)
# ----------------------------------------------------------------------------
def make_grid(group_edges):
    edges = np.unique(np.round(np.r_[DENSE_TOP, group_edges], 12))
    x, w = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        x.append((GL_X + 1) / 2 * (b - a) + a)
        w.append(GL_W * (b - a) / 2)
    x, w = np.concatenate(x), np.concatenate(w)
    Qx = np.stack([Q(k, x) for k in range(K)], 0)                    # K x n
    P = (Qx[:, None, :] * Qx[None, :, :]).reshape(K * K, -1)          # products Q_k Q_k'
    return dict(x=x, w=w, Qx=Qx, P=P, W2=np.outer(w, w))


# ----------------------------------------------------------------------------
# 2D maximum entropy
# ----------------------------------------------------------------------------
def maxent_density_2d(M, g, lam0=None, eps=1e-6, tol=1e-7, maxit=200):
    """Masses f (n x n) of exp(sum lam_kl Q_k(u) Q_l(v)) matching M (K x K, M[0,0]=1).
    Minimises the regularised dual  log Z(lam) - lam.M + eps/2 |lam|^2  by damped
    Newton (Armijo). Everything is separable: log f = Qx' Lam Qx, moments = Qx f Qx',
    Hessian = reshape(P f P'). Returns (f, lam, max moment error)."""
    free = np.ones((K, K), bool); free[0, 0] = False
    fi = free.ravel()
    lam = np.zeros((K, K)) if lam0 is None else lam0.copy()

    def state(L):
        e = g["Qx"].T @ L @ g["Qx"]
        mx = e.max()
        ew = np.exp(e - mx) * g["W2"]
        Z = ew.sum()
        f = ew / Z
        mu = g["Qx"] @ f @ g["Qx"].T
        D = np.log(Z) + mx - (L * M)[free].sum() + eps / 2 * (L[free] ** 2).sum()
        return f, mu, D

    f, mu, D = state(lam)
    grad = (mu - M)[free] + eps * lam[free]
    for _ in range(maxit):
        if np.abs(grad).max() < tol:
            break
        S4 = (g["P"] @ f @ g["P"].T).reshape(K, K, K, K)              # (k, k', l, l')
        H = S4.transpose(0, 2, 1, 3).reshape(K * K, K * K) - np.outer(mu.ravel(), mu.ravel())
        Hs = H[np.ix_(fi, fi)] + eps * np.eye(fi.sum())
        step = np.linalg.solve(Hs, grad)
        t, accepted = 1.0, False
        while t > 1e-8:
            Ln = lam.copy(); Ln[free] -= t * step
            fn, mun, Dn = state(Ln)
            if Dn <= D - 1e-4 * t * grad @ step:
                lam, f, mu, D, accepted = Ln, fn, mun, Dn, True
                break
            t /= 2
        if not accepted:
            break
        grad = (mu - M)[free] + eps * lam[free]
    err = np.abs(mu - M); err[0, 0] = 0.0
    return f, lam, float(err.max())


# ----------------------------------------------------------------------------
# Group means for one date
# ----------------------------------------------------------------------------
def pair_moments(kappa, cond, out):
    """K x K moments with rows = conditioning degree, cols = outcome degree."""
    ic, io = MEASURES.index(cond), MEASURES.index(out)
    idx = [0, 0, 0]
    idx[ic], idx[io] = slice(None), slice(None)
    S = kappa[tuple(idx)]
    return S if ic < io else S.T


def exact_mean(r, date, m):
    """∫_0^1 V(u) du with V the quantile function relative to the per-HH mean."""
    s = 0.0
    for a, b in zip(DENSE_TOP[:-1], DENSE_TOP[1:]):
        x = (GL_X + 1) / 2 * (b - a) + a
        s += (GL_W * (b - a) / 2 * r.quantile_at(date, m, x)).sum()
    return s


def group_means(r, date, cond, out, edges, g, method, warm):
    kappa = r._kappa(date)
    M = pair_moments(kappa, cond, out)
    labels = group_labels(edges)
    if method == "series":
        A = np.zeros(K)                                                # ∫ Q_l V
        for a, b in zip(DENSE_TOP[:-1], DENSE_TOP[1:]):
            x = (GL_X + 1) / 2 * (b - a) + a
            vals = r.quantile_at(date, out, x)
            A += np.array([(GL_W * (b - a) / 2 * Q(l, x) * vals).sum() for l in range(K)])
        res = []
        for lab, lo, hi in zip(labels, edges[:-1], edges[1:]):
            B = np.array([I(k, hi) - I(k, lo) for k in range(K)])
            res.append((lab, float(B @ M @ A / (hi - lo)), np.nan))
        return res
    f, lam, err = maxent_density_2d(M, g, lam0=warm.get((cond, out)))
    warm[(cond, out)] = lam
    V = r.quantile_at(date, out, g["x"])
    EV = f @ V                                                         # per conditioning node
    scale = exact_mean(r, date, out) / EV.sum()                       # groups add up to the exact mean
    res = []
    for lab, lo, hi in zip(labels, edges[:-1], edges[1:]):
        rows = (g["x"] > lo) & (g["x"] <= hi)
        res.append((lab, float(scale * EV[rows].sum() / (hi - lo)), err))
    return res


def group_labels(edges):
    n = len(edges) - 1
    out = []
    for i in range(n):
        width = round(100 * (edges[i + 1] - edges[i]))
        out.append(("bot" if i == 0 else "top" if i == n - 1 else "mid") + str(width))
    if len(set(out)) < n:  # e.g. quintiles: name groups by percentile range, p0-20 ... p80-100
        out = [f"p{round(100 * edges[i])}-{round(100 * edges[i + 1])}" for i in range(n)]
    return out


def has_measure(r, date, m):
    try:
        v = r.quantile_at(date, m, [0.5])
        return bool(np.isfinite(v).all())
    except Exception:
        return False


def load_anchors():
    p = os.path.join(REPO, "data", "aggregate_anchors.csv")
    if not os.path.exists(p):
        return None
    a = pd.read_csv(p)
    t = pd.to_datetime(a["time"].astype(str), errors="coerce")
    a["time"] = [f"{d.year}-Q{(d.month - 1) // 3 + 1}" if pd.notna(d) else s for d, s in zip(t, a["time"].astype(str))]
    return a.set_index("time")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dataset", default="PSID", choices=["PSID", "SCF", "CEX"])
    ap.add_argument("--kind", default="normal", choices=["normal", "average"])
    ap.add_argument("--date", default="2019-Q4", help="YYYY-Qq or 'all'")
    ap.add_argument("--groups", default="50,40,10", help="group sizes in percent, e.g. 50,40,10 or 40,40,20")
    ap.add_argument("--method", default="maxent", choices=["maxent", "series"])
    ap.add_argument("--out", default=os.path.join(HERE, "cross_conditional_means_output.csv"))
    args = ap.parse_args()

    sizes = np.array([float(s) for s in args.groups.split(",")]) / 100
    if abs(sizes.sum() - 1) > 1e-9:
        sys.exit("--groups must add up to 100")
    edges = np.r_[0.0, np.cumsum(sizes)]; edges[-1] = 1.0
    r = Reconstruction(os.path.join(REPO, "data", f"{args.dataset}_coefficients_{args.kind}.csv"))
    dates = r.available_dates() if args.date == "all" else [args.date]
    anchors = load_anchors()
    g = make_grid(edges)
    warm, rows = {}, []
    for date in dates:
        present = [m for m in MEASURES if has_measure(r, date, m)]
        for cond in present:
            for out in present:
                if cond == out:
                    continue
                for lab, rel, err in group_means(r, date, cond, out, edges, g, args.method, warm):
                    lvl = np.nan
                    col = f"{out}_per_hh"
                    if anchors is not None and col in anchors.columns and date in anchors.index:
                        lvl = rel * float(anchors.loc[date, col])
                    rows.append((date, out, cond, lab, rel, lvl, err))
    df = pd.DataFrame(rows, columns=["date", "outcome", "conditioning", "group", "rel_mean", "level", "moment_error"])
    df.to_csv(args.out, index=False)
    if len(dates) == 1:
        print(df.pivot_table(index=["outcome", "conditioning"], columns="group", values="rel_mean").round(3))
    print(f"wrote {args.out} ({len(df)} rows, method={args.method})")
    if args.method == "maxent" and len(df):
        print(f"max moment error of the max-entropy fits: {df['moment_error'].max():.2e}")


if __name__ == "__main__":
    main()
