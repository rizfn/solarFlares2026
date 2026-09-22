# The mean-field law with a running coupling.
#
# The rate equation fixes the spot exponent from one number, the branching ratio:
#
#     cos(pi tau_m) = (1 - q) / q .
#
# On the lattice q is not a number. The fraction of collisions that annihilate,
# 1 - q(s), varies over four decades (it is the U-shaped curve of plot_crossover).
# If the sub-leading balance is solved locally -- legitimate while q varies slowly
# on a log scale -- the lattice should carry a scale-dependent exponent
#
#     tau_loc(s) = 2 - arccos[(1 - q(s)) / q(s)] / pi ,
#
# and everything else follows from that one measured function:
#
#     P_coll(s) ~ s^(1 - 2 tau_loc)          min(i,j) over all collisions
#     P_emis(s)  =  P_coll(s) * [1 - q(s)]   the annihilating subset
#
# No fitted parameters anywhere: q(s) goes in, three spectra come out.
#
# It does not work well enough to build on. The mean residual tau_pred - tau_meas is
# -0.004 at p=0.5 and -0.018 at p=0.6, but climbs to +0.058 at p=0.9 and +0.062 at p=1:
# above the transition the running form predicts an upturn in the local exponent at
# large s that the spot spectrum does not show, and it is pinned near 3/2 there anyway,
# so the test has little leverage. At p=0.5 the predicted variation is about twice the
# measured one, in the other direction. Panel (d) still rebuilds the bent emission
# spectrum over eight decades, but that is mostly the (1 - q) factor, which is measured
# rather than predicted.
import os
import numpy as np
import matplotlib.pyplot as plt

DIR = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(DIR, "outputs")
PLOTS = os.path.join(DIR, "plots", "runningQ")

plt.rcParams.update({"font.size": 18, "axes.labelsize": 22, "xtick.labelsize": 16,
                     "ytick.labelsize": 16, "legend.fontsize": 14})

PERDEC = 4          # log bins per decade; coarse, because these are derivatives
MINC = 300          # a bin needs this many counts in every spectrum to be used
SMIN = 10.0         # below this the injection scale dominates and no continuum applies

COL = {0.3: "#7570b3", 0.5: "#1b9e77", 0.6: "#66a61e", 0.7: "#e6ab02",
       0.8: "#d95f02", 0.9: "#e7298a", 1.0: "#a6761d"}


def tag(L, rho, p):
    return "L_%d_rho_%g_p_%g" % (L, rho, p)


def load_hist(kind, L, rho, p, seeds=range(1, 9)):
    # sum the sparse (size, count) histograms over seeds
    acc = {}
    for s in seeds:
        f = os.path.join(OUT, "%s_%s_seed_%d.tsv" % (kind, tag(L, rho, p), s))
        if not os.path.exists(f) or os.path.getsize(f) < 20:
            continue
        d = np.loadtxt(f)
        if d.ndim == 1:
            d = d[None, :]
        for k, v in d:
            acc[k] = acc.get(k, 0.0) + v
    if not acc:
        raise OSError("no data for %s %s" % (kind, tag(L, rho, p)))
    k = np.array(sorted(acc))
    return k, np.array([acc[x] for x in k])


def tau_closed(q):
    # cos(pi tau) = (1 - q) / q, the root in 1 < tau < 2; undefined for q < 1/2
    a = (1.0 - q) / q
    ok = (a >= -1.0) & (a <= 1.0)
    return np.where(ok, 2.0 - np.arccos(np.clip(a, -1.0, 1.0)) / np.pi, np.nan)


def curves(L, rho, p):
    """The three spectra on one shared log grid, plus q(s) and the local slopes.

    Sharing the grid is what makes he/hc equal to 1 - q(s): the fraction of
    collisions at min(i,j) = s that annihilated rather than merged.
    """
    sc, cc = load_hist("chCollMin", L, rho, p)
    se, ce = load_hist("chEmisAll", L, rho, p)
    sm, cm = load_hist("chSpot", L, rho, p)
    top = max(sc.max(), se.max(), sm.max())
    e = np.geomspace(1, top + 1, int(PERDEC * np.log10(top)) + 2)
    hc = np.histogram(sc, e, weights=cc)[0]
    he = np.histogram(se, e, weights=ce)[0]
    hm = np.histogram(sm, e, weights=cm)[0]
    x, w = np.sqrt(e[:-1] * e[1:]), np.diff(e)
    k = (hc >= MINC) & (he >= MINC) & (hm >= MINC)
    x, lx = x[k], np.log(np.sqrt(e[:-1] * e[1:])[k])
    q = 1.0 - he[k] / hc[k]
    slope = lambda h: -np.gradient(np.log(h[k] / w[k]), lx)
    return dict(x=x, q=q, tau=tau_closed(q),
                spot=slope(hm), coll=slope(hc), emis=slope(he),
                Pemis=he[k] / w[k] / ce.sum())


def reconstruct(c):
    """Rebuild P_emis(s) from q(s) alone, up to one overall amplitude.

    ln P = -int (2 tau_loc - 1) dln s + ln(1 - q), integrated from the first
    usable bin. The amplitude is set by matching at that bin, so the comparison
    is entirely about shape.
    """
    lx = np.log(c["x"])
    a = 2.0 * c["tau"] - 1.0
    ok = np.isfinite(a)
    if ok.sum() < 3:
        return None, None
    lx, a = lx[ok], a[ok]
    integ = np.concatenate([[0.0], np.cumsum(np.diff(lx) * 0.5 * (a[1:] + a[:-1]))])
    lp = -integ + np.log(1.0 - c["q"][ok])
    P = np.exp(lp - lp[0]) * c["Pemis"][ok][0]
    return c["x"][ok], P


# --------------------------------------------------------------------- the figure

def figure(L=1024, rho=0.2, ps=(0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0),
           pool=((1024, 0.2), (512, 0.2), (256, 0.2), (512, 0.4), (512, 0.6),
                 (512, 0.1), (1024, 0.4), (1024, 0.6))):
    dat = {}
    for p in ps:
        try:
            dat[p] = curves(L, rho, p)
        except OSError:
            pass

    fig, ax = plt.subplots(2, 2, figsize=(15, 12))
    (a, b), (cx, d) = ax

    # (a) predicted vs measured local spot exponent, against s
    for p, c in dat.items():
        k = c["x"] > SMIN
        a.plot(c["x"][k], c["tau"][k], color=COL[p], lw=2.5)
        a.plot(c["x"][k], c["spot"][k], color=COL[p], lw=0, marker="o", ms=4.5,
               alpha=0.8)
    a.plot([], [], "k-", lw=2.5, label=r"$\tau_{\rm loc}(s)$ from $q(s)$")
    a.plot([], [], "ko", ms=4.5, label="measured spot slope")
    a.set_xscale("log")
    a.set_xlabel(r"$m$")
    a.set_ylabel(r"local spot exponent $\tau_m$")
    a.set_ylim(1.35, 1.95)
    a.legend(loc="upper right")
    a.text(0.03, 0.06, "(a)", transform=a.transAxes, fontsize=22)

    # (b) the same thing parametrically, pooled over L and rho
    lo, hi = 1.40, 1.90
    for LL, rr in pool:
        for p in ps:
            try:
                c = curves(LL, rr, p)
            except OSError:
                continue
            k = (c["x"] > SMIN) & np.isfinite(c["tau"])
            b.plot(c["tau"][k], c["spot"][k], "o", ms=4, alpha=0.45,
                   color=COL.get(p, "#888888"), mew=0)
    b.plot([lo, hi], [lo, hi], "k--", lw=2)
    for p in ps:
        if p in COL:
            b.plot([], [], "o", ms=7, color=COL[p], label="$p=%g$" % p)
    b.set_xlim(lo, hi); b.set_ylim(lo, hi)
    b.set_xlabel(r"predicted $\tau_{\rm loc}(s)$")
    b.set_ylabel(r"measured local slope")
    b.legend(loc="upper left", ncol=2, handletextpad=0.2, columnspacing=0.8)
    b.text(0.94, 0.06, "(b)", transform=b.transAxes, fontsize=22)

    # (c) the collision spectrum: 2 tau_loc - 1, no free parameters either
    for p, c in dat.items():
        k = c["x"] > SMIN
        cx.plot(c["x"][k], 2 * c["tau"][k] - 1, color=COL[p], lw=2.5)
        cx.plot(c["x"][k], c["coll"][k], color=COL[p], lw=0, marker="o", ms=4.5,
                alpha=0.8)
    cx.plot([], [], "k-", lw=2.5, label=r"$2\tau_{\rm loc}-1$")
    cx.plot([], [], "ko", ms=4.5, label="measured collision slope")
    cx.set_xscale("log")
    cx.set_xlabel(r"$s$")
    cx.set_ylabel(r"local slope of $P_{\rm coll}(s)$")
    cx.set_ylim(1.85, 2.85)
    cx.legend(loc="upper right")
    cx.text(0.03, 0.06, "(c)", transform=cx.transAxes, fontsize=22)

    # (d) the whole emission spectrum rebuilt from q(s), amplitude matched once
    for p, c in dat.items():
        k = c["x"] > SMIN
        cc = {kk: vv[k] if isinstance(vv, np.ndarray) else vv for kk, vv in c.items()}
        d.plot(cc["x"], cc["Pemis"], color=COL[p], lw=0, marker="o", ms=4.5,
               alpha=0.8)
        xr, Pr = reconstruct(cc)
        if xr is not None:
            d.plot(xr, Pr, color=COL[p], lw=2.5)
    d.plot([], [], "k-", lw=2.5, label="rebuilt from $q(s)$")
    d.plot([], [], "ko", ms=4.5, label="measured")
    d.set_xscale("log"); d.set_yscale("log")
    d.set_xlabel(r"$s$")
    d.set_ylabel(r"$P_{\rm emis}(s)$")
    d.legend(loc="upper right")
    d.text(0.03, 0.06, "(d)", transform=d.transAxes, fontsize=22)

    for k in (a, b, cx, d):
        k.grid(True, which="major", ls="-", lw=0.6, alpha=0.35)

    fig.tight_layout()
    os.makedirs(PLOTS, exist_ok=True)
    f = os.path.join(PLOTS, "runningQ_L_%d_rho_%g.png" % (L, rho))
    fig.savefig(f, dpi=300, bbox_inches="tight")
    print("wrote", f)


def table(L=1024, rho=0.2, ps=(0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)):
    # residual of the parametric collapse, decade by decade
    print("  p    bins   <pred-meas>   rms    range of tau_loc")
    for p in ps:
        try:
            c = curves(L, rho, p)
        except OSError:
            continue
        k = (c["x"] > SMIN) & np.isfinite(c["tau"])
        r = c["tau"][k] - c["spot"][k]
        print("%4g %6d %10.3f %8.3f    %.3f - %.3f"
              % (p, k.sum(), r.mean(), r.std(), c["tau"][k].min(), c["tau"][k].max()))


if __name__ == "__main__":
    table()
    figure()
