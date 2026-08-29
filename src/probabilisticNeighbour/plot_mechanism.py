# What sets the shape of the annihilation fraction.
#
# The flare spectrum is the collision spectrum times the fraction of collisions that
# annihilate. That fraction falls, turns, and rises again. Only the rising side is a
# straight line on log-log; the falling side is a hump, so it has no single slope.
#
# Counting a spot's collisions rather than its size explains the rising side: the
# chance of being destroyed climbs steadily with collisions survived, and size climbs
# with collisions too, so one slope divided by the other gives the rising slope.
import os
import glob
import numpy as np
import matplotlib.pyplot as plt

import plot_crossover as C

DIR = os.path.dirname(os.path.abspath(__file__))
PLOTS = os.path.join(DIR, "plots", "crossover")

plt.rcParams.update({"font.size": 18, "axes.labelsize": 22, "xtick.labelsize": 16,
                     "ytick.labelsize": 16, "legend.fontsize": 14})

PS = (0.7, 0.75, 0.8, 0.85, 0.9, 0.95)
COL = C.COL


def collisions(L, rho, p):
    """Pooled over seeds: collisions survived, the fraction of those collisions that
    annihilate, and the mean flare size produced at that count."""
    acc = {}
    for f in glob.glob(os.path.join(C.OUT, "chNcoll_%s_seed_*.tsv" % C.tag(L, rho, p))):
        d = np.loadtxt(f)
        if d.size == 0:
            continue
        if d.ndim == 1:
            d = d[None, :]
        for n, merge, ann, msize in d:
            k = round(np.log10(n), 4)
            acc.setdefault(k, np.zeros(4))
            acc[k] += [merge, ann, msize * ann, ann]
    ks = np.array(sorted(acc))
    v = np.array([acc[k] for k in ks])
    tot = v[:, 0] + v[:, 1]
    keep = tot > 2000
    return (10 ** ks[keep], (v[:, 1] / tot)[keep],
            (v[:, 2] / np.maximum(v[:, 3], 1))[keep])


def turning_point(x, y):
    ly = np.log(y)
    sm = np.convolve(ly, np.ones(3) / 3, mode="same")
    sm[0], sm[-1] = ly[0], ly[-1]
    j = int(np.argmin(sm))
    return j if 3 <= j <= len(y) - 4 else None


def figure(L=512, rho=0.2):
    os.makedirs(PLOTS, exist_ok=True)
    fig, ax = plt.subplots(2, 2, figsize=(16, 13))
    a, b, c, d = ax.ravel()

    # (a) slope measured point by point along the annihilation fraction. Flat means a
    # power law. The rising side is flat; the falling side is a hump with no one slope.
    for p in (0.8, 0.9):
        x, _, _, r = C.spectra(L, rho, p)
        j = list(x).index(C.crossover_scale(x, r))
        sl = np.gradient(np.log(r), np.log(x))
        k = x > 3
        a.plot(x[k], sl[k], "o-", color=COL[p], lw=2.2, ms=5, label="$p=%g$" % p)
        a.plot([x[j]], [sl[j]], "v", color=COL[p], ms=14, mec="k", mew=1.2)
    a.axhline(0, color="k", lw=0.8)
    a.axhline(0.28, color="k", ls="--", lw=1.4)
    a.text(2e4, 0.33, "0.28", fontsize=15)
    a.set(xscale="log", xlabel="flare size $s$", ylabel="slope of the annihilation fraction")
    a.legend(frameon=False, loc="lower right")

    # (b) the same fraction against collisions survived instead of size
    for p in PS:
        n, frac, _ = collisions(L, rho, p)
        b.plot(n, frac, "-o", color=COL[p], lw=2.2, ms=4, label="$p=%g$" % p)
        j = turning_point(n, frac)
        if j is not None:
            b.plot([n[j]], [frac[j]], "v", color=COL[p], ms=13, mec="k", mew=1.2)
    b.set(xscale="log", yscale="log", xlabel="collisions survived",
          ylabel="fraction that annihilate")
    b.legend(frameon=False, fontsize=13, ncol=2)

    # (c) how big a spot is by the time it has survived that many collisions
    for p in PS:
        n, _, ms = collisions(L, rho, p)
        k = ms > 1
        c.plot(n[k], ms[k], "-", color=COL[p], lw=2.4, label="$p=%g$" % p)
    c.set(xscale="log", yscale="log", xlabel="collisions survived",
          ylabel="mean flare size produced")
    c.legend(frameon=False, fontsize=13, ncol=2, loc="lower right")

    # (d) the two rates divided against the measured rising slope
    px, py = [], []
    for p in PS:
        n, frac, ms = collisions(L, rho, p)
        j = turning_point(n, frac)
        if j is None:
            continue
        rise = np.polyfit(np.log(n[j:]), np.log(frac[j:]), 1)[0]
        kz = (n > 3) & (ms > 3)
        growth = np.polyfit(np.log(n[kz]), np.log(ms[kz]), 1)[0]
        x, _, _, r = C.spectra(L, rho, p)
        i = list(x).index(C.crossover_scale(x, r))
        meas = np.polyfit(np.log(x[i:-3]), np.log(r[i:-3]), 1)[0]
        px.append(rise / growth); py.append(meas)
        d.plot(rise / growth, meas, "o", color=COL[p], ms=16, mec="k", mew=1.2,
               label="$p=%g$" % p)
    lim = (0.18, 0.32)
    d.plot(lim, lim, "k--", lw=1.4)
    d.set(xlim=lim, ylim=lim,
          xlabel="predicted: rise / growth, per collision",
          ylabel="measured rising slope")
    d.legend(frameon=False, fontsize=13, loc="lower right")

    for k, lab in zip((a, b, c, d), "abcd"):
        k.text(0.02, 0.97, "(%s)" % lab, transform=k.transAxes, va="top", fontsize=22)
        k.grid(True, which="major", ls="-", lw=0.6, alpha=0.35)
        k.grid(True, which="minor", ls=":", lw=0.4, alpha=0.25)
        k.set_axisbelow(True)
    fig.tight_layout()
    f = os.path.join(PLOTS, "mechanism_L_%d_rho_%g.png" % (L, rho))
    fig.savefig(f, dpi=300, bbox_inches="tight")
    print("wrote", f)


if __name__ == "__main__":
    figure()
