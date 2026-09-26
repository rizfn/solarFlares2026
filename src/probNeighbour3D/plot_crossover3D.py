# Does psi, the rising slope of the annihilation fraction 1 - q(s), depend on dimension?
#
# Same decomposition as ../probabilisticNeighbour/plot_crossover.py: the flare spectrum is
# the collision spectrum times the fraction of collisions that annihilate, and that
# fraction is U-shaped. In 2D the rising side has slope psi = 0.28 at every p, L and rho.
#
# In 3D it is the same number. The local slope sits on a plateau at 0.27-0.28 before a
# roll-off that starts at the same fraction of the system mass at L = 48 and 64, so the
# roll-off is finite size and the plateau is psi. What does change with dimension is the
# position of the U: s* is five to twenty times smaller in 3D, and the U survives at
# p = 0.5 and p = 1, where the 2D one is gone.
import os
import importlib.util
import numpy as np
import matplotlib.pyplot as plt

DIR = os.path.dirname(os.path.abspath(__file__))
OUT3 = os.path.join(DIR, "outputs")
OUT2 = os.path.join(DIR, "..", "probabilisticNeighbour", "outputs")
PLOTS = os.path.join(DIR, "plots", "crossover")

plt.rcParams.update({"font.size": 18, "axes.labelsize": 22, "xtick.labelsize": 16,
                     "ytick.labelsize": 16, "legend.fontsize": 14})

PERDEC, MINC = 6, 400
FIT_LO = 4.0      # psi window starts at FIT_LO * s*
FIT_TOP = 30.0    # 2D: ends at the largest collision / FIT_TOP, as in plot_mechanism
FIT_M = 0.02      # 3D: ends at FIT_M * rho L^3, below the finite-size roll-off
PS = (0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0)
COL = {0.5: "#1b9e77", 0.6: "#66a61e", 0.7: "#e6ab02", 0.75: "#e08214", 0.8: "#d95f02",
       0.85: "#c4451c", 0.9: "#e7298a", 0.95: "#a6329a", 1.0: "#1f78b4"}


def read(kind, root, L, rho, p, seed):
    f = os.path.join(root, "%s_L_%d_rho_%g_p_%g_seed_%d.tsv" % (kind, L, rho, p, seed))
    if not os.path.exists(f) or os.path.getsize(f) < 20:
        return None
    return np.loadtxt(f, ndmin=2)


def load_hist(kind, root, L, rho, p, seeds=range(1, 9)):
    # sum the sparse (size, count) histograms over seeds
    acc = {}
    for s in seeds:
        d = read(kind, root, L, rho, p, s)
        if d is None:
            continue
        for k, v in d:
            acc[k] = acc.get(k, 0.0) + v
    if not acc:
        raise OSError("no %s for L=%d rho=%g p=%g" % (kind, L, rho, p))
    k = np.array(sorted(acc))
    return k, np.array([acc[x] for x in k])


def seeds_of(root, L, rho, p):
    return [s for s in range(1, 9) if read("chEmisAll", root, L, rho, p, s) is not None]


def spectra(root, L, rho, p, seeds=range(1, 9), edges=None, minc=MINC):
    # 1 - q(s) is the ratio of the unnormalised flare and collision counts on shared bins
    sc, cc = load_hist("chCollMin", root, L, rho, p, seeds)
    se, ce = load_hist("chEmisAll", root, L, rho, p, seeds)
    if edges is None:
        top = max(sc.max(), se.max())
        edges = np.geomspace(1, top + 1, int(PERDEC * np.log10(top)) + 2)
    hc = np.histogram(sc, edges, weights=cc)[0]
    he = np.histogram(se, edges, weights=ce)[0]
    x = np.sqrt(edges[:-1] * edges[1:])
    k = hc >= minc
    return x[k], he[k] / np.maximum(hc[k], 1), edges, sc.max()


def crossover_scale(x, r):
    # minimum of 1 - q(s) on a 3-bin smoothed curve; nan if it sits at an edge
    if len(r) < 7:
        return np.nan
    lr = np.log(r)
    sm = np.convolve(lr, np.ones(3) / 3, mode="same")
    sm[0], sm[-1] = lr[0], lr[-1]
    j = int(np.argmin(sm))
    return np.nan if j in (0, len(r) - 1) else x[j]


def psi(root, L, rho, p, dim):
    """Rising slope of 1 - q(s) from FIT_LO s* to below the large-s roll-off. Window from
    the pooled data, slope refitted per seed on the same bins so the error is the seed
    spread. In 3D the roll-off sits at a fixed fraction of the total mass rho L^3 (panel
    c), well below the largest collision, so the 2D upper cut would fit into it."""
    x, r, e, top = spectra(root, L, rho, p)
    ss = crossover_scale(x, r)
    lo, hi = FIT_LO * ss, (FIT_M * rho * L ** 3 if dim == 3 else top / FIT_TOP)
    if not np.isfinite(ss) or hi < 5 * lo:   # under 0.7 decades is not a slope
        return np.nan, np.nan, ss, (lo, hi)
    fit = lambda x, r: np.polyfit(np.log(x[(x >= lo) & (x <= hi)]),
                                  np.log(r[(x >= lo) & (x <= hi)]), 1)[0]
    per = []
    for s in seeds_of(root, L, rho, p):
        xs, rs, _, _ = spectra(root, L, rho, p, seeds=[s], edges=e, minc=MINC // 4)
        per.append(fit(xs, rs))
    err = np.std(per, ddof=1) / np.sqrt(len(per)) if len(per) > 1 else np.nan
    return fit(x, r), err, ss, (lo, hi)


def plots2D():
    # the 2D plateau_fit, loaded by path since this directory has its own plots.py
    spec = importlib.util.spec_from_file_location(
        "plots2D", os.path.join(DIR, "..", "probabilisticNeighbour", "plots.py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def alpha_and_mf(root, L, rho, p, fit):
    # slope of the all-collision min(i,j) spectrum and the mean-field 2 tau_m - 1, both
    # from the same plateau estimator
    a = fit(*load_hist("chCollMin", root, L, rho, p))[0]
    return a, 2 * fit(*load_hist("chSpot", root, L, rho, p))[0] - 1


def local_slope(root, L, rho, p, perdec=4):
    # coarser bins than the fits use, since a derivative amplifies the noise
    x, r, _, top = spectra(root, L, rho, p)
    s = crossover_scale(x, r)
    e = np.geomspace(1, top + 1, int(perdec * np.log10(top)) + 2)
    x, r, _, _ = spectra(root, L, rho, p, edges=e)
    return x, np.gradient(np.log(r), np.log(x)), s


def figure(L3=64, L3b=48, L2=512, rho=0.2):
    os.makedirs(PLOTS, exist_ok=True)
    fig, ax = plt.subplots(2, 3, figsize=(23, 13))
    a, b, c, d, e, f = ax.ravel()
    runs = ((OUT3, L3, 3, "o", "#b2182b", "3D, $L=%d$" % L3),
            (OUT3, L3b, 3, "^", "#ef8a62", "3D, $L=%d$" % L3b),
            (OUT2, L2, 2, "s", "#2166ac", "2D, $L=%d$" % L2))

    # (a) the annihilation fraction in 3D, with its minimum marked
    for p in PS:
        x, r, _, _ = spectra(OUT3, L3, rho, p)
        a.plot(x, r, "-o", color=COL[p], lw=2.2, ms=4, label="$p=%g$" % p)
        s = crossover_scale(x, r)
        if np.isfinite(s):
            a.plot([s], [r[np.argmin(np.abs(x - s))]], "v", color=COL[p], ms=13,
                   mec="k", mew=1.2, zorder=5)
    a.set(xscale="log", yscale="log", xlabel="$s$", ylabel="$1-q(s)$")
    a.legend(frameon=False, fontsize=12, ncol=3, loc="lower right")

    # (b) local slope against s / s*, 3D solid, 2D dashed: the rising side plateaus at
    # the same height in both; the 3D box just holds fewer decades above s*
    for p in (0.7, 0.8, 0.9):
        for root, L, ls, mk in ((OUT3, L3, "-", "o"), (OUT2, L2, "--", "s")):
            x, sl, s = local_slope(root, L, rho, p)
            k = x > s
            b.plot(x[k] / s, sl[k], ls, marker=mk, color=COL[p], lw=2.2, ms=5)
        b.plot([], [], "-", color=COL[p], lw=3, label="$p=%g$" % p)
    b.plot([], [], "k-o", lw=2.2, ms=5, label="3D, $L=%d$" % L3)
    b.plot([], [], "k--s", lw=2.2, ms=5, label="2D, $L=%d$" % L2)
    b.axhline(0.28, color="k", ls=":", lw=1.6)
    b.text(1.2, 0.295, r"$\psi_{2D}=0.28$", fontsize=15)
    b.set(xscale="log", ylim=(-0.05, 0.4), xlabel="$s/s^*$",
          ylabel=r"local slope $d\ln(1-q)/d\ln s$")
    b.legend(frameon=False, fontsize=12, ncol=2, loc="lower right")

    # (c) the 3D roll-off collapses on s / (rho L^3): it is finite size, and the psi
    # window stops at the grey line
    for p in (0.8, 1.0):
        for L, ls, mk in ((L3, "-", "o"), (L3b, "--", "^")):
            x, sl, s = local_slope(OUT3, L, rho, p)
            k = x > s
            c.plot(x[k] / (rho * L ** 3), sl[k], ls, marker=mk, color=COL[p], lw=2.2, ms=5)
        c.plot([], [], "-", color=COL[p], lw=3, label="$p=%g$" % p)
    c.plot([], [], "k-o", lw=2.2, ms=5, label="$L=%d$" % L3)
    c.plot([], [], "k--^", lw=2.2, ms=5, label="$L=%d$" % L3b)
    c.axvline(FIT_M, color="0.5", lw=2)
    c.axhline(0.28, color="k", ls=":", lw=1.6)
    c.set(xscale="log", ylim=(-0.05, 0.4), xlabel=r"$s/\rho L^3$",
          ylabel=r"local slope $d\ln(1-q)/d\ln s$")
    c.legend(frameon=False, fontsize=12, ncol=2, loc="lower left")

    # (d) psi, and (e) s*, against p in both dimensions
    table = {}
    for root, L, dim, mk, col, lab in runs:
        rows = []
        for p in PS:
            try:
                v, err, s, win = psi(root, L, rho, p, dim)
            except OSError:
                continue
            rows.append((p, v, err, s, win))
        table[lab] = rows
        pp, vv, ee, sv = (np.array([r[i] for r in rows]) for i in range(4))
        d.errorbar(pp, vv, yerr=ee, fmt=mk + "-", color=col, ms=9, lw=2, capsize=3,
                   label=lab)
        e.plot(pp, sv, mk + "-", color=col, ms=9, lw=2, label=lab)
    d.axhline(0.28, color="k", ls=":", lw=1.6)
    d.set(xlabel="$p$", ylabel=r"$\psi$", ylim=(0.1, 0.35))
    d.legend(frameon=False, loc="lower right")
    e.set(yscale="log", xlabel="$p$", ylabel="$s^*$")
    e.legend(frameon=False)

    # (f) the baseline: all-collision slope alpha against the mean-field 2 tau_m - 1
    fit = plots2D().plateau_fit
    for root, L, dim, mk, col, lab in runs:
        al, mf = [], []
        for p in PS:
            try:
                v = alpha_and_mf(root, L, rho, p, fit)
            except OSError:
                continue
            al.append(v[0]); mf.append(v[1])
        f.plot(mf, al, mk, color=col, ms=11, mec="k", mew=1, label=lab)
    lim = (1.9, 2.5)
    f.plot(lim, lim, "k--", lw=1.4)
    f.set(xlim=lim, ylim=lim, xlabel=r"$2\tau_m-1$", ylabel=r"$\alpha$ (all collisions)")
    f.legend(frameon=False, loc="lower right")

    for k, lab in zip(ax.ravel(), "abcdef"):
        k.text(0.02, 0.97, "(%s)" % lab, transform=k.transAxes, va="top", fontsize=22)
        k.grid(True, which="major", ls="-", lw=0.6, alpha=0.35)
        k.grid(True, which="minor", ls=":", lw=0.4, alpha=0.25)
        k.set_axisbelow(True)
    fig.tight_layout()
    fn = os.path.join(PLOTS, "crossover3D_L_%d_%d_vs_2D_L_%d_rho_%g.png" % (L3, L3b, L2, rho))
    fig.savefig(fn, dpi=300, bbox_inches="tight")
    print("wrote", fn)

    for lab, rows in table.items():
        print(lab)
        for p, v, err, s, (lo, hi) in rows:
            print("  p=%-5g psi=%.3f +- %.3f  s*=%-7.3g window %.3g-%.3g (%.2f dec)"
                  % (p, v, err, s, lo, hi, np.log10(hi / lo) if hi > lo else np.nan))


if __name__ == "__main__":
    figure()
