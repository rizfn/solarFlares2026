# Is the small-s branch an injection artefact?
#
# probabilisticNeighbour has a U-shaped annihilating fraction 1 - q(s): falling at small
# s, rising at large s, with the minimum at the crossover scale s*. The falling branch is
# attributed to the (1-p) channel dropping lone minority monomers into wrong-sign domains.
# This model replaces that channel with a bipole -- the + and - land adjacent to each
# other -- and changes nothing else. If the attribution holds, the falling branch should
# weaken or vanish while the rising branch and the sliding exponent survive.
#
# It does not. The falling branch survives and steepens, because a bipole self-annihilates
# at s = 1 and supplies its own small-s excess. The useful result is the other one: at
# matched q, the two rules give the same tau_m, the same s* and the same psi, so the
# injection rule enters the cascade only by setting q. plot_collapse_vs_q is that figure.
#
# Every figure overlays the two models at matched p, L and rho.
import os
import numpy as np
import matplotlib.pyplot as plt

DIR = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(DIR, "outputs")
REF = os.path.join(DIR, "..", "probabilisticNeighbour", "outputs")
PLOTS = os.path.join(DIR, "plots")

plt.rcParams.update({"font.size": 18, "axes.labelsize": 22, "xtick.labelsize": 16,
                     "ytick.labelsize": 16, "legend.fontsize": 14})

PERDEC, MINC = 6, 400
COL = {0.0: "#8c8c8c", 0.3: "#7570b3", 0.5: "#1b9e77", 0.6: "#66a61e",
       0.7: "#e6ab02", 0.8: "#d95f02", 0.9: "#e7298a", 0.95: "#ce1256",
       1.0: "#a6761d"}


def load_hist(kind, L, rho, p, seeds=range(1, 9), root=None):
    # sum the sparse (size, count) histograms over seeds
    root = root or OUT
    acc = {}
    for s in seeds:
        f = os.path.join(root, "%s_L_%d_rho_%g_p_%g_seed_%d.tsv" % (kind, L, rho, p, s))
        if not os.path.exists(f) or os.path.getsize(f) < 20:
            continue
        d = np.loadtxt(f)
        if d.ndim == 1:
            d = d[None, :]
        for k, v in d:
            acc[k] = acc.get(k, 0.0) + v
    if not acc:
        raise OSError("no data for %s L=%d rho=%g p=%g" % (kind, L, rho, p))
    k = np.array(sorted(acc))
    return k, np.array([acc[x] for x in k])


def spectra(L, rho, p, bipole=True):
    """P_emis(s), P_coll(s) and 1 - q(s), all on one shared set of log bins.

    The ratio is only 1 - q(s) if the two spectra are binned identically and left
    unnormalised, so it is built here rather than from separately loaded curves.
    """
    pre, root = ("cb", OUT) if bipole else ("ch", REF)
    sc, cc = load_hist(pre + "CollMin", L, rho, p, root=root)
    se, ce = load_hist(pre + "EmisAll", L, rho, p, root=root)
    top = max(sc.max(), se.max())
    e = np.geomspace(1, top + 1, int(PERDEC * np.log10(top)) + 2)
    hc = np.histogram(sc, e, weights=cc)[0]
    he = np.histogram(se, e, weights=ce)[0]
    x, w = np.sqrt(e[:-1] * e[1:]), np.diff(e)
    k = hc >= MINC
    return x[k], he[k] / w[k] / ce.sum(), hc[k] / w[k] / cc.sum(), he[k] / hc[k]


def crossover_scale(x, r):
    # s* = the minimum of 1 - q(s), on a smoothed curve so one noisy bin cannot claim
    # it. nan if the minimum sits at an edge, i.e. the U is not resolved, which happens
    # only at p = 1 where there is no U in either model.
    if len(r) < 7:
        return np.nan
    lr = np.log(r)
    sm = np.convolve(lr, np.ones(3) / 3, mode="same")
    sm[0], sm[-1] = lr[0], lr[-1]
    j = int(np.argmin(sm))
    return np.nan if j in (0, len(r) - 1) else x[j]


def hill(sizes, counts, lo=10.0):
    # Hill estimator above lo, on the raw sparse histogram
    k = sizes >= lo
    s, c = sizes[k], counts[k]
    if c.sum() < 100:
        return np.nan
    return 1.0 + c.sum() / np.sum(c * np.log(s / (lo - 0.5)))


# ------------------------------------------------------------------ the decisive figure

def plot_qOfS(L=1024, rho=0.2, ps=(0.5, 0.6, 0.7, 0.8, 0.9)):
    fig, ax = plt.subplots(1, 2, figsize=(15, 6), sharey=True)
    for k, (bip, name) in enumerate([(False, "monomer injection"), (True, "bipole injection")]):
        for p in ps:
            try:
                x, _, _, r = spectra(L, rho, p, bipole=bip)
            except OSError:
                continue
            ax[k].plot(x, r, "-o", ms=3.5, lw=2, color=COL[p], label="$p=%g$" % p)
            ss = crossover_scale(x, r)
            if np.isfinite(ss):
                ax[k].plot([ss], [r[np.argmin(np.abs(x - ss))]], "v", ms=14,
                           color=COL[p], mec="k", mew=1.2, zorder=5)
        ax[k].set_xscale("log"); ax[k].set_yscale("log")
        ax[k].set_xlabel("$s$")
        ax[k].set_title(name, fontsize=20)
        ax[k].grid(True, which="major", ls="-", lw=0.6, alpha=0.35)
    ax[0].set_ylabel("$1 - q(s)$")
    ax[1].legend(loc="upper left", ncol=2)
    fig.tight_layout()
    os.makedirs(PLOTS, exist_ok=True)
    f = os.path.join(PLOTS, "qOfS_L_%d_rho_%g.png" % (L, rho))
    fig.savefig(f, dpi=300, bbox_inches="tight")
    print("wrote", f)


def plot_crossover_vs_p(rho=0.2, Ls=(256, 512, 1024),
                        ps=(0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95)):
    fig, a = plt.subplots(figsize=(9, 6.5))
    for L, mk in zip(Ls, ("^", "o", "s")):
        for bip, ls, lab in [(False, "--", "monomer"), (True, "-", "bipole")]:
            xs, ys = [], []
            for p in ps:
                try:
                    x, _, _, r = spectra(L, rho, p, bipole=bip)
                except OSError:
                    continue
                ss = crossover_scale(x, r)
                if np.isfinite(ss):
                    xs.append(p); ys.append(ss)
            if xs:
                a.plot(xs, ys, ls, marker=mk, ms=9, lw=2,
                       label="%s, $L=%d$" % (lab, L))
    a.set_yscale("log")
    a.set_xlabel("$p$"); a.set_ylabel("$s^*$")
    a.legend(fontsize=12, ncol=2)
    a.grid(True, which="major", ls="-", lw=0.6, alpha=0.35)
    fig.tight_layout()
    os.makedirs(PLOTS, exist_ok=True)
    f = os.path.join(PLOTS, "crossoverVsP_rho_%g.png" % rho)
    fig.savefig(f, dpi=300, bbox_inches="tight")
    print("wrote", f)


def global_q(L, rho, p, bipole=True, lo=10.0):
    """Branching ratio over all collisions above the injection scale.

    Excluding s < lo matters: both injection rules dump events at s = 1, by
    different mechanisms, and including them measures the injection rather than
    the cascade.
    """
    pre, root = ("cb", OUT) if bipole else ("ch", REF)
    sc, cc = load_hist(pre + "CollMin", L, rho, p, root=root)
    se, ce = load_hist(pre + "EmisAll", L, rho, p, root=root)
    return 1.0 - ce[se >= lo].sum() / cc[sc >= lo].sum()


def tau_closed(q):
    a = (1.0 - q) / q
    return np.where(np.abs(a) <= 1, 2.0 - np.arccos(np.clip(a, -1, 1)) / np.pi, np.nan)


def plot_collapse_vs_q(L=256, rho=0.2,
                       ps=(0.5, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0)):
    """The point of the whole comparison.

    Against p the two models look different. Against the measured q they do not:
    the injection rule enters the cascade only by setting q, and both tau_m and
    the crossover scale are functions of q alone.
    """
    fig, ax = plt.subplots(1, 2, figsize=(15, 6))
    for bip, mk, lab in [(False, "o", "monomer"), (True, "s", "bipole")]:
        pre, root = ("cb", OUT) if bip else ("ch", REF)
        qs, tm, ss = [], [], []
        for p in ps:
            try:
                q = global_q(L, rho, p, bipole=bip)
                sm, cm = load_hist(pre + "Spot", L, rho, p, root=root)
                x, _, _, r = spectra(L, rho, p, bipole=bip)
            except OSError:
                continue
            qs.append(q); tm.append(hill(sm, cm)); ss.append(crossover_scale(x, r))
        ax[0].plot(qs, tm, mk, ms=11, mew=1.2, mec="k", label=lab)
        ax[1].plot(qs, ss, mk, ms=11, mew=1.2, mec="k", label=lab)
    qq = np.linspace(0.55, 0.995, 300)
    ax[0].plot(qq, tau_closed(qq), "-", lw=2.5, color="#2ca25f",
               label=r"$\cos\pi\tau_m=\frac{1-q}{q}$")
    ax[0].set_xlabel("$q$"); ax[0].set_ylabel(r"$\tau_m$")
    ax[1].set_yscale("log")
    ax[1].set_xlabel("$q$"); ax[1].set_ylabel("$s^*$")
    for k in ax:
        k.legend(); k.grid(True, which="major", ls="-", lw=0.6, alpha=0.35)
    fig.tight_layout()
    os.makedirs(PLOTS, exist_ok=True)
    f = os.path.join(PLOTS, "collapseVsQ_L_%d_rho_%g.png" % (L, rho))
    fig.savefig(f, dpi=300, bbox_inches="tight")
    print("wrote", f)


def plot_exponents_vs_p(L=1024, rho=0.2,
                        ps=(0.0, 0.3, 0.5, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0)):
    fig, a = plt.subplots(figsize=(9, 6.5))
    for bip, ls, lab in [(False, "--", "monomer"), (True, "-", "bipole")]:
        pre, root = ("cb", OUT) if bip else ("ch", REF)
        pm, tm, ts = [], [], []
        for p in ps:
            try:
                sm, cm = load_hist(pre + "Spot", L, rho, p, root=root)
                se, ce = load_hist(pre + "EmisAll", L, rho, p, root=root)
            except OSError:
                continue
            pm.append(p); tm.append(hill(sm, cm)); ts.append(hill(se, ce))
        if not pm:
            continue
        a.plot(pm, ts, ls, marker="o", ms=8, lw=2, color="#b2182b",
               label=r"$\tau_s$, %s" % lab)
        a.plot(pm, tm, ls, marker="^", ms=8, lw=2, color="#ef8a62",
               label=r"$\tau_m$, %s" % lab)
    a.set_xlabel("$p$"); a.set_ylabel(r"power-law exponent $\tau$")
    a.legend(fontsize=13, ncol=2)
    a.grid(True, which="major", ls="-", lw=0.6, alpha=0.35)
    fig.tight_layout()
    os.makedirs(PLOTS, exist_ok=True)
    f = os.path.join(PLOTS, "exponentsVsP_L_%d_rho_%g.png" % (L, rho))
    fig.savefig(f, dpi=300, bbox_inches="tight")
    print("wrote", f)


def plot_spectra(L=1024, rho=0.2, ps=(0.5, 0.7, 0.9)):
    fig, ax = plt.subplots(1, 2, figsize=(15, 6))
    for p in ps:
        for bip, ls, lab in [(False, "--", "monomer"), (True, "-", "bipole")]:
            try:
                x, pe, pc, _ = spectra(L, rho, p, bipole=bip)
            except OSError:
                continue
            ax[0].plot(x, pe, ls, lw=2, color=COL[p],
                       label="$p=%g$, %s" % (p, lab))
            ax[1].plot(x, pe * x ** 2, ls, lw=2, color=COL[p])
    ax[0].set_xscale("log"); ax[0].set_yscale("log")
    ax[0].set_xlabel("$s$"); ax[0].set_ylabel("$P_{\\rm emis}(s)$")
    ax[0].legend(fontsize=12)
    ax[1].set_xscale("log"); ax[1].set_yscale("log")
    ax[1].set_xlabel("$s$"); ax[1].set_ylabel("$s^2 P_{\\rm emis}(s)$")
    for k in ax:
        k.grid(True, which="major", ls="-", lw=0.6, alpha=0.35)
    fig.tight_layout()
    os.makedirs(PLOTS, exist_ok=True)
    f = os.path.join(PLOTS, "spectra_L_%d_rho_%g.png" % (L, rho))
    fig.savefig(f, dpi=300, bbox_inches="tight")
    print("wrote", f)


if __name__ == "__main__":
    plot_qOfS()
    plot_collapse_vs_q()
    plot_crossover_vs_p()
    plot_exponents_vs_p()
    plot_spectra()
