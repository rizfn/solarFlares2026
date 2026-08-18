# The closed form against the well-mixed Monte Carlo.
#
# The sub-leading balance of the rate equation (renormalization.md S5) fixes the spot
# exponent from the branching ratio alone:  cos(pi tau_m) = (1 - q) / q,  on the branch
# tau_m in [3/2, 2]. It has no solution for q < 1/2. wellMixed.cpp runs the same process
# as an exact stochastic simulation with q imposed as a rule, so the two can be compared
# with no fitted parameter at all.
import os
import glob
import numpy as np
import matplotlib.pyplot as plt

DIR = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(DIR, "outputs")
PLOTS = os.path.join(DIR, "plots", "meanField")

plt.rcParams.update({"font.size": 18, "axes.labelsize": 22, "xtick.labelsize": 16,
                     "ytick.labelsize": 16, "legend.fontsize": 15})

WIN_LO, WIN_FTOP = 60.0, 3e-4   # scaling-window fit, as in plots.py


def tau_closed(q):
    r = (1 - np.asarray(q, float)) / np.asarray(q, float)
    return np.where(r <= 1, 2 - np.arccos(np.clip(r, -1, 1)) / np.pi, np.nan)


def window_tau(sizes, counts, lo=WIN_LO, ftop=WIN_FTOP, perdec=6, minc=20):
    # Log-log slope of the binned density over the scaling window.
    #
    # Not a Hill estimator: that averages everything above xmin, so it is biased up by
    # both the pre-asymptotic small-m region and the cutoff, and the bias survives any
    # amount of statistics (checked in plot_tauVsQ.py against cos(pi tau) = (1-q)/q).
    #
    # The upper edge is a quantile of the distribution, not the largest mass observed:
    # sizes.max() is a single extreme value that jumps around between parameter points,
    # and letting it set the window made tau(p) visibly jagged. Bins below minc counts
    # are dropped and the fit is weighted by counts, since var(log y) ~ 1/counts.
    # Tuned on the well-mixed data where the exponent is known: this leaves a residual
    # bias of 0.002 against the closed form.
    tot = counts.sum()
    cc = np.cumsum(counts[::-1])[::-1] / tot          # P(M >= m)
    hi = max(float(sizes[cc >= ftop].max()), 10 * lo)
    nb = max(5, int(round(perdec * np.log10(hi / lo))) + 1)
    edges = np.geomspace(lo, hi, nb)
    h, _ = np.histogram(sizes, bins=edges, weights=counts)
    x = np.sqrt(edges[:-1] * edges[1:])
    y = h / np.diff(edges)
    m = h >= minc
    if m.sum() < 4:
        return np.nan
    return -np.polyfit(np.log(x[m]), np.log(y[m]), 1, w=np.sqrt(h[m]))[0]


def per_seed_tau(q, N=200000, kind="wmSpotSize"):
    # one exponent per seed, so the spread across seeds gives an honest error bar
    out = []
    for f in sorted(glob.glob(os.path.join(OUT, f"{kind}_N_{N}_q_{q:g}_seed_*.tsv"))):
        if os.path.getsize(f) < 20:
            continue
        d = np.loadtxt(f, dtype=np.int64, ndmin=2)
        t = window_tau(d[:, 0].astype(float), d[:, 1].astype(float))
        if np.isfinite(t):
            out.append(t)
    return np.array(out)


def plot_tau_vs_q(N=200000):
    qs = sorted({float(f.split("_q_")[1].split("_seed")[0])
                 for f in glob.glob(os.path.join(OUT, f"wmSpotSize_N_{N}_q_*.tsv"))})
    x, y, e, n = [], [], [], []
    for q in qs:
        t = per_seed_tau(q, N)
        if len(t):
            x.append(q); y.append(t.mean()); e.append(t.std() / np.sqrt(len(t)))
            n.append(len(t))
    x, y, e = np.array(x), np.array(y), np.array(e)

    fig, ax = plt.subplots(figsize=(8.5, 6.5))
    # the closed form exists only for q >= 1/2; shade where it does not
    ax.axvspan(min(x.min(), 0.4) - 0.02, 0.5, color="#bdbdbd", alpha=0.35, lw=0)
    ax.text(0.5 - 0.012, 1.55, "no solution", rotation=90, ha="right", va="bottom",
            color="#606060", fontsize=15)

    qc = np.linspace(0.5, 0.999, 400)
    ax.plot(qc, tau_closed(qc), "-", color="#1a9850", lw=3.5, alpha=0.9,
            label=r"$\cos\pi\tau_m=\dfrac{1-q}{q}$", zorder=2)
    ok = x >= 0.5
    ax.errorbar(x[ok], y[ok], yerr=e[ok], fmt="o", color="k", ms=8, capsize=3,
                lw=1.5, zorder=3, label=f"well-mixed MC ({min(n)} seeds)")
    ax.errorbar(x[~ok], y[~ok], yerr=e[~ok], fmt="o", mfc="white", color="k", ms=8,
                capsize=3, lw=1.5, zorder=3)
    ax.plot(0.5, 2.0, "*", ms=20, color="#1a9850", mec="k", mew=0.8, zorder=4)

    ax.set_xlabel("branching ratio $q$")
    ax.set_ylabel(r"spot exponent $\tau_m$")
    ax.set_xlim(min(x.min(), 0.4) - 0.02, 1.005)
    ax.set_ylim(1.4, 2.9)
    ax.legend(frameon=False, loc="upper right")
    ax.grid(True, which="major", ls=":", lw=0.9, alpha=0.6)
    ax.set_axisbelow(True)
    fig.tight_layout()
    os.makedirs(PLOTS, exist_ok=True)
    fig.savefig(os.path.join(PLOTS, f"tauVsQ_N_{N}.png"), dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    plot_tau_vs_q()
