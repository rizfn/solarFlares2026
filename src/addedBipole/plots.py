import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

DIR = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(DIR, "outputs")
PLOTS = os.path.join(DIR, "plots")

plt.rcParams.update({"font.size": 18, "axes.labelsize": 22, "xtick.labelsize": 16,
                     "ytick.labelsize": 16, "legend.fontsize": 15})

XMIN = 10   # lower cutoff of the power-law fits, in flux units


def g(v):
    return f"{v:g}"

def grid(ax):
    ax.grid(True, which="major", ls=":", lw=0.9, alpha=0.6)
    ax.set_axisbelow(True)

def tag(L, rho, r):
    return f"L_{g(L)}_rho_{g(rho)}_r_{g(r)}"

def files(kind, L, rho, r):
    fs = sorted(glob.glob(os.path.join(OUT, f"{kind}_{tag(L, rho, r)}_seed_*.tsv")))
    fs = [f for f in fs if os.path.getsize(f) > 20]
    if not fs:
        raise OSError(f"no usable files for {kind} {tag(L, rho, r)}")
    return fs

def have_data(kind, L, rho, r):
    # a run still in flight has opened its files but written nothing
    return any(os.path.getsize(f) > 20
               for f in glob.glob(os.path.join(OUT, f"{kind}_{tag(L, rho, r)}_seed_*.tsv")))

def load_hist(kind, L, rho, r):
    # sum histograms over seeds
    total = {}
    for f in files(kind, L, rho, r):
        d = np.loadtxt(f, dtype=np.int64, ndmin=2)
        for s, c in d:
            total[s] = total.get(s, 0) + c
    sizes = np.array(sorted(total))
    return sizes.astype(float), np.array([total[s] for s in sizes], dtype=float)

def load_q(L, rho, r):
    # measured same-sign collision fraction, pooled over seeds
    coag = annih = 0
    for f in files("stats", L, rho, r):
        d = np.loadtxt(f, ndmin=2)
        coag += d[0, 0]; annih += d[0, 1]
    return coag / (coag + annih)

def load_snaps(L, rho, r):
    snaps = []
    for path in files("snapshots", L, rho, r):
        rows = []
        with open(path) as fh:
            for line in fh:
                fl = line.split()
                if len(fl) == L:
                    rows.append(np.fromiter(map(int, fl), dtype=np.int8, count=L))
        n = len(rows) // L
        if n:
            snaps.extend(np.array(rows[:n * L]).reshape(n, L, L))
    if not snaps:
        raise OSError(f"no usable snapshots for {tag(L, rho, r)}")
    return snaps

def logbin(sizes, counts, nb=26):
    m = sizes > 0
    sizes, counts = sizes[m], counts[m]
    edges = np.geomspace(sizes.min(), sizes.max() + 1, nb)
    h, _ = np.histogram(sizes, bins=edges, weights=counts)
    w = np.diff(edges)
    x = np.sqrt(edges[:-1] * edges[1:])
    ok = h > 0
    return x[ok], (h / w / counts.sum())[ok]

def mle(sizes, counts, xmin=XMIN):
    # Hill estimator for a discrete power law above xmin
    m = sizes >= xmin
    s, c = sizes[m], counts[m]
    if c.sum() < 100:
        return np.nan
    return 1 + c.sum() / np.sum(c * np.log(s / (xmin - 0.5)))

def local_slope(sizes, counts, nb=30):
    x, y = logbin(sizes, counts, nb)
    lx, ly = np.log(x), np.log(y)
    return x[1:-1], -(ly[2:] - ly[:-2]) / (lx[2:] - lx[:-2])

def plateau(sizes, counts, lo=30, hi_frac=50, nb=30):
    # median local log-slope over the clean interior: the honest test of whether the
    # distribution is a power law at all, since a curved log-log plot has no plateau
    x, sl = local_slope(sizes, counts, nb)
    m = (x >= lo) & (x <= sizes.max() / hi_frac)
    if m.sum() < 3:
        return np.nan, np.nan
    return np.median(sl[m]), sl[m].std()

def tau_closed(q):
    # cos(pi tau) = (1-q)/q on the branch tau in [3/2, 2]; no solution for q < 1/2
    q = np.asarray(q, dtype=float)
    ratio = (1 - q) / q
    return np.where(ratio <= 1, 2 - np.arccos(np.clip(ratio, -1, 1)) / np.pi, np.nan)


def plot_histograms(L=128, rho=0.2, rs=(0.0, 0.25, 0.5, 0.75, 1.0)):
    rs = [r for r in rs if have_data("spotSize", L, rho, r)]
    colors = plt.cm.plasma(np.linspace(0, 0.85, len(rs)))
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(14, 6))
    for r, col in zip(rs, colors):
        q = load_q(L, rho, r)
        for a, kind in ((a1, "spotSize"), (a2, "emission")):
            sizes, counts = load_hist(kind, L, rho, r)
            t, _ = plateau(sizes, counts)
            x, y = logbin(sizes, counts)
            a.plot(x, y, "o-", ms=4, lw=1, color=col,
                   label=rf"$r={g(r)}$, $q={q:.3f}$, $\tau={t:.2f}$")
    for a, xl, yl in ((a1, "spot size $m$", "$n(m)$"), (a2, "emission size $s$", "$P(s)$")):
        a.set_xscale("log"); a.set_yscale("log"); a.set_xlabel(xl); a.set_ylabel(yl)
        a.legend(frameon=False, fontsize=13); grid(a)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTS, "histograms",
                             f"histograms_L_{g(L)}_rho_{g(rho)}.png"), dpi=300)
    plt.close(fig)


def plot_local_slopes(L=128, rho=0.2, rs=(0.0, 0.25, 0.5, 0.75, 1.0)):
    # A power law is a horizontal line here. Curvature means a cutoff-limited fit
    # masquerading as an exponent, which is what mean field predicts for q < 1/2.
    rs = [r for r in rs if have_data("spotSize", L, rho, r)]
    colors = plt.cm.plasma(np.linspace(0, 0.85, len(rs)))
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(14, 6))
    for r, col in zip(rs, colors):
        for a, kind in ((a1, "spotSize"), (a2, "emission")):
            sizes, counts = load_hist(kind, L, rho, r)
            x, sl = local_slope(sizes, counts)
            a.plot(x, sl, "o-", ms=5, lw=1.6, color=col, label=rf"$r={g(r)}$")
    for a, xl, yl in ((a1, "spot size $m$", r"local slope $-d\log n/d\log m$"),
                      (a2, "emission size $s$", r"local slope $-d\log P/d\log s$")):
        a.set_xscale("log"); a.set_xlabel(xl); a.set_ylabel(yl)
        a.axhline(1.5, color="grey", ls="--", lw=1.5)
        a.axhline(2.0, color="grey", ls=":", lw=1.5)
        a.set_ylim(1.0, 4.0)
        a.legend(frameon=False, fontsize=13); grid(a)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTS, "histograms",
                             f"localSlopes_L_{g(L)}_rho_{g(rho)}.png"), dpi=300)
    plt.close(fig)


def plot_exponents_vs_r(L=128, rho=0.2,
                        rs=(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)):
    xs, tm, ts, qs, em = [], [], [], [], []
    for r in rs:
        try:
            a, ea = plateau(*load_hist("spotSize", L, rho, r))
            b, _ = plateau(*load_hist("emission", L, rho, r))
            q = load_q(L, rho, r)
        except OSError:
            continue
        xs.append(r); tm.append(a); ts.append(b); qs.append(q); em.append(ea)

    fig, ax = plt.subplots(figsize=(9, 6.5))
    ax.axhspan(1.5, 2.0, color="#ef8a62", alpha=0.15, lw=0)
    ax.errorbar(xs, tm, yerr=em, fmt="^-", color="#ef8a62", ms=8, lw=2.5,
                capsize=4, label=r"spot $\tau_m$")
    ax.plot(xs, ts, "o-", color="#b2182b", ms=8, lw=2.5, label=r"emission $\tau_s$")
    ax.plot(xs, 2 * np.array(tm) - 1, "k--", lw=2, label=r"$2\tau_m-1$")
    ax.set_xlabel("bipole probability $r$")
    ax.set_ylabel(r"power-law exponent $\tau$", color="#b2182b")
    ax.tick_params(axis="y", colors="#b2182b")
    ax.legend(frameon=False, loc="upper left")
    grid(ax)

    a2 = ax.twinx()
    a2.plot(xs, qs, "s--", color="#2166ac", ms=8, lw=2.5)
    a2.axhline(0.5, color="#2166ac", ls=":", lw=2)
    a2.text(0.02, 0.502, r"$q=1/2$: mean-field limit", color="#2166ac", fontsize=15)
    a2.set_ylabel("same-sign fraction $q$", color="#2166ac")
    a2.tick_params(axis="y", colors="#2166ac")
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTS, "exponents",
                             f"exponentsVsR_L_{g(L)}_rho_{g(rho)}.png"), dpi=300)
    plt.close(fig)


def plot_tau_vs_q(L=128, rhos=(0.2, 0.4, 0.6, 0.8),
                  rs=(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)):
    # The point of the model: it lives entirely at q <= 1/2, where the closed form
    # cos(pi tau) = (1-q)/q has no root at all. Anything measured here is 2D physics
    # the rate equation cannot see.
    colors = plt.cm.viridis(np.linspace(0, 0.85, len(rhos)))
    fig, ax = plt.subplots(figsize=(9, 6.5))
    qc = np.linspace(0.5, 1.0, 200)
    ax.plot(qc, tau_closed(qc), "-", color="#1a9850", lw=3, alpha=0.85,
            label=r"$\cos\pi\tau=(1-q)/q$")
    ax.axvline(0.5, color="grey", ls="--", lw=2)
    ax.axvspan(0.30, 0.5, color="grey", alpha=0.12, lw=0)
    for rho, col in zip(rhos, colors):
        q, t = [], []
        for r in rs:
            try:
                a, _ = plateau(*load_hist("spotSize", L, rho, r))
                q.append(load_q(L, rho, r)); t.append(a)
            except OSError:
                continue
        if q:
            ax.plot(q, t, "o", ms=9, color=col, label=rf"$\rho={g(rho)}$")
    ax.text(0.40, 2.45, "no mean-field\nsolution", color="grey", ha="center", fontsize=16)
    ax.set_xlabel("same-sign fraction $q$")
    ax.set_ylabel(r"spot exponent $\tau_m$")
    ax.set_xlim(0.32, 1.0); ax.set_ylim(1.4, 2.7)
    ax.legend(frameon=False, loc="lower left")
    grid(ax)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTS, "exponents", f"tauVsQ_L_{g(L)}.png"), dpi=300)
    plt.close(fig)


def plot_finite_size(rho=0.2, rs=(0.0, 0.5, 1.0), Ls=(64, 128, 256, 512)):
    # A cutoff-limited fit drifts with L; a genuine power law does not.
    fig, axes = plt.subplots(1, len(rs), figsize=(6 * len(rs), 5.6), squeeze=False)
    for a, r in zip(axes[0], rs):
        have = [L for L in Ls if have_data("spotSize", L, rho, r)]
        colors = plt.cm.viridis(np.linspace(0, 0.85, len(have)))
        for L, col in zip(have, colors):
            sizes, counts = load_hist("spotSize", L, rho, r)
            t, _ = plateau(sizes, counts)
            x, y = logbin(sizes, counts)
            a.plot(x, y, "o-", ms=4, lw=1, color=col, label=rf"$L={L}$, $\tau={t:.2f}$")
        a.set_xscale("log"); a.set_yscale("log")
        a.set_xlabel("spot size $m$"); a.set_ylabel("$n(m)$")
        a.set_title(rf"$r={g(r)}$", fontsize=20)
        a.legend(frameon=False, fontsize=13); grid(a)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTS, "finiteSize", f"finiteSize_rho_{g(rho)}.png"), dpi=300)
    plt.close(fig)


def plot_snapshots(L=128, rho=0.2, rs=(0.0, 0.5, 1.0)):
    cmap = ListedColormap(["#2166ac", "#f7f7f7", "#b2182b"])
    fig, axes = plt.subplots(1, len(rs), figsize=(5.2 * len(rs), 5.6), squeeze=False)
    for a, r in zip(axes[0], rs):
        try:
            snaps = load_snaps(L, rho, r)
        except OSError:
            continue
        a.imshow(snaps[-1], cmap=cmap, vmin=-1, vmax=1, interpolation="nearest")
        a.set_title(rf"$r={g(r)}$, $q={load_q(L, rho, r):.3f}$", fontsize=20)
        a.set_xticks([]); a.set_yticks([])
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTS, "snapshots",
                             f"snapshots_L_{g(L)}_rho_{g(rho)}.png"), dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    for sub in ("histograms", "exponents", "finiteSize", "snapshots"):
        os.makedirs(os.path.join(PLOTS, sub), exist_ok=True)
    plot_histograms()
    plot_local_slopes()
    plot_exponents_vs_r()
    plot_tau_vs_q()
    plot_finite_size()
    plot_snapshots()
