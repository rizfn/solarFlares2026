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

WIN_LO, WIN_FTOP = 60.0, 3e-4   # scaling-window fit: [lo, the mass above which a fraction ftop lies]


def g(v):
    return f"{v:g}"


def grid(ax):
    ax.grid(True, which="major", ls=":", lw=0.9, alpha=0.6)
    ax.set_axisbelow(True)


def tag(L, rho, p):
    return f"L_{g(L)}_rho_{g(rho)}_p_{g(p)}"


def files(kind, L, rho, p):
    fs = sorted(glob.glob(os.path.join(OUT, f"{kind}_{tag(L, rho, p)}_seed_*.tsv")))
    if not fs:
        raise OSError(f"no files for {kind} {tag(L, rho, p)}")
    return fs


def load_hist(kind, L, rho, p):
    total = {}
    for f in files(kind, L, rho, p):
        if os.path.getsize(f) < 20:
            continue
        d = np.loadtxt(f, dtype=np.int64, ndmin=2)
        for s, c in d:
            total[s] = total.get(s, 0) + c
    if not total:
        raise OSError(f"all {kind} files empty for {tag(L, rho, p)}")
    sizes = np.array(sorted(total))
    return sizes.astype(float), np.array([total[s] for s in sizes], dtype=float)


def logbin(sizes, counts, nb=26):
    m = sizes > 0
    sizes, counts = sizes[m], counts[m]
    edges = np.geomspace(sizes.min(), sizes.max() + 1, nb)
    h, _ = np.histogram(sizes, bins=edges, weights=counts)
    w = np.diff(edges)
    x = np.sqrt(edges[:-1] * edges[1:])
    ok = h > 0
    return x[ok], (h / w / counts.sum())[ok]


def load_snaps(L, rho, p, last=None):
    # each snapshot is L*L rows of L values, z-major -> reshape to (L, L, L) as [z][y][x]
    snaps = []
    for path in files("snapshots", L, rho, p):
        rows = []
        with open(path) as fh:
            for line in fh:
                f = line.split()
                if len(f) == L:
                    rows.append(np.fromiter(map(int, f), dtype=np.int8, count=L))
        n = len(rows) // (L * L)
        if n:
            a = np.array(rows[:n * L * L]).reshape(n, L, L, L)
            snaps.extend(a[-last:] if last else a)
    if not snaps:
        raise OSError(f"no usable snapshots for {tag(L, rho, p)}")
    return snaps


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


def available_ps(L, rho):
    return sorted({float(f.split("_p_")[1].split("_seed")[0])
                   for f in glob.glob(os.path.join(OUT, f"snapshots_L_{g(L)}_rho_{g(rho)}_p_*.tsv"))})


def voronoi_same_sign(snaps, seed=0):
    # Fraction of Voronoi-neighbour pairs sharing a sign, in 3D. Voronoi neighbours are
    # Delaunay edges, so the tetrahedra give them directly. Periodicity is handled by
    # replicating only the points within `margin` of a face, rather than all 26 image
    # copies -- a cell reaching further than a few mean spacings is rare enough to ignore,
    # and full replication is prohibitive in 3D. Centres get a tiny jitter because the
    # Delaunay of a perfect lattice is degenerate (cospherical points).
    from scipy.spatial import Delaunay
    same = tot = 0
    for k, a in enumerate(snaps):
        L = a.shape[0]
        sign = np.sign(a).astype(np.int8)
        z, y, x = np.nonzero(sign)
        n = len(x)
        rng = np.random.default_rng(seed + k)
        pts = np.column_stack([x, y, z]).astype(float) + 0.5
        pts += rng.uniform(-1e-6, 1e-6, (n, 3))
        margin = max(4.0 * (n / L ** 3) ** (-1 / 3.0), 4.0)

        allp = [pts]
        base = [np.arange(n)]
        for dx in (-L, 0, L):
            for dy in (-L, 0, L):
                for dz in (-L, 0, L):
                    if dx == dy == dz == 0:
                        continue
                    q = pts + np.array([dx, dy, dz])
                    keep = np.all((q > -margin) & (q < L + margin), axis=1)
                    if keep.any():
                        allp.append(q[keep]); base.append(np.arange(n)[keep])
        P = np.vstack(allp)
        B = np.concatenate(base)
        central = np.zeros(len(B), bool); central[:n] = True

        tri = Delaunay(P)
        s = tri.simplices
        e = np.vstack([s[:, [0, 1]], s[:, [0, 2]], s[:, [0, 3]],
                       s[:, [1, 2]], s[:, [1, 3]], s[:, [2, 3]]])
        e = e[central[e[:, 0]] | central[e[:, 1]]]
        i, j = B[e[:, 0]], B[e[:, 1]]
        m = i != j
        sgn = sign[z, y, x]
        pr = sgn[i[m]].astype(int) * sgn[j[m]].astype(int)
        same += (pr > 0).sum(); tot += len(pr)
    return same / tot


def plot_voronoi_same_sign(L=32, rhos=(0.2, 0.6), nsnap=3):
    fig, ax = plt.subplots(figsize=(8, 6))
    colors = plt.cm.viridis(np.linspace(0, 0.6, len(rhos)))
    for rho, col in zip(rhos, colors):
        xs, ys, es = [], [], []
        for p in available_ps(L, rho):
            snaps = load_snaps(L, rho, p, last=nsnap)
            vals = [voronoi_same_sign([s]) for s in snaps]
            xs.append(p); ys.append(np.mean(vals))
            es.append(np.std(vals) / np.sqrt(len(vals)))
        ax.errorbar(xs, ys, yerr=es, fmt="o-", color=col, ms=5, capsize=3,
                    label=rf"$\rho={g(rho)}$")
    ax.axhline(0.5, color="k", ls="--", lw=1)
    ax.set_xlabel("neighbour probability $p$")
    ax.set_ylabel(r"same-sign fraction $q_V$")
    ax.legend(frameon=False)
    grid(ax)
    fig.tight_layout()
    os.makedirs(os.path.join(PLOTS, "correlationLength"), exist_ok=True)
    fig.savefig(os.path.join(PLOTS, "correlationLength",
                             f"voronoiSameSign3D_L_{g(L)}.png"), dpi=300)
    plt.close(fig)


def plot_exponents_vs_p(L=32, rho=0.2, nsnap=3):
    xs, te, ts, qs = [], [], [], []
    for p in available_ps(L, rho):
        try:
            e = window_tau(*load_hist("emission", L, rho, p))
            m = window_tau(*load_hist("spotSize", L, rho, p))
            q = voronoi_same_sign(load_snaps(L, rho, p, last=nsnap))
        except OSError:
            continue
        xs.append(p); te.append(e); ts.append(m); qs.append(q)

    fig, ax = plt.subplots(figsize=(9, 6.5))
    ax.plot(xs, te, "o-", color="#b2182b", ms=8, lw=2.5, label=r"emission $\tau_s$")
    ax.plot(xs, ts, "^-", color="#ef8a62", ms=8, lw=2.5, label=r"spot $\tau_m$")
    # mean field ties the two together: a spot of mass m emits s ~ m, so tau_s = 2 tau_m - 1
    ax.plot(xs, 2 * np.array(ts) - 1, "--", color="k", lw=2,
            label=r"$2\tau_m - 1$ (predicted $\tau_s$)")
    ax.set_xlabel("neighbour probability $p$")
    ax.set_ylabel(r"power-law exponent $\tau$", color="#b2182b")
    ax.tick_params(axis="y", colors="#b2182b")
    ax.legend(frameon=False, loc="center left", bbox_to_anchor=(0.02, 0.35))
    grid(ax)

    a2 = ax.twinx()
    a2.plot(xs, qs, "s--", color="#2166ac", ms=8, lw=2.5)
    a2.axhline(0.5, color="#2166ac", ls=":", lw=1.5, alpha=0.7)
    a2.set_ylabel(r"same-sign fraction $q_V$", color="#2166ac")
    a2.tick_params(axis="y", colors="#2166ac")
    a2.set_ylim(0.45, 1.02)
    fig.tight_layout()
    os.makedirs(os.path.join(PLOTS, "exponents"), exist_ok=True)
    fig.savefig(os.path.join(PLOTS, "exponents",
                             f"exponents3D_L_{g(L)}_rho_{g(rho)}.png"), dpi=300)
    plt.close(fig)


def plot_histograms(L=32, rhos=(0.2, 0.6), ps=(0.0, 0.5, 1.0)):
    colors = plt.cm.viridis(np.linspace(0, 0.6, len(rhos)))
    for p in ps:
        fig, (a1, a2) = plt.subplots(1, 2, figsize=(14, 6))
        for a, kind in ((a1, "spotSize"), (a2, "emission")):
            for rho, col in zip(rhos, colors):
                try:
                    sizes, counts = load_hist(kind, L, rho, p)
                except OSError:
                    continue
                t = window_tau(sizes, counts)
                x, y = logbin(sizes, counts)
                a.plot(x, y, "o", ms=5, color=col,
                       label=rf"$\rho={g(rho)}$, $\tau={t:.2f}$")
        for a, xl, yl in [(a1, "spot size $m$", "$n(m)$"), (a2, "emission size $s$", "$P(s)$")]:
            a.set_xscale("log"); a.set_yscale("log"); a.set_xlabel(xl); a.set_ylabel(yl)
            a.legend(frameon=False); grid(a)
        fig.suptitle(f"3D, $p={g(p)}$", fontsize=20)
        fig.tight_layout()
        os.makedirs(os.path.join(PLOTS, "histograms"), exist_ok=True)
        fig.savefig(os.path.join(PLOTS, "histograms",
                                 f"histograms3D_L_{g(L)}_p_{g(p)}.png"), dpi=300)
        plt.close(fig)


def plot_slices(L=32, rhos=(0.2, 0.6), ps=(0.0, 0.5, 1.0)):
    # a 3D snapshot cannot be shown directly; take the z = L/2 plane through the cube
    cmap = ListedColormap(["#2166ac", "#f7f7f7", "#b2182b"])
    fig, axes = plt.subplots(len(rhos), len(ps), figsize=(3.3 * len(ps), 3.3 * len(rhos)))
    for i, rho in enumerate(rhos):
        for j, p in enumerate(ps):
            ax = np.atleast_2d(axes)[i][j]
            a = load_snaps(L, rho, p, last=1)[-1]
            ax.imshow(np.sign(a[L // 2]), cmap=cmap, vmin=-1, vmax=1,
                      interpolation="nearest")
            ax.set_xticks([]); ax.set_yticks([])
            if i == 0:
                ax.set_title(f"$p={g(p)}$")
            if j == 0:
                ax.set_ylabel(rf"$\rho={g(rho)}$")
    fig.tight_layout()
    os.makedirs(os.path.join(PLOTS, "slices"), exist_ok=True)
    name = f"slices3D_L_{g(L)}_rho_{'_'.join(g(r) for r in rhos)}_p_{'_'.join(g(p) for p in ps)}.png"
    fig.savefig(os.path.join(PLOTS, "slices", name), dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    plot_slices()
    plot_voronoi_same_sign()
    plot_exponents_vs_p()
    plot_exponents_vs_p(rho=0.6)
    plot_histograms()
