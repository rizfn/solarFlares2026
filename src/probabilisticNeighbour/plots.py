import os
import glob
import gzip
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

def opentext(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def files(kind, L, rho, p):
    fs = sorted(glob.glob(os.path.join(OUT, f"{kind}_{tag(L, rho, p)}_seed_*.tsv"))
                + glob.glob(os.path.join(OUT, f"{kind}_{tag(L, rho, p)}_seed_*.tsv.gz")))
    if not fs:
        raise OSError(f"no files for {kind} {tag(L, rho, p)}")
    return fs

def load_hist(kind, L, rho, p):
    # sum histograms over seeds; a run killed before it wrote leaves an empty file
    total = {}
    for f in files(kind, L, rho, p):
        if os.path.getsize(f) < 20:
            continue
        d = np.loadtxt(f, dtype=np.int64, ndmin=2)    # numpy reads .gz transparently
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

def load_snaps(L, rho, p):
    # all snapshots across seeds; blank separators are skipped, so reshape into LxL.
    # a run killed mid-write leaves a short final row, so keep only complete rows.
    snaps = []
    for path in files("snapshots", L, rho, p):
        rows = []
        with opentext(path) as fh:
            for line in fh:
                f = line.split()
                if len(f) == L:
                    rows.append(np.fromiter(map(int, f), dtype=np.int8, count=L))
        n = len(rows) // L
        if n:
            snaps.extend(np.array(rows[:n * L]).reshape(n, L, L))
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


def segregation(snaps, b=None):
    # Order parameter: mean |magnetization| of a b x b block, normalized by its
    # occupancy. A mixed surface still gives a nonzero value from sqrt(n) noise, so
    # subtract that floor in quadrature by re-measuring on sign-shuffled snapshots.
    L = snaps[0].shape[0]
    b = b or L // 4
    rng = np.random.default_rng(0)
    raw, floor = [], []
    for a in snaps:
        s = np.sign(a).astype(float)
        t = s.ravel().copy()
        occ = t != 0
        v = t[occ]; rng.shuffle(v); t[occ] = v
        for field, acc in ((s, raw), (t.reshape(L, L), floor)):
            blk = field.reshape(L // b, b, L // b, b)
            num = np.abs(blk.sum(axis=(1, 3)))
            den = np.abs(blk).sum(axis=(1, 3))
            acc.append((num[den > 0] / den[den > 0]).mean())
    m, m0 = np.mean(raw), np.mean(floor)
    return np.sqrt(max(m * m - m0 * m0, 0.0))

def voronoi_edges(sign, seed=0):
    # Periodic Voronoi (= Delaunay) neighbour pairs of the occupied sites of one
    # snapshot. Points are tiled 3x3 so that neighbours wrap; an edge of the periodic
    # graph then appears exactly twice among edges touching the central copy, which is a
    # uniform factor and drops out of any ratio. Site centres are jittered because the
    # Delaunay of a perfect lattice is degenerate (cocircular quadruples).
    from scipy.spatial import Delaunay
    L = sign.shape[0]
    y, x = np.nonzero(sign)
    n = len(x)
    rng = np.random.default_rng(seed)
    pts = np.column_stack([x, y]).astype(float) + 0.5 + rng.uniform(-1e-6, 1e-6, (n, 2))
    shifts = np.array([[i * L, j * L] for i in (-1, 0, 1) for j in (-1, 0, 1)])
    tiled = (pts[None, :, :] + shifts[:, None, :]).reshape(-1, 2)
    base = np.tile(np.arange(n), len(shifts))
    central = np.repeat(np.arange(len(shifts)), n) == 4   # shift (0,0)

    tri = Delaunay(tiled)
    s = tri.simplices
    e = np.vstack([s[:, [0, 1]], s[:, [1, 2]], s[:, [2, 0]]])
    keep = central[e[:, 0]] | central[e[:, 1]]
    e = e[keep]
    return sign[y, x], base[e[:, 0]], base[e[:, 1]]


def voronoi_same_sign(snaps, nmax=4):
    # Fraction of Voronoi-neighbour pairs sharing a sign. Unlike the lattice-neighbour
    # version this is defined for a dilute surface: every spot has neighbours regardless
    # of how many empty sites separate them.
    same = tot = 0
    for k, a in enumerate(snaps[-nmax:]):
        sgn, i, j = voronoi_edges(np.sign(a).astype(np.int8), seed=k)
        m = i != j                                        # a spot paired with its own image
        pr = sgn[i[m]].astype(int) * sgn[j[m]].astype(int)
        same += (pr > 0).sum(); tot += len(pr)
    return same / tot


def same_sign_frac(snaps):
    # q: fraction of occupied nearest-neighbour pairs that share a sign. A collision is
    # always between neighbours, so this is the branching ratio of the kinetic theory:
    # 1/2 when mixed, -> 1 when segregated.
    same = opp = 0
    for a in snaps:
        s = np.sign(a).astype(np.int8)
        for b in (np.roll(s, 1, 0), np.roll(s, 1, 1)):
            pr = s * b
            same += (pr > 0).sum(); opp += (pr < 0).sum()
    return same / (same + opp)


def available_ps(L, rho):
    # every p that has snapshot data at this (L, rho)
    return sorted({float(f.split("_p_")[1].split("_seed")[0])
                   for f in glob.glob(os.path.join(OUT, f"snapshots_L_{g(L)}_rho_{g(rho)}_p_*.tsv"))
                   + glob.glob(os.path.join(OUT, f"snapshots_L_{g(L)}_rho_{g(rho)}_p_*.tsv.gz"))})


def sign_correlation(snaps, nmax=8):
    # C(r): sign correlation of occupied pairs separated by r, conditioned on both sites
    # being occupied (the empty sites must not dilute it). Computed by FFT, radially
    # averaged. This is the honest correlation function -- q_V is only its value at the
    # nearest-neighbour distance, and converting that single number into a length
    # requires an assumed form, which is what made the earlier xi estimate unreliable.
    L = snaps[0].shape[0]
    num = np.zeros((L, L)); den = np.zeros((L, L))
    for a in snaps[-nmax:]:
        s = np.sign(a).astype(float)
        o = (a != 0).astype(float)
        num += np.fft.irfft2(np.abs(np.fft.rfft2(s)) ** 2, s=(L, L))
        den += np.fft.irfft2(np.abs(np.fft.rfft2(o)) ** 2, s=(L, L))
    num = np.fft.fftshift(num); den = np.fft.fftshift(den)
    c = L // 2
    yy, xx = np.indices((L, L))
    r = np.round(np.sqrt((yy - c) ** 2 + (xx - c) ** 2)).astype(int)
    C = np.bincount(r.ravel(), num.ravel()) / np.bincount(r.ravel(), den.ravel())
    return np.arange(len(C))[:c], C[:c]


def plot_correlation_function(rho=0.2, Ls=(128, 256, 512), ps=(0.6, 0.8, 0.9, 1.0)):
    # C(r) itself, one panel per p, all L overlaid. Read it as follows: if C(r) collapses
    # onto one curve as L grows, the length is intrinsic and finite. If instead the curve
    # keeps stretching with L, the length is set by the box -- which is what long-range
    # order looks like. Nothing is fitted here; the L-dependence is the whole message.
    # y is linear, not log: neutrality (sum_r C(r) = 0) forces C through zero at r ~ L/2,
    # and on log axes that sign change looks like a sharp cutoff, which it is not.
    # Top row vs r: collapse means an intrinsic, L-independent length.
    # Bottom row vs r/L: collapse means the only length is the box.
    fig, axes = plt.subplots(2, len(ps), figsize=(4.4 * len(ps), 8.5), sharey=True)
    colors = plt.cm.viridis(np.linspace(0, 0.75, len(Ls)))
    for j, p in enumerate(ps):
        for L, col in zip(Ls, colors):
            try:
                r, C = sign_correlation(load_snaps(L, rho, p))
            except OSError:
                continue
            m = r > 0
            axes[0][j].plot(r[m], C[m], "-", lw=2, color=col, label=rf"$L={L}$")
            axes[1][j].plot(r[m] / L, C[m], "-", lw=2, color=col, label=rf"$L={L}$")
        for i, xl in ((0, "$r$"), (1, "$r/L$")):
            a = axes[i][j]
            a.set_xscale("log"); a.axhline(0, color="k", lw=1)
            a.set_xlabel(xl); grid(a)
        axes[0][j].set_title(f"$p={g(p)}$")
        axes[0][j].legend(frameon=False, fontsize=13)
    for i in (0, 1):
        axes[i][0].set_ylabel("$C(r)$")
    fig.tight_layout()
    os.makedirs(os.path.join(PLOTS, "correlationLength"), exist_ok=True)
    fig.savefig(os.path.join(PLOTS, "correlationLength",
                             f"correlationFunction_rho_{g(rho)}.png"), dpi=300)
    plt.close(fig)


def plot_qv_finite_size(rho=0.2, Ls=(64, 128, 256, 512)):
    # The crossing test. A genuine critical point shows up as curves for different L
    # separating, and crossing at p_c. If the curves lie on top of each other all the way
    # to p=1, the only singular point is p=1 and the rest is a crossover.
    fig, ax = plt.subplots(figsize=(8.5, 6))
    colors = plt.cm.viridis(np.linspace(0, 0.75, len(Ls)))
    for L, col in zip(Ls, colors):
        xs, ys, es = [], [], []
        for p in available_ps(L, rho):
            try:
                snaps = load_snaps(L, rho, p)
            except OSError:
                continue
            vals = [voronoi_same_sign([s]) for s in snaps[-6:]]
            xs.append(p); ys.append(np.mean(vals))
            es.append(np.std(vals) / np.sqrt(len(vals)))
        if xs:
            ax.errorbar(xs, ys, yerr=es, fmt="o-", ms=6, lw=2, capsize=3,
                        color=col, label=rf"$L={L}$")
    ax.axhline(0.5, color="k", ls="--", lw=1)
    ax.set_xlabel("neighbour probability $p$")
    ax.set_ylabel(r"same-sign fraction $q_V$")
    ax.legend(frameon=False)
    grid(ax)
    fig.tight_layout()
    os.makedirs(os.path.join(PLOTS, "correlationLength"), exist_ok=True)
    fig.savefig(os.path.join(PLOTS, "correlationLength",
                             f"qvFiniteSize_rho_{g(rho)}.png"), dpi=300)
    plt.close(fig)


def load_max(L, rho, p):
    # largest-spot trajectories, one array per seed: step, maxAbs, sumAbs, nSpots
    out = []
    for f in files("maxSpot", L, rho, p):
        if os.path.getsize(f) < 40:
            continue
        d = np.loadtxt(f, ndmin=2)
        if len(d) > 1:
            out.append(d)
    if not out:
        raise OSError(f"no maxSpot data for {tag(L, rho, p)}")
    n = min(len(d) for d in out)
    return np.stack([d[:n] for d in out])


def plot_max_trajectories(L=128, rho=0.2, ps=(0.0, 0.2, 0.4, 0.6, 0.8, 1.0)):
    # Takayasu-style condensation would show the largest spot growing without bound;
    # annihilation may instead regulate it. Left: absolute mass of the largest spot.
    # Right: its share of the total mass, which is the condensation order parameter.
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(14, 6))
    colors = plt.cm.viridis(np.linspace(0, 0.85, len(ps)))
    for p, col in zip(ps, colors):
        try:
            d = load_max(L, rho, p)
        except OSError:
            continue
        t = d[0, :, 0]
        ok = t > 0
        # single seeds, not the seed average: the maximum grows and crashes, and
        # averaging over seeds smooths exactly the feature that answers the question
        a1.plot(t[ok], d[0, ok, 1], "-", lw=1.6, color=col, label=rf"$p={g(p)}$")
        share = (d[:, :, 1] / d[:, :, 2]).mean(axis=0)
        a2.plot(t[ok], share[ok], "-", lw=2, color=col, label=rf"$p={g(p)}$")
    a1.set_xscale("log"); a1.set_yscale("log")
    a1.set_xlabel("time $t$ (sweeps)"); a1.set_ylabel(r"largest spot $m_{\max}$")
    a2.set_xscale("log")
    a2.set_xlabel("time $t$ (sweeps)")
    a2.set_ylabel(r"mass share $m_{\max}/\sum m$")
    for a in (a1, a2):
        a.legend(frameon=False, fontsize=13, ncol=2); grid(a)
    fig.tight_layout()
    os.makedirs(os.path.join(PLOTS, "largestSpot"), exist_ok=True)
    fig.savefig(os.path.join(PLOTS, "largestSpot",
                             f"maxTrajectories_L_{g(L)}_rho_{g(rho)}.png"), dpi=300)
    plt.close(fig)


def plot_max_finite_size(rho=0.2, Ls=(64, 128, 256, 512), ps=(0.0, 0.5, 1.0)):
    # The decisive test. Steady-state maximum vs L: if the largest spot is regulated by
    # annihilation it is L-independent; if it condenses it grows with the system, and
    # m_max ~ L^2 (a finite share of all the mass) is the Takayasu-like case.
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(14, 6))
    colors = plt.cm.viridis(np.linspace(0, 0.7, len(ps)))
    for p, col in zip(ps, colors):
        xs, ys, ss = [], [], []
        for L in Ls:
            try:
                d = load_max(L, rho, p)
            except OSError:
                continue
            half = d.shape[1] // 2          # steady state: second half of the run
            xs.append(L)
            ys.append(d[:, half:, 1].mean())
            ss.append((d[:, half:, 1] / d[:, half:, 2]).mean())
        if not xs:
            continue
        a1.plot(xs, ys, "o-", ms=9, lw=2.5, color=col, label=rf"$p={g(p)}$")
        a2.plot(xs, ss, "o-", ms=9, lw=2.5, color=col, label=rf"$p={g(p)}$")
    xg = np.array(Ls, dtype=float)
    a1.plot(xg, xg ** 2 * 0.5, "k--", lw=2, label=r"$L^2$")
    a1.set_xscale("log"); a1.set_yscale("log")
    a1.set_xlabel("system size $L$"); a1.set_ylabel(r"steady-state $\langle m_{\max}\rangle$")
    a2.set_xscale("log")
    a2.set_xlabel("system size $L$"); a2.set_ylabel(r"mass share $m_{\max}/\sum m$")
    a2.set_ylim(0, 1)
    for a in (a1, a2):
        a.legend(frameon=False, fontsize=14); grid(a)
    fig.tight_layout()
    os.makedirs(os.path.join(PLOTS, "largestSpot"), exist_ok=True)
    fig.savefig(os.path.join(PLOTS, "largestSpot", f"maxVsL_rho_{g(rho)}.png"), dpi=300)
    plt.close(fig)


def plot_snapshots(L=128, rhos=(0.2, 0.6), ps=(0.0, 0.5, 1.0)):
    cmap = ListedColormap(["#2166ac", "#f7f7f7", "#b2182b"])
    fig, axes = plt.subplots(len(rhos), len(ps), figsize=(3.3 * len(ps), 3.3 * len(rhos)))
    for i, rho in enumerate(rhos):
        for j, p in enumerate(ps):
            ax = axes[i][j]
            ax.imshow(load_snaps(L, rho, p)[-1], cmap=cmap, vmin=-1, vmax=1, interpolation="nearest")
            ax.set_xticks([]); ax.set_yticks([])
            if i == 0:
                ax.set_title(f"$p={g(p)}$")
            if j == 0:
                ax.set_ylabel(rf"$\rho={g(rho)}$")
    fig.tight_layout()
    name = f"snapshots_L_{g(L)}_rho_{'_'.join(g(r) for r in rhos)}_p_{'_'.join(g(p) for p in ps)}.png"
    fig.savefig(os.path.join(PLOTS, "snapshots", name), dpi=300)
    plt.close(fig)


def plot_voronoi_same_sign(L=128, rhos=(0.2, 0.6)):
    fig, ax = plt.subplots(figsize=(8, 6))
    colors = plt.cm.viridis(np.linspace(0, 0.6, len(rhos)))
    for rho, col in zip(rhos, colors):
        ps = available_ps(L, rho)
        xs, ys, es = [], [], []
        for p in ps:
            vals = [voronoi_same_sign([s]) for s in load_snaps(L, rho, p)[-6:]]
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
                             f"voronoiSameSign_L_{g(L)}.png"), dpi=300)
    plt.close(fig)


def plot_snapshots_voronoi(L=128, rhos=(0.2, 0.6), ps=(0.0, 0.5, 1.0)):
    from scipy.spatial import Voronoi
    from matplotlib.collections import PolyCollection
    fig, axes = plt.subplots(len(rhos), len(ps), figsize=(3.3 * len(ps), 3.3 * len(rhos)))
    for i, rho in enumerate(rhos):
        for j, p in enumerate(ps):
            ax = np.atleast_2d(axes)[i][j]
            sign = np.sign(load_snaps(L, rho, p)[-1]).astype(np.int8)
            y, x = np.nonzero(sign)
            n = len(x)
            rng = np.random.default_rng(0)
            pts = np.column_stack([x, y]).astype(float) + 0.5 + rng.uniform(-1e-6, 1e-6, (n, 2))
            shifts = np.array([[a * L, b * L] for a in (-1, 0, 1) for b in (-1, 0, 1)])
            vor = Voronoi((pts[None] + shifts[:, None]).reshape(-1, 2))
            # draw every tile and let the axes clip: cells of image points still cover
            # part of the box, so drawing only the central copy leaves gaps at the edges
            polys, cols = [], []
            s = np.tile(sign[y, x], len(shifts))
            for k in range(len(shifts) * n):
                reg = vor.regions[vor.point_region[k]]
                if reg and -1 not in reg:
                    polys.append(vor.vertices[reg])
                    cols.append("#b2182b" if s[k] > 0 else "#2166ac")
            ax.add_collection(PolyCollection(polys, facecolors=cols, linewidths=0.15,
                                             edgecolors="white"))
            ax.set_xlim(0, L); ax.set_ylim(L, 0); ax.set_aspect("equal")
            ax.set_xticks([]); ax.set_yticks([])
            if i == 0:
                ax.set_title(f"$p={g(p)}$")
            if j == 0:
                ax.set_ylabel(rf"$\rho={g(rho)}$")
    fig.tight_layout()
    os.makedirs(os.path.join(PLOTS, "snapshotsVoronoi"), exist_ok=True)
    name = f"snapshotsVoronoi_L_{g(L)}_rho_{'_'.join(g(r) for r in rhos)}_p_{'_'.join(g(p) for p in ps)}.png"
    fig.savefig(os.path.join(PLOTS, "snapshotsVoronoi", name), dpi=300)
    plt.close(fig)


def plot_histograms(L=128, rhos=(0.2, 0.4, 0.6, 0.8), ps=(0.0, 0.5, 1.0), offset=12):
    colors = plt.cm.viridis(np.linspace(0, 0.85, len(rhos)))
    for p in ps:
        fig, (a1, a2) = plt.subplots(1, 2, figsize=(14, 6))
        for a, kind in ((a1, "spotSize"), (a2, "emission")):
            taus, amps, spans = [], [], []
            for rho, col in zip(rhos, colors):
                sizes, counts = load_hist(kind, L, rho, p)
                tau = window_tau(sizes, counts)
                x, y = logbin(sizes, counts)
                a.plot(x, y, "o", ms=5, color=col, label=rf"$\rho={g(rho)}$, $\tau={tau:.2f}$")
                f = (x >= WIN_LO) & (x <= sizes.max() / 20)
                if tau > 0 and f.sum() > 2:
                    taus.append(tau); amps.append(np.median(y[f] * x[f] ** tau))
                    spans.append((x[f].min(), x[f].max()))
            # one guideline for the whole panel, lifted clear of the data it describes
            if taus:
                tau = np.mean(taus)
                xg = np.geomspace(min(s[0] for s in spans), max(s[1] for s in spans), 50)
                a.plot(xg, offset * np.mean(amps) * xg ** -tau, "k--", lw=2,
                       label=rf"$\tau={tau:.2f}$")
        for a, xl, yl in [(a1, "spot size $m$", "$n(m)$"), (a2, "emission size $s$", "$P(s)$")]:
            a.set_xscale("log"); a.set_yscale("log"); a.set_xlabel(xl); a.set_ylabel(yl)
            a.legend(frameon=False); grid(a)
        fig.suptitle(f"$p={g(p)}$", fontsize=20)
        fig.tight_layout()
        fig.savefig(os.path.join(PLOTS, "histograms", f"histograms_L_{g(L)}_p_{g(p)}.png"), dpi=300)
        plt.close(fig)


def plot_exponents_vs_p(L=128, rho=0.2):
    # both cascade exponents and the Voronoi same-sign fraction share one x axis:
    # tau slides with p, while q_V says how ordered the surface is at that p
    xs, te, ts, qs = [], [], [], []
    for p in available_ps(L, rho):
        try:
            e = window_tau(*load_hist("emission", L, rho, p))
            m = window_tau(*load_hist("spotSize", L, rho, p))
            q = voronoi_same_sign(load_snaps(L, rho, p))
        except OSError:
            continue
        xs.append(p); te.append(e); ts.append(m); qs.append(q)

    fig, ax = plt.subplots(figsize=(9, 6.5))
    # mean-field windows: tau_m in [3/2, 2] as q runs 1 -> 1/2, and tau_s = 2 tau_m - 1
    ax.axhspan(1.5, 2.0, color="#ef8a62", alpha=0.15, lw=0)
    ax.axhspan(2.0, 3.0, color="#b2182b", alpha=0.10, lw=0)
    ax.plot(xs, te, "o-", color="#b2182b", ms=8, lw=2.5, label=r"emission $\tau_s$")
    ax.plot(xs, ts, "^-", color="#ef8a62", ms=8, lw=2.5, label=r"spot $\tau_m$")
    # mean field predicts the emission cascade from the spot one: a spot of mass m emits
    # s ~ m, so tau_s = 2 tau_m - 1. The gap to the measured tau_s is where that fails.
    ax.plot(xs, 2 * np.array(ts) - 1, "--", color="k", lw=2,
            label=r"$2\tau_m - 1$ (predicted $\tau_s$)")
    ax.set_xlabel("neighbour probability $p$")
    ax.set_ylabel(r"power-law exponent $\tau$", color="#b2182b")
    ax.tick_params(axis="y", colors="#b2182b")
    ax.set_ylim(1.3, 3.5)
    ax.legend(frameon=False, loc="center left", bbox_to_anchor=(0.02, 0.42))
    grid(ax)

    a2 = ax.twinx()
    a2.plot(xs, qs, "s--", color="#2166ac", ms=8, lw=2.5)
    a2.axhline(0.5, color="#2166ac", ls=":", lw=1.5, alpha=0.7)
    a2.set_ylabel(r"same-sign fraction $q_V$", color="#2166ac")
    a2.tick_params(axis="y", colors="#2166ac")
    a2.set_ylim(0.45, 1.02)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTS, "exponents", f"exponents_L_{g(L)}_rho_{g(rho)}.png"), dpi=300)
    plt.close(fig)


def window_slope(sizes, counts, lo, hi, nb=30):
    # least-squares log-log slope over a window fixed in absolute flux units, so the
    # same range of physical scales is fitted at every L
    x, y = logbin(sizes, counts, nb)
    m = (x >= lo) & (x <= hi)
    if m.sum() < 3:
        return np.nan
    return -np.polyfit(np.log(x[m]), np.log(y[m]), 1)[0]


# fit windows for the p=1 finite-size study, chosen to sit inside the plateau at every L
SPOT_WIN, EMIS_WIN = (100, 3e4), (50, 3e3)


def plot_finite_size(rho=0.2, p=1.0, Ls=(32, 48, 64, 96, 128, 256, 512, 1024)):
    # At p=1 the spot exponent is L-independent at the Takayasu value, while the
    # emission exponent drifts down with L, away from 2 tau_m - 1.
    Ls = [L for L in Ls
          if all(any(os.path.getsize(f) > 20
                     for f in glob.glob(os.path.join(OUT, f"{k}_{tag(L, rho, p)}_seed_*.tsv")))
                 for k in ("spotSize", "emission"))]
    colors = plt.cm.viridis(np.linspace(0, 0.88, len(Ls)))

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(14, 6))
    taus = {"spotSize": [], "emission": []}
    for L, col in zip(Ls, colors):
        for a, kind, win in ((a1, "spotSize", SPOT_WIN), (a2, "emission", EMIS_WIN)):
            sizes, counts = load_hist(kind, L, rho, p)
            t = window_slope(sizes, counts, *win)
            taus[kind].append(t)
            x, y = logbin(sizes, counts)
            a.plot(x, y, "o-", ms=4, lw=1, color=col, label=rf"$L={L}$, $\tau={t:.2f}$")
    xg = np.geomspace(*SPOT_WIN, 20)
    a1.plot(xg, 0.4 * xg ** -1.5, "k--", lw=2.5, label=r"$m^{-3/2}$ (Takayasu)")
    xg = np.geomspace(*EMIS_WIN, 20)
    a2.plot(xg, 4.0 * xg ** -2.0, "k--", lw=2.5, label=r"$s^{-2}=s^{-(2\tau_m-1)}$")
    a2.plot(xg, 0.35 * xg ** -1.5, "k:", lw=2.5, label=r"$s^{-3/2}=s^{-\tau_m}$")
    for a, xl, yl in ((a1, "spot size $m$", "$n(m)$"), (a2, "emission size $s$", "$P(s)$")):
        a.set_xscale("log"); a.set_yscale("log"); a.set_xlabel(xl); a.set_ylabel(yl)
        a.legend(frameon=False, fontsize=13); grid(a)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTS, "finiteSize",
                             f"finiteSize_rho_{g(rho)}_p_{g(p)}.png"), dpi=300)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 6.5))
    iL = 1.0 / np.array(Ls, dtype=float)
    tm, ts = np.array(taus["spotSize"]), np.array(taus["emission"])
    ax.axhline(1.5, color="grey", ls="--", lw=2)
    ax.axhline(2.0, color="grey", ls=":", lw=2)
    ax.text(0.0005, 1.52, r"$3/2$", color="grey", fontsize=17)
    ax.text(0.0005, 2.02, r"$2\tau_m-1=2$", color="grey", fontsize=17)
    ax.plot(iL, tm, "^-", color="#ef8a62", ms=9, lw=2.5, label=r"spot $\tau_m$")
    ax.plot(iL, ts, "o-", color="#b2182b", ms=9, lw=2.5, label=r"emission $\tau_s$")
    for a, y in ((iL, tm), (iL, ts)):  # 1/L extrapolation from the large-L points
        f = a <= 1 / 64.
        if f.sum() >= 3:
            b, c = np.polyfit(a[f], y[f], 1)
            xx = np.linspace(0, a[f].max(), 10)
            ax.plot(xx, c + b * xx, "-", color="k", lw=1, alpha=0.5)
            ax.plot([0], [c], "k*", ms=14)
    ax.set_xlabel("$1/L$"); ax.set_ylabel(r"power-law exponent $\tau$")
    ax.set_xlim(-0.001, 0.022); ax.set_ylim(1.35, 2.45)
    ax.legend(frameon=False, loc="lower right"); grid(ax)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTS, "finiteSize",
                             f"exponentsVsInvL_rho_{g(rho)}_p_{g(p)}.png"), dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    plot_finite_size()
    plot_snapshots()
    plot_snapshots(ps=(0.4, 0.6, 0.8))
    plot_snapshots_voronoi()
    plot_snapshots_voronoi(ps=(0.4, 0.6, 0.8))
    plot_max_trajectories()
    plot_max_trajectories(rho=0.6)
    plot_max_finite_size()
    plot_voronoi_same_sign()
    plot_qv_finite_size()
    plot_correlation_function()
    plot_histograms()
    plot_exponents_vs_p()
