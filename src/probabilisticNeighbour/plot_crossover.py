# Where the two-slope emission spectrum comes from.
#
# The emission spectrum P(s) bends at intermediate p, with a steep branch at small s
# and a shallow one at large s. It is not two populations and not two fixed points.
# It is one power law times one modulating factor:
#
#     P_emis(s)  =  P_coll(s) * [1 - q(s)]
#
# P_coll(s) is the distribution of min(i,j) over ALL collisions, which is a clean
# s^(1-2 tau_m) ~ s^-2 -- the mean-field result, exactly obeyed. [1 - q(s)] is the
# fraction of those collisions that actually annihilate, and it is U-shaped: falling
# at small s, rising at large s. Its minimum is the crossover scale.
#
# The U comes from the geometry. Box-counting the set where + and - sit close enough
# to react gives dimension 1 below a length xi(p) and dimension 2 above it: domains
# are compact, so reactions are confined to their 1D boundaries at short range, while
# at long range many domains tile the plane and the reactive set fills space again.
# xi(p) is independent of L for every p < 1, and only at p = 1 does it grow with L.
import os
import numpy as np
import matplotlib.pyplot as plt

DIR = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(DIR, "outputs")
PLOTS = os.path.join(DIR, "plots", "crossover")
SERVER = os.path.join(DIR, "..", "probabilisticNeighbourServer", "outputs")

plt.rcParams.update({"font.size": 18, "axes.labelsize": 22, "xtick.labelsize": 16,
                     "ytick.labelsize": 16, "legend.fontsize": 14})


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


def binned(sizes, counts, edges):
    # raw counts per log bin, on edges supplied by the caller
    h, _ = np.histogram(sizes, edges, weights=counts)
    return h


def edges_for(*hists, perdec=6):
    top = max(s.max() for s, _ in hists)
    return np.geomspace(1, top + 1, int(perdec * np.log10(top)) + 2)


def spectra(L, rho, p, perdec=6, minc=400):
    """P_coll(s), P_emis(s) and their ratio, all on ONE set of bin edges.

    The ratio is only 1 - q(s) if both spectra are binned identically and left
    unnormalised: it is the fraction of collisions with min(i,j) = s that were
    annihilations rather than coagulations, so it must lie in [0, 1].
    """
    sc, cc = load_hist("chCollMin", L, rho, p)
    se, ce = load_hist("chEmisAll", L, rho, p)
    e = edges_for((sc, cc), (se, ce), perdec=perdec)
    hc, he = binned(sc, cc, e), binned(se, ce, e)
    x, w = np.sqrt(e[:-1] * e[1:]), np.diff(e)
    k = hc >= minc
    return x[k], hc[k] / w[k] / cc.sum(), he[k] / w[k] / cc.sum(), he[k] / hc[k]


def crossover_scale(x, r):
    # s* = the minimum of 1 - q(s), located on a smoothed curve so a single noisy
    # bin cannot claim it. Returns nan if the minimum sits at either edge, which
    # means the U is not resolved inside the measured range.
    if len(r) < 7:
        return np.nan
    lr = np.log(r)
    sm = np.convolve(lr, np.ones(3) / 3, mode="same")
    sm[0], sm[-1] = lr[0], lr[-1]
    j = int(np.argmin(sm))
    return np.nan if j in (0, len(r) - 1) else x[j]


# ---------------------------------------------------------------- interface geometry

def wall_points(sign, seed=0):
    """Midpoints of the Voronoi bonds joining opposite-sign spots: the domain wall
    as a point set. Defined by the neighbour graph, not by lattice occupancy, so it
    does not thin out at low density -- counting boxes that merely *contain* both
    signs instead measures sparsity at small eps and returns negative dimensions."""
    from scipy.spatial import Delaunay
    L = sign.shape[0]
    y, x = np.nonzero(sign)
    n = len(x)
    rng = np.random.default_rng(seed)
    pts = np.column_stack([x, y]).astype(float) + 0.5 + rng.uniform(-1e-6, 1e-6, (n, 2))
    shifts = np.array([[i * L, j * L] for i in (-1, 0, 1) for j in (-1, 0, 1)])
    tiled = (pts[None, :, :] + shifts[:, None, :]).reshape(-1, 2)
    base = np.tile(np.arange(n), len(shifts))
    central = np.repeat(np.arange(len(shifts)), n) == 4
    tri = Delaunay(tiled)
    s = tri.simplices
    e = np.vstack([s[:, [0, 1]], s[:, [1, 2]], s[:, [2, 0]]])
    e = e[central[e[:, 0]] | central[e[:, 1]]]
    sg = sign[y, x]
    opp = sg[base[e[:, 0]]] * sg[base[e[:, 1]]] < 0
    return np.mod(0.5 * (tiled[e[opp, 0]] + tiled[e[opp, 1]]), L)


def box_dimension(snaps, L):
    es = np.array([e for e in (1, 2, 4, 8, 16, 32, 64, 128) if e <= L // 8])
    walls = [wall_points(s) for s in snaps]        # one Delaunay per snapshot
    N = []
    for e in es:
        N.append(np.mean([len(np.unique((w[:, 0] // e).astype(np.int64) * (L // e + 1)
                                        + (w[:, 1] // e).astype(np.int64))) for w in walls]))
    N = np.array(N)
    k = N > 4
    return es[k], -np.gradient(np.log(N[k]), np.log(es[k]))


def xi_from_dimension(es, d, thresh=1.5):
    # the length where the reactive set stops looking 1D and starts looking 2D
    for i in range(len(d) - 1):
        if d[i] < thresh <= d[i + 1]:
            f = (thresh - d[i]) / (d[i + 1] - d[i])
            return float(np.exp(np.log(es[i]) + f * (np.log(es[i + 1]) - np.log(es[i]))))
    return np.nan


# ------------------------------------------------------- mass vs distance to a wall

def load_mass_snaps(L, rho, p, seed=1):
    f = os.path.join(OUT, "massSnap_%s_seed_%d.tsv" % (tag(L, rho, p), seed))
    rows = []
    with open(f) as fh:
        for line in fh:
            v = line.split()
            if len(v) == L:
                rows.append(np.fromiter(map(int, v), dtype=np.int64, count=L))
    n = len(rows) // L
    return np.array(rows[:n * L]).reshape(n, L, L)


def wall_distance(ms, sigma=3.0):
    # coarse-grain the signs to find which domain each site belongs to, then measure
    # how far every site is from the nearest domain wall
    from scipy.ndimage import gaussian_filter, distance_transform_edt
    dom = np.sign(gaussian_filter(np.sign(ms).astype(float), sigma, mode="wrap"))
    wall = np.zeros_like(dom, dtype=bool)
    for ax in (0, 1):
        wall |= dom != np.roll(dom, 1, axis=ax)
        wall |= dom != np.roll(dom, -1, axis=ax)
    return distance_transform_edt(~wall)


def mass_vs_distance(L, rho, p, seed=1, nb=10):
    ms = load_mass_snaps(L, rho, p, seed)
    d_all, m_all = [], []
    for m in ms:
        d = wall_distance(m)
        k = m != 0
        d_all.append(d[k]); m_all.append(np.abs(m[k]))
    d = np.concatenate(d_all); m = np.concatenate(m_all).astype(float)
    e = np.geomspace(1, max(d.max(), 2), nb)
    out = []
    for i in range(len(e) - 1):
        k = (d >= e[i]) & (d < e[i + 1])
        if k.sum() > 200:
            out.append((np.sqrt(e[i] * e[i + 1]), m[k].mean(), k.sum()))
    return np.array(out)


# ------------------------------------------------------------------------- figure

PS_MAIN = (0.5, 0.6, 0.7, 0.8, 0.9)
COL = {0.0: "#8c8c8c", 0.3: "#7570b3", 0.5: "#1b9e77", 0.6: "#66a61e",
       0.65: "#a6a600", 0.7: "#e6ab02", 0.75: "#e08214", 0.8: "#d95f02",
       0.85: "#c4451c", 0.9: "#e7298a", 0.95: "#a6329a", 1.0: "#1f78b4"}


def snaps_for(L, rho, p):
    import plots as P
    old, P.OUT = P.OUT, SERVER
    try:
        return P.load_snaps(L, rho, p)[:6]
    finally:
        P.OUT = old


def figure(L=512, rho=0.2, Ls=(256, 512, 1024)):
    os.makedirs(PLOTS, exist_ok=True)
    fig, ax = plt.subplots(2, 2, figsize=(16, 13))
    a, b, c, d = ax.ravel()

    # (a) both spectra compensated by s^(2 tau_m - 1), with tau_m taken from the SPOT
    # distribution, not fitted here -- so a flat dashed line is a real test of the
    # mean-field law P(min = s) ~ s^(1 - 2 tau_m), not a tautology. Compensating by a
    # fixed s^2 instead only looks flat near p = 1, where tau_m happens to be ~3/2;
    # at p = 0.5, tau_m = 1.70 and the baseline genuinely falls as s^-0.4.
    import plots as P
    for p in PS_MAIN:
        x, yc, ye, _ = spectra(L, rho, p)
        tm, _, _ = P.plateau_fit(*load_hist("chSpot", L, rho, p))
        comp = x ** (2 * tm - 1)
        a.plot(x, yc * comp, "--", color=COL[p], lw=1.6)
        a.plot(x, ye * comp, "-", color=COL[p], lw=2.4,
               label=r"$p=%g\ (\tau_m=%.2f)$" % (p, tm))
    a.set(xscale="log", yscale="log", xlabel="$s$",
          ylabel=r"$s^{\,2\tau_m-1}P(s)$")
    a.plot([], [], "k--", lw=1.6, label="all collisions")
    a.plot([], [], "k-", lw=2.4, label="annihilations")
    a.legend(frameon=False, ncol=1, fontsize=11, loc="lower left",
             borderaxespad=0.4, labelspacing=0.25)

    # (b) their ratio: the annihilation fraction, and its minimum
    for p in PS_MAIN:
        x, _, _, r = spectra(L, rho, p)
        b.plot(x, r, "-o", color=COL[p], lw=2.2, ms=4, label="$p=%g$" % p)
        s = crossover_scale(x, r)
        if np.isfinite(s):
            b.plot([s], [r[list(x).index(s)]], "v", color=COL[p], ms=14, mec="k", mew=1.2)
    b.set(xscale="log", yscale="log", xlabel="$s$", ylabel="$1-q(s)$")
    b.legend(frameon=False, fontsize=13)

    # (c) dimension of the reactive set, and its independence of L
    mk = {256: "^", 512: "o", 1024: "s"}
    for p in (0.0, 0.6, 0.7, 0.8, 0.9, 1.0):
        try:
            es, dd = box_dimension(snaps_for(L, rho, p), L)
        except OSError:
            continue
        # p = 0 is the well-mixed baseline: even there the wall point set looks
        # low-dimensional at small eps, purely from the discreteness of the points
        c.plot(es, dd, "o-", color=COL[p], lw=2.4 if p else 1.6,
               ls="-" if p else "--", ms=6, label="$p=%g$" % p)
    c.axhline(1, color="k", ls=":", lw=1.2)
    c.axhline(2, color="k", ls=":", lw=1.2)
    c.set(xscale="log", xlabel=r"box size $\epsilon$", ylabel=r"$d_f(\epsilon)$", ylim=(0.5, 2.3))
    c.legend(frameon=False, fontsize=13, loc="lower right")

    # (d) does the crossover scale move with the system?
    for LL in Ls:
        ps, ss = [], []
        for p in (0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95):
            try:
                x, _, _, r = spectra(LL, rho, p)
            except OSError:
                continue
            s = crossover_scale(x, r)
            if np.isfinite(s):
                ps.append(p); ss.append(s)
        if ps:
            d.plot(ps, ss, mk[LL] + "-", lw=2.2, ms=9, label="$L=%d$" % LL)
    d.axvspan(0.70, 0.90, color="0.9", zorder=0)   # s* independent of L in here
    d.text(0.80, 0.93, "$s^*$ independent of $L$", transform=d.get_xaxis_transform(),
           ha="center", fontsize=14)
    d.set(yscale="log", xlabel="$p$", ylabel="$s^*$")
    d.legend(frameon=False, loc="lower left")

    for k, lab in zip((a, b, c, d), "abcd"):
        k.text(0.02, 0.97, "(%s)" % lab, transform=k.transAxes, va="top", fontsize=22)
        k.grid(True, which="major", ls="-", lw=0.6, alpha=0.35)
        k.grid(True, which="minor", ls=":", lw=0.4, alpha=0.25)
        k.set_axisbelow(True)
    fig.tight_layout()
    f = os.path.join(PLOTS, "crossover_L_%d_rho_%g.png" % (L, rho))
    fig.savefig(f, dpi=300, bbox_inches="tight")
    print("wrote", f)


if __name__ == "__main__":
    figure()
