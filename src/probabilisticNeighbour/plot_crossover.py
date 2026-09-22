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
# Why the U has that shape is not established. s* itself is well behaved: it is
# independent of L for 0.70 <= p <= 0.90, and correlatedOrBipole shows it is a function
# of q alone, the same under two very different injection rules. The depth is not -- the
# minimum keeps falling with L (roughly L^-0.7 over 256-1024), so the U deepens with
# system size while staying put.
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


# --------------------------------------------------------------- what sets the size

def age_profile(L, rho, p, seeds=range(1, 9)):
    """Mean age and distance walked by the emitting spot, against the size it emits.
    Pooled over seeds; bins are the fixed log grid written by emissionAge.cpp."""
    acc = {}
    for s in seeds:
        fn = os.path.join(OUT, "chAge_%s_seed_%d.tsv" % (tag(L, rho, p), s))
        if not os.path.exists(fn) or os.path.getsize(fn) < 20:
            continue
        d = np.loadtxt(fn)
        if d.ndim == 1:
            d = d[None, :]
        for size, n, age, disp, _rms in d:
            k = round(np.log10(size), 4)
            acc.setdefault(k, np.zeros(3))
            acc[k] += [n, n * age, n * disp]
    ks = np.array(sorted(acc))
    v = np.array([acc[k] for k in ks])
    keep = v[:, 0] > 500
    return 10 ** ks[keep], v[keep, 1] / v[keep, 0], v[keep, 2] / v[keep, 0]


def origin_fraction(L, rho, p, seeds=range(1, 9), perdec=3):
    """Share of emissions of size s produced by a spot that was injected at a random
    site rather than beside its own sign."""
    def acc(kind):
        a = {}
        for s in seeds:
            fn = os.path.join(OUT, "%s_%s_seed_%d.tsv" % (kind, tag(L, rho, p), s))
            if not os.path.exists(fn) or os.path.getsize(fn) < 20:
                continue
            d = np.loadtxt(fn)
            if d.ndim == 1:
                d = d[None, :]
            for k, v in d:
                a[k] = a.get(k, 0.0) + v
        return a
    ao, ar = acc("chEmisOwn"), acc("chEmisRand")
    ks = np.array(sorted(set(ao) | set(ar)))
    e = np.geomspace(1, ks.max() + 1, int(perdec * np.log10(ks.max())) + 2)
    ho, _ = np.histogram(ks, e, weights=np.array([ao.get(k, 0.0) for k in ks]))
    hr, _ = np.histogram(ks, e, weights=np.array([ar.get(k, 0.0) for k in ks]))
    tot = ho + hr
    k = tot > 2000
    return np.sqrt(e[:-1] * e[1:])[k], (hr / np.maximum(tot, 1))[k]


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

    # (c) what actually sets the emitted size: how long the spot lived and how far it
    # walked before it met an opposite sign. Small emissions are spots that die where
    # they were born; large ones have walked, and eaten, for a thousand times longer.
    c2 = c.twinx()
    for p in (0.6, 0.7, 0.8, 0.9):
        try:
            x, age, disp = age_profile(L, rho, p)
        except (OSError, ValueError):
            continue
        c.plot(x, age, "-", color=COL[p], lw=2.4, label="$p=%g$" % p)
        c2.plot(x, disp, "--", color=COL[p], lw=1.6)
    c.set(xscale="log", yscale="log", xlabel="$s$", ylabel="age at emission (sweeps)")
    c2.set_yscale("log"); c2.set_ylabel("distance walked", fontsize=18)
    c.plot([], [], "k-", lw=2.4, label="age")
    c.plot([], [], "k--", lw=1.6, label="distance")
    c.legend(frameon=False, fontsize=13, loc="upper left", bbox_to_anchor=(0.10, 0.98))

    # (d) does the crossover scale move with the system?
    mk = {256: "^", 512: "o", 1024: "s"}
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
