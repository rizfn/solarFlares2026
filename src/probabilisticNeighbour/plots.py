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

def tag(L, rho, p):
    return f"L_{g(L)}_rho_{g(rho)}_p_{g(p)}"

def files(kind, L, rho, p):
    fs = sorted(glob.glob(os.path.join(OUT, f"{kind}_{tag(L, rho, p)}_seed_*.tsv")))
    if not fs:
        raise OSError(f"no files for {kind} {tag(L, rho, p)}")
    return fs

def load_hist(kind, L, rho, p):
    # sum histograms over seeds
    total = {}
    for f in files(kind, L, rho, p):
        d = np.loadtxt(f, dtype=np.int64, ndmin=2)
        for s, c in d:
            total[s] = total.get(s, 0) + c
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
        with open(path) as fh:
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


def mle(sizes, counts, xmin=XMIN):
    # Hill estimator for a discrete power law above xmin (the -1/2 is the continuity
    # correction; without it tau is biased high when xmin is small)
    m = sizes >= xmin
    s, c = sizes[m], counts[m]
    if c.sum() < 100:
        return np.nan
    return 1 + c.sum() / np.sum(c * np.log(s / (xmin - 0.5)))


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


def correlation(snaps):
    # conditional sign correlation over occupied pairs, averaged over snapshots
    L = snaps[0].shape[0]
    num = np.zeros((L, L))
    den = np.zeros((L, L))
    for a in snaps:
        s = np.sign(a).astype(float)
        o = (a != 0).astype(float)
        num += np.fft.irfft2(np.abs(np.fft.rfft2(s)) ** 2, s=(L, L))
        den += np.fft.irfft2(np.abs(np.fft.rfft2(o)) ** 2, s=(L, L))
    num = np.fft.fftshift(num); den = np.fft.fftshift(den)
    c = L // 2
    yy, xx = np.indices((L, L))
    r = np.round(np.sqrt((yy - c) ** 2 + (xx - c) ** 2)).astype(int)
    C = np.bincount(r.ravel(), num.ravel()) / np.bincount(r.ravel(), den.ravel())
    return C

def xi_exp(snaps, rmax=6):
    # Decay rate of |C(r)| fitted from r=1 over a short window fixed in lattice units.
    # Neutrality (sum_r C(r) = 0) makes the tail of C an artifact of the box: at fixed p
    # it grows with L, while r<~4 is L-independent. Only the short range is physical, so
    # the window must not follow the tail (a self-consistent window chases it and gives
    # xi proportional to L).
    C = correlation(snaps)
    r = np.arange(1, rmax + 1)
    y = np.abs(C[1:rmax + 1])
    if (y <= 0).any() or abs(C[2]) > abs(C[1]):
        return np.nan
    slope = np.polyfit(r, np.log(y), 1)[0]
    return -1.0 / slope if slope < 0 else np.nan


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


def plot_correlation_length(L=128, rhos=(0.2, 0.4, 0.6, 0.8),
                            ps=(0.0, 0.1, 0.2, 0.3, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475,
                                0.5, 0.525, 0.55, 0.575, 0.6, 0.65, 0.7, 0.8, 0.9, 1.0)):
    fig, ax = plt.subplots(figsize=(8, 6))
    colors = plt.cm.viridis(np.linspace(0, 0.85, len(rhos)))
    for rho, col in zip(rhos, colors):
        xs, ys = [], []
        for p in ps:
            try:
                v = xi_exp(load_snaps(L, rho, p))
            except OSError:
                continue
            if np.isfinite(v):
                xs.append(p); ys.append(v)
        ax.plot(xs, ys, "o-", color=col, label=rf"$\rho={g(rho)}$")
    ax.set_xlabel("neighbour probability $p$")
    ax.set_ylabel(r"correlation length $\xi$")
    ax.legend(frameon=False)
    grid(ax)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTS, "correlationLength", f"correlationLength_L_{g(L)}.png"), dpi=300)
    plt.close(fig)


def plot_histograms(L=128, rhos=(0.2, 0.4, 0.6, 0.8), ps=(0.0, 0.5, 1.0), offset=12):
    colors = plt.cm.viridis(np.linspace(0, 0.85, len(rhos)))
    for p in ps:
        fig, (a1, a2) = plt.subplots(1, 2, figsize=(14, 6))
        for a, kind in ((a1, "spotSize"), (a2, "emission")):
            taus, amps, spans = [], [], []
            for rho, col in zip(rhos, colors):
                sizes, counts = load_hist(kind, L, rho, p)
                tau = mle(sizes, counts)
                x, y = logbin(sizes, counts)
                a.plot(x, y, "o", ms=5, color=col, label=rf"$\rho={g(rho)}$, $\tau={tau:.2f}$")
                f = (x >= XMIN) & (x <= sizes.max() / 20)
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


def plot_exponents_vs_p(L=128, rho=0.2, pc=0.6,
                        ps=(0.0, 0.1, 0.2, 0.3, 0.35, 0.4, 0.45, 0.5, 0.525, 0.55,
                            0.575, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0)):
    # the cascade exponent and the segregation order parameter share one x axis:
    # tau slides smoothly across the whole range, m switches on at p_c
    xs, te, ts, ms = [], [], [], []
    for p in ps:
        try:
            e = mle(*load_hist("emission", L, rho, p))
            s = mle(*load_hist("spotSize", L, rho, p))
            m = segregation(load_snaps(L, rho, p))
        except OSError:
            continue
        xs.append(p); te.append(e); ts.append(s); ms.append(m)

    fig, ax = plt.subplots(figsize=(9, 6.5))
    ax.axvline(pc, color="grey", ls=":", lw=2)
    ax.text(pc - 0.015, 3.35, rf"$p_c\approx{g(pc)}$", color="grey", ha="right", fontsize=17)
    # mean-field windows: tau_m in [3/2, 2] as q runs 1 -> 1/2, and tau_s = 2 tau_m - 1
    ax.axhspan(1.5, 2.0, color="#ef8a62", alpha=0.15, lw=0)
    ax.axhspan(2.0, 3.0, color="#b2182b", alpha=0.10, lw=0)
    ax.plot(xs, te, "o-", color="#b2182b", ms=8, lw=2.5, label=r"emission $\tau_s$")
    ax.plot(xs, ts, "^-", color="#ef8a62", ms=8, lw=2.5, label=r"spot $\tau_m$")
    ax.set_xlabel("neighbour probability $p$")
    ax.set_ylabel(r"power-law exponent $\tau$", color="#b2182b")
    ax.tick_params(axis="y", colors="#b2182b")
    ax.set_ylim(1.3, 3.5)
    ax.legend(frameon=False, loc="lower left")
    grid(ax)

    a2 = ax.twinx()
    a2.plot(xs, ms, "s--", color="#2166ac", ms=8, lw=2.5)
    a2.set_ylabel("segregation $m$", color="#2166ac")
    a2.tick_params(axis="y", colors="#2166ac")
    a2.set_ylim(-0.03, max(ms) * 1.15)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTS, "exponents", f"exponents_L_{g(L)}_rho_{g(rho)}.png"), dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    plot_snapshots()
    plot_correlation_length()
    plot_histograms()
    plot_exponents_vs_p()
