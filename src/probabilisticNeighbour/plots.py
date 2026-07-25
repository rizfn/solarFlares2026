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

def g(v):
    return f"{v:g}"

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
    # all snapshots across seeds; loadtxt skips the blank separators, so reshape into LxL
    snaps = []
    for path in files("snapshots", L, rho, p):
        snaps.extend(np.loadtxt(path, dtype=np.int8).reshape(-1, L, L))
    return snaps

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

def xi(snaps):
    # domain size = first zero crossing of the sign correlation, capped at L/2
    L = snaps[0].shape[0]
    lim = L // 2
    C = correlation(snaps)
    z = np.where(C[:lim] < 0)[0]
    return z[0] if len(z) else lim


def xi_exp(snaps, c=2.5):
    # decay rate of |C(r)| fitted from r=1; window is self-consistently ~c*xi, so it
    # tracks xi rather than L (a window tied to L makes different boxes disagree)
    C = correlation(snaps)
    L = snaps[0].shape[0]

    # neutrality forces C to change sign; fit only the leading constant-sign stretch
    rsign = L // 4
    for r in range(2, L // 4):
        if np.sign(C[r]) != np.sign(C[1]):
            rsign = r - 1
            break
    if rsign < 2:  # only C(1) survives; xi is unmeasurable at the sign-change crossover
        return np.nan
    if abs(C[2]) > abs(C[1]):  # C(1) sitting on its zero crossing: not a decaying tail
        return np.nan

    def fit(rmax):
        rmax = min(rmax, rsign)
        r = np.arange(1, rmax + 1)
        y = np.abs(C[1:rmax + 1])
        slope = np.polyfit(r, np.log(y), 1)[0]
        return -1.0 / slope if slope < 0 else np.nan

    v = fit(6)
    for _ in range(20):
        if not np.isfinite(v):
            return np.nan
        nv = fit(int(np.clip(round(c * v), 2, L // 4)))
        if not np.isfinite(nv):
            return np.nan
        converged = abs(nv - v) < 1e-3
        v = nv
        if converged:
            break
    # a length beyond L/4 is not resolvable in the box; don't report it
    return v if v < L / 4 else np.nan


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
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTS, "correlationLength", f"correlationLength_L_{g(L)}.png"), dpi=300)
    plt.close(fig)


def plot_critical_scaling(rho=0.2, Ls=(64, 128, 256, 512), fitL=512,
                          ps=(0.30, 0.36, 0.40, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48,
                              0.49, 0.50, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58,
                              0.60, 0.61, 0.64, 0.67, 0.70)):
    data = {}
    for L in Ls:
        xs, ys = [], []
        for p in ps:
            try:
                v = xi_exp(load_snaps(L, rho, p))
            except OSError:
                continue
            if np.isfinite(v):
                xs.append(p); ys.append(v)
        data[L] = (np.array(xs), np.array(ys))

    # xi ~ (p_c - p)^-nu: fit on the largest box, where xi has the most room to grow
    fx, fy = data[fitL]
    rise = fx >= fx[np.argmin(fy)]           # drop the flat disordered plateau
    rise &= fy < fitL / 8                    # and points pinned near the box ceiling
    fx, fy = fx[rise], fy[rise]
    # a power law cannot be fitted through a saturating curve; refuse rather than rail
    fitted = len(fx) >= 4
    if fitted:
        pcs, nus = np.arange(fx.max() + 0.005, 0.90, 0.002), np.arange(0.3, 3.0, 0.01)
        best = (1e18, 0.0, 0.0)
        for pc in pcs:
            for nu in nus:
                res = np.log(fy) + nu * np.log(pc - fx)
                c = np.sum((res - res.mean()) ** 2)
                if c < best[0]:
                    best = (c, pc, nu)
        _, pc, nu = best
        amp = np.exp(np.mean(np.log(fy) + nu * np.log(pc - fx)))
        if pc >= pcs[-1] or nu >= nus[-1]:   # railed: the data do not constrain the fit
            fitted = False
    if not fitted:
        pc = data[fitL][0][np.argmax(data[fitL][1] > fitL / 8)]  # crossing as a rough p_c

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(14, 6))
    colors = plt.cm.plasma(np.linspace(0, 0.8, len(Ls)))
    for L, col in zip(Ls, colors):
        x, y = data[L]
        a1.plot(x, y, "o-", color=col, label=f"$L={L}$")
        m = x < pc
        a2.plot(pc - x[m], y[m], "o-", color=col, label=f"$L={L}$")
    a1.axvline(pc, color="grey", ls=":", lw=1.5)
    a1.set_yscale("log")
    a1.set_xlabel("neighbour probability $p$"); a1.set_ylabel(r"$\xi$")
    a1.legend(frameon=False)
    a1.set_title(rf"$p_c\approx{pc:.2f}$")

    if fitted:
        t = np.array([(pc - fx).min(), (pc - fx).max()])
        a2.plot(t, amp * t ** -nu, "k--", lw=1.5, label=rf"$(p_c-p)^{{-{nu:.2f}}}$")
    a2.set_xscale("log"); a2.set_yscale("log")
    a2.set_xlabel(r"$p_c - p$"); a2.set_ylabel(r"$\xi$")
    a2.set_title(rf"$\nu\approx{nu:.2f}$" if fitted else "no fit: too few unsaturated points")
    a2.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTS, "criticalScaling", f"criticalScaling_rho_{g(rho)}.png"), dpi=300)
    plt.close(fig)


def plot_histograms(L=128, rhos=(0.2, 0.4, 0.6, 0.8), ps=(0.0, 0.5, 1.0)):
    colors = plt.cm.viridis(np.linspace(0, 0.85, len(rhos)))
    for p in ps:
        fig, (a1, a2) = plt.subplots(1, 2, figsize=(14, 6))
        for rho, col in zip(rhos, colors):
            x, y = logbin(*load_hist("spotSize", L, rho, p))
            a1.plot(x, y, "o", ms=5, color=col, label=rf"$\rho={g(rho)}$")
            x, y = logbin(*load_hist("emission", L, rho, p))
            a2.plot(x, y, "o", ms=5, color=col, label=rf"$\rho={g(rho)}$")
        for a, xl, yl in [(a1, "spot size $m$", "$n(m)$"), (a2, "emission size $s$", "$P(s)$")]:
            a.set_xscale("log"); a.set_yscale("log"); a.set_xlabel(xl); a.set_ylabel(yl)
            a.legend(frameon=False)
        fig.suptitle(f"$p={g(p)}$", fontsize=20)
        fig.tight_layout()
        fig.savefig(os.path.join(PLOTS, "histograms", f"histograms_L_{g(L)}_p_{g(p)}.png"), dpi=300)
        plt.close(fig)


if __name__ == "__main__":
    plot_snapshots()
    plot_correlation_length()
    plot_critical_scaling()
    plot_histograms()
