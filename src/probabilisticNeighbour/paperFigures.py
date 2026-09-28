# Publication figures for the probabilistic-neighbour model.
#
# Data: the server sweep in ../probabilisticNeighbourServer/outputs (7400 runs, L up to
# 1024) for the lattice, and ./outputs for the well-mixed Monte Carlo. Histograms are
# cached per seed on a fixed log grid (outputs/figcache), so every figure after the first
# run is a few seconds.
#
#   fig1  lattice snapshots at two densities, with the correlation length xi(p)
#   fig2  the spectra: n(m), P(s), both compensated, and the break s*(p)
#   fig3  exponents against the closed form cos(pi tau_m) = (1-q)/q
#   fig4  density scaling: rho = 0.05-0.9 on both sides of p_c, L = 512
#   fig5  the branching ratio q(p), which is what the closed form actually needs
#
# Colours. Four fixed ones, everything else built from them: red and blue are the two
# signs of a spot, exactly as on a lattice snapshot; C_SPOT carries every spot-size
# quantity and C_EMIS every flare quantity, so a panel's hue says which of the two
# distributions it is about.

import os
import glob
import gzip
import colorsys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

DIR = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(DIR, "outputs")                                  # well-mixed MC
SRV = os.path.join(DIR, "..", "probabilisticNeighbourServer", "outputs")   # lattice sweep
CACHE = os.path.join(OUT, "figcache")
PLOTS = os.path.join(DIR, "plots", "paper")

plt.rcParams.update({
    "font.size": 24, "axes.labelsize": 32, "xtick.labelsize": 25, "ytick.labelsize": 25,
    "legend.fontsize": 22, "axes.linewidth": 1.6, "xtick.major.width": 1.6,
    "ytick.major.width": 1.6, "xtick.minor.width": 1.1, "ytick.minor.width": 1.1,
    "xtick.major.size": 9, "ytick.major.size": 9, "xtick.minor.size": 5,
    "ytick.minor.size": 5, "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True, "figure.dpi": 110,
})

C_POS = "#b2182b"      # a + spot on the lattice
C_NEG = "#2166ac"      # a - spot
C_SPOT = "#ad8350"     # anything about spot sizes
C_EMIS = "#5f0f40"     # anything about flares
C_REF = "#666666"      # reference power laws, theory curves, identities


def tint(base, f):
    # f < 0 mixes base toward white, f > 0 toward black
    c = np.array(mcolors.to_rgb(base))
    return tuple(c + (1 - c) * (-f)) if f < 0 else tuple(c * (1 - f))


def hramp(base, stops, n=256):
    # Sequential colormap around one colour, stepping hue as well as lightness: stops are
    # (hue shift in degrees, lightness, saturation factor) in HLS. Lightness alone made
    # neighbouring curves too alike; a few tens of degrees of hue separates them.
    h, _, sat = colorsys.rgb_to_hls(*mcolors.to_rgb(base))
    cols = [colorsys.hls_to_rgb((h + dh / 360) % 1, lt, min(1.0, sat * fs)) for dh, lt, fs in stops]
    return mcolors.LinearSegmentedColormap.from_list("", cols, N=n)


CM_SPOT = hramp(C_SPOT, [(+24, 0.80, 1.25), (+8, 0.64, 1.1), (0, 0.50, 1.0), (-16, 0.37, 1.1)])   # gold -> tan -> sienna
CM_EMIS = hramp(C_EMIS, [(+38, 0.78, 0.85), (+15, 0.56, 0.9), (0, 0.36, 1.0), (-20, 0.24, 0.95)])  # salmon -> magenta -> violet

# Fixed binning grid shared by every cached histogram: exact integers to 10, then 20 bins
# per decade. Integer resolution at small s avoids the empty-bin comb of a pure log grid.
EDGES = np.unique(np.concatenate([np.arange(1, 11), np.geomspace(10, 1e10, 20 * 9 + 1)]))
CTR = np.sqrt(EDGES[:-1] * EDGES[1:])
WID = np.diff(EDGES)


def g(v):
    return f"{v:g}"


def grid(ax):
    ax.grid(True, which="major", ls=":", lw=0.9, alpha=0.5)
    ax.set_axisbelow(True)


def _bin_seeds(pattern):
    # one row of binned counts per seed file; numpy reads .tsv and .tsv.gz alike
    fs = sorted(glob.glob(pattern) + glob.glob(pattern + ".gz"))
    rows = []
    for f in fs:
        if os.path.getsize(f) < 20:
            continue
        d = np.loadtxt(f, dtype=np.int64, ndmin=2)
        if d.size == 0:
            continue
        h, _ = np.histogram(d[:, 0].astype(float), bins=EDGES, weights=d[:, 1].astype(float))
        rows.append(h)
    if not rows:
        raise OSError(f"no usable files for {pattern}")
    return np.array(rows)


def hist(kind, L, rho, p, root=None):
    # binned counts, shape (nseeds, nbins); `root` picks the server sweep (default) or
    # the local channel runs, which are the only ones covering rho = 0.05 to 0.6
    key = os.path.join(CACHE, f"{kind}_L_{g(L)}_rho_{g(rho)}_p_{g(p)}.npy")
    if os.path.exists(key):
        return np.load(key)
    h = _bin_seeds(os.path.join(root or SRV, f"{kind}_L_{g(L)}_rho_{g(rho)}_p_{g(p)}_seed_*.tsv"))
    os.makedirs(CACHE, exist_ok=True)
    np.save(key, h)
    return h


def wm_hist(kind, N, q):
    # binned counts for the well-mixed Monte Carlo
    key = os.path.join(CACHE, f"{kind}_N_{N:d}_q_{g(q)}.npy")
    if os.path.exists(key):
        return np.load(key)
    h = _bin_seeds(os.path.join(OUT, f"{kind}_N_{N:d}_q_{g(q)}_seed_*.tsv"))
    os.makedirs(CACHE, exist_ok=True)
    np.save(key, h)
    return h


def branching_ratio(L, rho, p, nseed=8):
    # q: the fraction of occupied nearest-neighbour pairs that share a sign. Every
    # collision is between neighbours, so this is exactly the branching ratio of the
    # kinetic theory -- 1/2 on a mixed surface, -> 1 once the signs segregate.
    key = os.path.join(CACHE, f"q_L_{g(L)}_rho_{g(rho)}_p_{g(p)}.npy")
    if os.path.exists(key):
        return float(np.load(key))
    fs = sorted(glob.glob(os.path.join(SRV, f"snapshots_L_{g(L)}_rho_{g(rho)}_p_{g(p)}_seed_*.tsv"))
                + glob.glob(os.path.join(SRV, f"snapshots_L_{g(L)}_rho_{g(rho)}_p_{g(p)}_seed_*.tsv.gz")))[:nseed]
    same = opp = 0
    for f in fs:
        rows = []
        with (gzip.open(f, "rt") if f.endswith(".gz") else open(f)) as fh:
            for line in fh:
                v = line.split()
                if len(v) == L:
                    rows.append(np.fromiter(map(int, v), dtype=np.int32, count=L))
        n = len(rows) // L
        if not n:
            continue
        a = np.sign(np.array(rows[:n * L], dtype=np.int32).reshape(n, L, L)).astype(np.int8)
        for b in (np.roll(a, 1, 1), np.roll(a, 1, 2)):
            pr = a.astype(np.int32) * b
            same += int((pr > 0).sum()); opp += int((pr < 0).sum())
    if same + opp == 0:
        raise OSError(f"no snapshots for L={L} rho={rho} p={p}")
    q = same / (same + opp)
    os.makedirs(CACHE, exist_ok=True)
    np.save(key, np.array(q))
    return q


NINT = 9                            # bins 0..8 are the exact integers 1..9


def rebin(h, k):
    # merge groups of k log-bins for plotting. The cache grid is 20 per decade, fine
    # enough for fitting but narrower than one integer near m ~ 10, which shows up as a
    # comb; k=4 gives 5 per decade and removes it without touching the integer bins.
    lo, hi = h[:NINT], h[NINT:]
    n = (len(hi) // k) * k
    e = np.concatenate([EDGES[:NINT + 1], EDGES[NINT + k::k][: n // k]])
    return np.concatenate([lo, hi[:n].reshape(-1, k).sum(axis=1)]), e


def density(h, minc=1, k=1):
    # normalized probability density from binned counts, dropping under-sampled bins
    tot = h.sum()
    if k > 1:
        h, e = rebin(h, k)
        ctr, wid = np.sqrt(e[:-1] * e[1:]), np.diff(e)
    else:
        ctr, wid = CTR, WID
    ok = h >= minc
    return ctr[ok], (h / wid / tot)[ok]


def available(kind, L, rho, root=None):
    root = root or SRV
    fs = (glob.glob(os.path.join(root, f"{kind}_L_{g(L)}_rho_{g(rho)}_p_*_seed_*.tsv"))
          + glob.glob(os.path.join(root, f"{kind}_L_{g(L)}_rho_{g(rho)}_p_*_seed_*.tsv.gz")))
    return sorted({float(f.split("_p_")[1].split("_seed")[0]) for f in fs})


def slope(h, lo, hi, minc=20):
    # weighted log-log slope of the binned density over [lo, hi]; var(log y) ~ 1/counts
    sel = (CTR >= lo) & (CTR <= hi) & (h >= minc)
    if sel.sum() < 4:
        return np.nan
    x, y, w = np.log(CTR[sel]), np.log((h / WID)[sel]), np.sqrt(h[sel])
    return -np.polyfit(x, y, 1, w=w)[0]


def window(h, lo=60.0, ftop=3e-4):
    # Scaling window [lo, hi]: hi is the quantile above which a fraction ftop of the mass
    # lies, not the largest event seen. A single extreme value jumps around between
    # parameter points and made tau(p) visibly jagged when it set the upper edge.
    cc = np.cumsum(h[::-1])[::-1] / h.sum()
    ok = cc >= ftop
    hi = CTR[ok].max() if ok.any() else CTR[-1]
    return lo, max(hi, 10 * lo)


def tau_seeds(H, lo=60.0, ftop=3e-4, **kw):
    # one exponent per seed, so the spread across seeds is an honest error bar
    a, b = window(H.sum(axis=0), lo, ftop)
    t = np.array([slope(h, a, b, **kw) for h in H])
    t = t[np.isfinite(t)]
    return t


def tau_pm(H, **kw):
    t = tau_seeds(H, **kw)
    if len(t) == 0:
        return np.nan, np.nan
    return t.mean(), (t.std(ddof=1) / np.sqrt(len(t)) if len(t) > 1 else 0.0)


def tau_closed(q):
    # Sub-leading balance of the rate equation: cos(pi tau_m) = (1-q)/q on tau_m in
    # [3/2, 2]. No power-law solution at all for q < 1/2 (renormalization_concise.md).
    r = (1 - np.asarray(q, float)) / np.asarray(q, float)
    return np.where(r <= 1, 2 - np.arccos(np.clip(r, -1, 1)) / np.pi, np.nan)


def kink(h, smooth=2):
    # s*: where the compensated density s^2 P(s) turns, i.e. where the local slope
    # crosses 2. Below p_c the curve falls monotonically and there is no turn.
    ok = h >= 50
    x, y = CTR[ok], (h / WID)[ok] * CTR[ok] ** 2
    if len(x) < 8:
        return np.nan
    ly = np.convolve(np.log(y), np.ones(2 * smooth + 1) / (2 * smooth + 1), mode="valid")
    lx = x[smooth:len(x) - smooth]
    i = int(np.argmin(ly))
    if i == 0 or i == len(ly) - 1:
        return np.nan                          # monotone: no interior minimum
    return float(lx[i])


def branches(h, smin=20.0, frac=0.3, ftop=1e-6, mindec=1.0):
    # Above p_c the emission density is two power laws joined at s*. Fit them on either
    # side of the turn, keeping a factor 1/frac clear of it so neither fit straddles the
    # bend, and return (tau_lo, tau_hi, s*). The turn is only accepted when at least
    # `mindec` decades of upper branch lie above it: below p_c the compensated curve does
    # have a shallow minimum, but it sits on the finite-size cutoff, not on a second
    # scaling regime, and calling that a kink would make p_c look smaller than it is.
    st = kink(h)
    if not np.isfinite(st):
        return np.nan, np.nan, np.nan
    cc = np.cumsum(h[::-1])[::-1] / h.sum()
    top = CTR[cc >= ftop].max()
    if top < st * 10 ** mindec:
        return np.nan, np.nan, np.nan
    lo = slope(h, smin, frac * st) if frac * st > 3 * smin else np.nan
    return lo, slope(h, st / frac, top), st




def upper_tau(H, mindec=1.0, frac=0.3, ftop=1e-6):
    # Mean and standard error of the upper-branch exponent. Every seed is fitted over the
    # same window, taken from the pooled histogram: re-detecting the break per seed makes
    # the window jump around and inflates the spread. Returns nan unless the branch spans
    # at least `mindec` decades -- a short fit still sits on the turn and reads high,
    # which is the artefact crossover.md warns about.
    h = H.sum(axis=0)
    st = kink(h)
    if not np.isfinite(st):
        return np.nan, np.nan
    cc = np.cumsum(h[::-1])[::-1] / h.sum()
    top = CTR[cc >= ftop].max()
    if top < st * 10 ** mindec:
        return np.nan, np.nan
    v = np.array([slope(x, st / frac, top) for x in H])
    v = v[np.isfinite(v)]
    if not len(v):
        return np.nan, np.nan
    return v.mean(), (v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else 0.0)


def mf_steady(q, M=4096, J=1.0, dt=0.05, steps=400000, tol=1e-11):
    # Steady state of the rate equation by forward Euler with FFT convolutions:
    #   dn_m/dt = (q/2) sum_{i+j=m} n_i n_j + (1-q) sum_j n_{m+j} n_j - n_m N + J d_{m,1}
    # Mass coagulating past the cutoff M leaves the system, which is what keeps the
    # cascade steady.
    n = np.zeros(M + 1)
    n[1] = 1.0
    nfft = 2 * (M + 1)
    for it in range(steps):
        N = n.sum()
        fn = np.fft.rfft(n, nfft)
        conv = np.fft.irfft(fn * fn, nfft)              # sum_{i+j=m} n_i n_j
        acf = np.fft.irfft(np.abs(fn) ** 2, nfft)       # sum_j n_j n_{j+m}
        dn = 0.5 * q * conv[:M + 1] + (1 - q) * acf[:M + 1] - n * N
        dn[1] += J
        nn = np.maximum(n + dt * dn, 0.0)
        nn[0] = 0.0
        if it % 5000 == 0 and it > 0 and np.abs(nn - n).sum() / (N + 1e-12) < tol:
            return nn
        n = nn
    return n


def mf_tau(q, lo=20, hi=512, **kw):
    # fitted well inside the cutoff, which bends the tail down
    n = mf_steady(q, **kw)
    m = np.arange(len(n))
    sel = (m >= lo) & (m <= hi) & (n > 0)
    return -np.polyfit(np.log(m[sel]), np.log(n[sel]), 1)[0]


def mf_curve(qs=np.arange(0.52, 0.99, 0.04)):
    key = os.path.join(CACHE, "rateEquation.npy")
    if os.path.exists(key):
        return np.load(key)
    out = np.array([[q, mf_tau(q)] for q in qs])
    os.makedirs(CACHE, exist_ok=True)
    np.save(key, out)
    return out


def wm_qs(N=200000, kind="wmSpotSize"):
    fs = glob.glob(os.path.join(OUT, f"{kind}_N_{N:d}_q_*_seed_*.tsv"))
    return sorted({float(f.split("_q_")[1].split("_seed")[0]) for f in fs})


def lattice_tau(L=512, rho=0.2):
    # (q, tau_m, err, tau_s, err, p) for every p of the sweep, with q measured from the
    # snapshots of the same runs. Above p_c the flare spectrum has two branches and no
    # single tau_s, so the upper-branch exponent is used there and flagged.
    rows = []
    for p in available("spotSize", L, rho):
        q = branching_ratio(L, rho, p)
        tm, em = tau_pm(hist("spotSize", L, rho, p))
        H = hist("emission", L, rho, p)
        broken = np.isfinite(branches(H.sum(axis=0))[2])
        ts, es = upper_tau(H) if broken else tau_pm(H, lo=20.0)
        rows.append((q, tm, em, ts, es, p, float(broken)))
    return np.array(rows)




RHOS_SW = (0.2, 0.4, 0.6, 0.8)                  # the server sweep: four densities at L = 128
RHOS_D = (0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9)   # run_density.sh, L = 512


def seq(vals, cmap, lo=0.0, hi=1.0):
    # one colour per value along a ramp
    v = np.asarray(vals, float)
    f = (v - v.min()) / max(v.max() - v.min(), 1e-12)
    return cmap(lo + f * (hi - lo))


def cbar(fig, ax, vals, cmap, label):
    # colorbar standing in for a legend with too many entries
    sm = plt.cm.ScalarMappable(norm=mcolors.Normalize(min(vals), max(vals)), cmap=cmap)
    cb = fig.colorbar(sm, ax=ax, pad=0.025, aspect=24)
    cb.ax.set_title(label, fontsize=30, pad=12)      # horizontal: a rotated one collides
    cb.ax.tick_params(labelsize=23)                  # with the next panel's y label
    return cb


def guide_on(ax, x, y, x0, dec, expo, off=0.12, label=None, lab_at=0.6, dy=0.16,
             color=C_REF, ls="--", lw=2.6):
    # Reference power law anchored a factor `off` from the curve (x, y) at x = x0, so it
    # sits in clear space beside the data rather than crossing it.
    y0 = np.exp(np.interp(np.log(x0), np.log(x), np.log(y))) * off
    x1 = x0 * 10 ** dec
    ax.plot([x0, x1], [y0, y0 * 10 ** (-dec * expo)], ls=ls, color=color, lw=lw, zorder=2)
    if label:
        xl = x0 * 10 ** (dec * lab_at)
        ax.text(xl, y0 * (xl / x0) ** (-expo) * dy, label, fontsize=27, ha="left",
                va="bottom" if dy > 1 else "top", color=color)   # dy > 1: above the line


def on_curve(x, y, s):
    return np.exp(np.interp(np.log(s), np.log(x), np.log(y)))


def save(fig, name):
    os.makedirs(PLOTS, exist_ok=True)
    fig.savefig(os.path.join(PLOTS, name), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("wrote", os.path.join("plots", "paper", name), flush=True)


def read_snaps(L, rho, p, nseed=8):
    # every complete L x L sign lattice in the first nseed snapshot files, as int8
    fs = sorted(glob.glob(os.path.join(SRV, f"snapshots_L_{g(L)}_rho_{g(rho)}_p_{g(p)}_seed_*.tsv"))
                + glob.glob(os.path.join(SRV, f"snapshots_L_{g(L)}_rho_{g(rho)}_p_{g(p)}_seed_*.tsv.gz")))[:nseed]
    out = []
    for f in fs:
        rows = []
        with (gzip.open(f, "rt") if f.endswith(".gz") else open(f)) as fh:
            for line in fh:
                v = line.split()
                if len(v) == L:
                    rows.append(np.fromiter(map(int, v), dtype=np.int64, count=L))
        n = len(rows) // L
        if n:
            out.extend(np.sign(np.array(rows[:n * L]).reshape(n, L, L)).astype(np.int8))
    if not out:
        raise OSError(f"no snapshots for L={L} rho={rho} p={p}")
    return out


def radial_corr(snaps, L):
    # C(r): sign correlation of occupied pairs at separation r, conditioned on both sites
    # being occupied, by FFT and radially averaged
    num = np.zeros((L, L)); den = np.zeros((L, L))
    for a in snaps:
        sg = a.astype(float); oc = (a != 0).astype(float)
        num += np.fft.irfft2(np.abs(np.fft.rfft2(sg)) ** 2, s=(L, L))
        den += np.fft.irfft2(np.abs(np.fft.rfft2(oc)) ** 2, s=(L, L))
    num = np.fft.fftshift(num); den = np.fft.fftshift(den)
    c = L // 2
    yy, xx = np.indices((L, L))
    r = np.round(np.hypot(yy - c, xx - c)).astype(int).ravel()
    return np.bincount(r, num.ravel()) / np.bincount(r, den.ravel())


def corr_length(L, rho, p, rmax=6):
    # xi from the short range of C(r): an exponential fit over r = 1..6.
    # Only the short range is physical. The lattice holds equal + and -, so sum_r C(r) = 0
    # and the tail is a finite-size artefact; any estimator that follows the tail returns
    # xi ~ L. In the segregated phase the true length is the box, and this estimator,
    # blind to r > 6, keeps growing without saturating at it -- read it as a monotone
    # measure of how ordered the surface is, not as a domain size.
    # Only the contiguous r >= 1 where C clears three times the noise enter the fit, the
    # noise being C measured on the same snapshots with the signs shuffled. If C(1) is
    # itself noise, xi = 0: at p = 0 there is no correlation to measure.
    key = os.path.join(CACHE, f"xi3_L_{g(L)}_rho_{g(rho)}_p_{g(p)}.npy")
    if os.path.exists(key):
        return float(np.load(key))
    snaps = read_snaps(L, rho, p)
    rng = np.random.default_rng(0)
    shuf = []
    for a in snaps:
        t = a.ravel().copy()
        occ = t != 0
        v = t[occ]; rng.shuffle(v); t[occ] = v
        shuf.append(t.reshape(L, L))
    C = radial_corr(snaps, L)
    floor = 3 * np.sqrt(np.mean(radial_corr(shuf, L)[1:rmax + 1] ** 2))
    rr = np.arange(1, rmax + 1)
    bad = np.nonzero(C[rr] <= floor)[0]
    rr = rr[:bad[0]] if len(bad) else rr
    if len(rr) == 0:
        xi = 0.0
    elif len(rr) == 1:
        xi = np.nan
    else:
        xi = -1.0 / np.polyfit(rr, np.log(C[rr]), 1)[0]
    os.makedirs(CACHE, exist_ok=True)
    np.save(key, np.array(xi))
    return xi


def fig1_snapshots(L=128, rhos=(0.2, 0.6), ps=(0.0, 0.5, 1.0)):
    # The surface itself: + spots red, - spots blue, empty sites white, one row per
    # density. Correlating the injection turns a salt-and-pepper mixture into domains of
    # one sign. The band below is the correlation length for both rows, with the snapshot
    # p values marked; the denser surface orders at slightly smaller p.
    cmap = mcolors.ListedColormap([C_NEG, "white", C_POS])
    nr, nc = len(rhos), len(ps)
    fig = plt.figure(figsize=(5.2 * nc + 1.2, 5.2 * nr + 3.8))
    gs = fig.add_gridspec(nr + 1, nc, height_ratios=[1] * nr + [0.62],
                          hspace=0.07, wspace=0.05)
    for i, r in enumerate(rhos):
        for j, p in enumerate(ps):
            ax = fig.add_subplot(gs[i, j])
            ax.imshow(read_snaps(L, r, p, nseed=1)[-1], cmap=cmap, vmin=-1, vmax=1,
                      interpolation="nearest")
            ax.set_xticks([]); ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_color(C_REF); sp.set_linewidth(1.2)
            if i == 0:
                ax.set_title(rf"$p={g(p)}$", fontsize=32, pad=10)
            if j == 0:
                ax.set_ylabel(rf"$\rho={g(r)}$", fontsize=32, labelpad=10)

    band = fig.add_subplot(gs[nr, :])
    for r, col, mk in zip(rhos, seq(rhos, CM_EMIS, lo=0.25, hi=0.95), ("o", "s")):
        pp = available("snapshots", L, r)
        xi = np.array([corr_length(L, r, p) for p in pp])
        ok = np.isfinite(xi)
        band.plot(np.array(pp)[ok], xi[ok], mk + "-", color=col, ms=11,
                  lw=2.8, mec="white", mew=1.2, label=rf"$\rho={g(r)}$")
    for p in ps:
        band.axvline(p, color=C_REF, ls=":", lw=1.8, zorder=0)
    band.set_xlim(-0.02, 1.02)
    band.set_xlabel("neighbour probability $p$")
    band.set_ylabel(r"$\xi$")
    band.legend(frameon=False, loc="upper left", ncol=2)
    grid(band)
    save(fig, f"fig1_snapshots_L_{g(L)}.png")


def fig2_spectra(L=512, rho=0.2, ps=None, Ls=(256, 512, 1024)):
    # The whole phenomenology in one figure, spot quantities in C_SPOT and flare
    # quantities in C_EMIS; the top row is the distributions, the bottom row the same
    # divided by a reference power law so that differences in slope become visible.
    #
    # Top left: the spot spectrum is a power law at every p whose exponent slides from 2
    #     down to the aggregation value 3/2, reached at p ~ 0.7 and then held.
    # Top right: the flare spectrum does the same below p_c, but above it breaks in two.
    # Bottom left: n(m) m^{3/2}. tau_m = 3/2 is now flat and tau_m = 2 has slope -1/2, so
    #     the slide becomes a fan and the p >= 0.7 curves plateaus.
    # Bottom right: P(s) s^2. A pure power law is a straight line and the break a minimum;
    #     above p_c a minimum opens at s* (circles) and the branch beyond rises as s^0.26
    #     at every p, i.e. tau_s^> = 1.74. Inset: s*(p) at three system sizes.
    ps = ps if ps is not None else available("emission", L, rho)
    fig = plt.figure(figsize=(19.5, 15.5))
    gs = fig.add_gridspec(2, 2, hspace=0.30, wspace=0.36)
    a, b, c, d = (fig.add_subplot(gs[i, j]) for i in (0, 1) for j in (0, 1))

    for p, cs, ce in zip(ps, seq(ps, CM_SPOT), seq(ps, CM_EMIS)):
        x, y = density(hist("spotSize", L, rho, p).sum(axis=0), minc=30, k=4)
        a.plot(x, y, "-", color=cs, lw=2.8)
        xc, yc = density(hist("spotSize", L, rho, p).sum(axis=0), minc=200, k=4)
        c.plot(xc, yc * xc ** 1.5, "-", color=cs, lw=2.8)
        h = hist("emission", L, rho, p).sum(axis=0)
        xe, ye = density(h, minc=30, k=4)
        b.plot(xe, ye, "-", color=ce, lw=2.8)
        xk, yk = density(h, minc=200, k=4)
        d.plot(xk, yk * xk ** 2, "-", color=ce, lw=2.8)
        s = branches(h)[2]
        if np.isfinite(s):
            d.plot(s, on_curve(xk, yk * xk ** 2, s), "o", color=ce, ms=16, mec="white",
                   mew=2.0, zorder=6)

    xs, ys = density(hist("spotSize", L, rho, 0.0).sum(axis=0), minc=30, k=4)
    guide_on(a, xs, ys, 2e1, 1.9, 2.0, off=0.035, label=r"$m^{-2}$", lab_at=0.75)
    xs, ys = density(hist("spotSize", L, rho, 1.0).sum(axis=0), minc=30, k=4)
    guide_on(a, xs, ys, 1e4, 2.2, 1.5, off=25, label=r"$m^{-3/2}$", lab_at=0.50, dy=1.8)
    a.set_xlabel("spot size $m$"); a.set_ylabel("$n(m)$")
    a.set_xscale("log"); a.set_yscale("log")
    a.set_xlim(0.7, 5e8); a.set_ylim(3e-15, 30)
    cbar(fig, a, ps, CM_SPOT, "$p$")

    xs, ys = density(hist("emission", L, rho, 0.0).sum(axis=0), minc=30, k=4)
    guide_on(b, xs, ys, 2e1, 1.5, 3.0, off=0.05, label=r"$s^{-3}$", lab_at=0.7)
    xs, ys = density(hist("emission", L, rho, 1.0).sum(axis=0), minc=30, k=4)
    guide_on(b, xs, ys, 1e3, 2.4, 1.74, off=25, label=r"$s^{-1.74}$", lab_at=0.50, dy=1.8)
    b.set_xlabel("flare size $s$"); b.set_ylabel("$P(s)$")
    b.set_xscale("log"); b.set_yscale("log")
    b.set_xlim(0.7, 5e6); b.set_ylim(3e-15, 30)
    cbar(fig, b, ps, CM_EMIS, "$p$")

    c.axhline(1.3, color=C_REF, ls="--", lw=2.6, zorder=2)
    c.text(3e5, 1.55, r"$\tau_m=3/2$", color=C_REF, fontsize=27)
    c.plot([6, 6e4], [0.11, 0.11 * 1e-2], "--", color=C_REF, lw=2.6, zorder=2)
    c.text(1.5, 2.5e-3, r"$\tau_m=2$", color=C_REF, fontsize=27)
    c.set_xlabel("spot size $m$"); c.set_ylabel(r"$m^{3/2}\,n(m)$")
    c.set_xscale("log"); c.set_yscale("log")
    c.set_xlim(0.7, 5e8); c.set_ylim(3e-4, 6)
    cbar(fig, c, ps, CM_SPOT, "$p$")

    d.plot([2e3, 2e5], [7.0, 7.0 * 100 ** 0.26], "--", color=C_REF, lw=2.8, zorder=2)
    d.text(1.2e4, 30, r"$s^{\,0.26}$", fontsize=27, color=C_REF, va="bottom")
    d.set_xlabel("flare size $s$"); d.set_ylabel(r"$s^{2}\,P(s)$")
    d.set_xscale("log"); d.set_yscale("log")
    d.set_xlim(0.7, 5e6); d.set_ylim(3e-8, 90)
    cbar(fig, d, ps, CM_EMIS, "$p$")

    ins = d.inset_axes([0.15, 0.10, 0.36, 0.29])
    for LL, mk, f in zip(Ls, ("o", "s", "^"), (0.25, 0.6, 0.95)):
        xx, yy = [], []
        for p in available("emission", LL, rho):
            s = branches(hist("emission", LL, rho, p).sum(axis=0))[2]
            if np.isfinite(s):
                xx.append(p); yy.append(s)
        ins.plot(xx, yy, mk + "-", color=CM_EMIS(f), ms=8, lw=2.0, mec="white", mew=1.0,
                 label=rf"${LL}$")
    ins.set_yscale("log")
    ins.set_xlabel("$p$", fontsize=21, labelpad=1)
    ins.set_ylabel("$s^*$", fontsize=21, labelpad=1)
    ins.set_xticks([0.7, 0.8, 0.9, 1.0])
    ins.tick_params(labelsize=17, top=False, right=False, direction="out", pad=2)
    ins.legend(frameon=False, fontsize=15, title="$L$", title_fontsize=15, loc="upper right",
               handlelength=1.0, labelspacing=0.12, borderpad=0.1, handletextpad=0.3)
    ins.set_facecolor("white")
    for ax in (a, b, c, d, ins):
        grid(ax)
    save(fig, f"fig2_spectra_L_{g(L)}_rho_{g(rho)}.png")


def fig3_mean_field(N=200000, L=512, rho=0.2):
    # The exponent is not universal, but it is not free either: the sub-leading balance of
    # the rate equation fixes it from the branching ratio alone,
    #     cos(pi tau_m) = (1 - q)/q,   tau_m in [3/2, 2],
    # with the two ends of the branch, (1/2, 2) and (1, 3/2), marked. The well-mixed Monte
    # Carlo shares every assumption of that derivation except the neglect of fluctuations
    # and sits on the curve to ~0.01; the lattice lies uniformly below, and that gap is
    # spatial correlation rather than a defect of the closed form.
    # Right: tau_s = 2 tau_m - 1, which follows from s = min(i, j) when both partners come
    # from the bulk: exact when well mixed, broken on the lattice in both directions.
    # Seed-to-seed errors are ~0.001, smaller than the markers, so none are drawn.
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(19, 8))

    qs = np.array(wm_qs(N))
    tm = np.array([tau_pm(wm_hist("wmSpotSize", N, q))[0] for q in qs])
    ts = np.array([tau_pm(wm_hist("wmEmission", N, q))[0] for q in qs])
    ok = qs >= 0.5
    lat = lattice_tau(L, rho)
    mix = lat[:, 6] == 0
    c_th = CM_EMIS(1.0)                     # the closed form: a result, so a colour of its own
    c_mc = tint(C_SPOT, 0.15)
    xs = dict(marker="x", ls="none", color=c_mc, ms=12, mew=3.0, zorder=4)
    dots = dict(marker="o", color=C_POS, ms=12, mec="white", mew=1.2)

    qc = np.linspace(0.5, 0.9995, 600)
    a1.plot(qc, tau_closed(qc), "-", color=c_th, lw=5, zorder=2,
            label=r"$\cos\pi\tau_m=\dfrac{1-q}{q}$")
    a1.plot(qs[ok], tm[ok], label="well-mixed MC", **xs)
    a1.plot(lat[:, 0], lat[:, 1], "-", lw=2.4, zorder=3, label=rf"2D lattice, $L={L}$", **dots)
    a1.plot([0.5, 1.0], [2.0, 1.5], "*", ms=28, color=c_th, mec="white", mew=1.5,
            zorder=5, ls="none")
    a1.set_xlabel("branching ratio $q$")
    a1.set_ylabel(r"spot exponent $\tau_m$")
    a1.set_xlim(0.495, 1.005); a1.set_ylim(1.42, 2.10)
    a1.legend(frameon=False, loc="upper right", labelspacing=0.35)
    grid(a1)

    lo, hi = 1.8, 4.2
    a2.plot([lo, hi], [lo, hi], "--", color=C_REF, lw=3.5, zorder=1, label=r"$\tau_s=2\tau_m-1$")
    a2.plot(2 * tm - 1, ts, label="well-mixed MC", **xs)
    a2.plot(2 * lat[mix, 1] - 1, lat[mix, 3], ls="none", zorder=4,
            label=r"2D lattice, $p<p_c$", **dots)
    a2.plot(2 * lat[~mix, 1] - 1, lat[~mix, 3], "o", ls="none", mfc="white", mec=C_POS,
            mew=2.2, ms=12, zorder=4, label=r"2D lattice, $p>p_c$")
    a2.set_xlabel(r"$2\tau_m-1$")
    a2.set_ylabel(r"flare exponent $\tau_s$")
    a2.set_xlim(lo, hi); a2.set_ylim(1.55, hi)
    a2.legend(frameon=False, loc="lower right", labelspacing=0.3, fontsize=20)
    grid(a2)

    fig.tight_layout()
    save(fig, f"fig3_meanField_N_{N:d}_L_{g(L)}.png")


def fig4_density_scaling(L=512, rhos=RHOS_D, p_lo=0.3, p_hi=0.9):
    # Density against the two distributions, on either side of the ordering transition,
    # over an eighteenfold range rho = 0.05 to 0.9 (run_density.sh, 8 seeds). Top row: a
    # mixed surface; bottom row: a segregated one. Left column spots, right column flares.
    # Everything lies on one curve; what density moves is the upper cutoff, which follows
    # the total mass rho L^2, and above p_c the position of the break s* (circles), which
    # drops about tenfold across the range. The slopes drift by at most ~0.1 when mixed
    # and ~0.05 when segregated.
    fig, axes = plt.subplots(2, 2, figsize=(18.5, 14.5))
    for i, pp in enumerate((p_lo, p_hi)):
        for j, (kind, cm, xl, yl) in enumerate(
                (("spotSize", CM_SPOT, "spot size $m$", "$n(m)$"),
                 ("emission", CM_EMIS, "flare size $s$", "$P(s)$"))):
            ax = axes[i][j]
            for r, col in zip(rhos, seq(rhos, cm)):
                h = hist(kind, L, r, pp).sum(axis=0)
                x, y = density(h, minc=30, k=4)
                ax.plot(x, y, "-", color=col, lw=3.0, label=rf"$\rho={g(r)}$")
                if kind == "emission":
                    s = branches(h)[2]
                    if np.isfinite(s):
                        ax.plot(s, on_curve(x, y, s), "o", color=col, ms=16, mec="white",
                                mew=2.0, zorder=6)
            ax.set_xscale("log"); ax.set_yscale("log")
            ax.set_xlabel(xl); ax.set_ylabel(yl)
            ax.set_xlim(0.7, 5e7); ax.set_ylim(1e-14, 60)   # common, so rows compare
            ax.set_title(rf"$p={g(pp)}$", fontsize=34, pad=12)
            # reference slope: the median exponent over the densities, below the curves
            Hs = [hist(kind, L, r, pp) for r in rhos]
            if kind == "spotSize":
                t = np.median([tau_pm(H)[0] for H in Hs])
                sym = "m"
            elif np.isfinite(branches(Hs[0].sum(axis=0))[2]):
                t = np.nanmedian([upper_tau(H)[0] for H in Hs])
                sym = "s"
            else:
                t = np.median([tau_pm(H, lo=20.0)[0] for H in Hs])
                sym = "s"
            x0, y0 = density(Hs[len(Hs) // 2].sum(axis=0), minc=30, k=4)
            xg = 3e2 if kind == "spotSize" or t > 2 else 1e3
            guide_on(ax, x0, y0, xg, 2.4 if t < 2 else 1.6, t, off=0.03)
            ax.text(xg / 1.6, on_curve(x0, y0, xg) * 0.03, rf"${sym}^{{-{t:.2f}}}$",
                    color=C_REF, fontsize=27, ha="right", va="center")
            if i == 0:
                ax.legend(frameon=False, loc="upper right", labelspacing=0.2,
                          handlelength=1.3, fontsize=21)
            grid(ax)
    fig.tight_layout()
    save(fig, f"fig4_densityScaling_L_{g(L)}_p_{g(p_lo)}_{g(p_hi)}.png")


def fig5_branching_ratio(L=128, rhos=RHOS_SW, Ls=(128, 256, 512, 1024), rho=0.2):
    # q is the only thing the closed form needs, so it is worth seeing on its own. It
    # leaves 1/2 -- the exactly uncorrelated value, which p = 0 reproduces to four digits
    # -- and climbs toward 1 as the signs segregate. Left: a denser surface reaches a
    # given q at slightly larger p, which is the whole of the density dependence in the
    # exponents. Right: at fixed density the curve is the same from L = 128 to 1024.
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(19, 8))
    for r, col in zip(rhos, seq(rhos, CM_SPOT, lo=0.1)):
        ps = available("snapshots", L, r)
        a1.plot(ps, [branching_ratio(L, r, p) for p in ps], "o-", color=col, ms=10,
                lw=2.6, mec="white", mew=1.2, label=rf"$\rho={g(r)}$")
    for LL, mk, col in zip(Ls, ("o", "s", "^", "D"), seq(Ls, CM_SPOT, lo=0.1)):
        ps = available("snapshots", LL, rho)
        a2.plot(ps, [branching_ratio(LL, rho, p) for p in ps], mk + "-", color=col,
                ms=10, lw=2.6, mec="white", mew=1.2, label=rf"$L={LL}$")
    for a, ttl in ((a1, rf"$L={L}$"), (a2, rf"$\rho={g(rho)}$")):
        a.axhline(0.5, color=C_NEG, ls=":", lw=2.5)
        a.axhline(1.0, color=C_POS, ls=":", lw=2.5)
        a.text(0.97, 0.515, "mixed", color=C_NEG, fontsize=23, va="bottom", ha="right")
        a.text(0.97, 0.985, "segregated", color=C_POS, fontsize=23, va="top", ha="right")
        a.set_xlabel("neighbour probability $p$")
        a.set_ylabel("branching ratio $q$")
        a.set_xlim(-0.03, 1.03); a.set_ylim(0.46, 1.04)
        a.set_title(ttl, fontsize=28, pad=12)
        a.legend(frameon=False, loc="upper left", bbox_to_anchor=(0.02, 0.9),
                 labelspacing=0.25)
        grid(a)
    fig.tight_layout()
    save(fig, f"fig5_branchingRatio_L_{g(L)}.png")


if __name__ == "__main__":
    fig1_snapshots()
    fig2_spectra()
    fig3_mean_field()
    fig4_density_scaling()
    fig5_branching_ratio()
