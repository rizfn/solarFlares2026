# Numerical steady state of the mean-field rate equation of renormalization.md, and its
# comparison with the simulation.
#
#   dn_m/dt = (q/2) sum_{i+j=m} n_i n_j + (1-q) sum_{j>=1} n_{m+j} n_j - n_m N + J d_{m,1}
#
# Forward Euler with FFT convolutions; mass coagulating past the cutoff M leaves the
# system, which is what keeps the cascade steady. The only input is the branching ratio
# q, which the simulation reports as the same-sign fraction of neighbouring spots.

import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import plots as P

CACHE = os.path.join(P.OUT, "meanField.npy")


def mf_steady(q, M=2048, J=1.0, dt=0.05, steps=250000, tol=1e-10):
    n = np.zeros(M + 1)
    n[1] = 1.0
    nfft = 2 * (M + 1)
    for it in range(steps):
        N = n.sum()
        fn = np.fft.rfft(n, nfft)
        conv = np.fft.irfft(fn * fn, nfft)             # sum_{i+j=m} n_i n_j
        acf = np.fft.irfft(np.abs(fn) ** 2, nfft)      # sum_j n_j n_{j+m}
        dn = 0.5 * q * conv[:M + 1] + (1 - q) * acf[:M + 1] - n * N
        dn[1] += J
        nn = np.maximum(n + dt * dn, 0.0)
        nn[0] = 0.0
        if it % 5000 == 0 and it > 0 and np.abs(nn - n).sum() / (N + 1e-12) < tol:
            return nn
        n = nn
    return n


def tau_closed(q):
    # Sub-leading balance of the rate equation (renormalization.md S5):
    # cos(pi tau) = (1-q)/q, on the branch tau in [3/2, 2]. No solution for q < 1/2.
    q = np.asarray(q, dtype=float)
    r = (1 - q) / q
    return np.where(r <= 1, 2 - np.arccos(np.clip(r, -1, 1)) / np.pi, np.nan)


def fit_tau(n, lo=10, hi=256):
    m = np.arange(len(n))
    sel = (m >= lo) & (m <= hi) & (n > 0)
    x, y = np.log(m[sel]), np.log(n[sel])
    return -np.polyfit(x, y, 1)[0]


def wm_tau(q, N=200000, kind="wmSpotSize"):
    # Hill fit to the well-mixed Monte Carlo (wellMixed.cpp), summed over seeds
    fs = sorted(glob.glob(os.path.join(P.OUT, f"{kind}_N_{N}_q_{P.g(q)}_seed_*.tsv")))
    if not fs:
        raise OSError(f"no well-mixed output for q={q}")
    total = {}
    for f in fs:
        for s, c in np.loadtxt(f, dtype=np.int64, ndmin=2):
            total[s] = total.get(s, 0) + c
    sizes = np.array(sorted(total), dtype=float)
    return P.mle(sizes, np.array([total[s] for s in sizes], dtype=float))


def measure(L=128, rho=0.2, ps=(0.0, 0.2, 0.4, 0.6, 0.8, 0.9, 1.0), **kw):
    # (p, q, tau_sim, tau_mf) for each p; cached, the solver takes ~a minute per q
    if os.path.exists(CACHE):
        return np.load(CACHE)
    rows = []
    for p in ps:
        try:
            q = P.same_sign_frac(P.load_snaps(L, rho, p))
            ts = P.mle(*P.load_hist("spotSize", L, rho, p))
        except OSError:
            continue
        tmf = fit_tau(mf_steady(q, **kw))
        rows.append((p, q, ts, tmf))
        print(f"p={p:.2f}  q={q:.3f}  tau_sim={ts:.3f}  tau_MF={tmf:.3f}", flush=True)
    rows = np.array(rows)
    np.save(CACHE, rows)
    return rows


def plot_mean_field(L=128, rho=0.2, **kw):
    p, q, ts, tmf = measure(L, rho, **kw).T
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(14, 6))

    a1.plot(p, q, "o-", color="#2166ac", ms=8, lw=2.5)
    a1.axhline(0.5, color="grey", ls=":", lw=1.5)
    a1.set_xlabel("neighbour probability $p$")
    a1.set_ylabel("same-sign fraction $q$")

    qc = np.linspace(0.5001, 0.999, 300)
    a2.plot(qc, tau_closed(qc), "-", color="#1a9850", lw=3, alpha=0.8,
            label=r"$\cos\pi\tau=\frac{1-q}{q}$")
    a2.plot(q, ts, "o-", color="#b2182b", ms=8, lw=2.5, label="simulation")
    a2.plot(q, tmf, "s--", color="#4d4d4d", ms=8, lw=2.5, label="mean field")
    a2.axhline(1.5, color="grey", ls=":", lw=1.5)
    a2.axhline(2.0, color="grey", ls=":", lw=1.5)
    a2.set_xlabel("same-sign fraction $q$")
    a2.set_ylabel(r"spot exponent $\tau_m$")
    # the mean field runs away once annihilation dominates; clip so the agreement above
    # q ~ 1/2 stays readable
    a2.set_ylim(1.4, 2.6)
    a2.text(q[0], 2.52, rf"$\tau_{{\rm MF}}={tmf[0]:.1f}$ at $q={q[0]:.2f}$",
            fontsize=15, color="#4d4d4d")
    a2.legend(frameon=False, loc="upper right")

    for a in (a1, a2):
        P.grid(a)
    fig.tight_layout()
    fig.savefig(os.path.join(P.PLOTS, "meanField", f"meanField_L_{P.g(L)}_rho_{P.g(rho)}.png"), dpi=300)
    plt.close(fig)


def plot_well_mixed(L=128, rho=0.2, N=200000,
                    qs=(0.386, 0.45, 0.50, 0.55, 0.619, 0.70, 0.760, 0.849, 0.952, 0.99)):
    # three routes to tau(q): the closed form, the stochastic well-mixed process, and
    # the 2D lattice. The first two share every assumption except fluctuations; the
    # third additionally has space, so the gap between them is spatial correlation.
    p2, q2, tm2, tmf2 = measure(L, rho).T
    ts2 = np.array([P.mle(*P.load_hist("emission", L, rho, p)) for p in p2])

    wq, wm, we = [], [], []
    for q in qs:
        try:
            m, e = wm_tau(q, N), wm_tau(q, N, "wmEmission")
        except OSError:
            continue
        wq.append(q); wm.append(m); we.append(e)
    wq, wm, we = map(np.array, (wq, wm, we))

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(14, 6))

    qc = np.linspace(0.5001, 0.999, 300)
    a1.plot(qc, tau_closed(qc), "-", color="#1a9850", lw=3.5, alpha=0.85,
            label=r"$\cos\pi\tau=\frac{1-q}{q}$")
    a1.plot(wq, wm, "o", color="black", ms=9, label="well-mixed MC")
    a1.plot(q2, tmf2, "s--", color="#4d4d4d", ms=7, lw=1.5, label="rate equation")
    a1.plot(q2, tm2, "^-", color="#b2182b", ms=9, lw=2.5, label="2D lattice")
    a1.axvline(0.5, color="grey", ls=":", lw=1.5)
    a1.set_xlabel("same-sign fraction $q$")
    a1.set_ylabel(r"spot exponent $\tau_m$")
    a1.set_ylim(1.4, 3.6)
    a1.legend(frameon=False)

    # tau_s = 2 tau_m - 1 assumes both partners come from the bulk: exact when mixed,
    # broken on the lattice where emissions only happen on domain interfaces
    lo, hi = 1.8, 4.0
    a2.plot([lo, hi], [lo, hi], "-", color="grey", lw=1.5, label=r"$\tau_s=2\tau_m-1$")
    a2.plot(2 * wm - 1, we, "o", color="black", ms=9, label="well-mixed MC")
    a2.plot(2 * tm2 - 1, ts2, "^", color="#b2182b", ms=9, label="2D lattice")
    a2.set_xlabel(r"$2\tau_m-1$")
    a2.set_ylabel(r"emission exponent $\tau_s$")
    a2.set_xlim(lo, hi); a2.set_ylim(lo, hi)
    a2.legend(frameon=False)

    for a in (a1, a2):
        P.grid(a)
    fig.tight_layout()
    fig.savefig(os.path.join(P.PLOTS, "meanField", f"wellMixed_N_{N}.png"), dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    plot_mean_field()
    plot_well_mixed()
