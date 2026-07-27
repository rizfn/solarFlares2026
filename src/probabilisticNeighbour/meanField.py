# Numerical steady state of the mean-field rate equation of renormalization.md, and its
# comparison with the simulation.
#
#   dn_m/dt = (q/2) sum_{i+j=m} n_i n_j + (1-q) sum_{j>=1} n_{m+j} n_j - n_m N + J d_{m,1}
#
# Forward Euler with FFT convolutions; mass coagulating past the cutoff M leaves the
# system, which is what keeps the cascade steady. The only input is the branching ratio
# q, which the simulation reports as the same-sign fraction of neighbouring spots.

import os
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


def fit_tau(n, lo=10, hi=256):
    m = np.arange(len(n))
    sel = (m >= lo) & (m <= hi) & (n > 0)
    x, y = np.log(m[sel]), np.log(n[sel])
    return -np.polyfit(x, y, 1)[0]


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


if __name__ == "__main__":
    plot_mean_field()
