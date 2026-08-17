# Collision kernel map: which pairs of masses actually meet.
# Data from collisionPairs.cpp (writes outputs/pairs_*.tsv with a C/A kind column).
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

DIR = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(DIR, "outputs")
PLOTS = os.path.join(DIR, "plots", "collisionKernel")

plt.rcParams.update({"font.size": 17, "axes.labelsize": 21, "xtick.labelsize": 15,
                     "ytick.labelsize": 15, "legend.fontsize": 14})


def load_pairs(L, rho, p, seed=1):
    f = os.path.join(OUT, f"pairs_L_{L:g}_rho_{rho:g}_p_{p:g}_seed_{seed}.tsv")
    m1, m2, kind = [], [], []
    with open(f) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            a, b, k = line.split()
            m1.append(int(a)); m2.append(int(b)); kind.append(k)
    return np.array(m1), np.array(m2), np.array(kind)


def collisionKernel_map(L=128, rho=0.2, ps=(0.5, 1.0), nb=70):
    for p in ps:
        m1, m2, kind = load_pairs(L, rho, p)
        fig, axes = plt.subplots(1, 2, figsize=(13.5, 6))
        for ax, k, name in ((axes[0], "C", "coagulation"), (axes[1], "A", "annihilation")):
            sel = kind == k
            a, b = m1[sel], m2[sel]
            # the kernel is symmetric in its two arguments, so mirror the pairs
            x = np.concatenate([a, b]); y = np.concatenate([b, a])
            ok = (x > 0) & (y > 0)
            # masses are integers, so log-spaced edges must land on integer boundaries or
            # the small-m bins alternate between one integer and none, and the map stripes
            edges = np.unique(np.round(np.geomspace(1, max(x.max(), y.max()) + 1, nb)))
            H, _, _ = np.histogram2d(x[ok], y[ok], bins=[edges, edges])
            im = ax.pcolormesh(edges, edges, np.ma.masked_equal(H.T, 0),
                               norm=LogNorm(vmin=1), cmap="magma")
            ax.set_facecolor("#000000")          # empty cells read as background, not holes
            ax.plot(edges, edges, "-", color="white", lw=1, alpha=0.35)
            ax.set_xscale("log"); ax.set_yscale("log"); ax.set_aspect("equal")
            ax.set_xlabel("$m_1$"); ax.set_title(f"{name}  ({sel.sum()} events)")
            fig.colorbar(im, ax=ax, label="collisions")
        axes[0].set_ylabel("$m_2$")
        fig.suptitle(f"$p={p:g}$", fontsize=20)
        fig.tight_layout()
        os.makedirs(PLOTS, exist_ok=True)
        fig.savefig(os.path.join(PLOTS, f"collisionKernel_L_{L:g}_rho_{rho:g}_p_{p:g}.png"), dpi=300)
        plt.close(fig)


if __name__ == "__main__":
    collisionKernel_map()
