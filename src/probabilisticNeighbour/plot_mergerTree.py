# Merger tree of the largest final spot: who ate whom, and when.
# Time runs downward. The x axis is layout only -- leaves are placed left to right in the
# order the recursion reaches them and each parent sits above its children; it carries no
# physical meaning. Every spot is born with mass 1, so each branch is drawn as a chain of
# segments whose width and colour follow its mass as it actually grew.
# Data from lineage.cpp (outputs/lineage_*.tsv).
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import LogNorm

DIR = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(DIR, "outputs")
PLOTS = os.path.join(DIR, "plots", "lineage")

plt.rcParams.update({"font.size": 17, "axes.labelsize": 21, "xtick.labelsize": 15,
                     "ytick.labelsize": 15})


def load_lineage(L, rho, p, seed):
    stem = f"L_{L:g}_rho_{rho:g}_p_{p:g}_seed_{seed}"
    births, events = {}, []
    with open(os.path.join(OUT, f"lineage_{stem}.tsv")) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split()
            if f[0] == "B":
                births[int(f[2])] = int(f[1])
            else:
                events.append((f[0], int(f[1]), int(f[2]), int(f[3]), int(f[4]), int(f[5])))
    final = {}
    with open(os.path.join(OUT, f"lineageFinal_{stem}.tsv")) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split()
            final[int(f[0])] = abs(int(f[3]))
    return births, events, final


def build_children(events):
    # the heavier spot survives and keeps its identity; record what it met and whether the
    # encounter added mass (coagulation) or removed it (annihilation)
    kids = {}
    for kind, step, tid, mid, mt, mm in events:
        if kind == "A" and mt == mm:
            continue                                  # mutual annihilation: no survivor
        win, lose, mlose = ((mid, tid, mt) if mm > mt else (tid, mid, mm))
        kids.setdefault(win, []).append((lose, step, mlose, kind))
    return kids


def mass_track(node, kids, births, tend, tol=0.05):
    # reconstruct mass over the branch's life: born at 1, gaining on C and losing on A.
    # Points are thinned to changes of more than `tol`, so the collection stays small.
    ev = sorted(kids.get(node, []), key=lambda e: e[1])
    t, m = [births.get(node, 0)], [1.0]
    cur = 1.0
    for _, step, ml, kind in ev:
        if step > tend:
            break
        cur = max(cur + (ml if kind == "C" else -ml), 1.0)
        if abs(cur - m[-1]) / m[-1] > tol:
            t.append(step); m.append(cur)
    t.append(tend); m.append(cur)
    return np.array(t), np.array(m)


def plot_tree(L=64, rho=0.2, p=1.0, seed=1, mmin=40):
    sys.setrecursionlimit(100000)
    births, events, final = load_lineage(L, rho, p, seed)
    kids = build_children(events)
    root = max(final, key=lambda k: final[k])
    tmax = max(births.values())

    # keep only branches that brought in at least mmin: the rest is dust being accreted
    def prune(n):
        out = []
        for c, t, m, kind in kids.get(n, []):
            if m >= mmin:
                out.append((c, t, m, prune(c)))
        out.sort(key=lambda e: e[1])
        return out

    tree = prune(root)
    segs, masses = [], []
    xpos, counter = {}, [0.0]

    def layout(node, sub, tmerge, mass):
        if sub:
            for c, t, m, s in sub:
                layout(c, s, t, m)
            x = float(np.mean([xpos[c] for c, _, _, _ in sub]))
        else:
            counter[0] += 1.0
            x = counter[0]
        xpos[node] = x
        tt, mm = mass_track(node, kids, births, tmerge)
        for i in range(len(tt) - 1):                    # the branch, thickening as it grows
            segs.append([(x, tt[i]), (x, tt[i + 1])])
            masses.append(mm[i])
        for c, t, m, s in sub:                          # the moment it is swallowed
            segs.append([(xpos[c], t), (x, t)])
            masses.append(m)
        return x

    layout(root, tree, tmax, final[root])

    masses = np.array(masses)
    order = np.argsort(masses)                          # heavy branches drawn on top
    segs = [segs[i] for i in order]
    masses = masses[order]
    lw = 0.3 + 3.4 * (np.log10(masses + 1) / np.log10(masses.max() + 1)) ** 2.5
    lc = LineCollection(segs, linewidths=lw, array=masses, cmap="YlGnBu",
                        norm=LogNorm(vmin=1, vmax=masses.max()), capstyle="round")

    fig, ax = plt.subplots(figsize=(11, 7.5))
    ax.add_collection(lc)
    ax.set_xlim(0, counter[0] + 1)
    ax.set_ylim(tmax, 0)
    ax.set_xticks([])
    ax.set_ylabel("time (sweeps)")
    ax.set_title(f"$p={p:g}$: ancestry of the largest spot "
                 f"(final mass {final[root]:,}; branches reaching $\\geq$ {mmin})",
                 fontsize=15)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    cb = fig.colorbar(lc, ax=ax, pad=0.02)
    cb.set_label("mass at that moment")
    fig.tight_layout()
    os.makedirs(PLOTS, exist_ok=True)
    fig.savefig(os.path.join(PLOTS, f"mergerTree_L_{L:g}_rho_{rho:g}_p_{p:g}.png"), dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    plot_tree()
