# Drainage basins of the cascade: colour every spot's birth site by which final spot its
# lineage ended up in. Data from lineage.cpp (outputs/lineage_*.tsv).
#
# "Basin" means ancestry, not conserved mass: annihilation destroys mass on the way, so a
# birth site belongs to the final spot that swallowed its line, not to a surviving gram.
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

DIR = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(DIR, "outputs")
PLOTS = os.path.join(DIR, "plots", "lineage")

# equal-width black and white rings. Round joins matter here: with the default miter the
# stroke overshoots at the star's sharp points and leaves the tips ragged.
STAR_PE = [pe.Stroke(linewidth=1.8, foreground="w", joinstyle="round", capstyle="round"),
           pe.Stroke(linewidth=0.9, foreground="k", joinstyle="round", capstyle="round"),
           pe.Normal()]

plt.rcParams.update({"font.size": 17, "axes.labelsize": 21, "xtick.labelsize": 15,
                     "ytick.labelsize": 15, "legend.fontsize": 14})


def load_lineage(L, rho, p, seed):
    stem = f"L_{L:g}_rho_{rho:g}_p_{p:g}_seed_{seed}"
    births, events = {}, []
    with open(os.path.join(OUT, f"lineage_{stem}.tsv")) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split()
            if f[0] == "B":
                births[int(f[2])] = (int(f[3]), int(f[4]), int(f[1]))
            else:
                # kind, step, survivorCandidate(target), mover, mTarget, mMover
                events.append((f[0], int(f[1]), int(f[2]), int(f[3]), int(f[4]), int(f[5])))
    final = {}
    with open(os.path.join(OUT, f"lineageFinal_{stem}.tsv")) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split()
            final[int(f[0])] = (int(f[1]), int(f[2]), int(f[3]))   # signed mass
    return births, events, final


def resolve(events):
    # who survived each collision: the heavier spot keeps its identity (ties go to the
    # target, matching the C++). Both die when opposite signs have equal mass.
    parent = {}
    for kind, step, tid, mid, mt, mm in events:
        if kind == "A" and mt == mm:
            parent[tid] = None; parent[mid] = None      # mutual annihilation ends both
            continue
        win, lose = (mid, tid) if mm > mt else (tid, mid)
        parent[lose] = win
    return parent


def root_of(i, parent, memo):
    path = []
    while True:
        if i in memo:
            r = memo[i]; break
        if i not in parent:
            r = i; break
        if parent[i] is None:
            r = None; break
        path.append(i); i = parent[i]
    for q in path:
        memo[q] = r
    return r


def plot_basins(L=64, rho=0.2, p=1.0, seed=1, nbasin=8, frac=0.25):
    # Births are aggregated per site, not scattered: over a long run every site is born on
    # many times. Only the last `frac` of the run is used -- a spot diffuses ~sqrt(t), so
    # ancestry from early times has wandered too far to show any catchment structure.
    births, events, final = load_lineage(L, rho, p, seed)
    parent = resolve(events)
    memo = {}
    top = sorted(final, key=lambda k: -abs(final[k][2]))[:nbasin]
    rank = {k: i for i, k in enumerate(top)}
    tmax = max(b[2] for b in births.values())
    t0 = tmax * (1 - frac)

    counts = np.zeros((L, L, nbasin))
    used = inwin = 0
    for bid, (x, y, t) in births.items():
        if t < t0:
            continue
        inwin += 1
        r = root_of(bid, parent, memo)
        if r in rank:
            counts[y, x, rank[r]] += 1; used += 1

    tot = counts.sum(axis=2)
    win = counts.argmax(axis=2)
    dom = np.divide(counts.max(axis=2), np.maximum(tot, 1))
    cmap = plt.cm.turbo(np.linspace(0.05, 0.95, nbasin))
    rgb = cmap[win][:, :, :3]
    # brightness = how strongly one basin dominates that site, rescaled off the 1/nbasin floor
    a = np.clip((dom - 1.0 / nbasin) / (1 - 1.0 / nbasin), 0, 1)[:, :, None]
    img = rgb * a ** 0.6
    # sites with no ancestry drawn are grey, so "no data" cannot be read as "contested"
    nodata = tot == 0
    img[nodata] = 0.0

    fig, ax = plt.subplots(figsize=(9.2, 7))
    ax.imshow(img, origin="lower", interpolation="nearest")
    handles, labels = [], []
    for i, k in enumerate(top):
        fx, fy, fm = final[k]
        ax.plot(fx, fy, "*", ms=19, color=cmap[i], mec="none", zorder=5,
                path_effects=STAR_PE)
        handles.append(plt.Line2D([], [], ls="", marker="*", ms=15, color=cmap[i],
                                  mec="none", path_effects=STAR_PE))
        labels.append(f"{fm:+,}")
    handles.append(plt.Line2D([], [], ls="", marker="s", ms=11, color="k"))
    labels.append("none of these 8")
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f"$p={p:g}$: which final spot each site's lineage feeds\n"
                 f"(last {frac:.0%} of the run; {used:,} of {inwin:,} births shown)",
                 fontsize=15)
    # a site is born on many times and those births need not share a destination, so the
    # hue is the most common one and the brightness is how large that majority was
    leg = ax.legend(handles, labels, frameon=False, loc="center left",
                    bbox_to_anchor=(1.01, 0.5), fontsize=13, handletextpad=0.4,
                    labelspacing=0.6, title="final spot mass\n\nhue: most common\n"
                                            "destination of a\nsite's births\n"
                                            "brightness: its\nshare of them")
    leg.get_title().set_fontsize(12)
    leg.get_title().set_multialignment("left")
    fig.tight_layout()
    os.makedirs(PLOTS, exist_ok=True)
    fig.savefig(os.path.join(PLOTS, f"basins_L_{L:g}_rho_{rho:g}_p_{p:g}.png"),
                dpi=300, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    plot_basins()
