# Showcase figures for talks, posters and pitches, into plots/paper/options/. White
# backgrounds and the paper's palette (paperFigures.py): red and blue are the two signs,
# C_SPOT and C_EMIS carry spot and flare quantities. The model itself is front and centre.
#
#   banner        order from chaos: the surface stitched across p, Phi riding on top
#   terrain       the signed mass field as a 3D landscape, at three values of p
#   magnetogram   the same field seen from above, at three values of p
#   hero          one large segregated magnetogram, nothing else on it

import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import LightSource
from scipy.ndimage import gaussian_filter, zoom
import paperFigures as F

OPT = os.path.join(F.PLOTS, "options")

# signed mass: heavy negative spots deep blue, light ones pale, empty white, then the
# same in red -- the snapshot colours of the paper, with mass as saturation
CM_SIGN = mcolors.LinearSegmentedColormap.from_list("", [
    F.tint(F.C_NEG, 0.45), F.C_NEG, F.tint(F.C_NEG, -0.55), F.tint(F.C_NEG, -0.9), "white",
    F.tint(F.C_POS, -0.9), F.tint(F.C_POS, -0.55), F.C_POS, F.tint(F.C_POS, 0.45)], N=1024)


def save(fig, name):
    os.makedirs(OPT, exist_ok=True)
    fig.savefig(os.path.join(OPT, name), dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("wrote plots/paper/options/" + name, flush=True)


def read_mass(L, rho, p):
    fs = glob.glob(os.path.join(F.OUT, "massSnapHero", f"massSnap_L_{L}_rho_{F.g(rho)}_p_{F.g(p)}_seed_*.tsv"))
    if not fs:
        raise OSError(f"no mass snapshot for L={L} rho={rho} p={p}")
    rows = []
    with open(fs[0]) as fh:
        for line in fh:
            v = line.split()
            if len(v) == L:
                rows.append(np.fromiter(map(int, v), dtype=np.int64, count=L))
    n = len(rows) // L
    return np.array(rows[:n * L]).reshape(n, L, L)[-1]


def fields(L, rho, ps, sigma=1.2, q=99.7):
    # Signed log-mass, lightly smoothed so heavy spots read as patches rather than single
    # pixels, all scaled by one high percentile taken over every panel: a common scale,
    # so the mixed surface really does look as light as it is next to the segregated one.
    fs = [gaussian_filter(np.sign(m) * np.log10(1 + np.abs(m)), sigma)
          for m in (read_mass(L, rho, p) for p in ps)]
    s = np.percentile(np.abs(np.concatenate([f.ravel() for f in fs])), q)
    return [np.clip(f / s, -1, 1) for f in fs]


def stretch(f, gamma=0.4):
    # colour only: |f|^gamma lifts light spots so their sign shows, while the heaviest
    # stay the deepest. Heights and the common scale are untouched.
    return np.sign(f) * np.abs(f) ** gamma


def opt_banner(L=128, rho=0.6, width=14, rows=3):
    # Order from chaos in one strip. Each vertical slice is a real snapshot at its own p,
    # from p = 0 on the left to p = 1 on the right, `rows` independent runs stacked. The
    # order parameter Phi rides along on top.
    ps = F.available("snapshots", L, rho)
    cols = []
    for p in ps:
        sn = F.read_snaps(L, rho, p, nseed=rows)
        picks = [sn[k * (len(sn) // rows)] for k in range(rows)]
        c0 = L // 2 - width // 2
        cols.append(np.vstack([s[:, c0:c0 + width] for s in picks]))
    img = np.hstack(cols).astype(float)
    sign = mcolors.ListedColormap([F.C_NEG, "white", F.C_POS])
    rgb = zoom(sign(mcolors.Normalize(-1, 1)(img))[..., :3], (3, 3, 1), order=0)
    fig = plt.figure(figsize=(26, 8.2))
    ax = fig.add_axes([0.04, 0.16, 0.92, 0.78])
    ax.imshow(rgb, extent=(0, 1, 0, 1), aspect="auto", interpolation="nearest")
    ph = np.array([F.order_parameter(L, rho, p) for p in ps])
    xs = np.linspace(0.5 / len(ps), 1 - 0.5 / len(ps), len(ps))
    ax.plot(xs, ph, color="white", lw=9, solid_capstyle="round", zorder=3)
    ax.plot(xs, ph, color=F.CM_EMIS(1.0), lw=4.2, solid_capstyle="round", zorder=4)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_yticks([0, 0.5, 1]); ax.set_ylabel(r"$\Phi$", fontsize=34)
    ax.set_xticks(np.linspace(0, 1, 6))
    ax.set_xticklabels([F.g(round(v, 1)) for v in np.linspace(0, 1, 6)])
    ax.set_xlabel("neighbour probability $p$", fontsize=30)
    ax.tick_params(top=False, right=False)
    for lab, x, ha in (("mixed", 0.012, "left"), ("segregated", 0.988, "right")):
        ax.text(x, 0.93, lab, fontsize=30, transform=ax.transAxes, va="top", ha=ha,
                bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none", alpha=0.85))
    save(fig, f"banner_L_{L}_rho_{F.g(rho)}.png")


def surface(ax, Z, elev=34, azim=-58):
    rgb = LightSource(azdeg=315, altdeg=38).shade_rgb(
        CM_SIGN(mcolors.Normalize(-1, 1)(stretch(Z)))[..., :3], Z, blend_mode="soft",
        vert_exag=3)
    X, Y = np.meshgrid(np.arange(Z.shape[1]), np.arange(Z.shape[0]))
    ax.plot_surface(X, Y, Z, facecolors=rgb, rstride=1, cstride=1, linewidth=0,
                    antialiased=False, shade=False)
    ax.view_init(elev=elev, azim=azim)
    ax.set_zlim(-1, 1)
    ax.set_box_aspect((1, 1, 0.38), zoom=1.15)
    ax.set_axis_off()


def crop_on_seam(f, crop):
    # centre the window on the column where the row-averaged sign flips, if there is one
    L = f.shape[0]
    prof = gaussian_filter(np.sign(f).mean(axis=0), 8)
    inner = prof[crop // 2:L - crop // 2]
    cx = int(np.argmin(np.abs(inner))) + crop // 2 if np.abs(inner).max() > 0.3 else L // 2
    return f[L // 2 - crop // 2:L // 2 + crop // 2, cx - crop // 2:cx + crop // 2]


def opt_terrain(L=512, rho=0.6, ps=(0.0, 0.5, 1.0), crop=200, smooth=1.6):
    # The mass field as a landscape at three values of p, on one common height scale:
    # positive spots are peaks, negative ones valleys (colour gamma-stretched, so light
    # spots still show their sign). Mixed, it is a low, choppy sea of
    # small spots of both signs; segregated, a red plateau of heavy spots faces a blue
    # basin across a cliff -- and the cliff is where the flares go off.
    fs = fields(L, rho, ps, sigma=smooth)
    fig = plt.figure(figsize=(10 * len(ps), 9))
    for k, (f, p) in enumerate(zip(fs, ps)):
        ax = fig.add_axes([k / len(ps) - 0.04, -0.12, 1 / len(ps) + 0.08, 1.12],
                          projection="3d", computed_zorder=False)
        surface(ax, crop_on_seam(f, crop))
        fig.text((k + 0.5) / len(ps), 0.9, rf"$p={F.g(p)}$", ha="center", fontsize=40)
    save(fig, f"terrain_L_{L}_rho_{F.g(rho)}.png")


def opt_terrain_single(L=512, rho=0.6, p=1.0, crop=200, smooth=1.6):
    f = fields(L, rho, (p,), sigma=smooth)[0]
    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_axes([-0.1, -0.25, 1.2, 1.5], projection="3d", computed_zorder=False)
    surface(ax, crop_on_seam(f, crop))
    save(fig, f"terrainSingle_L_{L}_rho_{F.g(rho)}_p_{F.g(p)}.png")


def opt_magnetogram(L=512, rho=0.6, ps=(0.0, 0.5, 1.0)):
    # Every pixel a spot, colour its sign, saturation its mass on a log scale. Mixed, the
    # surface is a fine speckle of light spots; segregated, two domains of heavy spots
    # face each other across a seam.
    fs = fields(L, rho, ps)
    fig, axes = plt.subplots(1, len(ps), figsize=(8 * len(ps), 8.6))
    for ax, f, p in zip(axes, fs, ps):
        im = ax.imshow(stretch(f), cmap=CM_SIGN, vmin=-1, vmax=1, interpolation="nearest")
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_color(F.C_REF); sp.set_linewidth(1.2)
        ax.set_title(rf"$p={F.g(p)}$", fontsize=34, pad=14)
    cb = fig.colorbar(im, ax=axes, orientation="horizontal", fraction=0.05, pad=0.03, aspect=60,
                      ticks=[-1, 0, 1])
    cb.set_ticklabels(["heavy $-$", "empty", "heavy $+$"])
    cb.ax.tick_params(labelsize=24)
    cb.set_label("signed spot mass (log scale, colour stretched)", fontsize=28)
    save(fig, f"magnetogram_L_{L}_rho_{F.g(rho)}.png")




def opt_hero(L=512, rho=0.6, p=1.0):
    # one big panel, nothing else on it: the poster image
    f = fields(L, rho, (p,), sigma=1.6)[0]
    fig = plt.figure(figsize=(12, 12))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.imshow(stretch(f), cmap=CM_SIGN, vmin=-1, vmax=1, interpolation="bicubic")
    ax.set_axis_off()
    save(fig, f"hero_L_{L}_rho_{F.g(rho)}_p_{F.g(p)}.png")

if __name__ == "__main__":
    opt_banner()
    opt_terrain()
    opt_terrain_single()
    opt_magnetogram()
    opt_hero()
