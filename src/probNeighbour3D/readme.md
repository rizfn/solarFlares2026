The `probabilisticNeighbour` model in three dimensions.

The rules are unchanged from 2D. Spots random-walk on a periodic cubic lattice, same
signs coagulate (`i,j -> i+j`), opposite signs partially annihilate and emit
`s = min(i,j)`, and the spot count is held fixed by injecting a `+1/-1` pair whenever one
is lost:

- with probability `p`, the `+` is placed next to an existing `+` and the `-` next to an
  existing `-` (co-localized addition),
- with probability `(1-p)`, the two land on independent random empty sites.

`p=0` is random addition, `p=1` is pure neighbour addition. The only change from the 2D
code is the lattice: six neighbours instead of four, and the site index `x + L*y + L*L*z`.

## Measuring the correlation

Same idea as in 2D: with empty sites everywhere, lattice adjacency is not a useful
notion of "neighbour", so the neighbour graph is built geometrically. In 3D the Voronoi
neighbours of a spot are its Delaunay edges, obtained from the tetrahedra of a Delaunay
triangulation of the spot positions. `q_V` is the fraction of Voronoi-neighbour pairs
sharing a sign: `1/2` for a mixed system, `-> 1` when segregated.

Two 3D-specific details:

- Periodicity is handled by replicating only the points within a few mean spacings of
  each face, not all 26 image copies. Full replication is affordable in 2D (9 copies) but
  not in 3D (27 copies of a much larger point set).
- A spot has ~15 Voronoi neighbours in 3D, against ~6 in 2D, so a given domain structure
  gives a *lower* `q_V` in 3D than in 2D: more of each spot's neighbours sit across a
  domain wall. A fully segregated test configuration (two slabs) at `L=24`, `rho=0.2`
  returns `q_V = 0.88`, not 1. Compare 3D numbers to 3D numbers, not to the 2D curve.

The estimator is checked against two configurations with known answers: random signs give
`q_V = 0.499` (should be 1/2), and two slabs give the 0.88 above.

## Contents

- `probNeighbour3D.cpp`, `run.sh` — the simulation and the sweep. Writes spot-size and
  emission histograms, snapshots, and largest-spot trajectories to `outputs/`
  (regenerable, not tracked).
- `emissionChannels3D.cpp`, `run_channels3D.sh`, `plot_crossover3D.py` — the psi-vs-dimension test below.
- `plots.py` — all figures, into `plots/`.
  - `plot_voronoi_same_sign` — `q_V(p)` at two densities.
  - `plot_exponents_vs_p` — both cascade exponents, the mean-field prediction
    `tau_s = 2 tau_m - 1`, and `q_V` on the right-hand axis.
  - `plot_histograms` — the spot-size and emission distributions.
  - `plot_slices` — the `z = L/2` plane through the cube, since a 3D snapshot cannot be
    shown directly.

## Is psi a function of dimension?

`emissionChannels3D.cpp` is the 3D copy of `../probabilisticNeighbour/emissionChannels.cpp`:
it records `min(i,j)` over all collisions and over annihilations, so their ratio is the
annihilating fraction `1 - q(s)` (see `../probabilisticNeighbour/crossover.md`). In 2D its
rising side has slope `psi = 0.28` everywhere. `run_channels3D.sh` sweeps `p` at `L=64`
(same site count as 2D `L=512`) and `L=48`, `rho=0.2`, 30000 sweeps.
`plot_crossover3D.py` compares the two dimensions, into `plots/crossover/`.

**psi does not depend on dimension.** The 3D local slope sits on a plateau at 0.27-0.28,
the 2D value. Fitted from `4 s*` to `0.02 rho L^3`:

| p | 0.7 | 0.75 | 0.8 | 0.85 | 0.9 | 0.95 | 1.0 |
|---|---|---|---|---|---|---|---|
| 3D, L=64 | 0.264 | 0.274 | 0.272 | 0.273 | 0.263 | 0.273 | 0.274 |
| 3D, L=48 | — | — | 0.269 | 0.270 | 0.275 | 0.274 | 0.269 |

Varying the lower cut from 3 to 6 `s*` and the upper cut from 0.01 to 0.04 `rho L^3`
moves these by at most 0.01. The upper cut matters: in 3D the slope rolls off above
~0.02 of the total mass `rho L^3`, at the same `s / rho L^3` for both `L`, so it is a
finite-size effect. Fitting up to the largest collision (the 2D recipe) runs into that
roll-off and reads 0.22-0.25. The box also leaves only 1-1.7 decades of plateau.

What does change with dimension is where the U sits. `s*` is 5 to 20 times smaller in
3D (17 against 120 at `p=0.8`), is the same at `L=48` and 64, and a U is still there at
`p=0.5` and `p=1`, where the 2D one is gone.

The mean-field baseline `alpha = 2 tau_m - 1` for the all-collision spectrum holds to
about 0.04 in 3D (same plateau estimator on both sides). On that estimator the 2D `L=512`
runs are off by about 0.09, more than the 0.02 quoted in `crossover.md`, which used a
different fit.
