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
- `plots.py` — all figures, into `plots/`.
  - `plot_voronoi_same_sign` — `q_V(p)` at two densities.
  - `plot_exponents_vs_p` — both cascade exponents, the mean-field prediction
    `tau_s = 2 tau_m - 1`, and `q_V` on the right-hand axis.
  - `plot_histograms` — the spot-size and emission distributions.
  - `plot_slices` — the `z = L/2` plane through the cube, since a 3D snapshot cannot be
    shown directly.
