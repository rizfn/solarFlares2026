A generalized version.

Now, when new spots are added, there's a probability `p` of being added to a neighbour of the same sign.

`p=1` regenerates the "neighbour addition" from earlier, `p=0` is the random addition model.

Neighbour addition leads to coarsening and phase separation at a critical `p`. The power law's exponent, however, varies continuously depending on `p`.

## Contents

- `probabilisticNeighbour.cpp`, `run.sh` — the simulation and the parameter sweeps. Writes histograms and snapshots to `outputs/` (regenerable, not tracked).
- `plots.py` — all figures, into `plots/`. `plot_exponents_vs_p` is the summary: the exponent slides smoothly with `p` while the segregation order parameter switches on at `p_c ≈ 0.6`.
- `wellMixed.cpp` — the same cascade on a complete graph (Takayasu with two signs and biased partner selection), with `q` imposed as a rule instead of emerging from correlations. This is the exact stochastic process behind the rate equation, so it separates fluctuation effects from spatial ones.
- `meanField.py` — numerical steady state of the mean-field rate equation, the closed-form `tau_closed(q)`, and the three-way comparison against `wellMixed` and the 2D lattice.
- `renormalization.md` — why the exponent varies continuously (a marginal cascade, a line of fixed points) while the spatial ordering has a single tuned critical point. The sub-leading balance picks the point on the line: `cos(pi tau) = (1-q)/q`, with no solution at all for `q < 1/2`. Written at tutorial pace, every step spelled out.
- `renormalization_concise.md` — the same derivation at manuscript length, for dropping into the paper.
- `collisionPairs.cpp` — instrumented copy of the simulation that dumps the mass pair of every annihilation, used in Section 5.1 of `renormalization.md` to find which mean-field assumption `tau_s = 2 tau_m - 1` loses above `p_c`.
- `tex/` — manuscript drafts. `main_revised.tex` is the current one; `main.tex` is the earlier version and its figures are missing.
