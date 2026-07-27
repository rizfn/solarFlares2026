A generalized version.

Now, when new spots are added, there's a probability `p` of being added to a neighbour of the same sign.

`p=1` regenerates the "neighbour addition" from earlier, `p=0` is the random addition model.

Neighbour addition leads to coarsening and phase separation at a critical `p`. The power law's exponent, however, varies continuously depending on `p`.

## Contents

- `probabilisticNeighbour.cpp`, `run.sh` — the simulation and the parameter sweeps. Writes histograms and snapshots to `outputs/` (regenerable, not tracked).
- `plots.py` — all figures, into `plots/`. `plot_exponents_vs_p` is the summary: the exponent slides smoothly with `p` while the segregation order parameter switches on at `p_c ≈ 0.6`.
- `meanField.py` — numerical steady state of the mean-field rate equation, compared against the simulation at the measured same-sign fraction `q`.
- `renormalization.md` — why the exponent varies continuously (a marginal cascade, a line of fixed points) while the spatial ordering has a single tuned critical point.
- `tex/` — manuscript drafts. `main_revised.tex` is the current one; `main.tex` is the earlier version and its figures are missing.
