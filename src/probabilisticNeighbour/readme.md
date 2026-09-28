A generalized version.

Now, when new spots are added, there's a probability `p` of being added to a neighbour of the same sign.

`p=1` regenerates the "neighbour addition" from earlier, `p=0` is the random addition model.

Neighbour addition leads to coarsening and phase separation at a critical `p`. The power law's exponent, however, varies continuously depending on `p`.

## Contents

- `probabilisticNeighbour.cpp`, `run.sh` — the simulation and the parameter sweeps. Writes histograms and snapshots to `outputs/` (regenerable, not tracked).
- `plots.py` — all figures, into `plots/`. `plot_exponents_vs_p` is the summary: the exponent slides smoothly with `p` while the surface segregates above `p ≈ 0.7`. That threshold is read off `plot_correlation_function` (C(r) collapses against `r` at `p=0.6`, against `r/L` from `p=0.8` up) and `plot_qv_finite_size` (curves separate above `p ≈ 0.65`); there is no finite-size-scaling determination of `p_c`, so the number is a bracket, not a measurement.
- `wellMixed.cpp` — the same cascade on a complete graph (Takayasu with two signs and biased partner selection), with `q` imposed as a rule instead of emerging from correlations. This is the exact stochastic process behind the rate equation, so it separates fluctuation effects from spatial ones.
- `meanField.py` — numerical steady state of the mean-field rate equation, the closed-form `tau_closed(q)`, and the three-way comparison against `wellMixed` and the 2D lattice.
- `renormalization.md` — why the exponent varies continuously (a marginal cascade, a line of fixed points) while the spatial ordering has a single tuned critical point. The sub-leading balance picks the point on the line: `cos(pi tau) = (1-q)/q`, with no solution at all for `q < 1/2`. Written at tutorial pace, every step spelled out.
- `renormalization_concise.md` — the same derivation at manuscript length, for dropping into the paper.
- `collisionPairs.cpp` — instrumented copy of the simulation that dumps the mass pair of every annihilation, for finding which mean-field assumption `tau_s = 2 tau_m - 1` loses above `p_c`.
- `emissionChannels.cpp`, `emissionAge.cpp`, `plot_crossover.py`, `plot_mechanism.py` — the crossover in the emission spectrum: `P_emis = P_coll * (1 - q(s))`, and what sets the scale where `1 - q(s)` turns. Written up in `crossover.md`.
- `plot_runningQ.py` — tests whether the closed form still holds locally once `q` is allowed to run with `s`. Partly: see its header.
- `paperFigures.py` — the five publication figures, into `plots/paper/`. It reads the big server sweep in `../probabilisticNeighbourServer/outputs` (7400 runs, up to `L=1024`) as well as this directory's well-mixed and channel runs, and caches every histogram per seed on a fixed log grid in `outputs/figcache/`, so only the first run is slow. Exponents are fitted per seed over a window taken from the pooled histogram, so the error bars are the seed-to-seed spread. Four fixed colours, everything else tinted from them: red and blue are the two signs of a spot as on a snapshot, `#ad8350` marks every spot-size quantity and `#5f0f40` every flare quantity, so a panel's hue says which distribution it is about; each ramp steps hue as well as lightness so neighbouring curves stay distinct. The density figure uses `../probabilisticNeighbourServer/run_density.sh` (`L=512`, `rho=0.05...0.9`).
