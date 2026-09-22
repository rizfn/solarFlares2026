Correlated injection, with the uncorrelated channel replaced by a bipole.

The cascade is identical to `probabilisticNeighbour`: spots random-walk on a periodic
square lattice, same signs coagulate (`i,j -> i+j`), opposite signs partially annihilate
and emit `s = min(i,j)`, and the spot count is held fixed by injecting a `+1/-1` pair
whenever one is lost. Only *where the pair goes* differs, and only in the `(1-p)` branch:

|  | probability `p` | probability `1-p` |
|---|---|---|
| `probabilisticNeighbour` | `+` next to a `+`, `-` next to a `-` | two **independent** random empty sites |
| this model | `+` next to a `+`, `-` next to a `-` | two **adjacent** empty sites (a bipole) |

`p=1` is the pure neighbour model in both cases, so the two only differ as `p` falls.

## The question

In `probabilisticNeighbour` the annihilating fraction

    1 - q(s) = P_emis(s) / P_coll(s)

is U-shaped: it falls at small `s`, turns at a crossover scale `s*`, and rises again
across three decades. `crossover.md` attributes the *falling* branch to the uncorrelated
channel dropping a lone minority monomer into a wrong-sign domain, where it annihilates
almost immediately and emits `s ~ 1`. Supporting numbers there: at `p=0.9` a randomly
dropped spot is three and a half times likelier to be destroyed than a placed one, and
the excess is concentrated at `s = 1` to `3`.

If that attribution is right, `s*` is an injection artefact and not a property of the
cascade — which matters, because it is the one feature of the 2D model with no law
behind it, and because the emission exponent is only clean above it.

A bipole cannot produce the artefact. The pair is its own nearest opposite sign, so it
either self-annihilates at once — emitting `s = 1`, restoring the count, a null event —
or the two halves separate and join the bulk as ordinary spots. Neither outcome drops a
lone minority charge into a hostile domain.

## What it does

The hypothesis is **wrong**, and the replacement is better.

The U in `1 - q(s)` does not weaken. It gets *deeper*, its falling branch gets *steeper*,
and `s*` moves *down* by a factor of 2 to 3 at every `p`. A bipole is its own nearest
opposite sign, so it self-annihilates at `s=1` and supplies its own small-`s` excess.
Both injection rules contaminate small `s`; neither is a clean control. So
`crossover.md`'s attribution of the falling branch to wrong-domain monomer drops does not
hold.

But comparing at fixed `p` was the wrong comparison, because the injection rule also
changes `q`. Against the measured `q` the two models **collapse onto one curve** — in the
exponent, the crossover scale, the depth of the U and `psi` alike (`L=1024`, `rho=0.2`):

| `q` | `tau_m` | `s*` | `(1-q)` at min | `psi` |
|---|---|---|---|---|
| 0.9439 (bip)  | 1.484 | 81.6 | 0.0241 | 0.275 |
| 0.9499 (mono) | 1.481 | 82.4 | 0.0254 | 0.264 |
| 0.9617 (bip)  | 1.477 | 52.4 | 0.0216 | 0.265 |
| 0.9704 (mono) | 1.473 | 53.6 | 0.0207 | 0.266 |

Two injection rules that differ as much as two rules can — one drops isolated charges into
random domains, the other places locally neutral dipoles — agree to a few percent whenever
they produce the same `q`. Tight for `q >~ 0.87`, loosening to a factor `~1.5` in `s*` by
`q ~ 0.65`, and the same at `L=256`, 512 and 1024. `p=1` is identical in both models by
construction (`q=0.9860`, `tau_m=1.466`, `tau_s=1.724`), which is the sweep's built-in
consistency check.

So `s*` is not an injection artefact and not a property of the injection rule at all: like
`tau_m`, it is a function of `q`.

One thing does **not** collapse across system size: the depth. The value of `1 - q` at the
minimum falls as roughly `L^-0.7` from `L=256` to 1024, at the same rate in both models,
while the small-`s` end and the position `s*` stay put — so the U deepens as the box grows.
Unexplained.

Note that the bipole also suppresses `q` on its own — see `addedBipole`, where pure
bipole injection drives `q` from `0.500` down to `0.388` and the exponent locks at the
`q = 1/2` endpoint `tau_m = 2`. So at small `p` this model is expected to sit at that
terminus, and the interesting window is `p >~ 0.6`, where the correlated channel
dominates and `q` is well above `1/2`.

## Contents

- `correlatedOrBipole.cpp`, `run.sh` — the simulation and the sweep. The sweep matches
  `probabilisticNeighbour/run_channels.sh` in `p`, `L`, `rho` and sweep count so the two
  models can be overlaid directly. Writes to `outputs/` (regenerable, not tracked).
- `plots.py` — all figures, into `plots/`. Every one overlays the two models.
  - `plot_qOfS` — `1 - q(s)` side by side, with `s*` marked.
  - `plot_collapse_vs_q` is the decisive one: `tau_m` and `s*` against the measured `q`,
    where the two models fall on a single curve.
  - `plot_crossover_vs_p` — `s*(p)` at three `L`, against the monomer model's curve.
  - `plot_exponents_vs_p` — `tau_m` and `tau_s`, to check the cascade is unchanged.
  - `plot_spectra` — the emission spectra, raw and compensated by `s^2`.
