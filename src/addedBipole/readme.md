Injection *correlation* instead of injection sign-sorting.

The cascade is unchanged from `probabilisticNeighbour`: spots random-walk, same signs
coagulate (`i,j -> i+j`), opposite signs partially annihilate and emit `s = min(i,j)`,
and the spot count is held fixed by injecting a `+1/-1` pair whenever one is lost. The
only parameter is *where that pair goes*:

- with probability `r`, the `+` and `-` are placed as lattice **neighbours** (a bipole),
- with probability `(1-r)`, they are placed at two **independent** random empty sites.

`r=1` is the `p=0` limit of `probabilisticNeighbour`; `r=0` is new.

## Why this parameter

The kinetic theory is controlled by `q`, the fraction of collisions that are same-sign.
The closed form derived in `probabilisticNeighbour/renormalization.md`,

    cos(pi tau) = (1 - q) / q,

has **no root at all for `q < 1/2`**: once annihilation outweighs coagulation, the
mean-field rate equation has no scale-free steady state. The well-mixed Monte Carlo
agrees — below `q=1/2` it returns a steep cutoff-limited fit that drifts with `q`,
which is the signature of no power law rather than of a large exponent.

An injected bipole is maximally anti-correlated: the two halves are born adjacent and
tend to annihilate each other before either finds anything else, which suppresses `q`.
Measured at `L=128`, `rho=0.2`:

    r = 0  ->  q = 0.500
    r = 1  ->  q = 0.388

So the whole sweep sits at or below the mean-field boundary. That is the point. In 2D
the two halves of a bipole can drift apart before annihilating and merge with other
spots instead, a channel the rate equation cannot represent, so a power law may survive
where mean field forbids one. Whether it does is what the sweep measures — and `q` is
measured from the collision counts rather than imposed, so the comparison is direct.

## What it does

A power law survives, over the whole forbidden window. Measured from the cumulative
`P(m > x) ~ x^(1-tau)`, which is far less noisy than the binned density:

| `r` | `q` | `tau_m` | `tau_s` |
|---|---|---|---|
| 0.0 | 0.500 | 1.96 | 3.00 |
| 0.5 | 0.449 | 1.98 | 3.03 |
| 1.0 | 0.389 | 1.99 | 3.14 |

Three things make this a real exponent rather than a cutoff-limited fit:

- the **local log-slope is flat** — `1.93`–`1.96` over three decades at `rho=0.8`, where
  the dynamic range is longest, with a scatter of `0.02`;
- it **does not drift with `L`** — `tau_m` = `1.978, 1.961, 1.972, 1.957` at
  `L = 64, 128, 256, 512` for `r=0`, while the cutoff `m_max` moves by a factor of 23;
- it **does not drift with density** either, `1.96`–`1.98` over `rho = 0.2`–`0.8`.

And it does not vary with `q`. That is the substantive point. For `q > 1/2` the exponent
slides continuously along the line `cos(pi tau) = (1-q)/q`; below `q = 1/2` it **locks**
at `tau_m = 2`, `tau_s = 3`, which are exactly that line's endpoint values (`q = 1/2`
gives `cos(pi tau) = 1`, so `tau = 2`). The line of fixed points terminates, and the
system sits at the terminus for every `q` below it rather than continuing to move.

So 2D drift does restore a scale-free steady state where the rate equation has none —
but it buys no new exponent. This is the annihilation-dominated regime, controlled by
two-species annihilation `A + B -> 0` with `d_c = 4`, so 2D is below the critical
dimension and the mean-field equation is simply the wrong description, not a nearly-right
one. `tau_s` tracks `2 tau_m - 1` at `r = 0` and drifts `~0.15` above it by `r = 1`; the
system is mixed at every `r` here (no segregation), so the identity is expected to hold,
and whether that residual is real or short-window bias is not settled.

## Contents

- `addedBipole.cpp`, `run.sh` — the simulation and the sweep. Writes spot-size and
  emission histograms, snapshots, and the measured `q` to `outputs/` (regenerable, not
  tracked).
- `plots.py` — all figures, into `plots/`.
  - `plot_local_slopes` is the decisive one: a power law is a horizontal line, and
    curvature is a cutoff-limited fit masquerading as an exponent.
  - `plot_tau_vs_q` places the measured points against the closed form, in the window
    where the closed form does not exist.
  - `plot_finite_size` is the other honest test — a cutoff-limited fit drifts with `L`,
    a real power law does not.
