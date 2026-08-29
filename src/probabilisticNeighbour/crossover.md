# Why the flare spectrum bends

## The model

A square lattice with periodic edges. Each occupied site holds a signed whole number
$\pm m$, a *spot*; the rest are empty. The number of spots is held fixed at $N=\rho L^2$.

One update picks an occupied site at random and moves its charge to a neighbour. If that
neighbour is empty, the spot moves. If it holds the same sign, the two merge into
$i+j$. If it holds the opposite sign, they cancel as far as they can, leaving $|i-j|$
and releasing a **flare** of size $s=\min(i,j)$. Whenever the spot count drops, one $+1$
and one $-1$ are injected: with probability $p$ each is placed beside a spot of its own
sign, otherwise both land on random empty sites.

Besides the density $\rho$ and the size $L$, the only parameter is $p$.

## The bend is not in the growth of spots

Flare sizes follow a broken power law — steep at small sizes, shallow at large ones,
with a bend at some size $s^*$. The bend is absent at $p\le0.5$ and at $p=1$, and
strongest near $p=0.8$.

It does not come from how spots grow. Record the smaller mass $\min(i,j)$ over *every*
collision, merging ones included, and that is a clean power law over four decades at
every $p$. Its slope is exactly what you would guess: spot masses are distributed as
$m^{-\tau_m}$ with $\tau_m\approx1.47$, and the smaller of two such masses is
distributed as $s^{-\alpha}$ with

$$\alpha = 2\tau_m-1 \approx 2 ,$$

which the data confirm to within $0.02$. Only the annihilating subset bends. Since the
flares are just that subset,

$$P_\text{flare}(s) = P_\text{collision}(s)\times f(s) ,$$

where $f(s)$ is the fraction of collisions at $\min(i,j)=s$ that annihilate rather than
merge. Everything about the bend lives in $f$.

## The shape of that fraction

$f(s)$ falls, turns at $s^*$, and rises again. **The two sides are not alike.** Measuring
the slope point by point along the curve (figure, panel a), the rising side sits at
$0.28$ across more than three decades, while the falling side is a hump: it steepens to
about $-0.7$ partway down and flattens again at both ends.

So the rising side is a power law and the falling side is not. Only the rising side has
an exponent; call it $\psi$, so that above the bend the flares fall as
$s^{-(\alpha-\psi)} \approx s^{-1.73}$. Quoting a matching exponent for the falling side
is meaningless, and doing so is what made that side appear to change with $p$: it was the
fitting range sliding along a curve.

**$\psi=0.28\pm0.01$ and it is universal.** The only thing it responds to is how wide a
stretch of the rising side is fitted. Pooling every run — sizes $L=256$ to $1024$,
densities $\rho=0.05$ to $0.6$, $p=0.6$ to $0.95$ — and sorting by the width $w$ of the
fit in decades:

| $w$ | 1–1.5 | 1.5–2 | 2–2.5 | 2.5–3 | 3–3.5 | 3.5–4.6 |
|---|---|---|---|---|---|---|
| $\psi$ | 0.180 | 0.215 | 0.240 | 0.263 | 0.261 | 0.282 |
| spread | ±.010 | ±.003 | ±.009 | ±.010 | ±.009 | ±.006 |

Each column mixes a twelvefold range of density, the whole range of $p$, and a factor
four in system size, yet $\psi$ agrees to $\pm0.01$ within a column. Narrow fits read low
because they sit near the turning point, where the falling side still contributes.
Taking only the widest fits gives $0.282\pm0.006$, which rules out $1/4$; $2/7=0.286$
sits comfortably inside.

## Why the fraction turns

> The crossover implies that collisions between opposite signs are likely when one is
> small, unlikely when both are medium, and likely again when both are large.

> The reason is the following. Space is phase separated into "red regions" and "blue
> regions". With probability ~0.5 * (1-p), a red spot is put in a blue region, leading to
> a small emission.

> However, with probaiblity `p`, the red is injected in the red region. How does it cause
> an emission? It needs to undergo a random walk, until it finally leaves the domain it's
> in. Once it arrives on the blue domain, it can annihilate and cause an emission.
> However, during the random walk, it's constantly growing, as it's merging with other
> reds, and newly injected reds (added with rate `p`) so that by the time it leaves the
> domain, it's large!

To test this, every spot is stamped at birth with the time, the place, and which rule
created it. The stamp travels with the spot and passes to the heavier partner when two
merge, so a spot keeps its identity as it grows. At $L=512$, $\rho=0.2$, $p=0.8$, the
spot that produces a flare of size $s$ has on average:

| $s$ | 1.8 | 12 | 121 | 1211 | 12106 |
|---|---|---|---|---|---|
| lived for (sweeps) | 0.6 | 8.7 | 45 | 169 | 604 |
| wandered (sites) | 1.1 | 5.5 | 13.3 | 25.4 | 46.6 |

Small flares come from spots that die where they were born; large ones from spots that
have lived a thousand times longer and wandered fifty sites. The distance goes as the
square root of the time at every $p$, so the wandering is an ordinary random walk.

Injecting into the wrong region does make the small flares. Spots dropped at random are
likelier to be destroyed than spots placed beside their own kind, and they release far
less when they are:

| $p$ | injections dropped at random | flares they cause | mean size released, random vs placed |
|---|---|---|---|
| 0.7 | 30% | 41.5% | 1.43 vs 1.90 |
| 0.8 | 20% | 34.1% | 1.49 vs 2.49 |
| 0.9 | 10% | 27.9% | 1.55 vs 3.74 |

At $p=0.9$ a randomly dropped spot is three and a half times likelier to be destroyed.
The excess is concentrated at sizes of one to three, so this starts the falling side but
does not account for all of it.

## What the picture predicts

Count a spot's collisions instead of its size. Each collision either merges it, so it
grows, or destroys it, so it flares. Two rates follow, both measured directly
(figure, panels b and c):

- how fast its chance of being destroyed climbs once it is established;
- how fast its size climbs, roughly size $\propto$ (collisions)$^{1.5\text{ to }1.7}$.

Dividing one by the other converts a slope in collisions into a slope in size, and gives
the rising side:

| $p$ | 0.70 | 0.75 | 0.80 | 0.85 | 0.90 | 0.95 |
|---|---|---|---|---|---|---|
| predicted $\psi$ | 0.225 | 0.221 | 0.261 | 0.243 | 0.236 | 0.254 |
| measured $\psi$ | 0.214 | 0.225 | 0.264 | 0.277 | 0.257 | 0.271 |

Agreement to about $0.02$ throughout. **$\psi$ is therefore not a quantity in its own
right** — it is one rate divided by another, and it is the same everywhere because both
of those are. This is not merely a change of units: spots that have survived the same
number of collisions have wildly different sizes, so there is no fixed conversion between
the two. That the slopes convert cleanly regardless is the result.

## What is still missing

Where the bend sits is not explained. The number of collisions a spot must survive to be
at its safest falls steadily with $p$ — 1211, 562, 383, 178, 121, 56, 38, 26 from $p=0.6$
to $0.95$ — but it is not a power of $1-p$: the fitted exponent slides from 2.5 to 1.5
depending on which points are used. The same is true of $s^*$ itself, so there may be no
clean law there to find.

Every attempt to bring in the size of the red and blue regions has failed in the same
direction. If a spot became vulnerable once it had crossed one region, the safest
collision count would grow as those regions grow; it shrinks instead, and the distance
covered at the bend shrinks too, from 30 sites at $p=0.6$ to 9 at $p=0.9$. The region
size moves opposite to everything it should explain, which suggests it is the wrong
variable rather than a badly measured right one.

$\tau_m$, the spread of spot masses, is not pinned down either. Over three and a half
decades $m^{-1.47}$ and $m^{-3/2}(\ln m)^{0.2}$ describe it equally well. At $p=1$ there
is no annihilation and the process is plain aggregation with injection, which gives
$3/2$; at $p=0.8$ the mean-field prediction is 1.61 against 1.47 measured, so the
two-dimensional value is neither.

$s^*$ is the same at $L=256$, 512 and 1024 for $0.70\le p\le0.90$, so the bend is real
and survives to infinite size. Outside that range it drifts with $L$, and at $p=1$ there
is no bend at all.
