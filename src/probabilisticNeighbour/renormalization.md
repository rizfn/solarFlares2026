# Why the exponent slides but the ordering snaps

Two different things in this model are scale-free, and they behave in opposite
ways when we turn the knob $p$:

- the **size distribution** of spots and emissions is a power law whose exponent
  drifts continuously with $p$;
- the **spatial arrangement of signs** is mixed at small $p$ and segregated at
  large $p$, with a sharp transition at some $p_c$.

A power law with a continuously varying exponent is suspicious — usually an
exponent is nailed down by a fixed point, and fixed points are isolated. So the
first job is to explain how the exponent can move at all, and the second is to
explain why the ordering does *not* move with it. Both answers are
renormalization arguments, but in different spaces: one in mass, one in real
space.

---

## 1. Strip the model down

Forget the lattice for a moment and look at what the dynamics actually does.
Every occupied site carries a signed integer flux $\pm m$. One update picks an
occupied site **uniformly at random** and moves it to a random neighbour. If the
target is occupied, then:

- **same sign**: the two merge, $i, j \to i+j$. Pure coagulation.
- **opposite sign**: they cancel as far as they can, releasing an emission of
  size $s=\min(i,j)$ and leaving $|i-j|$ behind. Annihilation.

Monomers are injected to keep the density fixed. So this is

> signed-mass coagulation + annihilation + monomer injection.

The important detail is that the update picks a *site*, not a unit of flux. A
big spot hops exactly as often as a small one. In mean field that makes the
collision **kernel constant**: two spots of masses $i$ and $j$ meet at a rate
$\propto n_i n_j$ with no extra mass dependence. No mass dependence in the
kernel means no mass scale in the dynamics, which is the whole reason a power
law shows up.

Everything $p$ does, kinetically, is set one number: the **branching ratio**

$$
q \;=\; \Pr(\text{a collision is same-sign}) .
$$

Segregated surface (large $p$) gives $q\to1$: with big domains almost every
collision happens inside a domain, between like signs. An uncorrelated surface
would give $q=\tfrac12$. The measured value runs monotonically from $q=0.39$ at
$p=0$ to $q=0.95$ at $p=1$ — note that it starts *below* $\tfrac12$, because
bipolar injection drops a $+$ and a $-$ right next to each other and so leaves
the surface slightly **anti**-correlated rather than merely random.

---

## 2. Build the master equation

Let $n_m$ be the number density of $+$ spots of mass $m$, and
$N=\sum_{m\ge1} n_m$ the total density. By $\pm$ symmetry the $-$ spots have the
same distribution. Set the kernel to $1$.

**Loss.** A spot of mass $m$ stops being a mass-$m$ spot the moment it collides
with *anything*, same sign or not. Collisions happen at rate $n_m N$, so

$$ \text{loss} = n_m N . $$

**Gain, from below.** Two smaller $+$ spots merge into an $m$. That needs
$i+j=m$, both positive, and the collision must be same-sign, which costs a
factor $q$:

$$ \text{gain}_\uparrow = \frac{q}{2}\sum_{i+j=m} n_i n_j . $$

The $\tfrac12$ is because the sum runs over ordered pairs $(i,j)$ and $(j,i)$,
which are the same event.

**Gain, from above.** A *bigger* $+$ spot gets cut down to $m$: a $+(m+j)$ meets
a $-j$, annihilates $j$ of it, and leaves $+m$. This needs an opposite-sign
collision, factor $1-q$:

$$ \text{gain}_\downarrow = (1-q)\sum_{j\ge1} n_{m+j}\, n_j . $$

In steady state gain $=$ loss, so for $m\ge2$

$$
\frac{q}{2}\sum_{i+j=m} n_i n_j \;+\; (1-q)\sum_{j\ge1} n_{m+j} n_j \;=\; n_m N .
\tag{$\star$}
$$

Injection only enters the $m=1$ equation, so it does not appear here. Read
$(\star)$ as a tug of war in mass space: coagulation pushes spots **up**,
annihilation drags them **down**, and $q$ sets who wins.

---

## 3. The equation has no scale — so it cannot fix $\tau$

Rescale mass, $m \to bm$. Nothing in $(\star)$ objects: the kernel is a
constant, the sums are convolutions with no built-in mass, and the only special
mass ($m=1$, injection) sits at the very bottom. A power law

$$ n_m = A\,m^{-\tau} $$

maps to $A (bm)^{-\tau} = (Ab^{-\tau})\, m^{-\tau}$: same shape, same exponent,
different amplitude. The *form* is a fixed point of the rescaling for **every**
$\tau$. That is the first hint that $\tau$ is not pinned.

Let us check it directly on $(\star)$. Put $n_m = Am^{-\tau}$ with $\tau>1$ and
keep the leading behaviour at large $m$. Because $n_m$ falls steeply, each sum
is dominated by the terms where one partner is *small*:

$$
\frac{q}{2}\sum_{i+j=m} n_i n_j
\;\simeq\; \frac{q}{2}\Big( \underbrace{\sum_{i\ \rm small} n_i\, n_{m-i}}_{\simeq\, N n_m}
+ \underbrace{\sum_{j\ \rm small} n_j\, n_{m-j}}_{\simeq\, N n_m} \Big)
= q\, n_m N ,
$$

where we used $n_{m-i}\simeq n_m$ for small $i$ and $\sum_i n_i = N$. The same
trick on the annihilation sum, with $n_{m+j}\simeq n_m$:

$$
(1-q)\sum_{j\ge1} n_{m+j} n_j \;\simeq\; (1-q)\, n_m N .
$$

Add them:

$$
q\,n_m N + (1-q)\, n_m N \;=\; n_m N ,
$$

which is exactly the right-hand side of $(\star)$ — **for any $\tau$, and for
any $q$**. The leading-order equation is satisfied identically. It does not
select an exponent.

This is what "line of fixed points" means here. The RG flow in mass space has a
whole one-parameter family of scale-invariant solutions, and $q$ is an
**exactly marginal** coordinate: changing it does not generate a flow, it just
slides the system along the line. The exponent is free to move continuously
because nothing is pushing it back.

---

## 4. What does select $\tau$: the mass flux

The leading order cancelled, so the answer hides in the sub-leading part — the
piece that carries a **net current of mass** toward large $m$.

Take pure coagulation first, $q=1$. Monomers go in at the bottom, mass cascades
upward, and gets removed at a cutoff (the largest spot the box can hold). In
steady state the rate at which mass crosses any intermediate scale $m$ must be
the *same at every scale* — otherwise mass would pile up somewhere.

A coagulation event carries mass across scale $m$ when a spot of mass $i<m$
merges with one of mass $j > m-i$, so the product lands above $m$; the mass
carried across is $i$. With a constant kernel the rate of such events is
$n_i n_j$, so

$$
\Phi(m) \;=\; \int_0^m \!\! di \int_{m-i}^\infty \!\! dj \;\; i\, n_i\, n_j .
$$

Do the inner integral with $n_j = Aj^{-\tau}$, valid for $\tau>1$:

$$
\int_{m-i}^\infty A j^{-\tau} dj \;=\; \frac{A}{\tau-1}\,(m-i)^{1-\tau} ,
$$

so

$$
\Phi(m) \;=\; \frac{A^2}{\tau-1}\int_0^m i\cdot i^{-\tau}\,(m-i)^{1-\tau}\, di
\;=\; \frac{A^2}{\tau-1}\int_0^m i^{\,1-\tau}(m-i)^{1-\tau}\, di .
$$

Substitute $i = mx$, $di = m\,dx$, and pull every power of $m$ out:

$$
\Phi(m) = \frac{A^2}{\tau-1}\; m^{\,(1-\tau)+(1-\tau)+1} \int_0^1 x^{1-\tau}(1-x)^{1-\tau} dx
= \frac{A^2\,B(2-\tau,\,2-\tau)}{\tau-1}\; m^{\,3-2\tau} ,
$$

with $B$ the Beta function (finite as long as $\tau<2$). So

$$ \Phi(m) \;\propto\; m^{\,3-2\tau} . $$

Demanding that the flux not depend on $m$ kills the exponent:

$$ 3-2\tau = 0 \quad\Longrightarrow\quad \boxed{\tau_m = \tfrac32} \qquad (q=1). $$

This is the standard constant-kernel, constant-flux (Smoluchowski /
Kolmogorov–Zakharov) exponent. Note what happened: the *shape* was fixed by
scale invariance, but the *number* came from a conservation law imposed on top
of it.

---

## 5. Turning annihilation back on

For $q<1$ every opposite-sign collision destroys $2\min(i,j)$ of flux. Mass is
now leaking out at **every** scale, not just at the cutoff, so $\Phi(m)$ is no
longer constant — it decays with $m$. A flux that decays as you go up means
fewer big spots than the constant-flux solution would give, i.e. a **steeper**
distribution. So lowering $q$ raises $\tau_m$ above $3/2$.

The other end is $q=\tfrac12$: coagulation-up and annihilation-down are exactly
equally likely, so the net current through mass space is **zero**. With no
cascade there is nothing to sustain the $3/2$ solution, and the balance sits at

$$ \boxed{\tau_m \to 2} \qquad (q=\tfrac12,\ \text{well mixed}) . $$

$\tau_m=2$ is the marginal value where $\sum_m m\, n_m \sim \sum_m m^{-1}$ is
logarithmic — the borderline between a steady state dominated by mass and one
dominated by number. (This end is a heuristic, not a derivation; only $3/2$
comes out cleanly.)

So the accessible range of $q$ maps to a continuum

$$ \tau_m(p): \quad 2 \;\;(p=0) \;\longrightarrow\; \tfrac32 \;\;(p=1), $$

and the measured curve does exactly this — see `plots/exponents/`.

### Solving $(\star)$ numerically

We do not have to stop at the two endpoints: `meanField.py` integrates the rate
equation forward in time to its steady state (FFT convolutions, mass past the
cutoff $M$ leaves the system) and fits $\tau_m$ to the result. Feeding it the
$q$ measured in the simulation:

| $p$ | $q$ | $\tau_m$ sim | $\tau_m$ mean field |
|-----|-----|-----|-----|
| 0.0 | 0.386 | 2.01 | 7.9 |
| 0.2 | 0.459 | 1.84 | 2.30 |
| 0.4 | 0.533 | 1.73 | 1.88 |
| 0.6 | 0.619 | 1.61 | 1.72 |
| 0.8 | 0.760 | 1.51 | 1.59 |
| 1.0 | 0.952 | 1.49 | 1.49 |

The agreement is essentially exact where coalescence dominates ($q\gtrsim0.7$)
and degrades steadily as $q$ drops toward $\tfrac12$, where the mean field runs
away entirely. That failure is expected rather than embarrassing: at $q<\tfrac12$
the steady state is controlled by two-species annihilation $A+B\to\emptyset$,
whose upper critical dimension is $d_c=4$. In $d=2$ that reaction is strongly
correlated — the two species segregate into anticorrelated patches on their own —
so a well-mixed rate equation has no chance. The simulation stays near $\tau_m=2$
because the real 2D dynamics keeps the cascade alive where mean field kills it.

---

## 6. Emissions inherit the exponent

An emission is $s=\min(i,j)$ of a colliding pair. In mean field the pair is
drawn $\propto n_i n_j$, so

$$
P(s) \;\simeq\; 2 \cdot \underbrace{\frac{n_s}{N}}_{\Pr(i=s)}
\cdot \underbrace{\frac{1}{N}\sum_{j\ge s} n_j}_{\Pr(j\ge s)} ,
$$

the factor $2$ for either partner being the smaller. For $\tau_m>1$ the tail sum
is $\sum_{j\ge s} n_j \simeq \frac{A}{\tau_m-1} s^{1-\tau_m}$, so

$$
P(s) \;\propto\; s^{-\tau_m}\cdot s^{1-\tau_m} = s^{-(2\tau_m-1)}
\qquad\Longrightarrow\qquad
\boxed{\ \tau_s = 2\tau_m - 1\ } .
$$

The two exponents are not independent: the spot distribution is primary. As
$\tau_m$ runs $2\to\tfrac32$, $\tau_s$ runs $3\to2$.

---

## 7. The other renormalization: the sign field

Now throw away the masses and keep only $\sigma_i = \mathrm{sign}(s_i) \in
\{+,-,0\}$. This field has a $\mathbb{Z}_2$ symmetry ($+\leftrightarrow-$) and
an order parameter — the coarse-grained magnetization $m$ — that is zero when
mixed and nonzero when segregated. Neighbour injection and same-sign
coalescence align neighbours; bipolar injection and diffusion scramble them. So
coarse-grained, the sign field is an Ising-like system with some effective
coupling $K(p)$ that grows with $p$.

Coarse-grain by a factor $b$: tile the lattice with $b\times b$ blocks and
replace each block by the sign of its total flux. We want the effective coupling
$K'$ between neighbouring blocks. Two neighbouring blocks touch along a face of
$b^{d-1}$ bonds, each of strength $K$. But a bond only transmits coupling
between the *block signs* to the extent that each of its two spins is correlated
with its own block's sign. Call that correlation $c$; then

$$ K' \;\approx\; b^{\,d-1}\, c^2\, K . $$

Two limits, and both are easy:

- **$K$ small.** Spins inside a block are nearly independent, so the block sign
  is decided by a $\sqrt{b^d}$ fluctuation and any single spin only knows about
  it at the level $c \sim b^{-d/2}$. Then
  $$ K' \approx b^{\,d-1}\, b^{-d}\, K = \frac{K}{b} \;<\; K . $$
  The coupling shrinks under coarse-graining: flow to $K=0$, the mixed phase.

- **$K$ large.** The block is fully aligned, every spin agrees with the block
  sign, $c=1$, so
  $$ K' \approx b^{\,d-1} K = b\,K \;>\; K \quad (d=2). $$
  The coupling grows: flow to $K=\infty$, the segregated phase.

The recursion $K' = R(K)$ therefore has $R(K)/K$ going continuously from $1/b$
at small $K$ to $b$ at large $K$. It must pass through $1$ somewhere, so there
is a $K_c$ with

$$ R(K_c) = K_c , \qquad R'(K_c) > 1 . $$

That fixed point is **unstable**: nudge $K$ above $K_c$ and it runs to infinity,
nudge it below and it runs to zero. Nothing sits at $K_c$ unless you put it
there. Linearizing, $R'(K_c) \equiv b^{1/\nu}$, which is the usual statement
that the correlation length diverges as

$$ \xi \sim |K-K_c|^{-\nu} \sim |p - p_c|^{-\nu} . $$

The contrast with Section 3 is the whole point. There, the marginal coordinate
$q$ generated no flow, so every value of it was a fixed point and the exponent
could slide. Here $K$ is **relevant**: it flows, only one value is fixed, and it
has to be tuned. Same knob $p$, two completely different roles — marginal for
the size cascade, relevant for the spatial order. That is why the exponent
varies smoothly right through $p_c$ while the order parameter switches on.

Note also that this section says nothing about *avalanche* criticality. The
power law does not need $p=p_c$; it is there for every $p$, because it comes
from a driven cascade (drive at $m=1$, dissipation at the cutoff) and not from a
tuned critical point. Widening the box widens the inertial range without
touching $\tau$.

---

## 8. What the simulation says

$L=128$, $\rho=0.2$, exponents from the Hill estimator above $m,s \ge 10$:

| $p$ | $\tau_m$ (spot) | $\tau_s$ (emission) | $2\tau_m - 1$ |
|-----|-----|-----|-----|
| 0.0 | 2.01 | 3.14 | 3.01 |
| 0.2 | 1.84 | 2.80 | 2.68 |
| 0.4 | 1.73 | 2.60 | 2.46 |
| 0.6 | 1.61 | 2.42 | 2.21 |
| 0.8 | 1.51 | 2.07 | 2.02 |
| 1.0 | 1.49 | 1.83 | 1.97 |

Both endpoints land on the predicted values: $\tau_m = 2.01$ at $p=0$ (theory
$2$) and $1.49$ at $p=1$ (theory $3/2$), sliding smoothly in between.

$\tau_s = 2\tau_m-1$ works at $p=0$ (3.14 vs 3.01) and fails at $p=1$ (1.83 vs
1.97). That failure is informative: when the surface is segregated, emissions
only happen on the *interfaces* between domains, so the colliding pair is not
drawn from the bulk distribution at all. The derivation in Section 6 assumed it
was.

The segregation order parameter is flat and near zero up to $p\approx0.6$ and
climbs steeply after, while $\tau$ does not notice $p_c$ at all.

---

## 9. Where this is shaky

- Sections 2–6 are mean field, and Section 5 shows exactly where that bites: the
  rate equation is quantitative for $q\gtrsim0.7$ and useless as $q\to\tfrac12$,
  because two-species annihilation in $d=2$ is far below its upper critical
  dimension. The segregated phase is also emphatically *not* well mixed, which is
  where the emission relation breaks.
- Only $\tau_m=3/2$ is derived in closed form. The $\tau_m\to2$ end is a
  heuristic and the interior of the line comes from integrating $(\star)$
  numerically. A controlled treatment of the annihilation-dominated end needs a
  Doi–Peliti field theory for coagulation + annihilation + injection.
- At $p$ near $1$ a few domain-spanning spots form a condensate that sits
  outside the scaling body, so the fitted exponent depends on where you put the
  lower cutoff.
- The block-spin estimate in Section 7 gets the *structure* right (two stable
  phases, one unstable fixed point in between, $K$ relevant) but not the numbers.
  It says nothing about which universality class, and the injection makes this a
  driven, non-equilibrium version of Ising rather than the equilibrium one.
