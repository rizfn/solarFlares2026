# Exponent selection in a driven coagulation–annihilation cascade

Manuscript-length companion to `renormalization.md`.

## Kinetics

Each occupied site carries a signed integer flux $\pm m$. An update selects an
occupied site uniformly and moves it to a random neighbour; like signs coalesce
($i,j\to i+j$), opposite signs partially annihilate, emitting $s=\min(i,j)$ and
leaving $|i-j|$. Monomers are injected to hold the density fixed.

Because the update selects a *site* rather than a unit of flux, all spots hop at
the same rate and the mean-field collision kernel is constant: masses $i,j$ meet
at a rate $\propto n_in_j$ with no further mass dependence. The injection
correlation $p$ enters the kinetics only through the branching ratio

$$ q = \Pr(\text{a collision is same-sign}) , $$

equal to $\tfrac12$ for an uncorrelated surface and $\to1$ when the signs
segregate. In simulation $q$ runs from $0.50$ at $p=0$ — exactly the
uncorrelated value, since injection there is at two independent random sites —
to $0.95$ at $p=1$.

## Rate equation

Let $n_m$ be the density of $+$ spots of mass $m$ and $N=\sum_{m\ge1}n_m$; the
$-$ spots are identical by symmetry. A mass-$m$ spot is lost at rate $n_mN$ on
any collision, and gained by same-sign coalescence of two smaller spots or by
annihilation of a larger one. In steady state, for $m\ge2$,

$$
\frac{q}{2}\sum_{i+j=m}n_in_j \;+\; (1-q)\sum_{j\ge1}n_{m+j}n_j \;=\; n_mN .
\tag{1}
$$

Injection enters only the $m=1$ equation.

## Marginality

Insert $n_m=Am^{-\tau}$, $\tau>1$. Since $n_m$ falls steeply, each convolution is
dominated by configurations in which one partner is $O(1)$ and the other carries
almost all of $m$. There $n_{m\mp i}\simeq n_m$ factors out, and summing the
small partner over $\sum_in_i=N$ leaves $n_mN$ per such configuration. The
coalescence sum has two of them — the small partner may be either member of the
pair — and the annihilation sum one, so

$$ \frac{q}{2}\,(2n_mN) + (1-q)\,(n_mN) = n_mN $$

identically, for every $\tau$ and every $q$. We refer to this leading
approximation as the **small-partner estimate**. It requires only that $n_m$
fall steeply, so it selects nothing. Equivalently, $m\to bm$ maps a power law to
a power law of the same exponent — $\tau$ is a marginal direction, the linearized
RG cannot fix it, and the steady states form a line. Selection resides in the
first sub-leading term.

## Sub-leading balance

With $n(m)=Am^{-\tau}$ for $m\ge1$, so $N=A/(\tau-1)$, write Eq. (1) as
$\tfrac{q}{2}C(\mu)+(1-q)K(\mu)=n(\mu)N$ with

$$
C(\mu)=\int_1^{\mu-1}\!\! n(i)n(\mu-i)\,di ,
\qquad
K(\mu)=\int_1^{\infty}\!\! n(\mu+j)n(j)\,dj .
$$

Rescaling $i=\mu x$ and $j=\mu y$ confines all $\mu$-dependence to a prefactor
and a cutoff $\delta=1/\mu$:

$$
C = A^2\mu^{1-2\tau}\!\!\int_\delta^{1-\delta}\!\! x^{-\tau}(1-x)^{-\tau}dx ,
\qquad
K = A^2\mu^{1-2\tau}\!\!\int_\delta^{\infty}\!\! y^{-\tau}(1+y)^{-\tau}dy .
\tag{2}
$$

Each integrand is non-integrable only where a bare $x^{-\tau}$ or $y^{-\tau}$
stands, and the divergence is confined to a single power: expanding the regular
factor,

$$ x^{-\tau}(1-x)^{-\tau} = x^{-\tau} + \tau x^{1-\tau} + O(x^{2-\tau}) , $$

only $\int_\delta x^{-\tau}dx=\delta^{1-\tau}/(\tau-1)+\mathrm{const}$ diverges,
the next term converging because $\tau<2$. $C$ has two such endpoints ($x=0,1$),
$K$ one ($y=0$; its tail $\sim y^{-2\tau}$ converges for $\tau>1$) — hence the
factor $2$ below. Subtracting the offending powers,

$$
x^{-\tau}(1-x)^{-\tau}-x^{-\tau}-(1-x)^{-\tau} \;\sim\; \tau x^{1-\tau}-1 ,
\qquad
y^{-\tau}\big[(1+y)^{-\tau}-1\big] \;\sim\; -\tau y^{1-\tau} ,
$$

leaves integrands integrable throughout for $1<\tau<2$, while the removed powers
integrate elementarily. Hence

$$
\int_\delta^{1-\delta}\!\!\! x^{-\tau}(1-x)^{-\tau}dx = \frac{2\delta^{1-\tau}}{\tau-1}+B_c+O(\delta^{2-\tau}),
\quad
\int_\delta^{\infty}\!\!\! y^{-\tau}(1+y)^{-\tau}dy = \frac{\delta^{1-\tau}}{\tau-1}+B_a+O(\delta^{2-\tau}),
$$

with

$$
B_c=\int_0^1\!\Big[x^{-\tau}(1-x)^{-\tau}-x^{-\tau}-(1-x)^{-\tau}\Big]dx-\frac{2}{\tau-1},
\quad
B_a=\int_0^\infty\!\! y^{-\tau}\Big[(1+y)^{-\tau}-1\Big]dy ,
$$

the extra constant in $B_c$ arising from the finite upper limits that $B_a$
lacks.

Both constants reduce to Beta functions, by different routes.

*Evaluating $B_c$.* For $\tau<1$ the two subtracted powers are separately
integrable, $\int_0^1x^{-\tau}dx=\int_0^1(1-x)^{-\tau}dx=1/(1-\tau)$, so together
they contribute $-2/(1-\tau)$ — exactly cancelling the explicit $-2/(\tau-1)$ in
the definition. There $B_c$ is thus the unsubtracted integral, which is the first
Euler form

$$
\int_0^1 x^{a-1}(1-x)^{b-1}dx = B(a,b) = \frac{\Gamma(a)\Gamma(b)}{\Gamma(a+b)}
$$

with $a-1=b-1=-\tau$, i.e. $a=b=1-\tau$, giving
$B_c=\Gamma(1-\tau)^2/\Gamma(2-2\tau)$. The subtracted definition converges for
all $\tau<2$ and this expression is analytic there; agreeing on $0<\tau<1$, they
agree on $1<\tau<2$.

*Evaluating $B_a$.* That route is unavailable here: as $y\to\infty$ the bracket
tends to $-1$ and the integrand to $-y^{-\tau}$, so $B_a$ exists only for
$\tau>1$ and never overlaps a convergent Euler integral. Integrate by parts
instead, with

$$
u=(1+y)^{-\tau}-1 , \quad du=-\tau(1+y)^{-\tau-1}dy ,
\qquad
dv=y^{-\tau}dy , \quad v=\frac{y^{1-\tau}}{1-\tau} .
$$

The boundary term $[uv]_0^\infty$ vanishes at both ends: as $y\to0$,
$u\simeq-\tau y$ and $uv\sim y^{2-\tau}\to0$ since $\tau<2$; as $y\to\infty$,
$u\to-1$ and $v\sim y^{1-\tau}\to0$ since $\tau>1$. Hence

$$
B_a = -\int_0^\infty v\,du
= \frac{\tau}{1-\tau}\int_0^\infty y^{1-\tau}(1+y)^{-\tau-1}dy .
$$

The remaining integral is the second Euler form
$\int_0^\infty y^{a-1}(1+y)^{-a-b}dy=B(a,b)$, with $a-1=1-\tau$ and
$a+b=\tau+1$, i.e. $a=2-\tau$ and $b=2\tau-1$ — both positive on $1<\tau<2$, so
it converges outright, with no continuation. Using
$\Gamma(2-\tau)=(1-\tau)\Gamma(1-\tau)$ and $\Gamma(\tau+1)=\tau\Gamma(\tau)$,
the prefactor cancels:

$$
B_a = \frac{\tau}{1-\tau}\cdot\frac{\Gamma(2-\tau)\Gamma(2\tau-1)}{\Gamma(\tau+1)}
= \frac{\tau}{1-\tau}\cdot\frac{(1-\tau)\Gamma(1-\tau)\Gamma(2\tau-1)}{\tau\,\Gamma(\tau)}
= \frac{\Gamma(1-\tau)\Gamma(2\tau-1)}{\Gamma(\tau)} .
$$

Collecting,

$$
B_c(\tau)=\frac{\Gamma(1-\tau)^2}{\Gamma(2-2\tau)} ,
\qquad
B_a(\tau)=\frac{\Gamma(1-\tau)\Gamma(2\tau-1)}{\Gamma(\tau)} .
$$

Using $\delta^{1-\tau}\mu^{1-2\tau}=\mu^{-\tau}$ and $N=A/(\tau-1)$,

$$
C(\mu)=2n(\mu)N + A^2B_c\,\mu^{1-2\tau} ,
\qquad
K(\mu)=n(\mu)N + A^2B_a\,\mu^{1-2\tau} .
\tag{3}
$$

The first terms reproduce the small-partner estimate; the second are what it
discards, smaller by $\mu^{1-\tau}$.

Substituting into Eq. (1), the $n(\mu)N$ terms cancel the right-hand side
identically and what survives is

$$
\frac{q}{2}B_c(\tau) + (1-q)B_a(\tau) = 0 .
\tag{4}
$$

Nothing else can absorb this residue: it carries $\mu^{1-2\tau}$ against the
leading $\mu^{-\tau}$; $A^2$ divides out, so no amplitude, density or injection
rate enters; and an added component $A'm^{-\tau'}$ with $\tau'>\tau$ cancels in
the leading balance and contributes only at $\mu^{1-\tau-\tau'}\ll\mu^{1-2\tau}$.
Equation (4) is a solvability condition in the single remaining unknown.

## Closed form

With $\Gamma(z)\Gamma(1-z)=\pi/\sin\pi z$ applied to $z=\tau$ and $z=2\tau-1$,
and $\sin(2\pi\tau-\pi)=-\sin2\pi\tau$, Eq. (4) becomes

$$
\frac{q}{1-q} = -2\,\frac{\Gamma(2\tau-1)\Gamma(2-2\tau)}{\Gamma(\tau)\Gamma(1-\tau)}
= \frac{2\sin\pi\tau}{\sin2\pi\tau} = \frac{1}{\cos\pi\tau} ,
$$

that is,

$$ \cos(\pi\tau_m) = \frac{1-q}{q} ,
\qquad
\tau_m(q) = 2-\frac{1}{\pi}\arccos\!\left(\frac{1-q}{q}\right) .
\tag{5}
$$

The root is unique: the derivation requires $1<\tau<2$ ($N$ diverges below,
$B_c$ and $B_a$ above), and $\cos\pi\tau$ increases monotonically from $-1$ to
$1$ across that interval. For $q\in[\tfrac12,1]$ Eq. (5) gives
$\tau_m\in[\tfrac32,2]$.

**Limits.** $q=1$ gives $\tau_m=\tfrac32$, the constant-flux (Takayasu /
Smoluchowski) exponent for constant-kernel aggregation with injection — obtained
here without invoking flux conservation, and reachable from Eq. (4) directly,
since $q=1$ requires $B_c=0$, i.e. $\Gamma(2-2\tau)=\infty$. $q=\tfrac12$ gives
$\tau_m=2$, where coalescence and annihilation carry equal and opposite mass
current.

**No solution for $q<\tfrac12$.** Eq. (5) would require $\cos\pi\tau_m>1$: no
power law solves Eq. (1) once annihilation outweighs coalescence. Direct
integration of Eq. (1) confirms this, relaxing to a cutoff-dominated rather than
scale-free steady state.

## Emissions

An emission is $s=\min(i,j)$ for a pair drawn $\propto n_in_j$, so
$P(s)\simeq2(n_s/N)(N^{-1}\sum_{j\ge s}n_j)\propto s^{-\tau_m}s^{1-\tau_m}$ and

$$
\tau_s = 2\tau_m-1 = 3-\frac{2}{\pi}\arccos\!\left(\frac{1-q}{q}\right) .
\tag{6}
$$

## Comparison

Eq. (5) is tested against three numerical routes: a well-mixed Monte Carlo of the
same process (complete graph, $q$ imposed as a rule, $N=2\times10^5$ spots),
direct integration of Eq. (1) to steady state, and the 2D lattice model
($L=128$, $\rho=0.2$). Exponents from the Hill estimator above $m\ge10$, which is biased
a few hundredths high by the pre-asymptotic region and the cutoff; `plots.py` uses a
windowed fit instead, and the figures carry those values, not these.

| $q$ | Eq. (5) | well-mixed MC | rate equation | 2D lattice |
|-----|-----|-----|-----|-----|
| 0.952 | 1.516 | 1.524 | 1.453 | 1.487 |
| 0.849 | 1.557 | 1.568 | 1.500 | 1.496 |
| 0.760 | 1.602 | 1.619 | 1.551 | 1.509 |
| 0.619 | 1.711 | 1.745 | 1.678 | 1.605 |
| 0.550 | 1.805 | 1.863 | 1.795 | — |
| 0.500 | 2.000 | 2.019 | — | — |
| 0.386 | — | 3.387 | 7.9 | 2.007 |

The Monte Carlo reproduces Eq. (5) to better than $1\%$ over most of the range
and $3\%$ at worst, the residual positive throughout and consistent with the
finite-$N$ cutoff, which steepens a Hill fit. That simulation shares every
assumption of the derivation except the neglect of fluctuations, so the agreement
establishes Eq. (5) as the exact well-mixed answer and identifies the remaining
gap to the lattice — uniformly downward, widening as $q$ falls — as spatial
correlation rather than a deficiency of the closed form.

Two further checks. Eq. (6) is satisfied by the Monte Carlo to three decimals
($\tau_s=2.324$ against $2\tau_m-1=2.324$ at $q=0.70$), as it must be when both
partners are genuinely drawn from the bulk; on the lattice it fails at large $q$,
where emissions occur only on domain interfaces. And below $q=\tfrac12$ the Monte
Carlo does not reproduce the rate equation's runaway, returning instead a steep
cutoff-limited fit that drifts with $q$ ($2.36$ at $q=0.45$, $3.39$ at
$q=0.386$) — the signature of no scale-free solution, rather than of a power law
with large exponent. The lattice nonetheless sustains a clean $\tau_m\approx2$
there, because the annihilation-dominated regime is controlled by two-species
annihilation $A+B\to\emptyset$, whose upper critical dimension is $d_c=4$; a
controlled treatment requires a Doi–Peliti field theory for coagulation,
annihilation and injection.

## Remark on the RG reading

Equation (4) states that the coefficient of the marginal operator vanishes, the
analogue of a beta function for $\tau$ set to zero. This is a correspondence
rather than a derivation: no coarse-graining transformation is iterated, and the
structure is that of a solvability condition removing a secular term. A
field-theoretic treatment should yield a genuine $\beta(\tau)$ whose zero
reproduces Eq. (4).
