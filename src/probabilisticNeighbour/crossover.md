# The two-slope emission spectrum in 2D

## Rules

A square lattice, periodic, $L^2$ sites. Each occupied site carries a signed integer
$\pm m$; empty sites carry $0$. A fixed number $N=\rho L^2$ of spots is maintained.

One update:

1. Pick an occupied site uniformly, move its charge to a random nearest neighbour.
2. If the target is empty, the spot moves there.
3. If the target holds the same sign, the two merge: $i,j\to i+j$.
4. If the target holds the opposite sign, they cancel partially, leaving $|i-j|$ and
   **emitting** $s=\min(i,j)$.
5. If the spot count has fallen below $N$, inject one $+1$ and one $-1$. With
   probability $p$ each is placed next to an existing spot of its own sign; otherwise
   both are placed at independent random empty sites.

$p$ is the only parameter beyond $\rho$ and $L$. It interpolates between uncorrelated
injection ($p=0$) and co-localised injection ($p=1$).

## Symbols

| | |
|---|---|
| $n(m)\sim m^{-\tau_m}$ | steady-state spot-mass distribution |
| $q$ | fraction of collisions that are same-sign (coagulations) |
| $P_\text{coll}(s)\sim s^{-\alpha}$ | distribution of $\min(i,j)$ over **all** collisions |
| $P_\text{emis}(s)$ | distribution of $s$ over annihilations only, i.e. the emissions |
| $\tau_s^<,\ \tau_s^>$ | emission slopes below and above the crossover |
| $s^*$ | crossover mass |
| $\varphi,\ \psi$ | slopes of $1-q(s)$ below and above $s^*$ |
| $\xi(p)$ | domain size |
| $d_f(\epsilon)$ | box dimension of the set where $+$ and $-$ are adjacent |

## Relations

The emission spectrum is the collision spectrum filtered by outcome:

$$P_\text{emis}(s) = P_\text{coll}(s)\,[1-q(s)] . \tag{1}$$

$\min$ of two independent draws from $n(m)$ gives $n(s)\int_s^\infty n$, so

$$\alpha = 2\tau_m - 1 . \tag{2}$$

$1-q(s)$ is not a scaling function of $s/s^*$: it is a sum of two channels with
independent amplitudes,

$$1-q(s) = A\,s^{-\varphi} + B\,s^{\psi} . \tag{3}$$

Substituting (3) into (1) gives the two slopes and their difference:

$$\tau_s^< = \alpha + \varphi, \qquad
  \tau_s^> = \alpha - \psi, \qquad
  \tau_s^< - \tau_s^> = \varphi + \psi . \tag{4}$$

Setting $\mathrm{d}(1-q)/\mathrm{d}s=0$ locates the crossover and fixes its depth:

$$s^{*\,(\varphi+\psi)} = \frac{\varphi A}{\psi B}, \qquad
  (1-q)_{\min} = B\,s^{*\psi}\left(1+\frac{\psi}{\varphi}\right) . \tag{5}$$

With $\alpha=2$, the mean emission $\langle s\rangle \propto \int s^{-1}[1-q(s)]\,\mathrm{d}s$
is logarithmically divergent at $\psi=0$ and cutoff-dominated for $\psi>0$: the
dissipated mass leaves through the rarest, largest events.

## Values

Equations (1), (2), (4), (5) are derived. No exponent *value* below is derived; all are
measured at $\rho=0.2$, $L=512$–$1024$.

| | | |
|---|---|---|
| $\tau_m$ | $1.465$–$1.52$ | drifts slightly with $p$; $1.465$ at $p=1$ |
| $\alpha$ | $2.01$–$2.03$ | flat under $s^{2\tau_m-1}$ compensation for all $p$ |
| $\psi$ | $0.26\pm0.02$ | independent of $p$ for $p\ge0.75$ |
| $\tau_s^>$ | $1.75$–$1.77$ | $=\alpha-\psi$ |
| $A$ | $\approx0.6$ | independent of $p$ |
| $B$ | $\approx0.011$ | roughly constant for $p\ge0.75$ |

$\varphi$ runs with $p$: $0.375,\,0.550,\,0.755,\,1.008,\,1.207$ at
$p=0.70,\,0.75,\,0.80,\,0.85,\,0.90$.

$\psi$ is a genuine power, not a logarithm read over too short a range: fitting the
rising branch to $A s^{\psi}$ against $A(\ln s)^{y}$ gives weighted rms residuals of
$0.024$–$0.038$ versus $0.063$–$0.092$ across $p=0.75$–$0.90$ and both $L$.

$\tau_m$ is *not* resolved. Over $\sim3.5$ decades a pure power $m^{-1.47}$ and
$m^{-3/2}(\ln m)^{0.2}$ fit equally well (rms $0.011$–$0.040$ against $0.012$–$0.033$),
with neither winning consistently. The measured value is consistent with $3/2$ but does
not establish it, and the small deficit cannot be attributed to a log correction on this
evidence.

## Geometry

$d_f(\epsilon)$ is measured by box-counting the midpoints of Voronoi bonds joining
opposite-sign spots. Counting boxes that merely *contain* both signs instead measures
density at small $\epsilon$ and returns negative dimensions.

$d_f = 1$ for $\epsilon<\xi(p)$ and $2$ for $\epsilon>\xi(p)$: domains are compact, so
reactions are confined to their boundaries at short range and fill space at long range.
The dimension crosses between two integers rather than taking a fractional value, which
is why (4) yields two clean slopes and not one that varies continuously.

$\xi \approx 4,\,16,\,32$ at $p=0.7,\,0.8,\,0.9$, independent of $L$. At $p=1$,
$d_f=1$ at every $\epsilon$ up to $L/8$ and $\xi$ grows with $L$; the crossover then
disappears and the spectrum is a single power law. Read $\xi$ as the departure from the
$p=0$ baseline, which is itself below $2$ at small $\epsilon$.

$s^*$ is identical at $L=256,512,1024$ for $0.70\le p\le0.90$, so the crossover survives
$L\to\infty$. Below $p\approx0.65$ it is cut off by $\xi\to L$, above $p\approx0.95$ by
the minimum leaving the measurable range.

## Open

No 2D exponent value is derived. Equation (2) fixes $\alpha$ *given* $\tau_m$, and (4)
fixes the emission slopes *given* $\varphi,\psi$; every chain terminates on a number
taken from simulation.

$\tau_m$ in 2D is the load-bearing unknown. At $p=1$ annihilation is absent and the
process is constant-kernel aggregation with injection, $\tau_m=3/2$. At intermediate $p$
the mean-field line $\cos\pi\tau_m=(1-q)/q$ predicts $1.61$ at $p=0.8$ against a measured
$1.47$, so the 2D value is neither the mean-field one nor demonstrably $3/2$.

$\psi\approx1/4$ would give $\tau_s^>=7/4$; this is numerology. $1/4=1/2d$ in $d=2$
predicts $\psi=1/6$, $\tau_s^>=11/6$ in three dimensions, which separates it from a
dimension-independent $\psi$.

$\xi$ grows with $p$ while $s^*$ falls, so $s^*$ is not a geometric length: by (5) the
geometry fixes $\psi$ and the injection rate fixes where the two channels of (3) cross.
No mass–length map bridges them — spot mass is uncorrelated with distance to a domain
wall (median mass $2$ at every distance, largest spots as often at a wall as deep
inside).

Two routes to $\psi$ that fail. Depletion: a spot of mass $m$ has swept a radius
$R\sim m^{1/d}$ clear of like-sign partners, so coagulation at contact is suppressed by
the pair correlation there. In $d=2$ that suppression is $1/\ln R$, giving
$1-q\sim\ln s$, and over $s=10^2$–$10^4$ at $p=0.8$ it predicts a rise of $1.9$ against
a measured $3.1$; the power form predicts $2.9$. Balance: equating injected and
annihilated mass gives $q/(1-q)=2(\langle s\rangle-f_c)$ with $f_c$ the fraction of
annihilations that cancel exactly, which misses the measured ratios by $30$–$50\%$ in
both directions, because injection restores the spot count in steps of two and $N$ is
fixed only on average.
