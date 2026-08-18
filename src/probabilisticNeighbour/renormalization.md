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
would give $q=\tfrac12$. The measured value runs monotonically from $q=0.50$ at
$p=0$ to $q=0.95$ at $p=1$ — at $p=0$ the two injected spots land on independent
random sites, so the surface is genuinely uncorrelated and the sweep starts
exactly at $\tfrac12$ and never enters $q<\tfrac12$.

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

Notice how little we actually used: only that $n_m$ falls fast enough for the
small-partner terms to dominate. The cancellation would have happened for *any*
steeply falling $n_m$, power law or not. That is a warning sign — an equation
that is satisfied by everything is telling you nothing, and we will have to look
harder.

This is what "line of fixed points" means here. The RG flow in mass space has a
whole one-parameter family of scale-invariant solutions, and $\tau$ is a
**marginal** direction: rescaling does not push it anywhere, so the linearized
RG has nothing to say about which member is realised.

But "nothing at leading order" is not "nothing". A marginal direction is exactly
the case where you have to go to the *next* order. Section 4 does that, and the
next order turns out to be computable in closed form.

---

## 4. Getting $\tau_m(q)$: keep the terms we threw away

Section 3 replaced both sums by $n_m N$. Those replacements are
*approximations*. Approximations have errors, and the errors did not cancel.
This section computes them.

The whole calculation is one idea repeated: **rescale the integration variable
by $\mu$, so that all the $\mu$-dependence sits either in a prefactor or in a
cutoff, then expand in the cutoff.** Nothing beyond first-year calculus is
needed, plus two standard integrals which are quoted where used.

### 4.1 Setting up

First go continuum. For $\mu\gg1$ the summand varies slowly from one integer to
the next, so sums become integrals with a relative error $O(1/\mu)$. We will be
keeping terms of relative size $\mu^{1-\tau}$, and since $\tau<2$ we have
$\mu^{1-\tau}\gg\mu^{-1}$, so the replacement is safe.

Take $n(m) = A m^{-\tau}$ for $m\ge1$ and $n(m)=0$ below — the lower cutoff is
the injection scale. Define the two collision integrals

$$
C(\mu) = \int_1^{\mu-1} n(i)\,n(\mu-i)\,di ,
\qquad
K(\mu) = \int_1^{\infty} n(\mu+j)\,n(j)\,dj ,
$$

so $(\star)$ becomes

$$ \frac{q}{2}\,C(\mu) + (1-q)\,K(\mu) = n(\mu)\,N . \tag{$\star'$} $$

We will need $N$ explicitly:

$$
N = \int_1^\infty A m^{-\tau} dm
= A\left[\frac{m^{1-\tau}}{1-\tau}\right]_1^\infty
= \frac{A}{\tau-1} ,
$$

which is finite because $\tau>1$. Note it is set by the *bottom* of the
distribution — there are far more monomers than big spots — so $N$ is just a
number, with no $\mu$ in it. Keep that in mind; it is why $N$ can appear on both
sides without spoiling anything.

Throughout, $1<\tau<2$. We will see at the end that this is exactly the range
where everything converges.

### 4.2 The coagulation integral, step by step

**Step 1. Substitute.** Put $i = \mu x$, so $di = \mu\,dx$, and the limits
$i=1$, $i=\mu-1$ become $x=1/\mu$, $x = 1-1/\mu$. Write $\delta \equiv 1/\mu$.
The two factors are

$$
n(i) = A(\mu x)^{-\tau} = A\mu^{-\tau}x^{-\tau} ,
\qquad
n(\mu-i) = A(\mu-\mu x)^{-\tau} = A\mu^{-\tau}(1-x)^{-\tau} .
$$

Multiplying, and collecting the $\mu$'s ($\mu^{-\tau}\cdot\mu^{-\tau}$ from the
two factors, $\mu$ from $di$):

$$
C(\mu) = A^2\,\mu^{1-2\tau}\int_{\delta}^{1-\delta} x^{-\tau}(1-x)^{-\tau}\,dx .
$$

**This is the whole trick.** Every trace of $\mu$ is now in the prefactor
$\mu^{1-2\tau}$ or in the cutoff $\delta = 1/\mu$. The integral itself is a pure
number once $\delta$ is fixed.

**Step 2. Locate the trouble.** Call the integral $I(\delta)$. Near $x=0$ the
factor $(1-x)^{-\tau}\to1$, so the integrand behaves like $x^{-\tau}$, and

$$ \int_\delta x^{-\tau}\,dx \;\sim\; \frac{\delta^{1-\tau}}{\tau-1} $$

blows up as $\delta\to0$, because $\tau>1$. The same happens at $x=1$ by the
symmetry $x\leftrightarrow1-x$. So $I(\delta)$ does not have a finite
$\delta\to0$ limit — but the divergence is a single, known power, and everything
else is finite. Let us separate the two.

**Step 3. Add and subtract the singular pieces.** Write the integrand as

$$
\underbrace{x^{-\tau}(1-x)^{-\tau} - x^{-\tau} - (1-x)^{-\tau}}_{\text{call this } g(x)}
\;+\; x^{-\tau} \;+\; (1-x)^{-\tau} .
$$

This is an identity — we have added and subtracted the same two terms. The point
is that $g(x)$ is *harmless at both ends*. Check $x\to0$: expand
$(1-x)^{-\tau} = 1+\tau x+O(x^2)$, so

$$
g(x) = x^{-\tau}\big[(1-x)^{-\tau}-1\big] - (1-x)^{-\tau}
\;\approx\; \tau\,x^{1-\tau} - 1 ,
$$

and $\int_0 x^{1-\tau}dx$ converges because $\tau<2$. By symmetry the same holds
at $x=1$. So $\int_0^1 g(x)\,dx$ is an ordinary finite number; call it $G(\tau)$.

**Step 4. Do the two elementary integrals.** These we can just evaluate:

$$
\int_\delta^{1-\delta} x^{-\tau}dx
= \left[\frac{x^{1-\tau}}{1-\tau}\right]_\delta^{1-\delta}
= \frac{\delta^{1-\tau} - (1-\delta)^{1-\tau}}{\tau-1}
\;\xrightarrow[\delta\to0]{}\; \frac{\delta^{1-\tau}}{\tau-1} - \frac{1}{\tau-1} .
$$

The $(1-x)^{-\tau}$ integral gives exactly the same by symmetry. Adding both:

$$
\frac{2\,\delta^{1-\tau}}{\tau-1} - \frac{2}{\tau-1} .
$$

**Step 5. Assemble.** Putting Steps 3 and 4 together,

$$
I(\delta) = \frac{2\,\delta^{1-\tau}}{\tau-1}
\;+\; \underbrace{\left[G(\tau) - \frac{2}{\tau-1}\right]}_{\equiv\, B_c(\tau)}
\;+\; O(\delta^{2-\tau}) .
$$

Three pieces, in decreasing order of size as $\delta\to0$: one that blows up like
$\delta^{1-\tau}$, one constant, one that vanishes like $\delta^{2-\tau}$. We
keep the first two and drop the third.

**Step 6. Translate back to $\mu$.** Substitute $\delta = 1/\mu$ into the
divergent piece:

$$
A^2\mu^{1-2\tau}\cdot\frac{2\mu^{\tau-1}}{\tau-1}
= \frac{2A^2}{\tau-1}\,\mu^{1-2\tau+\tau-1}
= \frac{2A^2}{\tau-1}\,\mu^{-\tau}
= 2\,\big(A\mu^{-\tau}\big)\cdot\frac{A}{\tau-1}
= 2\,n(\mu)\,N ,
$$

using $N = A/(\tau-1)$ from §4.1. So

$$ \boxed{\;C(\mu) = 2\,n(\mu)\,N \;+\; A^2\, B_c(\tau)\,\mu^{1-2\tau}\;} $$

The first term is *exactly* the "small partner dominates" estimate of Section 3
— which makes sense, since the divergence came from $x\approx0$, i.e. $i$ of
order $1$, i.e. one partner small. The second term is what Section 3 threw away.
Relative to the first it is a factor

$$ \frac{\mu^{1-2\tau}}{\mu^{-\tau}} = \mu^{1-\tau} \;\to\; 0 , $$

so it genuinely is subleading. Subleading is not zero.

**Step 7. Name $B_c$.** The subtraction in Step 3 is all we need — $B_c$ is
defined and finite. But it has a closed form worth knowing. If $\tau$ were less
than $1$ no subtraction would be necessary, and the integral would be the
standard Beta integral

$$
\int_0^1 x^{a-1}(1-x)^{b-1}dx = B(a,b) = \frac{\Gamma(a)\Gamma(b)}{\Gamma(a+b)} ,
$$

with $a=b=1-\tau$. Our subtracted $G(\tau)-2/(\tau-1)$ agrees with that whenever
both make sense, and it is a smooth function of $\tau$ on the whole range
$\tau<2$. Two smooth functions agreeing on an interval are the same function, so
the Beta formula carries over — the $\Gamma$'s handle negative arguments by
themselves:

$$ B_c(\tau) = \frac{\Gamma(1-\tau)^2}{\Gamma(2-2\tau)} . $$

(This is what "analytic continuation" means in practice: the subtraction *is*
the continuation. Verified numerically against direct quadrature of $G(\tau)$ to
six decimal places.)

### 4.3 The annihilation integral

Same technique, and slightly easier. Put $j = \mu y$, $dj = \mu\,dy$, with
$n(\mu+j) = A\mu^{-\tau}(1+y)^{-\tau}$ and $n(j)=A\mu^{-\tau}y^{-\tau}$:

$$
K(\mu) = A^2\mu^{1-2\tau}\int_\delta^\infty y^{-\tau}(1+y)^{-\tau}\,dy .
$$

Now there is only **one** bad end. At $y\to0$ the integrand is $\sim y^{-\tau}$,
divergent as before. At $y\to\infty$ it is $\sim y^{-2\tau}$, which converges
since $\tau>1$. So we subtract one singular piece instead of two:

$$
y^{-\tau}(1+y)^{-\tau} = \underbrace{y^{-\tau}\big[(1+y)^{-\tau}-1\big]}_{\text{finite at both ends}} + y^{-\tau} .
$$

The first bracket goes like $-\tau y^{1-\tau}$ at small $y$ (integrable, $\tau<2$)
and like $-y^{-\tau}$ at large $y$ (integrable, $\tau>1$), so its integral over
$(0,\infty)$ is a finite number. The second is elementary:

$$
\int_\delta^\infty y^{-\tau}dy = \frac{\delta^{1-\tau}}{\tau-1} .
$$

Hence

$$
\int_\delta^\infty y^{-\tau}(1+y)^{-\tau}dy
= \frac{\delta^{1-\tau}}{\tau-1} + B_a(\tau) + O(\delta^{2-\tau}) ,
\qquad
B_a(\tau) = \int_0^\infty y^{-\tau}\big[(1+y)^{-\tau}-1\big]dy .
$$

Translating back exactly as in Step 6 — one copy of the divergence instead of
two, so one copy of $n(\mu)N$ instead of two:

$$ \boxed{\;K(\mu) = n(\mu)\,N \;+\; A^2\, B_a(\tau)\,\mu^{1-2\tau}\;} $$

And by the same argument as Step 7, using the other standard Beta representation

$$
\int_0^\infty y^{a-1}(1+y)^{-(a+b)}dy = B(a,b) ,
$$

read off $a-1 = -\tau$ and $a+b = \tau$, so $a = 1-\tau$ and $b = 2\tau-1$:

$$ B_a(\tau) = \frac{\Gamma(1-\tau)\Gamma(2\tau-1)}{\Gamma(\tau)} . $$

### 4.4 The leftovers have to cancel by themselves

Substitute both boxed results into $(\star')$:

$$
\frac{q}{2}\Big[2n(\mu)N + A^2B_c\,\mu^{1-2\tau}\Big]
+ (1-q)\Big[n(\mu)N + A^2B_a\,\mu^{1-2\tau}\Big]
= n(\mu)N .
$$

Collect the two kinds of term. The $n(\mu)N$ terms give

$$
\frac{q}{2}\cdot2\,n(\mu)N + (1-q)\,n(\mu)N = \big[q + (1-q)\big]n(\mu)N = n(\mu)N ,
$$

which cancels the right-hand side exactly — that is Section 3's marginality,
reappearing, and this time we can see it is *exact* rather than approximate.
What remains is

$$
A^2\mu^{1-2\tau}\left[\frac{q}{2}B_c(\tau) + (1-q)B_a(\tau)\right] = 0 .
$$

Since $A\ne0$ and $\mu^{1-2\tau}\ne0$, the bracket must vanish:

$$ \frac{q}{2}\,B_c(\tau) + (1-q)\,B_a(\tau) = 0 . \tag{$\dagger$} $$

**Why nothing else can rescue a wrong $\tau$.** This is the crux, so it is worth
being explicit about what could have gone wrong.

1. *Could the leftover cancel against the leading term?* No. The leading terms
   carry $\mu^{-\tau}$ and the leftovers carry $\mu^{1-2\tau}$. Two different
   powers of $\mu$ cannot cancel each other at every $\mu$.
2. *Could we absorb it by choosing the amplitude $A$?* No. $A^2$ multiplies both
   leftovers, so it divides straight out of $(\dagger)$. This is also why the
   answer will depend on $q$ alone — not on the density, not on the injection
   rate.
3. *Could we patch it by adding another power to $n$,* say
   $n = Am^{-\tau} + A'm^{-\tau'}$ with $\tau'>\tau$? No. Recall the remark at
   the end of Section 3: the leading cancellation holds for *any* steeply falling
   $n$, so the $A'$ piece cancels up there along with everything else and
   contributes nothing at order $\mu^{-\tau}$. Down at the subleading order it
   contributes cross terms $\propto\mu^{1-\tau-\tau'}$ and $\mu^{1-2\tau'}$, and
   since $\tau'>\tau$ both are *smaller* than $\mu^{1-2\tau}$. They cannot cancel
   the term we need cancelled.

Every other dial is either fixed or divides out. The only free quantity left is
$\tau$, and $(\dagger)$ is one equation. One equation, one unknown.

### 4.5 Solving $(\dagger)$

Rearrange, then substitute the two $\Gamma$ expressions:

$$
\frac{q}{2}B_c = -(1-q)B_a
\quad\Longrightarrow\quad
\frac{q}{1-q} = -\frac{2B_a}{B_c}
= -2\,\frac{\Gamma(1-\tau)\Gamma(2\tau-1)}{\Gamma(\tau)}
\cdot\frac{\Gamma(2-2\tau)}{\Gamma(1-\tau)^2} .
$$

One factor of $\Gamma(1-\tau)$ cancels:

$$
\frac{q}{1-q} = -2\,\frac{\Gamma(2\tau-1)\,\Gamma(2-2\tau)}{\Gamma(\tau)\,\Gamma(1-\tau)} .
$$

Now use the reflection formula $\Gamma(z)\Gamma(1-z) = \pi/\sin(\pi z)$ on top
and bottom.

*Denominator*, with $z=\tau$:

$$ \Gamma(\tau)\Gamma(1-\tau) = \frac{\pi}{\sin\pi\tau} . $$

*Numerator*, with $z = 2\tau-1$. Then $1-z = 2-2\tau$, so the two $\Gamma$'s are
exactly of the reflection form:

$$ \Gamma(2\tau-1)\Gamma(2-2\tau) = \frac{\pi}{\sin\!\big(\pi(2\tau-1)\big)} . $$

Simplify that sine using $\sin(\theta-\pi) = \sin\theta\cos\pi - \cos\theta\sin\pi = -\sin\theta$:

$$ \sin\!\big(\pi(2\tau-1)\big) = \sin(2\pi\tau - \pi) = -\sin 2\pi\tau . $$

So the numerator is $-\pi/\sin2\pi\tau$, and

$$
\frac{q}{1-q}
= -2\cdot\frac{-\pi/\sin2\pi\tau}{\pi/\sin\pi\tau}
= \frac{2\sin\pi\tau}{\sin2\pi\tau} .
$$

Finally the double-angle identity $\sin2\pi\tau = 2\sin\pi\tau\cos\pi\tau$:

$$
\frac{q}{1-q} = \frac{2\sin\pi\tau}{2\sin\pi\tau\,\cos\pi\tau} = \frac{1}{\cos\pi\tau} .
$$

Everything cancels down to a bare cosine:

$$ \boxed{\;\cos(\pi\tau_m) = \frac{1-q}{q}\;} $$

**Which root?** $\cos(\pi\tau)=c$ has infinitely many solutions, but the
derivation is only valid for $1<\tau<2$ — below $1$ the density $N$ diverges,
above $2$ the finite parts $B_c,B_a$ diverge. On that interval $\pi\tau$ runs
over $(\pi,2\pi)$, where the cosine increases monotonically from $-1$ to $+1$.
So there is **exactly one** root, and we did not have to choose it: convergence
chose it for us.

To write it explicitly, put $\pi\tau = 2\pi - \theta$ with $\theta\in(0,\pi)$;
then $\cos\pi\tau = \cos\theta$, so $\theta = \arccos\big((1-q)/q\big)$ and

$$ \tau_m(q) = 2 - \frac{1}{\pi}\arccos\!\left(\frac{1-q}{q}\right) . $$

For $q\in[\tfrac12,1]$ the argument of $\arccos$ runs over $[0,1]$, so
$\theta\in[0,\tfrac{\pi}{2}]$ and $\tau_m\in[\tfrac32,2]$.

### 4.6 Does it work?

Check the two ends first:

- $q=1$: $\cos\pi\tau = 0 \Rightarrow \tau_m = \tfrac32$. This is the
  Takayasu/Smoluchowski constant-flux value, and the Appendix rederives it by a
  route that shares no steps with this one — a genuine consistency check.
- $q=\tfrac12$: $\cos\pi\tau = 1 \Rightarrow \tau_m = 2$. This end used to be a
  hand-wave; now it is a corollary.

And in between it is a genuine prediction. The sharpest test is `wellMixed.cpp`,
a Monte Carlo of the *same* process with the lattice thrown away: keep two pools
of signed masses, pick a spot at random, take its partner from the same sign with
probability $q$ and the opposite sign otherwise, coalesce or annihilate, and
refill with $\pm1$ pairs. That is the Takayasu model with two signs and biased
partner selection, and it is precisely the stochastic process whose deterministic
limit is $(\star)$ — it shares every assumption of the derivation **except** the
neglect of fluctuations. Alongside it, `meanField.py` integrates $(\star)$
directly, and the 2D lattice run has space as well:

| $q$ | closed form | well-mixed MC | rate equation | 2D lattice |
|-----|-----|-----|-----|-----|
| 0.952 | 1.516 | 1.524 | 1.453 | 1.487 |
| 0.849 | 1.557 | 1.568 | 1.500 | 1.496 |
| 0.760 | 1.602 | 1.619 | 1.551 | 1.509 |
| 0.700 | 1.641 | 1.662 | 1.596 | — |
| 0.619 | 1.711 | 1.745 | 1.678 | 1.605 |
| 0.550 | 1.805 | 1.863 | 1.795 | — |
| 0.500 | 2.000 | 2.019 | — | — |
| 0.386 | — | 3.387 | 7.9 | 2.007 |

The Monte Carlo lands on the closed form to better than $1\%$ over most of the
range and $3\%$ at worst, always slightly high — the direction in which a
finite-$N$ cutoff pushes a Hill fit. So the closed form is not merely a decent
approximation to the rate equation; it is the exact answer for the well-mixed
process.

That changes how the last column should be read. Since the Monte Carlo and the
derivation differ only by fluctuations and agree anyway, the whole remaining gap
to the lattice — uniformly downward, widening as $q$ falls — is **spatial
correlation**. Not the closed form being sloppy, and not the rate equation being
deterministic.

A second check comes free. The relation $\tau_s = 2\tau_m-1$ of Section 5 assumes
both colliding partners are drawn from the bulk distribution. In the well-mixed
Monte Carlo that assumption is exactly true, and the relation holds to three
decimals ($\tau_s = 2.324$ against $2\tau_m-1 = 2.324$ at $q=0.70$). On the
lattice it fails at large $q$. That confirms the explanation given in Section 7 —
emissions there occur only on domain interfaces — rather than leaving it a
plausible story.

**Why the table stops at $q=\tfrac12$.** For $q<\tfrac12$ we would need
$\cos\pi\tau_m = (1-q)/q > 1$: no solution. Not "a solution we cannot find" —
there is no power law that solves $(\star)$ at all. The rate equation agrees
emphatically, running off to $\tau=7.9$ at $q=0.386$. The Monte Carlo does not
run away, but neither does it find a power law: it returns a steep,
cutoff-limited fit that keeps drifting with $q$ ($2.36$ at $q=0.45$, $3.39$ at
$q=0.386$), which is what "no scale-free solution" looks like in a finite system.

The lattice, of course, does *not* die: it sits at a clean $\tau_m\approx2$ all
the way down to $q=0.386$. That is not a small correction to mean field, it is
mean field being wrong — see Section 8.

### 4.7 So could RG have told us this?

This is worth being precise about, because "the RG gives a line of fixed points"
sounds like a dead end and is not one.

The rescaling argument in Section 3 is the *linearized* RG about the fixed line,
and it says $\tau$ is a **marginal** direction. Marginal directions are, by
definition, the ones the linearized RG cannot resolve — that is what the word
means. When you have one, the position along it is fixed at the next order, and
finding it means computing the first non-vanishing term in the expansion.

That is exactly what §4.2–4.4 is. The "next order" here is the subleading term in
the small-$\delta$ expansion of the collision integrals, and it is clean enough
to do exactly: a single power $\mu^{1-2\tau}$ with a Beta-function coefficient.
Condition $(\dagger)$ is the statement that the coefficient of the marginal
operator vanishes — the analogue of a beta function for $\tau$, set to zero.

One caveat on that last sentence: it is a correspondence, not a derivation. We
never wrote a coarse-graining transformation and iterated it; what we have is a
**solvability condition** on a scaling ansatz, the same structure as removing
secular terms in Poincaré–Lindstedt, where the frequency of a perturbed
oscillator is undetermined at leading order and fixed by demanding the next
order stay bounded. A proper Doi–Peliti treatment should produce a genuine
$\beta(\tau)$ whose zero reproduces $(\dagger)$; that remains to be checked.

So RG does answer the question. It just cannot answer it at the order people
usually stop at.

---

## 5. Emissions inherit the exponent — in mean field

A spot of size $s$ has probability $P(s)\propto s^{-\tau_m}$, so a spot *larger*
than $s$ has the cumulative probability

$$ P(k>s) \;\propto\; s^{1-\tau_m} . $$

An emission is $s=\min(i,j)$ of a colliding pair. To emit more than $s$ you need
**two** things bigger than $s$ to merge, and in mean field the pair is drawn
$\propto n_i n_j$ — two independent draws — so the two probabilities multiply:

$$ P(\text{emission}>s) \;\propto\; s^{1-\tau_m}\cdot s^{1-\tau_m}
= s^{2-2\tau_m} . $$

Differentiating back from cumulative to density,

$$ P(s) \;\propto\; s^{1-2\tau_m}
\qquad\Longrightarrow\qquad
\boxed{\ \tau_s = 2\tau_m - 1\ } . $$

Hence $\tau_s>\tau_m$ always — a large emission is doubly rare. The spot
distribution is primary: as $\tau_m$ runs $2\to\tfrac32$, $\tau_s$ runs $3\to2$.
Combining with Section 4,

$$ \tau_s(q) = 3 - \frac{2}{\pi}\arccos\!\left(\frac{1-q}{q}\right) . $$

### 5.1 Why it fails above $p_c$

The only input was that the two partners are independent draws from $n(m)$. That
survives the mixed phase and dies at the transition: at $L=128$, $\rho=0.2$, the
measured $(\tau_s,\,2\tau_m-1)$ is $(2.60,2.51)$ at $p=0.3$ and $(2.18,2.18)$ at
$p=0.55$, but $(1.79,1.95)$ at $p=1$.

Dumping the mass pair of every annihilation (`collisionPairs.cpp`) shows which
half of the assumption goes. It is not independence but the distribution itself:
$n(m)$ counts each spot once per snapshot, whereas collisions count it once per
*event*, and a condensate pinned at a domain interface eats one arrival after
another. Its share of events is macroscopic and all at a single enormous mass, so
the collision marginal is much fatter than $n(m)$ — $\tau\approx1.13$–$1.31$
against $\tau_m\approx1.47$, at every $L$. Reweighting to that marginal shifts
the predicted exponent by $-0.3$ to $-0.6$; partner correlation pushes back by
$+0.1$ to $+0.3$.

Conditioning confirms the mechanism: restricted to events whose larger partner
exceeds $10^5$, $\min(i,j)$ is simply the arriving bulk spot, and the measured
exponent is $\tau_m$ itself ($1.48$–$1.54$ over $L=128$–$512$). Emissions are
thus a mixture of condensate collisions, carrying $\tau_m$, and ordinary bulk
pairs, carrying $2\tau_m-1$ — and the mixture settles at neither end:

| $L$ | 64 | 128 | 256 | 512 | 1024 |
|---|---|---|---|---|---|
| $\tau_m$ | 1.520 | 1.477 | 1.473 | 1.466 | 1.467 |
| $\tau_s$ | 1.864 | 1.774 | 1.757 | 1.679 | 1.676 |

$\tau_m$ is flat at $1.47$ — a clean plateau over five decades at $L=1024$, with
the early and late halves of the recording window agreeing to $0.003$ —
confirming Takayasu, the small residual below $3/2$ consistent with the
logarithmic corrections expected at $d_c=2$. $\tau_s$ falls and then
**saturates at $\approx1.68$**, strictly between $\tau_m$ and $2\tau_m-1$ and
fixed by neither.

*How* the mixture reaches that value is not resolved here. The condensate is a
single macroscopic object, so the event-ensemble marginal and the condensate's
share of the emitted flux are not self-averaging: single-seed estimates scatter
by $\sim0.2$ in exponent and $\sim0.15$ in flux share, non-monotone in $L$.
Settling that needs several seeds per size at the pair level. The saturation of
$\tau_s$ itself is solid — three seeds per size, and flat to $0.003$ over a
fourfold change in area.

The practical consequence is that $\tau_s$ is not a derived quantity above
$p_c$. Below the transition it is bound to $\tau_m$ by the identity; above it
the two must be measured independently.

---

## 6. The other renormalization: the sign field

Now throw away the masses and keep only $\sigma_i = \mathrm{sign}(s_i) \in
\{+,-,0\}$. This field has a $\mathbb{Z}_2$ symmetry ($+\leftrightarrow-$) and
an order parameter — the coarse-grained magnetization $m$ — that is zero when
mixed and nonzero when segregated. Neighbour injection and same-sign
coalescence align neighbours; random injection and diffusion scramble them. So
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

The contrast with Section 3 is the whole point. There, $\tau$ was marginal, so
every value was a fixed point and the exponent could slide. Here $K$ is
**relevant**: it flows, only one value is fixed, and it has to be tuned. Same
knob $p$, two completely different roles — marginal for the size cascade,
relevant for the spatial order. That is why the exponent varies smoothly right
through $p_c$ while the order parameter switches on.

Note also that this section says nothing about *avalanche* criticality. The
power law does not need $p=p_c$; it is there for every $p$, because it comes
from a driven cascade (drive at $m=1$, dissipation at the cutoff) and not from a
tuned critical point. Widening the box widens the inertial range without
touching $\tau$.

---

## 7. What the simulation says

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
drawn from the bulk distribution at all. The derivation in Section 5 assumed it
was.

The segregation order parameter is flat and near zero up to $p\approx0.6$ and
climbs steeply after, while $\tau$ does not notice $p_c$ at all.

---

## 8. Where this is shaky

- Sections 2–5 are mean field, and Section 4 shows exactly where that bites: the
  rate equation is quantitative for $q\gtrsim0.7$ and useless as $q\to\tfrac12$,
  because two-species annihilation in $d=2$ is far below its upper critical
  dimension. The segregated phase is also emphatically *not* well mixed, which is
  where the emission relation breaks.
- $\tau_m(q)$ is now derived in closed form for the whole line, but only *within*
  mean field — the closed form and the numerical integration agree with each
  other much better than either agrees with the simulation. A controlled
  treatment of the annihilation-dominated end still needs a Doi–Peliti field
  theory for coagulation + annihilation + injection.
- Condition $(\dagger)$ assumes a pure power law all the way from the injection
  scale to infinity. Real runs have a cutoff $m_c$, and $\mu \ll m_c$ is needed
  for the $\mu^{1-2\tau}$ bookkeeping to be the only correction that matters.
  The $q=1$ case is the awkward one: there $(\dagger)$ degenerates to $B_c=0$
  and the surviving solution carries a genuine flux to the cutoff, so the
  cutoff is not a spectator.
- At $p$ near $1$ a few domain-spanning spots form a condensate that sits
  outside the scaling body, so the fitted exponent depends on where you put the
  lower cutoff.
- The block-spin estimate in Section 6 gets the *structure* right (two stable
  phases, one unstable fixed point in between, $K$ relevant) but not the numbers.
  It says nothing about which universality class, and the injection makes this a
  driven, non-equilibrium version of Ising rather than the equilibrium one.

---

## Appendix. The $q=1$ corner by a different route

At $q=1$ there is no annihilation and the model reduces to constant-kernel
coagulation with monomer injection — the Takayasu aggregation model, whose
mean-field exponent is classic. It is a special case of Section 4, but it can
also be got by a conservation argument that never touches $(\dagger)$, and the
two agreeing is a useful check.

With no annihilation, mass enters at $m=1$ and cascades upward until it is
removed at a cutoff (the largest spot the box holds). In steady state the rate at
which mass crosses any intermediate scale $m$ must be the *same at every scale* —
otherwise mass would pile up somewhere.

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

Substitute $i = mx$, $di = m\,dx$ — the same trick as §4.2 — and pull every
power of $m$ out:

$$
\Phi(m) = \frac{A^2}{\tau-1}\; m^{\,(1-\tau)+(1-\tau)+1} \int_0^1 x^{1-\tau}(1-x)^{1-\tau} dx
= \frac{A^2\,B(2-\tau,\,2-\tau)}{\tau-1}\; m^{\,3-2\tau} ,
$$

with $B$ the Beta function (finite as long as $\tau<2$; note that here no
subtraction is needed, because the extra factor of $i$ makes the integrand
integrable at both ends). So

$$ \Phi(m) \;\propto\; m^{\,3-2\tau} . $$

Demanding that the flux not depend on $m$ kills the exponent:

$$ 3-2\tau = 0 \quad\Longrightarrow\quad \tau_m = \tfrac32 \qquad (q=1). $$

Compare with Section 4: setting $q=1$ in $(\dagger)$ leaves $\tfrac12 B_c(\tau)=0$,
i.e. $\Gamma(2-2\tau) = \infty$, i.e. $2-2\tau$ a non-positive integer. The root
in the valid range $1<\tau<2$ is $\tau=\tfrac32$. Same answer, and the two
arguments share nothing: one is a conservation law, the other a solvability
condition.

The reason this route stops at $q=1$ is that annihilation makes mass leak out at
*every* scale, so $\Phi(m)$ is no longer constant and there is no conservation
law to lean on. That is what forced the longer calculation.
