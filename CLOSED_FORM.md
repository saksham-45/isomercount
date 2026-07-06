# Closed-form counts for cyclic graphs and linear snakes

This report documents the exact closed forms for both enumeration problems
in this repository, and where the previous documentation was wrong.

Verification: both formulas were checked against 100 exact terms generated
by the Burnside enumerator (all-integer arithmetic, 0 mismatches), and the
linear-snake algebraic closed form was checked against those same terms
using 300-digit mpmath (agreement to ~223 digits — pure floating roundoff).

Production code: [src/analysis/closed_form.py](src/analysis/closed_form.py).

---

## 1. Cyclic graphs with adjacency (A053656)

### Formula

$$a(n) \;=\; \frac{1}{2n}\!\left[\;\sum_{d\mid n}\varphi\!\left(\frac{n}{d}\right)\!2^{d} \;+\; [n\text{ even}]\cdot n\cdot 2^{n/2-1}\right].$$

Equivalently, split by parity:

$$a(n)=
\begin{cases}
\displaystyle\frac{1}{2n}\sum_{d\mid n}\varphi(n/d)\,2^{d} & n \text{ odd}\\[6pt]
\displaystyle\frac{1}{2n}\!\left(\sum_{d\mid n}\varphi(n/d)\,2^{d}\right) + 2^{n/2-2} & n \text{ even, } n\ge 4\\[6pt]
2 & n = 2
\end{cases}$$

This is **OEIS A053656**. Because the count is a Burnside average over the
dihedral group $D_n$, it is intrinsically a divisor sum in $n$ — there is no
elementary closed form in $n$ alone. This is the simplest form possible.

### Why the adjacency constraint changes nothing

The 6-state transfer matrix

$$M = \begin{pmatrix}0&0&1&0&0&1\\0&0&0&1&0&0\\0&0&1&0&0&1\\1&0&0&1&0&0\\0&0&1&0&0&0\\1&0&0&1&0&0\end{pmatrix}$$

(states $=(\text{last vertex}, \text{current edge})\in\{-2,0,+2\}\times\{0,1\}$)
has characteristic polynomial

$$\chi_M(x) = x^{5}(x-2).$$

So $\text{tr}(M^n) = 2^n$ exactly — the same value you get with **no** transfer
matrix at all (all $2^n$ binary edge sequences). The adjacency rule redirects
the trace through the extra states but does not remove any closed walks.

The reflection-fixed count $\hat R(n) = 2^{n/2}$ (even $n$) and $0$ (odd $n$).

### Asymptotic

$$a(n) \sim \frac{2^{n}}{2n}, \qquad n\to\infty.$$

---

## 2. Linear snake polyominoes

### Sequence

$$N_1,\,N_2,\,\dots = 4,\;21,\;109,\;586,\;3326,\;19209,\;111871,\;653758,\;\dots$$

### Order-12 integer recurrence  (L ≥ 14)

$$N_L = 12N_{L-1} - 38N_{L-2} - 38N_{L-3} + 370N_{L-4} - 394N_{L-5}
- 556N_{L-6} + 1160N_{L-7} - 690N_{L-8} + 370N_{L-9}
+ 330N_{L-10} - 750N_{L-11} + 225N_{L-12}.$$

Seed values $N_2, \dots, N_{13}$ from the Burnside enumeration. $N_1 = 4$ is
an **exception** because at $L=1$ a single square is fixed by all 8 dihedral
symmetries, whereas at $L \ge 2$ only 4 of them (identity, $x$-reflection,
path-reverse composed with $180°$ rotation, path-reverse with $y$-reflection)
preserve the shape.

This recurrence was verified against exact ground truth for $L = 14, \dots, 100$
with zero mismatches.

### Characteristic polynomial

The minimal polynomial of the sequence for $L \ge 2$ is

$$P(x) = (x-3)(x-1)(x^{2}-3)(x^{2}-5x-5)(x^{2}-3x+1)(x^{4}-5x^{2}-5).$$

Degree 12. Roots:

| factor | roots | numeric |
|---|---|---|
| $x-3$ | $3$ | $3$ |
| $x-1$ | $1$ | $1$ |
| $x^{2}-3$ | $\pm\sqrt{3}$ | $\pm 1.732$ |
| $x^{2}-5x-5$ | $\mu_{1,2} = (5\pm 3\sqrt 5)/2$ | $5.854,\ -0.854$ |
| $x^{2}-3x+1$ | $\nu_{1,2} = (3\pm\sqrt 5)/2 = \varphi^{2},\ 1/\varphi^{2}$ | $2.618,\ 0.382$ |
| $x^{4}-5x^{2}-5$ | $\pm\sqrt{\mu_1},\ \pm i\sqrt{-\mu_2}$ | $\pm 2.420,\ \pm 0.924\,i$ |

### Algebraic closed form (L ≥ 2)

$$
N_L =\; -\tfrac{1}{4}\bigl(3^{L}+1+(\sqrt3)^{L}+(-\sqrt3)^{L}\bigr)
\;+\;\tfrac{5+2\sqrt5}{20}\,\mu_{1}^{L}\;+\;\tfrac{5-2\sqrt5}{20}\,\mu_{2}^{L}
\;+\;\tfrac{5+2\sqrt5}{20}\,\nu_{1}^{L}\;+\;\tfrac{5-2\sqrt5}{20}\,\nu_{2}^{L}
\;+\; G(L)
$$

where $G(L)$ is the contribution of the four roots of $x^{4}-5x^{2}-5$
(the residues are exact algebraic numbers but not pretty; the numerical
values are $A[\rho_1]\approx 0.9113$, $A[\rho_2]\approx 0.0359$, and
$A[\rho_{3,4}]\approx 0.02639 \pm 0.06386\,i$).

**Cleaner statement**: the $(x^2-3)$ term simplifies for even $L$ to
$-\tfrac{1}{4}\bigl(3^{L/2}+1\bigr)^{2}$, and for odd $L$ vanishes to
$-\tfrac{1}{4}(3^{L}+1)$.

### Asymptotic (the "clean" part of the closed form)

$$\boxed{\; N_L \;\sim\; A\cdot\mu^{L}, \qquad A=\tfrac{5+2\sqrt 5}{20},\quad \mu=\tfrac{5+3\sqrt 5}{2}\;}$$

Numerically:
$$A = 0.47360679774997896964091736687312762354\ldots$$
$$\mu = 5.85410196624968454461376626416967812779\ldots$$

Convergence is fast: $N_L / (A\mu^L) = 1 + O(\mu_2/\mu_1)^L \approx 1 + O(0.146^L)$.

### What the previous docs said, and why it was wrong

`rdme.md` and `testt.py` claimed
$$\mu = 3+\sqrt 5 \approx 5.236, \qquad A \approx 0.4736068.$$

* The $\mu$ value **5.236** is off by 12% from the true $\mu \approx 5.854$.
  It appears they identified the numeric value from a floating-point fit
  (which gave $\approx 5.854$) but then symbolically guessed the wrong
  algebraic form. The correct $\mu$ is the larger root of $x^2-5x-5$, not of
  $x^2-6x+4$ (which is what $3+\sqrt 5$ is a root of).
* The numeric $A \approx 0.4736068$ is right, but only to 6 decimals; the
  exact value is $A = (5+2\sqrt 5)/20 = 0.4736067977499789696\ldots$

`final_verification.py` also proposed a degree-15 polynomial
$$(x^2-5x-5)(x^2-3x+1)(x^4-5x^2-5)(x^2-2x-3)(x^4-2x^2-3)$$
which contains the correct three factors but adds extra spurious factors
$(x^{2}-2x-3) = (x-3)(x+1)$ and $(x^{4}-2x^{2}-3) = (x^{2}-3)(x^{2}+1)$.
The true minimal polynomial has $(x-3)$ and $(x^{2}-3)$ but **not** $(x+1)$
or $(x^{2}+1)$. The residual $324$ in that verification script is exactly the
constant contribution from those bogus factors — a hint the polynomial was
over-fitted.

---

## 3. What each result buys you

| Quantity | Best exact form | Cost | Range |
|---|---|---|---|
| $a(n)$ cyclic-adjacency | $\bigl(\text{rot sum} + n\cdot 2^{n/2-1}[n \text{ even}]\bigr)/(2n)$ | $O(\sqrt n \log n)$ arithmetic ops | any $n$ |
| $N_L$ linear snake | order-12 integer recurrence with the 12 coefficients above | $O(L)$ big-int adds | any $L$ |
| $N_L$ asymptotic | $\dfrac{5+2\sqrt 5}{20}\cdot\Bigl(\dfrac{5+3\sqrt 5}{2}\Bigr)^{L}$ | $O(1)$ mpf ops | approximate |

The recurrence is what's actually useful for computing big values exactly.
The algebraic closed form is what's useful for *understanding* the growth —
in particular the golden-ratio connection: the sub-dominant real eigenvalues
$\nu_{1,2} = (3\pm\sqrt 5)/2$ are exactly $\varphi^{\pm 2}$, and the dominant
$\mu_1 = (5+3\sqrt 5)/2$ can also be written as $\varphi^{2}\cdot\sqrt 5 \cdot (\text{small})$
or more usefully in terms of the Lucas numbers:

$$\mu_{1}^{L} + \mu_{2}^{L} \;=\; 5\bigl(\mu_{1}^{L-1}+\mu_{2}^{L-1}\bigr) + 5\bigl(\mu_{1}^{L-2}+\mu_{2}^{L-2}\bigr),$$

a Lucas-like sequence with characteristic $x^2 - 5x - 5$, initial pair $(2, 5)$.
