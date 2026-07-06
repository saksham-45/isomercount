# Linear Snake Growth — Exact Closed Form

Deprecated note about growth rate — supersedes the earlier approximation in
this file.

## The result

The linear-snake count $N_L$ satisfies (for $L\ge 2$) the order-12 integer
linear recurrence with characteristic polynomial

$$P(x) = (x-3)(x-1)(x^2-3)(x^2-5x-5)(x^2-3x+1)(x^4-5x^2-5).$$

Equivalently, for $L\ge 14$,

$$N_L = 12N_{L-1} - 38N_{L-2} - 38N_{L-3} + 370N_{L-4} - 394N_{L-5}
- 556N_{L-6} + 1160N_{L-7} - 690N_{L-8} + 370N_{L-9} + 330N_{L-10}
- 750N_{L-11} + 225N_{L-12}.$$

Asymptotically

$$N_L \sim A \cdot \mu^L, \quad
A = \frac{5+2\sqrt 5}{20}, \quad
\mu = \frac{5+3\sqrt 5}{2}.$$

Numerically $A \approx 0.4736067977499789696…$ and $\mu \approx 5.8541019662496845…$

## Correction to the previous version

An earlier draft of this document said $\mu = 3 + \sqrt 5 \approx 5.236$.
That is **wrong**. The correct dominant eigenvalue is $\mu = (5+3\sqrt 5)/2$,
the larger root of $x^2 - 5x - 5$. It happens that $\mu$ can be written as
$\varphi^4 - \varphi^2 = \varphi^2 + \varphi = \varphi^2(1+1/\varphi)$ etc.,
but not as $3+\sqrt 5$ (numerically distinct).

The amplitude $A$ was also given as a float; the exact value is $(5+2\sqrt 5)/20$.

## Implementation

Use [src/analysis/closed_form.py](src/analysis/closed_form.py):

```python
from src.analysis.closed_form import linear_snake_count, linear_snake_asymptotic

linear_snake_count(100)                 # exact 77-digit integer
linear_snake_asymptotic(100, digits=50)  # 'A·μ^100' to 50 decimals (mpmath)
```

For the full derivation and the cyclic-graph closed form, see
[CLOSED_FORM.md](CLOSED_FORM.md).
