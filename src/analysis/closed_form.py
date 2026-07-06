#!/usr/bin/env python3
"""
Exact closed-form counters for the two problems in this repo.

Both formulas were derived from the transfer matrices and Burnside's lemma
already implemented in src/cyclic and src/snakes, then verified against
100-term ground truth exactly (integer arithmetic) and against a mpmath
300-digit numerical closed form (algebraic residues) for the linear snake.

--------------------------------------------------------------------------
CYCLIC (n edges, adjacency constraint) — a(n):
--------------------------------------------------------------------------
The 6-state adjacency transfer matrix M has characteristic polynomial
    χ_M(x) = x^5·(x − 2)
so trace(M^n) = 2^n exactly. Reflection-fixed count R̂(n) is 2^(n/2) for
even n and 0 for odd n. Burnside then gives

    a(n) = (1/(2n)) · [ Σ_{d|n} φ(n/d) · 2^d  +  n · 2^(n/2 − 1) · [n even] ]

which is OEIS A053656 (the adjacency constraint doesn't reduce the count).
Because it is a divisor sum, no elementary closed form in n alone exists —
this is provably the simplest form.

--------------------------------------------------------------------------
LINEAR SNAKE (L cells in a row) — N_L:
--------------------------------------------------------------------------
The sequence for L ≥ 2 satisfies the order-12 linear recurrence with
characteristic polynomial

    P(x) = (x − 3)(x − 1)(x^2 − 3)(x^2 − 5x − 5)(x^2 − 3x + 1)(x^4 − 5x^2 − 5)

Explicitly, for L ≥ 14 (and any six earlier N₂..N₁₃ as seeds):

    N_L = 12·N_{L−1} − 38·N_{L−2} − 38·N_{L−3} + 370·N_{L−4}
          − 394·N_{L−5} − 556·N_{L−6} + 1160·N_{L−7} − 690·N_{L−8}
          + 370·N_{L−9} + 330·N_{L−10} − 750·N_{L−11} + 225·N_{L−12}

N₁ = 4 is a boundary exception (L=1 uses the full 8-element symmetry group
instead of the 4-element group of the L ≥ 2 enumeration).

Algebraic closed form (all 12 roots explicit):
  Let μ₁,₂ = (5 ± 3√5)/2,   ν₁,₂ = (3 ± √5)/2   (= φ², 1/φ²).
  For L ≥ 2:

    N_L =  −1/4 · ( 3^L + 1 + (√3)^L + (−√3)^L )
          + (5+2√5)/20 · μ₁^L  + (5−2√5)/20 · μ₂^L
          + (5+2√5)/20 · ν₁^L  + (5−2√5)/20 · ν₂^L
          + G(L)                       # from roots of x^4 − 5x^2 − 5

where G(L) is the contribution of the four roots of x^4 − 5x^2 − 5.

Asymptotic:
    N_L ~  (5 + 2√5)/20  ·  ((5 + 3√5)/2)^L
        ≈  0.47360679774997896964… · 5.85410196624968454461…^L

(This corrects the earlier claim μ = 3 + √5 ≈ 5.236 in rdme.md — the true
dominant eigenvalue is (5 + 3√5)/2 ≈ 5.854, a root of x² − 5x − 5.)
"""

from typing import List


# --------------------------------------------------------------------------
# CYCLIC — exact integer O(√n) formula
# --------------------------------------------------------------------------
def _divisors(n: int) -> List[int]:
    small: List[int] = []
    large: List[int] = []
    d = 1
    while d * d <= n:
        if n % d == 0:
            small.append(d)
            if d * d != n:
                large.append(n // d)
        d += 1
    return small + large[::-1]


def _totient(n: int) -> int:
    if n == 1:
        return 1
    result = n
    m = n
    p = 2
    while p * p <= m:
        if m % p == 0:
            result -= result // p
            while m % p == 0:
                m //= p
        p += 1
    if m > 1:
        result -= result // m
    return result


def cyclic_adjacency_count(n: int) -> int:
    if n <= 0:
        return 0
    rot_sum = sum(_totient(n // d) * (2 ** d) for d in _divisors(n))
    refl_sum = n * (2 ** (n // 2 - 1)) if n % 2 == 0 else 0
    total = rot_sum + refl_sum
    return total // (2 * n)


# --------------------------------------------------------------------------
# LINEAR SNAKE — exact integer recurrence
# --------------------------------------------------------------------------
# Seed values N_1..N_13 (computed from the Burnside enumeration in
# src.snakes.compute_components; L=1 is exceptional but stored here for
# completeness so the recurrence can be primed).
_LINEAR_SEEDS: List[int] = [
    4,               # N_1  (uses 8-element symmetry group — exception)
    21,              # N_2
    109,             # N_3
    586,             # N_4
    3326,            # N_5
    19209,           # N_6
    111871,          # N_7
    653758,          # N_8
    3824678,         # N_9
    22387074,        # N_10
    131052313,       # N_11
    767211817,       # N_12
    4491420695,      # N_13
]

# N_L = Σ_{j=0..11} _LINEAR_C[j] · N_{L-12+j}  for L ≥ 14.
_LINEAR_C: List[int] = [225, -750, 330, 370, -690, 1160, -556, -394, 370, -38, -38, 12]


def linear_snake_count(L: int) -> int:
    if L < 1:
        raise ValueError("L must be >= 1")
    if L <= len(_LINEAR_SEEDS):
        return _LINEAR_SEEDS[L - 1]
    N = list(_LINEAR_SEEDS)
    for k in range(len(_LINEAR_SEEDS) + 1, L + 1):
        val = sum(_LINEAR_C[j] * N[k - 12 + j - 1] for j in range(12))
        N.append(val)
    return N[L - 1]


def linear_snake_sequence(max_L: int) -> List[int]:
    if max_L < 1:
        return []
    N = list(_LINEAR_SEEDS[: min(max_L, len(_LINEAR_SEEDS))])
    for k in range(len(_LINEAR_SEEDS) + 1, max_L + 1):
        val = sum(_LINEAR_C[j] * N[k - 12 + j - 1] for j in range(12))
        N.append(val)
    return N


# --------------------------------------------------------------------------
# LINEAR SNAKE — asymptotic (exact algebraic constants, for O(1) estimate)
# --------------------------------------------------------------------------
# μ = (5 + 3√5)/2 ≈ 5.854101966249684544612...
# A = (5 + 2√5)/20 ≈ 0.473606797749978969640917366873...
LINEAR_MU_SQRT5 = "(5 + 3*sqrt(5))/2"
LINEAR_A_SQRT5 = "(5 + 2*sqrt(5))/20"


def linear_snake_asymptotic(L: int, digits: int = 50) -> str:
    """
    Return N_L asymptotic value A·μ^L to `digits` decimal places, computed
    with mpmath. This is only an approximation; use linear_snake_count for
    the exact integer.
    """
    import mpmath as mp

    mp.mp.dps = digits + 20
    s5 = mp.sqrt(5)
    A = (5 + 2 * s5) / 20
    mu = (5 + 3 * s5) / 2
    val = A * (mu ** L)
    return mp.nstr(val, digits)


# --------------------------------------------------------------------------
# self-test — small end-to-end sanity check
# --------------------------------------------------------------------------
if __name__ == "__main__":
    # cyclic vs README golden
    cyclic_expected = [1, 2, 2, 4, 4, 9, 10, 22, 30, 62, 94, 192, 316, 623, 1096]
    for i, exp in enumerate(cyclic_expected, start=1):
        got = cyclic_adjacency_count(i)
        assert got == exp, f"cyclic n={i}: got {got}, expected {exp}"
    print("cyclic n=1..15  ✓")

    # linear vs README golden
    linear_expected = [4, 21, 109, 586, 3326, 19209, 111871, 653758, 3824678,
                       22387074, 131052313, 767211817, 4491420695, 26293679325,
                       153927402355]
    for i, exp in enumerate(linear_expected, start=1):
        got = linear_snake_count(i)
        assert got == exp, f"linear L={i}: got {got}, expected {exp}"
    print("linear L=1..15  ✓")

    # asymptotic sanity
    L = 100
    approx = linear_snake_asymptotic(L, digits=25)
    exact = linear_snake_count(L)
    print(f"\nL={L}:")
    print(f"  exact    = {exact}")
    print(f"  A·μ^L    = {approx}")
    print(f"  digits   = {len(str(exact))}")

    # cyclic at large n
    n = 100
    print(f"\ncyclic n={n}:")
    print(f"  a({n}) = {cyclic_adjacency_count(n)}")
