#!/usr/bin/env python3
"""Regenerate the linear snake orbit sequence, print values, and verify recurrences.

Uses ComponentSolver (same logic as the linear-count fast path) to avoid copy/paste drift.

Notes:
- By default this script computes/prints L=1..n and checks the known recurrences.
- Building a full symbolic closed form is possible but slow; it is behind --closed-form.
"""

from __future__ import annotations

import argparse
import sympy as sp

from .compute_components import ComponentSolver

solver = ComponentSolver()


def linear_orbit(L: int) -> int:
    """Orbit count for linear snake of length L."""
    if L == 1:
        G0 = [(di, False) for di in range(8)]
        sv0 = sum(solver.count_fixed(L, di, rev, False) for di, rev in G0)
        sb0 = sum(solver.count_fixed(L, di, rev, False, True) for di, rev in G0)
        sbf = sum(solver.count_fixed(L, di, rev, True, True) for di, rev in G0)
        return int((2 * sv0 - sb0 + sbf) // (2 * 8))

    G0 = [(0, False), (4, False), (2, True), (5, True)]
    sv0 = sum(solver.count_fixed(L, di, rev, False) for di, rev in G0)
    sb0 = sum(solver.count_fixed(L, di, rev, False, True) for di, rev in G0)
    sbf = sum(solver.count_fixed(L, di, rev, True, True) for di, rev in G0)
    return int((2 * sv0 - sb0 + sbf) // (2 * 4))


def sequence(n: int):
    return [linear_orbit(L) for L in range(1, n + 1)]


def verify_recurrence(seq, coeff, start_L):
    d = len(coeff)
    for L in range(start_L, len(seq) + 1):
        rhs = sum(coeff[i] * seq[L - 2 - i] for i in range(d))
        if rhs != seq[L - 1]:
            return False, L
    return True, None


def verify_subsequence(arr, coeff, start_idx):
    """Verify recurrence on a plain 0-based subsequence array."""
    d = len(coeff)
    for n in range(start_idx, len(arr)):
        rhs = sum(coeff[i] * arr[n - 1 - i] for i in range(d))
        if rhs != arr[n]:
            return False, n
    return True, None


def exact_roots():
    """Return the 12 characteristic roots as exact algebraic SymPy expressions."""
    sqrt5 = sp.sqrt(5)
    sqrt3 = sp.sqrt(3)
    r1 = sp.Integer(1)
    r2 = sp.Integer(3)
    r3 = sqrt3
    r4 = -sqrt3
    r5 = (5 + 3 * sqrt5) / 2
    r6 = (5 - 3 * sqrt5) / 2
    r7 = (3 + sqrt5) / 2
    r8 = (3 - sqrt5) / 2
    # quartic roots: x = +/-sqrt((5 +/- 3*sqrt5)/2)
    r9 = sp.sqrt((5 + 3 * sqrt5) / 2)
    r10 = -r9
    r11 = sp.sqrt((5 - 3 * sqrt5) / 2)
    r12 = -r11
    return [r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12]


def exact_closed_form():
    """Build algebraic closed-form coefficients A_i for a(L) = sum A_i*r_i^(L-1)."""
    roots = exact_roots()
    init = sequence(12)  # a(1)..a(12)
    V = sp.Matrix([[r ** k for r in roots] for k in range(12)])  # k=0..11
    b = sp.Matrix(init)
    coeffs = list(sp.simplify(x) for x in V.LUsolve(b))
    return roots, coeffs


def eval_closed_form(L: int, roots, coeffs, prec: int = 100) -> int:
    """Evaluate the closed form at length L, return nearest int."""
    val = sum(coeffs[j] * (roots[j] ** (L - 1)) for j in range(len(roots)))
    return int(sp.N(sp.re(val), prec).round())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=120, help="Compute and print L=1..n (default: 120).")
    ap.add_argument(
        "--closed-form",
        action="store_true",
        help="Also build a (slow) algebraic closed form from characteristic roots and spot-check it.",
    )
    args = ap.parse_args()

    n = args.n
    seq = sequence(n)

    for L, val in enumerate(seq, start=1):
        print(f"L={L:3d}: {val}")

    # Full recurrence (order 12, valid L>=14)
    full_coeff = [12, -38, -38, 370, -394, -556, 1160, -690, 370, 330, -750, 225]
    ok, badL = verify_recurrence(seq, full_coeff, start_L=14)
    print("\nFull recurrence check:", "OK" if ok else f"FAIL at L={badL}")

    # Even subsequence (order 9, valid from L=20 onward)
    even = [seq[2 * k + 1] for k in range((len(seq) - 1) // 2)]  # L=2,4,...
    even_coeff = [60, -1126, 9403, -37168, 61132, -13515, -45260, 29850, -3375]
    ok, bad = verify_subsequence(even, even_coeff, start_idx=9)  # k index 9 => L=20
    print("Even-subsequence recurrence check:", "OK" if ok else f"FAIL at even k index {bad}")

    # Odd subsequence (order 8, valid from L=19 onward)
    odd = [seq[2 * k + 2] for k in range((len(seq) - 2) // 2)]  # L=3,5,...
    odd_coeff = [57, -955, 6538, -17554, 8470, 11895, -9575, 1125]
    ok, bad = verify_subsequence(odd, odd_coeff, start_idx=8)  # k index 8 => L=19
    print("Odd-subsequence recurrence check:", "OK" if ok else f"FAIL at odd k index {bad}")

    if args.closed_form:
        roots, coeffs = exact_closed_form()
        sample_L = [20, 40, 60, 80, 100]
        print("\nExact closed-form spot checks (rounded integer):")
        for L in sample_L:
            cf_val = eval_closed_form(L, roots, coeffs, prec=150)
            match = "match" if cf_val == seq[L - 1] else "DIFF"
            print(f"L={L:3d}: {cf_val} ({match})")


if __name__ == "__main__":
    main()
