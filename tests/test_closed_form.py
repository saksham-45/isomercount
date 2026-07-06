#!/usr/bin/env python3
"""
Tests for the exact closed-form counters in src.analysis.closed_form.

Every claim in CLOSED_FORM.md is verified here against the ground-truth
Burnside enumerator implementations in src.cyclic and src.snakes.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.analysis.closed_form import (
    cyclic_adjacency_count,
    linear_snake_count,
    linear_snake_sequence,
    linear_snake_asymptotic,
)
from src.cyclic.cyclic_sequences import count_adjacency_burnside
from src.snakes.compute_components import ComponentSolver


def test_cyclic_matches_burnside(n_max=50):
    fails = 0
    for n in range(1, n_max + 1):
        got = cyclic_adjacency_count(n)
        expected = count_adjacency_burnside(n)
        if got != expected:
            print(f"  cyclic n={n}: {got} vs Burnside {expected}")
            fails += 1
    assert fails == 0, f"{fails} mismatches"
    print(f"test_cyclic_matches_burnside n=1..{n_max}: OK")


def test_linear_snake_matches_burnside(L_max=40):
    solver = ComponentSolver()

    def truth(L):
        if L == 1:
            G0 = [(di, False) for di in range(8)]
            sv = sum(solver.count_fixed(L, di, rev, False) for di, rev in G0)
            sb = sum(solver.count_fixed(L, di, rev, False, True) for di, rev in G0)
            sbf = sum(solver.count_fixed(L, di, rev, True, True) for di, rev in G0)
            return int((2 * sv - sb + sbf) // (2 * 8))
        G0 = [(0, False), (4, False), (2, True), (5, True)]
        sv = sum(solver.count_fixed(L, di, rev, False) for di, rev in G0)
        sb = sum(solver.count_fixed(L, di, rev, False, True) for di, rev in G0)
        sbf = sum(solver.count_fixed(L, di, rev, True, True) for di, rev in G0)
        return int((2 * sv - sb + sbf) // (2 * 4))

    fails = 0
    for L in range(1, L_max + 1):
        got = linear_snake_count(L)
        expected = truth(L)
        if got != expected:
            print(f"  linear L={L}: {got} vs Burnside {expected}")
            fails += 1
    assert fails == 0, f"{fails} mismatches"
    print(f"test_linear_snake_matches_burnside L=1..{L_max}: OK")


def test_linear_snake_sequence_helper(L_max=30):
    seq = linear_snake_sequence(L_max)
    assert len(seq) == L_max
    for L in range(1, L_max + 1):
        assert seq[L - 1] == linear_snake_count(L), f"seq[{L-1}] != count({L})"
    print(f"test_linear_snake_sequence_helper L_max={L_max}: OK")


def test_asymptotic_convergence():
    """Check N_L / (A * μ^L) → 1 fast."""
    import mpmath as mp

    mp.mp.dps = 60
    for L in [30, 50, 100]:
        asy = mp.mpf(linear_snake_asymptotic(L, digits=40))
        actual = linear_snake_count(L)
        ratio = actual / asy
        assert abs(ratio - 1) < mp.mpf(10) ** (-8), f"L={L} ratio={ratio}"
    print("test_asymptotic_convergence L=30,50,100: OK")


def test_readme_golden_values():
    """The values shipped in README.md must match."""
    cyclic_expected = [1, 2, 2, 4, 4, 9, 10, 22, 30, 62, 94, 192, 316, 623, 1096]
    for i, exp in enumerate(cyclic_expected, start=1):
        assert cyclic_adjacency_count(i) == exp

    linear_expected = [4, 21, 109, 586, 3326, 19209, 111871, 653758, 3824678,
                       22387074, 131052313, 767211817, 4491420695,
                       26293679325, 153927402355]
    for i, exp in enumerate(linear_expected, start=1):
        assert linear_snake_count(i) == exp
    print("test_readme_golden_values: OK")


def test_scales_to_large_L():
    """The recurrence should handle L=500 without issue."""
    val = linear_snake_count(500)
    assert isinstance(val, int) and val > 0
    # Should have on the order of 500 * log10(μ) ≈ 500 * 0.767 ≈ 383 digits
    ndigits = len(str(val))
    assert 380 < ndigits < 386, f"L=500 gave {ndigits} digits"
    print(f"test_scales_to_large_L: N_500 has {ndigits} digits, OK")


if __name__ == "__main__":
    test_readme_golden_values()
    test_cyclic_matches_burnside()
    test_linear_snake_matches_burnside()
    test_linear_snake_sequence_helper()
    test_asymptotic_convergence()
    test_scales_to_large_L()
    print("\nAll closed-form tests passed.")
