#!/usr/bin/env python3
"""
A053656: Cyclic graphs with oriented edges on n nodes (up to dihedral symmetry)

This uses the BINARY edge-orientation encoding:
- Each of the n edges has 2 possible orientations (0 or 1)
- 2^n total configurations
- Count equivalence classes under dihedral symmetry (rotation + reflection)
"""

from typing import Iterator, Set, Tuple, Dict
from itertools import product
from math import gcd


Sequence = Tuple[int, ...]


def all_rotations(seq: Sequence) -> Iterator[Sequence]:
    """Generate all rotations of a sequence."""
    n = len(seq)
    for i in range(n):
        yield seq[i:] + seq[:i]


def canonical_form(seq: Sequence) -> Sequence:
    """
    Return the lexicographically smallest form under dihedral symmetry.
    For edge orientations, reflection reverses the sequence AND flips each bit.
    """
    if len(seq) == 0:
        return seq
    
    candidates = []
    
    # All rotations of original
    for rot in all_rotations(seq):
        candidates.append(rot)
    
    # Reflection: reverse AND flip each orientation
    # When you flip the cycle, edge i becomes edge (n-1-i) and its direction flips
    reflected = tuple(1 - x for x in reversed(seq))
    for rot in all_rotations(reflected):
        candidates.append(rot)
    
    return min(candidates)


def generate_oriented_cycles(n: int) -> Set[Sequence]:
    """
    Generate all oriented cycles on n edges under dihedral symmetry.
    """
    valid = set()
    
    for seq in product([0, 1], repeat=n):
        canon = canonical_form(seq)
        valid.add(canon)
    
    return valid


def a053656_formula(n: int) -> int:
    """
    Compute A053656(n) using the OEIS formula:
    a(n) = A000031(n)/2 + (if n even) 2^(n/2-2)
    
    where A000031(n) = (1/n) * sum_{d|n} phi(d) * 2^(n/d)
    """
    from sympy import divisors, totient
    
    # A000031(n) = number of binary necklaces of length n
    a000031 = sum(totient(d) * (2 ** (n // d)) for d in divisors(n)) // n
    
    result = a000031 // 2
    if n % 2 == 0:
        result += 2 ** (n // 2 - 2)
    
    return result


def convert_to_vertex_signature(edge_seq: Sequence) -> Sequence:
    """
    Convert an edge orientation sequence to vertex signature sequence.
    
    edge_seq[i] = orientation of edge between vertex i and vertex (i+1) mod n
    - 0 means edge points from i to i+1
    - 1 means edge points from i+1 to i
    
    Vertex signature = (edges pointing in) - (edges pointing out)
    """
    n = len(edge_seq)
    if n == 0:
        return ()
    
    signatures = []
    for v in range(n):
        # Edge from v-1 to v (edge index v-1)
        edge_left = edge_seq[(v - 1) % n]
        # If edge_left == 0: points from v-1 to v, so IN to v
        # If edge_left == 1: points from v to v-1, so OUT from v
        in_from_left = 1 if edge_left == 0 else 0
        out_to_left = 1 - in_from_left
        
        # Edge from v to v+1 (edge index v)
        edge_right = edge_seq[v]
        # If edge_right == 0: points from v to v+1, so OUT from v
        # If edge_right == 1: points from v+1 to v, so IN to v
        out_to_right = 1 if edge_right == 0 else 0
        in_from_right = 1 - out_to_right
        
        total_in = in_from_left + in_from_right
        total_out = out_to_left + out_to_right
        signature = (total_in - total_out) * 1  # Will be -2, 0, or +2
        
        # Wait, total_in + total_out = 2 always
        # signature = total_in - total_out = total_in - (2 - total_in) = 2*total_in - 2
        # So signature ∈ {-2, 0, 2}
        signatures.append(signature)
    
    return tuple(signatures)


def main():
    import argparse
    
    parser = argparse.ArgumentParser(description="Enumerate A053656")
    parser.add_argument("--max-n", type=int, default=15)
    parser.add_argument("--verbose", "-v", action="store_true")
    parser.add_argument("--show-formula", action="store_true", help="Also show formula-based values")
    args = parser.parse_args()
    
    print(f"Enumerating A053656 for n=1..{args.max_n}")
    print("=" * 60)
    
    counts = []
    for n in range(1, args.max_n + 1):
        cycles = generate_oriented_cycles(n)
        counts.append(len(cycles))
        
        if args.show_formula:
            formula_val = a053656_formula(n)
            match = "✓" if len(cycles) == formula_val else "✗"
            print(f"n={n}: enumerated={len(cycles)}, formula={formula_val} {match}")
        else:
            print(f"n={n}: {len(cycles)}")
        
        if args.verbose and len(cycles) <= 20:
            for seq in sorted(cycles):
                sig = convert_to_vertex_signature(seq)
                sig_str = "".join(f"{s:+d}" if s != 0 else "0" for s in sig)
                print(f"  edges={seq} -> vertex_sig=({sig_str})")
    
    print("=" * 60)
    print(f"a(1..{args.max_n}) = {counts}")
    
    # Expected A053656 values
    expected = [1, 2, 2, 4, 4, 9, 10, 22, 30, 62, 94, 192, 316, 623, 1096]
    if args.max_n <= len(expected):
        match = counts == expected[:args.max_n]
        print(f"Expected:        {expected[:args.max_n]}")
        print(f"Match: {'✓ YES' if match else '✗ NO'}")


if __name__ == "__main__":
    main()
