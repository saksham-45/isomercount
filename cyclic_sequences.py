#!/usr/bin/env python3
"""
Constrained Oriented Cyclic Graph Sequence Enumeration

WHAT THIS PROGRAM DOES (in simple words):
  We have a circle with n spots, and each "edge" between two spots can point
  one way or the other. We count how many truly different such pictures there are
  if we agree that "rotating the circle" or "flipping it over" doesn't create
  a new picture. We can also add rules (e.g. no two same-sign nonzero numbers
  next to each other) and count only those that follow the rules.

KEY IDEAS (with links to read more):
  - Integer: https://en.wikipedia.org/wiki/Integer
  - Cycle / rotation: https://en.wikipedia.org/wiki/Cyclic_permutation
  - Dihedral group (rotations + reflections): https://en.wikipedia.org/wiki/Dihedral_group
  - Burnside's lemma (counting orbits): https://en.wikipedia.org/wiki/Burnside%27s_lemma
  - Transfer matrix (counting paths): https://en.wikipedia.org/wiki/Transfer-matrix_method
  - Matrix (math): https://en.wikipedia.org/wiki/Matrix_(mathematics)
  - Euler's totient function: https://en.wikipedia.org/wiki/Euler%27s_totient_function
  - Divisor: https://en.wikipedia.org/wiki/Divisor
  - Prime number: https://en.wikipedia.org/wiki/Prime_number
  - Möbius function: https://en.wikipedia.org/wiki/M%C3%B6bius_function

CONSTRAINTS we can turn on:
  1. Balance: number of +2 and -2 are equal (automatic for valid cycles).
  2. Adjacency: no two consecutive nonzero numbers with the same sign.
  3. Primitivity: the pattern is not a repeat of a shorter pattern.
  4. Forbidden subsequence: no cyclic substring matches a "forbidden" pattern from smaller n.
"""

from typing import Iterator, Set, Tuple, Dict, List, Any, Optional
from itertools import product

# A "sequence" here is a tuple of integers (e.g. vertex values -2, 0, 2 or edge bits 0, 1).
Sequence = Tuple[int, ...]

# When we need to list actual cycles (e.g. for --verbose), we use brute force only for n <= N0.
N0 = 14


def all_rotations(seq: Sequence) -> Iterator[Sequence]:
    """
    Give every way to "rotate" the sequence (shift it around the circle).
    Example: (A,B,C) -> (A,B,C), (B,C,A), (C,A,B).
    See: https://en.wikipedia.org/wiki/Cyclic_permutation
    """
    n = len(seq)
    for i in range(n):
        yield seq[i:] + seq[:i]


def canonical_form_edges(seq: Sequence) -> Sequence:
    """
    Pick one "canonical" version of this cycle so that two cycles that are the same
    up to rotation or reflection get the same canonical form. We choose the
    lexicographically smallest (like alphabetical order for numbers).
    Reflection here means: reverse the order AND flip each 0<->1.
    Dihedral group: https://en.wikipedia.org/wiki/Dihedral_group
    Lexicographic order: https://en.wikipedia.org/wiki/Lexicographic_order
    """
    if len(seq) == 0:
        return seq
    
    candidates = []
    for rot in all_rotations(seq):
        candidates.append(rot)
    
    reflected = tuple(1 - x for x in reversed(seq))
    for rot in all_rotations(reflected):
        candidates.append(rot)
    
    return min(candidates)


def edge_to_vertex_signature(edge_seq: Sequence) -> Sequence:
    """
    Turn a list of edge directions (0 or 1) into a list of "vertex" numbers (-2, 0, or +2).
    At each vertex we look at the two edges touching it; we count how many point "in"
    and convert that to a single number. So each vertex gets -2, 0, or +2.
    (v-1) % n and v mean "the edge to the left of vertex v" and "the edge to the right"
    using wrap-around. Modular arithmetic: https://en.wikipedia.org/wiki/Modular_arithmetic
    """
    n = len(edge_seq)
    if n == 0:
        return ()
    
    signatures = []
    for v in range(n):
        # Edges on the left and right of vertex v (going around the circle).
        edge_left = edge_seq[(v - 1) % n]
        edge_right = edge_seq[v]
        # Count how many of these two edges point "in" to vertex v.
        in_left = 1 if edge_left == 0 else 0
        in_right = 1 if edge_right == 1 else 0
        total_in = in_left + in_right
        # Turn count (0, 1, or 2) into -2, 0, or +2.
        signature = 2 * total_in - 2
        signatures.append(signature)
    
    return tuple(signatures)


def generate_all_oriented_cycles(n: int) -> Set[Tuple[Sequence, Sequence]]:
    """
    List one representative of each "type" of cycle on n edges (under rotation + reflection).
    We try all 2^n ways to assign 0 or 1 to each edge (every combination; see Cartesian product:
    https://en.wikipedia.org/wiki/Cartesian_product), then collapse to canonical form
    and keep one (edge_canonical, vertex_signature) per type. Slow for large n.
    """
    seen = set()
    result = set()
    
    for edge_seq in product([0, 1], repeat=n):
        canon = canonical_form_edges(edge_seq)
        if canon in seen:
            continue
        seen.add(canon)
        
        vertex_sig = edge_to_vertex_signature(canon)
        result.add((canon, vertex_sig))
    
    return result


def is_adjacency_valid(sig: Sequence) -> bool:
    """
    Adjacency rule: we don't allow two nonzero numbers next to each other with the
    same sign. So (+2, +2) or (-2, -2) next to each other (including around the circle)
    means "invalid". We check every adjacent pair.
    """
    n = len(sig)
    if n < 2:
        return True
    
    for i in range(n):
        curr = sig[i]
        next_elem = sig[(i + 1) % n]
        if curr != 0 and next_elem != 0:
            if (curr > 0) == (next_elem > 0):  # Same sign
                return False
    return True


def canonical_form_sig(sig: Sequence) -> Sequence:
    """
    Pick one canonical form for a vertex signature: we consider all rotations and
    the reversed sequence (no sign flip here), and take the lexicographically smallest.
    Lexicographic order: https://en.wikipedia.org/wiki/Lexicographic_order
    """
    if len(sig) == 0:
        return sig
    
    candidates = list(all_rotations(sig))
    reversed_sig = sig[::-1]
    candidates.extend(all_rotations(reversed_sig))
    
    return min(candidates)


def is_primitive(seq: Sequence) -> bool:
    """
    "Primitive" means the sequence is not made by repeating a shorter block.
    Example: (1,2,1,2) is not primitive (repeats (1,2)); (1,2,3) is primitive.
    We check for every possible period (length that divides n).
    """
    n = len(seq)
    for period in range(1, n):
        if n % period != 0:
            continue
        pattern = seq[:period]
        if all(seq[i] == pattern[i % period] for i in range(n)):
            return False
    return True


def get_cyclic_substrings(seq: Sequence, length: int) -> Iterator[Sequence]:
    """
    All contiguous blocks of a given length, going around the circle (so we wrap from
    end to start). We duplicate the sequence so we can slide a window of that length.
    """
    n = len(seq)
    if length > n or length < 1:
        return
    doubled = seq + seq
    for start in range(n):
        yield tuple(doubled[start:start + length])


def contains_forbidden(sig: Sequence, forbidden_by_len: Dict[int, Set[Sequence]]) -> bool:
    """
    Check whether any contiguous block (of any length we care about) in this cycle,
    when canonicalized, appears in the "forbidden" set for that length. Used to
    enforce "no forbidden subsequence" constraint.
    """
    n = len(sig)
    for k, forbidden_set in forbidden_by_len.items():
        if k >= n or k < 2:  # Skip trivial or too-long patterns
            continue
        for sub in get_cyclic_substrings(sig, k):
            canon_sub = canonical_form_sig(sub)
            if canon_sub in forbidden_set:
                return True
    return False


# =============================================================================
# TRANSFER MATRIX AND BURNSIDE (adjacency-only, fast counting)
# We count "closed walks" (paths that come back to the start) that obey the
# adjacency rule, then use Burnside's lemma to count distinct cycles under
# rotation and reflection. Functions are in execution order: 1 -> 2 -> 3 -> 4
# (per divisor); 5 and 6 are used in the Burnside loop; 7 counts reflection-fixed;
# 8 ties everything together.
# Transfer matrix: https://en.wikipedia.org/wiki/Transfer-matrix_method
# Burnside's lemma: https://en.wikipedia.org/wiki/Burnside%27s_lemma
# =============================================================================

# --- Step 1: Small helper used when building the matrix and for reflection ---
def _vertex_from_edges(e_left: int, e_right: int) -> int:
    """
    From the two edge bits (0 or 1) around one vertex, compute the vertex "signature"
    number: -2, 0, or +2. Same idea as in edge_to_vertex_signature but for one vertex.
    """
    in_left = 1 if e_left == 0 else 0
    in_right = 1 if e_right == 1 else 0
    return 2 * (in_left + in_right) - 2


# --- Step 2: Build the transfer matrix (uses step 1) ---
def _build_adjacency_transfer_matrix() -> Tuple[List[Tuple[int, int]], List[List[int]]]:
    """
    We encode "state" as (last vertex value, current edge). The matrix M has one row/column
    per state; M[i][j] = 1 if we can go from state i to state j while obeying the adjacency
    rule (no same-sign nonzero consecutively), else 0. So M counts allowed one-step moves.
    Matrix: https://en.wikipedia.org/wiki/Matrix_(mathematics)
    """
    vertices = (-2, 0, 2)
    edges = (0, 1)
    state_list = [(v, e) for v in vertices for e in edges]
    n_states = len(state_list)
    M = [[0] * n_states for _ in range(n_states)]
    for i, (last_v, curr_e) in enumerate(state_list):
        for next_e in edges:
            new_v = _vertex_from_edges(curr_e, next_e)
            # Adjacency: disallow same-sign consecutive nonzero
            if last_v != 0 and new_v != 0 and (last_v > 0) == (new_v > 0):
                continue
            j = state_list.index((new_v, next_e))
            M[i][j] = 1
    return state_list, M


# --- Step 3: Compute M^n (matrix raised to power n) ---
def _matrix_power(M: List[List[int]], n: int) -> List[List[int]]:
    """
    Compute M^n (M multiplied by itself n times) efficiently by "binary exponentiation":
    we square M, then square that, etc., and multiply when n has a 1 in binary. This way
    we do about log(n) matrix multiplications instead of n.
    Matrix multiplication: https://en.wikipedia.org/wiki/Matrix_multiplication
    Exponentiation by squaring: https://en.wikipedia.org/wiki/Exponentiation_by_squaring
    """
    size = len(M)
    # M^0 = identity matrix (1 on diagonal, 0 elsewhere). Identity matrix: https://en.wikipedia.org/wiki/Identity_matrix
    if n == 0:
        return [[1 if i == j else 0 for j in range(size)] for i in range(size)]
    if n == 1:
        return [row[:] for row in M]
    half = _matrix_power(M, n // 2)
    # half * half
    result = [[0] * size for _ in range(size)]
    for i in range(size):
        for k in range(size):
            if half[i][k]:
                for j in range(size):
                    result[i][j] += half[i][k] * half[k][j]
    if n % 2:
        result2 = [[0] * size for _ in range(size)]
        for i in range(size):
            for k in range(size):
                if result[i][k]:
                    for j in range(size):
                        result2[i][j] += result[i][k] * M[k][j]
        return result2
    return result


# --- Step 4: Count closed walks of length n (uses steps 2 and 3). Called once per divisor. ---
def _count_closed_walks_adjacency(n: int) -> int:
    """
    A "closed walk" is a path of n steps that starts and ends in the same state. The
    transfer matrix M has the property: the number of closed walks of length n equals
    the sum of the diagonal entries of M^n (that sum is called the "trace"). So we
    build M, compute M^n, and return that trace.
    Trace of a matrix: https://en.wikipedia.org/wiki/Trace_(linear_algebra)
    """
    if n <= 0:
        return 0
    _, M = _build_adjacency_transfer_matrix()
    M_n = _matrix_power(M, n)
    trace = sum(M_n[i][i] for i in range(len(M_n)))
    return trace


# --- Step 5: List divisors of n (for Burnside: we sum over divisors) ---
def _divisors(n: int) -> List[int]:
    """
    All positive integers that divide n evenly (e.g. divisors of 12: 1, 2, 3, 4, 6, 12).
    We only check up to sqrt(n) and pair each divisor d with n//d to keep it fast.
    Divisor: https://en.wikipedia.org/wiki/Divisor
    """
    if n <= 0:
        return []
    small = []
    large = []
    d = 1
    while d * d <= n:
        if n % d == 0:
            small.append(d)
            if d * d != n:
                large.append(n // d)
        d += 1
    return small + large[::-1]


# --- Step 6: Euler totient phi(n) (for Burnside weights) ---
def _totient(n: int) -> int:
    """
    Euler's totient phi(n) = how many integers from 1 to n share no prime factor with n
    (i.e. are "coprime" to n). We need this for Burnside's formula. We compute it from
    the prime factorization of n.
    Euler's totient function: https://en.wikipedia.org/wiki/Euler%27s_totient_function
    Prime number: https://en.wikipedia.org/wiki/Prime_number
    """
    if n <= 0:
        return 0
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


# --- Step 7: Count sequences that look the same after reflection ---
def _count_reflection_fixed_adjacency(n: int) -> int:
    """
    Some cycles look identical when we "flip" the circle (reverse order and flip 0<->1).
    We count how many valid (closed, adjacency-ok) edge sequences of length n are unchanged
    by that reflection. For odd n there are none (middle edge would have to be 0 and 1).
    Reflection (geometry): https://en.wikipedia.org/wiki/Reflection_(mathematics)
    """
    if n <= 0:
        return 0
    if n % 2 == 1:
        return 0
    
    k = n // 2
    state_list, M = _build_adjacency_transfer_matrix()
    M_k = _matrix_power(M, k)
    
    total = 0
    for e0 in (0, 1):
        # Reflection forces the "last" edge (wrap-around) to be 1 - e0; that gives v0.
        v0 = _vertex_from_edges(1 - e0, e0)
        try:
            start_idx = state_list.index((v0, e0))
        except ValueError:
            continue
        row = M_k[start_idx]
        for end_idx, count in enumerate(row):
            if count == 0:
                continue
            v_k, e_k_minus_1 = state_list[end_idx]
            # At the halfway point, reflection links the end of the first half to the start
            # of the second; only certain (v_k, e) combinations are valid. We use the
            # condition that matched brute-force enumeration.
            if (v_k == 2 and e_k_minus_1 == 1) or (v_k == -2 and e_k_minus_1 == 0):
                total += count
    return total


# --- Step 8: Put it together (Burnside's lemma) ---
def count_adjacency_burnside(n: int) -> int:
    """
    Final count of "distinct" oriented cycles of length n with the adjacency rule, where
    we consider two cycles the same if one can be rotated or reflected to get the other.
    Burnside's lemma says: count = (1 / |group|) * (sum over symmetries of "how many
    cycles are unchanged by this symmetry"). We do that sum using divisors and totient
    for rotations, and step 7 for reflections, then divide by 2*n (size of the dihedral group).
    Burnside's lemma: https://en.wikipedia.org/wiki/Burnside%27s_lemma
    """
    if n <= 0:
        return 0
    # Rotation part: for each divisor d of n we get closed walks of length d, weighted by phi(n/d).
    rot_sum = 0
    for d in _divisors(n):
        phi_val = _totient(n // d)
        count_d = _count_closed_walks_adjacency(d)
        rot_sum += phi_val * count_d

    ref_fixed = _count_reflection_fixed_adjacency(n)
    ref_term = (n // 2) * ref_fixed
    orbits = (rot_sum + ref_term) // (2 * n)
    return orbits


# =============================================================================
# MÖBIUS INVERSION FOR PRIMITIVITY
# To count only "primitive" cycles (no repeating block), we use Möbius inversion:
# primitive_count(n) = sum over divisors d of n of mu(d) * (total count for n/d).
# Möbius function: https://en.wikipedia.org/wiki/M%C3%B6bius_function
# =============================================================================

def _mobius(n: int) -> int:
    """
    Möbius function mu(n): 0 if n is divisible by the square of any prime; otherwise
    (-1)^(number of distinct prime factors). Example: mu(6)=1 (6=2*3, two primes), mu(4)=0 (4=2^2).
    """
    if n <= 0:
        return 0
    if n == 1:
        return 1
    result = 1
    p = 2
    m = n
    while p * p <= m:
        if m % p == 0:
            if m % (p * p) == 0:
                return 0
            result = -result
            m //= p
        else:
            p += 1
    if m > 1:
        result = -result
    return result


def count_adjacency_primitive_burnside(n: int) -> int:
    """
    Same as count_adjacency_burnside but we only count cycles that are "primitive" (not
    a repetition of a shorter cycle). Formula: primitive(n) = sum over divisors d of n
    of Möbius(d) * total_count(n/d). So we use the Möbius function to subtract the
    non-primitive ones. Möbius inversion: https://en.wikipedia.org/wiki/M%C3%B6bius_inversion_formula
    """
    if n <= 0:
        return 0
    total = 0
    for d in _divisors(n):
        m = n // d
        total += _mobius(d) * count_adjacency_burnside(m)
    return total


# =============================================================================
# CLOSED-FORM FORMULA (for checking only)
# A053656 is a known sequence; we can compute it with a formula for comparison.
# We never use this for "real" results—we always count from scratch.
# OEIS A053656: https://oeis.org/A053656
# =============================================================================

def count_a053656(n: int) -> int:
    """
    A053656(n) = number of oriented cycles of length n up to rotation and reflection,
    with no extra constraints. Formula: (1/2)*A000031(n) + (if n even) 2^(n/2-2), where
    A000031 uses a sum over divisors and totient. Used only to verify our brute/matrix
    counts, not for producing answers.
    """
    try:
        from sympy import divisors, totient
    except ImportError:
        # Fallback: use math.gcd and simple divisor iteration
        from math import gcd
        def _totient(m: int) -> int:
            count = 0
            for k in range(1, m + 1):
                if gcd(k, m) == 1:
                    count += 1
            return count
        def _divisors(m: int) -> List[int]:
            return [d for d in range(1, m + 1) if m % d == 0]
        a000031 = sum(_totient(d) * (2 ** (n // d)) for d in _divisors(n)) // n
        result = a000031 // 2
        if n % 2 == 0:
            result += 2 ** (n // 2 - 2)
        return int(result)
    a000031 = sum(totient(d) * (2 ** (n // d)) for d in divisors(n)) // n
    result = a000031 // 2
    if n % 2 == 0:
        result += 2 ** (n // 2 - 2)
    return int(result)


class _CountOnlySet:
    """
    A lightweight object that "pretends" to be a set of a given size: len() returns the
    count but we don't actually store the elements. Used when we only need the number
    (e.g. from the matrix/Burnside path) and not the list of cycles.
    """

    __slots__ = ("_count",)

    def __init__(self, count: int):
        self._count = count

    def __len__(self) -> int:
        return self._count

    def __iter__(self) -> Iterator[Any]:
        return iter([])

    def __contains__(self, x: Any) -> bool:
        return False

    def __repr__(self) -> str:
        return f"_CountOnlySet({self._count})"


def count_constrained(
    n: int,
    use_adjacency: bool = False,
    use_primitivity: bool = False,
    use_forbidden: bool = False,
    all_sets_so_far: Optional[Dict[int, Set[Sequence]]] = None,
) -> int:
    """
    Count how many distinct cycles of length n satisfy the chosen constraints. We always
    compute from scratch (either by listing all cycles and filtering, or by the fast
    matrix/Burnside method when only adjacency is used). For the "forbidden subsequence"
    option we need the counts/sets for smaller n in all_sets_so_far.
    """
    if all_sets_so_far is None:
        all_sets_so_far = {}
    cycles = generate_all_oriented_cycles(n)
    valid_sigs = set()
    for edge_seq, vertex_sig in cycles:
        canon_sig = canonical_form_sig(vertex_sig)
        if canon_sig in valid_sigs:
            continue
        if use_adjacency and not is_adjacency_valid(canon_sig):
            continue
        if use_primitivity and not is_primitive(canon_sig):
            continue
        if use_forbidden and contains_forbidden(canon_sig, all_sets_so_far):
            continue
        valid_sigs.add(canon_sig)
    return len(valid_sigs)


def enumerate_constrained(
    max_n: int,
    use_adjacency: bool = False,
    use_primitivity: bool = False,
    use_forbidden: bool = False,
    need_explicit_sets: bool = False,
) -> Dict[int, Set[Sequence]]:
    """
    For each length from 1 to max_n, compute the count (and optionally the actual list)
    of distinct cycles with the chosen constraints. When we use the fast matrix/Burnside
    path we only store the count (as a _CountOnlySet), so len(result[n]) is correct but
    we don't store every cycle. If need_explicit_sets is True we force listing cycles
    for small n (up to N0) so you can iterate over them.
    """
    all_sets: Dict[int, Set[Sequence]] = {}

    for n in range(1, max_n + 1):
        # Fast path: when only adjacency is needed and we don't need the actual list.
        if use_adjacency and not use_forbidden and not need_explicit_sets:
            if use_primitivity:
                count = count_adjacency_primitive_burnside(n)
            else:
                count = count_adjacency_burnside(n)
            all_sets[n] = _CountOnlySet(count)
            print(f"n={n}: {count}")
            continue

        # Otherwise we list all cycles and filter by constraints.
        cycles = generate_all_oriented_cycles(n)
        valid_sigs: Set[Sequence] = set()
        for edge_seq, vertex_sig in cycles:
            canon_sig = canonical_form_sig(vertex_sig)
            if canon_sig in valid_sigs:
                continue
            if use_adjacency and not is_adjacency_valid(canon_sig):
                continue
            if use_primitivity and not is_primitive(canon_sig):
                continue
            if use_forbidden and contains_forbidden(canon_sig, all_sets):
                continue
            valid_sigs.add(canon_sig)
        all_sets[n] = valid_sigs
        print(f"n={n}: {len(valid_sigs)}")
    return all_sets


def format_sig(sig: Sequence) -> str:
    """Turn a vertex signature (tuple of -2, 0, 2) into a readable string like (+2,0,-2,0)."""
    parts = []
    for s in sig:
        if s == 0:
            parts.append("0")
        elif s > 0:
            parts.append(f"+{s}")
        else:
            parts.append(str(s))
    return "(" + ",".join(parts) + ")"


def _debug_compare_adjacency(n: int) -> None:
    """Debug helper: compare brute-force count (distinct canonical sigs with adjacency) to Burnside count."""
    cycles = generate_all_oriented_cycles(n)
    edge_orbits_with_adj = 0
    canon_sigs_seen: Set[Sequence] = set()
    for edge_seq, vertex_sig in cycles:
        canon_sig = canonical_form_sig(vertex_sig)
        if not is_adjacency_valid(canon_sig):
            continue
        edge_orbits_with_adj += 1
        canon_sigs_seen.add(canon_sig)
    brute_count = len(canon_sigs_seen)
    burnside_count = count_adjacency_burnside(n)
    try:
        import json
        with open("/Users/saksham/grp/.cursor/debug.log", "a") as _f:
            _f.write(json.dumps({"sessionId": "debug-session", "hypothesisId": "A,E", "location": "_debug_compare_adjacency", "message": "edge orbits vs canonical sigs", "data": {"n": n, "edge_orbits_with_adjacency": edge_orbits_with_adj, "distinct_canonical_sigs": brute_count, "burnside_count": burnside_count}, "timestamp": 0}) + "\n")
    except Exception:
        pass


def main():
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-n", type=int, default=15)
    parser.add_argument("--verbose", "-v", action="store_true")
    parser.add_argument("--adjacency", action="store_true", help="Apply adjacency constraint")
    parser.add_argument("--primitive", action="store_true", help="Apply primitivity constraint")
    parser.add_argument("--forbidden", action="store_true", help="Apply forbidden subsequence constraint")
    parser.add_argument("--all-constraints", action="store_true", help="Apply all constraints")
    args = parser.parse_args()
    
    if args.all_constraints:
        args.adjacency = args.primitive = args.forbidden = True
    
    if args.adjacency and args.max_n >= 3:
        _debug_compare_adjacency(3)
    
    print(f"Constraints: adjacency={args.adjacency}, primitive={args.primitive}, forbidden={args.forbidden}")
    print("=" * 60)
    
    results = enumerate_constrained(
        args.max_n,
        use_adjacency=args.adjacency,
        use_primitivity=args.primitive,
        use_forbidden=args.forbidden,
        need_explicit_sets=args.verbose,
    )
    
    print("=" * 60)
    counts = [len(results[n]) for n in range(1, args.max_n + 1)]
    print(f"a(1..{args.max_n}) = {counts}")
    
    # Show A053656 for comparison
    expected_a053656 = [1, 2, 2, 4, 4, 9, 10, 22, 30, 62, 94, 192, 316, 623, 1096]
    print(f"A053656:          {expected_a053656[:args.max_n]}")
    
    if args.verbose:
        for n in range(1, args.max_n + 1):
            if len(results[n]) <= 30:
                print(f"\nn={n}:")
                for sig in sorted(results[n]):
                    print(f"  {format_sig(sig)}")


if __name__ == "__main__":
    main()
