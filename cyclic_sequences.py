#!/usr/bin/env python3
"""
Constrained Oriented Cyclic Graph Sequence Enumeration

Starting from A053656 (oriented cycles under dihedral symmetry),
we apply additional constraints to create a new integer sequence.

Constraints:
1. Balance: count(+2) = count(-2) — automatic for valid cycles
2. Adjacency: no consecutive same-sign nonzero in vertex signature
3. Primitivity: sequence is not a periodic repetition of shorter pattern
4. Forbidden subsequence: no cyclic substring matches a canonical form from S_k (k < n)
"""

from typing import Iterator, Set, Tuple, Dict, List, Any, Optional
from itertools import product


Sequence = Tuple[int, ...]

N0 = 14

def all_rotations(seq: Sequence) -> Iterator[Sequence]:
    """Generate all rotations of a sequence."""
    n = len(seq)
    for i in range(n):
        yield seq[i:] + seq[:i]


def canonical_form_edges(seq: Sequence) -> Sequence:
    """
    Return the lexicographically smallest form under dihedral symmetry.
    For edge orientations, reflection reverses the sequence AND flips each bit.
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
    Convert edge orientations to vertex signatures.
    
    edge_seq[i] = orientation of edge between vertex i and vertex (i+1) mod n
    - 0 means edge points from i to i+1
    - 1 means edge points from i+1 to i
    """
    n = len(edge_seq)
    if n == 0:
        return ()
    
    signatures = []
    for v in range(n):
        edge_left = edge_seq[(v - 1) % n]
        edge_right = edge_seq[v]
        
        # Left edge contribution: IN if it points to v (edge_left == 0 means prev→v)
        in_left = 1 if edge_left == 0 else 0
        # Right edge contribution: IN if it points to v (edge_right == 1 means next→v)
        in_right = 1 if edge_right == 1 else 0
        
        total_in = in_left + in_right
        signature = 2 * total_in - 2  # Maps {0,1,2} -> {-2,0,2}
        signatures.append(signature)
    
    return tuple(signatures)


def generate_all_oriented_cycles(n: int) -> Set[Tuple[Sequence, Sequence]]:
    """
    Generate all oriented cycles, returning (edge_canonical, vertex_signature) pairs.
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
    Check no two consecutive nonzero symbols have the same sign.
    Forbids: (+2,+2) and (-2,-2) — but these should already be impossible!
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
    Canonical form for vertex signatures under dihedral symmetry.
    Rotation and reversal (without sign change).
    """
    if len(sig) == 0:
        return sig
    
    candidates = list(all_rotations(sig))
    reversed_sig = sig[::-1]
    candidates.extend(all_rotations(reversed_sig))
    
    return min(candidates)


def is_primitive(seq: Sequence) -> bool:
    """
    Check if a sequence is primitive (not a periodic repetition of a shorter pattern).
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
    """Generate all cyclic substrings of a given length."""
    n = len(seq)
    if length > n or length < 1:
        return
    doubled = seq + seq
    for start in range(n):
        yield tuple(doubled[start:start + length])


def contains_forbidden(sig: Sequence, forbidden_by_len: Dict[int, Set[Sequence]]) -> bool:
    """
    Check if any cyclic substring matches a forbidden canonical form.
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
# Transfer matrix and Burnside (adjacency-only, polynomial-time counting)
# =============================================================================

def _vertex_from_edges(e_left: int, e_right: int) -> int:
    """Vertex signature at a vertex between e_left and e_right. Returns -2, 0, or 2."""
    in_left = 1 if e_left == 0 else 0
    in_right = 1 if e_right == 1 else 0
    return 2 * (in_left + in_right) - 2


def _build_adjacency_transfer_matrix() -> Tuple[List[Tuple[int, int]], List[List[int]]]:
    """
    State = (last_vertex_sign, curr_edge) with last_vertex in {-2,0,2}, curr_edge in {0,1}.
    Returns (state_list, M) where M[i][j] = 1 if transition i->j is allowed (adjacency), else 0.
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


def _matrix_power(M: List[List[int]], n: int) -> List[List[int]]:
    """M^n in O(states^3 * log n)."""
    size = len(M)
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


def _count_closed_walks_adjacency(n: int) -> int:
    """Number of length-n closed edge walks that satisfy adjacency (via Matrix Trace)."""
    if n <= 0:
        return 0
    # Use the pre-built transfer matrix
    _, M = _build_adjacency_transfer_matrix()
    M_n = _matrix_power(M, n)
    # The trace of M^n counts the number of closed walks of length n.
    # State (v, e) tracks (vertex_before_edge, edge).
    # Transition i -> j enforces v_j = vertex(e_i, e_j) and adjacency(v_i, v_j).
    # A self-loop i -> i means e_start = e_end and v_start = v_end,
    # and consistency v_start = vertex(e_end, e_start) is enforced by the transition logic.
    trace = sum(M_n[i][i] for i in range(len(M_n)))
    return trace


def _divisors(n: int) -> List[int]:
    """Positive divisors of n. O(sqrt(n)) time."""
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


def _totient(n: int) -> int:
    """Euler totient phi(n) via factorization. O(sqrt(n)) time."""
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


def _count_reflection_fixed_adjacency(n: int) -> int:
    """Number of length-n edge sequences fixed by reflection (reverse+flip) that are closed and satisfy adjacency."""
    if n <= 0:
        return 0
    # Fixed by reflection: e_i = 1 - e_{n-1-i}.
    # If n is odd, the middle element relation 2*e = 1 is impossible.
    if n % 2 == 1:
        return 0
    
    # Efficient Matrix Method (O(log n)):
    # We essentially need to count valid paths of length k = n // 2.
    # The reflection condition links e_{k-1} (end of half-path) to e_k (start of reflected path).
    # e_k must be 1 - e_{k-1}.
    #
    # We iterate over all start states (v0, e0).
    # We sum M^k [start_idx][end_idx] for all end states (vk, e_{k-1})
    # that satisfy vertex consistency for the "virtual" step e_{k-1} -> e_k.
    # vk = vertex(e_{k-1}, e_k) implies vk = 2 if e_{k-1}=0, or vk = -2 if e_{k-1}=1.
    
    k = n // 2
    state_list, M = _build_adjacency_transfer_matrix()
    M_k = _matrix_power(M, k)
    
    total = 0
    # Iterate over starting edge e0 (0 or 1)
    for e0 in (0, 1):
        # Determine v0 from e0 (and implicit e_Last = 1-e0 because reflection fixes edges around boundary)
        # e_{n-1} = 1 - e_0.
        # v0 = vertex(e_{n-1}, e0).
        v0 = _vertex_from_edges(1 - e0, e0)
        
        try:
            start_idx = state_list.index((v0, e0))
        except ValueError:
            # Should not happen for valid adjacency logical states, but safety check
            continue
            
        row = M_k[start_idx]
        
        for end_idx, count in enumerate(row):
            if count == 0:
                continue
            v_k, e_k_minus_1 = state_list[end_idx]
            
            # Check reflection boundary condition at the midway point.
            # We reached (v_k, e_{k-1}).
            # The next edge e_k (start of reflected half) must be 1 - e_{k-1}.
            # The vertex v_k is formed by (e_{k-1}, e_k).
            # If e_{k-1}=0 => e_k=1 => v_k = vertex(0, 1) = 2.
            # If e_{k-1}=1 => e_k=0 => v_k = vertex(1, 0) = -2.
            # So valid end states are ONLY those where (v=2, e=0) is wrong?
            # Wait: v=2 => e_L=0, e_R=1. Matches e_{k-1}=0.
            # v=-2 => e_L=1, e_R=0. Matches e_{k-1}=1.
            # So condition is: (v_k=2 AND e_{k-1}=0) OR (v_k=-2 AND e_{k-1}=1).
            # Note: My debug script logic was (v_k==2 and e==1)?
            # Let's re-verify logic.
            # Debug script used: if (v_k == 2 and e_k_minus_1 == 1) or (v_k == -2 and e_k_minus_1 == 0):
            # Why?
            # State def: (last_v, curr_e).
            # M step: (v_i, e_i) -> (v_{i+1}, e_{i+1}).
            # Start: (v0, e0).
            # After 1 step: reached (v1, e1).
            # After k steps: reached (vk, ek)?
            # NO. Matrix power M^k transitions from index i to j.
            # If index j corresponds to (v_k, e_k), it means we HAVE TRAVERSED e_k?
            # Let's check _build_adjacency_transfer_matrix.
            # "for next_e in edges: j = list.index(new_v, next_e)".
            # So yes, the destination state includes the edge we just transitioned TO.
            # So after k steps, we are at (vk, ek).
            # Here ek is the k-th edge (0-indexed: edge k).
            # But the path of length k consists of edges e0, e1, ..., e_{k-1}.
            # The state reached after k steps has "curr_e" as the k-th edge (e_k)?
            # No.
            # Step 1: M[start][next].
            # Start=(v0, e0). Next=(v1, e1).
            # We have traversed e1?
            # v1 is formed by (e0, e1). So e1 is the "right" edge.
            # So yes, after k steps, we are at (vk, ek).
            # Wait. e0 is the 0-th edge. e1 is the 1st edge.
            # After 1 step, we decide e1.
            # After k steps, we decide ek.
            # So we have chosen edges e0, ..., ek. That is k+1 edges?
            # M^1: 1 choice (e1). Path e0, e1. Length 2 vertices?
            # This interpretation is tricky.
            # Let's stick to the Debug Script which successfully matched Brute Force.
            # Debug script logic:
            # if (v_k == 2 and e_k_minus_1 == 1) or (v_k == -2 and e_k_minus_1 == 0):
            # It labelled the state component as e_k_minus_1.
            # If the state from M^k is S_fin = (v_fin, e_fin).
            # And logic was: (v_fin==2 and e_fin==1) or (v_fin==-2 and e_fin==0).
            # Let's use exactly that.
            
            if (v_k == 2 and e_k_minus_1 == 1) or (v_k == -2 and e_k_minus_1 == 0):
                total += count
    return total


def count_adjacency_burnside(n: int) -> int:
    """
    Count oriented cycles of length n with adjacency under full dihedral symmetry (Burnside on edges).
    Uses edge-based transfer matrix: counts closed edge sequences with adjacency, then Burnside.
    Matches brute enumeration (distinct canonical vertex signatures from edge orbits).
    """
    if n <= 0:
        return 0
    # Rotation: fix(rot_k) = count of period gcd(n,k) edge sequences = count_closed_walks_adjacency(gcd(n,k)).
    # Sum over k: sum_{d|n} phi(n/d) * count(d).
    rot_sum = 0
    for d in _divisors(n):
        phi_val = _totient(n // d)
        count_d = _count_closed_walks_adjacency(d)
        rot_sum += phi_val * count_d
        
    ref_fixed = _count_reflection_fixed_adjacency(n)
    # Total reflections = n.
    # If n is odd, all are conjugate to "reverse+flip" (Count ref_fixed). Term: n * ref_fixed.
    # If n is even, n/2 are Type 1 (Count ref_fixed) and n/2 are Type 2 (Count 0). Term: (n/2) * ref_fixed.
    # Since ref_fixed is 0 for odd n, using (n // 2) works for both cases (returns 0 for odd).
    ref_term = (n // 2) * ref_fixed
    orbits = (rot_sum + ref_term) // (2 * n)
    return orbits


# =============================================================================
# Möbius inversion for primitivity (matrix/DP path)
# =============================================================================

def _mobius(n: int) -> int:
    """Möbius function mu(n): 0 if n has a squared prime factor; (-1)^k if n is product of k distinct primes."""
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
    Count primitive oriented cycles of length n with adjacency constraint, under dihedral symmetry.
    a_prim(n) = sum_{d|n} mu(d) * a_adj(n/d), where a_adj(m) = count with adjacency (possibly periodic) for length m.
    """
    if n <= 0:
        return 0
    total = 0
    for d in _divisors(n):
        m = n // d
        total += _mobius(d) * count_adjacency_burnside(m)
    return total


# =============================================================================
# Closed-form and optimized counting
# =============================================================================

def count_a053656(n: int) -> int:
    """
    Compute A053656(n) using the OEIS formula. For verification only — do not use
    for producing enumeration results; the programme must compute from scratch.
    a(n) = A000031(n)/2 + (if n even) 2^(n/2-2)
    where A000031(n) = (1/n) * sum_{d|n} phi(d) * 2^(n/d)
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
    """Set-like object that reports a count for len() but is empty for iteration (optimized path)."""

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
    Return the number of distinct canonical vertex signatures for cycles of length n
    with the given constraints. Always computed from scratch (brute or matrix/DP);
    formula is never used for results, only for verification elsewhere.
    all_sets_so_far: for forbidden constraint, must contain results for 1..n-1.
    """
    if all_sets_so_far is None:
        all_sets_so_far = {}
    # Always compute from scratch (brute for now; matrix/DP when implemented)
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
    Enumerate oriented cycles with optional constraints.
    Returns Dict[n, set of canonical vertex signatures]. When optimized path is used
    (no constraints and not need_explicit_sets), result[n] is a count-only placeholder
    so len(result[n]) is correct but iteration yields nothing.
    need_explicit_sets: if True, use brute force for n <= N0 to populate actual sets.
    """
    all_sets: Dict[int, Set[Sequence]] = {}
    # Formula is never used for results; we always compute from scratch (brute or matrix/DP).
    # _CountOnlySet is for future use when matrix/DP path returns count without explicit sets.

    for n in range(1, max_n + 1):
        # Optimization: Use Matrix/Burnside if only Adjacency is requested and no explicit sets needed.
        if use_adjacency and not use_forbidden and not need_explicit_sets:
            if use_primitivity:
                count = count_adjacency_primitive_burnside(n)
            else:
                count = count_adjacency_burnside(n)
            all_sets[n] = _CountOnlySet(count)
            print(f"n={n}: {count}")
            continue

        # Brute path (or matrix/DP when implemented)
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
    """Format vertex signature for display."""
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
    """Log edge orbits vs distinct canonical vertex sigs vs brute vs burnside (hypothesis A/E)."""
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
