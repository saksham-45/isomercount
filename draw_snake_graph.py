#!/usr/bin/env python3
"""
Enumerate snake graph shapes and distinct 4-edged colored snakes (efficient counting only).

Enumerates all distinct geometric arrangements (shapes) of L square cells on a grid
(self-avoiding snake polyominoes), and for each shape counts distinct edge colorings
up to the 32-element symmetry (8 dihedral × path reversal × global 0↔1 flip).

WHAT WE COUNT:
  - Each cell is a unit square with 4 edges (S,E,N,W), each colored 0 or 1.
  - Valid patterns: sum(edges) <= 2 per cell (11 patterns: all-0, four with one 1,
    six with two 1s). Consecutive cells share exactly one edge (hinge); values match.
  - Equivalence: two colorings are the same if related by (dihedral × reversal × flip).
  - For linear (L cells in a row): 1 shape; counts verified via brute-force for L=1,2,3.

Usage:
    python draw_snake_graph.py --table 10     # table: L, Shapes (A002013), Total for L=1..10
    python draw_snake_graph.py --table 10 --linear  # linear only: 1 shape per L
    python draw_snake_graph.py 5              # shapes count for L=5
    python draw_snake_graph.py 5 --linear     # linear only: 1 shape for L=5
    python draw_snake_graph.py 5 --full        # per-shape and total distinct snakes for L=5
    python draw_snake_graph.py 5 --all-shapes # include non-snake shapes (U-bends)
    python draw_snake_graph.py 3 --draw       # draw all unique snakes for L=3
    python draw_snake_graph.py 2 --draw --linear --save snakes_linear_L2.png
"""

import math
import argparse
import sys
import concurrent.futures
from typing import List, Tuple, Set, Dict, Optional

# Directions: S=0, E=1, N=2, W=3
S, E, N, W = 0, 1, 2, 3
DIR_DELTA = {S: (0, -1), E: (1, 0), N: (0, 1), W: (-1, 0)}
OPP = {S: N, E: W, N: S, W: E}

# ===================================================================
# 8 dihedral transformations of the plane (positions)
# ===================================================================
_DIHEDRAL = [
    lambda x, y: (x, y),        # 0  identity
    lambda x, y: (-y, x),       # 1  rot 90 CCW
    lambda x, y: (-x, -y),      # 2  rot 180
    lambda x, y: (y, -x),       # 3  rot 270 CCW
    lambda x, y: (x, -y),       # 4  reflect across x-axis
    lambda x, y: (-x, y),       # 5  reflect across y-axis
    lambda x, y: (y, x),        # 6  reflect across y=x
    lambda x, y: (-y, -x),      # 7  reflect across y=-x
]

# ===================================================================
# Edge direction permutations for each dihedral transform.
# After applying dihedral transform T to the plane, a cell's edge that
# was in world direction d now sits in direction d'. The permutation is:
#   new_edges[d'] = old_edges[d]   i.e.   new_edges[d] = old_edges[perm[d]]
#
# These were derived by tracking unit-square vertices through each transform.
# ===================================================================
_EDGE_PERMS = (
    (0, 1, 2, 3),  # 0  identity
    (3, 0, 1, 2),  # 1  rot 90 CCW:  S←W  E←S  N←E  W←N
    (2, 3, 0, 1),  # 2  rot 180:     S←N  E←W  N←S  W←E
    (1, 2, 3, 0),  # 3  rot 270 CCW: S←E  E←N  N←W  W←S
    (2, 1, 0, 3),  # 4  reflect x:   S↔N  E←E  W←W
    (0, 3, 2, 1),  # 5  reflect y:   S←S  E↔W  N←N
    (3, 2, 1, 0),  # 6  reflect y=x: S←W  E←N  N←E  W←S
    (1, 0, 3, 2),  # 7  reflect y=-x:S←E  E←S  N←W  W←N
)


# ===================================================================
# Self-avoiding walk enumeration
# ===================================================================

def _enumerate_walks(L: int) -> List[Tuple[Tuple[int, int], ...]]:
    """All self-avoiding walks of L cells starting at (0,0)."""
    if L <= 0:
        return []
    if L == 1:
        return [((0, 0),)]
    result: list = []

    def dfs(path, visited):
        if len(path) == L:
            result.append(tuple(path))
            return
        x, y = path[-1]
        for d in (S, E, N, W):
            dx, dy = DIR_DELTA[d]
            npos = (x + dx, y + dy)
            if npos not in visited:
                path.append(npos)
                visited.add(npos)
                dfs(path, visited)
                path.pop()
                visited.discard(npos)

    dfs([(0, 0)], {(0, 0)})
    return result


def _enumerate_snake_shapes_fast(L: int) -> List[Tuple[Tuple[int, int], ...]]:
    """
    Enumerate 2-sided snake polyomino shapes efficiently.
    Uses: (1) early snake pruning - reject non-snake branches immediately,
    (2) symmetry breaking - try all 4 first steps, deduplicate by canonical form.
    """
    if L <= 0:
        return []
    if L == 1:
        return [((0, 0),)]
    seen: Dict[Tuple, Tuple] = {}

    def dfs(path: List[Tuple[int, int]], visited: Set[Tuple[int, int]],
           tip: Tuple[int, int], prev: Tuple[int, int]) -> None:
        if len(path) == L:
            c = _canonical_shape(tuple(path))
            if c not in seen:
                seen[c] = tuple(path)
            return
        tx, ty = tip
        for d in (S, E, N, W):
            dx, dy = DIR_DELTA[d]
            npos = (tx + dx, ty + dy)
            if npos == prev or npos in visited:
                continue
            # Snake constraint: new cell's neighbors (except tip) must not be in path
            for nd in (S, E, N, W):
                ndx, ndy = DIR_DELTA[nd]
                nn = (npos[0] + ndx, npos[1] + ndy)
                if nn != tip and nn in visited:
                    break  # would create non-consecutive adjacency
            else:
                path.append(npos)
                visited.add(npos)
                dfs(path, visited, npos, tip)
                visited.discard(npos)
                path.pop()

    # Try all 4 first steps (E,N,W,S) - early snake pruning avoids most dead branches
    for dx, dy in [(1, 0), (0, 1), (-1, 0), (0, -1)]:
        path = [(0, 0), (dx, dy)]
        visited = {(0, 0), (dx, dy)}
        dfs(path, visited, (dx, dy), (0, 0))

    return [seen[c] for c in sorted(seen.keys())]


# ===================================================================
# Shape canonicalization (dihedral + path reversal)
# ===================================================================

def _normalize_walk(walk):
    """Translate so bounding box lower-left is at (0,0)."""
    mx = min(p[0] for p in walk)
    my = min(p[1] for p in walk)
    return tuple((x - mx, y - my) for x, y in walk)


def _canonical_shape(walk):
    """Canonical form under 8 dihedral × 2 reversal = 16 variants."""
    best = None
    for tf in _DIHEDRAL:
        tw = _normalize_walk(tuple(tf(x, y) for x, y in walk))
        rw = _normalize_walk(tuple(reversed(tw)))
        for w in (tw, rw):
            if best is None or w < best:
                best = w
    return best


def _is_snake_polyomino(walk) -> bool:
    """
    A snake polyomino requires that NO two non-consecutive cells in the path
    are grid-adjacent. E.g. a U-shape where cell 0 and cell 3 share a grid
    edge is NOT a snake polyomino (it would create a "hole" / extra adjacency).
    """
    L = len(walk)
    if L <= 2:
        return True
    pos_to_idx = {pos: i for i, pos in enumerate(walk)}
    for i, (x, y) in enumerate(walk):
        for d in (S, E, N, W):
            dx, dy = DIR_DELTA[d]
            nb = (x + dx, y + dy)
            if nb in pos_to_idx:
                j = pos_to_idx[nb]
                if abs(i - j) > 1:
                    return False
    return True


def enumerate_linear_shapes(L: int) -> List[Tuple[Tuple[int, int], ...]]:
    """
    Return the single linear (straight-line) shape: L cells in a row.
    Canonical form: ((0,0), (1,0), (2,0), ..., (L-1, 0)).
    """
    if L <= 0:
        return []
    return [tuple((i, 0) for i in range(L))]


def enumerate_shapes(L: int, snake_only: bool = True) -> List[Tuple[Tuple[int, int], ...]]:
    """
    Return one representative walk per distinct shape.
    If snake_only=True (default), only return shapes that are valid snake
    polyominoes (no non-consecutive grid adjacencies).
    Uses fast path with early pruning and symmetry breaking when snake_only=True.
    """
    if snake_only:
        return _enumerate_snake_shapes_fast(L)
    walks = _enumerate_walks(L)
    seen: Dict[Tuple, Tuple] = {}
    for walk in walks:
        if not _is_snake_polyomino(walk):
            continue
        c = _canonical_shape(walk)
        if c not in seen:
            seen[c] = walk
    return [seen[c] for c in sorted(seen.keys())]


# ===================================================================
# Shared edges between consecutive cells
# ===================================================================

def _shared_edges(walk):
    """Set of (cell_pos, direction) pairs that are shared (hinge) edges."""
    shared: set = set()
    for t in range(len(walk) - 1):
        dx = walk[t + 1][0] - walk[t][0]
        dy = walk[t + 1][1] - walk[t][1]
        for d in (S, E, N, W):
            if DIR_DELTA[d] == (dx, dy):
                shared.add((walk[t], d))
                shared.add((walk[t + 1], OPP[d]))
                break
    return shared


# ===================================================================
# Valid edge patterns (from cell types, n=4)
# ===================================================================
# For n=4, the 4 canonical cell types are:
#   (0,0,0,0), (0,0,0,1), (0,0,1,1), (0,1,0,1)
# Each can be rotated 4 ways, giving 1+4+4+2 = 11 distinct patterns.
# Equivalently: a pattern is valid iff sum(edges) <= 2.
# NOTE: We use sum <= 2 (not sum == 2). For "balanced" (sum == 2) only, use
# both_only=True in _fixed_count_burnside; that restricts to 6 patterns.

_VALID_PATTERNS_4: Tuple[Tuple[int, ...], ...] = tuple(
    (a, b, c, d)
    for a in (0, 1) for b in (0, 1) for c in (0, 1) for d in (0, 1)
    if a + b + c + d <= 2
)
# Precompute lookup by fixed edge: _VALID_BY_EDGE[dir][val] = patterns with edges[dir]==val
_VALID_BY_EDGE: Dict[int, Dict[int, List[Tuple[int, ...]]]] = {}
for _d in range(4):
    _VALID_BY_EDGE[_d] = {0: [], 1: []}
    for _p in _VALID_PATTERNS_4:
        _VALID_BY_EDGE[_d][_p[_d]].append(_p)


# ===================================================================
# Edge coloring enumeration (valid cell types, matching at hinges)
# ===================================================================

def _build_walk_adj(walk):
    """For each cell t>0, return (direction, predecessor_index) for matching."""
    L = len(walk)
    adj: List[List[Tuple[int, int]]] = [[] for _ in range(L)]
    for t in range(L - 1):
        dx = walk[t + 1][0] - walk[t][0]
        dy = walk[t + 1][1] - walk[t][1]
        for d in (S, E, N, W):
            if DIR_DELTA[d] == (dx, dy):
                adj[t + 1].append((OPP[d], t))  # cell t+1's entry edge
                break
    return adj


def _exit_dirs(walk):
    """For each cell t (0..L-2), return direction from cell t to cell t+1."""
    result = []
    for t in range(len(walk) - 1):
        dx = walk[t + 1][0] - walk[t][0]
        dy = walk[t + 1][1] - walk[t][1]
        for d in (S, E, N, W):
            if DIR_DELTA[d] == (dx, dy):
                result.append(d)
                break
    return result


# Precompute transfer blocks: (d_prev, d_next) -> 2x2 count matrix (v_prev, v_next)
# Count of patterns p with p[OPP[d_prev]]==v_prev and p[d_next]==v_next
_TRANSFER_BLOCK: Dict[Tuple[int, int], List[List[int]]] = {}
for _d_prev in range(4):
    _entry = OPP[_d_prev]
    for _d_next in range(4):
        _bl = [[0, 0], [0, 0]]
        for _p in _VALID_PATTERNS_4:
            vp = _p[_entry]
            vn = _p[_d_next]
            _bl[vp][vn] += 1
        _TRANSFER_BLOCK[(_d_prev, _d_next)] = _bl


def count_colorings_transfer_generic(walk: Tuple, valid_patterns: List[Tuple[int, ...]]) -> int:
    """
    Count valid edge colorings for this walk using transfer matrix with given valid patterns.
    """
    L = len(walk)
    if L == 0:
        return 0
    if L == 1:
        return len(valid_patterns)

    # Precompute transfer blocks for this specific valid_patterns set
    v_by_edge = {}
    for _d in range(4):
        v_by_edge[_d] = {0: 0, 1: 0}
        for _p in valid_patterns:
            v_by_edge[_d][_p[_d]] += 1

    trans_block = {}
    for _d_prev in range(4):
        _entry = OPP[_d_prev]
        for _d_next in range(4):
            _bl = [[0, 0], [0, 0]]
            for _p in valid_patterns:
                vp = _p[_entry]
                vn = _p[_d_next]
                _bl[vp][vn] += 1
            trans_block[(_d_prev, _d_next)] = _bl

    exit_dirs = _exit_dirs(walk)
    # Initial: cell 0, exit direction exit_dirs[0]
    d0 = exit_dirs[0]
    init = [0] * 8
    for p in valid_patterns:
        v = p[d0]
        init[d0 * 2 + v] += 1

    vec = init[:]
    for t in range(L - 2):
        d_prev = exit_dirs[t]
        d_next = exit_dirs[t + 1]
        bl = trans_block[(d_prev, d_next)]
        new_vec = [0] * 8
        for vp in (0, 1):
            if vec[d_prev * 2 + vp] == 0: continue
            for vn in (0, 1):
                cnt = bl[vp][vn]
                if cnt:
                    new_vec[d_next * 2 + vn] += vec[d_prev * 2 + vp] * cnt
        vec = new_vec

    # Last cell L-1
    d_last = exit_dirs[L - 2]
    entry_dir = OPP[d_last]
    total = 0
    for v in (0, 1):
        total += vec[d_last * 2 + v] * v_by_edge[entry_dir][v]
    return total

def count_colorings_transfer(walk: Tuple) -> int:
    """Count valid edge colorings for this walk using n=4 transfer matrix."""
    return count_colorings_transfer_generic(walk, _VALID_PATTERNS_4)


def _stabilizer_elements(walk):
    """
    Yield group elements that map this walk to its canonical form (for enumeration).
    Each element is (di, rev, perm, flip) with perm = _EDGE_PERMS[di].
    """
    canon_walk, autos = _compute_shape_autos(walk)
    for (di, rev, perm) in autos:
        for flip in (False, True):
            yield (di, rev, perm, flip)


def _stabilizer_elements_burnside(walk):
    """
    Yield all 32 group elements g that preserve the SHAPE (set of positions).
    Used for Burnside: we need Stab(shape) = {g : set(normalize(g(walk))) == set(normalize(walk))}.
    Each element is (di, rev, perm, flip).
    """
    L = len(walk)
    walk_set = set(_normalize_walk(walk))
    for di in range(8):
        tf = _DIHEDRAL[di]
        tw = tuple(tf(x, y) for x, y in walk)
        perm = _EDGE_PERMS[di]
        for rev in (False, True):
            w = tw if not rev else tuple(reversed(tw))
            nw = _normalize_walk(w)
            if set(nw) == walk_set:
                for flip in (False, True):
                    yield (di, rev, perm, flip)


def _apply_group_element(coloring, di, rev, perm, flip):
    """Apply symmetry (di, rev, perm, flip) to a coloring; returns new coloring tuple."""
    tc = tuple(
        (edges[perm[0]], edges[perm[1]], edges[perm[2]], edges[perm[3]])
        for edges in coloring
    )
    if rev:
        tc = tuple(reversed(tc))
    if flip:
        tc = tuple((1 - e[0], 1 - e[1], 1 - e[2], 1 - e[3]) for e in tc)
    return tc


def _cell_permutation(walk, di, rev):
    """
    Return permutation sigma: sigma[i] = old index k such that the cell at
    normalized transformed walk index i is the image of original cell k.
    """
    L = len(walk)
    tf = _DIHEDRAL[di]
    # Geometric image of each cell:
    image = [tf(x, y) for x, y in walk]
    if rev:
        image = list(reversed(image))
    
    # Normalize the image walk
    mx = min(p[0] for p in image)
    my = min(p[1] for p in image)
    nw = [ (x - mx, y - my) for x, y in image ]
    
    # sigma[i] is which old cell k ended up at nw[i].
    # NW[i] came from image[i]. 
    # image[i] is walk[i] transformed if not rev, else walk[L-1-i] transformed.
    # So if nw[i] == walk[i] for all i (which must hold if g is in Stab(walk)),
    # then the cell at position walk[i] in the new walk is old cell i if not rev, else L-1-i.
    return [i if not rev else L - 1 - i for i in range(L)]


def _fixed_count_burnside(walk, di, rev, perm, flip, both_only=False) -> int:
    """
    Count colorings c of `walk` such that g(c) = c (or fg(c) = c if flip=True).
    If both_only=True, only count colorings where all cells have sum == 2.
    """
    L = len(walk)
    
    # Optimization: if g is identity, use transfer matrix
    if di == 0 and not rev and not flip:
        v_set = _VALID_PATTERNS_4 if not both_only else [p for p in _VALID_PATTERNS_4 if sum(p) == 2]
        return count_colorings_transfer_generic(walk, v_set)

    sigma = _cell_permutation(walk, di, rev)

    def transform_cell(pat):
        t = (pat[perm[0]], pat[perm[1]], pat[perm[2]], pat[perm[3]])
        if flip:
            t = (1 - t[0], 1 - t[1], 1 - t[2], 1 - t[3])
        return t

    # Cycles of sigma
    seen = [False] * L
    cycles = []
    for start in range(L):
        if seen[start]:
            continue
        cyc = []
        i = start
        while not seen[i]:
            seen[i] = True
            cyc.append(i)
            i = sigma[i]
        cycles.append(cyc)

    # Patterns preserved by k applications of transform_cell
    valid_range = _VALID_PATTERNS_4 if not both_only else [p for p in _VALID_PATTERNS_4 if sum(p) == 2]

    def fix_set(k):
        result = []
        for p in valid_range:
            q = p
            valid_path = True
            for _ in range(k):
                q = transform_cell(q)
                if q not in valid_range:
                    valid_path = False
                    break
            if valid_path and q == p:
                result.append(p)
        return result

    fix_sets = [fix_set(len(c)) for c in cycles]
    cell_to_cyc_idx = {}
    for ci, cyc in enumerate(cycles):
        for cell in cyc:
            cell_to_cyc_idx[cell] = ci

    # Hinge adjacency
    hinges = []
    for t in range(L - 1):
        dx = walk[t + 1][0] - walk[t][0]
        dy = walk[t + 1][1] - walk[t][1]
        for d in (S, E, N, W):
            if DIR_DELTA[d] == (dx, dy):
                hinges.append((t, t + 1, d, OPP[d]))
                break

    cyc_assignment = [None] * len(cycles)
    
    # Pre-get patterns for a cycle assignment to avoid repeated work
    cyc_patterns = [None] * len(cycles)

    def get_pattern_local(cell, ci, pat):
        # find power of cell in cycle ci
        cyc = cycles[ci]
        pw = cyc.index(cell)
        p = pat
        for _ in range(pw):
            p = transform_cell(p)
        return p

    def count_assignments(c_idx):
        if c_idx == len(cycles):
            return 1
        
        c = 0
        for pat in fix_sets[c_idx]:
            # Before accepting pat, check all hinges between cells in this cycle
            # and already assigned cells.
            cyc_patterns[c_idx] = [get_pattern_local(cell, c_idx, pat) for cell in cycles[c_idx]]
            
            valid = True
            # For each cell in the current cycle
            for cell_in_cyc_local_idx, cell_idx in enumerate(cycles[c_idx]):
                p_current = cyc_patterns[c_idx][cell_in_cyc_local_idx]
                # check hinges involving cell_idx
                for t1, t2, d1, d2 in hinges:
                    if t1 == cell_idx:
                        other = t2
                        d_self, d_other = d1, d2
                    elif t2 == cell_idx:
                        other = t1
                        d_self, d_other = d2, d1
                    else: continue
                    
                    # check if 'other' is already assigned
                    other_cyc = cell_to_cyc_idx[other]
                    if other_cyc < c_idx:
                        # find local index of 'other' in its cycle
                        other_local_idx = cycles[other_cyc].index(other)
                        p_other = cyc_patterns[other_cyc][other_local_idx]
                        if p_current[d_self] != p_other[d_other]:
                            valid = False
                            break
                    elif other_cyc == c_idx:
                        # match within the same cycle
                        other_local_idx = cycles[c_idx].index(other)
                        p_other = cyc_patterns[c_idx][other_local_idx]
                        if p_current[d_self] != p_other[d_other]:
                            valid = False
                            break
                if not valid: break
            
            if valid:
                c += count_assignments(c_idx + 1)
        
        cyc_assignment[c_idx] = None
        cyc_patterns[c_idx] = None
        return c

    return count_assignments(0)


def count_orbits_burnside(walk) -> int:
    """
    Count orbits of valid colorings under the stabilizer of the walk × Flip
    using a modified Burnside formula.

    Group: G = G0 × {id, flip}, where G0 = stabilizer of the walk (geometric
    transforms that map the walk to its canonical form). The stabilizer is a
    subgroup of (8 dihedral × 2 reversal); flip is global 0↔1 on all edges.

    Formula: N = (2*sum_V0 - sum_Both0 + sum_BothFlip) / (2*|G0|)
    where sum_V0 = sum_g |Fix_V(g)|, Fix_V = all valid colorings (sum ≤ 2),
    sum_Both0 = sum_g |Fix_Vboth(g)|, Fix_Vboth = colorings with all cells sum=2,
    sum_BothFlip = sum_g |Fix_Vboth(f·g)| ( flip applied before g's fixed count).
    The "both" terms handle the flip involution in the orbit count.

    Matches the enumeration's "distinct canonical forms".
    Verified: brute-force orbit count agrees for L=1,2,3 (see test_linear_snake.py).
    """
    # 1. Standardize the walk: perform Burnside on the canonical representative.
    # This ensures the geometric stabilizer G0 always includes the identity.
    canon_walk, _ = _compute_shape_autos(walk)
    _, autos = _compute_shape_autos(canon_walk)
    G0 = [(di, rev, perm) for (di, rev, perm) in autos]

    # Formula: N = N0 - 0.5 * N_both_0 + 0.5 * N_both_flip
    # where N0 = sum |Fix_V(g)| / |G0|
    # N_both_0 = sum |Fix_Vboth(g)| / |G0|
    # N_both_flip = sum |Fix_Vboth(fg)| / |G0|

    sum_V0 = 0
    sum_Both0 = 0
    sum_BothFlip = 0

    for di, rev, perm in G0:
        sum_V0 += _fixed_count_burnside(canon_walk, di, rev, perm, flip=False, both_only=False)
        sum_Both0 += _fixed_count_burnside(canon_walk, di, rev, perm, flip=False, both_only=True)
        sum_BothFlip += _fixed_count_burnside(canon_walk, di, rev, perm, flip=True, both_only=True)

    g0_size = len(G0)
    # The formula can be rearranged: (sum_V0 - 0.5*sum_Both0 + 0.5*sum_BothFlip) / g0_size
    # = (2*sum_V0 - sum_Both0 + sum_BothFlip) / (2 * g0_size)
    total_num = 2 * sum_V0 - sum_Both0 + sum_BothFlip
    total_den = 2 * g0_size
    assert total_num % total_den == 0, (
        f"Burnside formula: total_num={total_num} not divisible by total_den={total_den}"
    )
    return total_num // total_den


def _partial_can_be_canonical(assignment: list, t: int, L: int, autos) -> bool:
    """
    Return True iff the partial coloring (assignment[0..t]) can be the prefix
    of the canonical representative. We pad with (0,0,0,0) and check that
    the canonical form of the padded coloring starts with our partial.
    """
    padded = tuple(assignment[i] for i in range(t + 1)) + ((0, 0, 0, 0),) * (L - t - 1)
    canon = _canonical_coloring(padded, autos)
    return tuple(assignment[i] for i in range(t + 1)) == canon[: t + 1]


def _enumerate_valid_colorings(walk, autos=None, canonical_first=False) -> List[Tuple[Tuple[int, ...], ...]]:
    """
    Enumerate all valid edge colorings for a walk (n=4).
    Valid = each cell's 4 edge values form a pattern with sum <= 2,
    AND shared edges between consecutive cells match.
    If canonical_first=True and autos is provided, only enumerate colorings
    that are the canonical representative (prune during backtrack).
    """
    L = len(walk)
    adj = _build_walk_adj(walk)
    if canonical_first and autos is None:
        _, autos = _compute_shape_autos(walk)

    result: list = []
    assignment: list = [None] * L

    def backtrack(t):
        if t == L:
            result.append(tuple(assignment))
            return
        if adj[t]:
            entry_dir, pred = adj[t][0]
            required_val = assignment[pred][OPP[entry_dir]]
            candidates = _VALID_BY_EDGE[entry_dir][required_val]
        else:
            candidates = _VALID_PATTERNS_4
        for pat in candidates:
            assignment[t] = pat
            if canonical_first and autos is not None:
                if not _partial_can_be_canonical(assignment, t, L, autos):
                    assignment[t] = None
                    continue
            backtrack(t + 1)
            assignment[t] = None

    backtrack(0)
    return result


# ===================================================================
# Canonical form for (walk, coloring) under full 32-element symmetry
# ===================================================================

def _compute_shape_autos(walk):
    """
    Precompute the relevant transforms for a given walk.

    For each of 8 dihedral × 2 reversal = 16 geometric transforms, check
    which ones map the walk to its canonical shape. Only those transforms
    can produce the canonical (walk, coloring) — because the walk component
    of the canonical form is always the shape's canonical walk (the minimum).

    Returns:
        canon_walk: the canonical walk for this shape
        autos: list of (di_index, reverse_flag, edge_perm) tuples
               that map `walk` to `canon_walk` after normalization
    """
    canon_walk = _canonical_shape(walk)
    autos = []
    for di in range(8):
        tf = _DIHEDRAL[di]
        tw = tuple(tf(x, y) for x, y in walk)
        perm = _EDGE_PERMS[di]
        for rev in (False, True):
            w = tw if not rev else tuple(reversed(tw))
            nw = _normalize_walk(w)
            if nw == canon_walk:
                autos.append((di, rev, perm))
    return canon_walk, autos


def _canonical_coloring(coloring, autos):
    """
    Find the canonical coloring using precomputed automorphisms.

    Instead of iterating over all 32 transforms, only iterates over the
    transforms that map this walk to its canonical walk (typically 2-8),
    times 2 for global value flip.

    Returns the canonical coloring tuple.
    """
    best = None
    for di, rev, perm in autos:
        # Apply edge permutation
        tc = tuple(
            (edges[perm[0]], edges[perm[1]], edges[perm[2]], edges[perm[3]])
            for edges in coloring
        )
        # Apply reversal
        c = tc if not rev else tuple(reversed(tc))

        # No flip
        if best is None or c < best:
            best = c
        # With global value flip (0↔1)
        cf = tuple((1 - e[0], 1 - e[1], 1 - e[2], 1 - e[3]) for e in c)
        if cf < best:
            best = cf

    return best


def _canonical_snake(walk, coloring):
    """
    Canonical form of a colored snake under the 32-element symmetry group:
      8 dihedral transforms × 2 path reversal × 2 global value flip.
    Returns (canonical_walk, canonical_coloring) as a comparable tuple.
    (Slow path — used only when no precomputed autos are available.)
    """
    canon_walk, autos = _compute_shape_autos(walk)
    return (canon_walk, _canonical_coloring(coloring, autos))


# ===================================================================
# Full colored snake enumeration
# ===================================================================

def _count_orbits_one_shape(walk, autos, colorings=None):
    """Count distinct orbits for one shape (no list storage)."""
    if colorings is None:
        colorings = _enumerate_valid_colorings(walk)
    seen = set()
    for col in colorings:
        cc = _canonical_coloring(col, autos)
        if cc not in seen:
            seen.add(cc)
    return len(seen)


def enumerate_colored_snakes(L, n=4, snake_only=True, linear_only=False, count_only=False, max_snakes=None):
    """
    Enumerate all distinct (shape, coloring) pairs for L-cell snake graphs.
    Returns:
      - all_snakes: list of (walk, coloring) representatives (empty if count_only)
      - shapes: list of walks (one per shape)
      - counts_per_shape: list of int (distinct colorings per shape)
    If linear_only=True, only the straight-line shape is used (L cells in a row).
    If count_only=True, only counts are computed (no list of snakes stored).
    If max_snakes is set, stop after collecting that many (for drawing large L).
    """
    if n != 4:
        raise NotImplementedError("Only n=4 supported for geometric enumeration")

    shapes = enumerate_linear_shapes(L) if linear_only else enumerate_shapes(L, snake_only=snake_only)
    all_snakes: list = []
    counts_per_shape: list = []

    for walk in shapes:
        if max_snakes is not None and len(all_snakes) >= max_snakes:
            break
        _canon_walk, autos = _compute_shape_autos(walk)
        colorings = _enumerate_valid_colorings(walk)
        seen: set = set()
        shape_reps: list = []
        for col in colorings:
            cc = _canonical_coloring(col, autos)
            if cc not in seen:
                seen.add(cc)
                if not count_only:
                    shape_reps.append((walk, col))
                    if max_snakes is not None and len(all_snakes) + len(shape_reps) >= max_snakes:
                        break
        counts_per_shape.append(len(seen))
        if not count_only:
            all_snakes.extend(shape_reps)
            if max_snakes is not None and len(all_snakes) >= max_snakes:
                all_snakes = all_snakes[:max_snakes]
                break

    return all_snakes, shapes, counts_per_shape


def _count_orbits_for_parallel(walk):
    """Helper for parallel processing."""
    return count_orbits_burnside(walk)

def count_all_colored_snakes(L, n=4, snake_only=True, linear_only=False, parallel=False):
    """
    Count total distinct colored snakes for L cells without storing the list.
    Returns (total, shapes, counts_per_shape).
    Uses the optimized Burnside counting logic.
    If linear_only=True, only the straight-line shape is used (L cells in a row).
    """
    if n != 4:
        raise NotImplementedError("Only n=4 supported for geometric enumeration")
    
    shapes = enumerate_linear_shapes(L) if linear_only else enumerate_shapes(L, snake_only=snake_only)
    
    if parallel:
        with concurrent.futures.ProcessPoolExecutor() as executor:
            counts = list(executor.map(_count_orbits_for_parallel, shapes))
    else:
        counts = [count_orbits_burnside(walk) for walk in shapes]
        
    return sum(counts), shapes, counts


# ===================================================================
# Drawing (matplotlib)
# ===================================================================

def _draw_snake(walk: Tuple, coloring: Tuple[Tuple[int, ...], ...], ax, cell_size: float = 1.0,
                edge_width: float = 2.5, color0: str = "#3498db", color1: str = "#e74c3c") -> None:
    """
    Draw one snake: walk = list of (x,y) cell positions, coloring = (S,E,N,W) per cell.
    Edge value 0 uses color0, value 1 uses color1.
    """
    try:
        from matplotlib.patches import Rectangle
    except ImportError:
        raise ImportError("matplotlib required for drawing. Install with: pip install matplotlib")

    nw = _normalize_walk(walk)
    for t, (x, y) in enumerate(nw):
        pat = coloring[t]  # (S, E, N, W)
        # Cell fill (light gray)
        ax.add_patch(Rectangle((x * cell_size, y * cell_size), cell_size, cell_size,
                               fill=True, facecolor="#f8f9fa", edgecolor="none", zorder=0))
        # Four edges: S=bottom, E=right, N=top, W=left
        edges = [
            ((x, y), (x + 1, y), pat[S]),           # S
            ((x + 1, y), (x + 1, y + 1), pat[E]),   # E
            ((x, y + 1), (x + 1, y + 1), pat[N]),  # N
            ((x, y), (x, y + 1), pat[W]),          # W
        ]
        for p1, p2, val in edges:
            ax.plot(
                [p1[0] * cell_size, p2[0] * cell_size],
                [p1[1] * cell_size, p2[1] * cell_size],
                color=color0 if val == 0 else color1,
                linewidth=edge_width,
                solid_capstyle="butt",
                zorder=1,
            )
    ax.set_aspect("equal")
    ax.axis("off")


def draw_all_snakes(snakes: List[Tuple[Tuple, Tuple]], cols: int = 5,
                   save_path: Optional[str] = None, cell_size: float = 0.9) -> None:
    """Draw all unique snakes in a grid. snakes = [(walk, coloring), ...]"""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib required for drawing. Install with: pip install matplotlib", file=sys.stderr)
        sys.exit(1)

    n = len(snakes)
    if n == 0:
        return
    ncols = min(cols, n)
    nrows = (n + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=(3 * ncols, 3 * nrows), squeeze=False)
    axes_flat = axes.flatten()

    for idx, (walk, coloring) in enumerate(snakes):
        ax = axes_flat[idx]
        _draw_snake(walk, coloring, ax, cell_size=cell_size)

        # Compute bounds for this snake
        nw = _normalize_walk(walk)
        xs = [p[0] for p in nw]
        ys = [p[1] for p in nw]
        pad = 0.3
        ax.set_xlim(min(xs) - pad, max(xs) + 1 + pad)
        ax.set_ylim(min(ys) - pad, max(ys) + 1 + pad)
        ax.set_title(f"#{idx + 1}", fontsize=10)

    for j in range(n, len(axes_flat)):
        axes_flat[j].axis("off")

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=120, bbox_inches="tight")
        plt.close()
        print(f"Saved to {save_path}")
    else:
        plt.show()


# ===================================================================
# CLI
# ===================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Enumerate snake polyomino shapes and distinct 4-edged colored snakes (counts only)."
    )
    parser.add_argument("L", type=int, nargs="?", default=None,
                        help="Number of cells (optional if --table used)")
    parser.add_argument("n", type=int, nargs="?", default=4,
                        help="Edges per cell (only 4 supported)")
    parser.add_argument("--full", action="store_true",
                        help="Print per-shape and total distinct snakes for given L")
    parser.add_argument("--all-shapes", action="store_true",
                        help="Include non-snake shapes (U-bends etc.)")
    parser.add_argument("--linear", action="store_true",
                        help="Only linear (straight-line) shapes: L cells in a row")
    parser.add_argument("--table", type=int, metavar="MAX_L", default=0,
                        help="Print count table for L=1..MAX_L (Shapes + Total)")
    parser.add_argument("--parallel", action="store_true",
                        help="Use parallel processing per shape (for --table)")
    parser.add_argument("--draw", action="store_true",
                        help="Draw all unique snakes (enumerates representatives)")
    parser.add_argument("--save", type=str, metavar="PATH", default=None,
                        help="Save drawing to file (use with --draw)")
    parser.add_argument("--cols", type=int, default=5,
                        help="Columns in snake drawing grid (default 5)")
    parser.add_argument("--limit", type=int, metavar="N", default=None,
                        help="Draw only first N snakes (use when total is large, e.g. L=6)")
    parser.add_argument("--export", type=str, metavar="PATH", default=None,
                        help="Export examples to JSON file (use with --limit for count)")
    args = parser.parse_args()

    if not args.table and args.L is None:
        parser.error("L is required unless --table is used")
    if args.table:
        import time
        max_l = max(1, args.table)
        label = "linear" if args.linear else "snake polyomino"
        print(f"{label} colored snakes (n=4), L=1..{max_l}")
        print(f"{'L':>4}  {'Shapes':>8}  {'Total':>16}")
        print("-" * 32)
        for L in range(1, max_l + 1):
            try:
                total, shapes, counts = count_all_colored_snakes(
                    L, n=4, snake_only=True, linear_only=args.linear, parallel=args.parallel
                )
                print(f"{L:>4}  {len(shapes):>8}  {total:>16,}")
            except Exception as e:
                print(f"{L:>4}  ERROR: {e}")
            sys.stdout.flush()
        return 0

    if args.n != 4:
        print("Error: only n=4 (squares) supported.")
        return 1

    snake_only = not args.all_shapes

    if args.draw or args.export:
        max_n = args.limit if args.limit is not None else (100 if args.export else None)
        all_snakes, shapes, counts = enumerate_colored_snakes(
            args.L, args.n, snake_only=snake_only, linear_only=args.linear, count_only=False,
            max_snakes=max_n,
        )
        total = len(all_snakes)
        print(f"L={args.L}: {total} unique snakes" + (f" (exporting first {args.limit})" if args.export and args.limit else ""))
        if args.export:
            import json
            to_export = all_snakes[: (args.limit or 100)]
            data = [{"index": i + 1, "coloring": [list(c) for c in walk_coloring[1]]}
                    for i, walk_coloring in enumerate(to_export)]
            with open(args.export, "w") as f:
                json.dump({"L": args.L, "linear": args.linear, "count": len(data), "examples": data}, f, indent=2)
            print(f"Exported {len(data)} examples to {args.export}")
        if args.draw:
            to_draw = all_snakes[: args.limit] if args.limit is not None else all_snakes
            if args.limit is not None and total > args.limit:
                print(f"Drawing first {args.limit} (use --limit to change)")
            if to_draw:
                draw_all_snakes(
                    to_draw,
                    cols=args.cols,
                    save_path=args.save,
                )
        return 0

    if args.full:
        _, shapes, counts = enumerate_colored_snakes(args.L, args.n,
                                                     snake_only=snake_only, linear_only=args.linear,
                                                     count_only=True)
        label = "linear" if args.linear else ("snake polyomino" if snake_only else "all self-avoiding")
        print(f"\nL={args.L}, n={args.n}  ({label} shapes)")
        print(f"{'Shape':>8}  {'Distinct colorings':>20}")
        print("-" * 32)
        for si, walk in enumerate(shapes):
            print(f"  #{si + 1:3d}    {counts[si]:>14d}")
        print("-" * 32)
        print(f"  Total:  {sum(counts):>14d}")
    else:
        shapes = enumerate_linear_shapes(args.L) if args.linear else enumerate_shapes(args.L, snake_only=snake_only)
        label = "linear" if args.linear else ("snake polyomino" if snake_only else "all")
        print(f"L={args.L}: {len(shapes)} distinct {label} shapes")
    return 0


if __name__ == "__main__":
    sys.exit(main())
