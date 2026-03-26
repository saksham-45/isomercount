import numpy as np
from sympy import symbols, expand, simplify, Matrix, solve

# Constants
S, E, N, W = 0, 1, 2, 3
OPP = {S: N, E: W, N: S, W: E}
EDGE_PERMS = (
    (0, 1, 2, 3), (3, 0, 1, 2), (2, 3, 0, 1), (1, 2, 3, 0),
    (2, 1, 0, 3), (0, 3, 2, 1), (3, 2, 1, 0), (1, 0, 3, 2)
)
VALID_PATTERNS = [(a,b,c,d) for a in (0,1) for b in (0,1) for c in (0,1) for d in (0,1) if a+b+c+d <= 2]
VALID_BOTH = [p for p in VALID_PATTERNS if sum(p) == 2]

class RigorousSolver:
    def __init__(self):
        self.patterns = np.array(VALID_PATTERNS)
        self.patterns_both = np.array(VALID_BOTH)

    def _get_transfer_matrix(self, pats):
        M = np.zeros((2, 2), dtype=object)
        for p in pats:
            M[p[W], p[E]] += 1
        return M

    def count_fixed(self, L, di, rev, flip, both_only=False):
        perm = EDGE_PERMS[di]
        pats = self.patterns_both if both_only else self.patterns
        if not rev:
            fixed_pats = []
            for p in pats:
                tp = [p[perm[0]], p[perm[1]], p[perm[2]], p[perm[3]]]
                if flip: tp = [1-x for x in tp]
                if list(p) == tp: fixed_pats.append(p)
            if not fixed_pats: return 0
            M = np.zeros((2, 2), dtype=object)
            for p in fixed_pats: M[p[W], p[E]] += 1
            if L == 1: return len(fixed_pats)
            v_init = np.zeros(2, dtype=object)
            for p in fixed_pats: v_init[p[E]] += 1
            v_end = np.zeros(2, dtype=object)
            for p in fixed_pats: v_end[p[W]] += 1
            return (v_init @ np.linalg.matrix_power(M, L - 2) @ v_end)
        else:
            k, rem = divmod(L, 2)
            M_full = self._get_transfer_matrix(pats)
            v_init = np.zeros(2, dtype=object)
            for p in pats: v_init[p[E]] += 1
            v_mid = (v_init @ np.linalg.matrix_power(M_full, k - 2)) if k > 1 else (v_init if k == 1 else None)
            count = 0
            if rem == 0:
                for p in pats:
                    gp = [p[perm[0]], p[perm[1]], p[perm[2]], p[perm[3]]]
                    if flip: gp = [1-x for x in gp]
                    if p[E] == gp[W]:
                        count += (1 if k == 1 else v_mid[p[W]])
            else:
                fixed_center = []
                for p in pats:
                    gp = [p[perm[0]], p[perm[1]], p[perm[2]], p[perm[3]]]
                    if flip: gp = [1-x for x in gp]
                    if list(p) == gp: fixed_center.append(p)
                if k == 0: return len(fixed_center)
                for pc in fixed_center:
                    for p_prev in pats:
                        if p_prev[E] == pc[W]:
                            gp_prev = [p_prev[perm[0]], p_prev[perm[1]], p_prev[perm[2]], p_prev[perm[3]]]
                            if flip: gp_prev = [1-x for x in gp_prev]
                            if gp_prev[W] == pc[E]:
                                count += (1 if k == 1 else v_mid[p_prev[W]])
            return count

    def solve_components(self, L):
        if L == 1: G0 = [(di, False) for di in range(8)]
        else: G0 = [(0, False), (4, False), (2, True), (5, True)]
        v0 = [self.count_fixed(L, di, rev, False) for di, rev in G0]
        b0 = [self.count_fixed(L, di, rev, False, True) for di, rev in G0]
        bf = [self.count_fixed(L, di, rev, True, True) for di, rev in G0]
        g_size = len(G0)
        nl = (2*sum(v0) - sum(b0) + sum(bf)) // (2 * g_size)
        return int(nl), v0, b0, bf

def berlekamp_massey(seq):
    # Minimal order using rank of Hankel matrix
    n = len(seq)
    for d in range(1, n // 2):
        H = []
        for i in range(d + 1):
            H.append(seq[i:i+d+1])
        if np.linalg.matrix_rank(np.array(H, dtype=float)) == d:
            # check if d is stable
            continue
    # Simple check: solve for coefficients and verify on remainder
    # Use d = 15 as target
    d = 15
    M = []
    for i in range(d): M.append(seq[i:i+d])
    M = Matrix(M)
    b = Matrix(seq[d:2*d])
    coeffs = M.LUsolve(b)
    for i in range(n - d):
        expected = sum(coeffs[j]*seq[i+j] for j in range(d))
        if seq[i+d] != expected: return None, None
    return d, coeffs

solver = RigorousSolver()
# Compute many terms to be sure
terms = [solver.solve_components(L)[0] for L in range(2, 60)]
order, coeffs = berlekamp_massey(terms)
print(f"Minimal Order (L>=2): {order}")
if order:
    x = symbols('x')
    poly = x**order - sum(coeffs[i]*x**i for i in range(order))
    print(f"Minimal Polynomial: {expand(poly)}")

# Check L=1 discontinuity
nl1, v1, b1, bf1 = solver.solve_components(1)
print(f"L=1 orbit count: {nl1} (Group size 8*2=16)")
nl2, v2, b2, bf2 = solver.solve_components(2)
print(f"L=2 orbit count: {nl2} (Group size 4*2=8)")

# Show term stabilization for (x-1) factor
# For L >= 2, the contribution of constant symmetries stabilizes.
comps = [solver.solve_components(L) for L in range(2, 10)]
# Look at sum(bf) - sum(b0)
offsets = [ (2*sum(c[1]) - sum(c[2]) + sum(c[3])) % 8 for c in comps]
print(f"Burnside Remainder mod 8: {offsets}")
