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

class MinimalOrderFinder:
    def __init__(self):
        self.patterns = np.array(VALID_PATTERNS)
        self.patterns_both = np.array([p for p in VALID_PATTERNS if sum(p) == 2])

    def _get_transfer_matrix(self, pats):
        M = np.zeros((2, 2), dtype=object)
        for p in pats: M[p[W], p[E]] += 1
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
                    if p[E] == gp[W]: count += (1 if k == 1 else v_mid[p[W]])
            else:
                fc = [p for p in pats if list(p) == [p[perm[0]], p[perm[1]], p[perm[2]], p[perm[3]] if not flip else 1-p[perm[0]], 1-p[perm[1]], 1-p[perm[2]], 1-p[perm[3]]]]
                # Fix for fc list comp
                center_fixed = []
                for p in pats:
                    gp = [p[perm[0]], p[perm[1]], p[perm[2]], p[perm[3]]]
                    if flip: gp = [1-x for x in gp]
                    if list(p) == gp: center_fixed.append(p)
                if k == 0: return len(center_fixed)
                for pc in center_fixed:
                    for p_prev in pats:
                        if p_prev[E] == pc[W]:
                            gp_prev = [p_prev[perm[0]], p_prev[perm[1]], p_prev[perm[2]], p_prev[perm[3]]]
                            if flip: gp_prev = [1-x for x in gp_prev]
                            if gp_prev[W] == pc[E]:
                                count += (1 if k == 1 else v_mid[p_prev[W]])
            return count

    def solve(self, L):
        if L == 1: G0 = [(di, False) for di in range(8)]
        else: G0 = [(0, False), (4, False), (2, True), (5, True)]
        v0 = sum(self.count_fixed(L, di, rev, False) for di, rev in G0)
        b0 = sum(self.count_fixed(L, di, rev, False, True) for di, rev in G0)
        bf = sum(self.count_fixed(L, di, rev, True, True) for di, rev in G0)
        return int((2*v0 - b0 + bf) // (2*len(G0)))

def find_minimal_order(seq):
    n = len(seq)
    for d in range(1, 20):
        H = []
        for i in range(d): H.append(seq[i:i+d])
        mat = np.array(H, dtype=float)
        if np.linalg.cond(mat) > 1e12: # Check for singularity
            # Minimal order found
            b = Matrix(seq[d:2*d])
            M = Matrix(H)
            # Find max d' such that M[:d', :d'] is invertible
            # Or just use the last stable d
            pass
    # Actual iterative search
    for d in range(1, n // 2):
        M = Matrix([seq[i:i+d] for i in range(d)])
        if M.det() != 0:
            # Check if it predicts next term
            coeffs = M.LUsolve(Matrix(seq[d:2*d]))
            if sum(coeffs[j]*seq[d+j] for j in range(d)) == seq[2*d]:
                # potentially minimal
                valid = True
                for i in range(n - d):
                    if sum(coeffs[j]*seq[i+j] for j in range(d)) != seq[i+d]:
                        valid = False; break
                if valid: return d, coeffs
    return None, None

finder = MinimalOrderFinder()
seq = [finder.solve(L) for L in range(2, 45)]
d, c = find_minimal_order(seq)
print(f"Minimal Order (L>=2): {d}")
if d:
    print(f"Coefficients: {[int(ci) for ci in c]}")
    x = symbols('x')
    poly = x**d - sum(c[i]*x**i for i in range(d))
    print(f"Minimal Polynomial: {expand(poly)}")
    
print("\n--- Structural Justifications ---")
print("11 Cell Patterns (sum <= 2):")
for p in VALID_PATTERNS: print(f"  {p}")

print("\nGroup Size Logic:")
print(f"L=1 geometric stabilizer size: {8} (Full square symmetry)")
print(f"L=2 horizontal path geometric stabilizer size: {4} (id, ref_x, rot180, ref_y)")

print("\nTransfer Matrix M (West -> East):")
M = finder._get_transfer_matrix(finder.patterns)
print(f"  M = {M.tolist()}")
print(f"  Determinant(M) = {M[0,0]*M[1,1] - M[0,1]*M[1,0]}")
print(f"  Characteristic: (4-x)(1-x) - 9 = x^2 - 5x - 5")
