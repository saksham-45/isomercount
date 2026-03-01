import numpy as np

# Directions: S=0, E=1, N=2, W=3
S, E, N, W = 0, 1, 2, 3
OPP = {S: N, E: W, N: S, W: E}

VALID_PATTERNS = [
    (a, b, c, d) for a in (0, 1) for b in (0, 1) for c in (0, 1) for d in (0, 1) 
    if a + b + c + d <= 2
]

EDGE_PERMS = (
    (0, 1, 2, 3), (3, 0, 1, 2), (2, 3, 0, 1), (1, 2, 3, 0),
    (2, 1, 0, 3), (0, 3, 2, 1), (3, 2, 1, 0), (1, 0, 3, 2)
)

VALID_BOTH = [p for p in VALID_PATTERNS if sum(p) == 2]

class FastSolver:
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
                if flip: tp = [1 - x for x in tp]
                if list(p) == tp: fixed_pats.append(p)
            if not fixed_pats: return 0
            
            M = np.zeros((2, 2), dtype=object)
            for p in fixed_pats:
                M[p[W], p[E]] += 1
            
            if L == 1: return len(fixed_pats)
            v_init = np.zeros(2, dtype=object)
            for p in fixed_pats: v_init[p[E]] += 1
            v_end = np.zeros(2, dtype=object)
            for p in fixed_pats: v_end[p[W]] += 1
            ML = np.linalg.matrix_power(M, L - 2)
            res = v_init @ ML
            return np.dot(res, v_end)
        else:
            k = L // 2
            M_full = self._get_transfer_matrix(pats)
            if L % 2 == 0:
                v_init = np.zeros(2, dtype=object)
                for p in pats: v_init[p[E]] += 1
                if k > 1:
                    ML = np.linalg.matrix_power(M_full, k - 2)
                    v_mid = v_init @ ML
                else: v_mid = None
                count = 0
                for p in pats:
                    gp = [p[perm[0]], p[perm[1]], p[perm[2]], p[perm[3]]]
                    if flip: gp = [1 - x for x in gp]
                    if p[E] == gp[W]:
                        if k == 1: count += 1
                        else: count += v_mid[p[W]]
                return count
            else:
                fixed_center = []
                for p in pats:
                    gp = [p[perm[0]], p[perm[1]], p[perm[2]], p[perm[3]]]
                    if flip: gp = [1 - x for x in gp]
                    if list(p) == gp: fixed_center.append(p)
                if not fixed_center: return 0
                if k == 0: return len(fixed_center)
                v_init = np.zeros(2, dtype=object)
                for p in pats: v_init[p[E]] += 1
                if k > 1:
                    ML = np.linalg.matrix_power(M_full, k - 2)
                    v_mid = v_init @ ML
                else: v_mid = None
                count = 0
                for pc in fixed_center:
                    for p_prev in pats:
                        if p_prev[E] == pc[W]:
                            gp_prev = [p_prev[perm[0]], p_prev[perm[1]], p_prev[perm[2]], p_prev[perm[3]]]
                            if flip: gp_prev = [1 - x for x in gp_prev]
                            if gp_prev[W] == pc[E]:
                                if k == 1: count += 1
                                else: count += v_mid[p_prev[W]]
                return count

    def solve(self, L):
        if L == 1: G0 = [(di, False) for di in range(8)]
        else: G0 = [(0, False), (4, False), (2, True), (5, True)]
        sum_V0 = sum(self.count_fixed(L, di, rev, False) for di, rev in G0)
        sum_Both0 = sum(self.count_fixed(L, di, rev, False, True) for di, rev in G0)
        sum_BothFlip = sum(self.count_fixed(L, di, rev, True, True) for di, rev in G0)
        total_num = 2 * sum_V0 - sum_Both0 + sum_BothFlip
        return int(total_num // (2 * len(G0)))

def solve_linear_recurrence(seq):
    from sympy import Matrix, symbols, Poly, solve
    n = len(seq)
    for d in range(1, n // 2):
        # We need d+1 equations to find d coeffs exactly (or d equations)
        # Solve c_0 s_0 + c_1 s_1 + ... + c_{d-1} s_{d-1} = s_d
        # ...
        A = []
        b = []
        for i in range(d):
            A.append(seq[i:i+d])
            b.append(seq[i+d])
        A = Matrix(A)
        b = Matrix(b)
        if A.det() == 0: continue
        coeffs = A.LUsolve(b)
        # Verify
        valid = True
        for i in range(d, n - d):
            actual = seq[i+d]
            expected = sum(coeffs[j] * seq[i+j] for j in range(d))
            if actual != expected:
                valid = False
                break
        if valid:
            return d, coeffs
    return None, None

solver = FastSolver()
seq = [solver.solve(L) for L in range(2, 45)]
print(f"Sequence L=2..44 computed.")

d, coeffs = solve_linear_recurrence(seq)
if d:
    print(f"Order: {d}")
    print("Coefficients [c_0, c_1, ..., c_{d-1}]:")
    print([int(c) for c in coeffs])
    x = symbols('x')
    poly = x**d - sum(coeffs[i] * x**i for i in range(d))
    print("Characteristic Polynomial P(x):")
    print(poly)
else:
    print("No recurrence found.")
