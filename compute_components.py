import numpy as np

# Directions: S=0, E=1, N=2, W=3
S, E, N, W = 0, 1, 2, 3
OPP = {S: N, E: W, N: S, W: E}

# 11 patterns with sum <= 2
VALID_PATTERNS = [
    (a, b, c, d) for a in (0, 1) for b in (0, 1) for c in (0, 1) for d in (0, 1) 
    if a + b + c + d <= 2
]

# 8 Dihedral transformations (edge permutations)
# 0: id, 1: r90, 2: r180, 3: r270, 4: ref_x, 5: ref_y, 6: ref_diag1, 7: ref_diag2
EDGE_PERMS = (
    (0, 1, 2, 3), # 0 identity
    (3, 0, 1, 2), # 1 rot 90 CCW
    (2, 3, 0, 1), # 2 rot 180
    (1, 2, 3, 0), # 3 rot 270 CCW
    (2, 1, 0, 3), # 4 reflect x (y -> -y)
    (0, 3, 2, 1), # 5 reflect y (x -> -x)
    (3, 2, 1, 0), # 6 reflect y=x
    (1, 0, 3, 2), # 7 reflect y=-x
)

VALID_BOTH = [p for p in VALID_PATTERNS if sum(p) == 2]

class ComponentSolver:
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
                else:
                    v_mid = None
                
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
                else:
                    v_mid = None
                
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

def main():
    solver = ComponentSolver()
    # G0 for L>1 is {0, 4, 2_rev, 5_rev}
    # For L=1 it's all 8.
    
    print("L | id | rx | r180 | ry | sum_V0 | sum_Both0 | sum_BothFlip | N_L")
    print("-" * 70)
    for L in range(1, 11):
        if L == 1:
            G0 = [(di, False) for di in range(8)]
            v_scores = [solver.count_fixed(L, di, rev, False) for di, rev in G0]
            sv0 = sum(v_scores)
            sb0 = sum(solver.count_fixed(L, di, rev, False, True) for di, rev in G0)
            sbf = sum(solver.count_fixed(L, di, rev, True, True) for di, rev in G0)
            # For table simplicity, just show the first 4 for L=1
            id_val, rx_val, r180_val, ry_val = v_scores[0], v_scores[4], v_scores[2], v_scores[5]
            nl = (2*sv0 - sb0 + sbf) // (2*8)
        else:
            G0 = [(0, False), (4, False), (2, True), (5, True)]
            id_val = solver.count_fixed(L, 0, False, False)
            rx_val = solver.count_fixed(L, 4, False, False)
            r180_val = solver.count_fixed(L, 2, True, False)
            ry_val = solver.count_fixed(L, 5, True, False)
            sv0 = id_val + rx_val + r180_val + ry_val
            sb0 = sum(solver.count_fixed(L, di, rev, False, True) for di, rev in G0)
            sbf = sum(solver.count_fixed(L, di, rev, True, True) for di, rev in G0)
            nl = (2*sv0 - sb0 + sbf) // (2*4)
        
        print(f"{L} | {id_val} | {rx_val} | {r180_val} | {ry_val} | {sv0} | {sb0} | {sbf} | {nl}")

if __name__ == "__main__":
    main()
