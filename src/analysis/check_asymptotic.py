import numpy as np
from sympy import sqrt, simplify

def solve_exact(L):
    # Simplified Burnside components for Identity
    # Transfer matrix matches West to East
    # M = [[4, 3], [3, 1]]
    # Dominant eigenvalue lambda = (5 + 3*sqrt(5))/2
    # For L cells, count is v_init @ M^(L-2) @ v_end
    # v_init = [7, 4] (Total patterns starting with W=0, W=1) -- wait
    # Actually, v_init[v_E] is count of patterns with p[E]=v_E
    # pats = [(0,0,0,0)...] total 11.
    # pats with E=0: (0,0,0,0), (0,0,0,1), (0,0,1,0), (1,0,0,0), (1,0,1,0), (0,1,0,0), (1,1,0,0) -> 7
    # pats with E=1: (0,0,0,1), (0,1,0,1), (1,0,0,1), (0,0,1,1) ... let's be precise.
    
    # Precise p[E] counts:
    # Patterns p = (S,E,N,W), sum <= 2
    # E=0: (0,0,0,0), (0,0,1,0), (1,0,0,0), (1,0,1,0), (0,0,0,1), (1,0,0,1), (0,0,1,1) -- wait
    # Let's just list them:
    pats = []
    for s in (0,1):
        for e in (0,1):
            for n in (0,1):
                for w in (0,1):
                    if s+e+n+w <= 2: pats.append((s,e,n,w))
    
    e0 = len([p for p in pats if p[1] == 0]) # E is index 1
    e1 = len([p for p in pats if p[1] == 1])
    # e0 = 7, e1 = 4. Total 11.
    
    w0 = len([p for p in pats if p[3] == 0]) # W is index 3
    w1 = len([p for p in pats if p[3] == 1])
    # w0 = 7, w1 = 4.
    
    M = np.array([[4, 3], [3, 1]], dtype=object)
    v_init = np.array([e0, e1], dtype=object)
    v_end = np.array([w0, w1], dtype=object)
    
    if L == 1: return 11
    res = v_init @ np.linalg.matrix_power(M, L - 2) @ v_end
    return int(res)

print("Checking Asymptotic Coefficient for Raw Count:")
L = 40
raw = solve_exact(L)
lam = (5 + 3*np.sqrt(5))/2
c_raw = (5 + 2*np.sqrt(5))/5
pred_raw = c_raw * (lam**L)
print(f"L={L}")
print(f"Actual Raw Code: {raw}")
print(f"Predicted Raw:   {pred_raw:.5e}")
print(f"Ratio: {raw/pred_raw}")

print("\nOrbit Count Asymptotic:")
# For L=40, G size is 8.
# Orbit count approx raw / 8?
# N_L ~ (raw + FixedPoints) / 8. FixedPoints grow slower.
# So N_L ~ raw / 8.
# If c_raw = (5+2sqrt(5))/5, then c_orbit = c_raw / 8 = (5+2sqrt(5))/40.
# The user claims it's (5+2sqrt(5))/20. This implies |G|=4?

# Let's check actual N_L for L=40
from sympy import Matrix
def solve_nl(L):
    # Use the logic from our previous solver
    # Simplified G0 size is 4 for L>1. Total group size is 2*4=8.
    # N_L = (2*sumV0 - sumB0 + sumBF) / (2 * 4) = (2*sumV0 - sumB0 + sumBF) / 8
    # sumV0 is mostly id_term + rx_term + r180_term + ry_term
    # id_term is the raw count.
    pass

# We already have results for L=20.
actual_l20 = 1058335309256578
pred_l20_40 = ((5 + 2*np.sqrt(5))/40) * (lam**20)
pred_l20_20 = ((5 + 2*np.sqrt(5))/20) * (lam**20)
print(f"L=20 actual: {actual_l20}")
print(f"Predicted (c/40): {pred_l20_40:.10e}")
print(f"Predicted (c/20): {pred_l20_20:.10e}")
print(f"Ratio (actual / c/40): {actual_l20 / pred_l20_40}")
print(f"Ratio (actual / c/20): {actual_l20 / pred_l20_20}")
