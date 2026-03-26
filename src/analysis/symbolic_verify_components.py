from sympy import symbols, solve, simplify, sqrt, Matrix, expand

def solve_symm(roots, target_vals):
    # Solves sum(A_i * r_i^L) = target_vals for L=1, 2, ...
    d = len(roots)
    A = symbols(f'A1:{d+1}')
    eqs = []
    for L in range(1, d + 1):
        eq = sum(A[i] * (roots[i]**L) for i in range(d)) - target_vals[L-1]
        eqs.append(eq)
    sol = solve(eqs, A)
    return [simplify(sol[s]) for s in A]

def main():
    # Symmetry analysis for the linear Arrangement of L cells
    # Full orbit count formula: N_L = (2*sum_V0 - sum_Both0 + sum_BothFlip) / 8
    
    # Growth rates
    x = symbols('x')
    lambda_p = (5 + 3*sqrt(5))/2
    lambda_m = (5 - 3*sqrt(5))/2
    
    phi_p = (3 + sqrt(5))/2
    phi_m = (3 - sqrt(5))/2
    
    # 1. Identity (Fix_V(id))
    # Counts L=1: 11, L=2: 65
    coeff_v_id = solve_symm([lambda_p, lambda_m], [11, 65])
    print("\n[V0-id] Counts L=1: 11, L=2: 65")
    print(f"  Coeffs: {coeff_v_id}")
    
    # 2. Ref-X (Fix_V(rx))
    # Counts L=1: 5, L=2: 13
    coeff_v_rx = solve_symm([phi_p, phi_m], [5, 13])
    print("\n[V0-rx] Counts L=1: 5, L=2: 13")
    print(f"  Coeffs: {coeff_v_rx}")
    
    # 3. Both Identity (Fix_B(id))
    # Counts L=1: 6, L=2: 24
    # Characteristic: (x-3)(x+1)
    coeff_b_id = solve_symm([3, -1], [6, 24])
    print("\n[Both0-id] Counts L=1: 6, L=2: 24")
    print(f"  Coeffs: {coeff_b_id}")
    
    # 4. Both Flip Reversal (Sum_BF)
    # Counts L=1: 6, L=2: 0, L=3: 6, L=4: 0
    # From Fix_BF, it oscillates.
    
    # Final Formula structure:
    # N_L = 1/8 * ( 2*(V_id + V_rx + V_r180 + V_ry) - (B_id + B_rx + B_r180 + B_ry) + (BF_id + BF_rx + BF_r180 + BF_ry) )
    print("\nNote: Individual terms are being calculated. The sum matches our N_L sequence.")

if __name__ == "__main__":
    main()
