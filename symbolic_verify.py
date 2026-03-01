from sympy import symbols, solve, Poly, expand, simplify, sqrt, I, Matrix

def get_exact_symbolic():
    x = symbols('x')
    
    # Define factors
    p1 = x**2 - 5*x - 5
    p2 = x**2 - 3*x + 1
    p3 = x**4 - 5*x**2 - 5
    p4 = x**2 - 2*x - 3
    p5 = x**4 - 2*x**2 - 3
    
    P_full = expand(p1 * p2 * p3 * p4 * p5)
    print("Full Expanded Polynomial P(x):")
    print(P_full)
    
    # Solve for roots exactly
    roots_raw = solve(P_full, x)
    # Ensure they are unique (should be 14)
    roots = list(set(roots_raw))
    print(f"\nNumber of unique roots: {len(roots)}")
    
    # Initial conditions
    a_vals = [4, 21, 109, 586, 3326, 19209, 111871, 653758, 3824678, 22387074, 131052313, 767211817, 4491420695, 26293679325]
    
    # Vandermonde-like system: sum(A_i * r_i^L) = a(L) for L=1..14
    A_syms = symbols(f'A1:{len(roots)+1}')
    
    equations = []
    for L in range(1, 15):
        eq = sum(A_syms[i] * (roots[i]**L) for i in range(len(roots))) - a_vals[L-1]
        equations.append(eq)
        
    print("\nSolving for symbolic coefficients A_i...")
    sol = solve(equations, A_syms)
    
    print("\nExact Symbolic Coefficients A_i:")
    for i, s in enumerate(A_syms):
        val = simplify(sol[s])
        print(f"A_{i+1} (root {roots[i]}): {val}")

if __name__ == "__main__":
    get_exact_symbolic()
