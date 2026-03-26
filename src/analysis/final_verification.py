from sympy import symbols, expand

def final_verify():
    x = symbols('x')
    P = (x**2-5*x-5)*(x**2-3*x+1)*(x**4-5*x**2-5)*(x**2-2*x-3)*(x**4-2*x**2-3)
    c = P.as_poly().all_coeffs()
    c = [int(val) for val in c]
    
    # Sequence L=1..20
    # Recalculating precisely for n=1..20 to ensure no errors
    # I already have the data from previous runs.
    seq_l1_to_l20 = [
        4, 21, 109, 586, 3326, 19209, 111871, 653758, 3824678, 22387074, 
        131052313, 767211817, 4491420695, 26293679325, 153927402355, 
        901112477374, 5275222867382, 30881755696845, 180785144397637, 1058335309256578
    ]
    
    print(f"Order: {len(c)-1}")
    print(f"Polynomial Coefficients (c_0 to c_{len(c)-1}):")
    print(c)
    
    # Check recurrence for L=15..20
    # Equation: sum_{j=0}^{14} c_j * N_{L-j} = 0
    # Note: seq is 0-indexed, so seq[L-1] is N_L
    for L in range(15, 21):
        # We need N_L down to N_{L-14}
        # indices L-1 down to L-15
        val = 0
        for j in range(15):
            val += c[j] * seq_l1_to_l20[L-1-j]
        print(f"L={L}: Residue = {val}")

if __name__ == "__main__":
    final_verify()
