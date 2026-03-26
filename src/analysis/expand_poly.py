from sympy import Matrix, symbols, Poly

def solve_linear_recurrence(seq):
    n = len(seq)
    for d in range(1, n // 2):
        M = []
        for i in range(d):
            M.append(seq[i:i+d])
        M = Matrix(M)
        if M.det() == 0: continue
        b = Matrix(seq[d:2*d])
        coeffs = M.LUsolve(b)
        
        # Verify on the rest
        valid = True
        for i in range(0, n - d):
            expected = sum(coeffs[j] * seq[i+j] for j in range(d))
            if seq[i+d] != expected:
                valid = False
                break
        if valid:
            return d, coeffs
    return None, None

# Recalculate or use sequence
# Using the data from the previous run
seq = [21, 109, 586, 3326, 19209, 111871, 653758, 3824678, 22387074, 
       131052313, 767211817, 4491420695, 26293679325, 153927402355, 
       901112477374, 5275222867382, 30881755696845, 180785144397637, 1058335309256578,
       6195548325603701, 36269389201990421, 212324203673727222, 1242967677843336330,
       7276332194600109593, 42597284411135431357, 249369904944837563502,
       1459846332029471191542, 8546124719082264639941, 50029525281488667503029,
       292878417937749440361274] # Seq starting from L=2

d, coeffs = solve_linear_recurrence(seq)
if d:
    x = symbols('x')
    poly = x**d - sum(coeffs[i] * x**i for i in range(d))
    print(f"Order: {d}")
    print("Coefficients [c_0, c_1, ..., c_{11}]:")
    print([int(c) for c in coeffs])
    print("Poly expanded:")
    print(poly)
