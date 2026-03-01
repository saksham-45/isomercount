import numpy as np

def find_recurrence(seq, max_order=15):
    n = len(seq)
    for d in range(1, max_order + 1):
        if 2 * d > n: break
        # Solve c_0*s_i + c_1*s_{i+1} + ... + c_{d-1}*s_{i+d-1} = s_{i+d}
        A = []
        b = []
        for i in range(n - d):
            A.append(seq[i:i+d])
            b.append(seq[i+d])
        A = np.array(A, dtype=object)
        b = np.array(b, dtype=object)
        
        # Try to solve exactly
        try:
            # Using float for faster check, then can refine
            Af = A.astype(float)
            bf = b.astype(float)
            if np.linalg.matrix_rank(Af) < d:
                continue
            coeffs = np.linalg.solve(Af.T @ Af, Af.T @ bf)
            # Check if it holds for the whole sequence
            predicted = []
            valid = True
            for i in range(n - d):
                val = sum(coeffs[j] * seq[i+j] for j in range(d))
                if abs(val - seq[i+d]) > 1e-2:
                    valid = False
                    break
            if valid:
                return d, coeffs
        except:
            continue
    return None, None

# Sequence from the L=1..20 run
nl_seq = [
    4, 21, 109, 586, 3326, 19209, 111871, 653758, 3824678, 22387074, 
    131052313, 767211817, 4491420695, 26293679325, 153927402355, 
    901112477374, 5275222867382, 30881755696845, 180785144397637, 1058335309256578
]

# Note: L=1 might be an outlier due to symmetry group size change (8 vs 4)
# Let's try L=2..20
order, coeffs = find_recurrence(nl_seq[1:], max_order=12)
if order:
    print(f"Recurrence found (L=2..20) of order {order}:")
    print("Coeffs:", coeffs)
    # Convert to integer if close
    int_coeffs = [round(c) for c in coeffs]
    print("Integer Coeffs:", int_coeffs)
    poly_str = "x^" + str(order)
    for i, c in enumerate(int_coeffs):
        sign = "+" if c >= 0 else "-"
        poly_str += f" {sign} {abs(c)}x^{i}"
    print("Polynomial:", poly_str)
else:
    print("No recurrence found up to order 12.")
