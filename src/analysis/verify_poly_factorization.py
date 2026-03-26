import sympy
from sympy import symbols, expand, simplify, solve, Matrix, sqrt

def verify():
    x = symbols('x')
    
    # User's proposed factors
    p1 = x**2 - 5*x - 5
    p2 = x**2 - 3*x + 1
    p3 = x**4 - 5*x**2 - 5
    p4 = x**2 - 2*x - 3
    p5 = x**4 - 2*x**2 - 3
    
    P_proposed = expand(p1 * p2 * p3 * p4 * p5)
    print("Expanded P(x):")
    print(P_proposed)
    
    # Recurrence coefficients (highest degree first)
    coeffs = [1, -10, 17, 82, -224, -124, 494, -86, 394, 110, -1060, -320, -495, -300, 225]
    P_target = sum(c * x**(14-i) for i, c in enumerate(coeffs))
    
    print("\nDegree of target:", P_target.as_poly().degree())
    print("Does expanded == target?", simplify(P_proposed - P_target) == 0)

    # Let's check the roots of p1
    r1 = solve(p1, x)
    print("\nRoots of p1 (Identity):", r1)
    # Roots are (5 +/- 3*sqrt(5))/2. No, solve(x^2 - 5x - 5) is (5 +/- sqrt(45))/2 = (5 +/- 3*sqrt(5))/2.
    # Dominant root is (5 + 3*sqrt(5))/2 approx 5.8541.
    
    # Roots of p2 (Vertical Flip rx)
    r2 = solve(p2, x)
    print("Roots of p2 (Vertical Flip rx):", r2)
    # (3 +/- sqrt(5))/2.
    
    # Roots of p3 (Reversal)
    r3 = solve(p3, x)
    print("Roots of p3 (Reversal):", r3)
    # These are ±sqrt((5 ± 3*sqrt(5))/2).

    # Roots of p4 (Both Identity)
    r4 = solve(p4, x)
    print("Roots of p4 (Both Identity):", r4)
    # Roots are 3 and -1.
    
    # Roots of p5 (Both Reversal)
    r5 = solve(p5, x)
    print("Roots of p5 (Both Reversal):", r5)
    # Roots are ±sqrt(3) and ±i.

if __name__ == "__main__":
    verify()
