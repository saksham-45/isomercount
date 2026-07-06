"""
Exact linear-snake count.

The previous version of this file used a float approximation

    predict_total(L) = 0.4736... * 5.8541...^L

which loses precision for L > ~15 and overflows for L > ~180. The values it
returned were WRONG once L got large enough for the exponent to matter.

The exact integer count for any L is available in
    src.analysis.closed_form.linear_snake_count(L)
which uses the order-12 integer recurrence derived from the transfer matrix.

For a fast asymptotic estimate (as a float / mpf), use
    src.analysis.closed_form.linear_snake_asymptotic(L, digits=50)
with the exact algebraic constants
    A = (5 + 2*sqrt(5))/20 ,   mu = (5 + 3*sqrt(5))/2 .
"""

from src.analysis.closed_form import linear_snake_count, linear_snake_asymptotic


def predict_total(L):
    """Exact integer count of distinct linear snakes with L cells."""
    return linear_snake_count(L)


def predict_total_asymptotic(L, digits=50):
    """Asymptotic estimate  A * mu^L  as a high-precision decimal string."""
    return linear_snake_asymptotic(L, digits=digits)
