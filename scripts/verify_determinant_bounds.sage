"""Exact verifier for determinant bounds on even rootless lattices using Poisson summation.

Expects as input a CSV with columns: rank, scale_factor, determinant_bound, coefficients

To run from the root of the repo: ``sage scripts/verify_determinant_bounds.sage data/determinant_bounds/certificates.csv''
"""

import csv
import sys
from pathlib import Path
from fractions import Fraction

def construct_polynomials(n, coefficients):
    """
    Construct `P` and `Q` using the normalized Laguerre basis.
    """
    R = PolynomialRing(QQ, 'u')
    u = R.gen()
    alpha = QQ(n)/2 - 1
    P = R(0)
    Q = R(0)
    for j, c in enumerate(coefficients):
        L = R(gen_laguerre(j, alpha, u)) # Cast to R needed for j=0
        L /= L(0)  # Normalize L(0) = 1

        Q += c * L
        P += (-1)^j * c * L
    return P, Q

def verify_positivity(Q):
    """
    Verify `Q(u) > 0` for all `u >= 0`.
    
    It suffices to check that `Q(0) > 0` and that `Q` has no nonnegative real roots.
    Recall that `AA` is predefined as the field of algebraic real numbers.
    """
    return Q(0) > 0 and not any(root >= 0 for root, _ in Q.roots(ring=AA))

def verify_negativity(P, scale_factor):
    """
    Verify `P(4*pi*m/scale_factor) < 0` for every integer `m >= 2`.
    
    It suffices to check that `P` has negative leading coefficient and that the bound holds
    for those values not exceeding the largest positive root of `P`.
    """
    roots = [r for r, _ in P.roots(ring=AA) if r > 0]
    if P(0) <= 0 or not roots or P.leading_coefficient() > 0:
        return False
    largest_root = max(roots)

    # Rigorous rational bounds for pi.
    pi_interval = RealIntervalField(128).pi()
    pi_lo = QQ(pi_interval.lower())
    pi_hi = QQ(pi_interval.upper())

    m = 2
    while True:
        lo = 4*m*pi_lo/scale_factor
        hi = 4*m*pi_hi/scale_factor
        if lo > largest_root:
            return True
        if any(lo <= r <= hi for r in roots):
            raise ValueError(f'Shell m={m} is too close to a root; increase pi precision')
        if P(lo) >= 0:
            return False
        m += 1

def verify_bound(rank, scale_factor, coefficients, determinant_bound):
    """
    Verify the claimed determinant bound using Poisson summation.
    """

    if rank < 1 or scale_factor <= 0 or sum(coefficients) != 1:
        raise ValueError('Require rank >= 1, scale_factor > 0, sum(coefficients) = 1')

    P, Q = construct_polynomials(rank, coefficients)
    assert verify_positivity(Q)
    assert verify_negativity(P, scale_factor)
        
    # Poisson summation gives det(L) >= (f_hat(0)/f(0))^2 = scale_factor^rank / P(0)^2.
    analytic_bound = scale_factor^rank / P(0)^2

    # Round up to an integer bound.
    minimum = analytic_bound.ceil()

    # In odd rank, the determinant must be even.
    if rank % 2 and minimum % 2:
        minimum += 1

    RIF = RealIntervalField(128)
    t = 2 * RIF.pi() / scale_factor
    if P(2*t) <= 0:
        raise ValueError('Require positive value at shell 1')
    root20bound = ZZ((analytic_bound / (1 + 20 * P(2*t) * exp(-t) / P(0))^2).ceil())
    if rank % 2 and root20bound % 2:
        root20bound += 1
    return minimum >= determinant_bound, root20bound

def main(filename):
    with open(filename, newline='') as f:
        for row in csv.DictReader(f):
            rank = ZZ(row['rank'])
            scale_factor = QQ(row['scale_factor'])
            coefficients = [QQ(Fraction(c.strip())) for c in row['coefficients'].strip('[]').split(',')]
            bound = ZZ(row['determinant_bound'])
            verify, root20bound = verify_bound(rank, scale_factor, coefficients, bound)
            if verify:
                print(f"Verified rank {rank} bound = {bound}; for <= 20 roots, bound = {root20bound}")               
            else:
                print(f"Failed rank {rank} bound = {bound}")

if __name__ == '__main__':
    default = Path(__file__).resolve().parents[1] / 'data/determinant_bounds/certificates.csv'
    main(sys.argv[1] if len(sys.argv) > 1 else default)
