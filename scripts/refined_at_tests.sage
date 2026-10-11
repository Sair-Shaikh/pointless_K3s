"""Keep the Weil polynomials that survive the refined Artin--Tate tests.

From ``pointless_K3s``, run ``sage scripts/refined_at_tests.sage``.
Pass field sizes, such as ``16 17``, to run only those files.
"""

import sys
from functools import lru_cache
from math import gcd
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1] / "data/weil_polynomials"
R = PolynomialRing(ZZ, "t")
t = R.gen()

# Lower bounds for rootless lattices of ranks 1 through 21.
DELTA = (0, 4, 12, 32, 64, 128, 192, 256, 256, 298, 370, 400,
         410, 406, 370, 318, 256, 194, 134, 86, 50, 26)
CYC = {n: R.cyclotomic_polynomial(n) for n in range(1, 201)
       if euler_phi(n) <= 21}
       
# For each base degree a and extension degree m, test
# H_{a,m}(t) = 1 + t^a + ... + t^((m-1)a), with am <= 72
# and no prime divisor of m above 19.
ARTIN_TATE_TESTS = [
    (base_degree, extension_degree,
     sum(t ** (base_degree * power) for power in range(extension_degree)))
    for extension_degree in range(2, 73)
    if all(prime <= 19 for prime, _ in ZZ(extension_degree).factor())
    for base_degree in range(1, 72 // extension_degree + 1)
]

def split_cyclotomic(f):
    """Separate the cyclotomic factors from the residual polynomial."""
    factors = []
    for order, cyclotomic in CYC.items():
        multiplicity = f.valuation(cyclotomic)
        if multiplicity:
            factors.append((order, multiplicity))
            f //= cyclotomic ** multiplicity
    return tuple(factors), f


@lru_cache(None)
def relevant_tests(factors):
    """Skip tests whose algebraic part already meets the determinant bound."""
    relevant = []
    for base_degree, extension_degree, test_polynomial in ARTIN_TATE_TESTS:
        invariant_rank, rank, u, norm = 1, 0, ZZ(1), ZZ(1)
        for n, exponent in factors:
            degree = CYC[n].degree() * exponent
            if base_degree % n == 0:
                invariant_rank += degree
            elif (base_degree * extension_degree) % n == 0:
                d = n // gcd(n, base_degree)
                rank += degree
                u *= CYC[d](1) ** (degree // CYC[d].degree())
            else:
                norm *= abs(CYC[n].resultant(test_polynomial)) ** exponent
        if rank:
            j0 = gcd(u, extension_degree ** invariant_rank)
            glue = j0 // (u // j0).squarefree_part()
            if glue * norm < DELTA[rank]:
                relevant.append((test_polynomial, rank, glue * norm))
    return relevant


def survives(f):
    factors, residual = split_cyclotomic(f)
    for h, rank, algebraic_bound in relevant_tests(factors):
        if algebraic_bound * abs(residual.resultant(h)) < DELTA[rank]:
            return False
    return True


def run(path):
    output = ROOT / "postAT" / path.name
    with path.open() as source, output.open("w") as target:
        for line in source:
            coefficients = [ZZ(x) for x in line.strip()[1:-1].split(",")]
            if coefficients[0] < 0:
                coefficients = [-coefficient for coefficient in coefficients]
                line = "[" + ", ".join(map(str, coefficients)) + "]\n"
            if survives(R(coefficients)):
                target.write(line)
    print(f"Wrote {output}")


if __name__ == "__main__":
    paths = ([ROOT / "preAT" / f"q{q}.txt" for q in sys.argv[1:]]
             or sorted((ROOT / "preAT").glob("q*.txt")))
    for path in paths:
        run(path)
