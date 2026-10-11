"""Enumerate Weil polynomials passing the Elsenhans--Jahnel test (based on Artin-Tate).

To run from the root of the repo: ``sage scripts/enumerate_weil_polys.sage``
"""

from pathlib import Path
from sage.rings.polynomial.weil.all import WeilPolynomials

PRE_AT = Path(__file__).resolve().parents[1] / "data/weil_polynomials/preAT"
field_sizes = [2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 17, 19]

def passes_elsenhans_jahnel(pol):
    """Test the residual value at -1 after removing all factors of t - 1."""
    t = pol.parent().gen()
    multiplicity = pol.valuation(t - 1)
    residual = pol // ((t - 1) ** multiplicity)
    return ZZ(residual(-1)).is_square()

for field_size in field_sizes:
    q = ZZ(field_size)
    with (PRE_AT / f"q{q}.txt").open("w") as output:
        for sign in (1, -1):
            for polynomial in WeilPolynomials(21, ZZ(1), lead=[q, q**2 + q + 1], sign=sign):
                if passes_elsenhans_jahnel(polynomial):
                    output.write(f"{list(polynomial)}\n")
