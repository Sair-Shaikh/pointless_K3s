"""Enumerate Weil polynomials passing the classical Artin--Tate square test.

From the repository root, run ``sage scripts/enumerate_weil_polys.sage``.
"""

from pathlib import Path
from sage.rings.polynomial.weil.all import WeilPolynomials

PRE_AT = Path(__file__).resolve().parents[1] / "data/weil_polynomials/preAT"
field_sizes = [2, 3, 4, 5, 7, 8, 9, 11, 13, 17, 17, 19]

def passes_artin_tate(polynomial):
    """Test the residual value at -1 after removing all factors t - 1."""
    t = polynomial.parent().gen()
    multiplicity = polynomial.valuation(t - 1)
    residual = polynomial // ((t - 1) ** multiplicity)
    return ZZ(residual(-1)).is_square()


for field_size in field_sizes:
    q = ZZ(field_size)
    with (PRE_AT / f"q{q}.txt").open("w") as output:
        for sign in (1, -1):
            for polynomial in WeilPolynomials(21, ZZ(1), lead=[q, q**2 + q + 1], sign=sign):
                if passes_artin_tate(polynomial):
                    output.write(f"{list(polynomial)}\n")
