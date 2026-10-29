"""
Exact verifier for Farkas certificates for nonexistence of Niemeier embeddings with rootless
orthogonal complement.

Expects as input a CSV with columns: rank, genus, matrix

To run from the root of the repo: ``sage scripts/verify_farkas_certificates.sage data/farkas_certificates/certificates.csv''
"""

import csv
import sys
from pathlib import Path

def matrix_from_upper_triangular(matrix_upper):
    n = len(matrix_upper)
    if not is_square(8*n+1):
        raise ValueError('input has wrong length')
    d = ZZ(sqrt(8*n + 1) - 1) // 2
    m = Matrix(ZZ, d, d)
    c = 0
    for i in range(d):
        for j in range(i, d):
            m[i, j] = matrix_upper[c]
            m[j, i] = matrix_upper[c]
            c += 1
    return(m)

def verify_genus(determinant, genus_str, gram_matrix):
    # Check determinant.
    if gram_matrix.det() != determinant:
        return False
        
    # Compute genus.
    L = IntegralLattice(gram_matrix)
    genus = L.genus()
    
    # Check evenness.
    if not L.is_even():
        return False

    # Break down genus string.
    genus_tmp = genus_str.strip('[]').split(',')
    genus_2 = []
    genus_odd = {}
    for i in genus_tmp:
        if i == "U":
            genus_2.append(i)
            continue
        f = sage_eval(i)
        num, denom = f.as_integer_ratio()
        if denom % 2 == 0:
            genus_2.append(f)
        else:
            p, k = denom.is_prime_power(get_data=True)
            if k == 0:
                raise ValueError('Bad local genus symbol')
            if p not in genus_odd:
                genus_odd[p] = []
            genus_odd[p].append(f)

    # Check odd part of genus.
    for p, l in genus_odd.items():
        loc = genus.local_symbol(p).symbol_tuple_list()
        if len(l) == 1:
            if loc[-1][-1] != legendre_symbol(l[0].numerator(), p):
                return False
        else:
            if len(loc) != 2 or loc[-1][0] != l[-1].denom().valuation(p) or len(l) != 2:
                return False
            if loc[-1][-1] != legendre_symbol(l[0].numerator() * l[1].numerator(), p):
                 return False

    # Check 2-part of genus.
    loc = genus.local_symbol(2).symbol_tuple_list()
    if genus_2 == ["U"]:
        if len(loc) != 2 or loc[-1] != [1, 2, 3, 0, 0]:
            return False
    elif len(genus_2) == 1:
        num, q = genus_2[0].as_integer_ratio()
        if len(loc) != 2 or loc[-1][0] != q.valuation(2):
            return False
        if (genus.local_symbol(2).gram_matrix()[-1,-1] // q) % q != num % q:
            return False
    else:
        if len(loc) == 2:
            if loc != [[0,6,3,0,0],[1,2,1,1,2]]:
                return False
        else:
            if len(loc) != 3 or loc[1][0] != 1 or loc[2][0] != 2 or genus_2[0].denom() != 4 or genus_2[1].denom() != 2:
                return False
            mat = genus.local_symbol(2).canonical_symbol()
            if mat[0][2] != (-1)^len(genus_odd):
                return False
            a = genus_2[0].numer()
            b = genus_2[1].numer()
            if mat[1][4] - a + b != 2:
                return False
    return True
    
def verify_certificate(gram_matrix):
    # Check symmetry.
    if not gram_matrix.is_symmetric():
        return False

    # Check positive definiteness.
    if not gram_matrix.is_positive_definite():
        return False

    # Check trace of matrix times T.
    E8 = matrix(ZZ, 8, 8, [2, -1, 0, 0, 0, 0, 0, 0, -1, 2, -1, 0, 0, 0, 0, 0, 0, -1, 2, -1, -1, 0, 0, 0, 0, 0, -1, 2, 0, 0, 0, 0, 0, 0, -1, 0, 2, -1, 0, 0, 0, 0, 0, 0, -1, 2, -1, 0, 0, 0, 0, 0, 0, -1, 2, -1, 0, 0, 0, 0, 0, 0, -1, 2])
    d = gram_matrix.dimensions()[0]
    T = E8.submatrix(0, 0, d, d)
    if (T * gram_matrix).trace() != 24:
        return False

    # Check that the last basis vector is a root.    
    if gram_matrix[d-1, d-1] != 2:
        return False
        
    # Check that the associated vector is not a root of T.
    v = gram_matrix.submatrix(0, d-1, d, 1)
    if (v.transpose() * T * v)[0, 0] == 2:
        return False

    return True

def main(filename):
    with open(filename, newline='') as f:
        for row in csv.DictReader(f):
            rank = ZZ(row['rank'])
            determinant = QQ(row['determinant'])
            genus_str = row['genus']
            gram_matrix = matrix_from_upper_triangular(sage_eval(row['matrix']))
            print()
            print(latex(gram_matrix))
            verify = verify_genus(determinant, genus_str, gram_matrix) and verify_certificate(gram_matrix)
            if verify:
                print(f'Verified determinant {determinant} genus {genus_str}')
            else:
                print(f'Failed to verify determinant {determinant} genus {genus_str}')
            

if __name__ == '__main__':
    default = Path(__file__).resolve().parents[1] / 'data/farkas_certificates/certificates.csv'
    main(sys.argv[1] if len(sys.argv) > 1 else default)
