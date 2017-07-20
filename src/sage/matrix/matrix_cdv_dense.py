"""
Helper methods for matrices over CDVR/CDVF
"""

from sage.rings.infinity import Infinity
from sage.rings.padics.precision_error import PrecisionError


def smith_normal_form(M, transformation):
    """
    Helper method for the computation of the Smith normal form

    See also :meth:`_matrix_smith_form`
    """
    n = M.nrows()
    m = M.ncols()
    S = M.parent()(M.list())
    smith = M.parent()(0)
    R = M.base_ring()
    if transformation:
        from sage.matrix.special import identity_matrix
        left = identity_matrix(R,n)
        right = identity_matrix(R,m)
    else:
        left = right = None
    val = -Infinity
    for piv in range(min(n,m)):
        curval = Infinity
        pivi = pivj = piv
        for i in range(piv,n):
            for j in range(piv,m):
                v = S[i,j].valuation()
                if v < curval:
                    pivi = i; pivj = j
                    curval = v
                    if v == val: break
            else:
                continue
            break
        val = curval

        if S[pivi,pivj] == 0:
            if curval is Infinity:
                break
            else:
                raise PrecisionError("Not enough precision to compute Smith normal form")

        S.swap_rows(pivi,piv)
        S.swap_columns(pivj,piv)
        if transformation:
            left.swap_rows(pivi,piv)
            right.swap_columns(pivj,piv)

        smith[piv,piv] = R(1) << curval
        inv = ~(S[piv,piv] >> curval)
        for i in range(piv+1,n):
            scalar = -inv * (S[i,piv] >> curval)
            scalar = scalar.lift_to_maximal_precision()
            S.add_multiple_of_row(i,piv,scalar,piv+1)
            if transformation:
                left.add_multiple_of_row(i,piv,scalar)
        if transformation:
            left.rescale_row(piv,inv)
            for j in range(piv+1,m):
                scalar = -inv * (S[piv,j] >> curval)
                scalar = scalar.lift_to_maximal_precision()
                right.add_multiple_of_column(j,piv,scalar)

    if transformation:
        prec = min([ x.precision_absolute() for x in M.list() ])
        if prec is not Infinity:
            prec -= curval
        left = left.apply_map(lambda x: x.add_bigoh(prec))
        return smith, left, right
    else:
        return smith
