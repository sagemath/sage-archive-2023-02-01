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
            if R.tracks_precision():
                scalar = scalar.lift_to_maximal_precision()
            S.add_multiple_of_row(i,piv,scalar,piv+1)
            if transformation:
                left.add_multiple_of_row(i,piv,scalar)
        if transformation:
            left.rescale_row(piv,inv)
            for j in range(piv+1,m):
                scalar = -inv * (S[piv,j] >> curval)
                if R.tracks_precision():
                    scalar = scalar.lift_to_maximal_precision()
                right.add_multiple_of_column(j,piv,scalar)

    if transformation:
        if R.tracks_precision():
            prec = min([ x.precision_absolute() for x in M.list() ])
            if prec is not Infinity:
                prec -= curval
            left = left.apply_map(lambda x: x.add_bigoh(prec))
        return smith, left, right
    else:
        return smith


def determinant(M):
    """
    Helper function for the computation of the determinant

    See also :meth:`_matrix_determinant`
    """
    n = M.nrows()
    S = M.parent()(M.list())
    R = M.base_ring()

    sign = 1; det = R(1)
    valdet = 0; val = -Infinity
    for piv in range(n):
        curval = Infinity
        for i in range(piv,n):
            for j in range(piv,n):
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
            if R.tracks_precision():
                return R(0, valdet + (n-piv)*val)
            else:
                return R(0)

        valdet += val
        S.swap_rows(pivi,piv)
        if pivi > piv: sign = -sign
        S.swap_columns(pivj,piv)
        if pivj > piv: sign = -sign

        det *= S[piv,piv]
        inv = ~(S[piv,piv] >> val)
        for i in range(piv+1,n):
            scalar = -inv * (S[i,piv] >> val)
            if R.tracks_precision():
                scalar = scalar.lift_to_maximal_precision()
            S.add_multiple_of_row(i,piv,scalar)

    if R.tracks_precision():
        relprec = +Infinity
        for i in range(n):
            prec = Infinity
            for j in range(n):
                p = S[i,j].precision_absolute()
                if p < prec: prec = p
            prec -= S[i,i].valuation()
            if prec < relprec: relprec = prec
        return (sign*det).add_bigoh(valdet+relprec)
    else:
        return sign*det
