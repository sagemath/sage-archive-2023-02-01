"""
Matrices over complete discrete valuation rings/fields
"""

from sage.matrix.matrix_generic_dense import Matrix_generic_dense

from sage.rings.infinity import Infinity
from sage.rings.padics.precision_error import PrecisionError


def smith_normal_form(M, transformation):
    """
    Helper method for the computation of the Smith normal form

    See also :meth:`_matrix_smith_form` and :meth:`integral_smith_form`
    """
    n = M.nrows()
    m = M.ncols()
    S = M.parent()(M.list())
    smith = M.parent()(0)
    R = M.base_ring()
    if R.tracks_precision():
        precM = min([ x.precision_absolute() for x in M.list() ])

    if transformation:
        from sage.matrix.special import identity_matrix
        left = identity_matrix(R,n)
        right = identity_matrix(R,m)

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

        if R.tracks_precision() and precM is not Infinity and val >= precM:
            raise PrecisionError("Not enough precision to compute Smith normal form")

        if val is Infinity:
            break

        S.swap_rows(pivi,piv)
        S.swap_columns(pivj,piv)
        if transformation:
            left.swap_rows(pivi,piv)
            right.swap_columns(pivj,piv)

        smith[piv,piv] = R(1) << val
        inv = ~(S[piv,piv] >> val)
        for i in range(piv+1,n):
            scalar = -inv * (S[i,piv] >> val)
            if R.tracks_precision():
                scalar = scalar.lift_to_maximal_precision()
            S.add_multiple_of_row(i,piv,scalar,piv+1)
            if transformation:
                left.add_multiple_of_row(i,piv,scalar)
        if transformation:
            left.rescale_row(piv,inv)
            for j in range(piv+1,m):
                scalar = -inv * (S[piv,j] >> val)
                if R.tracks_precision():
                    scalar = scalar.lift_to_maximal_precision()
                right.add_multiple_of_column(j,piv,scalar)

    if transformation:
        if R.tracks_precision() and precM is not Infinity:
            left = left.apply_map(lambda x: x.add_bigoh(precM-val))
        return smith, left, right
    else:
        return smith


def determinant(M):
    """
    Helper function for the computation of the determinant

    See also :meth:`_matrix_determinant`
    """
    n = M.nrows()

    # For 2x2 matrices, we use the formula
    if n == 2:
        return M[0,0]*M[1,1] - M[0,1]*M[1,0]

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
        relprec_neg = 0
        for i in range(n):
            prec = Infinity
            for j in range(n):
                p = S[i,j].precision_absolute()
            prec -= S[i,i].valuation()
            if prec < relprec: relprec = prec
            if prec < 0: relprec_neg += prec
        if relprec_neg < 0: relprec = relprec_neg
        return (sign*det).add_bigoh(valdet+relprec)
    else:
        return sign*det


class Matrix_cdvr_dense(Matrix_generic_dense):
    pass


class Matrix_cdvf_dense(Matrix_generic_dense):
    def integral_smith_form(self, transformation=True):
        """
        Return the integral Smith normal form of this matrix.

        INPUT:

        - ``transformation`` -- a boolean (default: True)
          Indicates whether the transformation matrices are returned

        NOTE:

        The integral Smith decomposition of a matrix `M`
        defined over a complete discrete valuation field
        is a writing of the form `L*M*R = S` where:

        - `L` and `R` are invertible matrices over the ring of
          integers

        - the only non-vanishing entries of `S` are located on
          the diagonal (through `S` might be not a square matrix)

        - if `d_i` denotes the `(i,i)` entry of `S`, then `d_i`
          divides `d_{i+1}` in the ring of integers for all `i`.

        The `d_i`'s are uniquely determined provided that they are
        normalized so that they are all either `0` or a power of the 
        distinguished uniformizer of the base ring.

        EXAMPLES::

            sage: A = Qp(5, prec=10, print_mode="digits")
            sage: M = matrix(A, 2, 2, [2, 7, 1, 6])

            sage: S, L, R = M.integral_smith_form()
            sage: S
            [ ...1     0]
            [    0 ...10]
            sage: L
            [...222222223          ...]
            [...444444444         ...2]
            sage: R
            [         ...1 ...2222222214]
            [            0          ...1]

        If not needed, it is possible to do not compute the
        transformations matrices L and R as follows::

            sage: M.integral_smith_form(transformation=False)
            [ ...1     0]
            [    0 ...10]

        This method works for rectangular matrices as well::

            sage: M = matrix(A, 3, 2, [2, 7, 1, 6, 3, 8])
            sage: S, L, R = M.integral_smith_form()
            sage: S
            [ ...1     0]
            [    0 ...10]
            [    0     0]
            sage: L
            [...222222223          ...          ...]
            [...444444444         ...2          ...]
            [...444444443         ...1         ...1]
            sage: R
            [         ...1 ...2222222214]
            [            0          ...1]

        An error is raised if the precision on the entries is
        not enough to determine the Smith normal form::

            sage: M = matrix(A, 2, 2, [1, 1, 1, 1])
            sage: M.integral_smith_form()
            Traceback (most recent call last):
            ...
            PrecisionError: Not enough precision to compute Smith normal form

        TESTS::

        We check that Smith decomposition works over various rings::

            sage: from sage.rings.padics.precision_error import PrecisionError
            sage: ring1 = Qp(7,10)
            sage: ring2 = Qq(7^2,names='a')
            sage: ring3 = Qp(7).extension(x^3-7, names='pi')
            sage: ring4 = LaurentSeriesRing(GF(7), name='t')
            sage: for A in [ ring1, ring2, ring4 ]:  # ring3 causes troubles (see ticket #23464)
            ....:     for _ in range(10):
            ....:         M = random_matrix(A,4)
            ....:         try:
            ....:             S, L, R = M.integral_smith_form()
            ....:         except PrecisionError:
            ....:             continue
            ....:         if L*M*R != S: raise RuntimeError
        """
        return smith_normal_form(self, transformation)
