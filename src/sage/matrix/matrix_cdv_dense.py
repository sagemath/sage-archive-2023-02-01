"""
Matrices over complete discrete valuation rings/fields
"""

from sage.matrix.matrix_generic_dense import Matrix_generic_dense

from sage.rings.infinity import Infinity
from sage.rings.padics.precision_error import PrecisionError




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
