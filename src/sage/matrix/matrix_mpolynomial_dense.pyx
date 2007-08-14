"""
Dense matrice over multivariate polynomials over fields.

This implementation inherits from Matrix_generic_dense, i.e. it is not
optimized for speed. Only some methods were added to
Matrix_generic_dense.

AUTHOR: Martin Albrecht <malb@informatik.uni-bremen.de>
"""

##############################################################################
#       Copyright (C) 2007  William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
##############################################################################

include "sage/ext/python.pxi"
cimport matrix_generic_dense

from sage.matrix.matrix2 cimport Matrix

cdef class Matrix_mpolynomial_dense(matrix_generic_dense.Matrix_generic_dense):   # dense or sparse
    """
    Dense matrix over a multivariate polynomial ring over a field.
    """

    def row_reduce(self):
        r"""
        Transform self to a matrix in row reduced form as far as this
        is possible, i.e. only perform division by constant elements.

        EXAMPLES:

            If all entries are constant, then this method performs the
            same operations as self.echelonize().

            sage: P.<x0,x1,y0,y1> = PolynomialRing(GF(127),4)
            sage: A = Matrix(P,4,4,[-14,0,45,-55,-61,-16,0,0,0,0,25,-62,-22,0,52,0]); A
            [-14   0  45 -55]
            [-61 -16   0   0]
            [  0   0  25 -62]
            [-22   0  52   0]

            sage: E = A.echelon_form()
            sage: A.row_reduce() # modifies self
            sage: E == A
            True

            If no entries are constant, nothing happens:

            sage: P.<x0,x1,y0,y1> = PolynomialRing(GF(2),4)
            sage: A = Matrix(P,2,2,[x0,y0,x0,y0]); A
            [x0 y0]
            [x0 y0]

            sage: B = A.copy()
            sage: A.row_reduce() # modifies self
            sage: B == A
            True

            A more interesting example:

            sage: P.<x0,x1,y0,y1> = PolynomialRing(GF(2),4)
            sage: l = [1, 1, 1, 1,     1, \
                       0, 1, 0, 1,    x0, \
                       0, 0, 1, 1,    x1, \
                       1, 1, 0, 0,    y0, \
                       0, 1, 0, 1,    y1, \
                       0, 1, 0, 0, x0*y0, \
                       0, 1, 0, 1, x0*y1, \
                       0, 0, 0, 0, x1*y0, \
                       0, 0, 0, 1, x1*y1]
            sage: A = Matrix(P,9,5,l)
            sage: B = A.copy()
            sage: B.row_reduce(); B
            [                 1                  0                  0                  0     x0*y0 + x1 + 1]
            [                 0                  1                  0                  0              x0*y0]
            [                 0                  0                  1                  0    x0*y0 + x0 + x1]
            [                 0                  0                  0                  1         x0*y0 + x0]
            [                 0                  0                  0                  0            x0 + y1]
            [                 0                  0                  0                  0        x1 + y0 + 1]
            [                 0                  0                  0                  0         x0*y1 + x0]
            [                 0                  0                  0                  0              x1*y0]
            [                 0                  0                  0                  0 x0*y0 + x1*y1 + x0]

            This is the same result SINGULAR's rowred command returns.

            sage: A._singular_().rowred()._sage_(P) == B
            True


        ALGORITHM: Gaussian elimination with division limited to
        constant entries.  Based on SINGULAR's rowred command.

        """
        from sage.matrix.constructor import matrix

        cdef int c, r, i, j, rc, start_row, nr, nc

        nr,nc = self.nrows(),self.ncols()
        F = self.base_ring().base_ring()
        if not F.is_field():
            raise TypeError, "row reduction only supported for polynomial rings over fields"
        cdef Matrix d = matrix(F,nr,nc)
        start_row = 0

        for r from 0 <= r < nr:
            for c from 0 <= c < nc:
                p = self.get_unsafe(r,c)
                if p.is_constant():
                    d.set_unsafe(r, c, p.constant_coefficient())

        for c from 0 <= c < nc:
            r = -1
            for rc from start_row <= rc < nr:
                if d.get_unsafe(rc, c):
                    r = rc
                    break
            if r!=-1:
                a_inverse = ~self.get_unsafe(r,c)
                self.rescale_row_c(r, a_inverse , c)
                self.swap_rows_c(r, start_row)

                for i from 0 <= i < nr:
                    if i != start_row:
                        minus_b = -self.get_unsafe(i,c)
                        self.add_multiple_of_row(i, start_row, minus_b, 0)

                start_row +=1

                d = 0 * d
                for i from start_row <= i < nr:
                    for j from c+1 <= j < nc:
                        if self.get_unsafe(i,j).is_constant():
                            d.set_unsafe(i,j, self.get_unsafe(i,j).constant_coefficient())


