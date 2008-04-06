"""
Dense matrice over multivariate polynomials over fields.

This implementation inherits from Matrix_generic_dense, i.e. it is not
optimized for speed only some methods were added.

AUTHOR: Martin Albrecht <malb@informatik.uni-bremen.de>
"""

#*****************************************************************************
#
#   SAGE: System for Algebra and Geometry Experimentation
#
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include "sage/ext/python.pxi"
include "sage/libs/singular/singular-cdefs.pxi"

from sage.libs.singular.singular cimport Conversion

cdef Conversion co
co = Conversion()

from sage.matrix.matrix_generic_dense cimport Matrix_generic_dense
from sage.matrix.matrix2 cimport Matrix

from sage.rings.polynomial.multi_polynomial_libsingular cimport MPolynomial_libsingular
from sage.rings.polynomial.multi_polynomial_libsingular cimport MPolynomialRing_libsingular
from sage.rings.polynomial.polynomial_singular_interface import can_convert_to_singular

cdef class Matrix_mpolynomial_dense(Matrix_generic_dense):
    """
    Dense matrix over a multivariate polynomial ring over a field.
    """

    def echelon_form(self, algorithm="row_reduction", **kwds):
        """
        Return an echelon form of self depending on chosen algorithm.

        By default only a usual row reduction with no divisions or
        column swaps is returned.

        If Gauss-Bareiss algorithm is chosen, column swaps are
        recorded and can be retrieved via E.swapped_columns().

        INPUT:
            algorithm -- string, which algorithm to use (default: 'row_reduction')
                         'row_reduction' (default) -- reduce as far as possible, only divide
                                   by constant entries
                         'frac' -- reduced echelon form over fraction field
                         'bareiss' -- fraction free Gauss-Bareiss algorithm with column swaps

        OUTPUT:
            matrix -- A row echelon form of A depending on the
            chosen algorithm, as an immutable matrix.  Note that self
            is *not* changed by this command.  Use A.echelonize() to
            change A in place.

        EXAMPLES:
            sage: P.<x,y> = MPolynomialRing(GF(127),2)
            sage: A = matrix(P,2,2,[1,x,1,y])
            sage: A
            [1 x]
            [1 y]


            sage: A.echelon_form()
            [     1      x]
            [     0 -x + y]


        The reduced row echelon form over the fraction field is as follows:
            sage: A.echelon_form('frac') # over fraction field
            [1 0]
            [0 1]

        Alternatively, the Gauss-Bareiss algorithm may be chosen.

            sage: E = A.echelon_form('bareiss'); E
            [    1     y]
            [    0 x - y]

        After the application of the Gauss-Bareiss algorithm the
        swapped columns may inspected.

            sage: E.swapped_columns(), E.pivots()
            ([0, 1], [0, 1])

            sage: A.swapped_columns(), A.pivots()
            (None, [0, 1])

            Another approach is to row reduce as far as possible.

            sage: A.echelon_form('row_reduction')
            [     1      x]
            [     0 -x + y]

        """
        x = self.fetch('echelon_form_'+algorithm)
        if x is not None: return x

        if  algorithm == "frac":
            E = self.matrix_over_field()
            E.echelonize(**kwds)
        else:
            E = self.copy()
            E.echelonize(algorithm=algorithm, **kwds)

        E.set_immutable()  # so we can cache the echelon form.
        self.cache('echelon_form_'+algorithm, E)

        if algorithm == "frac":
            self.cache('pivots', E.pivots())
        elif algorithm == "bareiss":
            l = E.swapped_columns()
            self.cache('pivots', sorted(l))
        elif algorithm == "row_reduction":
            pass

        return E

    def echelonize(self, algorithm="row_reduction", **kwds):
        r"""
        Transform self into a matrix in echelon form over the same
        base ring as self.

        If Gauss-Bareiss algorithm is chosen, column swaps are
        recorded and can be retrieved via self.swapped_columns().

        INPUT:
            algorithm -- string, which algorithm to use (default: 'row_reduction')
                   'row_reduction' -- reduce as far as possible,
                                      only divide by constant entries
                   'bareiss' -- fraction free Gauss-Bareiss algorithm with column swaps

        EXAMPLES:
            sage: P.<x,y> = MPolynomialRing(QQ,2)
            sage: A = matrix(P,2,2,[1/2,x,1,3/4*y+1])
            sage: A
            [      1/2         x]
            [        1 3/4*y + 1]

            sage: B = A.copy(); B.echelonize('bareiss'); B
            [              1       3/4*y + 1]
            [              0 x - 3/8*y - 1/2]

            sage: B = A.copy(); B.echelonize('row_reduction'); B
            [               1              2*x]
            [               0 -2*x + 3/4*y + 1]

            sage: P.<x,y> = MPolynomialRing(QQ,2)
            sage: A = matrix(P,2,3,[2,x,0,3,y,1]); A
            [2 x 0]
            [3 y 1]

            sage: E = A.echelon_form('bareiss'); E
            [1 3 y]
            [0 2 x]

            sage: E.swapped_columns()
            [2, 0, 1]

            sage: A.pivots()
            [0, 1, 2]

        """
        self.check_mutability()

        if self._nrows == 0 or self._ncols == 0:
            self.cache('in_echelon_form_'+algorithm, True)
            self.cache('rank', 0)
            self.cache('pivots', [])
            return

        x = self.fetch('in_echelon_form_'+algorithm)
        if not x is None: return  # already known to be in echelon form

        if algorithm == 'bareiss':
            self._echelonize_gauss_bareiss()
        elif algorithm == 'row_reduction':
            self._echelonize_row_reduction()
        else:
            raise ValueError, "Unknown algorithm '%s'"%algorithm

    def _echelonize_gauss_bareiss(self):
        """
        Transform self into a matrix in echelon form over the same
        base ring as self using the fraction free Gauss-Bareiss
        algorithm with column swaps.

        The performed column swaps can be accessed via
        self.swapped_columns().

        ALGORITHM: Uses libSINGULAR or SINGULAR
        """
        cdef int r,c
        cdef intvec *iv = NULL
        cdef ideal *res = NULL
        cdef matrix *i = NULL
        cdef ideal *ii = NULL
        cdef R = self.base_ring()
        cdef ring *_ring = NULL
        cdef int *ivv = NULL
        cdef MPolynomial_libsingular p

        x = self.fetch('in_echelon_form_bareiss')
        if not x is None: return  # already known to be in echelon form

        if PY_TYPE_CHECK(self.base_ring(), MPolynomialRing_libsingular):
            self.check_mutability()
            self.clear_cache()

            i = mpNew(self._ncols, self._nrows)
            _ring = (<MPolynomialRing_libsingular>R)._ring
            rChangeCurrRing(_ring)

            # we are transposing here also
            for r from 0 <= r < self._nrows:
                for c from 0 <= c < self._ncols:
                    p = <MPolynomial_libsingular>self.get_unsafe(r,c)
                    i.m[self._nrows * c + r] = p_Copy(p._poly, _ring)

            ii = idMatrix2Module(i) # kills i

            # this is actually a sparse implementation
            smCallNewBareiss(ii,0,0,res,&iv)

            ivv = iv.ivGetVec()
            l = []
            for r from 0 <= r < iv.rows():
                l.append(int(ivv[r]-1))

            # clear matrix
            for r from 0 <= r < self._nrows:
                for c from 0 <= c < self._ncols:
                    self.set_unsafe(r,c,R._zero_element)

            # we are transposing here also
            for r from 0 <= r < IDELEMS(res):
                for c from 0 <= c < res.rank:
                    p = co.new_MP(R, pTakeOutComp1(&res.m[r], c+1))
                    self.set_unsafe(r,c,p)

            self.cache('in_echelon_form_bareiss',True)
            self.cache('rank', IDELEMS(res))
            self.cache('pivots', range(IDELEMS(res)))
            self.cache('swapped_columns', l)

            id_Delete(&res,_ring)
            id_Delete(&ii, _ring)
            delete(iv)

        elif can_convert_to_singular(self.base_ring()):

            self.check_mutability()
            self.clear_cache()

            E,l = self._singular_().bareiss()._sage_(self.base_ring())

            # clear matrix
            for r from 0 <= r < self._nrows:
                for c from 0 <= c < self._ncols:
                    self.set_unsafe(r,c,R._zero_element)

            for r from 0 <= r < E.nrows():
                for c from 0 <= c < E.ncols():
                    self.set_unsafe(r,c, E[r,c])

            self.cache('in_echelon_form_bareiss',True)
            self.cache('rank', E.nrows())
            self.cache('pivots', range(E.nrows()))
            self.cache('swapped_columns', l)

        else:

            raise NotImplementedError, "cannot apply Gauss-Bareiss algorithm over this base ring"

    def _echelonize_row_reduction(self):
        r"""
        Transform self to a matrix in row reduced form as far as this
        is possible, i.e. only perform division by constant elements.

        EXAMPLES:

            If all entries are constant, then this method performs the
            same operations as self.echelon_form('default').

            sage: P.<x0,x1,y0,y1> = PolynomialRing(GF(127),4)
            sage: A = Matrix(P,4,4,[-14,0,45,-55,-61,-16,0,0,0,0,25,-62,-22,0,52,0]); A
            [-14   0  45 -55]
            [-61 -16   0   0]
            [  0   0  25 -62]
            [-22   0  52   0]

            sage: E1 = A.echelon_form()
            sage: E2 = A.echelon_form('row_reduction')
            sage: E1 == E2
            True

            If no entries are constant, nothing happens:

            sage: P.<x0,x1,y0,y1> = PolynomialRing(GF(2),4)
            sage: A = Matrix(P,2,2,[x0,y0,x0,y0]); A
            [x0 y0]
            [x0 y0]

            sage: B = A.copy()
            sage: A.echelonize('row_reduction') # modifies self
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
            sage: B.echelonize('row_reduction'); B
            [                 1                  0                  0                  0     x0*y0 + x1 + 1]
            [                 0                  1                  0                  0              x0*y0]
            [                 0                  0                  1                  0    x0*y0 + x0 + x1]
            [                 0                  0                  0                  1         x0*y0 + x0]
            [                 0                  0                  0                  0            x0 + y1]
            [                 0                  0                  0                  0        x1 + y0 + 1]
            [                 0                  0                  0                  0         x0*y1 + x0]
            [                 0                  0                  0                  0              x1*y0]
            [                 0                  0                  0                  0 x0*y0 + x1*y1 + x0]

            This is the same result as SINGULAR's rowred command returns.

            sage: E = A._singular_().rowred()._sage_(P)
            sage: E == B
            True


        ALGORITHM: Gaussian elimination with division limited to
        constant entries. Based on SINGULAR's rowred command.

        """
        from sage.matrix.constructor import matrix

        cdef int c, r, i, j, rc, start_row, nr, nc

        x = self.fetch('in_echelon_form_row_reduction')
        if not x is None: return  # already known to be in echelon form

        nr,nc = self.nrows(),self.ncols()
        F = self.base_ring().base_ring()
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

                d = d._parent(0)
                for i from start_row <= i < nr:
                    for j from c+1 <= j < nc:
                        if self.get_unsafe(i,j).is_constant():
                            d.set_unsafe(i,j, self.get_unsafe(i,j).constant_coefficient())

        self.cache('in_echelon_form_row_reduction',True)


    def swapped_columns(self):
        """
        Return a tuple representing the column swaps during the last
        application of the Gauss-Bareiss algorithms (see
        self.echelon_form() for details).

        The tuple as length equal to the rank of self and the value at
        the $i$-th position indicates the source column which was put
        as the $i$-th column.

        If no Gauss-Bareiss reduction was performed yet, None is
        returned.
        """
        return self.fetch('swapped_columns')

