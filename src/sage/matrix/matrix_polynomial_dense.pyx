"""
Dense matrices over univariate polynomials over fields

This implementation inherits from Matrix_generic_dense, i.e. it is not
optimized for speed only some methods were added.

AUTHOR:

* Kwankyu Lee <ekwankyu@gmail.com>
"""

#*****************************************************************************
#       Copyright (C) 2016 Kwankyu Lee <ekwankyu@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.polynomial.multi_polynomial_libsingular cimport new_MP

from sage.matrix.matrix_generic_dense cimport Matrix_generic_dense
from sage.matrix.matrix2 cimport Matrix

from sage.rings.polynomial.multi_polynomial_libsingular cimport MPolynomial_libsingular
from sage.rings.polynomial.multi_polynomial_libsingular cimport MPolynomialRing_libsingular
from sage.rings.polynomial.polynomial_singular_interface import can_convert_to_singular

from sage.libs.singular.function import singular_function

cdef class Matrix_polynomial_dense(Matrix_generic_dense):
    """
    Dense matrix over a univariate polynomial ring over a field.
    """

    def is_weak_popov(self):
        r"""
        Return ``True`` if the matrix is in weak Popov form.

        OUTPUT:

        A matrix over an ordered ring is in weak Popov form if all
        leading positions are different [MS2003]_. A leading position
        is the position `i` in a row with the highest order (for
        polynomials this is the degree), for multiple entries with
        equal but highest order the maximal `i` is chosen (which is
        the furthest to the right in the matrix).

        .. WARNING::

            This implementation only works for objects implementing a degree
            function. It is designed to work for polynomials.

        EXAMPLES:

        A matrix with the same leading position in two rows is not in weak
        Popov form. ::

            sage: PF = PolynomialRing(GF(2^12,'a'),'x')
            sage: A = matrix(PF,3,[x,x^2,x^3,x^2,x^2,x^2,x^3,x^2,x])
            sage: A.is_weak_popov()
            False

        If a matrix has different leading positions, it is in weak Popov
        form. ::

            sage: B = matrix(PF,3,[1,1,x^3,x^2,1,1,1,x^2,1])
            sage: B.is_weak_popov()
            True

        Weak Popov form is not restricted to square matrices. ::

            sage: PF = PolynomialRing(GF(7),'x')
            sage: D = matrix(PF,2,4,[x^2+1,1,2,x,3*x+2,0,0,0])
            sage: D.is_weak_popov()
            False

        Even a matrix with more rows than cols can still be in weak Popov
        form. ::

            sage: E = matrix(PF,4,2,[4*x^3+x,x^2+5*x+2,0,0,4,x,0,0])
            sage: E.is_weak_popov()
            True

        But a matrix with less cols than non zero rows is never in weak
        Popov form. ::

            sage: F = matrix(PF,3,2,[x^2,x,x^3+2,x,4,5])
            sage: F.is_weak_popov()
            False

        TESTS:

        A matrix to check if really the rightmost value is taken. ::

            sage: F = matrix(PF,2,2,[x^2,x^2,x,5])
            sage: F.is_weak_popov()
            True

        .. SEEALSO::

            - :meth:`weak_popov_form <sage.matrix.matrix2.weak_popov_form>`

        AUTHOR:

        - David Moedinger (2014-07-30)
        """
        t = set()
        for r in range(self.nrows()):
            max = -1
            for c in range(self.ncols()):
                if self[r, c].degree() >= max:
                    max = self[r, c].degree()
                    p = c
            if not max == -1:
                if p in t:
                    return False
                t.add(p)
        return True

    def weak_popov_form(self, transformation=None, ascend=None, old_call=True):
        """
        Return a matrix in weak Popov form which is row space-equivalent to
        the input matrix.

        A matrix is in weak Popov form if the leading positions of
        the non-zero rows are all different. The leading position of a row is
        the right-most position whose entry has maximal degree of the entries in
        that row.

        .. WARNING::

            This function currently does **not** compute the weak Popov form of a
            matrix, but rather a row reduced form (which is a slightly weaker
            requirement). See :meth:`row_reduced_form`.

        INPUT:

        - `transformation` - A boolean (default: `True`). If this is set to
          ``True``, the transformation matrix `U` will be returned as well: this
          is an invertible matrix over `k(x)` such that ``self`` equals `UW`,
          where `W` is the output matrix.

          Warning: the default of `transformation` will soon be set to ``False``,
          see :trac:`16896`. Until then, this function will print a deprecation
          warning as long as `transformation` was not explicitly set to ``True``
          or ``False``.

        - `ascend` - Deprecated and has no effect.

        - `old_call` - For backwards compatibility until the old calling
          convention will be deprecated (default: `True`). If `True`, then
          return `(W,U,d)`, where `U` is as when `transformation = True`, and
          `d` is a list of the degrees of the rows of `W`.
        """
        from sage.misc.superseded import deprecation
        deprecation(16888, "This function currently does *not* compute a weak Popov form, "
        "but rather a row reduced form. This function will soon be fixed (see Ticket #16742).")

        return self.row_reduced_form(transformation=transformation,
                ascend=ascend, old_call=old_call)

    def _weak_popov_form(self, transformation=False):
        """
        Return a matrix in weak Popov form which is row space-equivalent to
        the input matrix, if the input is over `k[x]` or `k(x)`.

        A matrix is in weak Popov form if the (row-wise) leading positions of
        the non-zero rows are all different. The leading position of a row is
        the right-most position whose entry has maximal degree of the entries in
        that row.

        INPUT:

        - ``transformation`` -- (default: `True`). If this is set to
          ``True``, the transformation matrix `U` will be returned as well: this
          is an invertible matrix over `k(x)` such that ``self`` equals `UW`,
          where `W` is the output matrix.
        """
        from sage.matrix.weak_popov import weak_popov_form_mulders_storjohann
        return weak_popov_form_mulders_storjohann(self, transformation)

    def row_reduced_form(self, transformation=None, ascend=None, old_call=True):
        r"""
        Return a row reduced form of this matrix.

        A matrix `M` is row reduced if the leading term matrix has the same rank
        as `M`. The leading term matrix of a polynomial matrix `M_0` is the matrix
        over `k` whose `(i,j)`'th entry is the `x^{d_i}` coefficient of `M_0[i,j]`,
        where `d_i` is the greatest degree among polynomials in the `i`'th row of `M_0`.

        INPUT:

        - `transformation` - A boolean (default: ``False``). If this is set to
          ``True``, the transformation matrix `U` will be returned as well: this
          is an invertible matrix over `k(x)` such that ``self`` equals `UW`,
          where `W` is the output matrix.

        - `ascend` - Deprecated and has no effect.

        - `old_call` - For backwards compatibility until the old calling
          convention will be deprecated (default: `True`). If `True`, then
          return `(W,U,d)`, where `U` is as when `transformation = True`, and
          `d` is a list of the degrees of the rows of `W`.

        OUTPUT:

        - `W` - a matrix over the same ring as `self` (i.e. either `k(x)` or
          `k[x]` for a field `k`) giving a row reduced form of ``self``.

        EXAMPLES::

            sage: R.<t> = GF(3)['t']
            sage: K = FractionField(R)
            sage: M = matrix([[(t-1)^2],[(t-1)]])
            sage: M.row_reduced_form(transformation=False, old_call=False)
            [    0]
            [t + 2]

        If ``self`` is an `n \times 1` matrix with at least one non-zero entry,
        `W` has a single non-zero entry and that entry is a scalar multiple of
        the greatest-common-divisor of the entries of ``self``.

        ::

            sage: M1 = matrix([[t*(t-1)*(t+1)],[t*(t-2)*(t+2)],[t]])
            sage: output1 = M1.row_reduced_form(transformation=False, old_call=False)
            sage: output1
            [0]
            [0]
            [t]

        The following is the first half of example 5 in [Hes2002]_ *except* that we
        have transposed ``self``; [Hes2002]_ uses column operations and we use row.

        ::

            sage: R.<t> = QQ['t']
            sage: M = matrix([[t^3 - t,t^2 - 2],[0,t]]).transpose()
            sage: M.row_reduced_form(transformation=False, old_call=False)
            [      t    -t^2]
            [t^2 - 2       t]

        The next example demonstrates what happens when ``self`` is a zero matrix.

        ::

            sage: R.<t> = GF(5)['t']
            sage: M = matrix(R, 2, [0,0,0,0])
            sage: M.row_reduced_form(transformation=False, old_call=False)
            [0 0]
            [0 0]

        In the following example, ``self`` has more rows than columns. Note also
        that the output is row reduced but not in weak Popov form (see
        :meth:`weak_popov_form`).

        ::

            sage: R.<t> = QQ['t']
            sage: M = matrix([[t,t,t],[0,0,t]])
            sage: M.row_reduced_form(transformation=False, old_call=False)
            [t t t]
            [0 0 t]

        The last example shows the usage of the transformation parameter.

        ::

            sage: Fq.<a> = GF(2^3)
            sage: Fx.<x> = Fq[]
            sage: A = matrix(Fx,[[x^2+a,x^4+a],[x^3,a*x^4]])
            sage: W,U = A.row_reduced_form(transformation=True,old_call=False);
            sage: W,U
            (
            [(a^2 + 1)*x^3 + x^2 + a                       a]  [      1 a^2 + 1]
            [                    x^3                   a*x^4], [      0                 1]
            )
            sage: U*W == A
            True
            sage: U.is_invertible()
            True

        NOTES:

         - For consistency with LLL and other algorithms in Sage, we have opted
           for row operations; however, references such as [Hes2002]_ transpose and use
           column operations.

         - There are multiple weak Popov forms of a matrix, so one may want to
           extend this code to produce a Popov form (see section 1.2 of [V]).  The
           latter is canonical, but more work to produce.

        .. SEEALSO::

            :meth:`is_weak_popov <sage.matrix.matrix0.is_weak_popov>`

        REFERENCES:

        - [Hes2002]_
        - [Kal1980]_

        """
        from sage.matrix.matrix_misc import row_reduced_form

        from sage.misc.superseded import deprecation
        if ascend is not None:
            deprecation(16888,
            "row_reduced_form: The `ascend` argument is deprecated "
            "and has no effect (see Ticket #16742).")
        if old_call == True:
            deprecation(16888,
            "row_reduced_form: The old calling convention is deprecated. "
            "In the future, only the matrix in row reduced form will be returned. "
            "Set `old_call = False` for that behaviour now, and to avoid this "
            "message (see Ticket #16742).")

        get_transformation = False
        if transformation is None:
            deprecation(16888,
            "row_reduced_form: The `transformation` argument will soon change to have "
            "default value to `False` from the current default value `True`. For now, "
            "explicitly setting the argument to `True` or `False` will avoid this message.")
            get_transformation = True
        elif old_call == True or transformation == True:
            get_transformation = True

        W_maybe_U = self._row_reduced_form(get_transformation)

        if not old_call:
            return W_maybe_U
        else:
            W, U = W_maybe_U
            row_deg = lambda r: max([e.degree() for e in r])
            d = []
            from sage.rings.all import infinity
            for r in W.rows():
                d.append(row_deg(r))
                if d[-1] < 0:
                    d[-1] = -infinity
            return (W,U,d)

    def _row_reduced_form(self, transformation=False):
        """
        Return a row reduced form of this matrix.

        INPUT:

        - `transformation` -- (default: `False`). If this is `True`,
           the transformation matrix is output.

        OUTPUT:

        If `transformation` is `False`, the output is `W`, a row reduced form of `M`.

        If `transformation` is `True`, this function will output a pair `(W,N)` of two matrices:

        1. `W` - a row reduced form of `M`.
        2. `N` - a unimodular matrix satisfying `N * W = M`.

        EXAMPLES::

            sage: Fq.<a> = GF(2^3)
            sage: Fx.<x> = Fq[]
            sage: A = matrix(Fx,[[x^2+a,x^4+a],[x^3,a*x^4]])
            sage: A._row_reduced_form(transformation=True)
            (
            [(a^2 + 1)*x^3 + x^2 + a                       a]  [      1 a^2 + 1]
            [                    x^3                   a*x^4], [      0                 1]
            )
        """

        # determine whether M has polynomial or rational function coefficients
        R = self.base_ring()
        t = R.gen()

        from sage.matrix.constructor import matrix

        num = self
        r = [list(v) for v in num.rows()]

        if transformation:
            N = matrix(num.nrows(), num.nrows(), R(1)).rows()

        rank = 0
        num_zero = 0
        if num.is_zero():
            num_zero = len(r)
        while rank != len(r) - num_zero:
            # construct matrix of leading coefficients
            v = []
            for w in map(list, r):
                # calculate degree of row (= max of degree of entries)
                d = max([e.numerator().degree() for e in w])

                # extract leading coefficients from current row
                x = []
                for y in w:
                    if y.degree() >= d and d >= 0:   x.append(y.coefficients(sparse=False)[d])
                    else:                            x.append(0)
                v.append(x)
            l = matrix(v)

            # count number of zero rows in leading coefficient matrix
            # because they do *not* contribute interesting relations
            num_zero = 0
            for v in l.rows():
                is_zero = 1
                for w in v:
                    if w != 0:
                        is_zero = 0
                if is_zero == 1:
                    num_zero += 1

            # find non-trivial relations among the columns of the
            # leading coefficient matrix
            kern = l.kernel().basis()
            rank = num.nrows() - len(kern)

            # do a row operation if there's a non-trivial relation
            if not rank == len(r) - num_zero:
                for rel in kern:
                    # find the row of num involved in the relation and of
                    # maximal degree
                    indices = []
                    degrees = []
                    for i in range(len(rel)):
                        if rel[i] != 0:
                            indices.append(i)
                            degrees.append(max([e.degree() for e in r[i]]))

                    # find maximum degree among rows involved in relation
                    max_deg = max(degrees)

                    # check if relation involves non-zero rows
                    if max_deg != -1:
                        i = degrees.index(max_deg)
                        rel /= rel[indices[i]]

                        for j in range(len(indices)):
                            if j != i:
                                # do the row operation
                                v = []
                                for k in range(len(r[indices[i]])):
                                    v.append(r[indices[i]][k] + rel[indices[j]] * t**(max_deg-degrees[j]) * r[indices[j]][k])
                                r[indices[i]] = v

                                if transformation:
                                    # If the user asked for it, record the row operation
                                    v = []
                                    for k in range(len(N[indices[i]])):
                                        v.append(N[indices[i]][k] + rel[indices[j]] * t**(max_deg-degrees[j]) * N[indices[j]][k])
                                    N[indices[i]] = v

                        # remaining relations (if any) are no longer valid,
                        # so continue onto next step of algorithm
                        break

        A = matrix(R, r)
        if transformation:
            return (A, matrix(N))
        else:
            return A

#    def echelon_form(self, algorithm='row_reduction', **kwds):
#        """
#        Return an echelon form of ``self`` using chosen algorithm.
#
#        By default only a usual row reduction with no divisions or column swaps
#        is returned.
#
#        If Gauss-Bareiss algorithm is chosen, column swaps are recorded and can
#        be retrieved via :meth:`swapped_columns`.
#
#        INPUT:
#
#        - ``algorithm`` -- string, which algorithm to use (default:
#          'row_reduction'). Valid options are:
#
#            - ``'row_reduction'`` (default) -- reduce as far as
#              possible, only divide by constant entries
#
#            - ``'frac'`` -- reduced echelon form over fraction field
#
#            - ``'bareiss'`` -- fraction free Gauss-Bareiss algorithm
#              with column swaps
#
#        OUTPUT:
#
#        The row echelon form of A depending on the chosen algorithm,
#        as an immutable matrix.  Note that ``self`` is *not* changed
#        by this command. Use ``A.echelonize()``` to change `A` in
#        place.
#
#        EXAMPLES::
#
#            sage: P.<x,y> = PolynomialRing(GF(127), 2)
#            sage: A = matrix(P, 2, 2, [1, x, 1, y])
#            sage: A
#            [1 x]
#            [1 y]
#            sage: A.echelon_form()
#            [     1      x]
#            [     0 -x + y]
#
#
#        The reduced row echelon form over the fraction field is as follows::
#
#            sage: A.echelon_form('frac') # over fraction field
#            [1 0]
#            [0 1]
#
#        Alternatively, the Gauss-Bareiss algorithm may be chosen::
#
#            sage: E = A.echelon_form('bareiss'); E
#            [    1     y]
#            [    0 x - y]
#
#        After the application of the Gauss-Bareiss algorithm the swapped columns
#        may inspected::
#
#            sage: E.swapped_columns(), E.pivots()
#            ((0, 1), (0, 1))
#            sage: A.swapped_columns(), A.pivots()
#            (None, (0, 1))
#
#        Another approach is to row reduce as far as possible::
#
#            sage: A.echelon_form('row_reduction')
#            [     1      x]
#            [     0 -x + y]
#        """
#        x = self.fetch('echelon_form_'+algorithm)
#        if x is not None: return x
#
#        if  algorithm == "frac":
#            E = self.matrix_over_field()
#            E.echelonize(**kwds)
#        else:
#            E = self.__copy__()
#            E.echelonize(algorithm=algorithm, **kwds)
#
#        E.set_immutable()  # so we can cache the echelon form.
#        self.cache('echelon_form_'+algorithm, E)
#
#        if algorithm == "frac":
#            self.cache('pivots', E.pivots())
#        elif algorithm == "bareiss":
#            l = E.swapped_columns()
#            self.cache('pivots', tuple(sorted(l)))
#        elif algorithm == "row_reduction":
#            pass
#
#        return E
#
#    def pivots(self):
#        """
#        Return the pivot column positions of this matrix as a list of integers.
#
#        This returns a list, of the position of the first nonzero entry in each
#        row of the echelon form.
#
#        OUTPUT:
#
#        A list of Python ints.
#
#        EXAMPLES::
#
#            sage: matrix([PolynomialRing(GF(2), 2, 'x').gen()]).pivots()
#            (0,)
#            sage: K = QQ['x,y']
#            sage: x, y = K.gens()
#            sage: m = matrix(K, [(-x, 1, y, x - y), (-x*y, y, y^2 - 1, x*y - y^2 + x), (-x*y + x, y - 1, y^2 - y - 2, x*y - y^2 + x + y)])
#            sage: m.pivots()
#            (0, 2)
#            sage: m.rank()
#            2
#        """
#        x = self.fetch('pivots')
#        if not x is None:
#            return x
#
#        self.echelon_form('frac')
#        x = self.fetch('pivots')
#
#        if x is None:
#            raise RuntimeError("BUG: matrix pivots should have been set but weren't, matrix parent = '%s'"%self.parent())
#        return x
#
#    def echelonize(self, algorithm='row_reduction', **kwds):
#        """
#        Transform self into a matrix in echelon form over the same base ring as
#        ``self``.
#
#        If Gauss-Bareiss algorithm is chosen, column swaps are recorded and can
#        be retrieved via :meth:`swapped_columns`.
#
#        INPUT:
#
#        - ``algorithm`` -- string, which algorithm to use. Valid options are:
#
#            - ``'row_reduction'`` -- reduce as far as possible, only
#              divide by constant entries
#
#            - ``'bareiss'`` -- fraction free Gauss-Bareiss algorithm
#              with column swaps
#
#        EXAMPLES::
#
#            sage: P.<x,y> = PolynomialRing(QQ, 2)
#            sage: A = matrix(P, 2, 2, [1/2, x, 1, 3/4*y+1])
#            sage: A
#            [      1/2         x]
#            [        1 3/4*y + 1]
#
#            sage: B = copy(A)
#            sage: B.echelonize('bareiss'); B
#            [              1       3/4*y + 1]
#            [              0 x - 3/8*y - 1/2]
#
#            sage: B = copy(A)
#            sage: B.echelonize('row_reduction'); B
#            [               1              2*x]
#            [               0 -2*x + 3/4*y + 1]
#
#            sage: P.<x,y> = PolynomialRing(QQ, 2)
#            sage: A = matrix(P,2,3,[2,x,0,3,y,1]); A
#            [2 x 0]
#            [3 y 1]
#
#            sage: E = A.echelon_form('bareiss'); E
#            [1 3 y]
#            [0 2 x]
#            sage: E.swapped_columns()
#            (2, 0, 1)
#            sage: A.pivots()
#            (0, 1, 2)
#        """
#        self.check_mutability()
#
#        if self._nrows == 0 or self._ncols == 0:
#            self.cache('in_echelon_form_'+algorithm, True)
#            self.cache('rank', 0)
#            self.cache('pivots', tuple())
#            return
#
#        x = self.fetch('in_echelon_form_'+algorithm)
#        if not x is None:
#            return  # already known to be in echelon form
#
#        if algorithm == 'bareiss':
#            self._echelonize_gauss_bareiss()
#        elif algorithm == 'row_reduction':
#            self._echelonize_row_reduction()
#        else:
#            raise ValueError("Unknown algorithm '%s'"%algorithm)
#
#    def _echelonize_gauss_bareiss(self):
#        """
#        Transform this matrix into a matrix in upper triangular form over the
#        same base ring as ``self`` using the fraction free Gauss-Bareiss
#        algorithm with column swaps.
#
#        The performed column swaps can be accessed via
#        :meth:`swapped_columns`.
#
#        EXAMPLE::
#
#            sage: R.<x,y> = QQ[]
#            sage: C = random_matrix(R, 2, 2, terms=2)
#            sage: C
#            [-6/5*x*y - y^2 -6*y^2 - 1/4*y]
#            [  -1/3*x*y - 3        x*y - x]
#
#            sage: E = C.echelon_form('bareiss')     # indirect doctest
#            sage: E
#            [ -1/3*x*y - 3                                                          x*y - x]
#            [            0 6/5*x^2*y^2 + 3*x*y^3 - 6/5*x^2*y - 11/12*x*y^2 + 18*y^2 + 3/4*y]
#            sage: E.swapped_columns()
#            (0, 1)
#
#        ALGORITHM:
#
#        Uses libSINGULAR or SINGULAR
#        """
#        cdef R = self.base_ring()
#
#        x = self.fetch('in_echelon_form_bareiss')
#        if not x is None:
#            return  # already known to be in echelon form
#
#        if isinstance(self.base_ring(), MPolynomialRing_libsingular):
#
#            self.check_mutability()
#            self.clear_cache()
#
#            singular_bareiss = singular_function("bareiss")
#            E, l = singular_bareiss(self.T)
#
#            m = len(E)
#            n = len(E[0])
#
#            for r in xrange(m):
#                for c in range(n):
#                    self.set_unsafe(r, c, E[r][c])
#                for c in xrange(n, self.ncols()):
#                    self.set_unsafe(r, c, R._zero_element)
#            for r in xrange(m, self.nrows()):
#                for c in xrange(self.ncols()):
#                    self.set_unsafe(r, c, R._zero_element)
#
#            from sage.rings.all import ZZ
#            l = [ZZ(e-1) for e in l]
#
#            self.cache('in_echelon_form_bareiss',True)
#            self.cache('rank', len(E))
#            self.cache('pivots', tuple(range(len(E))))
#            self.cache('swapped_columns', tuple(l))
#
#        elif can_convert_to_singular(self.base_ring()):
#
#            self.check_mutability()
#            self.clear_cache()
#
#            E,l = self.T._singular_().bareiss()._sage_(self.base_ring())
#
#            # clear matrix
#            for r from 0 <= r < self._nrows:
#                for c from 0 <= c < self._ncols:
#                    self.set_unsafe(r,c,R._zero_element)
#
#            for r from 0 <= r < E.nrows():
#                for c from 0 <= c < E.ncols():
#                    self.set_unsafe(c,r, E[r,c])
#
#            self.cache('in_echelon_form_bareiss',True)
#            self.cache('rank', E.nrows())
#            self.cache('pivots', tuple(range(E.nrows())))
#            self.cache('swapped_columns', l)
#
#        else:
#
#            raise NotImplementedError("cannot apply Gauss-Bareiss algorithm over this base ring")
#
#    def _echelonize_row_reduction(self):
#        """
#        Transform this matrix to a matrix in row reduced form as far as this is
#        possible, i.e., only perform division by constant elements.
#
#        EXAMPLES:
#
#        If all entries are constant, then this method performs the same
#        operations as ``echelon_form('default')`` ::
#
#            sage: P.<x0,x1,y0,y1> = PolynomialRing(GF(127), 4)
#            sage: A = Matrix(P, 4, 4, [-14,0,45,-55,-61,-16,0,0,0,0,25,-62,-22,0,52,0]); A
#            [-14   0  45 -55]
#            [-61 -16   0   0]
#            [  0   0  25 -62]
#            [-22   0  52   0]
#            sage: E1 = A.echelon_form()
#            sage: E2 = A.echelon_form('row_reduction')  # indirect doctest
#            sage: E1 == E2
#            True
#
#        If no entries are constant, nothing happens::
#
#            sage: P.<x0,x1,y0,y1> = PolynomialRing(GF(2), 4)
#            sage: A = Matrix(P, 2, 2, [x0, y0, x0, y0]); A
#            [x0 y0]
#            [x0 y0]
#
#            sage: B = copy(A)
#            sage: A.echelonize('row_reduction') # modifies self
#            sage: B == A
#            True
#
#        A more interesting example::
#
#            sage: P.<x0,x1,y0,y1> = PolynomialRing(GF(2), 4)
#            sage: l = [1, 1, 1, 1,     1, \
#                       0, 1, 0, 1,    x0, \
#                       0, 0, 1, 1,    x1, \
#                       1, 1, 0, 0,    y0, \
#                       0, 1, 0, 1,    y1, \
#                       0, 1, 0, 0, x0*y0, \
#                       0, 1, 0, 1, x0*y1, \
#                       0, 0, 0, 0, x1*y0, \
#                       0, 0, 0, 1, x1*y1]
#            sage: A = Matrix(P, 9, 5, l)
#            sage: B = A.__copy__()
#            sage: B.echelonize('row_reduction'); B
#            [                 1                  0                  0                  0     x0*y0 + x1 + 1]
#            [                 0                  1                  0                  0              x0*y0]
#            [                 0                  0                  1                  0    x0*y0 + x0 + x1]
#            [                 0                  0                  0                  1         x0*y0 + x0]
#            [                 0                  0                  0                  0            x0 + y1]
#            [                 0                  0                  0                  0        x1 + y0 + 1]
#            [                 0                  0                  0                  0         x0*y1 + x0]
#            [                 0                  0                  0                  0              x1*y0]
#            [                 0                  0                  0                  0 x0*y0 + x1*y1 + x0]
#
#        This is the same result as SINGULAR's ``rowred`` command which
#        returns::
#
#            sage: E = A._singular_().rowred()._sage_(P)
#            sage: E == B
#            True
#
#        ALGORITHM:
#
#        Gaussian elimination with division limited to constant
#        entries. Based on SINGULAR's rowred commamnd.
#        """
#        from sage.matrix.constructor import matrix
#
#        cdef int c, r, i, j, rc, start_row, nr, nc
#
#        x = self.fetch('in_echelon_form_row_reduction')
#        if not x is None: return  # already known to be in echelon form
#
#        nr,nc = self.nrows(),self.ncols()
#        F = self.base_ring().base_ring()
#        cdef Matrix d = matrix(F,nr,nc)
#        start_row = 0
#
#        for r from 0 <= r < nr:
#            for c from 0 <= c < nc:
#                p = self.get_unsafe(r,c)
#                if p.is_constant():
#                    d.set_unsafe(r, c, p.constant_coefficient())
#
#        for c from 0 <= c < nc:
#            r = -1
#            for rc from start_row <= rc < nr:
#                if d.get_unsafe(rc, c):
#                    r = rc
#                    break
#            if r!=-1:
#                a_inverse = ~self.get_unsafe(r,c)
#                self.rescale_row_c(r, a_inverse , c)
#                self.swap_rows_c(r, start_row)
#
#                for i from 0 <= i < nr:
#                    if i != start_row:
#                        minus_b = -self.get_unsafe(i,c)
#                        self.add_multiple_of_row(i, start_row, minus_b, 0)
#
#                start_row +=1
#
#                d = d._parent(0)
#                for i from start_row <= i < nr:
#                    for j from c+1 <= j < nc:
#                        if self.get_unsafe(i,j).is_constant():
#                            d.set_unsafe(i,j, self.get_unsafe(i,j).constant_coefficient())
#
#        self.cache('in_echelon_form_row_reduction',True)
#
#    def swapped_columns(self):
#        """
#        Return which columns were swapped during the Gauss-Bareiss reduction
#
#        OUTPUT:
#
#        Return a tuple representing the column swaps during the last application
#        of the Gauss-Bareiss algorithm (see :meth:`echelon_form` for details).
#
#        The tuple as length equal to the rank of self and the value at the
#        $i$-th position indicates the source column which was put as the $i$-th
#        column.
#
#        If no Gauss-Bareiss reduction was performed yet, None is
#        returned.
#
#        EXAMPLES::
#
#            sage: R.<x,y> = QQ[]
#            sage: C = random_matrix(R, 2, 2, terms=2)
#            sage: C.swapped_columns()
#            sage: E = C.echelon_form('bareiss')
#            sage: E.swapped_columns()
#            (0, 1)
#        """
#        return self.fetch('swapped_columns')



