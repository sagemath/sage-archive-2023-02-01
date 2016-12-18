"""
Dense matrices over univariate polynomials over fields

This implementation inherits from Matrix_generic_dense, i.e., it is not
optimized for speed but only some methods were added.

AUTHORS:

- Kwankyu Lee (2016-12-15): Initial version with code moved from other files.
"""

#*****************************************************************************
#       Copyright (C) 2016 Kwankyu Lee <ekwankyu@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.matrix.matrix_generic_dense cimport Matrix_generic_dense

cdef class Matrix_polynomial_dense(Matrix_generic_dense):
    """
    Dense matrix over a univariate polynomial ring over a field.
    """

    def is_weak_popov(self):
        r"""
        Return ``True`` if this matrix is in weak Popov form.

        OUTPUT:

        A matrix over a polynomial ring is in weak Popov form if all
        leading positions are different [MS2003]_. A leading position
        is the position `i` in a row with the largest degree, for multiple
        entries with largest degrees, the maximal `i` is chosen
        (which is the furthest to the right in the matrix).

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

        But a matrix with less cols than nonzero rows is never in weak
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

            - :meth:`weak_popov_form <sage.matrix.matrix_polynomial_dense.weak_popov_form>`

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

        A matrix is in weak Popov form if the leading positions of the nonzero
        rows are all different. The leading position of a row is the right-most
        position whose entry has maximal degree of the entries in that row.

        .. WARNING::

            This function currently does **not** compute the weak Popov form of a
            matrix, but rather a row reduced form (which is a slightly weaker
            requirement). See :meth:`row_reduced_form`.

        INPUT:

        - ``transformation`` -- A boolean (default: `True`). If this is set to
          ``True``, the transformation matrix `U` will be returned as well: this
          is an invertible matrix over `k(x)` such that ``self`` equals `UW`,
          where `W` is the output matrix.

          Warning: the default of ``transformation`` will soon be set to ``False``,
          see :trac:`16896`. Until then, this function will print a deprecation
          warning as long as ``transformation`` was not explicitly set to ``True``
          or ``False``.

        - ``ascend`` -- Deprecated and has no effect.

        - ``old_call`` -- For backwards compatibility until the old calling
          convention will be deprecated (default: `True`). If ``True``, then
          return `(W,U,d)`, where `U` is as when ``transformation = True``, and
          `d` is a list of the degrees of the rows of `W`.
        """
        from sage.misc.superseded import deprecation
        deprecation(16888, "This function currently does *not* compute a weak Popov form, "
        "but rather a row reduced form. This function will soon be fixed (see Ticket #16742).")

        return self.row_reduced_form(transformation=transformation,
                ascend=ascend, old_call=old_call)

    def _weak_popov_form(self, transformation=False, shifts=None):
        """
        Return a weak Popov form of this matrix.

        A matrix is in weak Popov form if the leading positions of the nonzero
        rows are all different. The leading position of a row is the right-most
        position whose entry has the maximal degree in the row.

        .. WARNING::

            This method will replace the current :meth:`weak_popov_form` when
            the current one is removed after deprecation period.

        INPUT:

        - ``transformation`` -- If this is ``True``, the transformation matrix
          is returned together with the weak Popov form.

        - ``shifts`` -- If this is a tuple of integers of length the number of
          columns of this matrix, then each integer of the tuple is added to
          the degree of the polynomial in the corresponding column to determine
          the leading position of the rows.

        ALGORITHM:

        This method implements the Mulders-Storjohann algorithm of [MS2003]_.

        EXAMPLES::

            sage: F.<a> = GF(2^4,'a')
            sage: PF.<x> = F[]
            sage: A = matrix(PF,[[1,a*x^17+1],[0,a*x^11+a^2*x^7+1]])
            sage: M, U = A._weak_popov_form(transformation=True)
            sage: U * A == M
            True
            sage: M.is_weak_popov()
            True
            sage: U.is_invertible()
            True

        A zero matrix will return itself::

            sage: Z = matrix(PF,5,3)
            sage: Z._weak_popov_form()
            [0 0 0]
            [0 0 0]
            [0 0 0]
            [0 0 0]
            [0 0 0]

        Matrices in weak popov form will just be returned untouched::

            sage: F.<a> = GF(17,'a')
            sage: PF.<x> = F[]
            sage: C = matrix(PF,[[1,7,x],[x^2,x,4],[2,x,11]]); C
            [  1   7   x]
            [x^2   x   4]
            [  2   x  11]
            sage: D = C._weak_popov_form()
            sage: D == C
            True

        And the transformation matrix will be the identity matrix::

            sage: C._weak_popov_form(transformation=True)
            (
            [  1   7   x]  [1 0 0]
            [x^2   x   4]  [0 1 0]
            [  2   x  11], [0 0 1]
            )

        .. SEEALSO::

            :meth:`is_weak_popov <sage.matrix.matrix_polynomial_dense.is_weak_popov>`
        """
        M = self.__copy__()
        x = M.base_ring().gen()

        if transformation:
            from sage.matrix.constructor import identity_matrix
            U = identity_matrix(M.base_ring(), M.nrows())

        def simple_transformation(i, j, pos):
            """
            Perform a simple transformation with row j on row i, reducing
            position pos.

            If M[i,pos] has lower degree than M[j,pos], nothing is changed.
            """
            pow = M[i,pos].degree() - M[j,pos].degree()
            if pow < 0:
                return
            coeff = -M[i, pos].leading_coefficient() / M[j, pos].leading_coefficient()

            multiple = x**pow * coeff
            M.add_multiple_of_row(i, j, multiple)
            if transformation:
                U.add_multiple_of_row(i, j, multiple)

        if shifts and len(shifts) != M.ncols():
            raise ValueError("The number of shifts must equal the number of columns")

        if shifts:
            def LP(v):
                best = None
                bestp = None
                for p in range(0,len(v)):
                    if not v[p].is_zero():
                        vpdeg = v[p].degree() + shifts[p]
                        if vpdeg >= best:
                            best = vpdeg
                            bestp = p
                if best == None:
                    return -1
                else:
                    return bestp
        else:
            def LP(v):
                bestp = 0
                best = v[0].degree()
                for p in range(1,len(v)):
                    vpdeg = v[p].degree()
                    if vpdeg >= best:
                        best = vpdeg
                        bestp = p
                if best < 0:
                    return -1
                else:
                    return bestp

        # initialise conflicts list and LP map
        LP_to_row = dict( (i,[]) for i in range(M.ncols()))
        conflicts = []
        for i in range(M.nrows()):
            lp = LP(M.row(i))
            if lp >= 0:
                ls = LP_to_row[lp]
                ls.append(i)
                if len(ls) > 1:
                    conflicts.append(lp)

        # while there is a conflict, do a simple transformation
        while conflicts:
            lp = conflicts.pop()
            ls = LP_to_row[lp]
            i, j = ls.pop(), ls.pop()
            if M[i,lp].degree() < M[j, lp].degree():
                j,i = i,j

            simple_transformation(i, j, lp)
            ls.append(j)
            lp_new = LP(M.row(i))
            if lp_new > -1:
                ls_new = LP_to_row[lp_new]
                ls_new.append(i)
                if len(ls_new) > 1:
                    conflicts.append(lp_new)

        if transformation:
            return M, U
        else:
            return M

    def row_reduced_form(self, transformation=None, ascend=None, old_call=True):
        r"""
        Return a row reduced form of this matrix.

        A matrix `M` is row reduced if the leading term matrix has the same rank
        as `M`. The leading term matrix of a polynomial matrix `M_0` is the matrix
        over `k` whose `(i,j)`'th entry is the `x^{d_i}` coefficient of `M_0[i,j]`,
        where `d_i` is the greatest degree among polynomials in the `i`'th row of `M_0`.

        INPUT:

        - ``transformation`` -- (default: ``False``). If this is set to
          ``True``, the transformation matrix `U` will be returned as well: this
          is an invertible matrix over `k(x)` such that ``self`` equals `UW`,
          where `W` is the output matrix.

        - ``ascend`` -- Deprecated and has no effect.

        - ``old_call`` -- For backwards compatibility until the old calling
          convention will be deprecated (default: ``True``). If ``True``, then
          return `(W,U,d)`, where `U` is as when ``transformation = True``, and
          `d` is a list of the degrees of the rows of `W`.

        OUTPUT:

        - `W` -- row reduced form of this matrix.

        EXAMPLES::

            sage: R.<t> = GF(3)['t']
            sage: K = FractionField(R)
            sage: M = matrix([[(t-1)^2],[(t-1)]])
            sage: M.row_reduced_form(transformation=False, old_call=False)
            [    0]
            [t + 2]

        If the matrix is an `n \times 1` matrix with at least one non-zero entry,
        `W` has a single non-zero entry and that entry is a scalar multiple of
        the greatest-common-divisor of the entries of the matrix.

        ::

            sage: M1 = matrix([[t*(t-1)*(t+1)],[t*(t-2)*(t+2)],[t]])
            sage: output1 = M1.row_reduced_form(transformation=False, old_call=False)
            sage: output1
            [0]
            [0]
            [t]

        The following is the first half of example 5 in [Hes2002]_ *except* that we
        have transposed the matrix; [Hes2002]_ uses column operations and we use row.

        ::

            sage: R.<t> = QQ['t']
            sage: M = matrix([[t^3 - t,t^2 - 2],[0,t]]).transpose()
            sage: M.row_reduced_form(transformation=False, old_call=False)
            [      t    -t^2]
            [t^2 - 2       t]

        The next example demonstrates what happens when the matrix is a zero matrix.

        ::

            sage: R.<t> = GF(5)['t']
            sage: M = matrix(R, 2, [0,0,0,0])
            sage: M.row_reduced_form(transformation=False, old_call=False)
            [0 0]
            [0 0]

        In the following example, the original matrix is already row reduced, but
        the output is a different matrix. This is because currently this method
        simply computes a weak Popov form, which is always also a row reduced matrix
        (see :meth:`weak_popov_form`). This behavior is likely to change when a faster
        algorithm designed specifically for row reduced form is implemented in Sage.

        ::

            sage: R.<t> = QQ['t']
            sage: M = matrix([[t,t,t],[0,0,t]]); M
            [t t t]
            [0 0 t]
            sage: M.row_reduced_form(transformation=False, old_call=False)
            [ t  t  t]
            [-t -t  0]

        The last example shows the usage of the transformation parameter.

        ::

            sage: Fq.<a> = GF(2^3)
            sage: Fx.<x> = Fq[]
            sage: A = matrix(Fx,[[x^2+a,x^4+a],[x^3,a*x^4]])
            sage: W,U = A.row_reduced_form(transformation=True,old_call=False);
            sage: W,U
            (
            [          x^2 + a           x^4 + a]  [1 0]
            [x^3 + a*x^2 + a^2               a^2], [a 1]
            )
            sage: U*W == A
            True
            sage: U.is_invertible()
            True

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

        - ``transformation`` -- (default: ``False``) If this is ``True``,
           the transformation matrix is output.

        OUTPUT:

        If ``transformation`` is ``True``, this function will output matrices ``W`` and ``N`` such that

        1. `W` -- a row reduced form of this matrix `M`.
        2. `N` -- a unimodular matrix satisfying `N * W = M`.

        If ``transformation`` is ``False``, the output is just `W`.

        EXAMPLES::

            sage: Fq.<a> = GF(2^3)
            sage: Fx.<x> = Fq[]
            sage: A = matrix(Fx,[[x^2+a,x^4+a],[x^3,a*x^4]])
            sage: A._row_reduced_form(transformation=True)
            (
            [          x^2 + a           x^4 + a]  [1 0]
            [x^3 + a*x^2 + a^2               a^2], [a 1]
            )
        """
        return self._weak_popov_form(transformation=transformation)
