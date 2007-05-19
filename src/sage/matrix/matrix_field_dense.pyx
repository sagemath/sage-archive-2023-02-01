"""nodoctest
Generic matrices over a field
"""

cimport matrix_pid_dense
import  matrix_pid_dense

import matrix_space
import sage.rings.polynomial.polynomial_ring
import sage.rings.number_field.number_field
import sage.misc.misc

from sage.rings.finite_field import is_FiniteField
from sage.rings.integer_mod_ring import is_IntegerModRing


cdef class Matrix_field_dense(matrix_pid_dense.Matrix_pid_dense):

    def __invert__(self):
        """
        Return this inverse of this matrix.

        Raises a ZeroDivisionError if the matrix has zero determinant, and
        raises an ArithmeticError, if the inverse doesn't exist
        because the matrix is nonsquare.

        EXAMPLES:

        """
        if not self.is_square():
            raise ArithmeticError, "self must be a square matrix"
        A = self.augment(self.parent().identity_matrix())
        B = A.echelon_form()
        if B[self.nrows()-1,self.ncols()-1] != 1:
            raise ZeroDivisionError, "self is not invertible"
        return B.matrix_from_columns(range(self.ncols(), 2*self.ncols()))

    def charpoly(self):
        """
        Return the characteristic polynomial of this matrix.

        ALGORITHM: Compute the Hessenberg form of the matrix and read
        off the characteristic polynomial from that.  The result is
        cached.

        EXAMPLES:
            sage: A = MatrixSpace(RationalField(),3)(range(9))
            sage: A.charpoly('x')
            x^3 - 12*x^2 - 18*x
            sage: A.trace()
            12
            sage: A.determinant()
            0
        """
        cdef int i

        try:
            return self._cache['charpoly']
        except KeyError:
            pass
        if self._nrows != self._ncols:
            raise ArithmeticError, "charpoly of non-square matrix not defined."

        R = sage.rings.polynomial.polynomial_ring.PolynomialRing(self.base_ring())
        zero = R(0)
        if self._nrows == 0:
            self._cache['charpoly'] = zero
            return zero
        time = sage.misc.misc.verbose(t=0)
        H = self.hessenberg_form()
        n = self.nrows()
        c = [zero]*(n+1)
        c[0] = R(1)
        X = R.gen()
        for m in range(1,n+1):
            c[m] = (X - R(H[m-1,m-1]))*c[m-1]
            t = 1
            for i in range(1,m):
                t = t*H[m-i, m-i-1]
                c[m] = c[m] - R(t*H[m-i-1,m-1])*c[m-i-1]
        sage.misc.misc.verbose('computed characteristic polynomial of %sx%s matrix'%
                     (self.nrows(), self.ncols()), time)
        f = c[n]
        self._cache['charpoly'] = f
        return f

    def column_space(self):
        """
        Return the vector space over the base ring spanned by the
        columns of this matrix.

        EXAMPLES:
            sage: M = MatrixSpace(QQ,3,3)
            sage: A = M([1,9,-7,4/5,4,3,6,4,3])
            sage: A.column_space()
            Vector space of degree 3 and dimension 3 over Rational Field
            Basis matrix:
            [1 0 0]
            [0 1 0]
            [0 0 1]
            sage: W = MatrixSpace(CC,2,2)
            sage: B = W([1, 2+3*I,4+5*I,9]); B
            [                       1.0000000000000000 2.0000000000000000 + 3.0000000000000000*I]
            [4.0000000000000000 + 5.0000000000000000*I                        9.0000000000000000]
            sage: B.column_space()
            Vector space of degree 2 and dimension 2 over Complex Field with 53 bits of precision
            Basis matrix:
            [                                                       1.0000000000000000 0.00000000000000044408920985006262 + 0.00000000000000088817841970012523*I]
            [                                                                        0               0.99999999999999978 - 0.000000000000000055511151231257827*I]
        """
        return self.column_module()

    def decomposition_of_subspace(self, M, is_diagonalizable=False):
        """
        Suppose the right action of self on M leaves M
        invariant. Return the decomposition of M as a list of pairs
        (W, is_irred) where is_irred is True if the charpoly of self
        acting on the factor W is irreducible.
        """
        if not self.is_square():
            raise ArithmeticError, "matrix must be square"
        if M.base_ring() != self.base_ring():
            raise ArithmeticError, "base rings are incompatible"
        if M.degree() != self.ncols():
            raise ArithmeticError, \
               "M must be a subspace of an %s-dimensional space"%self.ncols()

        time = sage.misc.misc.verbose(t=0)

        # 1. Restrict
        B = self.restrict(M)
        time0 = sage.misc.misc.verbose("restrict -- ", time)

        # 2. Decompose restriction
        D = B.decomposition(is_diagonalizable=is_diagonalizable, dual=False)

        cdef int dim
        dim = 0
        for A, _D in D:
            dim = dim + A.dimension()
        assert dim == M.dimension(), "bug in decomposition; " + \
               "the sum of the dimensions of the factors must equal the dimension of the space acted on"

        # 3. Lift decomposition to subspaces of ambient vector space.
        # Each basis vector for an element of D defines a linear combination
        # of the basis of W, and these linear combinations define the
        # corresponding subspaces of the ambient space M.

        sage.misc.misc.verbose("decomposition -- ", time0)
        C = M.basis_matrix()
        Z = M.ambient_vector_space()

        D = []
        for W, is_irred in D:
            tmp = []
            for x in W.basis():
                tmp.append(x*C)
            D.append((Z.subspace(tmp), is_irred))

        sage.misc.misc.verbose(t=time)
        return D

    def denominator(self):
        r"""
        Return the least common multiple of the denominators of the
        elements of self.

        If there is no denominator function for the base field, or no
        LCM function for the denominators, raise a TypeError.

        EXAMPLES:
            sage: A = MatrixSpace(RationalField(),2)(['1/2', '1/3', '1/5', '1/7'])
            sage: A.denominator()
            210

        Denominators are note defined for real numbers:
            sage: A = MatrixSpace(RealField(),2)([1,2,3,4])
            sage: A.denominator()
            Traceback (most recent call last):
            ...
            TypeError: denominator not defined for elements of the base ring

        We can even compute the denominator of matrix over the fraction field
        of $\Z[x]$.
            sage: K.<x> = FractionField(PolynomialRing(IntegerRing()))
            sage: A = MatrixSpace(K,2)([1/x, 2/(x+1), 1, 5/(x^3)])
            sage: A.denominator()
            x^4 + x^3

        Here's an example involving a cyclotomic field:
            sage: K.<z> = CyclotomicField(3)
            sage: M = MatrixSpace(K,3,sparse=True)
            sage: A = M([(1+z)/3,(2+z)/3,z/3,1,1+z,-2,1,5,-1+z])
            sage: print A
            [1/3*z + 1/3 1/3*z + 2/3       1/3*z]
            [          1       z + 1          -2]
            [          1           5       z - 1]
            sage: print A.denominator()
            3
        """
        if self.nrows() == 0 or self.ncols() == 0:
            return integer.Integer(1)
        R = self.base_ring()
        x = self.list()
        try:
            d = x[0].denominator()
        except AttributeError:
            raise TypeError, "denominator not defined for elements of the base ring"
        try:
            for y in x:
                d = d.lcm(y.denominator())
        except AttributeError:
            raise TypeError, "lcm function not defined for elements of the base ring"
        return d

    def echelon_form(self, include_zero_rows=True):
        """
        Returns the reduced row echelon form of self.

        INPUT:
            matrix -- an element A of a MatrixSpace

        OUTPUT:
            matrix -- The reduced row echelon form of A.
            Note that self is *not* changed by this command.

        EXAMPLES:
           sage: MS = MatrixSpace(RationalField(),2,3)
           sage: C = MS.matrix([1,2,3,4,5,6])
           sage: C.rank()
           2
           sage: C.nullity()
           1
           sage: C.echelon_form()
           [ 1  0 -1]
           [ 0  1  2]
        """
        try:
            return self.__echelon_form
        except AttributeError:
            pass

        R = self.base_ring()
        # Fix to work with finite fields and Z/nZ, which was
        # suggested by Dan Christensen <jdc@uwo.ca>.
        if (   (is_FiniteField(R) and R.is_prime_field()) or \
               is_IntegerModRing(R)  ) and R.characteristic() < 46340:
            p = R.characteristic()
            S = sage.ext.dense_matrix_pyx.Matrix_modint(p, self.nrows(), self.ncols(), self.list())
            S.echelon()
            A = self.parent()(S.list())    # most of time is spent here!?
            pivot_positions = S.pivots()

        else:

            t = sage.misc.misc.verbose("Generic echelon...")
            pivot_positions = []
            start_row = 0
            A = self.copy()
            nrows = A.nrows()
            ncols = A.ncols()
            cleared_a_column = False
            for c in range(ncols):
                sage.misc.misc.verbose("column %s of %s"%(c, ncols),t, level=2)
                for r in range(start_row, nrows):
                    if A.get((r,c)) != 0:
                        pivot_positions.append(c)
                        # Divide row r through by 1/A[r,c], so leading coefficient
                        # is 1.
                        z = A.get((r,c))
                        if z != 1:
                            A.rescale_row(r, ~z)
                        # Swap
                        if r != start_row:
                            A.swap_rows(r,start_row)
                        # Clear column
                        cleared_a_column = True
                        for i in range(nrows):
                            if i != start_row:
                                x = A.get((i,c))
                                if x != 0:
                                    # Add -x times start row to i
                                    A.add_multiple_of_row(i, start_row, -x)
                    if cleared_a_column:
                        start_row = start_row + 1
                        cleared_a_column = False
                        break
            # end for
            sage.misc.misc.verbose("Finished generic echelon.",t)
        #end if

        if not include_zero_rows:
            A = A.matrix_from_rows(range(len(pivot_positions)))
        A.__pivots = pivot_positions
        A.__rank = len(pivot_positions)
        A.set_immutable()
        return A

    def fcp(self, var='x'):
        """
        Return the factorization of the characteristic polynomial of self.

        EXAMPLES:
            sage: M = MatrixSpace(QQ,3,3)
            sage: A = M([1,9,-7,4/5,4,3,6,4,3])
            sage: A.fcp('x')
            (x^3 - 8*x^2 + 209/5*x - 286)
            sage: A = M([3, 0, -2, 0, -2, 0, 0, 0, 0])
            sage: A.fcp('x')
            (x - 3) * x * (x + 2)
        """
        return self.charpoly('x').factor()


    def hessenberg_form(self):
        """
        Return the Hessenberg form of self.

        The hessenberg form of a matrix $A$ is a matrix that is
        similar to $A$, so has the same characteristic polynomial as
        $A$, and is upper triangular except possible for entries right
        below the diagonal.

        ALGORITHM: See Henri Cohen's first book.

        EXAMPLES:
            sage: A = MatrixSpace(RationalField(),3)([2, 1, 1, -2, 2, 2, -1, -1, -1])
            sage: A.hessenberg_form()
            [  2 3/2   1]
            [ -2   3   2]
            [  0  -3  -2]

            sage: A = MatrixSpace(RationalField(),4)([2, 1, 1, -2, 2, 2, -1, -1, -1,1,2,3,4,5,6,7])
            sage: A.hessenberg_form()
            [    2  -7/2 -19/5    -2]
            [    2   1/2 -17/5    -1]
            [    0  25/4  15/2   5/2]
            [    0     0  58/5     3]
        """
        if not self.is_square():
            raise ArithmeticError, "self must be square"
        n = self.nrows()
        tm = sage.misc.misc.verbose("Computing Hessenberg Normal Form of %sx%s matrix"%(n,n))
        h = self.copy()
        for m in range(1,n-1):
            # Search for a non-zero entry in column m-1
            i = False
            for r in range(m+1,n):
                if h[r,m-1] != 0:
                    i = r
                    break
            if i:
                # Found a nonzero entry in column m-1 that is strictly below row m
                # Now set i to be the first nonzero position >= m in column m-1
                if h[m,m-1] != 0:
                    i = m
                t = h[i,m-1]
                if i>m:
                    h.swap_rows(i,m)
                    # We must do the corresponding column swap to
                    # maintain the characteristic polynomial (which is
                    # an invariant of Hessenberg form)
                    h.swap_columns(i,m)
                # Now the nonzero entry in position (m,m-1) is t.
                # Use t to clear the entries in column m-1 below m.
                for j in range(m+1,n):
                    if h[j,m-1] != 0:
                        u = h[j,m-1]/t
                        h.add_multiple_of_row(j, m, -u)
                        # To maintain charpoly, do the corresponding column operation,
                        # which doesn't mess up the matrix, since it only changes
                        # column m, and we're only worried about column m-1 right now.
                        # Add u*column_j to column_m.
                        h.add_multiple_of_column(m, j, u)
        sage.misc.misc.verbose("Finished Hessenberg Normal Form of %sx%s matrix"%(n,n),tm)
        return h

    def kernel(self, *args, **kwds):
        r"""
        Return the kernel of this matrix, as a vector space.

        INPUT:
            -- all additional arguments to the kernel function
               are passed directly onto the echelon call.

        \algorithm{Elementary row operations don't change the kernel,
        since they are just right multiplication by an invertible
        matrix, so we instead compute kernel of the column echelon
        form.  More precisely, there is a basis vector of the kernel
        that corresponds to each non-pivot row.  That vector has a 1
        at the non-pivot row, 0's at all other non-pivot rows, and for
        each pivot row, the negative of the entry at the non-pivot row
        in the column with that pivot element.}

        \note{Since we view matrices as acting on the right, but have
        functions for reduced \emph{row} echelon forms, we instead
        compute the reduced row echelon form of the transpose of this
        matrix, which is the reduced column echelon form.}

        EXAMPLES:

        A kernel of dimension one over $\Q$:x
            sage: A = MatrixSpace(QQ, 3)(range(9))
            sage: A.kernel()
            Vector space of degree 3 and dimension 1 over Rational Field
            Basis matrix:
            [ 1 -2  1]

        A trivial kernel:
            sage: A = MatrixSpace(QQ, 2)([1,2,3,4])
            sage: A.kernel()
            Vector space of degree 2 and dimension 0 over Rational Field
            Basis matrix:
            []

        Kernel of a zero matrix:
            sage: A = MatrixSpace(QQ, 2)(0)
            sage: A.kernel()
            Vector space of degree 2 and dimension 2 over Rational Field
            Basis matrix:
            [1 0]
            [0 1]

        Kernel of a non-square matrix:
            sage: A = MatrixSpace(QQ,3,2)(range(6))
            sage: A.kernel()
            Vector space of degree 3 and dimension 1 over Rational Field
            Basis matrix:
            [ 1 -2  1]

        The 2-dimensional kernel of a matrix over a cyclotomic field:
            sage: K = CyclotomicField(12); a=K.0
            sage: M = MatrixSpace(K,4,2)([1,-1, 0,-2, 0,-a**2-1, 0,a**2-1])
            sage: M
            [             1             -1]
            [             0             -2]
            [             0 -zeta12^2 - 1]
            [             0  zeta12^2 - 1]
            sage: M.kernel()
            Vector space of degree 4 and dimension 2 over Cyclotomic Field of order 12 and degree 4
            Basis matrix:
            [               0                1                0     -2*zeta12^2]
            [               0                0                1 -2*zeta12^2 + 1]

        A nontrivial kernel over a complicated base field.
            sage: K = FractionField(MPolynomialRing(QQ, 2))
            sage: M = MatrixSpace(K, 2)([[K.1, K.0], [K.1, K.0]])
            sage: M
            [x1 x0]
            [x1 x0]
            sage: M.kernel()
            Vector space of degree 2 and dimension 1 over Fraction Field of Polynomial Ring in x0, x1 over Rational Field
            Basis matrix:
            [ 1 -1]
        """

        R = self.base_ring()

        if self.nrows() == 0:    # from a 0 space
            V = sage.modules.free_module.VectorSpace(R, self.nrows())
            return V.zero_subspace()

        elif self.ncols() == 0:  # to a 0 space
            return sage.modules.free_module.VectorSpace(R, self.nrows())

        if isinstance(R, sage.rings.number_field.number_field.NumberField_generic):
            A = self._pari_().mattranspose()
            B = A.matker()
            n = self.nrows()
            V = sage.modules.free_module.VectorSpace(R, n)

            # This used to be:
            # basis = [V([R(x) for x in b]) for b in B]

            basis = []
            for b in B:
                z = []
                for x in b:
                    z.append(R(x))
                basis.append(V(z))

            return V.subspace(basis)

        E = self.transpose().echelon_form(*args, **kwds)
        pivots = E.pivots()
        pivots_set = set(pivots)
        basis = []
        VS = sage.modules.free_module.VectorSpace
        V = VS(R, self.nrows())
        ONE = R(1)
        for i in xrange(self.nrows()):
            if not (i in pivots_set):
                v = V(0)
                v[i] = ONE
                for r in range(len(pivots)):
                    v[pivots[r]] = -E[r,i]
                basis.append(v)
        return V.subspace(basis)


    def is_invertible(self):
        """
        Return True if this matrix is invertible.

        EXAMPLES:

            sage: M = MatrixSpace(QQ,3,3)
            sage: A = M([1,9,-7,4/5,4,3,6,4,3])
            sage: A.is_invertible()
            True
            sage: W = MatrixSpace(CC,2,2)
            sage: i = CC.0
            sage: B = W([1, 2+3*i, 4+5*i, 9])
            sage: B.is_invertible()
            True
            sage: N = MatrixSpace(QQ,2,2)
            sage: C = N([2,1,2,1])
            sage: C.is_invertible()
            False
        """
        return self.is_square() and self.rank() == self.nrows()

    def restrict(self, V, check=True):
        """
        Returns the matrix that defines the action of self on the
        chosen basis for the invariant subspace V.  If V is an
        ambient, returns self (not a copy of self).

        INPUT:
            V -- vector subspace
            check -- (optional) default: True; if False may not check
                     that V is invariant (hence can be faster).
        OUTPUT:
            a matrix

        WARNING:
        This function returns an nxn matrix, where V has dimension n.
        It does \emph{not} check that V is in fact invariant under
        self, unless check is True (not the default).

        EXAMPLES:
            sage: V = VectorSpace(QQ, 3)
            sage: M = MatrixSpace(QQ, 3)
            sage: A = M([1,2,0, 3,4,0, 0,0,0])
            sage: W = V.subspace([[1,0,0], [0,1,0]])
            sage: A.restrict(W)
            [1 2]
            [3 4]
            sage: A.restrict(W, check=True)
            [1 2]
            [3 4]

        We illustrate the warning about invariance not being checked
        by default, by giving a non-invariant subspace.  With the default
        check=False this function returns the 'restriction' matrix, which
        is meaningless as check=True reveals.
            sage: W2 = V.subspace([[1,0,0], [0,1,1]])
            sage: A.restrict(W2, check=False)
            [1 2]
            [3 4]
            sage: A.restrict(W2, check=True)
            Traceback (most recent call last):
            ...
            ArithmeticError: subspace is not invariant under matrix
        """
        if not isinstance(V, sage.modules.free_module.FreeModule_generic):
            raise TypeError, "V must be a Vector Space"
        if V.base_field() != self.base_ring():
            raise TypeError, "base rings must be the same"
        if V.degree() != self.nrows():
            raise IndexError, "degree of V (=%s) must equal number of rows of self (=%s)"%(V.degree(), self.nrows())
        if V.rank() == 0:
            return self.new_matrix(nrows=0, ncols=0)

        if not check and V.base_ring().is_field() and not V.has_user_basis():
            B = V.echelonized_basis_matrix()
            P = B.pivots()
            return B*self.matrix_from_columns(P)
        else:
            n = V.rank()
            try:
                # todo optimize so only involves matrix multiplies ?
                C = []
                for b in V.basis():
                    C.append(V.coordinate_vector(b*self))
                # Previous list comp version: C = [V.coordinate_vector(b*self) for b in V.basis()]
            except ArithmeticError:
                raise ArithmeticError, "subspace is not invariant under matrix"
            return self.new_matrix(n, n, C, sparse=False)

    def restrict_domain(self, V):
        """
        Compute the matrix relative to the basis for V on the domain
        obtained by restricting self to V, but not changing the
        codomain of the matrix.  This is the matrix whose rows are the
        images of the basis for V.

        INPUT:
            V -- vector space (subspace of ambient space on which self acts)

        SEE ALSO: restrict()

        EXAMPLES:
            sage: V = VectorSpace(QQ, 3)
            sage: M = MatrixSpace(QQ, 3)
            sage: A = M([1,2,0, 3,4,0, 0,0,0])
            sage: W = V.subspace([[1,0,0], [1,2,3]])
            sage: A.restrict_domain(W)
            [1 2 0]
            [3 4 0]
            sage: W2 = V.subspace_with_basis([[1,0,0], [1,2,3]])
            sage: A.restrict_domain(W2)
            [ 1  2  0]
            [ 7 10  0]
        """
        # Previous version: D = [b*self for b in V.basis()]
        D = []
        for b in V.basis():
            D.append(b*self)
        return self.new_matrix(V.dimension(), self.ncols(), D)


    def row_space(self):
        """
        Return the vector space over the base field spanned by the
        rows of self.

        EXAMPLES:
            sage: A = MatrixSpace(QQ, 2,2)([1,2,3,4])
            sage: A.row_space()
            Vector space of degree 2 and dimension 2 over Rational Field
            Basis matrix:
            [1 0]
            [0 1]
            sage: A.row_span(IntegerRing())
            Free module of degree 2 and rank 2 over Integer Ring
            Echelon basis matrix:
            [1 0]
            [0 2]
        """
        return self.row_module()

    def row_span(self, R=None):
        r"""
        Return the R-module spanned by the rows of self.

        INPUT:
            R -- (optional) principal ideal ring, such that the
                 entries of self coerce into R.  The default
                 is the base ring.
        OUTPUT:
            a free R-module.

        EXAMPLES:

        We define a $2\times 3$ matrix over $\Q$, then consider its row span
        over both $\Q$ and $\Z$.

            sage: A = MatrixSpace(QQ, 2,3)([1,2,3, '1/3', 4, 2])
            sage: A
            [  1   2   3]
            [1/3   4   2]
            sage: M1 = A.row_span(); M1
            Vector space of degree 3 and dimension 2 over Rational Field
            Basis matrix:
            [   1    0 12/5]
            [   0    1 3/10]
            sage: M2 = A.row_span(IntegerRing()); M2
            Free module of degree 3 and rank 2 over Integer Ring
            Echelon basis matrix:
            [1/3   4   2]
            [  0  10   3]

        Note that the determinants of the inner product matrices are
        different, though their discriminants differ by a square:

            sage: M1.inner_product_matrix()
            [ 169/25   18/25]
            [  18/25 109/100]
            sage: d1 = M1.discriminant(); d1
            137/20
            sage: d2 = M2.discriminant(); d2
            685/9
            sage: d2/d1
            100/9
        """
        if R is None:
            return self.row_module()
        M = sage.modules.free_module.FreeModule(R, self.ncols())
        return M.span(self.rows())

    def wiedemann(self, i, t=0):
        """
        Application of Wiedemann's algorithm to the i-th standard
        basis vector.

        If the optimal parameter t is nonzero, use only the first t
        linear recurrence relations.
        """
        i = int(i); t=int(t)
        if self.nrows() != self.ncols():
            raise ArithmeticError, "matrix must be square."
        n = self.nrows()
        v = sage.modules.free_module.VectorSpace(self.base_ring(), n).gen(i)
        tm = sage.misc.misc.verbose('computing iterates...')
        cols = self.iterates(v, 2*n).columns()
        tm = sage.misc.misc.verbose('computed iterates', tm)
        f = None
        # Compute the minimal polynomial of the linear recurrence
        # sequence corresponding to the 0-th entries of the iterates,
        # then the 1-th entries, etc.
        if t == 0:
            R = range(n)
        else:
            R = [t]
        for i in R:
            tm = sage.misc.misc.verbose('applying berlekamp-massey')
            g = berlekamp_massey.berlekamp_massey(cols[i].list())
            sage.misc.misc.verbose('berlekamp-massey done', tm)
            if f is None:
                f = g
            else:
                f = f.lcm(g)
            if f.degree() == n:
                break
        return f

