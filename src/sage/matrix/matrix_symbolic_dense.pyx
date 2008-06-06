"""
Symbolic matrices

Matrices with symbolic entries.  The underlying representation is a
pointer to a Maxima object.
"""


from sage.structure.element cimport ModuleElement, RingElement
from sage.structure.factorization import Factorization

cimport matrix_dense
cimport matrix

cdef maxima
from sage.calculus.calculus import maxima

from sage.calculus.calculus import symbolic_expression_from_maxima_string, SymbolicVariable, var_cmp

cdef class Matrix_symbolic_dense(matrix_dense.Matrix_dense):
    r"""
    The \class{Matrix_generic_dense} class derives from \class{Matrix}, and
    defines functionality for dense matrices over any base ring.
    Matrices are represented by a list of elements in the base ring,
    and element access operations are implemented in this class.
    """
    ########################################################################
    # LEVEL 1 functionality
    #   * __new__   (not needed)
    # x * __init__
    #   * __dealloc__   (not needed)
    # x * set_unsafe
    # x * get_unsafe
    # x * def _pickle
    # x * def _unpickle
    ########################################################################
    def __init__(self, parent, entries, copy, bint coerce):
        """
        Create a symbolic matrix, which is internally represented
        by a maxima matrix.

        EXAMPES:
            sage: matrix(SR, 2, 2, range(4))
            [0 1]
            [2 3]
            sage: matrix(SR, 2, 2, var('t'))
            [t 0]
            [0 t]
        """
        matrix.Matrix.__init__(self, parent)
        self.__variables = None
        self._simp = None

        cdef Py_ssize_t i, n

        if entries is None:
            entries = 0

        if not isinstance(entries, list):
            try:
                x = parent.base_ring()(entries)
                is_list = 0
                coerce = 0
            except TypeError:
                try:
                    entries = list(entries)
                    is_list = 1
                except TypeError:
                    raise TypeError, "entries must be coercible to a list or the basering"

        else:
            is_list = 1

        if is_list:
            n = self._nrows * self._ncols
            if len(entries) != n:
                raise TypeError, "entries has the wrong length"

            if coerce:
                R = parent.base_ring()
                for i from 0 <= i < n:
                    entries[i] = R(entries[i])
            self.set_from_list(entries)

        else:
            if x != parent.base_ring()(0):
                if self._nrows != self._ncols:
                    raise TypeError, "nonzero scalar matrix must be square"
                self._maxima = maxima(x) * maxima.ident(self._nrows)
            else:
                self._maxima = maxima.zeromatrix(self._nrows, self._ncols)

    def _new_c(self):
        """
        Called when creating a new matrix.

        EXAMPLES:
            sage: matrix(SR,0)   # this implicitly calls _new_c
            []
        """
        cdef Matrix_symbolic_dense M = Matrix_symbolic_dense.__new__(Matrix_symbolic_dense, 0, 0, 0)
        matrix.Matrix.__init__(M, self._parent)
        return M

    cdef set_from_list(self, entries):
        """
        Set a matrix from a list of entries, all known to be in the symbolic ring.
        """
        s = ','.join(['[' + ','.join([z._maxima_init_() for z in entries[i:i+self._ncols]]) + ']' for
                      i from 0 <= i < self._ncols*self._nrows by self._ncols])
        self._maxima = maxima('matrix(%s)'%s)

    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, value):
        """
        Set the i,j entry without bounds checking.
        """
        maxima.setelmx(value, i+1, j+1, self._maxima)

    cdef get_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        """
        Get the i,j entry without bounds checking.
        """
        entry = symbolic_expression_from_maxima_string(maxima.eval("%s[%s,%s]" % (self._maxima._name, i+1, j+1)))
        return self._parent.base_ring()(entry)

    def _pickle(self):
        """
        Used when serializing a symbolic matrix.

        EXAMPLES:
            sage: M = matrix(SR, 2, 2, range(4))
            sage: type(M)
            <type 'sage.matrix.matrix_symbolic_dense.Matrix_symbolic_dense'>
            sage: loads(dumps(M)) == M         # implicitly calls _pickle
            True
        """
        return self._list(), 0

    def _unpickle(self, data, int version):
        """
        Used when de-serializing a symbolic matrix.

        EXAMPLES:
            sage: m = matrix(SR, 2, [sqrt(2), 3, pi, e]); m
            [sqrt(2)       3]
            [     pi       e]
            sage: loads(dumps(m))   # implicitly calls _unpickle
            [sqrt(2)       3]
            [     pi       e]
        """
        if version == 0:
            self.set_from_list(data)
        else:
            raise RuntimeError, "unknown matrix version"

    def __richcmp__(matrix.Matrix self, right, int op):  # always need for mysterious reasons.
        """
        Called implicitly when doing comparisons.

        EXAMPLES:
            sage: m = matrix(SR, 2, [sqrt(2), 3, pi, e])
            sage: cmp(m,m)
            0
            sage: cmp(m,3)
            -1
        """
        return self._richcmp(right, op)

    def __hash__(self):
        """
        Return hash of this matrix.

        NOTE: The hash of a symbolic matrix agrees with its hash in Maxima.

        EXAMPLES:
            sage: m = matrix(SR, 2, [sqrt(2), 3, pi, e]); m
            [sqrt(2)       3]
            [     pi       e]
            sage: hash(m) #random due to architecture dependence
            1532248127            # 32-bit
            1653238849131003967   # 64-bit
            sage: m.__hash__() #random due to architecture dependence
            1532248127            # 32-bit
            1653238849131003967   # 64-bit
            sage: hash(maxima(m)) #random due to architecture dependence
            1532248127            # 32-bit
            1653238849131003967   # 64-bit
        """
        return hash(self._maxima)

    ########################################################################
    # LEVEL 2 functionality
    # x * cdef _add_c_impl
    #   * cdef _mul_c_impl
    #   * cdef _cmp_c_impl
    # x * __neg__
    # x * __invert__
    # x * __copy__
    # x * _multiply_classical
    # x * _list -- copy of the list of underlying elements
    #   * _dict -- copy of the sparse dictionary of underlying elements
    ########################################################################

    cdef ModuleElement _add_c_impl(self, ModuleElement right):
        cdef Matrix_symbolic_dense M = self._new_c()
        M._maxima = self._maxima + (<Matrix_symbolic_dense>right)._maxima
        return M

    def __neg__(self):
        """
        Return the negative of this matrix.

        EXAMPLES:
            sage: -matrix(SR, 2, range(4))
            [ 0 -1]
            [-2 -3]
        """
        cdef Matrix_symbolic_dense M = self._new_c()
        M._maxima = -self._maxima
        return M

    def __invert__(self):
        """
        Return the inverse of this matrix.

        EXAMPLES:
            sage: M = matrix(SR, 2, var('a,b,c,d'))
            sage: ~M
            [ d/(a*d - b*c) -b/(a*d - b*c)]
            [-c/(a*d - b*c)  a/(a*d - b*c)]
            sage: M = matrix(SR, 3, 3, range(9)) - var('t')
            sage: (~M * M).simplify_rational()
            [1 0 0]
            [0 1 0]
            [0 0 1]
        """
        cdef Matrix_symbolic_dense M = self._new_c()
        M._maxima = self._maxima.invert()
        return M

    def __copy__(self):
        """
        Make a copy of this matrix.

        EXAMPLES:
        We make a copy of a matrix and change one entry:
            sage: m = matrix(SR, 2, [sqrt(2), 3, pi, e])
            sage: n = copy(m)
            sage: n[0,0] = sin(1)
            sage: m
            [sqrt(2)       3]
            [     pi       e]
            sage: n
            [sin(1)      3]
            [    pi      e]
        """
        cdef Matrix_symbolic_dense M = self._new_c()
        M._maxima = self._maxima.copymatrix()
        if self.subdivisions is not None:
            M.subdivide(*self.get_subdivisions())
        return M

    def _multiply_classical(left, Matrix_symbolic_dense right):
        """
        Multiply this matrix by another matrix.

        EXAMPLES:
            sage: m = matrix(SR, 4, [1..4^2])
            sage: m * m              # _multiply_classical is called implicitly
            [ 90 100 110 120]
            [202 228 254 280]
            [314 356 398 440]
            [426 484 542 600]

            sage: m = matrix(SR, 3, [1, 2, 3]); m
            [1]
            [2]
            [3]
            sage: m.transpose() * m
            [14]
        """
        if left._ncols != right._nrows:
            raise IndexError, "Number of columns of left must equal number of rows of other."
        cdef Matrix_symbolic_dense M = Matrix_symbolic_dense.__new__(Matrix_symbolic_dense, 0, 0, 0)
        P = left.matrix_space(left._nrows, right._ncols)
        matrix.Matrix.__init__(M, P)

        res = left._maxima.dot(right._maxima)

        #Since maxima doesn't return a 1x1 matrix when multiplying
        #a 1xn and a nx1 matrix, we need to wrap it in a matrix.
        if left._nrows == 1 and right._ncols == 1:
            res = res.parent().matrix([res])

        M._maxima = res
        return M

    def _list(self):
        r"""
        Return cached underlying list of entries of self.

        It is dangerous to use this if you change it! -- instead use
        \code{self.list()} which returns a copy.

        EXAMPLES:
            sage: m = matrix(SR, 2, [sqrt(2), 3, pi, e]); m
            [sqrt(2)       3]
            [     pi       e]
            sage: v = m._list(); v
            [sqrt(2), 3, pi, e]
            sage: w = m.list(); w
            [sqrt(2), 3, pi, e]

        We now illustrate the caching behavior of \code{list} versus \code{_list}.
            sage: v[0] = 10
            sage: m._list()
            [10, 3, pi, e]
            sage: m.list()
            [10, 3, pi, e]
            sage: w[0] = 100
            sage: m.list()
            [10, 3, pi, e]
            sage: m._list()
            [10, 3, pi, e]
        """
        l = self.fetch('_list')
        if not l is None: return l
        v = [self.get_unsafe(i,j) for i from 0 <= i < self._nrows for j from 0 <= j < self._ncols]
        self.cache('_list', v)
        return v

    cdef ModuleElement _lmul_c_impl(self, RingElement right):
        """
        Multiply a scalar times a matrix efficiently.

        EXAMPLES:
            sage: m = matrix(SR, 2, [1..4]); sqrt(2)*m
            [  sqrt(2) 2*sqrt(2)]
            [3*sqrt(2) 4*sqrt(2)]
        """
        cdef Matrix_symbolic_dense M = self._new_c()
        M._maxima = right._maxima_() * self._maxima
        return M


    ########################################################################
    # LEVEL 3 functionality (Optional)
    #  x * cdef _sub_c_impl
    #    * __deepcopy__
    #    * __invert__
    #  x * _multiply_classical
    #    * Matrix windows -- only if you need strassen for that base
    #    * Other functions (list them here):
    ########################################################################

    cdef ModuleElement _sub_c_impl(self, ModuleElement right):
        cdef Matrix_symbolic_dense M = self._new_c()
        M._maxima = self._maxima - (<Matrix_symbolic_dense>right)._maxima
        return M

    def _maxima_(self, session=None):
        """
        Return maxima version of this matrix.

        INPUT:
            session -- Maxima session.  If None, returns the session sage.calculus.calculus.maxima.

        OUTPUT:
            a maxima object with parent session

        EXAMPLES:
            sage: m = matrix(SR, 2, [sqrt(2), 3, pi, e])
            sage: m._maxima_()
            matrix([sqrt(2),3],[%pi,%e])

        Now put the same matrix in a completely different Maxima sesssion:
            sage: m._maxima_(maxima)
            matrix([sqrt(2),3],[%pi,%e])
        """
        if session is None:
            return self._maxima
        else:
            return session(self._maxima)

    def determinant(self):
        """
        Returns the determinant of self, as calculated by
        maxima using a Gaussian elimination-like algorithm.

        EXAMPLES:
            sage: M = matrix(SR, 2, 2, [x,2,3,4])
            sage: M.determinant()
            4*x - 6
            sage: M = matrix(SR, 3,3,range(9))
            sage: M.det()
            0
            sage: t = var('t')
            sage: M = matrix(SR, 2, 2, [cos(t), sin(t), -sin(t), cos(t)])
            sage: M.det()
            sin(t)^2 + cos(t)^2
        """
        det = repr(self._maxima.determinant())
        return self._parent.base_ring()(symbolic_expression_from_maxima_string(det))

    def permanent(self):
        """
        Returns the permanent of self.

        EXAMPLES:
            sage: M = matrix(SR, 2, 2, [x,2,3,4])
            sage: M.permanent()
            4*x + 6
        """
        perm = repr(self._maxima.permanent())
        return self._parent.base_ring()(symbolic_expression_from_maxima_string(perm))

    def rank(self):
        """
        EXAMPLES:
            sage: M = matrix(SR, 5, 5, range(25))
            sage: M.rank()
            2
            sage: M = matrix(SR, 5, 5, range(25)) - var('t')
            sage: M.rank()
            5

        WARNING: rank may return the wrong answer if it cannot determine that
                 a matrix element that is equivalent to zero is indeed so.
        """
        rank = repr(self._maxima.rank())
        return self._parent.base_ring()(symbolic_expression_from_maxima_string(rank))

    def charpoly(self, var='x'):
        """
        Compute the characteristic polynomial of self, using maxima.

        EXAMPLES:
            sage: M = matrix(SR, 2, 2, var('a,b,c,d'))
            sage: M.charpoly('t')
            (a - t)*(d - t) - b*c
            sage: matrix(SR, 5, [1..5^2]).charpoly().expand()
            x^5 - 65*x^4 - 250*x^3

        """
        # Maxima has the definition det(matrix-xI) instead of
        # det(xI-matrix), which is what Sage uses elsewhere.  We
        # correct for the discrepancy.
        if self.nrows() % 2 == 0:
            res = repr(self._maxima.charpoly(var))
        else:
            res = repr(-self._maxima.charpoly(var))

        return self._parent.base_ring()(symbolic_expression_from_maxima_string(res))

#    def echelonize(self, **kwds):
#        self._maxima = self._maxima.echelon()

    def eigenvalues(self, solution_set=False):
        """
        Compute the eigenvalues by solving the characteristic
        polynomial in maxima

        EXAMPLE:
            sage: a=matrix(SR,[[1,2],[3,4]])
            sage: a.eigenvalues()
            [(5 - sqrt(33))/2, (sqrt(33) + 5)/2]

        """
        tmp = SymbolicVariable('tmp_var')
        sols = self.charpoly(tmp).expand().solve(tmp)
        if solution_set:
            return sols
        else:
            values = []
            for sol in sols:
                if sol.left() != tmp or tmp in sol.right().variables():
                    raise ValueError, "Unable to symbolically extract eigenvalues, use solution_set=True option."
                values.append(sol.right())
            return values

    def exp(self):
        r"""
        Return the matrix exponential of this matrix $X$, which is the matrix
        $$
             e^X = \sum_{k=0}^{\infty} \frac{X^k}{k!}.
        $$

        This function depends on maxima's matrix exponentiation
        function, which does not deal well with floating point
        numbers.  If the matrix has floating point numbers, they will
        be rounded automatically to rational numbers during the
        computation.

        EXAMPLES:
            sage: m = matrix(SR,2, [0,x,x,0]); m
            [0 x]
            [x 0]
            sage: m.exp()
            [e^(-x)*(e^(2*x) + 1)/2 e^(-x)*(e^(2*x) - 1)/2]
            [e^(-x)*(e^(2*x) - 1)/2 e^(-x)*(e^(2*x) + 1)/2]
            sage: exp(m)
            [e^(-x)*(e^(2*x) + 1)/2 e^(-x)*(e^(2*x) - 1)/2]
            [e^(-x)*(e^(2*x) - 1)/2 e^(-x)*(e^(2*x) + 1)/2]
            sage: m.exp().expand()
            [e^x/2 + e^(-x)/2 e^x/2 - e^(-x)/2]
            [e^x/2 - e^(-x)/2 e^x/2 + e^(-x)/2]

        Exp works on 0x0 and 1x1 matrices:
            sage: m = matrix(SR,0,[]); m
            []
            sage: m.exp()
            []
            sage: m = matrix(SR,1,[2]); m
            [2]
            sage: m.exp()
            [e^2]

        Commuting matrices $m, n$ have the property that $e^{m+n} =
        e^m e^n$ (but non-commuting matrices need not):
            sage: m = matrix(SR,2,[1..4]); n = m^2
            sage: a = exp(m+n) - exp(m)*exp(n)
            sage: bool(a == 0)
            True

        The input matrix must be square:
            sage: m = matrix(SR,2,3,[1..6]); exp(m)
            Traceback (most recent call last):
            ...
            ValueError: exp only defined on square matrices

        In this example we take the symbolic answer and make it numerical at the end:
            sage: exp(matrix(SR, [[1.2, 5.6], [3,4]])).change_ring(RDF)
            [346.557487298 661.734590934]
            [354.500673715 677.424782765]

        Another example involving the reversed identity matrix, which we clumsily create.
            sage: m = identity_matrix(SR,4); m = matrix(list(reversed(m.rows()))) * x
            sage: exp(m)
            [e^(-x)*(e^(2*x) + 1)/2                      0                      0 e^(-x)*(e^(2*x) - 1)/2]
            [                     0 e^(-x)*(e^(2*x) + 1)/2 e^(-x)*(e^(2*x) - 1)/2                      0]
            [                     0 e^(-x)*(e^(2*x) - 1)/2 e^(-x)*(e^(2*x) + 1)/2                      0]
            [e^(-x)*(e^(2*x) - 1)/2                      0                      0 e^(-x)*(e^(2*x) + 1)/2]

        """
        if not self.is_square():
            raise ValueError, "exp only defined on square matrices"
        if self.nrows() == 0:
            return self
        cdef Matrix_symbolic_dense M = self._new_c()
        # Maxima's matrixexp function chokes on floating point numbers
        # so we automatically convert floats to rationals by passing
        # keepfloat: false
        z = maxima('matrixexp(%s), keepfloat: false'%self._maxima.name())
        if self.nrows() == 1:
            # We do the following, because Maxima stupidly exp's 1x1 matrices into non-matrices!
            M._maxima = maxima('matrix([%s])'%z.name())
        else:
            M._maxima = z
        return M


    #############################################
    # Simplification commands
    #############################################
    def is_simplified(self):
        """
        Return True if self is the result of running simplify() on a symbolic
        matrix.  This has the semantics of 'has_been_simplified'.

        EXAMPLES:
            sage: var('x,y,z')
            (x, y, z)
            sage: m = matrix([[z, (x+y)/(x+y)], [x^2, y^2+2]]); m
            [      z       1]
            [    x^2 y^2 + 2]
            sage: m.is_simplified()
            False
            sage: ms = m.simplify(); ms
            [      z       1]
            [    x^2 y^2 + 2]

            sage: m.is_simplified()
            False
            sage: ms.is_simplified()
            True
        """

        return self._simp is self

    def simplify(self):
        """
        Simplifies self.

        EXAMPLES:
            sage: var('x,y,z')
            (x, y, z)
            sage: m = matrix([[z, (x+y)/(x+y)], [x^2, y^2+2]]); m
            [      z       1]
            [    x^2 y^2 + 2]
            sage: m.simplify()
            [      z       1]
            [    x^2 y^2 + 2]
        """
        if self._simp is not None:
            return self._simp

        S = self.parent()([x.simplify() for x in self.list()])
        S._simp = S
        self._simp = S
        return S


    def simplify_trig(self):
        """
        EXAMPLES:
            sage: theta = var('theta')
            sage: M = matrix(SR, 2, 2, [cos(theta), sin(theta), -sin(theta), cos(theta)])
            sage: ~M
            [ cos(theta)/(sin(theta)^2 + cos(theta)^2) -sin(theta)/(sin(theta)^2 + cos(theta)^2)]
            [ sin(theta)/(sin(theta)^2 + cos(theta)^2)  cos(theta)/(sin(theta)^2 + cos(theta)^2)]
            sage: (~M).simplify_trig()
            [ cos(theta) -sin(theta)]
            [ sin(theta)  cos(theta)]
        """
        cdef Matrix_symbolic_dense M = self._new_c()
        M._maxima = self._maxima.trigexpand().trigsimp()
        return M

    def simplify_rational(self):
        """
        EXAMPLES:
            sage: M = matrix(SR, 4, 4, range(16)) - var('t')
            sage: (~M*M)[0,0]
            (-((10 - t)*(15 - t) - 154)*(5 - t) + 6*(9*(15 - t) - 143) - 7*(126 - 13*(10 - t)))*t/((-((10 - t)*(15 - t) - 154)*(5 - t) + 6*(9*(15 - t) - 143) - 7*(126 - 13*(10 - t)))*t - 4*((10 - t)*(15 - t) - 154) + 6*(8*(15 - t) - 132) - 7*(112 - 12*(10 - t)) + 2*((132 - 8*(15 - t))*(5 - t) + 4*(9*(15 - t) - 143) - 28) + 3*((112 - 12*(10 - t))*(5 - t) - 4*(126 - 13*(10 - t)) + 24)) + 4*((15 - t)*(t - 10) + 2*(9*(15 - t) - 143) - 3*(126 - 13*(10 - t)) + 154)/((-((10 - t)*(15 - t) - 154)*(5 - t) + 6*(9*(15 - t) - 143) - 7*(126 - 13*(10 - t)))*t - 4*((10 - t)*(15 - t) - 154) + 6*(8*(15 - t) - 132) - 7*(112 - 12*(10 - t)) + 2*((132 - 8*(15 - t))*(5 - t) + 4*(9*(15 - t) - 143) - 28) + 3*((112 - 12*(10 - t))*(5 - t) - 4*(126 - 13*(10 - t)) + 24)) + 8*(6*(15 - t) - 2*((5 - t)*(15 - t) - 91) + 3*(14*(5 - t) - 78) - 98)/((-((10 - t)*(15 - t) - 154)*(5 - t) + 6*(9*(15 - t) - 143) - 7*(126 - 13*(10 - t)))*t - 4*((10 - t)*(15 - t) - 154) + 6*(8*(15 - t) - 132) - 7*(112 - 12*(10 - t)) + 2*((132 - 8*(15 - t))*(5 - t) + 4*(9*(15 - t) - 143) - 28) + 3*((112 - 12*(10 - t))*(5 - t) - 4*(126 - 13*(10 - t)) + 24)) + 12*(7*(10 - t) - 3*((5 - t)*(10 - t) - 54) + 2*(11*(5 - t) - 63) - 66)/((-((10 - t)*(15 - t) - 154)*(5 - t) + 6*(9*(15 - t) - 143) - 7*(126 - 13*(10 - t)))*t - 4*((10 - t)*(15 - t) - 154) + 6*(8*(15 - t) - 132) - 7*(112 - 12*(10 - t)) + 2*((132 - 8*(15 - t))*(5 - t) + 4*(9*(15 - t) - 143) - 28) + 3*((112 - 12*(10 - t))*(5 - t) - 4*(126 - 13*(10 - t)) + 24))
            sage: (~M * M).simplify_rational()
            [1 0 0 0]
            [0 1 0 0]
            [0 0 1 0]
            [0 0 0 1]
        """
        cdef Matrix_symbolic_dense M = self._new_c()
        M._maxima = self._maxima.fullratsimp()
        return M

    def factor(self):
        """
        Operates pointwise on each element.

        EXAMPLES:
            sage: M = matrix(SR, 2, 2, x^2 - 2*x + 1); M
            [x^2 - 2*x + 1             0]
            [            0 x^2 - 2*x + 1]
            sage: M.factor()
            [(x - 1)^2         0]
            [        0 (x - 1)^2]
        """
        cdef Matrix_symbolic_dense M = self._new_c()
        M._maxima = self._maxima.factor()
        return M

    def expand(self):
        """
        Operates pointwise on each element.

        EXAMPLES:
            sage: M = matrix(2, 2, range(4)) - var('x')
            sage: M*M
            [        x^2 + 2         3 - 2*x]
            [2*(3 - x) - 2*x   (3 - x)^2 + 2]
            sage: (M*M).expand()
            [       x^2 + 2        3 - 2*x]
            [       6 - 4*x x^2 - 6*x + 11]
        """
        cdef Matrix_symbolic_dense M = self._new_c()
        M._maxima = self._maxima.expand()
        return M



    def variables(self, vars=tuple([])):
        """
        Returns the variables of self.

        EXAMPLES:
            sage: var('a,b,c,x,y')
            (a, b, c, x, y)
            sage: m = matrix([[x, x+2], [x^2, x^2+2]]); m
            [      x   x + 2]
            [    x^2 x^2 + 2]
            sage: m.variables()
            (x,)
            sage: m = matrix([[a, b+c], [x^2, y^2+2]]); m
            [      a   c + b]
            [    x^2 y^2 + 2]
            sage: m.variables()
            (a, b, c, x, y)
        """
        if self.__variables is not None:
            return self.__variables

        vars = list(set(sum([list(op.variables()) for op in self.list()], list(vars))))

        vars.sort(var_cmp)
        vars = tuple(vars)
        self.__variables = vars
        return vars

    def number_of_arguments(self):
        """
        Returns the number of arguments that self can take.

        EXAMPLES:
            sage: var('a,b,c,x,y')
            (a, b, c, x, y)
            sage: m = matrix([[a, (x+y)/(x+y)], [x^2, y^2+2]]); m
            [      a       1]
            [    x^2 y^2 + 2]
            sage: m.number_of_arguments()
            3

            sage: M = MatrixSpace(SR,2,2)
            sage: m = M(sin+1)
            sage: m.variables()
            ()
            sage: m.number_of_arguments()
            1

        """
        if self.__number_of_args is not None:
            return self.__number_of_args

        variables = self.variables()
        if not self.is_simplified():
            n = self.simplify().number_of_arguments()
        else:
            # We need to do this maximum to correctly handle the case where
            # self is something like (sin+1)
            n = max( max([i.number_of_arguments() for i in self.list()]+[0]), len(variables) )
        self.__number_of_args = n
        return n

    def arguments(self):
        """
        Returns a tuple of the arguments that self can take.

        EXAMPLES:
            sage: var('x,y,z')
            (x, y, z)
            sage: M = MatrixSpace(SR,2,2)
            sage: M(x).arguments()
            (x,)
            sage: M(x+sin).arguments()
            (x,)

            sage: M(sin+1).arguments()
            ()
            sage: M(sin+1).number_of_arguments()
            1
        """
        return self.variables()

    def __call__(self, *args, **kwargs):
        """
        EXAMPLES:
            sage: var('x,y,z')
            (x, y, z)
            sage: M = MatrixSpace(SR,2,2)
            sage: h = M(sin+cos)
            sage: h
            [sin + cos         0]
            [        0 sin + cos]
            sage: h(1)
            [sin(1) + cos(1)               0]
            [              0 sin(1) + cos(1)]
            sage: h(x)
            [sin(x) + cos(x)               0]
            [              0 sin(x) + cos(x)]
            sage: h = 3*M(sin)
            sage: h(1)
            [3*sin(1)        0]
            [       0 3*sin(1)]
            sage: h(x)
            [3*sin(x)        0]
            [       0 3*sin(x)]

            sage: M(sin+1)(1)
            [sin(1) + 1          0]
            [         0 sin(1) + 1]
            sage: M(x+sin)(5)
            [sin(5) + 5          0]
            [         0 sin(5) + 5]

            sage: f = M([0,x,y,z]); f
            [0 x]
            [y z]
            sage: f.arguments()
            (x, y, z)
            sage: f()
            [0 x]
            [y z]
            sage: f(1)
            [0 1]
            [y z]
            sage: f(1,2)
            [0 1]
            [2 z]
            sage: f(1,2,3)
            [0 1]
            [2 3]
            sage: f(x=1)
            [0 1]
            [y z]
            sage: f(x=1,y=2)
            [0 1]
            [2 z]
            sage: f(x=1,y=2,z=3)
            [0 1]
            [2 3]
            sage: f({x:1,y:2,z:3})
            [0 1]
            [2 3]

            sage: f(1, x=2)
            Traceback (most recent call last):
            ...
            ValueError: args and kwargs cannot both be specified
            sage: f(1,2,3,4)
            Traceback (most recent call last):
            ...
            ValueError: the number of arguments must be less than or equal to 3
        """
        if kwargs and args:
            raise ValueError, "args and kwargs cannot both be specified"

        if len(args) == 1 and isinstance(args[0], dict):
            kwargs = dict([(repr(x[0]), x[1]) for x in args[0].iteritems()])

        if kwargs:
            #Handle the case where kwargs are specified
            new_entries = []
            for entry in self.list():
                try:
                    new_entries.append( entry(**kwargs) )
                except ValueError:
                    new_entries.append(entry)
        else:
            #Handle the case where args are specified

            #Get all the variables
            variables = list( self.arguments() )

            if len(args) > self.number_of_arguments():
                raise ValueError, "the number of arguments must be less than or equal to %s"%self.number_of_arguments()

            new_entries = []
            for entry in self.list():
                try:
                    entry_vars = entry.variables()
                    if len(entry_vars) == 0:
                        if len(args) != 0:
                            new_entries.append( entry(args[0]) )
                        else:
                            new_entries.append( entry )
                        continue
                    else:
                        indices = [i for i in map(variables.index, entry_vars) if i < len(args)]
                        if len(indices) == 0:
                            new_entries.append( entry )
                        else:
                            new_entries.append( entry(*[args[i] for i in indices]) )
                except ValueError:
                    new_entries.append( entry )

        return self.parent(new_entries)

    def fcp(self, var='x'):
        """
        Return the factorization of the characteristic polynomial of self.

        INPUT:
            var -- (default: 'x') name of variable of charpoly

        EXAMPLES:
            sage: a = matrix(SR,[[1,2],[3,4]])
            sage: a.fcp()
            x^2 - 5*x - 2
            sage: [i for i in a.fcp()]
            [(x^2 - 5*x - 2, 1)]
            sage: a = matrix(SR,[[1,0],[0,2]])
            sage: a.fcp()
            (x - 2) * (x - 1)
            sage: [i for i in a.fcp()]
            [(x - 2, 1), (x - 1, 1)]
            sage: a = matrix(SR, 5, [1..5^2])
            sage: a.fcp()
            x^3 * (x^2 - 65*x - 250)
            sage: list(a.fcp())
            [(x, 3), (x^2 - 65*x - 250, 1)]

        """
        return Factorization(self.charpoly(var).factor_list())
