r"""
Sparse Matrices over a general ring

EXAMPLES::

    sage: R.<x> = PolynomialRing(QQ)
    sage: M = MatrixSpace(QQ['x'],2,3,sparse=True); M
    Full MatrixSpace of 2 by 3 sparse matrices over Univariate Polynomial Ring in x over Rational Field
    sage: a = M(range(6)); a
    [0 1 2]
    [3 4 5]
    sage: b = M([x^n for n in range(6)]); b
    [  1   x x^2]
    [x^3 x^4 x^5]
    sage: a * b.transpose()
    [            2*x^2 + x           2*x^5 + x^4]
    [      5*x^2 + 4*x + 3 5*x^5 + 4*x^4 + 3*x^3]
    sage: pari(a)*pari(b.transpose())
    [2*x^2 + x, 2*x^5 + x^4; 5*x^2 + 4*x + 3, 5*x^5 + 4*x^4 + 3*x^3]
    sage: c = copy(b); c
    [  1   x x^2]
    [x^3 x^4 x^5]
    sage: c[0,0] = 5; c
    [  5   x x^2]
    [x^3 x^4 x^5]
    sage: b[0,0]
    1
    sage: c.dict()
    {(0, 0): 5, (0, 1): x, (0, 2): x^2, (1, 0): x^3, (1, 1): x^4, (1, 2): x^5}
    sage: c.list()
    [5, x, x^2, x^3, x^4, x^5]
    sage: c.rows()
    [(5, x, x^2), (x^3, x^4, x^5)]
    sage: TestSuite(c).run()
    sage: d = c.change_ring(CC['x']); d
    [5.00000000000000                x              x^2]
    [             x^3              x^4              x^5]
    sage: latex(c)
    \left(\begin{array}{rrr}
    5 & x & x^{2} \\
    x^{3} & x^{4} & x^{5}
    \end{array}\right)
    sage: c.sparse_rows()
    [(5, x, x^2), (x^3, x^4, x^5)]
    sage: d = c.dense_matrix(); d
    [  5   x x^2]
    [x^3 x^4 x^5]
    sage: parent(d)
    Full MatrixSpace of 2 by 3 dense matrices over Univariate Polynomial Ring in x over Rational Field
    sage: c.sparse_matrix() is c
    True
    sage: c.is_sparse()
    True
"""

cimport matrix
cimport matrix_sparse
cimport sage.structure.element
from sage.structure.element cimport ModuleElement

import sage.misc.misc as misc

cdef class Matrix_generic_sparse(matrix_sparse.Matrix_sparse):
    r"""
    The ``Matrix_generic_sparse`` class derives from
    ``Matrix``, and defines functionality for sparse
    matrices over any base ring. A generic sparse matrix is represented
    using a dictionary with keys pairs `(i,j)` and values in
    the base ring.

    The values of the dictionary must never be zero.
    """
    ########################################################################
    # LEVEL 1 functionality
    #   * __cinit__
    #   * __init__
    #   * __dealloc__
    #   * set_unsafe
    #   * get_unsafe
    #   * def _pickle
    #   * def _unpickle
    ########################################################################
    def __cinit__(self, parent, entries=0, coerce=True, copy=True):
        self._entries = {}  # crucial so that pickling works

    def __init__(self, parent, entries=None, coerce=True, copy=True):
        cdef Py_ssize_t i, j
        matrix.Matrix.__init__(self, parent)

        R = self._base_ring
        self._zero = R(0)
        if not isinstance(entries, (list, dict)):
            if entries is None:
                x = R.zero_element()
            else:
                x = R(entries)
            entries = {}
            if x != self._zero:
                if self._nrows != self._ncols:
                    raise TypeError, "scalar matrix must be square"
                for i from 0 <= i < self._nrows:
                    entries[(int(i),int(i))] = x

        if isinstance(entries, list) and len(entries) > 0 and \
                hasattr(entries[0], "is_vector"):
            entries = _convert_sparse_entries_to_dict(entries)

        if isinstance(entries, list):
            if len(entries) != self.nrows() * self.ncols():
                raise TypeError, "entries has the wrong length"
            x = entries
            entries = {}
            k = 0
            for i from 0 <= i < self._nrows:
                for j from 0 <= j < self._ncols:
                    if x[k] != 0:
                        entries[(int(i),int(j))] = x[k]
                    k = k + 1
            copy = False

        if not isinstance(entries, dict):
            raise TypeError, "entries must be a dict"

        if coerce:
            try:
                v = {}
                for k, x in entries.iteritems():
                    w = R(x)
                    if w != 0:
                        v[k] = w
                entries = v
            except TypeError:
                raise TypeError, "Unable to coerce entries to %s"%R
        else:
            if copy:
                # Make a copy
                entries = dict(entries)
            for k in entries.keys():
                if entries[k].is_zero():
                    del entries[k]

        self._entries = entries

    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, value):
        # TODO: maybe make faster with Python/C API
        k = (int(i),int(j))
        if value.is_zero():
            try:
                del self._entries[k]
            except KeyError:
                pass
        else:
            self._entries[k] = value

    cdef get_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        # TODO: maybe make faster with Python/C API
        try:
            return self._entries[(int(i),int(j))]
        except KeyError:
            return self._zero

    def _pickle(self):
        version = 0
        return self._entries, version

    def _unpickle(self, data, int version):
        """
        EXAMPLES::

            sage: a = matrix([[1,10],[3,4]],sparse=True); a
            [ 1 10]
            [ 3  4]
            sage: loads(dumps(a)) == a
            True
        """
        cdef Py_ssize_t i, j, k

        if version == 0:
            self._entries = data
            self._zero = self._base_ring(0)
        else:
            raise RuntimeError, "unknown matrix version (=%s)"%version

    def __richcmp__(matrix.Matrix self, right, int op):  # always need for mysterious reasons.
        return self._richcmp(right, op)
    def __hash__(self):
        return self._hash()

    ########################################################################
    # LEVEL 2 functionality
    # x  * cdef _add_
    #    * cdef _mul_
    #    * cdef _cmp_c_impl
    #    * __neg__
    #    * __invert__
    # x  * __copy__
    #    * _multiply_classical
    # x  * _list -- copy of the list of underlying elements
    # x  * _dict -- copy of the sparse dictionary of underlying elements
    ########################################################################

    cpdef ModuleElement _add_(self, ModuleElement _other):
        """
        EXAMPLES::

            sage: a = matrix([[1,10],[3,4]],sparse=True); a
            [ 1 10]
            [ 3  4]
            sage: a+a
            [ 2 20]
            [ 6  8]

        ::

            sage: a = matrix([[1,10,-5/3],[2/8,3,4]],sparse=True); a
            [   1   10 -5/3]
            [ 1/4    3    4]
            sage: a+a
            [    2    20 -10/3]
            [  1/2     6     8]
        """
        # Compute the sum of two sparse matrices.
        # This is complicated because of how we represent sparse matrices.
        # Algorithm:
        #   1. Sort both entry coefficient lists.
        #   2. March through building a new list, adding when the two i,j are the same.
        cdef Py_ssize_t i, j, len_v, len_w
        cdef Matrix_generic_sparse other
        other = <Matrix_generic_sparse> _other
        v = self._entries.items()
        v.sort()
        w = other._entries.items()
        w.sort()
        s = {}
        i = 0  # pointer into self
        j = 0  # pointer into other
        len_v = len(v)   # could be sped up with Python/C API??
        len_w = len(w)
        while i < len_v and j < len_w:
            vi = v[i][0]
            wj = w[j][0]
            if vi < wj:
                s[vi] = v[i][1]
                i = i + 1
            elif vi > wj:
                s[wj] = w[j][1]
                j = j + 1
            else:  # equal
                sm = v[i][1] + w[j][1]
                if not sm.is_zero():
                    s[vi] = sm
                i = i + 1
                j = j + 1
        while i < len(v):
            s[v[i][0]] = v[i][1]
            i = i + 1
        while j < len(w):
            s[w[j][0]] = w[j][1]
            j = j + 1

        cdef Matrix_generic_sparse A
        A = Matrix_generic_sparse.__new__(Matrix_generic_sparse, self._parent, 0,0,0)
        matrix.Matrix.__init__(A, self._parent)
        A._entries = s
        A._zero = self._zero
        A._base_ring = self._base_ring
        return A

    def __copy__(self):
        A = self.__class__(self._parent, self._entries, copy = True, coerce=False)
        if self._subdivisions is not None:
            A.subdivide(*self.subdivisions())
        return A


    def _list(self):
        """
        Return all entries of self as a list of numbers of rows times
        number of columns entries.
        """
        x = self.fetch('list')
        if not x is None:
            return x
        v = [self._base_ring(0)]*(self._nrows * self._ncols)
        for ij, x in self._entries.iteritems():
            i, j = ij          # todo: optimize
            v[i*self._ncols + j] = x
        self.cache('list', v)
        return v

    def _dict(self):
        """
        Return the underlying dictionary of self.
        """
        return self._entries

    ########################################################################
    # LEVEL 3 functionality -- matrix windows, etc.
    ########################################################################

    def _nonzero_positions_by_row(self, copy=True):
        x = self.fetch('nonzero_positions')
        if x is not None:
            if copy:
                return list(x)
            return x

        v = self._entries.keys()
        v.sort()

        self.cache('nonzero_positions', v)
        if copy:
            return list(v)

        return v

    def _nonzero_positions_by_column(self, copy=True):
        x = self.fetch('nonzero_positions_by_column')
        if x is not None:
            if copy:
                return list(x)
            return x

        v = self._entries.keys()
        v.sort(_cmp_backward)

        self.cache('nonzero_positions_by_column', v)
        if copy:
            return list(v)
        return v

    ######################
    # Echelon support
    ######################


##     def dense_matrix(self):
##         import sage.matrix.matrix
##         M = sage.matrix.matrix.MatrixSpace(self.base_ring(),
##                                            self._nrows, self._ncols,
##                                            sparse = False)(0)
##         for i, j, x in self._entries
##             M[i,j] = x
##         return M

##     def nonpivots(self):
##         # We convert the pivots to a set so we have fast
##         # inclusion testing
##         X = set(self.pivots())
##         # [j for j in xrange(self.ncols()) if not (j in X)]
##         np = []
##         for j in xrange(self.ncols()):
##             if not (j in X):
##                 np.append(j)
##         return np

##     def matrix_from_nonpivot_columns(self):
##         """
##         The sparse matrix got by deleted all pivot columns.
##         """
##         return self.matrix_from_columns(self.nonpivots())

##     def matrix_from_columns(self, columns):
##         """
##         Returns the sparse submatrix of self composed of the given
##         list of columns.

##         INPUT:
##             columns -- list of int's
##         OUTPUT:
##             a sparse matrix.
##         """
##         # Let A be this matrix and let B denote this matrix with
##         # the columns deleted.
##         # ALGORITHM:
##         # 1. Make a table that encodes the function
##         #          f : Z --> Z,
##         #    f(j) = column of B where the j-th column of A goes
##         # 2. Build new list of entries and return resulting matrix.
##         C = set(columns)
##         X = []
##         j = 0
##         for i in xrange(self.ncols()):
##             if i in C:
##                 X.append(j)
##                 j = j + 1
##             else:
##                 X.append(-1)    # column to be deleted.
##         entries2 = []
##         for i, j, x in self.entries():
##             if j in C:
##                 entries2.append((i,X[j],x))
##         return SparseMatrix(self.base_ring(), self.nrows(),
##                             len(C), entries2, coerce=False, sort=False)

##     def matrix_from_rows(self, rows):
##         """
##         Returns the sparse submatrix of self composed of the given
##         list of rows.

##         INPUT:
##             rows -- list of int's
##         OUTPUT:
##             a sparse matrix.
##         """
##         R = set(rows)
##         if not R.issubset(set(xrange(self.nrows()))):
##             raise ArithmeticError, "Invalid rows."
##         X = []
##         i = 0
##         for j in xrange(self.nrows()):
##             if j in R:
##                 X.append(i)
##                 i = i + 1
##             else:
##                 X.append(-1)    # row to be deleted.
##         entries2 = []
##         for i, j, x in self.entries():
##             if i in R:
##                 entries2.append((X[i],j,x))
##         return SparseMatrix(self.base_ring(), len(R),
##                             self.ncols(), entries2, coerce=False, sort=False)


##     def transpose(self):
##         """
##         Returns the transpose of self.
##         """
##         entries2 = [] # [(j,i,x) for i,j,x in self.entries()]
##         for i,j,x in self.entries():
##             entries2.append((j,i,x))
##         return SparseMatrix(self.base_ring(), self.ncols(),
##                             self.nrows(), entries2, coerce=False, sort=True)

##     def __rmul__(self, left):
##         return self.scalar_multiple(left)

##     def scalar_multiple(self, left):
##         R = self.base_ring()
##         left = R(left)
##         if left == R(1):
##             return self
##         if left == R(0):
##             return SparseMatrix(R, self.nrows(), self.ncols(), coerce=False, sort=False)

##         X = []
##         for i, j, x in self.list():
##             X.append((i,j,left*x))
##         return SparseMatrix(self.base_ring(), self.nrows(),
##                             self.ncols(), X, coerce=False, sort=False)


####################################################################################
# Various helper functions
####################################################################################
import matrix_space

def Matrix_sparse_from_rows(X):
    """
    INPUT:


    -  ``X`` - nonempty list of SparseVector rows


    OUTPUT: Sparse_matrix with those rows.

    EXAMPLES::

        sage: V = VectorSpace(QQ,20,sparse=True)
        sage: v = V(0)
        sage: v[9] = 4
        sage: from sage.matrix.matrix_generic_sparse import Matrix_sparse_from_rows
        sage: Matrix_sparse_from_rows([v])
        [0 0 0 0 0 0 0 0 0 4 0 0 0 0 0 0 0 0 0 0]
        sage: Matrix_sparse_from_rows([v, v, v, V(0)])
        [0 0 0 0 0 0 0 0 0 4 0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 4 0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 4 0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
    """
    cdef Py_ssize_t i, j

    if not isinstance(X, (list, tuple)):
        raise TypeError, "X (=%s) must be a list or tuple"%X
    if len(X) == 0:
        raise ArithmeticError, "X must be nonempty"
    entries = {}
    R = X[0].base_ring()
    ncols = X[0].degree()
    for i from 0 <= i < len(X):
        for j, x in X[i].iteritems():
            entries[(i,j)] = x
    M = matrix_space.MatrixSpace(R, len(X), ncols, sparse=True)
    return M(entries, coerce=False, copy=False)


## def _sparse_dot_product(v, w):
##     """
##     INPUT:
##         v and w are dictionaries with integer keys.
##     """
##     x = set(v.keys()).intersection(set(w.keys()))
##     a = 0
##     for k in x:
##         a = a + v[k]*w[k]
##     return a

cdef _convert_sparse_entries_to_dict(entries):
    e = {}
    for i in xrange(len(entries)):
        for j, x in (entries[i].dict()).iteritems():
            e[(i,j)] = x
    return e



def _cmp_backward(x,y):  # todo: speed up via Python/C API
    cdef int c
    # compare two 2-tuples, but in reverse order, i.e., second entry than first
    c = cmp(x[1],y[1])
    if c: return c
    c = cmp(x[0],y[0])
    if c: return c
    return 0
