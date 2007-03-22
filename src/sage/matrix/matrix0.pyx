"""
Base class for matrices, part 0

NOTE: For design documentation see matrix/docs.py.

EXAMPLES:
    sage: matrix(2,[1,2,3,4])
    [1 2]
    [3 4]
"""

################################################################################
#       Copyright (C) 2005, 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL), version 2.
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
################################################################################

include "../ext/stdsage.pxi"
include "../ext/cdefs.pxi"
include "../ext/python.pxi"

import sage.modules.free_module
import sage.misc.latex
import sage.structure.coerce
import sage.rings.integer

from   sage.misc.misc import verbose, get_verbose
from   sage.structure.sequence import _combinations, Sequence

cimport sage.structure.element
from   sage.structure.element    cimport ModuleElement, Element, RingElement, Vector
from   sage.structure.mutability cimport Mutability

from sage.rings.ring cimport CommutativeRing

import sage.modules.free_module

import matrix_misc



cdef class Matrix(sage.structure.element.Matrix):
    r"""
    A generic matrix.

    The \class{Matrix} class is the base class for all matrix classes.
    To create a \class{Matrix}, first create a \class{MatrixSpace},
    then coerce a list of elements into the \class{MatrixSpace}.  See
    the documentation of \class{MatrixSpace} for more details.

    EXAMPLES:

    We illustrate matrices and matrix spaces.  Note that no actual
    matrix that you make should have class Matrix; the class should
    always be derived from Matrix.

        sage: M = MatrixSpace(CDF,2,3); M
        Full MatrixSpace of 2 by 3 dense matrices over Complex Double Field
        sage: a = M([1,2,3,  4,5,6]); a
        [1.0 2.0 3.0]
        [4.0 5.0 6.0]
        sage: type(a)
        <type 'sage.matrix.matrix_complex_double_dense.Matrix_complex_double_dense'>
        sage: parent(a)
        Full MatrixSpace of 2 by 3 dense matrices over Complex Double Field

        sage: matrix(CDF, 2,3, [1,2,3, 4,5,6])
        [1.0 2.0 3.0]
        [4.0 5.0 6.0]
        sage: Mat(CDF,2,3)(range(1,7))
        [1.0 2.0 3.0]
        [4.0 5.0 6.0]

        sage: Q.<i,j,k> = QuaternionAlgebra(QQ, -1,-1)
        sage: matrix(Q,2,1,[1,2])
        [1]
        [2]
    """
    def __init__(self, parent):
        """
        EXAMPLES:
            sage: import sage.matrix.matrix0
            sage: A = sage.matrix.matrix0.Matrix(MatrixSpace(QQ,2))
            sage: type(A)
            <type 'sage.matrix.matrix0.Matrix'>
        """
        self._parent = parent
        self._base_ring = parent.base_ring()
        self._nrows = parent.nrows()
        self._ncols = parent.ncols()
        self._mutability = Mutability(False)
        self._cache = {}

    def copy(self):
        """
        Make a copy of self.  If self is immutable, the copy will be mutable.

        WARNING: The individual elements aren't themselves copied
        (though the list is copied).    This shouldn't matter, since ring
        elements are (almost!) always immutable in SAGE.

        EXAMPLES:
            sage: R.<x> = QQ['x']
            sage: a = matrix(R,2,[x+1,2/3,  x^2/2, 1+x^3]); a
            [  x + 1     2/3]
            [1/2*x^2 x^3 + 1]
            sage: b = a.copy()
            sage: b[0,0] = 5
            sage: b
            [      5     2/3]
            [1/2*x^2 x^3 + 1]
            sage: a
            [  x + 1     2/3]
            [1/2*x^2 x^3 + 1]

            sage: b = copy(a)
            sage: f = b[0,0]; f[0] = 10
            Traceback (most recent call last):
            ...
            IndexError: polynomials are immutable
        """
        return self.__copy__()

    def list(self):
        """
        List of elements of self.  It is safe to change the returned list.

        EXAMPLES:
            sage: R.<x,y> = QQ[]
            sage: a = matrix(R,2,[x,y,x*y, y,x,2*x+y]); a
            [      x       y     x*y]
            [      y       x y + 2*x]
            sage: v = a.list(); v
             [x, y, x*y, y, x, y + 2*x]

        Notice that changing the returned list does not change a (the
        list is a copy):
            sage: v[0] = 25
            sage: a
            [      x       y     x*y]
            [      y       x y + 2*x]
        """
        return list(self._list())

    def _list(self):
        """
        Unsafe version of the list method, mainly for internal use.
        This may return the list of elements, but as an *unsafe*
        reference to the underlying list of the object.  It is might
        be dangerous if you change entries of the returned list.

        EXAMPLES:
        Using _list is potentially fast and memory efficient, but
        very dangerous (at least for generic dense matrices).

            sage: a = matrix(QQ['x,y'],2,range(6)); a
            [0 1 2]
            [3 4 5]
            sage: v = a._list(); v
            [0, 1, 2, 3, 4, 5]

        If you change an entry of the list, the corresponding entry
        of the matrix will be changed (but without clearing any caches
        of computing information about the matrix):
            sage: v[0] = -2/3; v
            [-2/3, 1, 2, 3, 4, 5]
            sage: a._list()
            [-2/3, 1, 2, 3, 4, 5]

        Now the 0,0 entry of the matrix is $-2/3$, which is weird.
            sage: a[0,0]
            -2/3

        See:
            sage: a
            [-2/3    1    2]
            [   3    4    5]
        """
        cdef Py_ssize_t i, j

        x = self.fetch('list')
        if not x is None:
            return x
        x = []
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                x.append(self.get_unsafe(i, j))
#        self.cache('list', x)
        return x

    def dict(self):
        """
        Dictionary of the elements of self with keys pairs (i,j) and
        values the nonzero entries of self.

        It is safe to change the returned dictionary.

        EXAMPLES:
            sage: R.<x,y> = QQ[]
            sage: a = matrix(R,2,[x,y,0, 0,0,2*x+y]); a
            [      x       y       0]
            [      0       0 y + 2*x]
            sage: d = a.dict(); d
            {(0, 1): y, (1, 2): y + 2*x, (0, 0): x}

        Notice that changing the returned list does not change a (the
        list is a copy):
            sage: d[0,0] = 25
            sage: a
            [      x       y       0]
            [      0       0 y + 2*x]
        """
        return dict(self._dict())

    def _dict(self):
        """
        Unsafe version of the dict method, mainly for internal use.
        This may return the dict of elements, but as an *unsafe*
        reference to the underlying dict of the object.  It is might
        be dangerous if you change entries of the returned dict.

        EXAMPLES:
        Using _dict is potentially fast and memory efficient, but
        very dangerous (at least for generic sparse matrices).

            sage: a = matrix(QQ['x,y'],2,range(6), sparse=True); a
            [0 1 2]
            [3 4 5]
            sage: v = a._dict(); v
            {(0, 1): 1, (1, 2): 5, (1, 0): 3, (0, 2): 2, (1, 1): 4}

        If you change a key of the dictionary, the corresponding entry
        of the matrix will be changed (but without clearing any caches
        of computing information about the matrix):
            sage: v[0,1] = -2/3; v
            {(0, 1): -2/3, (1, 2): 5, (1, 0): 3, (0, 2): 2, (1, 1): 4}
            sage: a._dict()
            {(0, 1): -2/3, (1, 2): 5, (1, 0): 3, (0, 2): 2, (1, 1): 4}
            sage: a[0,1]
            -2/3

        But the matrix doesn't know the entry changed, so it returns the cached
        version of its print representation:

            sage: a
            [0 1 2]
            [3 4 5]

        If we change an entry, the cache is cleared, and the correct print
        representation appears:
            sage: a[1,2]=10
            sage: a
            [   0 -2/3    2]
            [   3    4   10]
        """
        d = self.fetch('dict')
        if not d is None:
            return d

        cdef Py_ssize_t i, j
        d = {}
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                x = self.get_unsafe(i, j)
                if x != 0:
                    d[(int(i),int(j))] = x
        self.cache('dict', d)
        return d

    ###########################################################
    # Cache
    ###########################################################
    def _clear_cache(self):
        """
        Clear anything cached about this matrix.
        """
        self.clear_cache()

    cdef clear_cache(self):
        """
        Clear the properties cache.
        """
        self._cache = None

    cdef fetch(self, key):
        """
        Try to get an element from the cache; if there isn't anything there, return None.
        """
        if self._cache is None:
            return None
        try:
            return self._cache[key]
        except KeyError:
            return None

    cdef cache(self, key, x):
        """
        Record x in the cache with given key.
        """
        if self._cache is None:
            self._cache = {}
        self._cache[key] = x

    def _get_cache(self):
        return self._cache

    ###########################################################
    # Mutability and bounds checking
    ###########################################################

    cdef check_bounds(self, Py_ssize_t i, Py_ssize_t j):
        """
        This function gets called when you're about to access the i,j
        entry of this matrix.  If i, j are out of range, and
        IndexError is raised.
        """
        if i<0 or i >= self._nrows or j<0 or j >= self._ncols:
            raise IndexError, "matrix index out of range"

    cdef check_mutability(self):
        """
        This function gets called when you're about to change this matrix.

        If self is immutable, a ValueError is raised, since you should
        never change a mutable matrix.

        If self is mutable, the cache of results about self is deleted.
        """
        if self._mutability._is_immutable:
            raise ValueError, "matrix is immutable; please change a copy instead (use self.copy())."
        else:
            self._cache = {}

    cdef check_bounds_and_mutability(self, Py_ssize_t i, Py_ssize_t j):
        """
        This function gets called when you're about to set the i,j
        entry of this matrix.  If i or j is out of range, an
        IndexError exception is raised.

        If self is immutable, a ValueError is raised, since you should
        never change a mutable matrix.

        If self is mutable, the cache of results about self is deleted.
        """
        if self._mutability._is_immutable:
            raise ValueError, "matrix is immutable; please change a copy instead (use self.copy())."
        else:
            self._cache = {}

        if i<0 or i >= self._nrows or j<0 or j >= self._ncols:
            raise IndexError, "matrix index out of range"


    def set_immutable(self):
        r"""
        Call this function to matrix a matrix immutable.

        Matrices are always mutable by default, i.e., you can change
        their entries using \code{A[i,j] = x}.  However, mutable
        matrices aren't hasheable, so can't be used as keys in
        dictionaries, etc.  Also, often when implementing a class, you
        might compute a matrix associated to it, e.g., the matrix of a
        Hecke operator.  If you return this matrix to the user you're
        really returning a reference and the user could then change an
        entry; this could be confusing. Thus you shoulds set such a
        matrix immutable.

        EXAMPLES:
            sage: A = Matrix(QQ, 2, 2, range(4))
            sage: A.is_mutable()
            True
            sage: A[0,0] = 10
            sage: A
            [10   1]
            [ 2   3]

        Mutable matrices are not hasheable, so can't be used as keys for dictionaries:
            sage: hash(A)
            Traceback (most recent call last):
            ...
            TypeError: mutable matrices are unhashable
            sage: v = {A:1}
            Traceback (most recent call last):
            ...
            TypeError: mutable matrices are unhashable


        If we make A immutable it suddenly is hasheable.
            sage: A.set_immutable()
            sage: A.is_mutable()
            False
            sage: A[0,0] = 10
            Traceback (most recent call last):
            ...
            ValueError: matrix is immutable; please change a copy instead (use self.copy()).
            sage: hash(A)
            12
            sage: v = {A:1}; v
            {[10  1]
             [ 2  3]: 1}
        """
        self._mutability.set_immutable()

    def is_immutable(self):
        """
        Return True if this matrix is immutable.

        See the documentation for self.set_immutable for more details
        about mutability.

        EXAMPLES:
            sage: A = Matrix(QQ['t','s'], 2, 2, range(4))
            sage: A.is_immutable()
            False
            sage: A.set_immutable()
            sage: A.is_immutable()
            True
        """
        return bool(self._mutability._is_immutable)

    def is_mutable(self):
        """
        Return True if this matrix is mutable.

        See the documentation for self.set_immutable for more details
        about mutability.

        EXAMPLES:
            sage: A = Matrix(QQ['t','s'], 2, 2, range(4))
            sage: A.is_mutable()
            True
            sage: A.set_immutable()
            sage: A.is_mutable()
            False
        """
        return bool(self._mutability.is_mutable())

    ###########################################################
    # Entry access
    #    The first two must be overloaded in the derived class
    ###########################################################
    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, object x):
        """
        Set entry quickly without doing any bounds checking.  Calling
        this with invalid arguments is allowed to produce a
        segmentation fault.

        This is fast since it is a cdef function and there is no bounds checking.
        """
        raise NotImplementedError, "this must be defined in the derived class (type=%s)"%type(self)

    cdef get_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        """
        Entry access, but fast since it might be without
        bounds checking.

        This is fast since it is a cdef function and there is no bounds checking.
        """
        raise NotImplementedError, "this must be defined in the derived type."

    def __iter__(self):
        return matrix_misc.row_iterator(self)

    def __getitem__(self, ij):
        """
        Return element or row of self.

        INPUT:
            ij -- tuple (i,j) with i, j integers
        or
            ij -- integer

        USAGE:
            A[i, j] -- the i,j of A, and
            A[i]    -- the i-th row of A.

        EXAMPLES:
            sage: A = Matrix(Integers(2006),2,2,[-1,2,3,4])
            sage: A[0,0]
            2005
            sage: A[0]
            (2005, 2)

            sage: a = MatrixSpace(ZZ,3)(range(9)); a
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: a[1,2]
            5
            sage: a[0]
            (0, 1, 2)
            sage: a[4,7]
            Traceback (most recent call last):
            ...
            IndexError: matrix index out of range
            sage: a[-1,0]
            Traceback (most recent call last):
            ...
            IndexError: matrix index out of range
        """
        cdef Py_ssize_t i, j
        cdef object x

        if PyTuple_Check(ij):
            # ij is a tuple, so we get i and j efficiently, construct corresponding integer entry.
            if PyTuple_Size(ij) != 2:
                raise IndexError, "index must be an integer or pair of integers"
            i = <object> PyTuple_GET_ITEM(ij, 0)
            j = <object> PyTuple_GET_ITEM(ij, 1)
            self.check_bounds(i, j)
            return self.get_unsafe(i, j)
        else:
            # If ij is not a tuple, coerce to an integer and get the row.
            i = ij
            return self.row(i)

    def __getslice__(self,  Py_ssize_t i,  Py_ssize_t j):
        """
        Returns a submatrix of this matrix got from the i-th through j-th rows.

        USAGE:
            A[i:j] -- return submatrix from those rows.

        EXAMPLES:
            sage: n=10;a=matrix(QQ,n,range(n^2))
            sage: a[0:3]
            [ 0  1  2  3  4  5  6  7  8  9]
            [10 11 12 13 14 15 16 17 18 19]
            [20 21 22 23 24 25 26 27 28 29]
        """
        if i < 0: i = self._nrows + i
        if j < 0: j = self._nrows + j
        if i >= self._nrows:
            i = self._nrows - 1
        if j >= self._nrows:
            j = self._nrows - 1
        return self.matrix_from_rows(range(i,j))

    def __setitem__(self, ij, x):
        """
        Set position i,j of this matrix to x.

        INPUT:
            ij -- tuple (i,j), where i is the row and j the column
                  Alternatively, ij can be an integer, and the ij-th row is set.
            x -- something that can be coerced to the base ring of this matrix.

        USAGE:
            A[i, j] = x -- set the (i,j) entry of A
            A[i]    = x -- set the ith row of A

        EXAMPLES:
            sage: A = Matrix(Integers(2006),2,2,[-1,2,3,4]); A
            [2005    2]
            [   3    4]
            sage: A[0,0] = 5; A
            [5 2]
            [3 4]
            sage: A[0] = [2,3]
            sage: A
            [2 3]
            [3 4]
            sage: A.set_immutable()
            sage: A[0,0] = 7
            Traceback (most recent call last):
            ...
            ValueError: matrix is immutable; please change a copy instead (use self.copy()).

            sage: a = matrix(ZZ,2,3, range(6)); a
            [0 1 2]
            [3 4 5]
            sage: a[0,0] = 10
            sage: a
            [10  1  2]
            [ 3  4  5]
        """
        cdef Py_ssize_t i, j

        if PyTuple_Check(ij):
            # ij is a tuple, so we get i and j efficiently, construct corresponding integer entry.
            if PyTuple_Size(ij) != 2:
                raise IndexError, "index must be an integer or pair of integers"
            i = <object> PyTuple_GET_ITEM(ij, 0)
            j = <object> PyTuple_GET_ITEM(ij, 1)
            self.check_bounds_and_mutability(i, j)
            self.set_unsafe(i, j, self._coerce_element(x))
        else:
            # If ij is not a tuple, coerce to an integer and set the whole row.
            i = ij
            if self._ncols == 0: # degenerate case
                return
            self.check_bounds_and_mutability(i, 0)
            for j from 0 <= j < self._ncols:
                self.set_unsafe(i, j, self._coerce_element(x[j]))

    cdef _coerce_element(self, x):
        """
        Return coercion of x into the base ring of self.
        """
        if PY_TYPE_CHECK(x, Element) and (<Element> x)._parent is self._base_ring:
            return x
        return self._base_ring(x)

    ###########################################################
    # Pickling
    ###########################################################
    def __reduce__(self):
        """
        EXAMPLES:
            sage: a = matrix(Integers(8),3,range(9))
            sage: a == loads(dumps(a))
            True
        """
        data, version = self._pickle()
        return unpickle, (self.__class__, self._parent, self._mutability,
                                          self._cache, data, version)

    def _pickle(self):
        raise NotImplementedError

    def _test_pickle(self):
        a, b = self.__reduce__()
        print "__reduce__:"
        print a
        print b
        print "Now call a with b:"
        c = a(*b)
        print "Got c = ", c
        if self == c:
            print "Pickle success."
        else:
            print "Pickle failure."

    ###########################################################
    # Base Change
    ###########################################################
    def base_ring(self):
        return self._base_ring

    def change_ring(self, ring):
        """
        Return the matrix obtained by coercing the entries of this
        matrix into the given ring.

        Always returns a copy (unless self is immutable, in which case
        returns self).

        EXAMPLES:
            sage: A = Matrix(QQ, 2, 2, [1/2, 1/3, 1/3, 1/4])
            sage: A.parent()
             Full MatrixSpace of 2 by 2 dense matrices over Rational Field
            sage: A.change_ring(GF(25,'a'))
            [3 2]
            [2 4]
            sage: A.change_ring(GF(25,'a')).parent()
             Full MatrixSpace of 2 by 2 dense matrices over Finite Field in a of size 5^2
            sage: A.change_ring(ZZ)
            Traceback (most recent call last):
            ...
            TypeError: matrix has denominators so can't change to ZZ.
        """
        if ring is self._base_ring:
            if self._mutability._is_immutable:
                return self
            return self.copy()

        M = sage.matrix.matrix_space.MatrixSpace(ring, self._nrows, self._ncols, sparse=self.is_sparse())
        return M(self.list(), coerce=True, copy=False)

    def _matrix_(self, R):
        """
        EXAMPLES:
            sage: A = Matrix(ZZ[['t']], 2, 2, range(4))
            sage: A.parent()
            Full MatrixSpace of 2 by 2 dense matrices over Power Series Ring in t over Integer Ring
            sage: A._matrix_(QQ[['t']])
            [0 1]
            [2 3]
            sage: A._matrix_(QQ[['t']]).parent()
            Full MatrixSpace of 2 by 2 dense matrices over Power Series Ring in t over Rational Field
        """
        return self.change_ring(R)

    ###########################################################
    # Representation -- string, latex, etc.
    ###########################################################
    def __repr__(self):
        if self._nrows < max_rows and self._ncols < max_cols:
            return self.str()
        if self.is_sparse():
            s = 'sparse'
        else:
            s = 'dense'
        return "%s x %s %s matrix over %s"%(self._nrows, self._ncols, s, self.base_ring())

    def str(self):
        r"""
        EXAMPLES:
            sage: R = PolynomialRing(QQ,6,'z')
            sage: a = matrix(2,3, R.gens())
            sage: a.__repr__()
            '[z0 z1 z2]\n[z3 z4 z5]'
        """
        #x = self.fetch('repr')  # too confusing!!
        #if not x is None: return x
        cdef Py_ssize_t nr, nc, r, c
        nr = self._nrows
        nc = self._ncols

        if nr == 0 or nc == 0:
            return "[]"

        # compute column widths
        S = []
        for x in self.list():
            S.append(str(x))

        tmp = []
        for x in S:
            tmp.append(len(x))

        width = max(tmp)
        rows = []
        m = 0

        # compute rows
        for r from 0 <= r < nr:
            s = ""
            for c from 0 <= c < nc:
                if c == nc - 1:
                    sep=""
                else:
                    sep=" "
                entry = S[r*nc+c]
                if c == 0:
                    m = max(m, len(entry))
                entry = " "*(width-len(entry)) + entry
                s = s + entry + sep
            rows.append(s)

        # Put brackets around in a single string
        tmp = []
        for row in rows:
            tmp.append("[%s]"%row)
        s = "\n".join(tmp)
        #self.cache('repr',s)
        return s

##     def _latex_sparse(self, variable="x"):
##         r"""
##         Return a latex string that represents this matrix as a sparse
##         matrix.  The rows are printed as sums $\sum a_i x_i$, where
##         $x$ is the variable.

##         EXAMPLES:

##         """
##         cdef Py_ssize_t nr, nc, i, j
##         nr = self._nrows
##         nc = self._ncols
##         s = "\\left(\\begin{align*}\n"
##         for i from 0 <= i < nr:
##             v = []
##             for j from 0 <= j < nc:
##                 x = self.get_unsafe(i, j)
##                 if x != 0:
##                     v.append((j, x))
##             for j from 0 <= j < len(v):
##                 s  = s + "%s*%s_{%s}"%(v[j][1], variable, v[j][0])
##                 if j == 0:
##                     s = s + "& + "
##                 elif j < len(v) - 1:
##                     s =  s + " + "
##                 else:
##                     s =  s + "\\\\\n"
##         s = s + "\n\\end{align*}"
##         return s

    def _latex_(self):
        r"""
        Return latex representation of this matrix.

        EXAMPLES:
            sage: R = PolynomialRing(QQ,4,'z')
            sage: a = matrix(2,2, R.gens())
            sage: b = a*a
            sage: latex(b)
            \left(\begin{array}{rr}
            z_{1}z_{2} + z_{0}^{2}&z_{1}z_{3} + z_{0}z_{1}\\
            z_{2}z_{3} + z_{0}z_{2}&z_{3}^{2} + z_{1}z_{2}
            \end{array}\right)
        """
        cdef Py_ssize_t nr, nc, r, c
        nr = self._nrows
        nc = self._ncols
        if nr == 0 or nc == 0:
            return "()"

        # compute column widths
        S = []
        for x in self.list():
            S.append(str(x))

        tmp = []
        for x in S:
            tmp.append(len(x))

        width = max(tmp)

        S = self.list()
        rows = []
        m = 0

        # compute rows
        latex = sage.misc.latex.latex
        for r from 0 <= r < nr:
            s = ""
            for c from 0 <= c < nc:
                if c == nc-1:
                    sep=""
                else:
                    sep="&"
                entry = latex(S[r*nc+c])
                if c == 0:
                    m = max(m, len(entry))
                s = s + entry + sep
            rows.append(s)

        # Put brackets around in a single string
        tmp = []
        for row in rows:
            tmp.append("%s"%row)
        s = "\\\\\n".join(tmp)
        return "\\left(\\begin{array}{%s}\n"%('r'*nc) + s + "\n\\end{array}\\right)"



    ###################################################
    ## Basic Properties
    ###################################################
    def ncols(self):
        """
        Return the number of rows of this matrix.

        EXAMPLES:
            sage: M = MatrixSpace(QQ, 2, 3)
            sage: A = M([1,2,3, 4,5,6])
            sage: A
            [1 2 3]
            [4 5 6]
            sage: A.ncols()
            3
            sage: A.nrows()
            2

        AUTHORS:
            -- Naqi Jaffery (2006-01-24): examples
        """
        return self._ncols

    def nrows(self):
        r"""
        Return the number of rows of this matrix.

        EXAMPLES:
            sage: M = MatrixSpace(QQ,6,7)
            sage: A = M([1,2,3,4,5,6,7, 22,3/4,34,11,7,5,3, 99,65,1/2,2/3,3/5,4/5,5/6, 9,8/9, 9/8,7/6,6/7,76,4, 0,9,8,7,6,5,4, 123,99,91,28,6,1024,1])
            sage: A
            [   1    2    3    4    5    6    7]
            [  22  3/4   34   11    7    5    3]
            [  99   65  1/2  2/3  3/5  4/5  5/6]
            [   9  8/9  9/8  7/6  6/7   76    4]
            [   0    9    8    7    6    5    4]
            [ 123   99   91   28    6 1024    1]
            sage: A.ncols()
            7
            sage: A.nrows()
            6

        AUTHORS:
            -- Naqi Jaffery (2006-01-24): examples
        """
        return self._nrows


    ###################################################
    # Functions
    ###################################################
    def act_on_polynomial(self, f):
        """
        Returns the polynomial f(self*x).

        INPUT:
            self -- an nxn matrix
            f -- a polynomial in n variables x=(x1,...,xn)

        OUTPUT:
            The polynomial f(self*x).

        EXAMPLES:
            sage: R.<x,y> = QQ[]
            sage: x, y = R.gens()
            sage: f = x**2 - y**2
            sage: M = MatrixSpace(QQ, 2)
            sage: A = M([1,2,3,4])
            sage: A.act_on_polynomial(f)
            -12*y^2 - 20*x*y - 8*x^2
        """
        cdef Py_ssize_t i, j, n

        if self._nrows != self._ncols:
            raise ArithmeticError, "self must be square"

        F = f.base_ring()
        vars = f.parent().gens()
        n = len(self.rows())
        ans = []
        for i from 0 <= i < n:
            tmp = []
            for j from 0 <= j < n:
                tmp.append(self.get_unsafe(i,j)*vars[j])
            ans.append( sum(tmp) )
        return f(tuple(ans))

    ###################################################
    # Arithmetic
    ###################################################
    def commutator(self, other):
        """
        Return the commutator self*other - other*self.

        EXAMPLES:
            sage: A = Matrix(QQ[['t']], 2, 2, range(4))
        """
        return self*other - other*self

    ###################################################
    # Row and column operations
    # The _c versions do no bounds checking and all
    # assume input values have parent that is self._base_ring.
    ###################################################
    cdef check_row_bounds_and_mutability(self, Py_ssize_t r1, Py_ssize_t r2):
        if self._mutability._is_immutable:
            raise ValueError, "matrix is immutable; please change a copy instead (use self.copy())."
        else:
            self._cache = {}
        if r1<0 or r1 >= self._nrows or r2<0 or r2 >= self._nrows:
            raise IndexError, "matrix row index out of range"

    cdef check_column_bounds_and_mutability(self, Py_ssize_t c1, Py_ssize_t c2):
        if self._mutability._is_immutable:
            raise ValueError, "matrix is immutable; please change a copy instead (use self.copy())."
        else:
            self._cache = {}
        if c1<0 or c1 >= self._ncols or c2<0 or c2 >= self._ncols:
            raise IndexError, "matrix column index out of range"

    def swap_columns(self, Py_ssize_t c1, Py_ssize_t c2):
        """
        Swap columns c1 and c2 of self.

        EXAMPLES:
            sage: M = MatrixSpace(QQ,3,3)
            sage: A = M([1,9,-7,4/5,4,3,6,4,3])
            sage: A
            [  1   9  -7]
            [4/5   4   3]
            [  6   4   3]
            sage: A.swap_columns(1,2); A
            [  1  -7   9]
            [4/5   3   4]
            [  6   3   4]
        """
        self.check_column_bounds_and_mutability(c1, c2)
        if c1 != c2:
            self.swap_columns_c(c1, c2)

    cdef swap_columns_c(self, Py_ssize_t c1, Py_ssize_t c2):
        cdef Py_ssize_t r
        for r from 0 <= r < self._nrows:
            a = self.get_unsafe(r, c2)
            self.set_unsafe(r, c2, self.get_unsafe(r,c1))
            self.set_unsafe(r, c1, a)

    def swap_rows(self, r1, r2):
        """
        Swap rows r1 and r2 of self.

        EXAMPLES:
            sage: M = MatrixSpace(QQ,3,3)
            sage: A = M([1,9,-7,4/5,4,3,6,4,3])
            sage: A
            [  1   9  -7]
            [4/5   4   3]
            [  6   4   3]
            sage: A.swap_rows(0,2); A
            [  6   4   3]
            [4/5   4   3]
            [  1   9  -7]
        """
        self.check_row_bounds_and_mutability(r1, r2)
        if r1 != r2:
            self.swap_rows_c(r1, r2)

    cdef swap_rows_c(self, Py_ssize_t r1, Py_ssize_t r2):
        cdef Py_ssize_t c
        for c from 0 <= c < self._ncols:
            a = self.get_unsafe(r2, c)
            self.set_unsafe(r2, c, self.get_unsafe(r1, c))
            self.set_unsafe(r1, c, a)


    def add_multiple_of_row(self, Py_ssize_t i, Py_ssize_t j,    s,   Py_ssize_t col_start=0):
        """
        Add s times row j to row i.

        EXAMPLES:
            sage: a = matrix(2,3,range(6)); a
            [0 1 2]
            [3 4 5]
            sage: a.add_multiple_of_row(1,0,-3)
            sage: a
            [ 0  1  2]
            [ 3  1 -1]
        """
        self.check_row_bounds_and_mutability(i,j)
        s = self._coerce_element(s)
        self.add_multiple_of_row_c(i, j, s, col_start)

    cdef add_multiple_of_row_c(self, Py_ssize_t i, Py_ssize_t j,    s,   Py_ssize_t col_start):
        cdef Py_ssize_t c
        for c from col_start <= c < self._ncols:
            self.set_unsafe(i, c, self.get_unsafe(i, c) + s*self.get_unsafe(j, c))

    def add_multiple_of_column(self, Py_ssize_t i, Py_ssize_t j, s, Py_ssize_t row_start=0):
        """
        Add s times column j to column i.

        EXAMPLES:
            sage: a = matrix(QQ,2,3,range(6)); a
            [0 1 2]
            [3 4 5]
            sage: a.add_multiple_of_column(1,2,-1/2)
            sage: a
            [  0   0   2]
            [  3 3/2   5]
        """
        self.check_column_bounds_and_mutability(i,j)
        s = self._coerce_element(s)
        self.add_multiple_of_column_c(i, j, s, row_start)

    cdef add_multiple_of_column_c(self, Py_ssize_t i, Py_ssize_t j, s, Py_ssize_t row_start):
        cdef Py_ssize_t r
        for r from row_start <= r < self._nrows:
            self.set_unsafe(r, i, self.get_unsafe(r, i) + s*self.get_unsafe(r, j))

    def rescale_row(self, Py_ssize_t i, s, Py_ssize_t start_col=0):
        """
        Replace i-th row of self by s times i-th row of self.

        start_row -- only rescale enries at that column and to the right

        EXAMPLES:
            sage: a = matrix(2,3,range(6)); a
            [0 1 2]
            [3 4 5]
            sage: a.add_multiple_of_row(1,0,-3); a
            [ 0  1  2]
            [ 3  1 -1]
        """
        self.check_row_bounds_and_mutability(i, i)
        s = self._coerce_element(s)
        self.rescale_row_c(i, s, start_col)

    cdef rescale_row_c(self, Py_ssize_t i, s, Py_ssize_t start_col):
        cdef Py_ssize_t j
        for j from start_col <= j < self._ncols:
            self.set_unsafe(i, j, self.get_unsafe(i, j)*s)


    def rescale_col(self, Py_ssize_t i, s, Py_ssize_t start_row=0):
        """
        Replace i-th col of self by s times i-th col of self.

        INPUT:
            i -- ith column
            s -- scalar
            start_row -- only rescale enries at this row and lower

        EXAMPLES:
        We rescale the last column of a matrix over the rational numbers:

            sage: a = matrix(QQ,2,3,range(6)); a
            [0 1 2]
            [3 4 5]
            sage: a.rescale_col(2,1/2); a
            [  0   1   1]
            [  3   4 5/2]
            sage: R.<x> = QQ[]

        We rescale a matrix over a polynomial ring.
            sage: a = matrix(R,2,3,[1,x,x^2,x^3,x^4,x^5]); a
            [  1   x x^2]
            [x^3 x^4 x^5]
            sage: a.rescale_col(2,1/2); a
            [      1       x 1/2*x^2]
            [    x^3     x^4 1/2*x^5]

        We try and fail to rescale a matrix over the integers by a non-integer:
            sage: a = matrix(ZZ,2,3,[0,1,2, 3,4,4]); a
            [0 1 2]
            [3 4 4]
            sage: a.rescale_col(2,1/2); a
            Traceback (most recent call last):
            ...
            TypeError: no coercion of this rational to integer

        We rescale the integer matrix's column 2 column by $-1$, which is an integer.
            sage: a.rescale_col(2,-1); a
            [ 0  1 -2]
            [ 3  4 -4]

        To rescale the matrix by 1/2, you must change the base ring to the rationals:
            sage: b = a.change_ring(QQ); b
            [ 0  1 -2]
            [ 3  4 -4]
            sage: b.rescale_col(2, 1/2); b
            [ 0  1 -1]
            [ 3  4 -2]
        """
        self.check_column_bounds_and_mutability(i, i)
        s = self._coerce_element(s)
        self.rescale_col_c(i, s, start_row)

    cdef rescale_col_c(self, Py_ssize_t i, s, Py_ssize_t start_row):
        cdef Py_ssize_t j
        for j from start_row <= j < self._nrows:
            self.set_unsafe(j, i, self.get_unsafe(j, i)*s)

    def set_row_to_multiple_of_row(self, Py_ssize_t i, Py_ssize_t j, s):
        """
        Set row i equal to s times row j.

        EXAMPLES:
            sage: a = matrix(QQ,2,3,range(6)); a
            [0 1 2]
            [3 4 5]
            sage: a.set_row_to_multiple_of_row(1,0,-3)
            sage: a
            [ 0  1  2]
            [ 0 -3 -6]
        """
        self.check_row_bounds_and_mutability(i,j)
        cdef Py_ssize_t n
        s = self._base_ring(s)
        for n from 0 <= n < self._ncols:
            self.set_unsafe(i, n, s * self.get_unsafe(j, n))  # self[i] = s*self[j]

    def _set_row_to_negative_of_row_of_A_using_subset_of_columns(self, Py_ssize_t i, Matrix A,
                                                                 Py_ssize_t r, cols,
                                                                 cols_index=None):
        """
        Set row i of self to -(row r of A), but where we only take
        the given column positions in that row of A.  We do not
        zero out the other entries of self's row i either.


        INPUT:
            i -- integer, index into the rows of self
            A -- a matrix
            r -- integer, index into rows of A
            cols -- a *sorted* list of integers.
            (cols_index -- ignored)

        EXAMPLES:
            sage: a = matrix(QQ,2,3,range(6)); a
            [0 1 2]
            [3 4 5]
            sage: a._set_row_to_negative_of_row_of_A_using_subset_of_columns(0,a,1,[1,2])
            sage: a
            [-4 -5  2]
            [ 3  4  5]
        """
        self.check_row_bounds_and_mutability(i,i)
        if r < 0 or r >= A.nrows():
            raise IndexError, "invalid row"
        # this function exists just because it is useful for modular symbols presentations.
        cdef Py_ssize_t l
        l = 0
        for k in cols:
            self.set_unsafe(i,l,-A.get_unsafe(r,k))               #self[i,l] = -A[r,k]
            l += 1

    ###################################################
    # Matrix-vector multiply
    ###################################################
    def linear_combination_of_rows(self, v):
        """
        Return the linear combination of the rows of self given by the
        coefficients in the list v.

        INPUT:
            v -- list of length at most the number of rows of self.

        EXAMPLES:
            sage: a = matrix(QQ,2,3,range(6)); a
            [0 1 2]
            [3 4 5]
            sage: a.linear_combination_of_rows([1,2])
            (6, 9, 12)
            sage: a.linear_combination_of_rows([0,0])
            (0, 0, 0)
        """
        cdef Py_ssize_t i, n
        R = self.rows()
        n = len(R)
        if len(v) > n:
            raise ValueError, "length of v (=%s) must be at most the number (=%s) of rows."%(len(v), n)
        zero = self._base_ring(0)
        s = None
        for i from 0 <= i < len(v):
            if v[i] != zero:
                a = v[i] * R[i]
                if s is None:
                    s = a
                else:
                    s += a
        if s is None:
            return self.parent().row_space()(0)
        return s

    def linear_combination_of_columns(self, v):
        """
        Return the linear combination of the columns of self given by the coefficients in
        the list v.

        INPUT:
            v -- list

        EXAMPLES:
            sage: a = matrix(QQ,2,3,range(6)); a
            [0 1 2]
            [3 4 5]
            sage: a.linear_combination_of_columns([1,1,1])
            (3, 12)
            sage: a.linear_combination_of_columns([0,0,0])
            (0, 0)
        """
        cdef Py_ssize_t i, n
        C = self.columns()
        n = len(C)
        if len(v) != n:
            raise ValueError, "length of v must equal number of columns."
        zero = self._base_ring(0)
        s = None
        for i from 0 <= i < n:
            if v[i] != zero:
                a = v[i] * C[i]
                if s is None:
                    s = a
                else:
                    s = s + a
        if s is None:
            return self.parent().column_space()(0)
        return s


    ###################################################
    # Predicates
    ###################################################

    def is_dense(self):
        """
        Returns True if this is a dense matrix.

        In SAGE, being dense is a property of the underlying
        representation, not the number of nonzero entries.

        EXAMPLES:
            sage: matrix(QQ,2,2,range(4)).is_dense()
            True
            sage: matrix(QQ,2,2,range(4),sparse=True).is_dense()
            False
        """
        return bool(self.is_dense_c())

    def is_sparse(self):
        """
        Return True if this is a sparse matrix.

        In SAGE, being sparse is a property of the underlying
        representation, not the number of nonzero entries.

        EXAMPLES:
            sage: matrix(QQ,2,2,range(4)).is_sparse()
            False
            sage: matrix(QQ,2,2,range(4),sparse=True).is_sparse()
            True
        """
        return bool(self.is_sparse_c())

    def is_square(self):
        """
        Return True precisely if this matrix is square, i.e., has the same
        number of rows and columns.

        EXAMPLES:
            sage: matrix(QQ,2,2,range(4)).is_square()
            True
            sage: matrix(QQ,2,3,range(6)).is_square()
            False
        """
        return bool(self._nrows == self._ncols)

    def is_invertible(self):
        r"""
        Return True if this matrix is invertible.

        EXAMPLES:
        The following matrix is invertible over $\Q$ but not over $\Z$.
            sage: A = MatrixSpace(ZZ, 2)(range(4))
            sage: A.is_invertible()
            False
            sage: A.matrix_over_field().is_invertible()
            True

        The inverse function is a constructor for matrices over the
        fraction field, so it can work even if A is not invertible.
            sage: ~A   # inverse of A
            [-3/2  1/2]
            [   1    0]

        The next matrix is invertible over $\Z$.
            sage: A = MatrixSpace(IntegerRing(),2)([1,10,0,-1])
            sage: A.is_invertible()
            True
            sage: ~A                # compute the inverse
            [ 1 10]
            [ 0 -1]

        The following nontrivial matrix is invertible over $\Z[x]$.
            sage: R.<x> = PolynomialRing(IntegerRing())
            sage: A = MatrixSpace(R,2)([1,x,0,-1])
            sage: A.is_invertible()
            True
            sage: ~A
            [ 1  x]
            [ 0 -1]
        """
        return self.is_square() and self.determinant().is_unit()

    ###################################################
    # Invariants of a matrix
    ###################################################

    def pivots(self):
        """
        Return the pivot column positions of this matrix as a list of Python integers.

        This returns a list, of the position of the first nonzero entry in each row
        of the echelon form.

        OUTPUT:
             list -- a list of Python ints

        EXAMPLES:
        """
        x = self.fetch('pivots')
        if not x is None: return x
        self.echelon_form()
        x = self.fetch('pivots')
        if x is None:
            print self
            print self.nrows()
            print self.dict()
            self.save('/home/was/a')
            raise RuntimeError, "BUG: matrix pivots should have been set but weren't, matrix parent = '%s'"%self.parent()
        return x

    def rank(self):
        x = self.fetch('rank')
        if not x is None: return x
        r = len(self.pivots())
        self.cache('rank', r)
        return r

    cdef _set_pivots(self, X):
        self.cache('pivots', X)

    def nonpivots(self):
        """
        Return the list of i such that the i-th column of self
        is NOT a pivot column of the reduced row echelon form
        of self.

        OUTPUT:
            list -- sorted list of (Python) integers

        EXAMPLES:
            sage: a = matrix(QQ,3,3,range(9)); a
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: a.echelon_form()
            [ 1  0 -1]
            [ 0  1  2]
            [ 0  0  0]
            sage: a.nonpivots()
            [2]
        """
        x = self.fetch('nonpivots')
        if not x is None: return x

        X = set(self.pivots())
        np = []
        for j in xrange(self.ncols()):
            if not (j in X):
                np.append(j)
        self.cache('nonpivots',np)
        return np

    def nonzero_positions(self, copy=True, column_order=False):
        """
        Returns the sorted list of pairs (i,j) such that self[i,j] != 0.

        INPUT:
            copy -- (default: True) It is safe to change the resulting
                    list (unless you give the option copy=False).
            column_order -- (default: False)  If true, returns the list of pairs (i,j)
                    such that self[i,j] != 0, but sorted by columns, i.e.,
                    column j=0 entries occur first, then column j=1 entries, etc.

        EXAMPLES:
            sage: a = matrix(QQ, 2,3, [1,2,0,2,0,0]); a
            [1 2 0]
            [2 0 0]
            sage: a.nonzero_positions()
            [(0, 0), (0, 1), (1, 0)]
            sage: a.nonzero_positions(copy=False)
            [(0, 0), (0, 1), (1, 0)]
            sage: a.nonzero_positions(column_order=True)
            [(0, 0), (1, 0), (0, 1)]
            sage: a = matrix(QQ, 2,3, [1,2,0,2,0,0], sparse=True); a
            [1 2 0]
            [2 0 0]
            sage: a.nonzero_positions()
            [(0, 0), (0, 1), (1, 0)]
            sage: a.nonzero_positions(copy=False)
            [(0, 0), (0, 1), (1, 0)]
            sage: a.nonzero_positions(column_order=True)
            [(0, 0), (1, 0), (0, 1)]
        """
        if column_order:
            return self._nonzero_positions_by_column(copy)
        else:
            return self._nonzero_positions_by_row(copy)

    def _nonzero_positions_by_row(self, copy=True):
        x = self.fetch('nonzero_positions')
        if not x is None:
            if copy:
                return list(x)
            return x
        cdef Py_ssize_t i, j
        z = self._base_ring(0)
        nzp = []
        for i from 0 <= i < self._nrows:
           for j from 0 <= j < self._ncols:
                if self.get_unsafe(i,j) != z:
                    nzp.append((i,j))
        self.cache('nonzero_positions', nzp)
        if copy:
            return list(nzp)
        return nzp

    def _nonzero_positions_by_column(self, copy=True):
        """
        Returns the list of pairs (i,j) such that self[i,j] != 0, but sorted by columns, i.e.,
        column j=0 entries occur first, then column j=1 entries, etc.

        It is safe to change the resulting list (unless you give the option copy=False).
        """
        x = self.fetch('nonzero_positions_by_column')
        if not x is None:
            if copy:
                return list(x)
            return x
        cdef Py_ssize_t i, j
        z = self._base_ring(0)
        nzp = []
        for j from 0 <= j < self._ncols:
            for i from 0 <= i < self._nrows:
                if self.get_unsafe(i,j) != z:
                    nzp.append((i,j))
        self.cache('nonzero_positions_by_column', nzp)
        if copy:
            return list(nzp)
        return nzp

    def nonzero_positions_in_column(self, Py_ssize_t i):
        """
        Return a sorted list of the integers j such that self[j,i] is
        nonzero, i.e., such that the j-th position of the i-th column
        is nonzero.

        INPUT:
            i -- an integer

        OUTPUT:
            list

        EXAMPLES:
            sage: a = matrix(QQ, 3,2, [1,2,0,2,0,0]); a
            [1 2]
            [0 2]
            [0 0]
            sage: a.nonzero_positions_in_column(0)
            [0]
            sage: a.nonzero_positions_in_column(1)
            [0, 1]

        You'll get an IndexError, if you select an invalid column:
            sage: a.nonzero_positions_in_column(2)
            Traceback (most recent call last):
            ...
            IndexError: matrix column index out of range
        """
        cdef Py_ssize_t j
        z = self._base_ring(0)
        tmp = []

        if i<0 or i >= self._ncols:
            raise IndexError, "matrix column index out of range"
        for j from 0 <= j < self._nrows:
            if self.get_unsafe(j,i) != z:
                tmp.append(j)
        return tmp

    def nonzero_positions_in_row(self, Py_ssize_t i):
        """
        Return the integers j such that self[i,j] is nonzero, i.e.,
        such that the j-th position of the i-th row is nonzero.

        INPUT:
            i -- an integer

        OUTPUT:
            list

        EXAMPLES:
            sage: a = matrix(QQ, 3,2, [1,2,0,2,0,0]); a
            [1 2]
            [0 2]
            [0 0]
            sage: a.nonzero_positions_in_row(0)
            [0, 1]
            sage: a.nonzero_positions_in_row(1)
            [1]
            sage: a.nonzero_positions_in_row(2)
            []
        """
        cdef Py_ssize_t j

        if i<0 or i >= self._nrows:
            raise IndexError, "matrix row index out of range"

        z = self._base_ring(0)
        tmp = []

        for j from 0 <= j < self._ncols:
            if self.get_unsafe(i,j) != z:
                tmp.append(j)
        return tmp

    ###################################################
    # Arithmetic
    ###################################################
    cdef Vector _vector_times_matrix_c_impl(self, Vector v):
        """
        Returns the vector times matrix product.

        INPUT:
             v -- a free module element.

        OUTPUT:
            The the vector times matrix product v*A.

        EXAMPLES:
            sage: B = matrix(QQ,2, [1,2,3,4])
            sage: V = VectorSpace(QQ, 2)
            sage: v = V([-1,5])
            sage: v*B
            (14, 18)
            sage: -1*B.row(0) + 5*B.row(1)
            (14, 18)
            sage: B*v    # computes B*v, where v is interpreted as a column vector.
            (9, 17)
            sage: -1*B.column(0) + 5*B.column(1)
            (9, 17)
        """
        M = sage.modules.free_module.FreeModule(self._base_ring, self.ncols(), sparse=self.is_sparse())
        if not PY_TYPE_CHECK(v, sage.modules.free_module_element.FreeModuleElement):
            v = M(v)
        if self.nrows() != v.degree():
            raise ArithmeticError, "number of rows of matrix must equal degree of vector"
        s = M(0)
        for i in xrange(self.nrows()):
            if v[i] != 0:
                s = s + v[i]*self.row(i)
        return s

    cdef Vector _matrix_times_vector_c_impl(self, Vector v):
        M = sage.modules.free_module.FreeModule(self._base_ring, self.nrows(), sparse=self.is_sparse())
        if not PY_TYPE_CHECK(v, sage.modules.free_module_element.FreeModuleElement):
            v = M(v)
        if self.ncols() != v.degree():
            raise ArithmeticError, "number of columns of matrix must equal degree of vector"
        s = M(0)
        for i in xrange(self.ncols()):
            if v[i] != 0:
                s = s + self.column(i, from_list=True)*v[i]
        return s

    def iterates(self, v, n, rows=True):
        r"""
        Let $A$ be this matrix and $v$ be a free module element.  If
        rows is True, return a matrix whose rows are the entries of
        the following vectors:
        $$
          v, v A, v A^2, \ldots, v A^{n-1}.
        $$

        If rows is False, return a matrix whose columns are the entries
        of the following vectors:
        $$
          v, Av, A^2 v, \ldots, A^{n-1} v.
        $$

        INPUT:
            v -- free module element
            n -- nonnegative integer

        EXAMPLES:
            sage: A = matrix(ZZ,2, [1,1,3,5]); A
            [1 1]
            [3 5]
            sage: v = vector([1,0])
            sage: A.iterates(v,0)
            []
            sage: A.iterates(v,5)
            [  1   0]
            [  1   1]
            [  4   6]
            [ 22  34]
            [124 192]

        Another example:
            sage: a = matrix(ZZ,3,range(9)); a
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: v = vector([1,0,0])
            sage: a.iterates(v,4)
            [  1   0   0]
            [  0   1   2]
            [ 15  18  21]
            [180 234 288]
            sage: a.iterates(v,4,rows=False)
            [  1   0  15 180]
            [  0   3  42 558]
            [  0   6  69 936]
        """
        n = int(n)
        if n >= 2 and self.nrows() != self.ncols():
            raise ArithmeticError, "matrix must be square if n >= 2."
        if n == 0:
            return self.matrix_space(n, self.ncols())(0)
        m = self.nrows()
        M = sage.modules.free_module.FreeModule(self._base_ring, m, sparse=self.is_sparse())
        v = M(v)
        X = [v]

        if rows:
            for _ in range(n-1):
                X.append(X[len(X)-1]*self)
            MS = self.matrix_space(n, m)
            return MS(X)
        else:
            for _ in range(n-1):
                X.append(self*X[len(X)-1])
            MS = self.matrix_space(n, m)
            return MS(X).transpose()

    cdef ModuleElement _add_c_impl(self, ModuleElement right):
        """
        Add two matrices with the same parent.
        """
        cdef Py_ssize_t i, j
        cdef Matrix A
        A = self.new_matrix()
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                A.set_unsafe(i,j, self.get_unsafe(i,j) + (<Matrix>right).get_unsafe(i,j))
        return A

    cdef ModuleElement _sub_c_impl(self, ModuleElement right):
        """
        Subtract two matrices with the same parent.

        EXAMPLES:
            sage: R.<x,y> = FreeAlgebra(QQ,2)
            sage: a = matrix(2,2, [1,2,x*y,y*x])
            sage: b = matrix(2,2, [1,2,y*x,y*x])
        """
        cdef Py_ssize_t i, j
        cdef Matrix A
        A = self.new_matrix()
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                A.set_unsafe(i,j, self.get_unsafe(i,j) - (<Matrix>right).get_unsafe(i,j))
        return A


    def __mod__(self, p):
        r"""
        Return matrix mod $p$, returning again a matrix over the same
        base ring.

        \note{Use \code{A.Mod(p)} to obtain a matrix over the residue
        class ring modulo $p$.}

        EXAMPLES:
            sage: M = Matrix(ZZ, 2, 2, [5, 9, 13, 15])
            sage: M % 7
            [5 2]
            [6 1]
            sage: parent(M % 7)
            Full MatrixSpace of 2 by 2 dense matrices over Integer Ring
        """
        cdef Py_ssize_t i
        v = self.list()
        for i from 0 <= i < len(v):
            v[i] = v[i] % p
        return self.new_matrix(entries = v, copy=False, coerce=True)

    def mod(self, p):
        """
        Return matrix mod $p$, over the reduced ring.

        EXAMPLES:
            sage: M = matrix(ZZ, 2, 2, [5, 9, 13, 15])
            sage: M.mod(7)
            [5 2]
            [6 1]
            sage: parent(M.mod(7))
            Full MatrixSpace of 2 by 2 dense matrices over Ring of integers modulo 7
        """
        return self.change_ring(self._base_ring.quotient_ring(p))


    cdef ModuleElement _rmul_c_impl(self, RingElement left):
        # derived classes over a commutative base *just* overload _lmul_c_impl (!!)
        """
        EXAMPLES:
            sage: a = matrix(QQ['x'],2,range(6))
            sage: (3/4) * a
            [   0  3/4  3/2]
            [ 9/4    3 15/4]

            sage: R.<x,y> = QQ[]
            sage: a = matrix(R,2,3,[1,x,y,-x*y,x+y,x-y]); a
            [       1        x        y]
            [  -1*x*y    y + x -1*y + x]
            sage: (x*y) * a
            [             x*y            x^2*y            x*y^2]
            [      -1*x^2*y^2    x*y^2 + x^2*y -1*x*y^2 + x^2*y]

            sage: R.<x,y> = FreeAlgebra(ZZ,2)
            sage: a = matrix(R,2,3,[1,x,y,-x*y,x+y,x-y]); a
            [    1     x     y]
            [ -x*y x + y x - y]
            sage: (x*y) * a
            [          x*y         x*y*x         x*y^2]
            [     -x*y*x*y x*y*x + x*y^2 x*y*x - x*y^2]
        """
        if PY_TYPE_CHECK(self._base_ring, CommutativeRing):
            return self._lmul_c_impl(left)
        cdef Py_ssize_t r,c
        cdef RingElement x
        x = self._base_ring(left)
        cdef Matrix ans
        ans = self._parent.zero_matrix()
        for r from 0 <= r < self._nrows:
            for c from 0 <= c < self._ncols:
                ans.set_unsafe(r, c, x._mul_c(<RingElement>self.get_unsafe(r, c)))
        return ans

    cdef ModuleElement _lmul_c_impl(self, RingElement right):
        # derived classes over a commutative base *just* overload this and not _rmul_c_impl
        """
        EXAMPLES:
        A simple example in which the base ring is commutative:
            sage: a = matrix(QQ['x'],2,range(6))
            sage: a*(3/4)
            [   0  3/4  3/2]
            [ 9/4    3 15/4]

        An example in which the base ring is not commutative:
            sage: F.<x,y> = FreeAlgebra(QQ,2)
            sage: a = matrix(2,[x,y,x^2,y^2]); a
            [  x   y]
            [x^2 y^2]
            sage: x * a
            [  x^2   x*y]
            [  x^3 x*y^2]
            sage: a * y
            [  x*y   y^2]
            [x^2*y   y^3]

            sage: R.<x,y> = FreeAlgebra(ZZ,2)
            sage: a = matrix(R,2,3,[1,x,y,-x*y,x+y,x-y]); a
            [    1     x     y]
            [ -x*y x + y x - y]
            sage: a * (x*y)
            [          x*y         x^2*y         y*x*y]
            [     -x*y*x*y x^2*y + y*x*y x^2*y - y*x*y]
        """
        cdef Py_ssize_t r,c
        cdef RingElement x
        x = self._base_ring(right)
        cdef Matrix ans
        ans = self._parent.zero_matrix()
        for r from 0 <= r < self._nrows:
            for c from 0 <= c < self._ncols:
                ans.set_unsafe(r, c, (<RingElement>self.get_unsafe(r, c))._mul_c(x))
        return ans

    cdef sage.structure.element.Matrix _matrix_times_matrix_c_impl(self, sage.structure.element.Matrix right):
        r"""
        Return the product of two matrices, a vector*matrix and
        matrix*vector product, or a scalar*matrix or matrix*scalar
        product.

        EXAMPLE of matrix times matrix over same base ring:
        We multiply matrices over $\QQ[x,y]$.
            sage: R.<x,y> = QQ[]
            sage: a = matrix(R,2,3,[1,x,y,-x*y,x+y,x-y]); a
            [       1        x        y]
            [  -1*x*y    y + x -1*y + x]
            sage: b = a.transpose(); b
            [       1   -1*x*y]
            [       x    y + x]
            [       y -1*y + x]
            sage: a*b
            [          1 + y^2 + x^2      -1*y^2 + x*y + x^2]
            [     -1*y^2 + x*y + x^2 2*y^2 + 2*x^2 + x^2*y^2]
            sage: b*a
            [        1 + x^2*y^2   x - x*y^2 - x^2*y   y + x*y^2 - x^2*y]
            [  x - x*y^2 - x^2*y y^2 + 2*x*y + 2*x^2  -1*y^2 + x*y + x^2]
            [  y + x*y^2 - x^2*y  -1*y^2 + x*y + x^2 2*y^2 - 2*x*y + x^2]

        We verify that the matrix multiplies are correct by comparing them
        with what PARI gets:
            sage: gp(a)*gp(b) - gp(a*b)
            [0, 0; 0, 0]
            sage: gp(b)*gp(a) - gp(b*a)
            [0, 0, 0; 0, 0, 0; 0, 0, 0]

        EXAMPLE of matrix times matrix over different base rings:
            sage: a = matrix(ZZ,2,2,range(4))
            sage: b = matrix(GF(7),2,2,range(4))
            sage: c = matrix(QQ,2,2,range(4))
            sage: d = a*b; d
            [2 3]
            [6 4]
            sage: parent(d)
            Full MatrixSpace of 2 by 2 dense matrices over Finite Field of size 7
            sage: parent(b*a)
            Full MatrixSpace of 2 by 2 dense matrices over Finite Field of size 7
            sage: d = a*c; d
            [ 2  3]
            [ 6 11]
            sage: parent(d)
            Full MatrixSpace of 2 by 2 dense matrices over Rational Field
            sage: d = b*c; d
            Traceback (most recent call last):
            ...
            TypeError: Base rings must be the same.
            sage: d = b*c.change_ring(GF(7)); d
            [2 3]
            [6 4]


        EXAMPLE of matrix times matrix where one matrix is sparse and
        the other is dense (in such mixed cases, the result is always
        dense):
            sage: a = matrix(ZZ,2,2,range(4),sparse=True)
            sage: b = matrix(GF(7),2,2,range(4),sparse=False)
            sage: c = a*b; c
            [2 3]
            [6 4]
            sage: parent(c)
            Full MatrixSpace of 2 by 2 dense matrices over Finite Field of size 7
            sage: c = b*a; c
            [2 3]
            [6 4]
            sage: parent(c)
            Full MatrixSpace of 2 by 2 dense matrices over Finite Field of size 7

        EXAMPLE of matrix multiplication over a noncommutative base ring:
            sage: R.<x,y> = FreeAlgebra(QQ,2)
            sage: x*y - y*x
            x*y - y*x
            sage: a = matrix(2,2, [1,2,x,y])
            sage: b = matrix(2,2, [x,y,x^2,y^2])
            sage: a*b
            [  x + 2*x^2   y + 2*y^2]
            [x^2 + y*x^2   x*y + y^3]
            sage: b*a
            [    x + y*x   2*x + y^2]
            [x^2 + y^2*x 2*x^2 + y^3]

        EXAMPLE of row vector times matrix (vectors are row vectors, so
        matrices act from the right):
            sage: a = matrix(2,3, range(6)); a
            [0 1 2]
            [3 4 5]
            sage: V = ZZ^2
            sage: v = V([-2,3]); v
            (-2, 3)
            sage: v*a
            (9, 10, 11)

        This is not allowed, since v is a {\em row} vector:
            sage: a*v
            Traceback (most recent call last):
            ...
            TypeError: incompatible dimensions

        This illustrates how coercion works:
            sage: V = QQ^2
            sage: v = V([-2,3]); v
            (-2, 3)
            sage: parent(v*a)
            Vector space of dimension 3 over Rational Field

        EXAMPLE of matrix times column vector:
          (column vectors are not implemented yet)  TODO TODO

        EXAMPLE of scalar times matrix:
            sage: a = matrix(2,3, range(6)); a
            [0 1 2]
            [3 4 5]
            sage: b = 3*a; b
            [ 0  3  6]
            [ 9 12 15]
            sage: parent(b)
            Full MatrixSpace of 2 by 3 dense matrices over Integer Ring
            sage: b = (2/3)*a; b
            [   0  2/3  4/3]
            [   2  8/3 10/3]
            sage: parent(b)
            Full MatrixSpace of 2 by 3 dense matrices over Rational Field

        EXAMPLE of matrix times scalar:
            sage: a = matrix(2,3, range(6)); a
            [0 1 2]
            [3 4 5]
            sage: b = a*3; b
            [ 0  3  6]
            [ 9 12 15]
            sage: parent(b)
            Full MatrixSpace of 2 by 3 dense matrices over Integer Ring
            sage: b = a*(2/3); b
            [   0  2/3  4/3]
            [   2  8/3 10/3]
            sage: parent(b)
            Full MatrixSpace of 2 by 3 dense matrices over Rational Field

        EXAMPLE of scalar multiplication in the noncommutative case:
            sage: R.<x,y> = FreeAlgebra(ZZ,2)
            sage: a = matrix(2,[x,y,x^2,y^2])
            sage: a * x
            [  x^2   y*x]
            [  x^3 y^2*x]
            sage: x * a
            [  x^2   x*y]
            [  x^3 x*y^2]
            sage: a*x - x*a
            [             0     -x*y + y*x]
            [             0 -x*y^2 + y^2*x]
        """
        # Both self and right are matrices with compatible dimensions and base ring.
        if self._will_use_strassen(right):
            return self._multiply_strassen(right)
        else:
            return self._multiply_classical(right)

    cdef int _will_use_strassen(self, Matrix right) except -2:
        """
        Whether or not matrix multiplication of self by right should
        be done using Strassen.

        Overload _strassen_default_cutoff to return -1 to not use
        Strassen by default.
        """
        cdef int n
        n = self._strassen_default_cutoff(right)
        if n == -1:
            return 0  # do not use Strassen
        if self._nrows > n and self._ncols > n and \
               right._nrows > n and right._ncols > n:
            return 1
        return 0

    cdef int _will_use_strassen_echelon(self) except -2:
        """
        Whether or not matrix multiplication of self by right should
        be done using Strassen.

        Overload this in derived classes to not use Strassen by default.
        """
        cdef int n
        n = self._strassen_default_echelon_cutoff()
        if n == -1:
            return 0  # do not use Strassen
        if self._nrows > n and self._ncols > n:
            return 1
        return 0

    def __neg__(self):
        """
        Return the negative of self.

        EXAMPLES:
            sage: a = matrix(ZZ,2,range(4))
            sage: a.__neg__()
            [ 0 -1]
            [-2 -3]
            sage: -a
            [ 0 -1]
            [-2 -3]
        """
        return self._lmul_c_impl(self._base_ring(-1))

    def __invert__(self):
        r"""
        Return this inverse of this matrix, as a matrix over the
        fraction field.

        Raises a \code{ZeroDivisionError} if the matrix has zero
        determinant, and raises an \code{ArithmeticError}, if the
        inverse doesn't exist because the matrix is nonsquare.  Also,
        note, e.g., that the inverse of a matrix over $\ZZ$ is always
        a matrix defined over $\Q$ (even if the entries are integers).

        EXAMPLES:
            sage: A = MatrixSpace(ZZ, 2)([1,1,3,5])
            sage: ~A
            [ 5/2 -1/2]
            [-3/2  1/2]
            sage: A.__invert__()
            [ 5/2 -1/2]
            [-3/2  1/2]

        Even if the inverse lies in the base field, the result is still a matrix
        over the fraction field.
            sage: I = MatrixSpace(ZZ,2)(1)  # identity matrix
            sage: ~I
            [1 0]
            [0 1]
            sage: (~I).parent()
            Full MatrixSpace of 2 by 2 dense matrices over Rational Field

        This is analogous to the situation for ring elements, e.g., for $\ZZ$ we have:
            sage: parent(~1)
            Rational Field
        """
        if not self.base_ring().is_field():
            return ~self.matrix_over_field()
        if not self.is_square():
            raise ArithmeticError, "self must be a square matrix"
        A = self.augment(self.parent().identity_matrix())
        B = A.echelon_form()
        if B[self._nrows - 1,  self._ncols - 1] != 1:
            raise ZeroDivisionError, "self is not invertible"
        return B.matrix_from_columns(range(self._ncols, 2*self._ncols))

    def __pos__(self):
        """
        Return +self, which is just self, of course.

        EXAMPLES:
            sage: a = matrix(ZZ,2,range(4))
            sage: +a
            [0 1]
            [2 3]
            sage: a.__pos__()
            [0 1]
            [2 3]
        """
        return self

    def __pow__(self, n, ignored):
        """
        EXAMPLES:
            sage: MS = MatrixSpace(QQ, 3, 3)
            sage: A = MS([0, 0, 1, 1, 0, '-2/11', 0, 1, '-3/11'])
            sage: A * A^(-1) == 1
            True
            sage: A^4
            [      -3/11     -13/121   1436/1331]
            [    127/121   -337/1331 -4445/14641]
            [    -13/121   1436/1331 -8015/14641]
            sage: A.__pow__(4)
            [      -3/11     -13/121   1436/1331]
            [    127/121   -337/1331 -4445/14641]
            [    -13/121   1436/1331 -8015/14641]
        """
        if not self.is_square():
            raise ArithmeticError, "self must be square"

        return RingElement.__pow__(self, n, ignored)

    ###################################################
    # Comparison
    ###################################################
    cdef long _hash(self) except -1:
        raise NotImplementedError

    cdef int _cmp_c_impl(left,Element right) except -2:
        """
        Compare two matrices.

        Matrices are compared in lexicographic order on the underlying
        list of coefficients.   A dense matrix and a sparse matrix
        are equal if their coefficients are the same.

        EXAMPLES:
        EXAMPLE cmparing sparse and dense matrices:
            sage: matrix(QQ,2,range(4)) == matrix(QQ,2,range(4),sparse=True)
            True
            sage: matrix(QQ,2,range(4)) == matrix(QQ,2,range(4),sparse=True)
            True

        Dictionary order:
            sage: matrix(ZZ,2,[1,2,3,4]) < matrix(ZZ,2,[3,2,3,4])
            True
            sage: matrix(ZZ,2,[1,2,3,4]) > matrix(ZZ,2,[3,2,3,4])
            False
            sage: matrix(ZZ,2,[0,2,3,4]) < matrix(ZZ,2,[0,3,3,4], sparse=True)
            True
        """
        raise NotImplementedError  # this is defined in the derived classes

    def __nonzero__(self):
        """
        EXAMPLES:
            sage: M = Matrix(ZZ, 2, 2, [0, 0, 0, 0])
            sage: bool(M)
            False
            sage: M = Matrix(ZZ, 2, 2, [1, 2, 3, 5])
            sage: bool(M)
            True
        """
        z = self._base_ring(0)
        cdef Py_ssize_t i, j
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                if self.get_unsafe(i,j) != z:
                    return True
        return False

    cdef int _strassen_default_cutoff(self, Matrix right) except -2:
        return -1

    cdef int _strassen_default_echelon_cutoff(self) except -2:
        return -1




#######################
# Unpickling
#######################

def unpickle(cls, parent, mutability, cache, data, version):
    r"""
    Unpickle a matrix.  This is only used internally by SAGE.
    Users should never call this function directly.

    EXAMPLES:
    We illustrating saving and loading several different types of matrices.

    OVER $\ZZ$:
        sage: A = matrix(ZZ,2,range(4))
        sage: loads(dumps(A))
        [0 1]
        [2 3]

    Sparse OVER $\QQ$:

    Dense over $\QQ[x,y]$:

    Dense over finite field.

    """
    cdef Matrix A
    A = cls.__new__(cls, parent, 0,0,0)
    A._parent = parent  # make sure -- __new__ doesn't have to set it, but unpickle may need to know.
    A._nrows = parent.nrows()
    A._ncols = parent.ncols()
    A._mutability = mutability
    A._base_ring = parent.base_ring()
    A._cache = cache
    if version >= 0:
        A._unpickle(data, version)
    else:
        A._unpickle_generic(data, version)
    return A


max_rows = 20
max_cols = 50

def set_max_rows(n):
    global max_rows
    max_rows = n

def set_max_cols(n):
    global max_cols
    max_cols = n
