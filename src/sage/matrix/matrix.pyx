"""
Abstract base class for matrices.

For design documentation see matrix/docs.py.
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
from   sage.rings.number_field.all import is_NumberField

cimport sage.structure.element
from   sage.structure.element    cimport ModuleElement, Element, RingElement, Vector
from   sage.structure.mutability cimport Mutability

import sage.modules.free_module

import matrix_misc

def is_Matrix(x):
    """
    EXAMPLES:
        sage: is_Matrix(0)
        False
        sage: is_Matrix(matrix([[1,2],[3,4]]))
        True
    """
    return isinstance(x, Matrix)

cdef class MatrixWindow  # forward declare

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
        <type 'sage.matrix.matrix_field.Matrix_field'>
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
            sage: import sage.matrix.matrix
            sage: A = sage.matrix.matrix.Matrix(MatrixSpace(QQ,2))
            sage: type(A)
            <type 'matrix.Matrix'>
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
            error

            sage: b
            [ x + 10     2/3]
            [1/2*x^2 x^3 + 1]
        """
        return self.__copy__()

    def list(self):
        """
        List of elements of self.  It is safe to change the returned list.

        EXAMPLES:
            sage: R.<x,y> = QQ[]
            sage: a = matrix(R,2,[x,y,x*y, y,x,2*x+y]); a
            sage: v = a.list(); v

        Notice that changing the returned list does not change a (the
        list is a copy):
            sage: v[0] = 25
            sage: a
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

        But the matrix doesn't know the entry changed, so it returns
        the cached version of its print representation:
            sage: a
            [0 1 2]
            [3 4 5]

        If we change an entry, the cache is cleared, and the correct print
        representation appears:
            sage: a[1,2]=10
            sage: a
            [-2/3    1    2]
            [   8    4   10]
        """
        cdef Py_ssize_t i, j

        x = self.fetch('list')
        if not x is None:
            return x
        x = []
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                x.append(self.get_unsafe(i, j))
        self.cache('list', x)
        return x

    def dict(self):
        """
        Dictionary of the elements of self with keys pairs (i,j) and
        values the nonzero entries of self.

        It is safe to change the returned dictionary.

        EXAMPLES:
            sage: R.<x,y> = QQ[]
            sage: a = matrix(R,2,[x,y,0, 0,0,2*x+y]); a
            sage: d = a.dict(); d

        Notice that changing the returned list does not change a (the
        list is a copy):
            sage: d[0,0] = 25
            sage: a
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
    cdef clear_cache(self):
        """
        Clear the properties cache.
        """
        self._cache = {}

    cdef fetch(self, key):
        """
        Try to get an element from the cache; if there isn't anything there, return None.
        """
        try:
            return self._cache[key]
        except KeyError:
            return None

    cdef cache(self, key, x):
        """
        Record x in the cache with given key.
        """
        self._cache[key] = x

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
            boom
            sage: v = {A:1}
            boom

        If we make A immutable it suddenly is hasheable.
            sage: A.set_immutable()
            sage: A.is_mutable()
            False
            sage: A[0,0] = 10
            Traceback (most recent call last):
            ...
            <type 'exceptions.ValueError'>: matrix is immutable; please change a copy instead (use self.copy())
            sage: hash(A)
            ...
            sage: v = {A:1}; v
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
            IndexError: index out of bounds
            sage: a[-1,0]
            Traceback (most recent call last):
            ...
            IndexError: index out of bounds
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
            ValueError: object is immutable; please change a copy instead.

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
            <type 'exceptions.TypeError'>: Unable to coerce rational (=1/2) to an Integer.
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
        r"""
        EXAMPLES:
            sage: R = PolynomialRing(QQ,6,'z')
            sage: a = matrix(2,3, R.gens())
            sage: a.__repr__()
            '[z0 z1 z2]\n[z3 z4 z5]'
        """
        x = self.fetch('repr')
        if not x is None: return x
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
        self.cache('repr',s)
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
        """
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
    # Coercion to Various Systems
    ###################################################

    def _pari_init_(self):
        """
        EXAMPLES:
            sage: R.<x> = QQ['x']
            sage: a = matrix(R,2,[x+1,2/3,  x^2/2, 1+x^3]); a
            [  x + 1     2/3]
            [1/2*x^2 x^3 + 1]
            sage: b = pari(a); b
            [x + 1, 2/3; 1/2*x^2, x^3 + 1]
            sage: a.determinant()
            sage: b.determinant()

        """
        w = self.list()
        cdef Py_ssize_t nr, nc, i, j
        nr = self._nrows
        nc = self._ncols
        v = []
        for i from 0 <= i < nr:
            tmp = []
            for j from 0 <= j < nc:
                tmp.append(w[i*nc + j]._pari_init_())
            v.append( ','.join(tmp))
        return 'Mat([%s])'%(';'.join(v))

    def _gap_init_(self):
        """
        EXAMPLES:
            sage: A = MatrixSpace(QQ,3)([1,2,3,4/3,5/3,6/4,7,8,9])
            sage: g = gap(A); g
            [ [ 1, 2, 3 ], [ 4/3, 5/3, 3/2 ], [ 7, 8, 9 ] ]
            sage: g.IsMatrix()
            true
        """
        cdef Py_ssize_t i, j
        v = []
        for i from 0 <= i < self._nrows:
            tmp = []
            for j from 0 <= j < self._ncols:
                tmp.append(self.get_unsafe(i,j)._gap_init_())
            v.append( '[%s]'%(','.join(tmp)) )
        return '[%s]'%(','.join(v))

    def _mathematica_init_(self):
       """
       EXAMPLES:
           sage: A = MatrixSpace(QQ,3)([1,2,3,4/3,5/3,6/4,7,8,9])
           sage: g = mathematica(A); g
       """
       cdef Py_ssize_t i, j
       v = []
       for i from 0 <= i < self._nrows:
           w = []
           for j from 0 <= j < self._ncols:
               w.append('{%s}'%self.get_unsafe(i, j))
           v.append(','.join(w))
       return '{%s}'%(','.join(v))

    def _magma_init_(self):
        r"""
        EXAMPLES:
        We first coerce a square matrix.
            sage: A = MatrixSpace(QQ,3)([1,2,3,4/3,5/3,6/4,7,8,9])
            sage: B = magma(A); B                       # optional
            [  1   2   3]
            [4/3 5/3 3/2]
            [  7   8   9]
            sage: B.Type()                              # optional
            AlgMatElt
            sage: B.Parent()                            # optional
            Full Matrix Algebra of degree 3 over Rational Field

        We coerce a non-square matrix over $\Z/8\Z$.
            sage: A = MatrixSpace(Integers(8),2,3)([-1,2,3,4,4,-2])
            sage: B = magma(A); B                       # optional
            [7 2 3]
            [4 4 6]
            sage: B.Type()                              # optional
            ModMatRngElt
            sage: B.Parent()                            # optional
            Full RMatrixSpace of 2 by 3 matrices over IntegerRing(8)
        """
        K = self._base_ring._magma_init_()
        if self._nrows != self._ncols:
            s = 'MatrixAlgebra(%s, %s)'%(K, self.nrows())
        else:
            s = 'RMatrixSpace(%s, %s, %s)'%(K, self.nrows(), self.ncols())
        v = []
        for x in self.list():
            v.append(x._magma_init_())
        return s + '![%s]'%(','.join(v))

    def _singular_(self, singular=None):
        """
        Tries to coerce this matrix to a singular matrix.
        """
        if singular is None:
            from sage.interfaces.all import singular as singular_default
            singular = singular_default
        try:
            self.base_ring()._singular_(singular)
        except (NotImplementedError, AttributeError):
            raise TypeError, "Cannot coerce to Singular"

        return singular.matrix(self.nrows(),self.ncols(),singular(self.list()))

    def numeric(self, typecode=None):
        """
        Return the Numeric array associated to this matrix (if possible).
        All entries must be coercible to the given typecode.

        INPUT:
            typecode -- optional (default: Numeric.Float64)
        """
        import Numeric
        if typecode is None:
            typecode = Numeric.Float64
        A = Numeric.array(self.list(), typecode=typecode)
        return Numeric.resize(A,(self._nrows, self._ncols))


    ###################################################
    # Construction functions
    ###################################################

    def matrix_over_field(self):
        """
        Return this matrix, but with entries viewed as elements of the
        fraction field of the base ring (assuming it is defined).

        EXAMPLES:
            sage: A = MatrixSpace(IntegerRing(),2)([1,2,3,4])
            sage: B = A.matrix_over_field()
            sage: B
            [1 2]
            [3 4]
            sage: B.parent()
            Full MatrixSpace of 2 by 2 dense matrices over Rational Field
        """
        return self.change_ring(self.base_ring().fraction_field())


    def lift(self):
        """
        EXAMPLES:
            sage: M = Matrix(ZZ/7, 2, 2, [5, 9, 13, 15]) ; M
            [5 2]
            [6 1]
            sage: M.lift()
            [5 2]
            [6 1]
            sage: parent(M.lift())
            Full MatrixSpace of 2 by 2 dense matrices over Integer Ring
        """
        return self.change_ring(self._base_ring.cover_ring())

    #############################################################################################
    # rows, columns, sparse_rows, sparse_columns, dense_rows, dense_columns, row, column
    #############################################################################################

    def columns(self, copy=True):
        """
        Return a list of the columns of self.

        INPUT:
            copy -- (default: True) if True, return a copy of the list of
                    columns, which is safe to change.

        If self is sparse, returns columns as sparse vectors, and if self
        is dense returns them as dense vectors.
        """
        x = self.fetch('columns')
        if not x is None:
            if copy: return list(x)
            return x
        if self.is_sparse():
            columns = self.sparse_columns(copy=copy)
        else:
            columns = self.dense_columns(copy=copy)
        self.cache('columns', columns)
        if copy: return list(columns)
        return columns

    def rows(self, copy=True):
        """
        Return a list of the rows of self.

        INPUT:
            copy -- (default: True) if True, return a copy of the list of rows, which is safe to change.

        If self is sparse, returns rows as sparse vectors, and if self
        is dense returns them as dense vectors.
        """
        x = self.fetch('rows')
        if not x is None:
            if copy: return list(x)
            return x
        if self.is_sparse():
            rows = self.sparse_rows(copy=copy)
        else:
            rows = self.dense_rows(copy=copy)
        self.cache('rows', rows)
        if copy: return list(rows)
        return rows

    def dense_columns(self, copy=True):
        """
        Return list of the dense columns of self.

        INPUT:
            copy -- (default: True) if True, return a copy so you can modify it safely

        EXAMPLES:
        An example over the integers:
            sage: a = matrix(3,3,range(9)); a
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: a.dense_columns()
            [(0, 3, 6), (1, 4, 7), (2, 5, 8)]

        We do an example over a polynomial ring:
            sage: R.<x> = QQ[ ]
            sage: a = matrix(R, 2, [x,x^2, 2/3*x,1+x^5]); a
            [      x     x^2]
            [  2/3*x x^5 + 1]
            sage: a.dense_columns()
            [(x, 2/3*x), (x^2, x^5 + 1)]
            sage: a = matrix(R, 2, [x,x^2, 2/3*x,1+x^5], sparse=True)
            sage: c = a.dense_columns(); c
            [(x, 2/3*x), (x^2, x^5 + 1)]
            sage: parent(c[1])
            Ambient free module of rank 2 over the principal ideal domain Univariate Polynomial Ring in x over Rational Field
        """
        x = self.fetch('dense_columns')
        if not x is None:
            if copy: return list(x)
            return x

        F = sage.modules.free_module.FreeModule(self._base_ring, self._nrows, sparse=False)
        C = []
        cdef Py_ssize_t i, j
        cdef object v
        for j from 0 <= j < self._ncols:
            v = PyList_New(self._nrows)
            for i from 0 <= i < self._nrows:
                o = self.get_unsafe(i,j)
                Py_INCREF(o)  # since we are about to set it, which doesn't increment the ref count.
                PyList_SET_ITEM(v, i, o)
            C.append(F(v, coerce=False,copy=False,check=False))
        # cache result
        self.cache('dense_columns', C)
        if copy:
            return list(C)
        else:
            return C

    def dense_rows(self, copy=True):
        """
        Return list of the dense rows of self.

        INPUT:
            copy -- (default: True) if True, return a copy so you can modify it safely

        EXAMPLES:
            sage: m = Mat(ZZ,3,3,dense=True)(range(9)); m
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: v = m.dense_rows(); v
            [(0, 1, 2), (3, 4, 5), (6, 7, 8)]
            sage: v is m.dense_rows()
            True
            sage: m[0,0] = 10
            sage: m.dense_rows()
            [(10, 1, 2), (3, 4, 5), (6, 7, 8)]
        """
        x = self.fetch('dense_rows')
        if not x is None:
            if copy: return list(x)
            return x

        F = sage.modules.free_module.FreeModule(self._base_ring, self._ncols, sparse=False)
        R = []
        cdef Py_ssize_t i, j
        cdef object o
        cdef object v
        for i from 0 <= i < self._nrows:
            v = PyList_New(self._ncols)
            for j from 0 <= j < self._ncols:
                o = self.get_unsafe(i,j)
                Py_INCREF(o)  # since we are about to set it.
                PyList_SET_ITEM(v, j, o)
            R.append(F(v, coerce=False,copy=False,check=False))
        # cache result
        self.cache('dense_rows', R)
        if copy:
            return list(R)
        else:
            return R


    def sparse_columns(self, copy=True):
        """
        Return list of the sparse columns of self.

        INPUT:
             copy -- (default: True) if True, return a copy so you can modify it safely

        EXAMPLES:
            sage: a = matrix(2,3,range(6)); a
            [0 1 2]
            [3 4 5]
            sage: v = a.sparse_columns(); v
            [(0, 3), (1, 4), (2, 5)]
            sage: v[1].is_sparse()
            True
        """
        x = self.fetch('sparse_columns')
        if not x is None:
            if copy: return list(x)
            return x

        F = sage.modules.free_module.FreeModule(self._base_ring, self._nrows, sparse=True)

        C = []
        k = 0
        entries = {}
        cdef Py_ssize_t i, j

        for i, j in self.nonzero_positions(copy=False, column_order=True):
            if j > k:
                # new column -- emit vector
                while len(C) < k:
                    C.append(F(0))
                C.append(F(entries, coerce=False, copy=False, check=False))
                entries = {}
                k = j
            entries[i] = self.get_unsafe(i, j)

        # finish up
        while len(C) < k:
            C.append(F(0))
        C.append(F(entries, coerce=False, copy=False, check=False))
        while len(C) < self._ncols:
            C.append(F(0))

        # cache result
        self.cache('sparse_columns', C)
        if copy:
            return list(C)
        else:
            return C

    def sparse_rows(self, copy=True):
        """
        Return list of the sparse rows of self.

        INPUT:
            copy -- (default: True) if True, return a copy so you can modify it safely

        EXAMPLES:
            sage: m = Mat(ZZ,3,3,sparse=True)(range(9)); m
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: v = m.sparse_rows(); v
            [(0, 1, 2), (3, 4, 5), (6, 7, 8)]
            sage: v is m.sparse_rows()
            True
            sage: v[1].is_sparse()
            True
            sage: m[0,0] = 10
            sage: m.sparse_rows()
            [(10, 1, 2), (3, 4, 5), (6, 7, 8)]
        """
        x = self.fetch('sparse_rows')
        if not x is None:
            if copy: return list(x)
            return x

        F = sage.modules.free_module.FreeModule(self._base_ring, self._ncols, sparse=True)

        R = []
        k = 0
        entries = {}
        cdef Py_ssize_t i, j

        for i, j in self.nonzero_positions(copy=False):
            if i > k:
                # new row -- emit vector
                while len(R) < k:
                    R.append(F(0))
                R.append(F(entries, coerce=False, copy=False, check=False))
                entries = {}
                k = i
            entries[j] = self.get_unsafe(i, j)

        # finish up
        while len(R) < k:
            R.append(F(0))
        R.append(F(entries, coerce=False, copy=False, check=False))
        while len(R) < self._ncols:
            R.append(F(0))

        # cache result
        self.cache('sparse_rows', R)
        if copy:
            return list(R)
        else:
            return R

    def column(self, Py_ssize_t i, from_list=False):
        """
        Return the i-th column of this matrix as a vector.

        This column is a dense vector if and only if the matrix is a
        dense matrix.

        INPUT:
            i -- integer
            from_list -- bool (default: False); if true, returns the
                         ith element of self.columns(), which may be
                         faster, but requires building a list of all
                         columns the first time it is called after an
                         entry of the matrix is changed.

        EXAMPLES:
            sage: a = matrix(2,3,range(6)); a
            [0 1 2]
            [3 4 5]
            sage: a.column(1)
            (1, 4)

        If the column is negative, it wraps around, just like with list
        indexing, e.g., -1 gives the right-most column:
            sage: a.column(-1)
            (2, 5)
        """
        i = i % self._ncols
        if i < 0:
            i = i + self._ncols
        if from_list:
            return self.columns(copy=False)[i]
        cdef Py_ssize_t j
        V = sage.modules.free_module.FreeModule(self._base_ring,
                                     self._nrows, sparse=self.is_sparse())
        tmp = []
        for j from 0 <= j < self._nrows:
            tmp.append(self.get_unsafe(j,i))
        return V(tmp, coerce=False, copy=False, check=False)

    def row(self, Py_ssize_t i, from_list=False):
        """
        Return the i-th row of this matrix as a vector.

        This row is a dense vector if and only if the matrix is a
        dense matrix.

        INPUT:
            i -- integer
            from_list -- bool (default: False); if true, returns the
                         ith element of self.rows(), which may be
                         faster, but requires building a list of all
                         rows the first time it is called after an
                         entry of the matrix is changed.

        EXAMPLES:
            sage: a = matrix(2,3,range(6)); a
            [0 1 2]
            [3 4 5]
            sage: a.row(0)
            (0, 1, 2)
            sage: a.row(1)
            (3, 4, 5)
            sage: a.row(-1)  # last row
            (3, 4, 5)
        """
        i = i % self._nrows
        if i < 0:
            i = i + self._nrows
        if from_list:
            return self.rows(copy=False)[i]
        cdef Py_ssize_t j
        V = sage.modules.free_module.FreeModule(self._base_ring,
                                      self._ncols, sparse=self.is_sparse())
        tmp = []
        for j from 0 <= j < self._ncols:
            tmp.append(self.get_unsafe(i,j))
        return V(tmp, coerce=False, copy=False, check=False)


    ############################################################################################
    # Building matrices out of other matrices, rows, or columns
    ############################################################################################
    def stack(self, Matrix other):
        """
        Return the augmented matrix self on top of other:
           [ self  ]
           [ other ]

        EXAMPLES:
            sage: M = Matrix(QQ, 2, 3, range(6))
            sage: N = Matrix(QQ, 1, 3, [10,11,12])
            sage: M.stack(N)
            [ 0  1  2]
            [ 3  4  5]
            [10 11 12]
        """
        cdef Py_ssize_t r, c

        if not (self._base_ring is other.base_ring()):
            raise TypeError, "base rings must be the same"
        if self._ncols != other._ncols:
            raise TypeError, "number of columns must be the same"

        v = self.list() + other.list()
        Z = self.new_matrix(nrows = self._nrows + other._nrows, entries=v, coerce=False, copy=False)
        return Z

    def matrix_from_columns(self, columns):
        """
        Return the matrix constructed from self using columns with
        indices in the columns list.


        EXAMPLES:
            sage: M = MatrixSpace(Integers(8),3,3)
            sage: A = M(range(9)); A
            [0 1 2]
            [3 4 5]
            [6 7 0]
            sage: A.matrix_from_columns([2,1])
            [2 1]
            [5 4]
            [0 7]
        """
        if not isinstance(columns, (list, tuple)):
            raise TypeError, "columns (=%s) must be a list of integers"%columns
        cdef Matrix A
        A = self.new_matrix(ncols = len(columns))
        k = 0
        for i in columns:
            i = int(i)
            if i < 0 or i >= self.ncols():
                raise IndexError, "column %s out of range"%i
            for r in xrange(self.nrows()):
                A.set_unsafe(r,k, self.get_unsafe(r,i))
            k = k + 1
        return A

    def matrix_from_rows(self, rows):
        """
        Return the matrix constructed from self using rows with indices
        in the rows list.

        EXAMPLES:
            sage: M = MatrixSpace(Integers(8),3,3)
            sage: A = M(range(9)); A
            [0 1 2]
            [3 4 5]
            [6 7 0]
            sage: A.matrix_from_rows([2,1])
            [6 7 0]
            [3 4 5]
        """
        if not isinstance(rows, (list, tuple)):
            raise TypeError, "rows must be a list of integers"
        cdef Matrix A
        A = self.new_matrix(nrows = len(rows))
        k = 0
        for i in rows:
            i = int(i)
            if i < 0 or i >= self.nrows():
                raise IndexError, "row %s out of range"%i
            for c in xrange(self.ncols()):
                A.set_unsafe(k,c, self.get_unsafe(i,c))
            k = k + 1
        return A

    def matrix_from_rows_and_columns(self, rows, columns):
        """
        Return the matrix constructed from self from the given
        rows and columns.

        EXAMPLES:
            sage: M = MatrixSpace(Integers(8),3,3)
            sage: A = M(range(9)); A
            [0 1 2]
            [3 4 5]
            [6 7 0]
            sage: A.matrix_from_rows_and_columns([1], [0,2])
            [3 5]
            sage: A.matrix_from_rows_and_columns([1,2], [1,2])
            [4 5]
            [7 0]

        Note that row and column indices can be reordered or repeated:
            sage: A.matrix_from_rows_and_columns([2,1], [2,1])
            [0 7]
            [5 4]

        For example here we take from row 1 columns 2 then 0 twice,
        and do this 3 times.
            sage: A.matrix_from_rows_and_columns([1,1,1],[2,0,0])
            [5 3 3]
            [5 3 3]
            [5 3 3]

        AUTHOR:
            -- Jaap Spies (2006-02-18)
        """
        if not isinstance(rows, list):
            raise TypeError, "rows must be a list of integers"
        if not isinstance(columns, list):
            raise TypeError, "columns must be a list of integers"
        cdef Matrix A
        A = self.new_matrix(nrows = len(rows), ncols = len(columns))
        r = 0
        c = len(columns)
        tmp = []
        for j in columns:
            if j >= 0 and j < self.ncols():
                tmp.append(int(j))
        columns = tmp
        if c != len(columns):
            raise IndexError, "column index out of range"
        for i in rows:
            i = int(i)
            if i < 0 or i >= self.nrows():
                raise IndexError, "row %s out of range"%i
            k = 0
            for j in columns:
                A.set_unsafe(r,k, self.get_unsafe(i,j))
                k = k + 1
            r = r + 1
        return A

    ####################################################################################
    # Change of representation between dense and sparse.
    ####################################################################################

    def dense_matrix(self):
        """
        If this matrix is sparse, return a dense matrix with the same
        entries.  If this matrix is dense, return this matrix (not a
        copy).

        NOTE: The definition of"dense" and "sparse" in SAGE have
        nothing to do with the number of nonzero entries.  Sparse
        and dense are properties of the underlying representation
        of the matrix.

        EXAMPLES:
            sage: A = MatrixSpace(QQ,2, sparse=True)([1,2,0,1])
            sage: A.is_sparse()
            True
            sage: B = A.dense_matrix()
            sage: B.is_sparse()
            False
            sage: A*B
            [1 4]
            [0 1]
            sage: A.parent()
            Full MatrixSpace of 2 by 2 sparse matrices over Rational Field
            sage: B.parent()
            Full MatrixSpace of 2 by 2 dense matrices over Rational Field
            sage: (A*B).parent()
            Full MatrixSpace of 2 by 2 sparse matrices over Rational Field
            sage: (B*A).parent()
            Full MatrixSpace of 2 by 2 dense matrices over Rational Field
        """
        if self.is_dense():
            return self
        return self.new_matrix(self._nrows, self._ncols, entries = self.list(), coerce=False,
                               copy = False, sparse=False)

    def sparse_matrix(self):
        """
        If this matrix is dense, return a sparse matrix with
        the same entries.  If this matrix is sparse, return this
        matrix (not a copy).

        NOTE: The definition of "dense" and "sparse" in SAGE have
        nothing to do with the number of nonzero entries.  Sparse
        and dense are properties of the underlying representation
        of the matrix.

        EXAMPLES:
            sage: A = MatrixSpace(QQ,2, sparse=False)([1,2,0,1])
            sage: A.is_sparse()
            False
            sage: B = A.sparse_matrix()
            sage: B.is_sparse()
            True
            sage: A
            [1 2]
            [0 1]
            sage: B
            [1 2]
            [0 1]
            sage: A*B
            [1 4]
            [0 1]
            sage: A.parent()
            Full MatrixSpace of 2 by 2 dense matrices over Rational Field
            sage: B.parent()
            Full MatrixSpace of 2 by 2 sparse matrices over Rational Field
            sage: (A*B).parent()
            Full MatrixSpace of 2 by 2 dense matrices over Rational Field
            sage: (B*A).parent()
            Full MatrixSpace of 2 by 2 sparse matrices over Rational Field
        """
        if self.is_sparse():
            return self
        return self.new_matrix(self._nrows, self._ncols, entries = self.dict(), coerce=False,
                               copy = False, sparse=True)

    def matrix_space(self, nrows=None, ncols=None, sparse=None):
        if nrows is None:
            nrows = self._nrows
        if ncols is None:
            ncols = self._ncols
        if sparse is None:
            sparse = self.is_sparse()
        return self.parent().matrix_space(nrows, ncols, sparse=sparse)

    def new_matrix(self, nrows=None, ncols=None, entries=0,
                   coerce=True, copy=True, sparse=None):
        """
        Create a matrix in the parent of this space with the given
        number of rows, columns, etc.  The default parameters are
        the same as for self.

        WARNING: This function called with no arguments returns the 0
        matrix by default.
        """
        return self.matrix_space(nrows, ncols, sparse=sparse)(entries=entries,
                                             coerce=coerce, copy=copy)

    def augment(self, Matrix other):
        """
        Return the augmented matrix of the form [self | other].

        EXAMPLES:
            sage: M = MatrixSpace(QQ,2,2)
            sage: A = M([1,2, 3,4])
            sage: A
            [1 2]
            [3 4]
            sage: N = MatrixSpace(QQ,2,1)
            sage: B = N([9,8])
            sage: B
            [9]
            [8]
            sage: A.augment(B)
            [1 2 9]
            [3 4 8]
            sage: B.augment(A)
            [9 1 2]
            [8 3 4]
            sage: M = MatrixSpace(QQ,3,4)
            sage: A = M([1,2,3,4, 0,9,8,7, 2/3,3/4,4/5,9/8])
            sage: A
            [  1   2   3   4]
            [  0   9   8   7]
            [2/3 3/4 4/5 9/8]
            sage: N = MatrixSpace(QQ,3,2)
            sage: B = N([1,2, 3,4, 4,5])
            sage: B
            [1 2]
            [3 4]
            [4 5]
            sage: A.augment(B)
            [  1   2   3   4   1   2]
            [  0   9   8   7   3   4]
            [2/3 3/4 4/5 9/8   4   5]
            sage: B.augment(A)
            [  1   2   1   2   3   4]
            [  3   4   0   9   8   7]
            [  4   5 2/3 3/4 4/5 9/8]

        AUTHORS:
            -- Naqi Jaffery (2006-01-24): examples
        """
        if not (self._base_ring is other.base_ring()):
            raise TypeError, "base rings must be the same"
        if self._nrows != other._nrows:
            raise TypeError, "number of rows must be the same"

        cdef Matrix Z
        Z = self.new_matrix(ncols = self._ncols + other._ncols)

        cdef Py_ssize_t r, c
        for r from 0 <= r < self._nrows:
            for c from 0 <= c < self._ncols:
                Z.set_unsafe(r,c, self.get_unsafe(r,c))
        nc = self.ncols()

        for r from 0 <= r < other._nrows:
            for c from 0 <= c < other._ncols:
                Z.set_unsafe(r, c+nc, other.get_unsafe(r,c))

        return Z

    def block_sum(self, Matrix other):
        """
        Return the block matrix that has self and other on the diagonal:
        [self |    0  ]
        [  0  | other ]


        EXAMPLES:
            sage: A = Matrix(QQ[['t']], 2, 2, range(4))
            sage: A.block_sum(100*A)
            [  1   2   0   0]
            [  3   4   0   0]
            [  0   0 100 200]
            [  0   0 300 400]
        """
        if not isinstance(other, Matrix):
            raise TypeError, "other must be a Matrix"
        top = self.augment(self.new_matrix(ncols=other._ncols))
        bottom = other.new_matrix(ncols=self._ncols).augment(other)
        return top.stack(bottom)

    def adjoint(self):
        """
        Returns the adjoint matrix of self (matrix of cofactors).

        INPUT:
            M -- a square matrix

        OUTPUT:
            N -- the adjoint matrix, such that
                N * M = M * N = M.parent(M.det())

        ALGORITHM:
            Use PARI

        EXAMPLES:
            sage: M = Matrix(ZZ,2,2,[5,2,3,4]) ; M
            [5 2]
            [3 4]
            sage: N = M.adjoint() ; N
            [ 4 -2]
            [-3  5]
            sage: M * N
            [14  0]
            [ 0 14]
            sage: N * M
            [14  0]
            [ 0 14]
            sage: M = Matrix(QQ,2,2,[5/3,2/56,33/13,41/10]) ; M
            [  5/3  1/28]
            [33/13 41/10]
            sage: N = M.adjoint() ; N
            [ 41/10  -1/28]
            [-33/13    5/3]
            sage: M * N
            [7363/1092         0]
            [        0 7363/1092]

        BUGS:
            only implemented for matrices over ZZ or QQ
            PARI can deal with more general base rings
        """
        if self._nrows != self._ncols:
            raise ArithmeticError, "matrix must be square"
        raise NotImplementedError

    ####################################################################################
    # LLL
    ####################################################################################

    def lllgram(self):
        """
        LLL reduction of the lattice whose gram matrix is self.

        INPUT:
            M -- gram matrix of a definite quadratic form

        OUTPUT:
            U -- unimodular transformation matrix such that

                U.transpose() * M * U

            is LLL-reduced

        ALGORITHM:
            Use PARI

        EXAMPLES:
            sage: M = Matrix(ZZ, 2, 2, [-5,3,3,-2]) ; M
            [-5  3]
            [ 3 -2]
            sage: U = M.lllgram() ; U
            [1 1]
            [1 2]
            sage: U.transpose() * M * U
            [-1  0]
            [ 0 -1]
            sage: M = Matrix(QQ,2,2,[269468, -199019/2, -199019/2, 36747]) ; M
            [   269468 -199019/2]
            [-199019/2     36747]
            sage: U = M.lllgram() ; U
            [-113  -24]
            [-306  -65]
            sage: U.transpose() * M * U
            [  2 1/2]
            [1/2   3]

        Semidefinite and indefinite forms raise a ValueError:

            sage: Matrix(ZZ,2,2,[2,6,6,3]).lllgram()
            Traceback (most recent call last):
            ...
            ValueError: not a definite matrix
            sage: Matrix(ZZ,2,2,[1,0,0,-1]).lllgram()
            Traceback (most recent call last):
            ...
            ValueError: not a definite matrix

        BUGS:
            should work for semidefinite forms (PARI is ok)
        """
        if self._nrows != self._ncols:
            raise ArithmeticError, "matrix must be square"
        raise NotImplementedError


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
            [      1       x 1/4*x^2]
            [    x^3     x^4 1/4*x^5]

        We try and fail to rescale a matrix over the integers by a non-integer:
            sage: a = matrix(ZZ,2,3,[0,1,2, 3,4,4]); a
            [0 1 2]
            [3 4 4]
            sage: a.rescale_col(2,1/2); a
            Traceback (most recent call last):
            ...
            TypeError: Unable to coerce rational (=1/2) to an Integer.

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

    def set_row_to_multiple_of_row(self, i, j, s):
        """
        Set row i equal to s times row j.

        EXAMPLES:
            sage: a = matrix(QQ,2,3,range(6)); a
            [0 1 2]
            [3 4 5]
            sage: a.set_row_to_multiple_of_row(1,2,-3)
            sage: a
            [ 0  1  2]
            [ 0 -3 -6]
        """
        self[i] = s*self[j]

    def _set_row_to_negative_of_row_of_A_using_subset_of_columns(self, Py_ssize_t i, Matrix A,
                                                                 Py_ssize_t r, cols):
        """
        Set row i of self to -(row r of A), but where we only take
        the given column positions in that row of A:

        INPUT:
            i -- integer, index into the rows of self
            A -- a matrix
            r -- integer, index into rows of A
            cols -- a *sorted* list of integers.

        EXAMPLES:
            sage: a = matrix(QQ,2,3,range(6)); a
            [0 1 2]
            [3 4 5]
            sage: a._set_row_to_negative_of_row_of_A_using_subset_of_columns(0,a,1,[1,2])
            sage: a
            [-4 -5  2]
            [ 3  4  5]
        """
        # this function exists just because it is useful for modular symbols presentations.
        cdef Py_ssize_t l
        z = self._base_ring(0)
        l = 0
        for k in cols:
            self[i,l] = -A[r,k]
            l = l + 1

    ###################################################
    # Matrix-vector multiply
    ###################################################
    def linear_combination_of_rows(self, v):
        """
        Return the linear combination of the rows of self given by the coefficients in
        the list v.

        INPUT:
            v -- list

        EXAMPLES:
            sage: a = matrix(QQ,2,3,range(6)); a
            [0 1 2]
            [3 4 5]
            sage: a.linear_combination_of_rows([1,2])
            (6, 9, 12)
        """
        cdef Py_ssize_t i
        R = self.rows()
        s = 0
        for i from 0 <= i < len(v):
            s = s + v[i] * R[i]
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
        """
        cdef Py_ssize_t i

        C = self.columns()
        s = 0
        for i from 0 <= i < len(v):
            s = s + v[i]*C[i]
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
        return self._parent.is_dense()

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
        return self._parent.is_sparse()

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
            sage: A = MatrixSpace(IntegerRing(), 2)(range(4))
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
        x = self.fetch('pivots')
        if not x is None: return x
        self.echelon_form()
        return self.fetch('pivots')

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

    def prod_of_row_sums(self, cols):
        r"""
        Calculate the product of all row sums of a submatrix of $A$ for a
        list of selected columns \code{cols}.

        EXAMPLES:
            sage: a = matrix(QQ, 2,2, [1,2,3,2]); a
            [1 2]
            [3 2]
            sage: a.prod_of_row_sums([0,1])
            15

        Another example:
            sage: a = matrix(QQ, 2,3, [1,2,3,2,5,6]); a
            [1 2 3]
            [2 5 6]
            sage: a.prod_of_row_sums([1,2])
            55

        AUTHOR:
            -- Jaap Spies (2006-02-18)
        """
        cdef Py_ssize_t c, row
        pr = 1
        for row from 0 <= row < self._nrows:
            tmp = []
            for c in cols:
                if c<0 or c >= self._ncols:
                    raise IndexError, "matrix column index out of range"
                tmp.append(self.get_unsafe(row, c))
            pr = pr * sum(tmp)
        return pr

    def permanent(self):
        r"""
        Calculate and return the permanent of this $m \times n$ matrix using
        Ryser's algorithm.

        Let $A = (a_{i,j})$ be an $m \times n$ matrix over any
        commutative ring, with $m \le n$.   The permanent of $A$ is
        \[
        \text{per}(A) = \sum_\pi a_{1,\pi(1)}a_{2,\pi(2)} \cdots a_{m\,pi(m)}
        \]
        where the summation extends over all one-to-one functions $\pi$ from
        $\{1, \ldots, m\}$ to $\{1, \ldots, n\}$.

        The product $ a_{1,\pi(1)}a_{2,\pi(2)} \cdots a_{m,\pi(m)}$ is called
        diagonal product. So the permanent of an $m \times n$ matrix $A$ is the
        sum of all the diagonal products of $A$.

        Modification of theorem 7.1.1. from Brualdi and Ryser:
        Combinatorial Matrix Theory.
        Instead of deleting columns from $A$, we choose columns from $A$ and
        calculate the product of the row sums of the selected submatrix.

        INPUT:
            A -- matrix of size m x n with m <= n

        OUTPUT:
            permanent of matrix A

        EXAMPLES:
            sage: M = MatrixSpace(ZZ,4,4)
            sage: A = M([1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1])
            sage: A.permanent()
            24

            sage: M = MatrixSpace(QQ,3,6)
            sage: A = M([1,1,1,1,0,0,0,1,1,1,1,0,0,0,1,1,1,1])
            sage: A.permanent()
            36

            sage: M = MatrixSpace(RR,3,6)
            sage: A = M([1.0,1.0,1.0,1.0,0,0,0,1.0,1.0,1.0,1.0,0,0,0,1.0,1.0,1.0,1.0])
            sage: A.permanent()
            36.000000000000000

        See Sloane's sequence OEIS A079908(3) = 36, "The Dancing School Problems"

            sage: print sloane_sequence(79908)                # optional (internet connection)
            Searching Sloane's online database...
            [79908, 'Solution to the Dancing School Problem with 3 girls: f(3,n).', [1, 4, 14, 36, 76, 140, 234, 364, 536, 756, 1030, 1364, 1764, 2236, 2786, 3420, 4144, 4964, 5886, 6916, 8060, 9324, 10714, 12236, 13896, 15700, 17654, 19764, 22036, 24476, 27090, 29884, 32864, 36036, 39406, 42980, 46764, 50764, 54986, 59436]]

            sage: M = MatrixSpace(ZZ,4,5)
            sage: A = M([1,1,0,1,1,0,1,1,1,1,1,0,1,0,1,1,1,0,1,0])
            sage: A.permanent()
            32

            See Minc: Permanents, Example 2.1, p. 5.

            sage: M = MatrixSpace(QQ,2,2)
            sage: A = M([1/5,2/7,3/2,4/5])
            sage: A.permanent()
            103/175

            sage: R.<a> = PolynomialRing(IntegerRing())
            sage: A = MatrixSpace(R,2)([[a,1], [a,a+1]])
            sage: A.permanent()
            a^2 + 2*a

            sage: R.<x,y> = MPolynomialRing(IntegerRing(),2)
            sage: A = MatrixSpace(R,2)([x, y, x^2, y^2])
            sage: A.permanent()
            x0*x1^2 + x0^2*x1


        AUTHOR:
            -- Jaap Spies (2006-02-16)
                Copyright (C) 2006 Jaap Spies <j.spies@hccnet.nl>
                Copyright (C) 2006 William Stein <wstein@gmail.com>
            -- Jaap Spies (2006-02-21): added definition of permanent

        NOTES:
            -- Currently optimized for dense matrices over QQ.
        """
        cdef Py_ssize_t m, n, r

        perm = 0
        m = self._nrows
        n = self._ncols
        if not m <= n:
            raise ValueError, "must have m <= n, but m (=%s) and n (=%s)"%(m,n)

        from   sage.structure.sequence import _combinations
        from sage.rings.arith import binomial
        for r from 1 <= r < m+1:
            lst = _combinations(range(n), r)
            tmp = []
            for cols in lst:
                tmp.append(self.prod_of_row_sums(cols))
            s = sum(tmp)
            perm = perm + (-1)**(m-r) * binomial(n-r, m-r) * s
        return perm

    def permanental_minor(self, Py_ssize_t k):
        r"""
        Calculates the permanental $k$-minor of a $m \times n$ matrix.

        This is the sum of the permanents of all possible $k$ by $k$
        submatices of $A$.

        See Brualdi and Ryser: Combinatorial Matrix Theory, p. 203.
        Note the typo $p_0(A) = 0$ in that reference!  For
        applications see Theorem 7.2.1 and Theorem 7.2.4.

        Note that the permanental $m$-minor equals $per(A)$.

        For a (0,1)-matrix $A$ the permanental $k$-minor counts the
        number of different selections of $k$ 1's of $A$ with no two
        of the 1's on the same line.

        INPUT:
            self -- matrix of size m x n with m <= n

        OUTPUT:
            permanental k-minor of matrix A

        EXAMPLES:
            sage: M = MatrixSpace(ZZ,4,4)
            sage: A = M([1,0,1,0,1,0,1,0,1,0,10,10,1,0,1,1])
            sage: A.permanental_minor(2)
            114

            sage: M = MatrixSpace(ZZ,3,6)
            sage: A = M([1,1,1,1,0,0,0,1,1,1,1,0,0,0,1,1,1,1])
            sage: A.permanental_minor(0)
            1
            sage: A.permanental_minor(1)
            12
            sage: A.permanental_minor(2)
            40
            sage: A.permanental_minor(3)
            36

        Note that if k == m the permanental k-minor equals per(A)

            sage: A.permanent()
            36

            sage: A.permanental_minor(5)
            0

        For C the "complement" of A:

            sage: M = MatrixSpace(ZZ,3,6)
            sage: C = M([0,0,0,0,1,1,1,0,0,0,0,1,1,1,0,0,0,0])
            sage: m, n = 3, 6
            sage: sum([(-1)^k * C.permanental_minor(k)*factorial(n-k)/factorial(n-m) for k in range(m+1)])
            36

            See Theorem 7.2.1 of Brualdi: and Ryser: Combinatorial Matrix Theory: per(A)

        AUTHOR:
            - Jaap Spies (2006-02-19)
        """
        m = self._nrows
        n = self._ncols
        if not m <= n:
            raise ValueError, "must have m <= n, but m (=%s) and n (=%s)"%(m,n)

        R = self._base_ring
        if k == 0:
            return R(1)
        if k > m:
            return R(0)

        k = int(k)
        pm = 0
        for cols in _combinations(range(n),k):
            for rows in _combinations(range(m),k):
                pm = pm + self.matrix_from_rows_and_columns(rows, cols).permanent()
        return pm

    def rook_vector(self, check = False):
        r"""
        Returns rook vector of this matrix.

        Let $A$ be a general $m$ by $n$ (0,1)-matrix with $m \le n$.
        We identify $A$ with a chessboard where rooks can be placed on
        the fields corresponding with $a_{ij} = 1$. The number $r_k =
        p_k(A)$ (the permanental $k$-minor) counts the number of ways
        to place $k$ rooks on this board so that no two rooks can
        attack another.

        The rook vector is the list consisting of $r_0, r_1, \ldots, r_m$.

        The rook polynomial is defined by $r(x) = \sum_{k=0}^m r_k x^k$.

        INPUT:
            self -- m by n matrix with m <= n
            check -- True or False (default), optional

        OUTPUT:
            rook vector

        EXAMPLES:
            sage: M = MatrixSpace(ZZ,3,6)
            sage: A = M([1,1,1,1,0,0,0,1,1,1,1,0,0,0,1,1,1,1])
            sage: A.rook_vector()
            [1, 12, 40, 36]

            sage: R.<x> = PolynomialRing(ZZ)
            sage: rv = A.rook_vector()
            sage: rook_polynomial = sum([rv[k] * x^k for k in range(len(rv))])
            sage: rook_polynomial
            36*x^3 + 40*x^2 + 12*x + 1

        AUTHOR:
            - Jaap Spies (2006-02-24)
        """
        m = self._nrows
        n = self._ncols
        if not m <= n:
            raise ValueError, "must have m <= n, but m (=%s) and n (=%s)"%(m,n)

        if check:
            # verify that self[i, j] in {0, 1}
            for i in range(m):
                for j in range(n):
                    x = self.get_unsafe(i, j)
                    if not (x == 0 or x == 1):
                        raise ValueError, "must have zero or one, but we have (=%s)"%x

        tmp = []
        for k in range(m+1):
            tmp.append(self.permanental_minor(k))
        return tmp


    def determinant(self):
        r"""
        Return the determinant of self.

        ALGORITHM: This is computed using the very stupid expansion by
        minors stupid \emph{naive generic algorithm}.  For matrices
        over more most rings more sophisticated algorithms can be
        used.  (Type \code{A.determinant?} to see what is done for a
        specific matrix A.)

        EXAMPLES:
            sage: A = MatrixSpace(Integers(8),3)([1,7,3, 1,1,1, 3,4,5])
            sage: A.determinant()
            6
            sage: A.determinant() is A.determinant()
            False
            sage: A.set_immutable()
            sage: A.determinant() is A.determinant()
            True

        We compute the determinant of the arbitrary 3x3 matrix:
            sage: R = PolynomialRing(QQ,9,'x')
            sage: A = matrix(R,3,R.gens())
            sage: A
            [x0 x1 x2]
            [x3 x4 x5]
            [x6 x7 x8]
            sage: A.determinant()
            -1*x2*x4*x6 + x2*x3*x7 + x1*x5*x6 - x1*x3*x8 - x0*x5*x7 + x0*x4*x8

        We create a matrix over $\Z[x,y]$ and compute its determinant.
            sage: R.<x,y> = MPolynomialRing(IntegerRing(),2)
            sage: A = MatrixSpace(R,2)([x, y, x**2, y**2])
            sage: A.determinant()
            x0*x1^2 - x0^2*x1
        """
        if self._nrows != self._ncols:
            raise ValueError, "self must be square"

        d = self.fetch('det')
        if not d is None: return d

        cdef Py_ssize_t i, n

        # if charpoly known, then det is easy.
        D = self.fetch('charpoly')
        if not D is None:
            c = D.iteritems()[0][1][0]
            if self._nrows % 2:
                c = -c
            d = self._coerce_element(c)
            self.cache('det', d)
            return d

        # if over an integral domain, get the det by computing charpoly.
        if self._base_ring.is_integral_domain():
            c = self.charpoly('x')[0]
            if self._nrows % 2:
                c = -c
            d = self._coerce_element(c)
            self.cache('det', d)
            return d

        # fall back to very very stupid algorithm
        # TODO: surely there is something much better, even in total generality...
        # this is ridiculous.
        n = self._ncols
        R = self.parent().base_ring()
        if n == 0:
            d = R(1)
            self.cache('det', d)
            return d

        elif n == 1:
            d = self.get_unsafe(0,0)
            self.cache('det', d)
            return d

        d = R(0)
        s = R(1)
        A = self.matrix_from_rows(range(1, n))
        sgn = R(-1)
        for i from 0 <= i < n:
            v = range(0,i) + range(i+1,n)
            B = A.matrix_from_columns(v)
            d = d + s*self.get_unsafe(0,i) * B.determinant()
            s = s*sgn

        self.cache('det', d)
        return d


    # shortcuts
    def det(self, *args, **kwds):
        """
        Synonym for self.determinant(...).

        EXAMPLES:
            sage: A = MatrixSpace(Integers(8),3)([1,7,3, 1,1,1, 3,4,5])
            sage: A.det()
            6
        """
        return self.determinant(*args, **kwds)

    def __abs__(self):
        """
        Synonym for self.determinant(...).

        EXAMPLES:
            sage: a = matrix(QQ, 2,2, [1,2,3,4]); a
            [1 2]
            [3 4]
            sage: abs(a)
            -2
        """
        return self.determinant()

    def characteristic_polynomial(self, *args, **kwds):
        """
        Synonym for self.charpoly(...).

        EXAMPLES:
            sage: a = matrix(QQ, 2,2, [1,2,3,4]); a
            [1 2]
            [3 4]
            sage: a.characteristic_polynomial('T')
            T^2 - 5*T - 2
        """
        return self.charpoly(*args, **kwds)

    def charpoly(self, var, algorithm="hessenberg"):
        r"""
        Return the characteristic polynomial of self, as a polynomial
        over the base ring.

        ALGORITHM: Compute the Hessenberg form of the matrix and read
        off the characteristic polynomial from that.  The result is
        cached.

        INPUT:
            var -- a variable name
            algorithm -- string:
                  'hessenberg' -- default (use Hessenberg form of matrix)

        EXAMPLES:
        First a matrix over $\Z$:
            sage: A = MatrixSpace(ZZ,2)( [1,2,  3,4] )
            sage: f = A.charpoly('x')
            sage: f
            x^2 - 5*x - 2
            sage: f.parent()
            Univariate Polynomial Ring in x over Integer Ring
            sage: f(A)
            [0 0]
            [0 0]

        An example over $\Q$:
            sage: A = MatrixSpace(QQ,3)(range(9))
            sage: A.charpoly('x')
            x^3 - 12*x^2 - 18*x
            sage: A.trace()
            12
            sage: A.determinant()
            0

        We compute the characteristic polynomial of a matrix over
        the polynomial ring $\Z[a]$:
            sage: R.<a> = PolynomialRing(ZZ)
            sage: M = MatrixSpace(R,2)([a,1,  a,a+1]); M
            [    a     1]
            [    a a + 1]
            sage: f = M.charpoly('x'); f
            x^2 + (-2*a - 1)*x + a^2
            sage: f.parent()
            Univariate Polynomial Ring in x over Univariate Polynomial Ring in a over Integer Ring
            sage: M.trace()
            2*a + 1
            sage: M.determinant()
            a^2

        We compute the characteristic polynomial of a matrix over the
        multi-variate polynomial ring $\Z[x,y]$:
            sage: R.<x,y> = MPolynomialRing(ZZ,2)
            sage: A = MatrixSpace(R,2)([x, y, x^2, y^2])
            sage: f = A.charpoly('x'); f
            x^2 + (-1*x1^2 - x0)*x + x0*x1^2 - x0^2*x1

        It's a little difficult to distinguish the variables.  To fix this,
        we temporarily view the indeterminate as $Z$:
            sage: with localvars(f.parent(), 'Z'): print f
            Z^2 + (-1*x1^2 - x0)*Z + x0*x1^2 - x0^2*x1

        We could also compute f in terms of Z from the start:
            sage: A.charpoly('Z')
            Z^2 + (-1*x1^2 - x0)*Z + x0*x1^2 - x0^2*x1
        """
        D = self.fetch('charpoly')
        if not D is None:
            if D.has_key(var):
                return D[var]
        else:
            D = {}
            self.cache('charpoly',D)

        f = self._charpoly_hessenberg(var)
        D[var] = f   # this caches it.
        return f


    def fcp(self, var):
        """
        Return the factorization of the characteristic polynomial of
        self.

        INPUT:
            var -- name of variable of charpoly

        EXAMPLES:
            sage: M = MatrixSpace(QQ,3,3)
            sage: A = M([1,9,-7,4/5,4,3,6,4,3])
            sage: A.fcp('x')
            (x^3 - 8*x^2 + 209/5*x - 286)
            sage: A = M([3, 0, -2, 0, -2, 0, 0, 0, 0])
            sage: A.fcp('x')
            (x - 3) * x * (x + 2)
        """
        return self.charpoly(var).factor()

##     def minimal_polynomial(self, var, algorithm=''):
##         """
##         Synonym for self.charpoly(...).

##         EXAMPLES:
##             sage: ???
##         """
##         return self.minpoly(*args, **kwds)

##     def minpoly(self, *args, **kwds):
##         """
##         EXAMPLES:
##             sage: ???
##         """
##         raise NotImplementedError

    def denominator(self):
        r"""
        Return the least common multiple of the denominators of the
        elements of self.

        If there is no denominator function for the base field, or no
        LCM function for the denominators, raise a TypeError.

        EXAMPLES:
            sage: A = MatrixSpace(QQ,2)(['1/2', '1/3', '1/5', '1/7'])
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


    def trace(self):
        """
        Return the trace of self, which is the sum of the
        diagonal entries of self.

        INPUT:
            self -- a square matrix
        OUTPUT:
            element of the base ring of self
        EXAMPLES:
        """
        if self._nrows != self._ncols:
            raise ArithmeticError, "matrix must be square"
        R = self._base_ring
        cdef Py_ssize_t i
        cdef object s
        s = R(0)
        for i from 0 <= i < self._nrows:
            s = s + self.get_unsafe(i,i)
        return s

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
            sage: MS = MatrixSpace(QQ, 2,2)
            sage: B = MS.matrix([1,2,1,2])
            sage: V = VectorSpace(QQ, 2)
            sage: v = V([1,2])
            sage: B.vector_matrix_multiply(v)     # computes v*B
            (3, 6)
            sage: Bt = B.transpose()
            sage: Bt.vector_matrix_multiply(v)    # computes B*v
            (5, 5)
        """
        M = sage.modules.free_module.FreeModule(self._base_ring, self.ncols())
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
        M = sage.modules.free_module.FreeModule(self._base_ring, self.ncols())
        if not PY_TYPE_CHECK(v, sage.modules.free_module_element.FreeModuleElement):
            v = M(v)
        if self.nrows() != v.degree():
            raise ArithmeticError, "number of rows of matrix must equal degree of vector"
        s = M(0)
        for i in xrange(self.nrows()):
            if v[i] != 0:
                s = s + self.column(i, from_list=True)*v[i]
        return s


    def iterates(self, v, n):
        """
        Let $A$ be this matrix and $v$ be a free module element.
        Return a vector whose rows are the entries of the following
        vectors:
        $$
           v, v A, v A^2, \ldots, v A^{n-1}.
        $$

        INPUT:
            v -- free module element
            n -- nonnegative integer

        EXAMPLES:
            sage: A = MatrixSpace(IntegerRing(), 2)([1,1,3,5]); A
            [1 1]
            [3 5]
            sage: v = FreeModule(IntegerRing(), 2)([1,0])
            sage: A.iterates(v,0)
            []
            sage: A.iterates(v,5)
            [  1   0]
            [  1   1]
            [  4   6]
            [ 22  34]
            [124 192]
        """
        n = int(n)
        if n >= 2 and self.nrows() != self.ncols():
            raise ArithmeticError, "matrix must be square if n >= 2."
        if n == 0:
            return self.matrix_space(n, self.ncols())(0)
        M = sage.modules.free_module.FreeModule(self._base_ring, self.ncols())
        if not PY_TYPE_CHECK(v, sage.modules.free_module_element.FreeModuleElement):
            v = M(v)
        X = [v]
        for _ in range(n-1):
            X.append(X[len(X)-1]*self)
        MS = self.matrix_space(n, self.ncols())
        return MS(X)

    cdef ModuleElement _add_c_impl(self, ModuleElement right):
        """
        Add two matrices with the same parent.

        EXAMPLES:
            sage:
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
        return self._add_c_impl(right._left_scalar_multiply(-1))

##         cdef Py_ssize_t i, j
##         cdef Matrix A
##         A = self.new_matrix()
##         for i from 0 <= i < self._nrows:
##             for j from 0 <= j < self._ncols:
##                 A.set_unsafe(i,j, self.get_unsafe(i,j) - (<Matrix>right).get_unsafe(i,j))
##         return A

##     def _div_(self, right):
##         """
##         EXAMPLES:
##             sage: ???
##         """
##         # TODO TODO
##         raise NotImplementedError


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
            Full MatrixSpace of 2 by 2 dense matrices over Finite Field of size 7
        """
        return self.change_ring(self._base_ring.quotient_ring(p))


    def _scalar_multiply(self, x):
        """
        EXAMPLES:
            sage: a = matrix(QQ,2,range(6))
            sage: a._scalar_multiply(3/4)
            [   0  3/4  3/2]
            [ 9/4    3 15/4]
        """
        cdef Py_ssize_t i

        if not x in self._base_ring:
            x = self._base_ring(x)

        v = self.list()
        for i from 0 <= i < len(v):
            v[i] = x * v[i]
        return self.new_matrix(entries = v, copy=False, coerce=False)

    def _right_scalar_multiply(self, x):
        """
        EXAMPLES:
        A simple example in which the base ring is commutative:
            sage: a = matrix(QQ,2,range(6))
            sage: a._right_scalar_multiply(3/4)
            [   0  3/4  3/2]
            [ 9/4    3 15/4]

        An example in which the base ring is not commutative:
            sage: F.<x,y> = FreeAlgebra(QQ,2)
            sage: a = matrix(2,[x,y,x^2,y^2]); a
            [  x   y]
            [x^2 y^2]
            sage: a._scalar_multiply(x)
            [  x^2   x*y]
            [  x^3 x*y^2]
            sage: a._right_scalar_multiply(x)
            [  x^2   y*x]
            [  x^3 y^2*x]
            sage: a._right_scalar_multiply(y)
            [  x*y   y^2]
            [x^2*y   y^3]
        """
        cdef Py_ssize_t i

        if not x in self._base_ring:
            x = self._base_ring(x)

        v = self.list()
        for i from 0 <= i < len(v):
            v[i] = v[i]*x
        return self.new_matrix(entries = v, copy=False, coerce=False)

    def _left_scalar_multiply(self, x):
        """
        EXAMPLES:
            sage: R.<x,y> = QQ[]
            sage: a = matrix(R,2,3,[1,x,y,-x*y,x+y,x-y]); a
            [       1        x        y]
            [  -1*x*y    y + x -1*y + x]
            sage: a._left_scalar_multiply(x*y)
            [             x*y            x^2*y            x*y^2]
            [      -1*x^2*y^2    x*y^2 + x^2*y -1*x*y^2 + x^2*y]
            sage: R.<x,y> = FreeAlgebra(ZZ,2)
            sage: a = matrix(R,2,3,[1,x,y,-x*y,x+y,x-y]); a
            [    1     x     y]
            [ -x*y x + y x - y]
            sage: a._left_scalar_multiply(x*y)
            [          x*y         x*y*x         x*y^2]
            [     -x*y*x*y x*y*x + x*y^2 x*y*x - x*y^2]
            sage: a._right_scalar_multiply(x*y)
            [          x*y         x^2*y         y*x*y]
            [     -x*y*x*y x^2*y + y*x*y x^2*y - y*x*y]
        """
        return self._scalar_multiply(x)

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
            TypeError: base rings must be compatible
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
            TypeError: cannot multiply matrix times row vector -- instead compute row vector times matrix

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
            ???
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
            ???
            sage: parent(b)
            Full MatrixSpace of 2 by 3 dense matrices over Rational Field

        EXAMPLE of scalar multiplication in the noncommutative case:
            sage: todo
        """
        # Both self and right are matrices with compatible dimensions and base ring.
        if self._will_use_strassen(right):
            return self._multiply_strassen(right)
        else:
            return self._multiply_classical(right)

    cdef ModuleElement _rmul_c_impl(right, RingElement left):
        """
        The product left * right, where right is a matrix and left is
        guaranteed to be in the base ring.
        """
        return right._scalar_multiply(left)

    cdef ModuleElement _lmul_c_impl(left, RingElement right):
        """
        The product left * right, where left is a matrix and right is
        guaranteed to be in the base ring.
        """
        return left._right_scalar_multiply(right)

    cdef int _will_use_strassen(self, Matrix right) except -2:
        """
        Whether or not matrix multiplication of self by right should
        be done using Strassen.

        Overload this in derived classes to not use Strassen by default.
        """
        cdef int n
        n = self._strassen_default_cutoff(right)
        if n == -1:
            return 0  # do not use Strassen
        if self._nrows > n and self._ncols > n and right._nrows > n and right._ncols > n:
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
        return self._left_scalar_multiply(self._base_ring(-1))

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
        n = int(sage.rings.integer.Integer(n))    # coerce to integer so fractions give error.
        if n == 0:
            return self.parent()(1)
        if n < 0:
            return (~self)**(-n)
        ans = self.parent()(1)
        apow = self
        while n != 0:
            if n%2 != 0:
                ans = ans * apow
            n = n/2
            if n == 0:  # to not waste time doing an extra multiplication/increment
                break
            apow = apow * apow
        return ans


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
            sage:
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

    #####################################################################################
    # Generic Hessenberg Form and charpoly algorithm
    #####################################################################################
    def hessenberg_form(self):
        """
        Return Hessenberg form of self.

        If the base ring is merely an integral domain (and not a
        field), the Hessenberg form will (in general) only be defined
        over the fraction field of the base ring.

        EXAMPLES:
            sage: A = matrix(ZZ,4,[2, 1, 1, -2, 2, 2, -1, -1, -1,1,2,3,4,5,6,7])
            sage: h = A.hessenberg_form(); h
            [    2  -7/2 -19/5    -2]
            [    2   1/2 -17/5    -1]
            [    0  25/4  15/2   5/2]
            [    0     0  58/5     3]
            sage: parent(h)
            Full MatrixSpace of 4 by 4 dense matrices over Rational Field
            sage: A.hessenbergize()
            Traceback (most recent call last):
            ...
            TypeError: Hessenbergize only possible for matrices over a field
        """
        X = self.fetch('hessenberg_form')
        if not X is None:
            return X
        R = self._base_ring
        if not R.is_field():
            try:
                K = self._base_ring.fraction_field()
                H = self.change_ring(K)
                H.hessenbergize()
            except TypeError, msg:
                raise TypeError, "%s\nHessenberg form only possible for matrices over a field"%msg
        else:
            H = self.copy()
            H.hessenbergize()
        #end if
        self.cache('hessenberg_form', H)
        return H

    def hessenbergize(self):
        """
        Tranform self to Hessenberg form.

        The hessenberg form of a matrix $A$ is a matrix that is
        similar to $A$, so has the same characteristic polynomial as
        $A$, and is upper triangular except possible for entries right
        below the diagonal.

        ALGORITHM: See Henri Cohen's first book.

        EXAMPLES:
            sage: A = MatrixSpace(QQ,3)([2, 1, 1, -2, 2, 2, -1, -1, -1])
            sage: A.hessenberg_form()
            [  2 3/2   1]
            [ -2   3   2]
            [  0  -3  -2]

            sage: A = MatrixSpace(QQ,4)([2, 1, 1, -2, 2, 2, -1, -1, -1,1,2,3,4,5,6,7])
            sage: A.hessenberg_form()
            [    2  -7/2 -19/5    -2]
            [    2   1/2 -17/5    -1]
            [    0  25/4  15/2   5/2]
            [    0     0  58/5     3]
        """
        tm = verbose("Computing Hessenberg Normal Form of %sx%s matrix"%(n,n))

        if not self.is_square():
            raise TypeError, "self must be square"

        if not self._base_ring.is_field():
            raise TypeError, "Hessenbergize only possible for matrices over a field"

        self.check_mutability()

        cdef Py_ssize_t i, j, m, n, r
        zero = self._base_ring(0)
        one = self._base_ring(1)
        n = self._nrows
        for m from 1 <= m < n-1:
            # Search for a non-zero entry in column m-1
            i = -1
            for r from m+1 <= r < n:
                if self.get_unsafe(r, m-1) != zero:
                    i = r
                    break
            if i != -1:
                # Found a nonzero entry in column m-1 that is strictly below row m
                # Now set i to be the first nonzero position >= m in column m-1
                if self.get_unsafe(m,m-1) != zero:
                    i = m
                t = self.get_unsafe(i,m-1)
                t_inv = None
                if i > m:
                    self.swap_rows_c(i,m)
                    # We must do the corresponding column swap to
                    # maintain the characteristic polynomial (which is
                    # an invariant of Hessenberg form)
                    self.swap_columns_c(i,m)
                # Now the nonzero entry in position (m,m-1) is t.
                # Use t to clear the entries in column m-1 below m.
                for j from m+1 <= j < n:
                    x = self.get_unsafe(j, m-1)
                    if x != zero:
                        if t_inv is None:
                            t_inv = one / t
                        u = x * t_inv
                        self.add_multiple_of_row_c(j, m, -u, 0)
                        # To maintain charpoly, do the corresponding column operation,
                        # which doesn't mess up the matrix, since it only changes
                        # column m, and we're only worried about column m-1 right now.
                        # Add u*column_j to column_m.
                        self.add_multiple_of_column_c(m, j, u, 0)
        verbose("Finished Hessenberg Normal Form of %sx%s matrix"%(n,n),tm)


    def _charpoly_hessenberg(self, var):
        """
        Transforms self in place to its Hessenberg form then computes
        and returns the coefficients of the characteristic polynomial of
        this matrix.

        INPUT:
            var -- name of the indeterminate of the charpoly.

        The characteristic polynomial is represented as a vector of
        ints, where the constant term of the characteristic polynomial
        is the 0th coefficient of the vector.

        EXAMPLES:
            sage: matrix(QQ,3,range(9))._charpoly_hessenberg('Z')
            Z^3 - 12*Z^2 - 18*Z
            sage: matrix(ZZ,3,range(9))._charpoly_hessenberg('Z')
            Z^3 - 12*Z^2 - 18*Z
            sage: matrix(GF(7),3,range(9))._charpoly_hessenberg('Z')
            Z^3 + 2*Z^2 + 3*Z
            sage: matrix(QQ['x'],3,range(9))._charpoly_hessenberg('Z')
            Z^3 + (-12)*Z^2 + (-18)*Z
            sage: matrix(ZZ['ZZ'],3,range(9))._charpoly_hessenberg('Z')
            Z^3 + (-12)*Z^2 + (-18)*Z
        """
        if self._nrows != self._ncols:
            raise ArithmeticError, "charpoly not defined for non-square matrix."

        # Replace self by its Hessenberg form
        cdef Matrix H
        H = self.hessenberg_form()

        # We represent the intermediate polynomials that come up in
        # the calculations as rows of an (n+1)x(n+1) matrix, since
        # we've implemented basic arithmetic with such a matrix.
        # Please see the generic implementation of charpoly in
        # matrix.py to see more clearly how the following algorithm
        # actually works.  (The implementation is clearer (but slower)
        # if one uses polynomials to represent polynomials instead of
        # using the rows of a matrix.)  Also see Cohen's first GTM,
        # Algorithm 2.2.9.

        cdef Py_ssize_t i, m, n,
        n = self._nrows

        cdef Matrix c
        c = H.new_matrix(nrows=n+1,ncols=n+1)    # the 0 matrix
        one = H._coerce_element(1)
        c.set_unsafe(0,0,one)

        for m from 1 <= m <= n:
            # Set the m-th row of c to (x - H[m-1,m-1])*c[m-1] = x*c[m-1] - H[m-1,m-1]*c[m-1]
            # We do this by hand by setting the m-th row to c[m-1]
            # shifted to the right by one.  We then add
            # -H[m-1,m-1]*c[m-1] to the resulting m-th row.
            for i from 1 <= i <= n:
                c.set_unsafe(m, i, c.get_unsafe(m-1,i-1))
            c.add_multiple_of_row_c(m, m-1, -H.get_unsafe(m-1, m-1), 0)
            t = one
            for i from 1 <= i < m:
                t = t * H.get_unsafe(m-i,m-i-1)
                # Set the m-th row of c to c[m] - t*H[m-i-1,m-1]*c[m-i-1]
                c.add_multiple_of_row_c(m, m-i-1, - t*H.get_unsafe(m-i-1,m-1), 0)

        # The answer is now the n-th row of c.
        v = PyList_New(n+1)     # this is really sort of v = []..."
        for i from 0 <= i <= n:
            # Finally, set v[i] = c[n,i]
            o = c.get_unsafe(n,i)
            Py_INCREF(o); PyList_SET_ITEM(v, i, o)

        R = self._base_ring[var]    # polynomial ring over the base ring
        return R(v)

    #####################################################################################
    # Decomposition: kernel, image, decomposition
    #####################################################################################

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

        A kernel of dimension one over $\Q$:
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
        K = self.fetch('kernel')
        if not K is None:
            return K
        R = self._base_ring

        if self._nrows == 0:    # from a 0 space
            V = sage.modules.free_module.VectorSpace(R, self._nrows)
            Z = V.zero_subspace()
            self.cache('kernel', Z)
            return Z

        elif self._ncols == 0:  # to a 0 space
            Z = sage.modules.free_module.VectorSpace(R, self._nrows)
            self.cache('kernel', Z)
            return Z

        if is_NumberField(R):
            A = self._pari_().mattranspose()
            B = A.matker()
            n = self._nrows
            V = sage.modules.free_module.VectorSpace(R, n)
            basis = eval('[V([R(x) for x in b]) for b in B]', {'V':V, 'B':B, 'R':R})
            Z = V.subspace(basis)
            self.cache('kernel', Z)
            return Z

        E = self.transpose().echelon_form(*args, **kwds)
        pivots = E.pivots()
        pivots_set = set(pivots)
        basis = []
        VS = sage.modules.free_module.VectorSpace
        V = VS(R, self.nrows())
        ONE = R(1)
        for i in xrange(self._nrows):
            if not (i in pivots_set):
                v = V(0)
                v[i] = ONE
                for r in range(len(pivots)):
                    v[pivots[r]] = -E[r,i]
                basis.append(v)
        W = V.subspace(basis)
        if W.dimension() != len(basis):
            raise RuntimeError, "bug in kernel function in matrix.pyx -- basis got from echelon form not a basis."
        self.cache('kernel', W)
        return W



    def decomposition(self, is_diagonalizable=False, dual=False):
        """
        Returns the decomposition of the free module on which this
        matrix acts from the right, along with whether this matrix
        acts irreducibly on each factor.  The factors are guaranteed
        to be sorted in the same way as the corresponding factors of
        the characteristic polynomial.

        Let A be the matrix acting from the on the vector space V of
        column vectors.  Assume that A is square.  This function
        computes maximal subspaces W_1, ..., W_n corresponding to
        Galois conjugacy classes of eigenvalues of A.  More precisely,
        let f(X) be the characteristic polynomial of A.  This function
        computes the subspace $W_i = ker(g_(A)^n)$, where g_i(X) is an
        irreducible factor of f(X) and g_i(X) exactly divides f(X).
        If the optional parameter is_diagonalizable is True, then we
        let W_i = ker(g(A)), since then we know that ker(g(A)) =
        $ker(g(A)^n)$.

        If dual is True, also returns the corresponding decomposition
        of V under the action of the transpose of A.  The factors are
        guarenteed to correspond.

        OUTPUT:
            Sequence -- list of pairs (V,t), where V is a vector spaces
                    and t is a bool, and t is True exactly when the
                    charpoly of self on V is irreducible.

            (optional) list -- list of pairs (W,t), where W is a vector
                    space and t is a bool, and t is True exactly
                    when the charpoly of the transpose of self on W
                    is irreducible.

        EXAMPLES:
            sage: MS1 = MatrixSpace(ZZ,4)
            sage: MS2 = MatrixSpace(QQ,6)
            sage: A = MS1.matrix([3,4,5,6,7,3,8,10,14,5,6,7,2,2,10,9])
            sage: B = MS2(range(36))
            sage: B*11   # random output
            [-11  22 -11 -11 -11 -11]
            [ 11 -22 -11 -22  11  11]
            [-11 -11 -11 -22 -22 -11]
            [-22  22 -22  22 -11  11]
            [ 22 -11  11 -22  11  22]
            [ 11  11  11 -22  22  22]
            sage: decomposition(A)
            [(Ambient free module of rank 4 over the principal ideal domain Integer Ring, True)]
            sage: decomposition(B)
            [(Vector space of degree 6 and dimension 4 over Rational Field
            Basis matrix:
            [ 1  0  0  0 -5  4]
            [ 0  1  0  0 -4  3]
            [ 0  0  1  0 -3  2]
            [ 0  0  0  1 -2  1],
              False),
             (Vector space of degree 6 and dimension 2 over Rational Field
            Basis matrix:
            [ 1  0 -1 -2 -3 -4]
            [ 0  1  2  3  4  5],
              True)]
        """
        if not self.is_square():
            raise ArithmeticError, "self must be a square matrix"

        if self.nrows() == 0:
            return decomp_seq([])

        f = self.charpoly('x')
        E = decomp_seq([])

        # Idea: For optimization, could compute powers of self
        #       up to max degree of any factor.  Then get g(self)
        #       by taking a linear combination.

        if dual:
            Edual = decomp_seq([])
        F = f.factor()
        if len(F) == 1:
            V = sage.modules.free_module.FreeModule(
                              self.base_ring(), self.nrows())
            m = F[0][1]
            if dual:
                return decomp_seq([(V,m==1)]), decomp_seq([(V,m==1)])
            else:
                return decomp_seq([(V,m==1)])
        F.sort()
        for g, m in f.factor():
            if is_diagonalizable:
                B = g(self)
            else:
                B = g(self) ** m
            E.append((B.kernel(), m==1))
            if dual:
                Edual.append((B.transpose().kernel(), m==1))
        if dual:
            return E, Edual
        return E

    def decomposition_of_subspace(self, M, is_diagonalizable=False):
        """
        Suppose the right action of self on M leaves M
        invariant. Return the decomposition of M as a list of pairs
        (W, is_irred) where is_irred is True if the charpoly of self
        acting on the factor W is irreducible.

        EXAMPLES:
            sage: t = matrix(QQ, 3, [3, 0, -2, 0, -2, 0, 0, 0, 0]); t
            [ 3  0 -2]
            [ 0 -2  0]
            [ 0  0  0]
            sage: t.fcp('X')   # factored charpoly
            (X - 3) * X * (X + 2)
            sage: v = kernel(t*(t+2)); v   # an invariant subspace
            Vector space of degree 3 and dimension 2 over Rational Field
            Basis matrix:
            [0 1 0]
            [0 0 1]
            sage: D = t.decomposition_of_subspace(v); D
            [
            (Vector space of degree 3 and dimension 1 over Rational Field
            Basis matrix:
            [0 0 1], 1),
            (Vector space of degree 3 and dimension 1 over Rational Field
            Basis matrix:
            [0 1 0], 1)
            ]
            sage: t.restrict(D[0][0])
            [0]
            sage: t.restrict(D[1][0])
            [-2]
        """
        if not sage.modules.free_module.is_FreeModule(M):
            raise TypeError, "M must be a free module."
        if not self.is_square():
            raise ArithmeticError, "matrix must be square"
        if M.base_ring() != self.base_ring():
            raise ArithmeticError, "base rings are incompatible"
        if M.degree() != self.ncols():
            raise ArithmeticError, \
               "M must be a subspace of an %s-dimensional space"%self.ncols()

        time = verbose(t=0)

        # 1. Restrict
        B = self.restrict(M)
        time0 = verbose("restrict -- ", time)

        # 2. Decompose restriction
        D = B.decomposition(is_diagonalizable=is_diagonalizable, dual=False)

        assert sum(eval('[A.dimension() for A,_ in D]',{'D':D})) == M.dimension(), \
               "bug in decomposition; " + \
               "the sum of the dimensions of the factors must equal the dimension of the acted on space."

        # 3. Lift decomposition to subspaces of ambient vector space.
        # Each basis vector for an element of D defines a linear combination
        # of the basis of W, and these linear combinations define the
        # corresponding subspaces of the ambient space M.

        verbose("decomposition -- ", time0)
        C = M.basis_matrix()
        Z = M.ambient_vector_space()

        D = eval('[(Z.subspace([x*C for x in W.basis()]), is_irred) for W, is_irred in D]',\
                 {'C':C, 'D':D, 'Z':Z})
        D = decomp_seq(D)

        verbose(t=time)
        return D

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
            raise IndexError, "degree of V (=%s) must equal number of rows of self (=%s)"%(\
                V.degree(), self.nrows())
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
                C = eval('[V.coordinate_vector(b*self) for b in V.basis()]',{'V':V, 'self':self})
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
        e = eval('[b*self for b in V.basis()]', {'self':self, 'V':V})
        return self.new_matrix(V.dimension(), self.ncols(), e)

    def eigenspaces(self):
        """
        Return a list of pairs
             (e, V)
        where e runs through all eigenvalues (up to Galois conjugation)
        of this matrix, and V is the corresponding eigenspace.

        WARNING: Uses a somewhat naive algorithm (simply factors the
        characteristic polynomial and computes kernels directly over
        the extension field).  TODO: Implement the better algorithm
        that is in dual_eigenvector in sage/hecke/module.py.

        EXAMPLES:
        We compute the eigenspaces of the matrix of the Hecke operator
        $T_2$ on a space:

            sage: A = ModularSymbols(43).T(2).matrix()
            sage: A.eigenspaces()
            [
            (3, [
            (1, 0, 1/7, 0, -1/7, 0, -2/7)
            ]),
            (-2, [
            (0, 1, 0, 1, -1, 1, -1),
            (0, 0, 1, 0, -1, 2, -1)
            ]),
            (a, [
            (0, 1, 0, -1, -a - 1, 1, -1),
            (0, 0, 1, 0, -1, 0, -a + 1)
            ])
            ]

        Next we compute the eigenspaces over the finite field
        of order 11:

            sage: A = ModularSymbols(43, base_ring=GF(11), sign=1).T(2).matrix()
            sage: A.eigenspaces()
            [
            (9, [
            (0, 0, 1, 5)
            ]),
            (3, [
            (1, 6, 0, 6)
            ]),
            (x, [
            (0, 1, 0, 5*x + 10)
            ])
            ]

        Finally, we compute the eigenspaces of a $3\times 3$ matrix.

            sage: A = Matrix(QQ,3,3,range(9))
            sage: A.eigenspaces()
            [
            (0, [
            (1, -2, 1)
            ]),
            (a, [
            (1, 1/15*a + 2/5, 2/15*a - 1/5)
            ])
            ]
        """
        x = self.fetch('eigenvectors')
        if not x is None:
            return x
        f = self.charpoly('x')
        G = f.factor()
        V = []
        for h, e in G:
            F = h.root_field()
            W = (self.change_ring(F) - F.gen(0)).kernel()
            V.append((F.gen(0), W.basis()))
        V = sage.structure.sequence.Sequence(V, cr=True)
        self.cache('eigenvectors', V)
        return V

    #####################################################################################
    # Generic Echelon Form
    ###################################################################################
    def echelonize(self, algorithm="default", cutoff=0):
        r"""
        Transform self into a matrix in echelon form.

        INPUT:
            algorithm -- string, which algorithm to use (default: 'default')
                   'default' -- use a default algorithm, chosen by SAGE
                   'strassen' -- use a Strassen divide and conquer algorithm
            cutoff -- integer; only used if the Strassen algorithm is selected.

        EXAMPLES:
            sage: a = matrix(QQ,3,range(9)); a
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: a.echelonize()
            sage: a
            [ 1  0 -1]
            [ 0  1  2]
            [ 0  0  0]

        An immutable matrix cannot be transformed into echelon form.
        Use \code{self.echelon_form()} instead:

            sage: a = matrix(QQ,3,range(9)); a.set_immutable()
            sage: a.echelonize()
            Traceback (most recent call last):
            ...
            ValueError: matrix is immutable; please change a copy instead (use self.copy()).
            sage: a.echelon_form()
            [ 1  0 -1]
            [ 0  1  2]
            [ 0  0  0]

        Echelon form over the integers is what is also classically often known as
        Hermite normal form:
            sage: a = matrix(ZZ,3,range(9))
            sage: a.echelonize(); a
            [ 3  0 -3]
            [ 0  1  2]
            [ 0  0  0]

        Echelon form is not defined for any integral domain; you may have to explicitly
        base extend to the fraction field, if that is what you want.
            sage: R.<x,y> = QQ[]
            sage: a = matrix(R, 2, [x,y,x,y])
            sage: a.echelonize()
            Traceback (most recent call last):
            ...
            ValueError: Echelon form not defined over this base ring.
            sage: b = a.change_ring(R.fraction_field())
            sage: b.echelonize()
            [  1 y/x]
            [  0   0]

        Echelon form is not defined over arbitrary rings:
            sage: a = matrix(Integers(9),3,range(9))
            sage: a.echelonize()
            Traceback (most recent call last):
            ...
            ValueError: Echelon form not defined over this base ring.
        """
        self.check_mutability()
        if algorithm == 'default':
            if self._will_use_strassen_echelon():
                algorithm = 'strassen'
            else:
                algorithm = 'classical'
        try:
            if algorithm == 'classical':
                self._echelon_in_place_classical()
            elif algorithm == 'strassen':
                self._echelon_strassen(cutoff)
            else:
                raise ValueError, "Unknown algorithm '%s'"%algorithm
        except ArithmeticError, msg:
            print msg
            raise ValueError, "Echelon form not defined over this base ring."

    def echelon_form(self, algorithm="default", cutoff=0):
        """
        Return the echelon form of self.

        INPUT:
            matrix -- an element A of a MatrixSpace

        OUTPUT:
            matrix -- The reduced row echelon form of A, as an
            immutable matrix.  Note that self is *not* changed by this
            command.  Use A.echelonize() to change A in place.

        EXAMPLES:
           sage: MS = MatrixSpace(QQ,2,3)
           sage: C = MS.matrix([1,2,3,4,5,6])
           sage: C.rank()
           2
           sage: C.nullity()
           1
           sage: C.echelon_form()
           [ 1  0 -1]
           [ 0  1  2]



        """
        x = self.fetch('echelon_form')
        if not x is None:
            return x
        E = self.copy()
        E.echelonize(algorithm = algorithm, cutoff=cutoff)
        E.set_immutable()  # so we can cache the echelon form.
        self.cache('echelon_form', E)
        self.cache('pivots', E.pivots())
        return E

    def _echelon_classical(self):
        """
        Return the echelon form of self.
        """
        E = self.fetch('echelon_classical')
        if not E is None:
            return E
        E = self.copy()
        E._echelon_in_place_classical()
        self.cache('echelon_classical', E)
        return E

    def _echelon_in_place_classical(self):
        """
        Return the echelon form of self and set the pivots of self.
        """
        cdef Py_ssize_t start_row, c, r, nr, nc, i
        if self.fetch('in_echelon_form'):
            return

        if not self._base_ring.is_field():
            raise ValueError, "Echelon form not defined over this base ring."

        self.check_mutability()

        start_row = 0
        nr = self._nrows
        nc = self._ncols
        pivots = []

        for c from 0 <= c < nc:
            if PyErr_CheckSignals(): raise KeyboardInterrupt
            for r from start_row <= r < nr:
                if self.get_unsafe(r, c) != 0:
                    pivots.append(c)
                    a_inverse = ~self.get_unsafe(r,c)
                    self.rescale_row(r, a_inverse, c)
                    self.swap_rows(r, start_row)
                    for i from 0 <= i < nr:
                        if i != start_row:
                            if self.get_unsafe(i,c) != 0:
                                minus_b = -self.get_unsafe(i, c)
                                self.add_multiple_of_row(i, start_row, minus_b, c)
                    start_row = start_row + 1
                    break
        self.cache('pivots', pivots)
        self.cache('in_echelon_form', True)

    #####################################################################################
    # Windowed Strassen Matrix Multiplication and Echelon
    # Precise algorithms invented and implemented by David Harvey and Robert Bradshaw
    # at William Stein's MSRI 2006 Summer Workshop on Modular Forms.
    #####################################################################################
    cdef int _strassen_default_cutoff(self, Matrix right) except -2:
        return -1

    cdef int _strassen_default_echelon_cutoff(self) except -2:
        return -1

    def _multiply_strassen(self, Matrix right, int cutoff=0):
        """
        Multiply self by the matrix right using a Strassen-based
        asymptotically fast arithmetic algorithm.

        ALGORITHM: Custom algorithm for arbitrary size matrices
        designed by David Harvey and Robert Bradshaw, based on
        Strassen's algorithm.

        INPUT:
            cutoff -- integer (default: 0 -- let class decide).

        EXAMPLES:
            sage: ?
        """
        if self._ncols != right._nrows:
            raise ArithmeticError, "Number of columns of self must equal number of rows of right."
        if not self._base_ring is right.base_ring():
            raise TypeError, "Base rings must be the same."

        if cutoff == 0:
            cutoff = self._strassen_default_cutoff(right)

        if cutoff <= 1:
            raise ValueError, "cutoff must be at least 2"

        output = self.new_matrix(self._nrows, right._ncols)

        self_window   = self.matrix_window()
        right_window  = right.matrix_window()
        output_window = output.matrix_window()

        self._strassen_window_multiply(output_window, self_window, right_window, cutoff)
        return output

    def _echelon_strassen(self, int cutoff=0):
        """
        In place Strassen echelon of self, and sets the pivots.

        ALGORITHM: Custom algorithm for arbitrary size matrices
        designed by David Harvey and Robert Bradshaw, based on
        Strassen's algorithm.

        EXAMPLES:
            sage: ?
        """
        self.check_mutability()

        if not self._base_ring.is_field():
            raise ValueError, "Echelon form not defined over this base ring."

        if cutoff == 0:
            cutoff = self._strassen_default_echelon_cutoff()

        if cutoff <= 1:
            raise ValueError, "cutoff must be at least 2"

        if self._nrows < cutoff or self._ncols < cutoff:
            self._echelon_in_place_classical()
            return

        pivots = self._strassen_echelon(self.matrix_window(), cutoff)
        self._set_pivots(pivots)


    def matrix_window(self, Py_ssize_t row=0, Py_ssize_t col=0,
                      Py_ssize_t nrows=-1, Py_ssize_t ncols=-1):
        """
        Return the requested matrix window.

        EXAMPLES:

        """
        if nrows == -1:
            nrows = self._nrows
            ncols = self._ncols
        return MatrixWindow(self, row, col, nrows, ncols)

    def _strassen_window_multiply(self, MatrixWindow C, MatrixWindow A,
                                  MatrixWindow B, int cutoff):
        """
        Multiplies the submatrices specified by A and B, places result
        in C.  Assumes that A and B have compatible dimensions to be
        multiplied, and that C is the correct size to receive the
        product, and that they are all defined over the same ring.

        Uses strassen multiplication at high levels and then uses MatrixWindow
        methods at low levels.

        todo: doc cutoff parameter as soon as I work out what it really means

        EXAMPLES:
            The following matrix dimensions are chosen especially to exercise the
            eight possible parity combinations that ocould ccur while subdividing
            the matrix in the strassen recursion. The base case in both cases will
            be a (4x5) matrix times a (5x6) matrix.

            TODO -- the doctests below are currently not
            tested/enabled/working -- enable them when linear algebra
            restructing gets going.

            sage: dim1 = 64; dim2 = 83; dim3 = 101
            sage: R = MatrixSpace(QQ, dim1, dim2)
            sage: S = MatrixSpace(QQ, dim2, dim3)
            sage: T = MatrixSpace(QQ, dim1, dim3)


            sage: A = R.random_element(range(-30, 30))
            sage: B = S.random_element(range(-30, 30))
            sage: C = T(0)
            sage: D = T(0)

            sage: A_window = A.matrix_window(0, 0, dim1, dim2)
            sage: B_window = B.matrix_window(0, 0, dim2, dim3)
            sage: C_window = C.matrix_window(0, 0, dim1, dim3)
            sage: D_window = D.matrix_window(0, 0, dim1, dim3)

            sage: strassen_window_multiply(C, A, B, 2)   # use strassen method
            sage: D_window.set_to_prod(A, B)             # use naive method
            sage: C == D
            True

            sage: dim1 = 79; dim2 = 83; dim3 = 101
            sage: R = MatrixSpace(QQ, dim1, dim2)
            sage: S = MatrixSpace(QQ, dim2, dim3)
            sage: T = MatrixSpace(QQ, dim1, dim3)

            sage: A = R.random_element(range(30))
            sage: B = S.random_element(range(30))
            sage: C = T(0)
            sage: D = T(0)

            sage: A_window = A.matrix_window(0, 0, dim1, dim2)
            sage: B_window = B.matrix_window(0, 0, dim2, dim3)
            sage: C_window = C.matrix_window(0, 0, dim1, dim3)

            sage: strassen_window_multiply(C, A, B, 2)   # use strassen method
            sage: D.set_to_prod(A, B)                    # use naive method

            sage: C == D
            True

        AUTHOR: David Harvey
        """
        # todo -- I'm not sure how to interpret "cutoff". Should it be...
        # (a) the minimum side length of the matrices (currently implemented below)
        # (b) the maximum side length of the matrices
        # (c) the total number of entries being multiplied
        # (d) something else entirely?

        cdef int A_nrows, A_ncols, B_ncols
        A_nrows = A.nrows()
        A_ncols = A.ncols()        # this should also be the number of rows of B
        B_ncols = B.ncols()

        if (A_nrows <= cutoff) or (A_ncols <= cutoff) or (B_ncols <= cutoff):
            # note: this code is only reached if the TOP level is already beneath
            # the cutoff. In a typical large multiplication, the base case is
            # handled directly (see below).
            C.set_to_prod(A, B)
            return

        # Construct windows for the four quadrants of each matrix.
        # Note that if the side lengths are odd we're ignoring the
        # final row/column for the moment.

        cdef int A_sub_nrows, A_sub_ncols, B_sub_ncols
        A_sub_nrows = A_nrows >> 1
        A_sub_ncols = A_ncols >> 1     # this is also like B_sub_nrows
        B_sub_ncols = B_ncols >> 1

        A00 = A.matrix_window(0,           0,           A_sub_nrows, A_sub_ncols)
        A01 = A.matrix_window(0,           A_sub_ncols, A_sub_nrows, A_sub_ncols)
        A10 = A.matrix_window(A_sub_nrows, 0,           A_sub_nrows, A_sub_ncols)
        A11 = A.matrix_window(A_sub_nrows, A_sub_ncols, A_sub_nrows, A_sub_ncols)
        B00 = B.matrix_window(0,           0,           A_sub_ncols, B_sub_ncols)
        B01 = B.matrix_window(0,           B_sub_ncols, A_sub_ncols, B_sub_ncols)
        B10 = B.matrix_window(A_sub_ncols, 0,           A_sub_ncols, B_sub_ncols)
        B11 = B.matrix_window(A_sub_ncols, B_sub_ncols, A_sub_ncols, B_sub_ncols)

        # Allocate temp space.

        # todo: can I cache the bound A.new_empty_window method?

        S0 = A.new_empty_window(A_sub_nrows, A_sub_ncols)
        S1 = A.new_empty_window(A_sub_nrows, A_sub_ncols)
        S2 = A.new_empty_window(A_sub_nrows, A_sub_ncols)
        S3 = A.new_empty_window(A_sub_nrows, A_sub_ncols)

        T0 = A.new_empty_window(A_sub_ncols, B_sub_ncols)
        T1 = A.new_empty_window(A_sub_ncols, B_sub_ncols)
        T2 = A.new_empty_window(A_sub_ncols, B_sub_ncols)
        T3 = A.new_empty_window(A_sub_ncols, B_sub_ncols)

        Q0 = A.new_empty_window(A_sub_nrows, B_sub_ncols)
        Q1 = A.new_empty_window(A_sub_nrows, B_sub_ncols)
        Q2 = A.new_empty_window(A_sub_nrows, B_sub_ncols)


        # Preparatory matrix additions/subtractions.

        # todo: we can probably save some memory in these
        # operations by reusing some of the buffers, if we interleave
        # these additions with the multiplications (below)

        # (I believe we can save on one S buffer and one T buffer)

        # S0 = A10 + A11,  T0 = B01 - B00
        # S1 = S0 - A00,   T1 = B11 - T0
        # S2 = A00 - A10,  T2 = B11 - B01
        # S3 = A01 - S1,   T3 = B10 - T1

        S0.set_to_sum(A10, A11)
        S1.set_to_diff(S0, A00)
        S2.set_to_diff(A00, A10)
        S3.set_to_diff(A01, S1)

        T0.set_to_diff(B01, B00)
        T1.set_to_diff(B11, T0)
        T2.set_to_diff(B11, B01)
        T3.set_to_diff(B10, T1)


        # The relations we need now are:

        # P0 = A00*B00
        # P1 = A01*B10
        # P2 =  S0*T0
        # P3 =  S1*T1
        # P4 =  S2*T2
        # P5 =  S3*B11
        # P6 = A11*T3

        # U0 = P0 + P1
        # U1 = P0 + P3
        # U2 = U1 + P4
        # U3 = U2 + P6
        # U4 = U2 + P2
        # U5 = U1 + P2
        # U6 = U5 + P5

        # We place the final answer into the matrix:
        # U0 U6
        # U3 U4

        U0 = C.matrix_window(0,           0,           A_sub_nrows, B_sub_ncols)
        U6 = C.matrix_window(0,           B_sub_ncols, A_sub_nrows, B_sub_ncols)
        U3 = C.matrix_window(A_sub_nrows, 0,           A_sub_nrows, B_sub_ncols)
        U4 = C.matrix_window(A_sub_nrows, B_sub_ncols, A_sub_nrows, B_sub_ncols)

        if (A_sub_nrows <= cutoff) or (A_sub_ncols <= cutoff) or (B_sub_ncols <= cutoff):
            # This is the base case, so we use MatrixWindow methods directly.

            # This next chunk is arranged so that each output cell gets written
            # to exactly once. This is important because the output blocks might
            # be quite fragmented in memory, whereas our temporary buffers
            # (Q0, Q1, Q2) will be quite localised, so we can afford to do a bit
            # of arithmetic in them.

            Q0.set_to_prod(A00, B00)         # now Q0 holds P0
            Q1.set_to_prod(A01, B10)         # now Q1 holds P1
            U0.set_to_sum(Q0, Q1)            # now U0 is correct
            Q0.add_prod(S1, T1)              # now Q0 holds U1
            Q1.set_to_prod(S2, T2)           # now Q1 holds P4
            Q1.add(Q0)                       # now Q1 holds U2
            Q2.set_to_prod(A11, T3)          # now Q2 holds P6
            U3.set_to_sum(Q1, Q2)            # now U3 is correct
            Q2.set_to_prod(S0, T0)           # now Q2 holds P2
            U4.set_to_sum(Q2, Q1)            # now U4 is correct
            Q0.add(Q2)                       # now Q0 holds U5
            Q2.set_to_prod(S3, B11)          # now Q2 holds P5
            U6.set_to_sum(Q0, Q2)            # now U6 is correct

        else:
            # Recurse into sub-products.

            self._strassen_window_multiply(Q0, A00, B00, cutoff)   # now Q0 holds P0
            self._strassen_window_multiply(Q1, A01, B10, cutoff)   # now Q1 holds P1
            U0.set_to_sum(Q0, Q1)                            # now U0 is correct
            self._strassen_window_multiply(Q1, S1, T1, cutoff)     # now Q1 holds P3
            Q0.add(Q1)                                       # now Q0 holds U1
            self._strassen_window_multiply(Q1, S2, T2, cutoff)     # now Q1 holds P4
            Q1.add(Q0)                                       # now Q1 holds U2
            self._strassen_window_multiply(Q2, A11, T3, cutoff)    # now Q2 holds P6
            U3.set_to_sum(Q1, Q2)                            # now U3 is correct
            self._strassen_window_multiply(Q2, S0, T0, cutoff)     # now Q2 holds P2
            U4.set_to_sum(Q2, Q1)                            # now U4 is correct
            Q0.add(Q2)                                       # now Q0 holds U5
            self._strassen_window_multiply(Q2, S3, B11, cutoff)    # now Q2 holds P5
            U6.set_to_sum(Q0, Q2)                            # now U6 is correct


        # Now deal with the leftover row and/or column (if they exist).

        if B_ncols & 1:
            B_last_col = B.matrix_window(0, B_ncols-1, A_ncols, 1)
            C_last_col = C.matrix_window(0, B_ncols-1, A_nrows, 1)
            C_last_col.set_to_prod(A, B_last_col)

        if A_nrows & 1:
            A_last_row = A.matrix_window(A_nrows-1, 0, 1, A_ncols)
            if B_ncols & 1:
                B_bulk = B.matrix_window(0, 0, A_ncols, B_ncols-1)
                C_last_row = C.matrix_window(A_nrows-1, 0, 1, B_ncols-1)
            else:
                B_bulk = B
                C_last_row = C.matrix_window(A_nrows-1, 0, 1, B_ncols)
            C_last_row.set_to_prod(A_last_row, B_bulk)

        if A_ncols & 1:
            A_last_col = A.matrix_window(0, A_ncols-1, A_sub_nrows << 1, 1)
            B_last_row = B.matrix_window(A_ncols-1, 0, 1, B_sub_ncols << 1)
            C_bulk = C.matrix_window(0, 0, A_sub_nrows << 1, B_sub_ncols << 1)
            C_bulk.add_prod(A_last_col, B_last_row)

    cdef subtract_strassen_product(self, result, A, B, int cutoff):
        """
        EXAMPLES:
            sage: ?
        """
        cutoff = 1000000 # for testing
        if (result.ncols() < cutoff or result.nrows() < cutoff):
            result.subtract_prod(A, B)
        else:
            to_sub = self._strassen_window_multiply(A, B, cutoff)
            result.subtract(to_sub)


    def _strassen_echelon(self, A, cutoff):
        """
        Compute echelon form, in place.
        Internal function, call with M.echelonize(algorithm="strassen")
        Based on work of Robert Bradshaw and David Harvey at MSRI workshop in 2006.

        INPUT:
            A -- matrix window
            cutoff -- size at which algorithm reverts to naive gaussian elemination and multiplication

        OUTPUT:
            The list of pivot columns

        EXAMPLE:
            sage: import sage.matrix.matrix_rational_dense as m
            sage: parent = MatrixSpace(QQ, 5, 30)
            sage: data = parent.random_element(range(18), prob=.2).list() # test lots of non-pivots
            sage: A = m.Matrix_rational_dense(parent, data)
            sage: T = A.echelonize(alg="gauss")
            sage: E = A.echelon_strassen(4)
            sage: A.copy() == T.copy()  # fix when parents are changed
            True

        AUTHORS:
            -- Robert Bradshaw
        """
        # The following notation will be used in the comments below, which should be understood to give
        # the general idea of what's going on, as if there were no inconvenient non-pivot columns.
        # The original matrix is given by [ A B ]
        #                                 [ C D ]
        # For compactness, let A' denote the inverse of A
        # top_left, top_right, bottom_left, and bottom_right loosely correspond to A, B, C, and D respectively,
        # however, the "cut" between the top and bottom rows need not be the same.

        cdef int nrows, ncols
        nrows = A.nrows()
        ncols = A.ncols()

        if (nrows < cutoff or ncols < cutoff):
            return A.echelon_in_place()

        cdef int top_h, bottom_cut, bottom_h, bottom_start, top_cut
        cdef int prev_pivot_count
        cdef int split
        split = nrows / 2

        top = A.matrix_window(0, 0, split, ncols)
        bottom = A.matrix_window(split, 0, nrows-split, ncols)

        top_pivots = self._strassen_echelon(top, cutoff)
        # effectively "multiplied" top row by A^{-1}
        #     [  I  A'B ]
        #     [  C   D  ]

        top_pivot_intervals = int_range(top_pivots)
        top_h = len(top_pivots)

        if top_h == 0:
            #                                 [ 0 0 ]
            # the whole top is a zero matrix, [ C D ]. Run echelon on the bottom
            bottom_pivots = self._strassen_echelon(bottom, cutoff)
            #             [  0   0  ]
            # we now have [  I  C'D ], proceed to sorting

        else:
            bottom_cut = max(top_pivots) + 1
            bottom_left = bottom.matrix_window(0, 0, nrows-split, bottom_cut)

            if top_h == ncols:
                bottom.set_to_zero()
                # [ I ]
                # [ 0 ]
                # proceed to sorting

            else:
                if bottom_cut == top_h:
                    clear = bottom_left
                else:
                    clear = bottom_left.to_matrix().matrix_from_cols(top_pivots).matrix_window() # TODO: read only, can I do this faster? Also below
                # Subtract off C time top from the bottom_right
                if bottom_cut < ncols:
                    bottom_right = bottom.matrix_window(0, bottom_cut, nrows-split, ncols-bottom_cut)
                    self.subtract_strassen_product(bottom_right, clear, top.matrix_window(0, bottom_cut, top_h, ncols-bottom_cut), cutoff);
                # [  I      A'B   ]
                # [  *   D - CA'B ]

                # Now subtract off C times the top from the bottom_left (pivots -> 0)
                if bottom_cut == top_h:
                    bottom_left.set_to_zero()
                    bottom_start = bottom_cut

                else:
                    for cols in top_pivot_intervals:
                        bottom_left.matrix_window(0, cols[0], nrows-split, cols[1]).set_to_zero()
                    non_pivots = int_range(0, bottom_cut) - top_pivot_intervals
                    for cols in non_pivots:
                        if cols[0] == 0: continue
                        prev_pivot_count = len(top_pivot_intervals - int_range(cols[0]+cols[1], bottom_cut - cols[0]+cols[1]))
                        self.subtract_strassen_product(bottom_left.matrix_window(0, cols[0], nrows-split, cols[1]),
                                                         clear.matrix_window(0, 0, nrows-split, prev_pivot_count),
                                                         top.matrix_window(0, cols[0], prev_pivot_count, cols[1]),
                                                         cutoff)
                    bottom_start = non_pivots._intervals[0][0]
                # [  I      A'B   ]
                # [  0   D - CA'B ]

                # Now recursively do echelon form on the bottom
                bottom_pivots_rel = self._strassen_echelon(bottom.matrix_window(0, bottom_start, nrows-split, ncols-bottom_start), cutoff)
                # [  I  A'B ]
                # [  0  I F ]
                bottom_pivots = []
                for pivot in bottom_pivots_rel:
                    bottom_pivots.append(pivot + bottom_start)
                bottom_h = len(bottom_pivots)

                if bottom_h + top_h == ncols:
                    top.matrix_window(0, bottom_cut, split, ncols-bottom_cut).set_to_zero()
                    # [ I 0 ]
                    # [ 0 I ]
                    # proceed to sorting

                elif bottom_h == 0:
                    pass
                    # [  I  A'B ]
                    # [  0   0  ]
                    # proceed to sorting

                else:
                    #     [  I  A'B ]  =  [  I  E  G  ]
                    # let [  0  I F ]  =  [  0  I  F  ]
                    top_cut = max(max(bottom_pivots) + 1, bottom_cut)

                    # Note: left with respect to leftmost non-zero column of bottom
                    top_left = top.matrix_window(0, bottom_start, top_h, top_cut - bottom_start)

                    if top_cut - bottom_start == bottom_h:
                        clear = top_left
                    else:
                        clear = top_left.to_matrix().matrix_from_cols(bottom_pivots_rel).matrix_window()

                    # subtract off E times bottom from top right
                    if top_cut < ncols:
                        top_right = top.matrix_window(0, top_cut, top_h, ncols - top_cut)
                        self.subtract_strassen_product(top_right, clear, bottom.matrix_window(0, top_cut, top_h, ncols - top_cut), cutoff);
                    # [  I  *  G - EF ]
                    # [  0  I     F   ]

                    # Now subtract of E times bottom from top left
                    if top_cut - bottom_start == bottom_h:
                        top_left.set_to_zero()

                    else:
                        bottom_pivot_intervals = int_range(bottom_pivots)
                        for cols in bottom_pivot_intervals:
                            top.matrix_window(0, cols[0], split, cols[1]).set_to_zero()
                        non_pivots = int_range(bottom_start, top_cut - bottom_start) - bottom_pivot_intervals - top_pivot_intervals
                        for cols in non_pivots:
                            if cols[0] == 0: continue
                            prev_pivot_count = len(bottom_pivot_intervals - int_range(cols[0]+cols[1], top_cut - cols[0]+cols[1]))
                            self.subtract_strassen_product(top.matrix_window(0, cols[0], split, cols[1]),
                                                             clear.matrix_window(0, 0, split, prev_pivot_count),
                                                             bottom.matrix_window(0, cols[0], prev_pivot_count, cols[1]),
                                                             cutoff)
                    # [  I  0  G - EF ]
                    # [  0  I     F   ]
                    # proceed to sorting

        # subrows already sorted...maybe I could do this more efficiently in cases with few pivot columns (e.g. merge sort)

        pivots = top_pivots
        pivots.extend(bottom_pivots)
        pivots.sort()

        cdef int i, cur_row
        for cur_row from 0 <= cur_row < len(pivots):
            pivot = pivots[cur_row]
            for i from cur_row <= i < nrows:
                if not A.element_is_zero(i, pivot):
                    break
            if i > cur_row:
                A.swap_rows(i, cur_row)

        return pivots

############################# end matrix class ####################


#########################################################################
# Generic matrix windows, which are used for block echelon and strassen #
# algorithms.                                                           #
#########################################################################
cdef class MatrixWindow:

    def __init__(MatrixWindow self, Matrix matrix, int row, int col, int nrows, int ncols):

        self._matrix = matrix
        self._row = row
        self._col = col
        self._nrows = nrows
        self._ncols = ncols

    def __repr__(self):
        return "Matrix window of size %s x %s at (%s,%s):\n%s"%(
            self._nrows, self._ncols, self._row, self._col, self._matrix)

    def set_unsafe(self, Py_ssize_t i, Py_ssize_t j, x):
        self._matrix.set_unsafe(i + self._row, j + self._col, x)

    def get_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        return self._matrix.get_unsafe(i + self._row, j + self._col)

    def __setitem__(self, ij, x):
        cdef Py_ssize_t i, j
        if PyTuple_Check(ij):
            # ij is a tuple, so we get i and j efficiently, construct corresponding integer entry.
            if PyTuple_Size(ij) != 2:
                raise IndexError, "index must be an integer or pair of integers"
            i = <object> PyTuple_GET_ITEM(ij, 0)
            j = <object> PyTuple_GET_ITEM(ij, 1)
            if i<0 or i >= self._nrows or j<0 or j >= self._ncols:
                raise IndexError, "matrix index out of range"
            self.set_unsafe(i, j, x)
        else:
            # If ij is not a tuple, coerce to an integer and set the row.
            i = ij
            for j from 0 <= j < self._ncols:
                self.set_unsafe(i, j, x)

    def __getitem__(self, ij):
        cdef Py_ssize_t i, j
        cdef object x

        if PyTuple_Check(ij):
            # ij is a tuple, so we get i and j efficiently, construct corresponding integer entry.
            if PyTuple_Size(ij) != 2:
                raise IndexError, "index must be an integer or pair of integers"
            i = <object> PyTuple_GET_ITEM(ij, 0)
            j = <object> PyTuple_GET_ITEM(ij, 1)
            if i<0 or i >= self._nrows or j<0 or j >= self._ncols:
                raise IndexError, "matrix index out of range"
            return self.get_unsafe(i, j)
        else:
            # If ij is not a tuple, coerce to an integer and get the row.
            i = ij
            return self.row(i)

    def matrix(MatrixWindow self):
        """
        Returns the underlying matrix that this window is a view of.
        """
        return self._matrix


    def to_matrix(MatrixWindow self):
        """
        Returns an actual matrix object representing this view.
        """
        a = self._matrix.new_matrix(self._nrows, self._ncols)
        a.matrix_window().set_to(self)
        return a

    def matrix_window(MatrixWindow self, int row=0, int col=0, int n_rows=-1, int n_cols=-1):
        """
        Returns a matrix window relative to this window of the underlying matrix.
        """
        if row == 0 and col == 0 and n_rows == self._nrows and n_cols == self._ncols:
            return self
        return self._matrix.matrix_window(self._row + row, self._col + col, n_rows, n_cols)

    def nrows(MatrixWindow self):
        return self._nrows

    def ncols(MatrixWindow self):
        return self._ncols

    def set_to(MatrixWindow self, MatrixWindow A):
        """
        Change self, making it equal A.
        """
        cdef Py_ssize_t i, j
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                self.set_unsafe(i, j, A.get_unsafe(i, j))

    def set_to_zero(MatrixWindow self):
        cdef Py_ssize_t i, j
        z = self._base_ring(0)
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                self._matrix.matrix.set_unsafe(i, j, z)

    def add(MatrixWindow self, MatrixWindow A):
        cdef Py_ssize_t i, j
        if self._nrows != A._nrows or self._ncols != A._ncols:
            raise ArithmeticError, "incompatible dimensions"
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                self.set_unsafe(i, j, self.get_unsafe(i, j) + A.get_unsafe(i, j))

    def subtract(MatrixWindow self, MatrixWindow A):
        cdef Py_ssize_t i, j
        if self._nrows != A._nrows or self._ncols != A._ncols:
            raise ArithmeticError, "incompatible dimensions"
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                self.set_unsafe(i, j, self.get_unsafe(i, j) - A.get_unsafe(i, j))

    def set_to_sum(MatrixWindow self, MatrixWindow A, MatrixWindow B):
        cdef Py_ssize_t i, j
        if self._nrows != A._nrows or self._ncols != A._ncols:
            raise ArithmeticError, "incompatible dimensions"
        if self._nrows != B._nrows or self._ncols != B._ncols:
            raise ArithmeticError, "incompatible dimensions"
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                self.set_unsafe(i, j, A.get_unsafe(i, j) + B.get_unsafe(i, j))

    def set_to_diff(MatrixWindow self, MatrixWindow A, MatrixWindow B):
        cdef Py_ssize_t i, j
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                self.set_unsafe(i, j, A.get_unsafe(i, j) - B.get_unsafe(i, j))

    def set_to_prod(MatrixWindow self, MatrixWindow A, MatrixWindow B):
        cdef Py_ssize_t i, j, k
        if A._ncols != B._nrows or self._nrows != A._nrows or self._ncols != B._ncols:
            raise ArithmeticError, "incompatible dimensions"
        for i from 0 <= i < A._nrows:
            for j from 0 <= j < B._ncols:
                s = 0
                for k from 0 <= k < A._ncols:
                    s = s + A.get_unsafe(i, k) * B.get_unsafe(k, j)
                self.set_unsafe(i, j, s)

    def add_prod(MatrixWindow self, MatrixWindow A, MatrixWindow B):
        cdef Py_ssize_t i, j, k
        if A._ncols != B._nrows or self._nrows != A._nrows or self._ncols != B._ncols:
            raise ArithmeticError, "incompatible dimensions"
        for i from 0 <= i < A._nrows:
            for j from 0 <= j < B._ncols:
                s = self.get_unsafe(i, j)
                for k from 0 <= k < A._ncols:
                    s = s + A.get_unsafe(i, k) * B.get_unsafe(k, j)
                self.set_unsafe(i, j, s)

    def subtract_prod(MatrixWindow self, MatrixWindow A, MatrixWindow B):
        cdef Py_ssize_t i, j, k
        if A._ncols != B._nrows or self._nrows != A._nrows or self._ncols != B._ncols:
            raise ArithmeticError, "incompatible dimensions"
        for i from 0 <= i < A._nrows:
            for j from 0 <= j < B._ncols:
                s = self.get_unsafe(i, j)
                for k from 0 <= k < A._ncols:
                    s = s - A.get_unsafe(i, k) * B.get_unsafe(k, j)
                self.set_unsafe(i, j, s)

    def swap_rows(MatrixWindow self, int a, int b):
        self._matrix.swap_rows(self._row + a, self._row + b)

    def echelon_in_place(MatrixWindow self):
        """
        Calculate the echelon form of this matrix, returning the list of pivot columns
        """
        echelon = self.to_matrix()
        echelon.echelonize() # TODO: read only, only need to copy pointers
        self.set_to(echelon.matrix_window())
        return echelon.pivots()

    def element_is_zero(MatrixWindow self, int i, int j):
        return self._matrix.matrix[i+self._row][j+self._col] == 0

    def new_empty_window(MatrixWindow self, int nrows, int ncols):
        return self._matrix.new_matrix(nrows, ncols).matrix_window()


################################

# lots of room for optimization....
# eventually, should I just pass these around rather than lists of ints for pivots?
# would need new from_cols
class int_range:
    r"""
    Useful class for dealing with pivots in the strassen echelon, could have much more general application
    AUTHORS:
      -- Robert Bradshaw

    """
    def __init__(self, indices=None, range=None):
        if indices is None:
            self._intervals = []
            return
        elif not range is None:
            self._intervals = [(int(indices), int(range))]
        else:
            self._intervals = []
            if len(indices) == 0:
                return
            indices.sort()
            start = None
            last = None
            for ix in indices:
                if last is None:
                    start = ix
                elif ix-last > 1:
                    self._intervals.append((start, last-start+1))
                    start = ix
                last = ix
            self._intervals.append((start, last-start+1))

    def __repr__(self):
        return str(self._intervals)

    def intervals(self):
        return self._intervals

    def to_list(self):
        all = []
        for iv in self._intervals:
            for i in range(iv[0], iv[0]+iv[1]):
                all.append(i)
        return all

    def __iter__(self):
        return self._intervals.__iter__()

    def __len__(self):
        len = 0
        for iv in self._intervals:
            len = len + iv[1]
        return len

    # Yes, these two could be a lot faster...
    # Basically, this class is for abstracting away what I was trying to do by hand in several places
    def __add__(self, right):
        all = self.to_list()
        for i in right.to_list():
            all.append(i)
        return int_range(all)

    def __sub__(self, right):
        all = self.to_list()
        for i in right.to_list():
            if i in all:
                all.remove(i)
        return int_range(all)

    def __mul__(self, right):
        intersection = []
        all = self.to_list()
        for i in right.to_list():
            if i in all:
                intersection.append(i)
        return int_range(intersection)

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


cdef decomp_seq(v):
    return Sequence(v, universe=tuple, check=False, cr=True)
