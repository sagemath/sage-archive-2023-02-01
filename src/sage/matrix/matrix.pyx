r"""
Matrices over an arbitrary ring

AUTHORS:
    -- William Stein
    -- Martin Albrecht: conversion to Pyrex
    -- Jaap Spies: various functions
    -- Gary Zablackis: fixed a sign bug in generic determinant.

Elements of matrix spaces are of class \code{Matrix} (or a class
derived from Matrix).  They can be either sparse or dense, and can be
defined over any base ring.


EXAMPLES:

We create the $2\times 3$ matrix
$$
   \left(\begin{matrix} 1&2&3\\4&5&6 \end{matrix}\right)
$$
as an element of a matrix space over $\Q$:

    sage: M = MatrixSpace(RationalField(),2,3)
    sage: A = M([1,2,3, 4,5,6])
    sage: A
    [1 2 3]
    [4 5 6]
    sage: A.parent()
    Full MatrixSpace of 2 by 3 dense matrices over Rational Field

We next change the top-right entry of $A$.  Note that matrix indexing
is $0$-based in SAGE, so the top right entry is $(0,2)$, which should
be thought of as ``row number $0$, column number 2''.

    sage: A[0,2] = 389
    sage: A
    [  1   2 389]
    [  4   5   6]

Also notice how matrices print.  All columns have the same width and
entries in a given column are right justified.   Next we compute the
reduced row echelon form of $A$.

    sage: A.echelon_form()
    [      1       0 -1933/3]
    [      0       1  1550/3]

We save and load a matrix:
    sage: A = Matrix(Integers(8),3,range(9))
    sage: loads(dumps(A)) == A
    True


MUTABILITY: Matrices are either immutable or not.  When initially
created, matrices are typically mutable, so one can change their
entries.  Once a matrix $A$ is made immutable using
\code{A.set_immutable()} the entries of $A$ cannot be changed, and $A$
can never be made mutable again.  However, properies of $A$ such as
its rank, characteristic polynomial, etc., are all cached so
computations involving $A$ may be more efficient.  Once $A$ is made
immutable it cannot be changed back.  However, one can obtain a
mutable copy of $A$ using \code{A.copy()}.

EXAMPLES:
    sage: A = Matrix(RR,2,[1,10,3.5,2])
    sage: A.set_immutable()
    sage: A.copy() is A
    False

The echelon form method always returns immutable matrices with known
rank.

EXAMPLES:
    sage: A = Matrix(Integers(8),3,range(9))
    sage: A.determinant()
    0
    sage: A[0,0] = 5
    sage: A.determinant()
    1
    sage: A.set_immutable()
    sage: A[0,0] = 5
    Traceback (most recent call last):
    ...
    ValueError: object is immutable; please change a copy instead.


\subsection{Implementation Discussion}
Class Diagram:
\begin{verbatim}
Matrix (*) -- abstract base
    Matrix_generic_dense
    Matrix_generic_sparse
    Matrix_integer_dense
    Matrix_integer_sparse
    Matrix_rational_dense
    Matrix_rational_sparse
    Matrix_cyclo_dense
    Matrix_cyclo_sparse
    Matrix_modn_dense
    Matrix_modn_sparse
    Matrix_RR_dense
    Matrix_RR_sparse
    Matrix_CC_dense
    Matrix_CC_sparse
    Matrix_RDF_dense
    Matrix_RDF_sparse
    Matrix_CDF_dense
    Matrix_CDF_sparse
\end{verbatim}

The corresponding files in the sage/matrix library code
directory are named

          [matrix] [base ring] [dense or sparse].

See the files \code{matrix_template.pxd} and \code{matrix_template.pyx}.

LEVEL 1. For each base field it is necessary to completely implement the following
   functionality for that base ring:
   * __new__       -- should use sage_malloc from ext/stdsage.pxi  (only needed if allocate memory)
   * __init__      -- this signature: 'def __init__(self, parent, entries, copy, coerce)'
   * __dealloc__   -- use sage_free (only needed if allocate memory)
   * set_unsafe(self, size_t i, size_t j, x) -- a cdef'd method that doesn't do bounds or any other checks; can core dump if given bad i,j
   * get_unsafe(self, size_t i, size_t j) -- a cdef'd method that doesn't do checks

LEVEL 2. After getting the special class with the above code to work, implement any subset
         of the following purely for speed reasons:
   * cdef pickle(self):
          return data, version
   * cdef unpickle(self, data, int version):
          reconstruct matrix from given data and version; may assume _parent, _nrows, and _ncols are set.
          Use version numbers so if you change the pickle strategy then
          old objects still unpickle.
   * _add_sibling_cdef
   * _sub_sibling_cdef
   * _mul_cousin_cdef
   * _cmp_sibling_cdef
   * __neg__
   * __invert__
   * __copy__
   * __deepcopy__
   * multiply_classical
   * list -- copy of the list of underlying elements
   * dict -- copy of the sparse dictionary of underlying elements

LEVEL 3
   * Matrix windows -- only if you need strassen for that base

   * Other functions, e.g., transpose, for which knowing the
     specific representation can be helpful.

NOTES:
   * For cacheing, use self.fetch and self.cache.
   * Any method that can change the matrix should call check_mutability() first.

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

import sage.modules.free_module
import sage.misc.latex
import sage.structure.coerce
from   sage.structure.sequence import _combinations
import sage.rings.integer

cimport sage.structure.element

from sage.structure.mutability cimport Mutability

def is_Matrix(x):
    """
    EXAMPLES:
        sage: ???
    """
    return isinstance(x, Matrix)

cdef class MatrixWindow  # forward declare

cdef class Matrix(sage.structure.element.ModuleElement):
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

    Note that matrices in SAGE are (currently) only implemented over commutative rings:
        sage: Q = QuaternionAlgebra(QQ, -1,-1)
        sage: matrix(Q,2,1,[1,2])
        Traceback (most recent call last):
        ...
        TypeError: base_ring must be a commutative ring

    (This is just because this assumption is implicitly made in some functions.  It
    could be fixed if there were interest.)
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
        self._nrows = parent.nrows()
        self._ncols = parent.ncols()
        self._mutability = Mutability(False)
        self._cache = {}

    def __hash__(self):
        """
        EXAMPLES:
            sage: A = Matrix(ZZ[['t']], 2, 2, range(4))
            sage: B = A.copy()
            sage: A.__hash__() == B.__hash__()
            True
            sage: A[0,0] = -1
            sage: A.__hash__() == B.__hash__()
            False
        """
        h = self.fetch('hash')
        if not h is None: return h

        if not self._mutability._is_immutable:
            raise TypeError, "mutable matrices are unhashable"

        h = hash(str(self))
        self.cache('hash', h)
        return h

    def __copy__(self):
        """
        Return a copy of this matrix.  Changing the entries of the
        copy will not change the entries of this matrix.

        EXAMPLES:
            sage: ???
        """
        return self.new_matrix(entries=self.list(), coerce=False, copy=False)

    def copy(self):
        """
        Make a copy of self.

        WARNING: The individual elements aren't themselves copied
        (though the list is copied).    This shouldn't matter, since ring
        elements are (almost!) always immutable in SAGE.

        EXAMPLES:
            sage: R.<x> = QQ['x']
            sage: a = matrix(R,2,[x+1,2/3,  x^2/2, 1+x^3]); a
            [  x + 1     2/3]
            [1/2*x^2 x^3 + 1]

            sage: b = copy(a)

            sage: b[0,0] = 5

            sage: b
            [      5     2/3]
            [1/2*x^2 x^3 + 1]

            sage: a
            [  x + 1     2/3]
            [1/2*x^2 x^3 + 1]

            sage: b = copy(a)

            sage: f = b[0,0]; f[0] = 10

            sage: b
            [ x + 10     2/3]
            [1/2*x^2 x^3 + 1]

            sage: a
            [ x + 10     2/3]
            [1/2*x^2 x^3 + 1]
        """
        return self.__copy__()

    def list(self):
        """
        EXAMPLES:
            sage: ???
        """
        cdef size_t i, j

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
        EXAMPLES:
            sage: ???
        """
        cdef size_t i, j

        x = self.fetch('dict')
        if not x is None:
            return x
        x = {}
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < nself._ncols:
                x = self.get_unsafe(i, j)
                if x != 0:
                    x[(i,j)] = x
        self.cache('dict', x)
        return x

    ###########################################################
    # Cache
    ###########################################################
    cdef clear_cache(self):
        self._cache = {}

    cdef fetch(self, key):
        try:
            return self._cache[key]
        except KeyError:
            return None

    cdef cache(self, key, x):
        self._cache[key] = x

    ###########################################################
    # Mutability and bounds checking
    ###########################################################

    cdef check_bounds(self, size_t i, size_t j):
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

    cdef check_bounds_and_mutability(self, size_t i, size_t j):
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
        """
        EXAMPLES:
            sage: A = Matrix(QQ['x','y'], 2, 2, range(4))
            sage: A.is_mutable()
            True
            sage: A[0,0] = 10
            sage: A
            [10   1]
            [ 2   3]
            sage: A.set_immutable()
            sage: A.is_mutable()
            False
            sage: A[0,0] = 10
            Traceback (most recent call last):
            ...
            <type 'exceptions.ValueError'>: matrix is immutable; please change a copy instead (use self.copy())
        """
        self._mutability.set_immutable()

    def is_immutable(self):
        """
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
    cdef set_unsafe(self, size_t i, size_t j, object x):
        """
        Set entry, but potentially faster because it might be without
        bounds checking.
        """
        raise NotImplementedError, "this must be defined in the derived type."

    cdef get_unsafe(self, size_t i, size_t j):
        """
        Entry access, but potentially faster since it might be without
        bounds checking.
        """
        raise NotImplementedError, "this must be defined in the derived type."

    def __getitem__(self, ij):
        """
        Return element or row of self.

        INPUT:
            ij -- tuple (i,j) with i, j ints
        or
            ij -- int

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
        cdef size_t i, j
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
        cdef size_t i, j

        if PyTuple_Check(ij):
            # ij is a tuple, so we get i and j efficiently, construct corresponding integer entry.
            if PyTuple_Size(ij) != 2:
                raise IndexError, "index must be an integer or pair of integers"
            i = <object> PyTuple_GET_ITEM(ij, 0)
            j = <object> PyTuple_GET_ITEM(ij, 1)
            self.check_bounds_and_mutability(i, j)
            self.set_unsafe(i, j, x)
        else:
            # If ij is not a tuple, coerce to an integer and set the row.
            i = ij
            for j from 0 <= j < self._ncols:
                self.set_unsafe(i, j, x)

    def column(self, size_t i):
        """
        Return the i-th column of this matrix.

        EXAMPLES:
            sage: ???
        """
        cdef size_t j
        V = sage.modules.free_module.FreeModule(self.base_ring(), self._ncols)
        tmp = []
        for j from 0 <= j < self._nrows:
            tmp.append(self.get_unsafe(j,i))
        return V(tmp)

    def columns(self):
        """
        Returns the list of columns of self, as vectors.

        EXAMPLES:
            sage: A = MatrixSpace(RationalField(), 2,3) (xrange(6))
            sage: A
            [0 1 2]
            [3 4 5]
            sage: A.columns()
            [(0, 3), (1, 4), (2, 5)]
        """
        cdef size_t j
        x = self.fetch('columns')
        if not x is None: return x
        tmp = []
        for j from 0 <= j < self._nrows:
            tmp.append(self.column(j))
        self.cache('columns', tmp)
        return tmp

    def row(self, size_t i):
        """
        Return the i-th row of this matrix.

        EXAMPLES:
            sage: ???
        """
        cdef size_t j
        V = sage.modules.free_module.FreeModule(self.base_ring(), self.ncols())
        tmp = []
        for j from 0 <= j < self._ncols:
            tmp.append(self.get_unsafe(i,j))
        return V(tmp)

    def rows(self):
        """
        Return a list of the rows of this matrix.

        EXAMPLES:
            sage: ???
        """
        cdef size_t i
        x = self.fetch('rows')
        if not x is None: return x
        tmp = []
        for i from 0 <= i < self._nrows:
            tmp.append(self.row(i))
        self.cache('rows', tmp)
        return tmp

    ###########################################################
    # Pickling
    ###########################################################
    def __reduce__(self):
        data, version = self.pickle()
        return unpickle, (self.__class__, self._parent, data, version)

    cdef pickle(self):
        version = 0
        data = self.list()  # linear list of all elements

    cdef unpickle(self, data, int version):
        cdef size_t i, j, k
        if version == 0:
            # data is a *list* of the entries of the matrix.
            # TODO: Change the data[k] below to use the fast list access macros from the Python/C API
            k = 0
            for i from 0 <= i < self._nrows:
                for j from 0 <= j < self._ncols:
                    self.set_unsafe(i, j, data[k])
                    k = k + 1
        else:
            raise RuntimeError, "unknown version"


    ###########################################################
    # Base Change
    ###########################################################

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
            sage: A.change_ring(GF(25))
            [3 2]
            [2 4]
            sage: A.change_ring(GF(25)).parent()
             Full MatrixSpace of 2 by 2 dense matrices over Finite Field in a of size 5^2
            sage: A.change_ring(ZZ)
            Traceback (most recent call last):
            ...
            <type 'exceptions.TypeError'>: Unable to coerce rational (=1/2) to an Integer.
        """
        if ring is self.base_ring():
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
        """
        EXAMPLES:
            sage: ???
        """
        cdef size_t nr, nc, r, c
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

        return s

    def _latex_sparse(self, variable="x"):
        r"""
        Return a latex string that represents this matrix as a sparse
        matrix.  The rows are printed as sums $\sum a_i x_i$, where
        $x$ is the variable.

        EXAMPLES:
            sage: ???
        """
        cdef size_t nr, nc, i, j
        nr = self._nrows
        nc = self._ncols
        s = "\\left(\\begin{align*}\n"
        for i from 0 <= i < nr:
            v = []
            for j in 0 <= j < nc:
                x = self.get_unsafe(i, j)
                if x != 0:
                    v.append((j, x))
            for j in xrange(len(v)):
                s  = s + "%s*%s_{%s}"%(v[j][1], variable, v[j][0])
                if j == 0:
                    s = s + "& + "
                elif j < len(v) - 1:
                    s =  s + " + "
                else:
                    s =  s + "\\\\\n"
        s = s + "\n\\end{align*}"
        return s

    def _latex_(self):
        """
        EXAMPLES:
            sage: ???
        """
        cdef size_t nr, nc, r, c
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
        cdef size_t nr, nc, i, j
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
        cdef size_t i, j
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
       cdef size_t i, j
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
        K = self.base_ring()._magma_init_()
        if self._nrows != self._ncols:
            s = 'MatrixAlgebra(%s, %s)'%(K, self.nrows())
        else:
            s = 'RMatrixSpace(%s, %s, %s)'%(K, self.nrows(), self.ncols())
        v = []
        for x in self.list():
            v.append(x._magma_init_())
        return s + '![%s]'%(','.join(v))


    ###################################################
    # Construction functions
    ###################################################
    def transpose(self):
        """
        Returns the transpose of self, without changing self.

        EXAMPLES:
        We create a matrix, compute its transpose, and note that the
        original matrix is not changed.
            sage: M = MatrixSpace(QQ,  2)
            sage: A = M([1,2,3,4])
            sage: B = A.transpose()
            sage: print B
            [1 3]
            [2 4]
            sage: print A
            [1 2]
            [3 4]
        """
        f = []
        e = self.list()
        (nc, nr) = (self.ncols(), self.nrows())
        for j in xrange(nc):
            for i in xrange(nr):
                f.append(e[i*nc + j])
        return self.new_matrix(nrows = nc, ncols = nr,
                               entries = f, copy=False,
                               coerce=False)

    def antitranspose(self):
        f = []
        e = self.list()
        (nc, nr) = (self.ncols(), self.nrows())
        for j in reversed(xrange(nc)):
            for i in reversed(xrange(nr)):
                f.append(e[i*nc + j])
        return self.new_matrix(nrows = nc, ncols = nr,
                               entries = f, copy=False, coerce=False)

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
        return self.change_ring(self.base_ring().cover_ring())

    def sparse_columns(self):
        """
        EXAMPLES:
            sage: ???
        """
        if self.__sparse_columns!=None:
            return self.__sparse_columns
        else:
            C =[]
            for _ in xrange(self.ncols()):
                C.append({})
            for i, j in self.nonzero_positions():
                C[i][j] = self.get(i,j)
            if self.is_immutable():
                self.__sparse_columns = C
            return C

    def sparse_rows(self):
        """
        EXAMPLES:
            sage: ???
        """
        if self.__sparse_rows != None:
            return self.__sparse_rows
        else:
            R = []
            for _ in xrange(self.ncols()):
                R.append({})
            for i, j in self.nonzero_positions():
                R[i][j] = self.get((i,j))
            if self.is_immutable():
                self.__sparse_rows = R
            return R

    def stack(self, Matrix other):
        """
        Return the augmented matrix self on top of other:
           [ self  ]
           [ other ]

        EXAMPLES:
            sage: M = Matrix(RationalField(), 2, 3, range(6))
            sage: N = Matrix(RationalField(), 1, 3, [10,11,12])
            sage: M.stack(N)
            [ 0  1  2]
            [ 3  4  5]
            [10 11 12]
        """
        cdef size_t r, c

        if not (self.base_ring() is other.base_ring()):
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
            sage: A = MatrixSpace(RationalField(),2, sparse=True)([1,2,0,1])
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
            sage: A = MatrixSpace(RationalField(),2, sparse=False)([1,2,0,1])
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

    def dense_row(self, size_t n):
        """
        Return the n-th row of self as a dense vector.

        EXAMPLES:
            sage: ???
        """
        cdef size_t i
        if n < 0 or n >= self._nrows:
            raise IndexError, "n must be between 0 and the number of rows"
        V = sage.modules.free_module.VectorSpace(self.base_ring(), self.ncols(), sparse=False)
        # TODO -- better way to build a list using Python API??
        tmp = []
        for i from 0 <= i < self._ncols:
            tmp.append(self.get_unsafe(n, i))

        return V(tmp)

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
        return self.matrix_space(nrows, ncols, sparse=sparse).matrix(entries, coerce, copy)

    def augment(self, Matrix other):
        """
        Return the augmented matrix of the form [self | other].

        EXAMPLES:
            sage: M = MatrixSpace(RationalField(),2,2)
            sage: A = M([1,2, 3,4])
            sage: A
            [1 2]
            [3 4]
            sage: N = MatrixSpace(RationalField(),2,1)
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
            sage: M = MatrixSpace(RationalField(),3,4)
            sage: A = M([1,2,3,4, 0,9,8,7, 2/3,3/4,4/5,9/8])
            sage: A
            [  1   2   3   4]
            [  0   9   8   7]
            [2/3 3/4 4/5 9/8]
            sage: N = MatrixSpace(RationalField(),3,2)
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
        if not (self.base_ring() is other.base_ring()):
            raise TypeError, "base rings must be the same"
        if self._nrows != other._nrows:
            raise TypeError, "number of rows must be the same"

        cdef Matrix Z
        Z = self.new_matrix(ncols = self._ncols + other._ncols)

        cdef size_t r, c
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
    ## Properties
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
            sage: M = MatrixSpace(RationalField(),6,7)
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
            sage: R = MPolynomialRing(RationalField(), 2, ['x','y'])
            sage: x, y = R.gens()
            sage: f = x**2 - y**2
            sage: M = MatrixSpace(RationalField(),2)
            sage: A = M([1,2,3,4])
            sage: A.act_on_polynomial(f)
            -12*y^2 - 20*x*y - 8*x^2
        """
        cdef size_t i, j, n

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
    ###################################################

    def swap_columns(self, size_t c1, size_t c2):
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
        cdef size_t nc, r

        nc = self.ncols()
        if c1 < 0 or c1 >= nc:
            raise IndexError, "c1 invalid column"
        if c2 < 0 or c2 >= nc:
            raise IndexError, "c2 invalid column"
        if c1 == c2:
            return

        self.check_mutability()

        for r from 0 <= r < self._nrows:
            a = self.get_unsafe(r, c2)
            b = self.get_unsafe(r, c1)
            self.set_unsafe(r, c1, a)
            self.set_unsafe(r, c2, b)


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
        cdef size_t nc, r

        nr = self.nrows()
        if r1 < 0 or r1 >= nc:
            raise IndexError, "invalid row"
        if r2 < 0 or r2 >= nc:
            raise IndexError, "invalid row"
        if r1 == c2:
            return

        self.check_mutability()

        for c from 0 <= c < self._ncols:
            a = self.get_unsafe(r2, c)
            b = self.get_unsafe(r1, c)
            self.set_unsafe(r1, c, a)
            self.set_unsafe(r2, c, b)

    def linear_combination_of_rows(self, v):
        """
        EXAMPLES:
            sage: ???
        """
        cdef size_t i
        R = self.rows()
        s = 0
        for i from 0 <= i < len(v):
            s = s + v[i] * R[i]
        return s

    def linear_combination_of_columns(self, v):
        """
        EXAMPLES:
            sage: ???
        """
        cdef size_t i

        C = self.columns()
        s = 0
        for i from 0 <= i < len(v):
            s = s + v[i]*C[i]
        return s

    def set_row_to_multiple_of_row(self, i, j, s):
        """
        Set row i equal to s times row j.

        EXAMPLES:
            sage: ???
        """
        self[i] = s*self[j]

    def _set_row_to_negative_of_row_of_A_using_subset_of_columns(self, i, A, size_t r, cols):
        """
        EXAMPLES:
            sage: ???
        """
        # this function exists just because it is useful for modular symbols presentations.
        cdef size_t l
        l = 0
        rows = A.sparse_rows()
        v = rows[r]
        for k in cols:
            if k in v.keys():
                self.set_unsafe(i, l, -v[k])
            l = l + 1


    def add_multiple_of_row(self, size_t i, size_t j, s):
        """
        Add s times row j to row i.

        EXAMPLES:
            sage: ???
        """
        cdef size_t c
        self.check_mutability()
        if i<0 or i >= self._nrows or j<0 or j >= self._nrows:
            raise IndexError, "matrix row index out of range"

        s = self.base_ring()(s)
        for c from 0 <= c < self._ncols:
            self.set_unsafe(i, c, self.get_unsafe(i, c) + s*self.get_unsafe(j, c))

    def add_multiple_of_column(self, size_t i, size_t j, s):
        """
        Add s times column j to column i.

        EXAMPLES:
            sage: ???
        """
        cdef size_t r
        self.check_mutability()
        if i<0 or i >= self._ncols or j<0 or j >= self._ncols:
            raise IndexError, "matrix column index out of range"
        s = self.base_ring()(s)
        for r from 0 <= r < self._nrows:
            self.set_unsafe(r, i, self.get_unsafe(r, i) + s*self.get_unsafe(r, j))

    def rescale_row(self, size_t i, s):
        """
        Replace i-th row of self by s times i-th row of self.

        EXAMPLES:
            sage: ???
        """
        if s == 1:
            return
        self.check_mutability()

        if i<0 or i >= self._nrows:
            raise IndexError, "matrix row index out of range"

        for j from 0 <= j < self._ncols:
            self.set_unsafe(i, j, self.get_unsafe(i, j)*s)


    ###################################################
    # Predicates
    ###################################################

    def is_dense(self):
        """
        EXAMPLES:
            sage: ???
        """
        return self._parent.is_dense()

    def is_sparse(self):
        """
        EXAMPLES:
            sage: ???
        """
        return self._parent.is_sparse()

    def is_square(self):
        """
        EXAMPLES:
            sage: ???
        """
        return self._nrows == self._ncols

    def nonpivots(self):
        """
        Return the list of i such that the i-th column of self
        is NOT a pivot column of the reduced row echelon form
        of self.

        OUTPUT:
            list -- sorted list of integers
        EXAMPLES:
            sage: ???
        """
        X = set(self.pivots())
        tmp = []
        for j in xrange(self.ncols()):
            if not (j in X):
                tmp.append(j)
        return tmp

    def nonzero_positions(self):
        """
        Returns the set of pairs (i,j) such that self[i,j] != 0.

        EXAMPLES:
            sage: ???
        """
        cdef size_t i, j
        z = self.base_ring()(0)
        tmp = set()
        for i from 0 <= i < self._nrows:
           for j from 0 <= j < self._ncols:
                if self.get_unsafe(i,j) != z:
                    tmp.add((i,j))
        return tmp

    def nonzero_positions_in_column(self, size_t i):
        """
        Return the integers j such that self[j,i] is nonzero, i.e.,
        such that the j-th position of the i-th column is nonzero.

        EXAMPLES:
            sage: ???
        """
        cdef size_t j
        z = self.base_ring()(0)
        tmp = set()

        if i<0 or i >= self._ncols:
            raise IndexError, "matrix column index out of range"
        for j from 0 <= j < self._nrows:
            if self.get_unsafe(j,i) != z:
                tmp.add(j)
        return tmp

    def nonzero_positions_in_row(self, size_t i):
        """
        Return the integers j such that self[i,j] is nonzero, i.e.,
        such that the j-th position of the i-th row is nonzero.

        EXAMPLES:
            sage: ???
        """
        cdef size_t j

        if i<0 or i >= self._nrows:
            raise IndexError, "matrix row index out of range"

        z = self.base_ring()(0)
        tmp = set()

        for j from 0 <= j < self._nrows:
            if self.get_unsafe(i,j) != z:
                tmp.add(j)
        return tmp

    def prod_of_row_sums(self, cols):
        r"""
        Calculate the product of all row sums of a submatrix of $A$ for a
        list of selected columns \code{cols}.

        EXAMPLES:
            sage: ???

        AUTHOR:
            -- Jaap Spies (2006-02-18)
        """
        cdef size_t c, row
        pr = 1
        for row from 0 <= row < self._nrows:
            tmp = []
            for c in cols:
                if c<0 or c >= self._ncols:
                    raise IndexError, "matrix column index out of range"
                tmp.append(self.get_unsafe(row, c))
            pr = pr * sum(tmp)
        return pr

    ###################################################
    # Invariants of a matrix
    ###################################################
    def pivots(self):
        try:
            return self._cache['pivots']
        except KeyError:
            self.echelon_form()
        return self._cache['pivots']

    cdef _set_pivots(self, X):
        self.cache('pivots', X)

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

            sage: R = PolynomialRing(IntegerRing(),'a'); a = R.gen()
            sage: A = MatrixSpace(R,2)([[a,1], [a,a+1]])
            sage: A.permanent()
            a^2 + 2*a

            sage: R = MPolynomialRing(IntegerRing(),2); x,y = R.gens()
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
        cdef size_t m, n, r

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

    def permanental_minor(self, size_t k):
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

        R = self.base_ring()
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

            sage: x = PolynomialRing(IntegerRing(),'x').gen()
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
        """
        cdef size_t i, n

        try:
            return self._cache['determinant']
        except KeyError:
            pass

        if self._nrows != self._ncols:
            raise ValueError, "self must be square"

        n = self._nrows
        R = self.parent().base_ring()
        if n == 0: return R(1)

        d = R(0)
        s = R(1)
        A = self.matrix_from_rows(range(1, n))
        sgn = R(-1)
        for i from 0 <= i < n:
            v = range(n)
            del v[i]
            B = A.matrix_from_columns(v)
            d = d + s*self.get_unsafe(0,i) * B.determinant()
            s = s*sgn

        self._cache['determinant'] = d

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
        EXAMPLES:
            sage: ???
        """
        return self.determinant()

    def characteristic_polynomial(self, *args, **kwds):
        """
        Synonym for self.charpoly(...).

        EXAMPLES:
            sage: ???
        """
        return self.charpoly(*args, **kwds)

    def charpoly(self, *args, **kwds):
        raise NotImplementedError

    def minimal_polynomial(self, *args, **kwds):
        """
        Synonym for self.charpoly(...).

        EXAMPLES:
            sage: ???
        """
        return self.minpoly(*args, **kwds)

    def minpoly(self, *args, **kwds):
        """
        EXAMPLES:
            sage: ???
        """
        raise NotImplementedError

    def trace(self):
        """
        Return the trace of self, which is the sum of the
        diagonal entries of self.
        INPUT:
            self -- a square matrix
        OUTPUT:
            element of the base ring of self
        """
        if self.nrows() != self.ncols():
            raise ArithmeticError, "matrix must be square"
        R = self.base_ring()
        tmp =  []
        for i in xrange(self.nrows()):
            tmp.append(self[i,i])

        return sum(tmp,R(0))

    ###################################################
    # Arithmetic
    ###################################################
    def vector_matrix_multiply(self, v):
        """
        Returns the vector times matrix product.

        INPUT:
             v -- a free module element.

        OUTPUT:
            The the vector times matrix product v*A.

        EXAMPLES:
            sage: MS = MatrixSpace(RationalField(), 2,2)
            sage: B = MS.matrix([1,2,1,2])
            sage: V = VectorSpace(RationalField(), 2)
            sage: v = V([1,2])
            sage: B.vector_matrix_multiply(v)     # computes v*B
            (3, 6)
            sage: Bt = B.transpose()
            sage: Bt.vector_matrix_multiply(v)    # computes B*v
            (5, 5)
        """
        M = sage.modules.free_module.FreeModule(self.base_ring(), self.ncols())
        if not PY_TYPE_CHECK(v, sage.modules.free_module_element.FreeModuleElement):
            v = M(v)
        if self.nrows() != v.degree():
            raise ArithmeticError, "number of rows of matrix must equal degree of vector"
        s = M(0)
        for i in xrange(self.nrows()):
            if v[i] != 0:
                s = s + v[i]*self.row(i)
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
        M = sage.modules.free_module.FreeModule(self.base_ring(), self.ncols())
        if not PY_TYPE_CHECK(v, sage.modules.free_module_element.FreeModuleElement):
            v = M(v)
        X = [v]
        for _ in range(n-1):
            X.append(X[len(X)-1]*self)
        MS = self.matrix_space(n, self.ncols())
        return MS(X)

    cdef sage.structure.element.ModuleElement _add_sibling_cdef(self, sage.structure.element.ModuleElement right):
        """
        Add two matrices with the same parent.

        EXAMPLES:
            sage:
        """
        cdef size_t i, j
        cdef Matrix A
        A = self.new_matrix()
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                A.set_unsafe(i,j, self.get_unsafe(i,j) + (<Matrix>right).get_unsafe(i,j))
        return A

##     cdef sage.structure.element.ModuleElement _sub_sibling_cdef(self,
##                                                 sage.structure.element.ModuleElement right):
##         cdef size_t i, j
##         A = self.new_matrix()
##         for i from 0 <= i < self._nrows:
##             for j from 0 <= j < self._ncols:
##                 A.set_unsafe(i,j, self.get_unsafe(i,j) - (<Matrix>right).get_unsafe(i,j))
##         return A

    def _div_(self, right):
        """
        EXAMPLES:
            sage: ???
        """
        if not isinstance(self,Matrix):
            self,right = right,self

        return self*(~right)

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
        cdef size_t i
        v = self.list()
        for i from 0 <= i < len(v):
            v[i] = v[i] % p
        return self.new_matrix(entries = v, copy=False, coerce=True)

    def Mod(self, p):
        """
        Return matrix mod $p$, over the reduced ring.

        EXAMPLES:
            sage: M = Matrix(ZZ, 2, 2, [5, 9, 13, 15])
            sage: M.Mod(7)
            [5 2]
            [6 1]
            sage: parent(M.Mod(7))
            Full MatrixSpace of 2 by 2 dense matrices over Finite Field of size 7
        """
        return self.change_ring(self.base_ring().quotient_ring(p))


    def _scalar_multiply(self, x):
        """
        EXAMPLES:
            sage: ???
        """
        cdef size_t i

        if not x in self.base_ring():
            x = self.base_ring()(x)

        v = self.list()
        for i from 0 <= i < len(v):
            v[i] = x * v[i]
        return self.new_matrix(entries = v, copy=False, coerce=False)

    def _right_scalar_multiply(self, x):
        """
        EXAMPLES:
            sage: ?
        """
        cdef size_t i

        if not x in self.base_ring():
            x = self.base_ring()(x)

        v = self.list()
        for i from 0 <= i < len(v):
            v[i] = v[i]*x
        return self.new_matrix(entries = v, copy=False, coerce=False)

    def _left_scalar_multiply(self, x):
        """
        EXAMPLES:
            sage: ?
        """
        # same as _right, because base ring is always commutative.
        return self._scalar_multiply(x)

    def __mul__(self, right):
        """
        EXAMPLES:
            sage: ???
        """
        if not PY_TYPE_CHECK(self, Matrix):
            # it is not a vector, since if it were the vector __mul__ would be called.
            return right._left_scalar_multiply(self)

        if PY_TYPE_CHECK(right, sage.modules.free_module_element.FreeModuleElement):
            raise TypeError, "cannot multiply matrix times row vector -- instead compute row vector times matrix"

        if not PY_TYPE_CHECK(right, Matrix):
            # the only possibility is to coerce to a scalar.
            # TODO: what if right is not in the
            return self._right_scalar_multiply(right)

        # Now both self and right are matrices.

        # Now we have to use coercion rules.
        # If parents are compatible, multiply. Compatible means that
        # they have the same base ring and are both either dense or sparse.
        #
        #  If not, use the following rules:
        #
        #   * error if the matrix dimensions don't match up,
        #
        #   * check for canonical coercion of the base rings in
        #      exactly one direction -- if not error
        #
        #   * if one matrix is sparse and the other is dense, the
        #     product is dense.


        # First we check that matrix multiplication is defined.
        if (<Matrix> self)._ncols != (<Matrix> right)._nrows:
            raise ArithmeticError, "number of columns of self must equal number of rows of right."

        # check parents directly for speed on square matrices (e.g. 2x2 matrices over ZZ)
        if not ( (<Matrix> self)._parent is (<Matrix> right)._parent
                or (<Matrix> self)._parent.base_ring() is (<Matrix> right)._parent.base_ring() ):
            # The base rings are not the same
            try:
                self, right = sage.structure.coerce.canonical_base_coercion(self, right)
            except TypeError:
                raise TypeError, "base rings must be compatible"
            # Either an error was just raised, or self and right now have the same base ring.

        # Next we deal with the possiblity that one could be sparse and the other dense.
        if self.is_sparse() and not right.is_sparse():
            self = self.dense_matrix()
        elif right.is_sparse() and not self.is_dense():
            right = right.dense_matrix()

        return (<Matrix> self)._mul_cousin_cdef(right)

    cdef _mul_cousin_cdef(self, Matrix right):
        """
        Multiply two matrices that are assumed to be compatable and
        defined over the same base ring.

        EXAMPLES:
            sage: ???
        """
        if self._will_use_strassen(right):
            return self.multiply_strassen(right)
        else:
            return self.multiply_classical(right)

    def multiply_classical(self, Matrix right):
        """
        Multiply the matrices self and right using the classical $O(n^3)$
        algorithm.

        This method assumes that self and right have the same parent and
        compatable dimensions.

        EXAMPLES
            sage: include the 0 rows and 0 columns cases
        """
        cdef size_t i, j, k
        if self._ncols != right._nrows:
            raise IndexError, "Number of columns of self must equal number of rows of other."
        v = []
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < right._ncols:
                s = 0
                for k from 0 <= k < self._ncols:
                    s = s + self.get_unsafe(i,k) * right.get_unsafe(k,j)
                v.append(s)
        return self.new_matrix(self._nrows, right._ncols, entries = v, coerce=False, copy=False)

    cdef int _will_use_strassen(self, Matrix right) except -1:
        """
        Whether or not matrix multiplication of self by right should
        be done using Strassen.

        Overload this in derived classes to not use Strassen by default.

        EXAMPLES:
            sage: ???
        """
        cdef int n
        n = self._strassen_default_cutoff(right)
        if self._nrows > n and self._ncols > n and right._nrows > n and right._ncols > n:
            return 1
        return 0

    def __neg__(self):
        """
        EXAMPLES:
            sage: ?
        """
        return self._left_scalar_multiply(self.base_ring()(-1))

    def __pos__(self):
        """
        EXAMPLES:
            sage: ?
        """
        return self

    def __pow__(self, n, ignored):
        """
        EXAMPLES:
            sage: MS = MatrixSpace(RationalField(), 3, 3)
            sage: A = MS([0, 0, 1, 1, 0, '-2/11', 0, 1, '-3/11'])
            sage: A * A**(-1) == 1
            True
            sage: A**4
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

    def __rmul__(self, left):
        """
        EXAMPLES:
            sage: ?
        """
        return self._left_scalar_multiply(left)

    def _sub_(self, right):
        """
        EXAMPLES:
            sage: ?
        """
        return self + right._left_scalar_multiply(-1)



    ###################################################
    # Comparison
    ###################################################
    def __richcmp__(self, right, int op):
        """
        EXAMPLES:
            sage: ???
        """
        # TODO: change the last "==" here to "is"
        if not PY_TYPE_CHECK(right, Matrix) or not PY_TYPE_CHECK(self, Matrix) or not ((<Matrix>right)._parent == (<Matrix>self)._parent):
            # todo: can make faster using the cdef interface to coerce
            return sage.structure.coerce.cmp(self, right)

        cdef int r
        r = (<Matrix>self)._cmp_sibling_cdef(right)

        if op == 0:  #<
            return bool(r  < 0)
        elif op == 2: #==
            return bool(r == 0)
        elif op == 4: #>
            return bool(r  > 0)
        elif op == 1: #<=
            return bool(r <= 0)
        elif op == 3: #!=
            return bool(r != 0)
        elif op == 5: #>=
            return bool(r >= 0)

    cdef _cmp_sibling_cdef(self, Matrix right):
        return cmp(self.list(), right.list())

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
        return self == self.parent()(0)


    #####################################################################################
    # Windowed Strassen Matrix Multiplication and Echelon
    # Precise algorithms invented and implemented by David Harvey and Robert Bradshaw
    # at William Stein's MSRI 2006 Summer Workshop on Modular Forms.
    #####################################################################################
    cdef int _strassen_default_cutoff(self, Matrix right) except -1:
        """
        EXAMPLES:
            sage: ?
        """
        return 128

    def strassen_multiply(self, Matrix right, int cutoff=0):
        """
        EXAMPLES:
            sage: ?
        """
        if self._ncols != right._nrows:
            raise ArithmeticError, "Number of columns of self must equal number of rows of right."
        if not self.base_ring() is right.base_ring():
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

    def strassen_echelon(self, cutoff):
        """
        EXAMPLES:
            sage: ?
        """
        pivots = self._strassen_echelon(self.matrix_window(), cutoff)
        self._set_pivots(pivots)
        return self

    def matrix_window(self, size_t row=0, size_t col=0, size_t nrows=-1, size_t ncols=-1):
        if nrows == -1:
            nrows = self._nrows
            ncols = self._ncols
        return MatrixWindow(self, row, col, nrows, ncols)

    def _strassen_window_multiply(self, MatrixWindow C, MatrixWindow A, MatrixWindow B, int cutoff):
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

            sage.: A_window = A.matrix_window(0, 0, dim1, dim2)
            sage.: B_window = B.matrix_window(0, 0, dim2, dim3)
            sage.: C_window = C.matrix_window(0, 0, dim1, dim3)
            sage.: D_window = D.matrix_window(0, 0, dim1, dim3)

            sage.: strassen_window_multiply(C, A, B, 2)   # use strassen method
            sage.: D_window.set_to_prod(A, B)             # use naive method
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

            sage.: A_window = A.matrix_window(0, 0, dim1, dim2)
            sage.: B_window = B.matrix_window(0, 0, dim2, dim3)
            sage.: C_window = C.matrix_window(0, 0, dim1, dim3)

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

    cdef subtract_strassen_product(result, A, B, int cutoff):
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


    def _strassen_echelon(A, cutoff):
        """
        Compute echelon form, in place.
        Internal function, call with M.echelon(alg="block")
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
            sage: T = A.echelon(alg="gauss")
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

        top_pivots = strassen_echelon(top, cutoff)
        # effectively "multiplied" top row by A^{-1}
        #     [  I  A'B ]
        #     [  C   D  ]

        top_pivot_intervals = int_range(top_pivots)
        top_h = len(top_pivots)

        if top_h == 0:
            #                                 [ 0 0 ]
            # the whole top is a zero matrix, [ C D ]. Run echelon on the bottom
            bottom_pivots = strassen_echelon(bottom, cutoff)
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
                bottom_pivots_rel = strassen_echelon(bottom.matrix_window(0, bottom_start, nrows-split, ncols-bottom_start), cutoff)
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
############################# end class ####################

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

    def set_unsafe(self, size_t i, size_t j, x):
        self._matrix.set_unsafe(i + self._row, j + self._col, x)

    def get_unsafe(self, size_t i, size_t j):
        return self._matrix.get_unsafe(i + self._row, j + self._col)

    def __setitem__(self, ij, x):
        cdef size_t i, j
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
        cdef size_t i, j
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
        cdef size_t i, j
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                self.set_unsafe(i, j, A.get_unsafe(i, j))

    def set_to_zero(MatrixWindow self):
        cdef size_t i, j
        z = self.base_ring()(0)
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                self._matrix.matrix.set_unsafe(i, j, z)

    def add(MatrixWindow self, MatrixWindow A):
        cdef size_t i, j
        if self._nrows != A._nrows or self._ncols != A._ncols:
            raise ArithmeticError, "incompatible dimensions"
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                self.set_unsafe(i, j, self.get_unsafe(i, j) + A.get_unsafe(i, j))

    def subtract(MatrixWindow self, MatrixWindow A):
        cdef size_t i, j
        if self._nrows != A._nrows or self._ncols != A._ncols:
            raise ArithmeticError, "incompatible dimensions"
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                self.set_unsafe(i, j, self.get_unsafe(i, j) - A.get_unsafe(i, j))

    def set_to_sum(MatrixWindow self, MatrixWindow A, MatrixWindow B):
        cdef size_t i, j
        if self._nrows != A._nrows or self._ncols != A._ncols:
            raise ArithmeticError, "incompatible dimensions"
        if self._nrows != B._nrows or self._ncols != B._ncols:
            raise ArithmeticError, "incompatible dimensions"
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                self.set_unsafe(i, j, A.get_unsafe(i, j) + B.get_unsafe(i, j))

    def set_to_diff(MatrixWindow self, MatrixWindow A, MatrixWindow B):
        cdef size_t i, j
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                self.set_unsafe(i, j, A.get_unsafe(i, j) - B.get_unsafe(i, j))

    def set_to_prod(MatrixWindow self, MatrixWindow A, MatrixWindow B):
        cdef size_t i, j, k
        if A._ncols != B._nrows or self._nrows != A._nrows or self._ncols != B._ncols:
            raise ArithmeticError, "incompatible dimensions"
        for i from 0 <= i < A._nrows:
            for j from 0 <= j < B._ncols:
                s = 0
                for k from 0 <= k < A._ncols:
                    s = s + A.get_unsafe(i, k) * B.get_unsafe(k, j)
                self.set_unsafe(i, j, s)

    def add_prod(MatrixWindow self, MatrixWindow A, MatrixWindow B):
        cdef size_t i, j, k
        if A._ncols != B._nrows or self._nrows != A._nrows or self._ncols != B._ncols:
            raise ArithmeticError, "incompatible dimensions"
        for i from 0 <= i < A._nrows:
            for j from 0 <= j < B._ncols:
                s = self.get_unsafe(i, j)
                for k from 0 <= k < A._ncols:
                    s = s + A.get_unsafe(i, k) * B.get_unsafe(k, j)
                self.set_unsafe(i, j, s)

    def subtract_prod(MatrixWindow self, MatrixWindow A, MatrixWindow B):
        cdef size_t i, j, k
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
        calculate the echelon form of this matrix, returning the list of pivot columns
        """
        echelon = self.to_matrix()
        echelon.echelon() # TODO: read only, only need to copy pointers
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

    EXAMPLES:
        sage: ?

    AUTHORS:
      -- Robert Bradshaw

    """
    def __init__(self, indices=None, range=None):
        """
        EXAMPLES:
            sage: ?
        """
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
        """
        EXAMPLES:
            sage: ?
        """
        return str(self._intervals)

    def intervals(self):
        """
        EXAMPLES:
            sage: ?
        """
        return self._intervals

    def to_list(self):
        """
        EXAMPLES:
            sage: ?
        """
        all = []
        for iv in self._intervals:
            for i in range(iv[0], iv[0]+iv[1]):
                all.append(i)
        return all

    def __iter__(self):
        """
        EXAMPLES:
            sage: ?
        """
        return self._intervals.__iter__()

    def __len__(self):
        """
        EXAMPLES:
            sage: ?
        """
        len = 0
        for iv in self._intervals:
            len = len + iv[1]
        return len

    # Yes, these two could be a lot faster...
    # Basically, this class is for abstracting away what I was trying to do by hand in several places
    def __add__(self, right):
        """
        EXAMPLES:
            sage: ?
        """
        all = self.to_list()
        for i in right.to_list():
            all.append(i)
        return int_range(all)

    def __sub__(self, right):
        """
        EXAMPLES:
        sage: ?
        """
        all = self.to_list()
        for i in right.to_list():
            if i in all:
                all.remove(i)
        return int_range(all)

    def __mul__(self, right):
        """
        In the boolean sense, i.e. intersection

        EXAMPLES:
            sage: ?
        """
        intersection = []
        all = self.to_list()
        for i in right.to_list():
            if i in all:
                intersection.append(i)
        return int_range(intersection)

#######################
# Unpickling
#######################

def unpickle(cls, parent, data, version):
    """
    Unpickle a matrix.
    """
    cdef Matrix A
    try:
        A = cls.__new__(cls, parent, 0,0)
        A._parent = parent  # make sure -- __new__ doesn't have to set it, but unpickle may need to know.
        A._nrows = parent.nrows()
        A._ncols = parent.ncols()
        A.unpickle(data, version)
    except Exception, msg:
        print msg
    return A

