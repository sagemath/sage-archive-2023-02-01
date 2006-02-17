"""
Sparse matrices

EXAMPLES:
    sage: M = MatrixSpace(QQ, 2, 3, sparse=True)
    sage: A = M([1,2,3, 1,1,1])
    sage: A
    [1 2 3]
    [1 1 1]
    sage: A.echelon_form()
    [ 1  0 -1]
    [ 0  1  2]

    sage: M = MatrixSpace(QQ, 1000,1000, sparse=True)
    sage: A = M(0)
    sage: A[1,1] = 5


    sage: from sage.matrix.sparse_matrix import SparseMatrix
    sage: x = SparseMatrix(QQ, 5,10)
    sage: x.randomize(5)
    sage: x.echelon_form()       # random output
    [
    1, 0, 0, 0, 0, 0, 0, 0, 0, 10/29,
    0, 1, 0, 0, 0, 0, 0, 0, 0, -4/29,
    0, 0, 1, 0, 0, 0, 0, 0, 0, -12/29,
    0, 0, 0, 0, 0, 0, 1, 0, 0, 24/29,
    0, 0, 0, 0, 0, 0, 0, 0, 1, 4/29
    ]
"""
import random, weakref

import sage.rings.arith as arith
import sage.misc.misc as misc
import sage.matrix.dense_matrix_pyx
import sage.rings.integer_ring as integer_ring
import sage.rings.rational_field as rational_field
import sage.rings.rational as rational
import sage.rings.integer as integer
import sage.rings.rational_field as rational_field
import sage.rings.integer_ring as integer_ring
import sage.rings.ring as ring
import sparse_matrix_pyx

from sage.misc.proof import proof

START_PRIME = 20011  # used for modular algorithms

#################################################################
#
# Generic sparse vectors
#
#################################################################

_objsSVS = {}
class SparseVectorSpace(object):
    def __new__(cls, base_ring, degree):
        key = (base_ring, degree)
        if _objsSVS.has_key(key):
            x = _objsSVS[key]()
            if x != None: return x
        O = object.__new__(Sparse_vector_space_generic)
        _objsSVS[key] = weakref.ref(O)
        return O

class Sparse_vector_space_generic(SparseVectorSpace):
    """
    Sparse vector space.
    """
    def __init__(self, base_ring, degree):
        self.__base_ring = base_ring
        self.__degree = degree

    def __repr__(self):
        return "Ambient sparse vector space of dimension %s over %s."%(
            self.__base_ring, self.__degree)

    def __call__(self, x):
        if (isinstance(x, (int, ring.Integer)) and x == 0) or x in self.__base_ring:
            return SparseVector(self.__base_ring, self.__degree)
        if x in self:
            return x
        raise TypeError, "Cannot coerce %s into %s."%(x, self)

    def __cmp__(self, other):
        if isinstance(other, SparseVectorSpace) and \
           self.degree() == other.degree() and \
           self.base_ring() == other.base_ring():
            return 0
        return -1

    def base_ring(self):
        return self.__base_ring

    def degree(self):
        return self.__degree

    def dimension(self):
        """
        Return the dimension of this vector space.

        EXAMPLES:
            sage: V = VectorSpace(QQ,30,sparse=True)
            sage: V.dimension()
            30
            sage: v = V(0)
            sage: v[12]=2; v[9]=1; v[29]=8
            sage: v
            (0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8)
            sage: W = V.subspace([v, 2*v])
            sage: W.dimension()
            1
        """
        return self.__degree


class SparseVector:
    """
    A generic sparse vector.
    """
    def __init__(self, base_ring, degree, entries=[],
                 coerce=True, sort=True):
        self.__degree = degree
        self.__base_ring = base_ring
        self.set_entries(entries, coerce, sort)

    def set_entries(self, entries, coerce=True, sort=True):
        if not coerce:
            self.__entries = list(entries)  # copy
            if sort:
                self.__entries.sort()
            return
        self.__entries = []
        R = self.base_ring()
        ZERO = R(0)
        for i, x in entries:
            if i < 0 or i >= self.degree():
                raise ArithmeticError, "Invalid index."
            if x != ZERO:
                self.__entries.append((i,R(x)))
        if sort:
            self.__entries.sort()

    def entries(self):
        return self.__entries

    def list(self):
        return self.__entries

    def __repr__(self):
        return str(self.dense_vector())

    def dense_list(self):
        R = self.__base_ring
        v = [R(0) for _ in xrange(self.__degree)]
        for i, x in self.list():
            v[i] = x
        return v

    def dense_vector(self):
        import sage.modules.free_module
        V = sage.modules.free_module.VectorSpace(self.__base_ring,
                                        self.degree(), sparse=False)
        return V(entries = self.dense_list(),
                 coerce_entries=False)

    def degree(self):
        return self.__degree

    def base_ring(self):
        return self.__base_ring

    def __add__(self, right):
        return self.add(right)

    def add(self, right, multiple=1):
        if not isinstance(right, SparseVector):
            raise TypeError
        multiple = self.base_ring()(multiple)
        if self.degree() != right.degree():
            raise ArithmeticError, "Degrees must be the same."
        if self.base_ring() != right.base_ring():
            raise ArithmeticError, "Base rings must be equal."
        v = self.entries()
        w = right.entries()
        z = []
        i = 0  # index into entries of v
        j = 0  # index into entries of w
        while i < len(v) or j < len(w):
            if i >= len(v):   # just copy w in
                z.append((w[j][0],multiple*w[j][1]))
                j += 1
            elif j >= len(w):   # just copy v in
                z.append(v[i])
                i += 1
            elif v[i][0] < w[j][0]:  # copy entry from v in
                z.append(v[i])
                i += 1
            elif v[i][0] > w[j][0]: # copy entry from w in
                s = multiple*w[j][1]
                if s != 0:
                    z.append((w[j][0], s))
                j += 1
            else:                                 # equal, so add and copy
                s = v[i][1] + multiple*w[j][1]
                if s != 0:
                    z.append((v[i][0],s))
                i += 1
                j += 1
            #end if
        # end while
        return SparseVector(self.base_ring(), self.degree(), z, coerce=False)

    def __rmul__(self, left):
        return self.scalar_multiple(left)

    def parent(self):
        return SparseVectorSpace(self.base_ring(), self.degree())

    def scalar_multiple(self, left):
        left = self.base_ring()(left)
        X = []
        for i, x in self.entries():
            X.append((i,left*x))
        return SparseVector(self.base_ring(), self.degree(),
                             X, coerce=False, sort=False)

    def rescale(self, scalar):
        """
        Rescale self by multiplying each entry by scalar.

        NOTE: Sparse vectors are stored as lists of 2-tuples, and
        2-tuples are immutable, so this function is not more efficient
        than simply using scalar*self.
        """
        scalar = self.base_ring()(scalar)
        self.__entries = [(i, scalar*x) for i, x, in self.__entries]

    def __sub__(self, right):
        return self + (-1)*right

    def __neg__(self):
        return (-1)*self

    def randomize(self, sparcity=4):
        entries = []
        R = self.base_ring()
        X = []
        for _ in xrange(sparcity):
            i = random.randrange(self.degree())
            if i in X:
                continue
            X.append(i)
            x = R.random()
            if x != 0:
                entries.append((i, x))
        self.set_entries(entries, coerce=False, sort=True)

    def __getitem__(self, i):
        if i < 0:
            i = self.degree() + i
        for i2, x in self.__entries:
            if i == i2:
                return x
        if i < 0 or i >= self.degree():
            raise IndexError
        return self.base_ring()(0)

    def num_nonzero(self):
        return len(self.__entries)

    def first_nonzero_position(self):
        if len(self.__entries) == 0:
            return self.degree()
        return self.__entries[0][0]

    def first_nonzero_entry(self):
        if len(self.__entries) == 0:
            return self.base_ring()(0)
        return self.__entries[0][1]


#################################################################
#
# Generic sparse matrices
#
#################################################################

class SparseMatrix(object):
    def __new__(cls, base_ring, nrows, ncols, entries=[],
                coerce=True, sort=True, copy=True):
        nrows = int(nrows)
        #if not isinstance(nrows, int):
        #    raise TypeError, "nrows must be an int"
        ncols = int(ncols)
        #if not isinstance(ncols, int):
        #    raise TypeError, "ncols must be an int"
        if not isinstance(entries, list):
            raise TypeError, "entries must be a list"
        if base_ring is rational_field.RationalField():
            return object.__new__(Sparse_matrix_rational)
        return object.__new__(Sparse_matrix_generic)

def SparseMatrix_from_rows(X):
    """
    INPUT:
        X -- nonempty list of SparseVector rows

    OUTPUT:
        Sparse_matrix with those rows.

    EXAMPLES:
        sage: V = VectorSpace(QQ,20,sparse=True)
        sage: v = V(0)
        sage: v[9] = 4
        sage: from sage.matrix.sparse_matrix import *
        sage: SparseMatrix_from_rows([v])
        [
        0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
        ]
        sage: SparseMatrix_from_rows([v, v, v, V(0)])
        [
        0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
        ]
    """
    if not isinstance(X, (list, tuple)):
        raise TypeError, "X (=%s) must be a list or tuple"%X
    if len(X) == 0:
        raise ArithmeticError, "X must be nonempty"
    entries = []
    R = X[0].base_ring()
    ncols = X[0].degree()
    for i in range(len(X)):
        for j, x in X[i].entries().iteritems():
            entries.append((i,j,x))
    return SparseMatrix(R, len(X), ncols, entries, coerce=False, sort=False)

class Sparse_matrix_generic(SparseMatrix):
    """
    A generic sparse matrix.
    """
    def __init__(self, base_ring, nrows, ncols, entries=[],
                 coerce=True, sort=True, copy=True):
        self.__nrows = nrows
        self.__ncols = ncols
        self.__base_ring = base_ring
        self.set_entries(entries, coerce, sort, copy)

    def __repr__(self):
        return str(self.dense_matrix())

    def copy(self):
        E = list(self.entries())
        return SparseMatrix(self.__base_ring, self.__nrows, \
                self.__ncols, E, coerce=False, sort=False, copy=False)

    def __add__(self, other):
        if not isinstance(other, Sparse_matrix_generic):
            raise TypeError
        if self.nrows() != other.nrows() or self.ncols() != other.ncols():
            raise ArithmeticError, "Matrices must have the same number of rows and columns."
        if self.base_ring() != other.base_ring():
            raise ArithmeticError, "Matrices must have the same base ring."

        # Compute the sum of two sparse matrices.
        # This is complicated because of how we represent sparse matrices as tuples (i,j,x).
        # Algorithm:
        #   1. Sort both entry lists.
        #   2. March through building a new list, adding when the two i,j are the same.
        v = self.__entries
        v.sort()
        w = other.__entries
        w.sort()
        s = []
        i = 0  # pointer into self
        j = 0  # pointer into other
        while i < len(v) and j < len(w):
            vi = (v[i][0], v[i][1])
            wj = (w[j][0], w[j][1])
            if vi < wj:
                s.append(tuple(v[i]))
                i += 1
            elif vi > wj:
                s.append(tuple(w[j]))
                j += 1
            else:  # equal
                s.append((vi[0],vi[1],v[i][2] + w[j][2]))
                i += 1
                j += 1
        if i < len(v):
            while i < len(v):
                s.append(tuple(v[i]))
                i += 1
        if j < len(w):
            while j < len(w):
                s.append(tuple(w[j]))
                j += 1
        A = SparseMatrix(self.base_ring(), self.nrows(),
                         self.ncols(), s, coerce=False)
        return A

    def __sub__(self, right):
        return self + (-1)*right

    def __neg__(self):
        return (-1)*self

    def __getitem__(self, i):
        return self.row(i)

    def row(self, i):
        """
        Return the ith row of this matrix as a sparse vector.

        WARNING: Sparse matrices are stored as a single list of triples (i,j,x),
        so extracting the i-th row is expensive.  This command builds a redundant
        representation of the matrix as a list of sparse vectors, thus doubling
        the memory requirement.
        """
        i = int(i)
        #if not isinstance(i, int):
        #    raise TypeError, "i must be an int"
        if i < 0: i = self.nrows() + i
        return self.rows()[i]

    def rows(self):
        """
        Return a list of the sparse rows of this matrix.  The ith
        element of this list is the i-th row of this matrix.
        """
        try:
            return self.__rows
        except AttributeError:
            X = []
            entries = []
            k = 0
            R = self.base_ring()
            n = self.ncols()
            for i, j, x in self.entries():
                if i > k:
                    while len(X) <= k-1:
                        X.append(SparseVector(R, n))
                    X.append(SparseVector(R, n, entries, coerce=False, sort=False))
                    entries = []
                    k = i
                entries.append((j,x))
            while len(X) <= k-1:
                X.append(SparseVector(R, n))
            X.append(SparseVector(R, n, entries, coerce=False, sort=False))
            while len(X) < self.nrows():
                X.append(SparseVector(R, n))
            self.__rows = X
        return self.__rows

    def __getslice__(self, i, j):
        i = int(i)
        #if not  isinstance(i, int):
        #    raise TypeError, "i must be an int"
        j = int(j)
        #if not  isinstance(j, int):
        #    raise TypeError, "j must be an int"
        if i < 0: i = self.nrows() + i
        if j < 0: j = self.nrows() + j
        if i >= self.nrows():
            i = self.nrows() - 1
        if j >= self.nrows():
            j = self.nrows() - 1
        return self.matrix_from_rows(xrange(i,j))

    def nrows(self):
        return self.__nrows

    def ncols(self):
        return self.__ncols

    def base_ring(self):
        return self.__base_ring

    def entries(self):
        return self.__entries

    def list(self):
        return self.__entries

    def dict(self):
        X = {}
        for i, j, x in self.entries():
            X[(i,j)] = x
        return X

    def dense_matrix(self):
        """
        todo
        """
        import sage.matrix.matrix
        M = sage.matrix.matrix.MatrixSpace(self.base_ring(),
                                     self.nrows(),
                                     self.ncols(),
                                     sparse = False)(0)
        for i, j, x in self.list():
            M[i,j] = x
        return M

    def nonpivots(self):
        # We convert the pivots to a set so we have VERY FAST
        # inclusion testing
        X = set(self.pivots())
        return [j for j in xrange(self.ncols()) if not (j in X)]

    def pivots(self):
        try:
            return self.__pivots
        except AttributeError:
            self.echelon_form()
        return self.__pivots

    def _set_pivots(self, X):
        self.__pivots = X

    def matrix_from_nonpivot_columns(self):
        """
        The sparse matrix got by deleted all pivot columns.
        """
        return self.matrix_from_columns(self.nonpivots())

    def matrix_from_columns(self, columns):
        """
        Returns the sparse submatrix of self composed of the given
        list of columns.

        INPUT:
            columns -- list of int's
        OUTPUT:
            a sparse matrix.
        """
        # Let A be this matrix and let B denote this matrix with
        # the columns deleted.
        # ALGORITHM:
        # 1. Make a table that encodes the function
        #          f : Z --> Z,
        #    f(j) = column of B where the j-th column of A goes
        # 2. Build new list of entries and return resulting matrix.
        C = set(columns)
        X = []
        j = 0
        for i in xrange(self.ncols()):
            if i in C:
                X.append(j)
                j += 1
            else:
                X.append(-1)    # column to be deleted.
        entries2 = []
        for i, j, x in self.entries():
            if j in C:
                entries2.append((i,X[j],x))
        return SparseMatrix(self.base_ring(), self.nrows(),
                            len(C), entries2, coerce=False, sort=False)

    def matrix_from_rows(self, rows):
        """
        Returns the sparse submatrix of self composed of the given
        list of rows.

        INPUT:
            rows -- list of int's
        OUTPUT:
            a sparse matrix.
        """
        R = set(rows)
        if not R.issubset(set(xrange(self.nrows()))):
            raise ArithmeticError, "Invalid rows."
        X = []
        i = 0
        for j in xrange(self.nrows()):
            if j in R:
                X.append(i)
                i += 1
            else:
                X.append(-1)    # row to be deleted.
        entries2 = []
        for i, j, x in self.entries():
            if i in R:
                entries2.append((X[i],j,x))
        return SparseMatrix(self.base_ring(), len(R),
                            self.ncols(), entries2, coerce=False, sort=False)


    def set_entries(self, entries, coerce=True, sort=True, copy=True):
        """
        entries is a list of triples (i,j,x) and the
        x must be nonzero.

        This function does *not* check that i and j are in bounds
        or that the x are all nonzero.
        """
        try:
            del self.__rows
        except AttributeError:
            pass
        try:
            del self.__echelon_form
        except AttributeError:
            pass
        if not coerce:
            if copy:
                self.__entries = list(entries)   # copy
            else:
                self.__entries = entries
            if sort:
                self.__entries.sort()
            return
        self.__entries = []
        R = self.base_ring()
        self.__entries = [(i,j,R(z)) for i,j,z in entries]
        if sort:
            self.__entries.sort()

    def randomize(self, sparcity=4, exact=True):
        entries = []
        R = self.base_ring()
        for i in range(self.nrows()):
            if exact:
                r = sparcity
            else:
                r = random.randrange(sparcity)
            X = []
            for j in range(0,r+1):
                x = R.random()
                if x != R(0):
                    k = random.randrange(0,self.ncols())
                    if not (k in X):
                        entries.append((i,k, x))
                        X.append(k)
        self.set_entries(entries, coerce=False, sort=True)

    def transpose(self):
        """
        Returns the transpose of self.
        """
        entries2 = [(j,i,x) for i,j,x in self.entries()]
        return SparseMatrix(self.base_ring(), self.ncols(),
                            self.nrows(), entries2, coerce=False, sort=True)

    def __rmul__(self, left):
        return self.scalar_multiple(left)

    def scalar_multiple(self, left):
        R = self.base_ring()
        left = R(left)
        if left == R(1):
            return self
        if left == R(0):
            return SparseMatrix(R, self.nrows(), self.ncols(), coerce=False, sort=False)

        X = []
        for i, j, x in self.list():
            X.append((i,j,left*x))
        return SparseMatrix(self.base_ring(), self.nrows(),
                            self.ncols(), X, coerce=False, sort=False)

    def echelon_form(self, params=None):
        """
        Returns the echelon form of this matrix.

        INPUT:
           params -- ignored.
        """
        # ALGORITHM:
        # Since we know nothing about the base field, we use a generic
        # algorithm.  Since sparse matrices are stored as triples
        # (i,j,x), which is not a suitable format for row operations,
        # we first convert to a list of sparse rows, then directly
        # perform a generic echelon algorithm on that list of rows.
        try:
            return self.__echelon_form
        except AttributeError:
            pass
        t = misc.verbose("Started generic sparse echelon.")
        K = self.base_ring()
        ONE = K(1)
        if not K.is_field():
            raise ArithmeticError, "The base ring must be a field."
        X = self.rows()
        nrows = self.nrows()
        ncols = self.ncols()
        pivot_positions = []
        start_row = 0
        nrows = self.nrows()
        ncols = self.ncols()
        for c in range(ncols):
            N = [(X[r].num_nonzero(),r) for r in xrange(start_row, nrows) \
                 if X[r].first_nonzero_position() == c]
            if len(N) == 0:
                continue
            N.sort()
            r = N[0][1]
            leading = X[r].first_nonzero_entry()
            if leading != 0:
                pivot_positions.append(c)
                # 1. Rescale
                X[r].rescale(ONE/leading)
                # 2. Swap
                X[r], X[start_row] = X[start_row], X[r]
                # 3. Clear column
                for i in range(nrows):
                    if i != start_row:
                        s = X[i][c]
                        if s != 0:
                            X[i] = X[i].add(X[start_row], -s)
            # endif
            start_row += 1
        #endfor
        self.__pivots = pivot_positions
        E = SparseMatrix_from_rows(X)
        E.__pivots = pivot_positions
        self.__echelon_form = E
        misc.verbose("Finished generic echelon.",t)
        return E


#################################################################
#
# Sparse matrices over the rational numbers.
#
# This is a kick-ass tool for computing echelon forms of
# very sparse matrices using a multi-modular method.
# E.g., on my 1.8Ghz Pentium-M laptop, it can compute
# the reduced row-echelon form of a 50000x100000 with
# up to 3 nonzero entries per row (between -2 and 2)
# in under 10 minutes.  And this computation uses only
# about 50MB RAM!
#
# sage: from sage.matrix.sparse_matrix import SparseMatrix
# sage: x = SparseMatrix(Q, 50000,100000)
# sage: x.randomize(3)
# sage: set_verbose(2)
# sage: time v=x.echelon_form(height_guess=1)
# echelon_modular: p=20011
# echelon_modular: p=20021
# echelon_modular: start rr
# echelon_modular: done (time = 3.08)
# Time: 541.28 seconds
#
# With proof=False, it takes only 7 minutes (on slower 1.6Ghz laptop):
# sage: time v=x.echelon_form(height_guess=1)
# echelon_form: height_guess = 1
# echelon_form: p=20011
# echelon_form: time to reduce: (time = 0.6)
# echelon_form: time to echelon: (time = 418.94)
# echelon_form: time for pivot compare (time = 0.0)
# echelon_form: start crt and rr
# echelon_form: crt time is (time = 0.49)
# echelon_form: rr time is (time = 0.27)
# echelon_form: Not checking validity of result (since proof=False).
# echelon_form: total time (time = 420.3)
# Time: 421.0 seconds
#
# On this same slower laptop with proof=False, it takes 30 minutes
# to do a random 100000x200000 matrix.
#
#################################################################

class Sparse_matrix_rational(Sparse_matrix_generic):
    """
    A sparse matrix over the rational numbers.

    Once you create the matrix you can not change specific entries.
    However you can compute the echelon form, which uses a sparse
    multi-modular method.
     Create with
        Sparse_matrix_rational(base_ring, nrows, ncols, entries)
     INPUT:
        nrows -- integer, the number of rows
        ncols -- integer, the number of columns
        entries-- list of tuples (i,j,x), where 0 <= i < nrows, 0
                  <= j < ncols, and x is a Python object that has a
                  string representation that can be parsed as a gmp mpq.
                  All the x's are assumed *nonzero*!
    """
    def __init__(self, base_ring, nrows, ncols, entries=[],
                 coerce=True, sort=True, copy=True):
        Sparse_matrix_generic.__init__(self, base_ring,
                                       nrows, ncols, entries,
                                       coerce, sort, copy)

    def dense_matrix(self):
        A = sage.matrix.dense_matrix_pyx.Matrix_rational(self.nrows(), self.ncols())
        for i, j, x in self.list():
            A[i,j] = x
        return A

    def denom(self):
        """
        Return the least common multiple of the denominators of
        the entries of self.
        """
        d = integer.Integer(1)
        for _, _, x in self.list():
            d = d.lcm(x.denom())
        return d

    def clear_denom(self):
        """
        Replace self by n*self, where n is the least common multiple
        of the denominators of all entries of self.
        """
        d = self.denom()
        if d == integer.Integer(1):
            return
        d = rational.Rational(d)
        entries = self.list()
        for k in range(len(entries)):
            i, j, x = entries[k]
            entries[k] = (i,j,x*d)

    def echelon_form(self, height_guess=None):
        """
        Returns echelon form of self, possibly modifying self by rescaling.
        Uses a sparse multi-modular method.
        ALGORITHM:
        The following is a modular algorithm for computing the echelon
        form.  Define the height of a matrix to be the max of the
        absolute values of the entries.
        Input: Matrix A with n columns (this).
         0. Rescale input matrix A to have integer entries.  This does
            not change echelon form and makes reduction modulo lots of
            primes significantly easier if there were denominators.
            Henceforth we assume A has integer entries.
         1. Let c be a guess for the height of the echelon form.  E.g.,
            c=1000, since matrix is sparse and application is modular symbols.
         2. Let M = n * c * H(A) + 1,
            where n is the number of columns of A.
         3. List primes p_1, p_2, ..., such that the product of
            the p_i is at least M.
         4. Try to compute the rational reconstruction CRT echelon form
            of A mod the product of the p_i.  If rational
            reconstruction fails, compute 1 more echelon forms mod the
            next prime, and attempt again.  Make sure to keep the
            result of CRT on the primes from before, so we don't have
            to do that computation again.  Let E be this matrix.
         5. Compute the denominator d of E.
            Try to prove that result is correct by checking that
                  H(d*E)*ncols(A)*H(A) < (prod of reduction primes)
            where H denotes the height.   If this fails, do step 4 with
            a few more primes
        """
        try:
            # If self.__echelon_form is true, this means the matrix self is already
            # in echelon form, so we return self.
            # Otherwise, we return what is stored in the __echelon_form variable.
            if isinstance(self.__echelon_form, bool) and self.__echelon_form:
                return self
            return self.__echelon_form
        except AttributeError:
            pass
        if self.nrows() == 0 or self.ncols() == 0:
            return self
        self.clear_denom()
        hA = long(self.height())
        if height_guess == None:
            height_guess = 100000*hA**4
        tm = misc.verbose("height_guess = %s"%height_guess, level=2)
        if proof:
            M = self.ncols() * height_guess * hA  +  1
        else:
            M = height_guess + 1
        p = START_PRIME
        X = []
        best_pivots = []
        prod = 1
        problem = 0
        while True:
            while prod < M:
                problem += 1
                if problem > 20:
                    misc.verbose("stuck in sparse_matrix reduce")
                t = misc.verbose("p=%s"%p, level=2)
                A = self.matrix_modint(p)
                t = misc.verbose("time to reduce matrix mod p:",t, level=2)
                A.echelon()
                t = misc.verbose("time to put reduced matrix in echelon form:",t, level=2)
                c = sage.matrix.dense_matrix_pyx.cmp_pivots(best_pivots, A.pivots())
                if c <= 0:
                    best_pivots = A.pivots()
                    X.append(A)
                    prod = prod * p
                else:
                    # do not save A since it is bad.
                    if misc.LEVEL > 1:
                        misc.verbose("Excluding this prime (bad pivots).")
                p = arith.next_prime(p)
                t = misc.verbose("time for pivot compare", t, level=2)
            # Find set of best matrices.
            Y = []
            # recompute product, since may drop bad matrices
            prod = 1
            for i in range(len(X)):
                if sage.matrix.dense_matrix_pyx.cmp_pivots(
                                   best_pivots, X[i].pivots()) <= 0:
                    Y.append(X[i])
                    prod = prod * X[i].prime()
            try:
                t = misc.verbose("start crt and rr", level=2)
                A = SparseMatrix_using_crt(Y)
                t = misc.verbose("crt time is", t, level=2)
                E = rational_reconstruction(A, integer.Integer(prod))
                misc.verbose("rr time is",t, level=2)
            except ValueError:
                misc.verbose("Redoing with several more primes", level=2)
                for i in range(3):
                    M = M * START_PRIME
                continue
            if not proof:
                misc.verbose("Not checking validity of result (since proof=False).", level=2)
                break
            d = E.denom()
            dE = E.scalar_multiple(d)
            hdE = long(dE.height())
            if hdE * self.ncols() * hA < prod:
                break
            for i in range(3):
                M = M * START_PRIME
        #end while
        misc.verbose("total time",tm, level=2)
        self._set_pivots(best_pivots)
        E._set_pivots(best_pivots)
        E.__echelon_form = True
        self.__echelon_form = E
        return E

    def height(self):
        """
        Returns the height of self, which is the maximum of the absolute
        values of all numerators and denominators of the elements of self.
         Since 0 = 0/1 has denominator 1, the height is at least 1.
        """
        h = integer.Integer(1)
        for _, _, x in self.list():
            a = x.height()
            if a > h:
                h = a
        return h

    def matrix_modint(self, n):
        X = []
        for i, j, x in self.list():
            a = x.mod_ui(n)
            if a:
                X.append((i,j,a))
        return sparse_matrix_pyx.Matrix_modint(n, self.nrows(), self.ncols(), X)

    def randomize(self, sparcity=4, bound=2, bound_denom=2, exact=False):
        """
        The sparcity is a bound on the number of nonzeros per row.
        """
        try:
            del self.__pivots   # since randomizing changes the pivots.
        except AttributeError:
            pass
        entries = []
        for i in range(self.nrows()):
            if exact:
                r = sparcity
            else:
                r = random.randrange(sparcity)
            X = []
            for j in range(0,r+1):
                x = rational.Rational("%s/%s"%(
                    random.randrange(-bound,bound), random.randrange(1,bound_denom)))
                if x != rational.Rational(0):
                    k = random.randrange(0,self.ncols())
                    if not (k in X):
                        entries.append((i,k, x))
                        X.append(k)
        self.set_entries(entries, coerce=False, sort=False)

    def set_entries(self, entries, coerce=True, sort=True, copy=True):
        try:
            del self.__echelon_form
        except AttributeError:
            pass
        Sparse_matrix_generic.set_entries(self, entries, coerce, sort, copy)


def crt_dict(X, M, M2):
    """
    Given a list X of dictionaries and a corresponding list M of
    Integer moduli (of the same length), use the CRT to make
    a single dict.

    M2 - should be the partial products of elements of M, so
         M2[i] = M[0] * ... * M[i-1],
         and M2[0] = 1.
         This is an optimization.
    """
    if len(X) != len(M):
        raise RuntimeError, "The number of dictionaries must be the same as the number of moduli."
    n = len(X)
    if n == 0:
        return []
    if n == 1:
        Y = []
        for x in X[0].keys():
            Y.append((x,X[0][x]))
        return Y

    # 1. take the set-theoretic union of the keys
    K = set(X[0].keys())
    for i in range(1,n):
        K.union(set(X[i].keys()))

    Y = []
    zero = integer.Integer(0)
    one = integer.Integer(1)
    # 3. for each key in the union, set up and solve a crt
    for x in K:
        a = zero
        for i in range(n):
            A = X[i]
            if A.has_key(x):
                b = A[x]
            else:
                b = zero
            a = a.crt(b, M2[i], M[i])
        Y.append((x, a))
    return Y

# TODO: The following function is *really* slow, and is the
# bottleneck in a lot of computations.  It must be moved to
# the right place (probably sparse_matrix_pyx).  It could
# be very fast if done right!
def SparseMatrix_using_crt(X):
    """
    Given Matrix_modint's X, this function uses the chinese remainder
    theorem to create a SparseMatrix over the integers.
    INPUT:
       X -- a list of sparse_matrix_pyx.Matrix_modint's modulo coprime moduli.
    OUTPUT:
       A single SparseMatrix with integer entries that reduces to all
       the X[i].
    """
    if not isinstance(X,list) or len(X) == 0:
        raise TypeError

    P = [integer.Integer(A.prime()) for A in X]
    B = integer_ring.crt_basis(P)

    # Now we compute the linear combination
    #    sum B[i] * X[i].
    A = SparseMatrix(integer_ring.IntegerRing(),
                     X[0].nrows(), X[0].ncols(), X[0].list(),
                     coerce=True, sort=False)

    A = A.scalar_multiple(B[0])
    for i in range(1,len(X)):
        V = SparseMatrix(integer_ring.IntegerRing(),
                X[i].nrows(), X[i].ncols(), X[i].list(),
                coerce=True, sort=False)
        A += V.scalar_multiple(B[i])
    return A

def rational_reconstruction(A, m, denom_optimization=False):
    """
    Lift A mod m using rational reconstruction.

    Given a sparse matrix A over the integers and an integer m, this
    function uses rational reconstruction element-by-element to try
    and find a matrix B over the rational numbers that reduces
    to A modulo m.  The entries of B will be uniquely determined by
    the condition that the numerator and denominator have absolute
    value at most sqrt(m/2).  If no such B exists, this function raises
    a ValueError.

    INPUT:
        A -- SparseMatrix over the integers
        m -- an integer

    OUTPUT:
        B -- SparseMatrix over the rationals
    """
    # NOTE: There is a "denominator optimization that I should try, namely
    # to multiply entry under consideration by denominators so far, and try
    # lifting with denom 1 before trying rr.  ...
    if denom_optimization:
        pass
    else:
        entries = [(i, j, rational.rational_reconstruction(x,m))
                   for i, j, x in A.entries()]
    return SparseMatrix(rational_field.RationalField(), A.nrows(),
                        A.ncols(), entries, coerce=False, sort=False)


