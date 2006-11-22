r"""
Spaces of matrices over a ring or field

You can create any space $\text{Mat}_{n\times m}(R)$ of either dense
or sparse matrices with given number of rows and columns over any
commutative or noncommutative ring.

EXAMPLES:
    sage: MS = MatrixSpace(QQ,6,6,sparse=True); MS
    Full MatrixSpace of 6 by 6 sparse matrices over Rational Field
    sage: MS.base_ring()
    Rational Field
    sage: MS = MatrixSpace(ZZ,3,5,sparse=False); MS
    Full MatrixSpace of 3 by 5 dense matrices over Integer Ring
"""

# System imports
import random
import weakref

# SAGE matrix imports
import matrix
import matrix_generic_dense
import matrix_generic_sparse

## import matrix_domain_dense
## import matrix_domain_sparse

## import matrix_pid_dense
## import matrix_pid_sparse

## import matrix_field_dense
## import matrix_field_sparse

import matrix_modn_dense
import matrix_modn_sparse

import matrix_integer_dense
## import matrix_integer_sparse

import matrix_rational_dense
## import matrix_rational_sparse

## import matrix_cyclo_dense
## import matrix_cyclo_sparse


# SAGE imports
import sage.structure.parent_gens as parent_gens
import sage.rings.ring as ring
import sage.rings.rational_field as rational_field
import sage.rings.integer_ring as integer_ring
import sage.rings.integer as integer
import sage.rings.field as field
import sage.rings.finite_field as finite_field
import sage.rings.principal_ideal_domain as principal_ideal_domain
import sage.rings.integral_domain as integral_domain
import sage.rings.number_field.all
import sage.rings.integer_mod_ring
import sage.misc.latex as latex
from sage.misc.misc import xsrange

import sage.modules.free_module_element

from sage.structure.sequence import Sequence

def is_MatrixSpace(x):
    """
    returns true if self is an instance of MatrixSpace
    returns false if self is not an instance of MatrixSpace

    EXAMPLES:
	sage: MS = MatrixSpace(QQ,2)
	sage: A = MS.random_element()
	sage: is_MatrixSpace(MS)
	True
	sage: is_MatrixSpace(A)
	False
	sage: is_MatrixSpace(5)
	False
    """

    return isinstance(x, MatrixSpace_generic)

_cache = {}
def MatrixSpace(base_ring, nrows, ncols=None, sparse=False):
    """
    Create with the command

          MatrixSpace(base_ring , nrows [, ncols] [, sparse])

    The default value of the optional argument sparse is False.
    The default value of the optional argument ncols is nrows.

    INPUT:
         base_ring -- a ring
         nrows -- int, the number of rows
         ncols -- (default nrows) int, the number of columns
         sparse -- (default false) whether or not matrices are given
                   a sparse representation
    OUTPUT:
        The unique space of all nrows x ncols matrices over base_ring.

    EXAMPLES:
        sage: MS = MatrixSpace(RationalField(),2)
        sage: MS.base_ring()
        Rational Field
        sage: MS.dimension()
        4
        sage: MS.dims()
        (2, 2)
        sage: B = MS.basis()
        sage: B
        [
        [1 0]
        [0 0],
        [0 1]
        [0 0],
        [0 0]
        [1 0],
        [0 0]
        [0 1]
        ]
        sage: B[0]
        [1 0]
        [0 0]
        sage: B[1]
        [0 1]
        [0 0]
        sage: B[2]
        [0 0]
        [1 0]
        sage: B[3]
        [0 0]
        [0 1]
        sage: A = MS.matrix([1,2,3,4])
        sage: A
        [1 2]
        [3 4]
        sage: MS2 = MatrixSpace(RationalField(),2,3)
        sage: B = MS2.matrix([1,2,3,4,5,6])
        sage: A*B
        [ 9 12 15]
        [19 26 33]

        sage: M = MatrixSpace(ZZ, 10)
        sage: M
        Full MatrixSpace of 10 by 10 dense matrices over Integer Ring
        sage: loads(M.dumps()) == M
        True

    EXAMPLES OF EACH TYPE OF MATRIX SPACE:
    In this section, we create sparse and dense matrix spaces and
    matrices in those spaces for a wide range of base rings.

    Dense and sparse generic matrices:
        sage: ??

    Dense and sparse matrices over a domain:

    Dense and sparse matrices over a principal ideal domain:

    Dense and sparse matrices over a field:

    Dense and sparse matrices over the rational numbers:

    Dense and sparse matrices over a cyclotomic field:

    Dense and sparse matrices over the ring ZZ of integers:

    Dense and sparse matrices over the ring of integers:

    Dense and sparse matrices modulo $n$:

    """
    if ncols is None: ncols = nrows
    nrows = int(nrows); ncols = int(ncols); sparse=bool(sparse)
    key = (base_ring, nrows, ncols, sparse)
    if _cache.has_key(key):
        M = _cache[key]()
        if not M is None: return M

    if not sage.rings.ring.is_Ring(base_ring):
        raise TypeError, "base_ring (=%s) must be a ring"%base_ring

    M = MatrixSpace_generic(base_ring, nrows, ncols, sparse)

    _cache[key] = weakref.ref(M)
    return M



class MatrixSpace_generic(parent_gens.ParentWithGens):
    """
    The space of all nrows x ncols matrices over base_ring.

    EXAMPLES:
        sage: MatrixSpace(ZZ,10,5)
        Full MatrixSpace of 10 by 5 dense matrices over Integer Ring
        sage: MatrixSpace(ZZ,10,2^33)
        Traceback (most recent call last):                                   # 32-bit
        ...                                                                  # 32-bit
        ValueError: number of rows and columns must be less than 2^32 (on a 32-bit computer -- use a 64-bit computer for bigger matrices)    # 32-bit
        Full MatrixSpace of 10 by 8589934592 dense matrices over Integer Ring   # 64-bit
    """
    def __init__(self,  base_ring,
                        nrows,
                        ncols=None,
                        sparse=False):
        parent_gens.ParentWithGens.__init__(self, base_ring)
        if not isinstance(base_ring, ring.Ring):
            raise TypeError, "base_ring must be a ring"
        if ncols == None: ncols = nrows
        nrows = int(nrows)
        #if not isinstance(nrows, int):
        #    raise TypeError, "nrows must be an int"
        ncols = int(ncols)
        #if not isinstance(ncols, int):
        #    raise TypeError, "ncols must be an int"
        if nrows < 0:
            raise ArithmeticError, "nrows must be nonnegative"
        if ncols < 0:
            raise ArithmeticError, "ncols must be nonnegative"

        if sage.misc.misc.is_64bit():
            if nrows >= 2**64 or ncols >= 2**64:
                raise ValueError, "number of rows and columns must be less than 2^64"
        else:
            if nrows >= 2**32 or ncols >= 2**32:
                raise ValueError, "number of rows and columns must be less than 2^32 (on a 32-bit computer -- use a 64-bit computer for bigger matrices)"

        self.__nrows = nrows
        self.__is_sparse = sparse
        if ncols == None:
            self.__ncols = nrows
        else:
            self.__ncols = ncols
        self.__matrix_class = self._get_matrix_class()

    def __call__(self, entries=0, coerce=True, copy=True):
        if isinstance(entries, list) and len(entries) > 0 and \
           sage.modules.free_module_element.is_FreeModuleElement(entries[0]):
            if self.__is_sparse:
                e = {}
                for i in xrange(len(entries)):
                    for j, x in entries[i].iteritems():
                        e[(i,j)] = x
                entries = e
            else:
                entries = sum([v.list() for v in entries],[])
        if not self.__is_sparse and isinstance(entries, dict):
            entries = dict_to_list(entries, self.__nrows, self.__ncols)
            coerce = True
            copy = False
        elif self.__is_sparse and isinstance(entries, (list, tuple)):
            entries = list_to_dict(entries, self.__nrows, self.__ncols)
            coerce = True
            copy = False
        return self.matrix(entries, copy=copy, coerce=coerce)

    def base_extend(self, R):
        """
        INPUT:
            R -- ring
        """
        return MatrixSpace(R, self.__nrows, self.__ncols, self.__is_sparse)

    def _coerce_impl(self, x):
        """
        EXAMPLES:
            sage: MS1 = MatrixSpace(QQ,3)
            sage: MS2 = MatrixSpace(ZZ,4,5,true)
            sage: A = MS1(range(9))
            sage: D = MS2(range(20))
            sage: MS1._coerce_(A)
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: MS2._coerce_(D)
            [ 0  1  2  3  4]
            [ 5  6  7  8  9]
            [10 11 12 13 14]
            [15 16 17 18 19]
        """
        if isinstance(x, matrix.Matrix):
            if self.is_sparse() and x.is_dense():
                raise TypeError, "cannot coerce dense matrix into sparse space for arithmetic"
            if x.nrows() == self.nrows() and x.ncols() == self.ncols():
                if self.base_ring().has_coerce_map_from(x.base_ring()):
                    return self(x)
                raise TypeError, "no canonical coercion"
        return self._coerce_try(x, self.base_ring())

    def __cmp__(self, other):
        """
        Compare this matrix space with other.

        If other is not a matrix space, return something arbitrary but
        deterministic.  Otherwise, compare based on base ring, then on
        number of rows and columns.

        EXAMPLES:

        """
        if isinstance(other, MatrixSpace_generic):
            return cmp((self.base_ring(), self.__nrows, self.__ncols),
                       (other.base_ring(), other.__nrows, other.__ncols))
        return cmp(type(self), type(other))

    def _repr_(self):
        """
        Returns the string representation of a MatrixSpace

        EXAMPLES:
    	sage: MS = MatrixSpace(ZZ,2,4,true)
    	sage: repr(MS)
    	'Full MatrixSpace of 2 by 4 sparse matrices over Integer Ring'
    	sage: MS
    	Full MatrixSpace of 2 by 4 sparse matrices over Integer Ring

        """
        if self.is_sparse():
            s = "sparse"
        else:
            s = "dense"
        return "Full MatrixSpace of %s by %s %s matrices over %s"%(
                    self.__nrows, self.__ncols, s, self.base_ring())

    def _latex_(self):
        r"""
        Returns the latex representation of a MatrixSpace

        EXAMPLES:
    	sage: MS3 = MatrixSpace(QQ,6,6,true)
    	sage: latex(MS3)
        \mbox{\rm Mat}_{6\times 6}(\mathbf{Q})
        """
        return "\\mbox{\\rm Mat}_{%s\\times %s}(%s)"%(self.nrows(), self.ncols(),
                                                      latex.latex(self.base_ring()))

    def _get_matrix_class(self):
        """
        Returns the class of self

        EXAMPLES:
            sage: MS1 = MatrixSpace(QQ,4)
            sage: MS2 = MatrixSpace(ZZ,4,5,true)
            sage: MS1._get_matrix_class()
            <type 'sage.matrix.matrix_rational_dense.Matrix_rational_dense'>
            sage: MS2._get_matrix_class()
            <type 'sage.matrix.matrix_generic_sparse.Matrix_generic_sparse'>
        """
        R = self.base_ring()
        if self.is_dense():
            if sage.rings.integer_ring.is_IntegerRing(R):
                return matrix_integer_dense.Matrix_integer_dense
            elif sage.rings.rational_field.is_RationalField(R):
                return matrix_rational_dense.Matrix_rational_dense
            elif sage.rings.integer_mod_ring.is_IntegerModRing(R) and R.order() < 46340:
                return matrix_modn_dense.Matrix_modn_dense
            # the default
            return matrix_generic_dense.Matrix_generic_dense

        else:
            if sage.rings.integer_mod_ring.is_IntegerModRing(R) and R.order() < 46340:
                return matrix_modn_sparse.Matrix_modn_sparse
            # the default
            return matrix_generic_sparse.Matrix_generic_sparse


    def basis(self):
        try:
            return self.__basis
        except AttributeError:
            v = [self() for _ in range(self.dimension())]
            one = self.base_ring()(1)
            i = 0
            for r in range(self.__nrows):
                for c in range(self.__ncols):
                    v[i][r,c] = one
                    v[i].set_immutable()
                    i += 1
            B = Sequence(v, universe=self, check=False, immutable=True, cr=True)
            self.__basis = B
            return B

    def dimension(self):
        """
        Returns (m rows) * (n cols) of self as Integer

        EXAMPLES:
     	sage: MS = MatrixSpace(ZZ,4,6)
    	sage: u = MS.dimension()
    	sage: u - 24 == 0
    	True
        """
        return self.__nrows * self.__ncols

    def dims(self):
        """
        Returns (m row, n col) representation of self dimension

        EXAMPLES:
    	sage: MS = MatrixSpace(ZZ,4,6)
    	sage: MS.dims()
    	(4, 6)
        """
        return (self.__nrows, self.__ncols)

    def identity_matrix(self):
        """
        Create an identity matrix in self.  (Must be a space of square matrices).

        EXAMPLES:
    	sage: MS1 = MatrixSpace(ZZ,4)
    	sage: MS2 = MatrixSpace(QQ,3,4)
    	sage: I = MS1.identity_matrix()
    	sage: I
    	[1 0 0 0]
    	[0 1 0 0]
    	[0 0 1 0]
    	[0 0 0 1]
    	sage: Er = MS2.identity_matrix()
        Traceback (most recent call last):
        ...
    	TypeError: self must be a space of square matrices
        """

        if self.__nrows != self.__ncols:
            raise TypeError, "self must be a space of square matrices"
        A = self(0)
        for i in xrange(self.__nrows):
            A[i,i] = 1
        return A

    def is_dense(self):
        """
        Returns true if self is dense
        Returns false if self is sparse
        """
        return not self.__is_sparse

    def is_sparse(self):
        """
        Returns True if self is sparse
        Returns False if self is dense
        """
        return self.__is_sparse

    def gen(self, n):
        return self.basis()[n]

    def ngens(self):
        return self.dimension()

    def matrix(self, x=0, coerce=True, copy=True):
        """
        Create a matrix in self.  The entries can be specified either
        as a single list of length nrows*ncols, or as a list of
        lists.

        EXAMPLES:
            sage: M = MatrixSpace(ZZ, 2)
            sage: M.matrix([[1,0],[0,-1]])
            [ 1  0]
            [ 0 -1]
            sage: M.matrix([1,0,0,-1])
            [ 1  0]
            [ 0 -1]
        """
        if isinstance(x, (xrange,xsrange)):
            x = list(x)
        elif isinstance(x, (int, integer.Integer)) and x==1:
            return self.identity_matrix()
        if matrix.is_Matrix(x):
            if x.parent() is self:
                if x.is_immutable():
                    return x
                else:
                    return x.copy()
            x = x.list()
        if isinstance(x, list) and len(x) > 0:
            if isinstance(x[0], list):
                x = sum(x,[])
            elif hasattr(x[0], "is_vector"): # TODO: is this the best way to test that?
                e = []
                for v in x:
                    e = e + v.list()
                copy = False # deep copy?
                x = e
            elif isinstance(x[0], tuple):
                x = list(sum(x,()))
        return self.__matrix_class(self, entries=x, copy=copy, coerce=coerce)

    def matrix_space(self, nrows, ncols, sparse=False):
        if nrows is None:
            nrows = self.__nrows
        if ncols is None:
            ncols = self.__ncols
        return MatrixSpace(self.base_ring(), nrows, ncols,
                        sparse=sparse)

    def ncols(self):
        return self.__ncols

    def nrows(self):
        return self.__nrows

    def random_element(self, X=[-2,-1,1,2], prob=1.0, coerce=True):
        """
        Returns a random element of self.
        """
        prob=float(prob)
        if not isinstance(X, list):
            raise TypeError, "X must be a list"
        if not isinstance(coerce, bool):
            raise TypeError, "coerce must be a bool"
        R = self.base_ring()
        if coerce:
            X = [R(a) for a in X]
        zero = R(0)

        if self.is_sparse():
            nc = self.ncols()
            num_per_row = int(prob * nc) + 1
            z = range(num_per_row)
            v = {}
            for i in xrange(self.nrows()):
                for k in z:
                    v[(i,random.randint(0,nc-1))] = random.choice(X)
        else:
            def f():
                if random.random() < prob:
                    return random.choice(X)
                else:
                    return zero
            v = [f() for _ in xrange(self.nrows()*self.ncols())]
        return self(v, coerce=False, copy=False)

_random = 1

def dict_to_list(entries, nrows, ncols):
    v = [0]*(nrows*ncols)
    for ij, y in entries:
        i,j = ij
        v[i*ncols + j] = y
    return v

def list_to_dict(entries, nrows, ncols):
    d = {}
    if ncols == 0 or nrows == 0:
        return d
    for i in range(len(entries)):
        x = entries[i]
        if x != 0:
            col = i % ncols
            row = i // ncols
            d[(row,col)] = x
    return d


