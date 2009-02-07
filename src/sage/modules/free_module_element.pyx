r"""
Elements of free modules

AUTHOR:
    -- William Stein
    -- Josh Kantor

TODO:
    Change to use a get_unsafe / set_unsafe, etc., structure exactly
    like with matrices, since we'll have to define a bunch of special
    purpose implementations of vectors easily and systematically.

EXAMPLES:
We create a vector space over $\QQ$ and a subspace of this space.
    sage: V = QQ^5
    sage: W = V.span([V.1, V.2])

Arithmetic operations always return something in the ambient space,
since there is a canonical map from $W$ to $V$ but not from $V$ to $W$.

    sage: parent(W.0 + V.1)
    Vector space of dimension 5 over Rational Field
    sage: parent(V.1 + W.0)
    Vector space of dimension 5 over Rational Field
    sage: W.0 + V.1
    (0, 2, 0, 0, 0)
    sage: W.0 - V.0
    (-1, 1, 0, 0, 0)

Next we define modules over $\ZZ$ and a finite field.
    sage: K = ZZ^5
    sage: M = GF(7)^5

Arithmetic between the $\QQ$  and $\ZZ$ modules is defined, and
the result is always over $\QQ$, since there is a canonical coercion
map to $\QQ$.
    sage: K.0 + V.1
    (1, 1, 0, 0, 0)
    sage: parent(K.0 + V.1)
    Vector space of dimension 5 over Rational Field

Since there is no canonical coercion map to the finite field from $\QQ$
the following arithmetic is not defined:
    sage: V.0 + M.0
    Traceback (most recent call last):
    ...
    TypeError: unsupported operand parent(s) for '+': 'Vector space of dimension 5 over Rational Field' and 'Vector space of dimension 5 over Finite Field of size 7'

However, there is a map from $\ZZ$ to the finite field, so the
following is defined, and the result is in the finite field.
    sage: w = K.0 + M.0; w
    (2, 0, 0, 0, 0)
    sage: parent(w)
    Vector space of dimension 5 over Finite Field of size 7
    sage: parent(M.0 + K.0)
    Vector space of dimension 5 over Finite Field of size 7

Matrix vector multiply:
    sage: MS = MatrixSpace(QQ,3)
    sage: A = MS([0,1,0,1,0,0,0,0,1])
    sage: V = QQ^3
    sage: v = V([1,2,3])
    sage: v * A
    (2, 1, 3)

TESTS:
    sage: D = 46341
    sage: u = 7
    sage: R = Integers(D)
    sage: p = matrix(R,[[84, 97, 55, 58, 51]])
    sage: 2*p.row(0)
    (168, 194, 110, 116, 102)
"""

import math
import operator

include '../ext/cdefs.pxi'
include '../ext/stdsage.pxi'
import sage.misc.misc as misc
import sage.misc.latex

from sage.structure.sequence import Sequence

from sage.structure.element cimport Element, ModuleElement, RingElement, Vector as element_Vector
from sage.structure.element import canonical_coercion

import sage.rings.arith

from sage.rings.ring import is_Ring
from sage.rings.infinity import Infinity
import sage.rings.integer_ring
import sage.rings.integer
from sage.rings.real_double import RDF
from sage.rings.complex_double import CDF
from sage.misc.derivative import multi_derivative

#from sage.matrix.matrix cimport Matrix

def is_FreeModuleElement(x):
    return isinstance(x, FreeModuleElement)

def vector(arg0, arg1=None, arg2=None, sparse=None):
    r"""
    Return a vector over R with given entries.

    CALL FORMATS:
        1. vector(object)
        2. vector(ring, object)
        3. vector(object, ring)
        4. vector(numpy_array)

    In each case, give sparse=[True|False] as an option.

    INPUT:
        elts -- entries of a vector (either a list or dict).
        R -- ring
        sparse -- optional

    OUTPUT:
        An element of the free module over R of rank len(elts).

    EXAMPLES:
        sage: v = vector([1,2,3]); v
        (1, 2, 3)
        sage: v.parent()
        Ambient free module of rank 3 over the principal ideal domain Integer Ring
        sage: v = vector([1,2,3/5]); v
        (1, 2, 3/5)
        sage: v.parent()
        Vector space of dimension 3 over Rational Field

    All entries must \emph{canonically} coerce to some common ring:
        sage: v = vector([17, GF(11)(5), 19/3]); v
        Traceback (most recent call last):
        ...
        TypeError: unable to find a common ring for all elements

        sage: v = vector([17, GF(11)(5), 19]); v
        (6, 5, 8)
        sage: v.parent()
        Vector space of dimension 3 over Finite Field of size 11
        sage: v = vector([17, GF(11)(5), 19], QQ); v
        (17, 5, 19)
        sage: v.parent()
        Vector space of dimension 3 over Rational Field
        sage: v = vector((1,2,3), QQ); v
        (1, 2, 3)
        sage: v.parent()
        Vector space of dimension 3 over Rational Field
        sage: v = vector(QQ, (1,2,3)); v
        (1, 2, 3)
        sage: v.parent()
        Vector space of dimension 3 over Rational Field
        sage: v = vector(vector([1,2,3])); v
        (1, 2, 3)
        sage: v.parent()
        Ambient free module of rank 3 over the principal ideal domain Integer Ring

    You can also use \code{free_module_element}, which is the same as \code{vector}.
        sage: free_module_element([1/3, -4/5])
        (1/3, -4/5)

    Make a vector mod 3 out of a vector over ZZ:
        sage: vector(vector([1,2,3]), GF(3))
        (1, 2, 0)

    Here we illustrate the creation of sparse vectors by using a dictionary:
        sage: vector({1:1.1, 3:3.14})
        (0.000000000000000, 1.10000000000000, 0.000000000000000, 3.14000000000000)

    Any 1 dimensional numpy array of type float or complex may be passed to vector. The result
    will vector in the appropriate dimensional vector space over the real double field or the
    complex double field. The data in the array must be contiguous so columnwise slices of numpy matrices
    will rase an exception.

        sage: import numpy
        sage: x=numpy.random.randn(10)
        sage: y=vector(x)
        sage: v=numpy.random.randn(10)*numpy.complex(0,1)
        sage: w=vector(v)

    If any of the arguments to vector have Python type int, long,
    real, or complex, they will first be coerced to the appropriate
    Sage objects.  This fixes trac \#3847:
        sage: v = vector([int(0)]); v
        (0)
        sage: v[0].parent()
        Integer Ring
        sage: v = vector(range(10)); v
        (0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
        sage: v[3].parent()
        Integer Ring
        sage: v = vector([float(23.4), int(2), complex(2+7*I), long(1)]); v
        (23.4, 2.0, 2.0 + 7.0*I, 1.0)
        sage: v[1].parent()
        Complex Double Field
    """
    if hasattr(arg0, '_vector_'):
        if arg1 is None:
            try:
                arg1 = sage.rings.integer_ring.ZZ
                return arg0._vector_(arg1)
            except TypeError:
                return arg0._vector_()
        return arg0._vector_(arg1)

    if hasattr(arg1, '_vector_'):
        return arg1._vector_(arg0)

    if sage.rings.integer.is_Integer(arg1) or isinstance(arg1,(int,long)):
        if arg2 is None:
            arg1 = [0]*arg1
        else:
            if len(arg2) != arg1:
                raise ValueError, "incompatible degrees in vector constructor"
            arg1 = arg2

    if is_Ring(arg0):
        R = arg0
        v = arg1
    elif is_Ring(arg1):
        R = arg1
        v = arg0
    else:
        v = arg0
        R = None
    if isinstance(v, dict):
        if sparse is None:
            sparse = True
        v, R = prepare_dict(v, R)

    from numpy import ndarray
    from free_module import VectorSpace
    if isinstance(v,ndarray):
        if len(v.shape)==1:
            if str(v.dtype).count('float')==1:
                V=VectorSpace(RDF,v.shape[0])
                import vector_real_double_dense
                _v=vector_real_double_dense.Vector_real_double_dense(V, v)
                return _v
            if str(v.dtype).count('complex')==1:
                V=VectorSpace(CDF,v.shape[0])
                import vector_complex_double_dense
                _v=vector_complex_double_dense.Vector_complex_double_dense(V, v)
                return _v

    else:
        if sparse is None:
            sparse = False
        v, R = prepare(v, R)
    if sparse:
        import free_module  # slow -- can we improve
        return free_module.FreeModule(R, len(v), sparse=True)(v)
    else:
        return (R**len(v))(v)

free_module_element = vector

def prepare(v, R):
    v = Sequence(v, universe=R, use_sage_types=True)
    ring = v.universe()
    if not is_Ring(ring):
        raise TypeError, "unable to find a common ring for all elements"
    return v, ring

def prepare_dict(w, R):
    """
    EXAMPLES:
        sage: from sage.modules.free_module_element import prepare_dict
        sage: prepare_dict({3:1 , 5:3}, QQ)
        ([0, 0, 0, 1, 0, 3], Rational Field)
        sage: prepare_dict({},QQ)
        ([], Rational Field)
    """
    Z = w.items()
    cdef Py_ssize_t n
    n = max([-1]+[key for key,value in Z])+1
    X = [0]*n
    for key, value in Z:
        X[key] = value
    return prepare(X, R)

cdef class FreeModuleElement(element_Vector):   # abstract base class
    """
    An element of a generic free module.
    """
    def __init__(self, parent):
        self._parent = parent
        self._degree = parent.degree()
        self._is_mutable = 1

    def _magma_init_(self, magma):
        r"""
        Convert self to Magma.

        EXAMPLES:
            sage: F = FreeModule(ZZ, 2, inner_product_matrix=matrix(ZZ, 2, 2, [1, 0, 0, -1]))
            sage: v = F([1, 2])
            sage: M = magma(v); M # optional - magma
            (1 2)
            sage: M.Type() # optional - magma
            ModTupRngElt
            sage: M.Parent() # optional - magma
            Full RSpace of degree 2 over Integer Ring
            Inner Product Matrix:
            [ 1  0]
            [ 0 -1]
            sage: M.sage() # optional - magma
            (1, 2)
            sage: M.sage() == v # optional - magma
            True
            sage: M.sage().parent() is v.parent() # optional - magma
            True

            sage: v = vector(QQ, [1, 2, 5/6])
            sage: M = magma(v); M # optional - magma
            (  1   2 5/6)
            sage: M.Type() # optional - magma
            ModTupFldElt
            sage: M.Parent() # optional - magma
            Full Vector space of degree 3 over Rational Field
            sage: M.sage() # optional - magma
            (1, 2, 5/6)
            sage: M.sage() == v # optional - magma
            True
            sage: M.sage().parent() is v.parent() # optional - magma
            True
        """
        # Get a reference to Magma version of parent.
        R = magma(self.parent())
        # Get list of coefficients.
        v = ','.join([a._magma_init_(magma) for a in self.list()])
        return '%s![%s]' % (R.name(), v)

    def __hash__(self):
        if self._is_mutable:
            raise TypeError, "mutable vectors are unhasheable"
        return hash(tuple(self))

    def _vector_(self, R=None):
        r"""Return self as a vector.

        EXAMPLES:
            sage: v = vector(ZZ, [2, 12, 22])
            sage: vector(v)
            (2, 12, 22)
            sage: vector(GF(7), v)
            (2, 5, 1)
            sage: vector(v, ZZ['x', 'y'])
            (2, 12, 22)

	    sage: vector(vector((1, 6.8)))
	    (1.00000000000000, 6.80000000000000)
	    sage: vector(vector(SR, (1, sqrt(2)) ) )
	    (1, sqrt(2))
        """
        if R is None:
            R = self.base_ring()
        return self.change_ring(R)

    def _matrix_(self, R=None):
        r"""Return self as a row matrix.

        EXAMPLES:
            sage: v = vector(ZZ, [2, 12, 22])
            sage: vector(v)
            (2, 12, 22)
            sage: vector(GF(7), v)
            (2, 5, 1)
            sage: vector(v, ZZ['x', 'y'])
            (2, 12, 22)
        """
        if R is None:
            R = self.base_ring()
        from sage.matrix.constructor import matrix
        return matrix(R, [list(self)])

    def transpose(self):
        r"""Return self as a column matrix.

        EXAMPLES:
            sage: v = vector(ZZ, [2, 12, 22])
            sage: transpose(vector(v))
            [ 2]
            [12]
            [22]

            sage: transpose(vector(GF(7), v))
            [2]
            [5]
            [1]

            sage: transpose(vector(v, ZZ['x', 'y']))
            [ 2]
            [12]
            [22]
        """
        return self._matrix_().transpose()

    def _hash(self):
        return hash(tuple(list(self)))

    def copy(self):
        """
        Make a copy of this vector.

        EXAMPLES:
            sage: v = vector([1..5]); v
            (1, 2, 3, 4, 5)
            sage: w = v.copy()
            sage: v == w
            True
            sage: v is w
            False

            sage: v = vector([1..5], sparse=True); v
            (1, 2, 3, 4, 5)
            sage: v.copy()
            (1, 2, 3, 4, 5)
        """
        return self.__copy__()

    def __copy__(self):
        if self.is_sparse():
            return self.parent()(self.dict())
        else:
            return self.parent()(self.list())

    def set_immutable(self):
        """
        Make this vector immutable.  This operation can't be undone.

        EXAMPLES:
            sage: v = vector([1..5]); v
            (1, 2, 3, 4, 5)
            sage: v[1] = 10
            sage: v.set_immutable()
            sage: v[1] = 10
            Traceback (most recent call last):
            ...
            ValueError: vector is immutable; please change a copy instead (use self.copy())
        """
        self._is_mutable = 0

    def is_mutable(self):
        """
        Return True if this vector is mutable, i.e., the entries can be changed.

        EXAMPLES:
            sage: v = vector(QQ['x,y'], [1..5]); v.is_mutable()
            True
            sage: v.set_immutable()
            sage: v.is_mutable()
            False
        """
        return self._is_mutable

    def is_immutable(self):
        """
        Return True if this vector is immutable, i.e., the entries cannot be changed.

        EXAMPLES:
            sage: v = vector(QQ['x,y'], [1..5]); v.is_immutable()
            False
            sage: v.set_immutable()
            sage: v.is_immutable()
            True
        """
        return not self._is_mutable

    def change_ring(self, R):
        """
        Change the base ring of this vector, by coercing each element
        of this vector into R.

        EXAMPLES:
            sage: v = vector(QQ['x,y'], [1..5]); v.change_ring(GF(3))
            (1, 2, 0, 1, 2)
        """
        P = self.parent()
        if P.base_ring() is R:
            return self
        return P.change_ring(R)(self)

    def additive_order(self):
        """
        Return the additive order of self.

        EXAMPLES:
            sage: v = vector(Integers(4), [1,2])
            sage: v.additive_order()
            4

            sage: v = vector([1,2,3])
            sage: v.additive_order()
            +Infinity

            sage: v = vector(Integers(30), [6, 15]); v
            (6, 15)
            sage: v.additive_order()
            10
            sage: 10*v
            (0, 0)
        """
        v = [None]*self.degree()
        cdef int i
        for i from 0 <= i < self.degree():
            v[i] = self[i].additive_order()
            if v[i] == +Infinity:
               return +Infinity
        return sage.rings.arith.LCM(v)

    def iteritems(self):
        return self.dict(copy=False).iteritems()

    def __abs__(self):
        """
        Return the square root of the sum of the squares of the entries of this vector.

        EXAMPLES:
            sage: v = vector([1..5]); abs(v)
            sqrt(55)
            sage: v = vector(RDF, [1..5]); abs(v)
            7.4161984871
        """
        return sum([x**2 for x in self.list()]).sqrt()

    def norm(self, p=sage.rings.integer.Integer(2)):
        """
        Return the p-norm of this vector, where p can be a real
        number >= 1, Infinity, or a symbolic expression.
        If p=2 (default), this is the usual Euclidean norm;
        if p=Infinity, this is the maximum norm; if p=1, this is
        the taxicab (Manhattan) norm.

        EXAMPLES:
            sage: v = vector([1,2,-3])
            sage: v.norm(5)
            276^(1/5)

        The default is the usual Euclidean norm:
            sage: v.norm()
            sqrt(14)
            sage: v.norm(2)
            sqrt(14)

        The infinity norm is the maximum size of any entry:
            sage: v.norm(Infinity)
            3

        Any real or symbolic value works:
            sage: v=vector(RDF,[1,2,3])
            sage: v.norm(5)
            3.07738488539
            sage: v.norm(pi/2)
            4.2165958647
            sage: var('a b c d p')
            (a, b, c, d, p)
            sage: v=vector([a, b, c, d])
            sage: v.norm(p)
            (abs(d)^p + abs(c)^p + abs(b)^p + abs(a)^p)^(1/p)
        """
        abs_self = [abs(x) for x in self]
        if p == Infinity:
            return max(abs_self)
        try:
            pr = RDF(p)
            if pr < 1:
                raise ValueError, "%f is not greater than or equal to 1" %(pr)
        except TypeError:
            pass

        s = sum([a**p for a in abs_self])
        return s**(1/p)

    cdef int _cmp_c_impl(left, Element right) except -2:
        """
        EXAMPLES:
            sage: v = vector(QQ, [0,0,0,0])
            sage: v == 0
            True
            sage: v == 1
            False
            sage: v == v
            True
            sage: w = vector(QQ, [-1,0,0,0])
            sage: w < v
            True
            sage: w > v
            False
        """
        cdef Py_ssize_t i
        cdef int c
        for i from 0 <= i < left.degree():
            c = cmp(left[i], right[i])
            if c: return c
        return 0

    def __nonzero__(self):
        """
        EXAMPLES:
            sage: V = vector(ZZ, [0, 0, 0, 0])
            sage: bool(V)
            False
            sage: V = vector(ZZ, [1, 2, 3, 5])
            sage: bool(V)
            True
        """
        return self != 0

    def __getitem__(self, i):
        raise NotImplementedError

    def __invert__(self):
        raise NotImplementedError

    def __len__(self):
        return self.parent().degree()

    def __mod__(self, p):
        """
        EXAMPLES:
            sage: V = vector(ZZ, [5, 9, 13, 15])
            sage: V % 7
            (5, 2, 6, 1)
            sage: parent(V % 7)
            Ambient free module of rank 4 over the principal ideal domain Integer Ring
        """
        return self.parent()([x % p for x in self.list()], copy=False, coerce=False, check=False)

    def Mod(self, p):
        """
        EXAMPLES:
            sage: V = vector(ZZ, [5, 9, 13, 15])
            sage: V.Mod(7)
            (5, 2, 6, 1)
            sage: parent(V.Mod(7))
            Vector space of dimension 4 over Ring of integers modulo 7
        """
        return self.change_ring(self.base_ring().quotient_ring(p))

    def list(self, copy=True):
        d = self.degree()
        v = [0]*d
        for i in range(d):
            v[i] = self[i]
        return v

    def list_from_positions(self, positions):
        cdef Py_ssize_t j
        v = [0]*len(positions)
        j = 0
        for i in positions:
            v[j] = self[i]
            j = j + 1
        return v

    def lift(self):
        """
        EXAMPLES:
            sage: V = vector(Integers(7), [5, 9, 13, 15]) ; V
            (5, 2, 6, 1)
            sage: V.lift()
            (5, 2, 6, 1)
            sage: parent(V.lift())
            Ambient free module of rank 4 over the principal ideal domain Integer Ring
        """
        return self.change_ring(self.base_ring().cover_ring())

    def __pos__(self):
        return self

    def __pow__(self, n, dummy):
        raise NotImplementedError

    def _repr_(self):
        """
        String representation of a vector.

        EXAMPLES:
            sage: vector(QQ, [])._repr_()
            '()'
            sage: vector(QQ, range(5))._repr_()
            '(0, 1, 2, 3, 4)'

        Symbolic are not displayed using ASCII art.
            sage: x = var('x')
            sage: v = vector([x/(2*x)+sqrt(2)+var('theta')^3,x/(2*x)]); v
            (theta^3 + sqrt(2) + 1/2, 1/2)
            sage: v._repr_()
            '(theta^3 + sqrt(2) + 1/2, 1/2)'
        """
        d = self.degree()
        if d == 0: return "()"
        # compute column widths
        S = [repr(x) for x in self.list(copy=False)]
        #width = max([len(x) for x in S])
        s = "("
        for i in xrange(d):
            if i == d-1:
                sep = ""
            else:
                sep=", "
            entry = S[i]
            #if i > 0:
            #    entry = " "*(width-len(entry)) + entry
            s = s + entry + sep
        s = s + ")"
        return s

    def _maple_init_(self):
        """
        EXAMPLES:
            sage: v = vector(ZZ, 4, range(4))               #optional
            sage: maple(v)                                  #optional
            Vector[row](4, [0,1,2,3])

            sage: v = vector(QQ, 3, [2/3, 0, 5/4])          #optional
            sage: maple(v)                                  #optional
            Vector[row](3, [2/3,0,5/4])

            sage: P.<x> = ZZ[]                                       #optional
            sage: v = vector(P, 3, [x^2 + 2, 2*x + 1, -2*x^2 + 4*x]) #optional
            sage: maple(v)                                           #optional
            Vector[row](3, [x^2+2,2*x+1,-2*x^2+4*x])
        """
        return "Vector[row](%s)"%(str(self.list()))

    def __setitem__(self, i, x):
        raise NotImplementedError

    def __getslice__(self, Py_ssize_t i, Py_ssize_t j):
        """
        EXAMPLES:
            sage: v = vector(QQ['x,y'], [1,2, 'x*y', 'x^2-y^2']); v
            (1, 2, x*y, x^2 - y^2)
            sage: v[1:]
            (2, x*y, x^2 - y^2)
            sage: v[:2]
            (1, 2)
            sage: type(v[1:])
            <type 'sage.modules.free_module_element.FreeModuleElement_generic_dense'>
            sage: v = vector(CDF,[1,2,(3,4)]); v
            (1.0, 2.0, 3.0 + 4.0*I)
            sage: w = v[1:]; w
            (2.0, 3.0 + 4.0*I)
            sage: parent(w)
            Vector space of dimension 2 over Complex Double Field
        """
        return vector(self.base_ring(), self.list()[i:j])

    def __setslice__(self, i, j, value):
        """
        EXAMPLES:
            sage: v = vector(CDF,[1,2,(3,4)]); v
            (1.0, 2.0, 3.0 + 4.0*I)
            sage: v[1:] = (1,3); v
            (1.0, 1.0, 3.0)

            sage: v.set_immutable()
            sage: v[1:2] = [3,5]
            Traceback (most recent call last):
            ...
            ValueError: vector is immutable; please change a copy instead (use self.copy())
        """
        if not self._is_mutable:
            raise ValueError, "vector is immutable; please change a copy instead (use self.copy())"
        cdef Py_ssize_t k, d, n
        d = self.degree()
        R = self.base_ring()
        n = 0
        for k from i <= k < j:
            if k >= d:
                return
            if k >= 0:
                self[k] = R(value[n])
                n = n + 1



    def __richcmp__(left, right, int op):
        cdef int ld, rd
        if not isinstance(left, FreeModuleElement) or not isinstance(right, FreeModuleElement):
            # use the generic compare
            return (<Element>left)._richcmp(right, op)
        ld = (<FreeModuleElement>left)._degree
        rd = (<FreeModuleElement>right)._degree
        if ld < rd:
            return (<Element>left)._rich_to_bool(op, -1)
        elif ld > rd:
            return (<Element>left)._rich_to_bool(op, 1)
        if (<FreeModuleElement>left)._parent.base_ring() is (<FreeModuleElement>right)._parent.base_ring():
            return (<Element>left)._rich_to_bool(op, (
                    <FreeModuleElement>left)._cmp_same_ambient_c(right))
        return (<Element>left)._richcmp(right, op)


    cdef int _cmp_same_ambient_c(left, FreeModuleElement right):
        return cmp(left.list(copy=False), right.list(copy=False))

    def degree(self):
        return self._degree

    def denominator(self):
        R = self.base_ring()
        if self.degree() == 0: return 1
        x = self.list()
        d = x[0].denominator()
        for y in x:
            d = d.lcm(y.denominator())
        return d

    def dict(self, copy=True):
        e = {}
        for i in xrange(self.degree()):
            c = self[i]
            if c != 0:
                e[i] = c
        return e

    #############################
    # Plotting
    #############################
    def plot(self, plot_type=None, **kwds):
        """
        INPUT:

            plot_type -- (default: 'arrow' if v has 3 or fewer components,
            otherwise 'step') type of plot.  Options are 'arrow' to draw an
            arrow; 'point' to draw a point at the coordinates
            specified by the vector; 'step' to draw a step function
            representing the coordinates of the vector.  Both 'arrow'
            and 'point' raise exceptions if the vector has more than 3
            dimensions.

        EXAMPLES:
            sage: v = vector(RDF, (1,2))
            sage: eps = 0.1
            sage: plot(v, plot_type='arrow')
            sage: plot(v, plot_type='point')
            sage: plot(v, plot_type='step') # calls v.plot_step()
            sage: plot(v, plot_type='step', eps=eps, xmax=5, hue=0)
            sage: v = vector(RDF, (1,2,1))
            sage: plot(v) # defaults to an arrow plot
            sage: plot(v, plot_type='arrow')
            sage: from sage.plot.plot3d.shapes2 import frame3d
            sage: plot(v, plot_type='point')+frame3d((0,0,0), v.list())
            sage: plot(v, plot_type='step') # calls v.plot_step()
            sage: plot(v, plot_type='step', eps=eps, xmax=5, hue=0)
            sage: v = vector(RDF, (1,2,3,4))
            sage: plot(v) # defaults to a step plot


        """
        # Give sensible defaults based on the vector length
        if plot_type is None:
            if len(self)<=3:
                plot_type='arrow'
            else:
                plot_type='step'

        if plot_type == 'arrow' or plot_type == 'point':
            dimension = len(self)
            if dimension == 3:
                from sage.plot.plot3d.shapes import arrow3d, Sphere
                # Sphere complains if the radius is given twice,
                # so we have to delete it from kwds if it is given.
                radius = kwds.pop('radius', .02)

                if plot_type == 'arrow':
                    return arrow3d((0,0,0), self, radius=radius, **kwds)
                else:
                    return Sphere(radius, **kwds).translate(self.list())
            elif dimension < 3:
                vectorlist = self.list()
                if dimension < 2:
                    # pad to make 2-dimensional
                    vectorlist.extend([0]*(2-dimension))

                from sage.plot.all import arrow, point
                if plot_type == 'arrow':
                    return arrow((0,0), vectorlist, **kwds)
                else:
                    return point(vectorlist, **kwds)
            else:
                raise ValueError, "arrow and point plots require vectors with 3 or fewer components"

        elif plot_type == 'step':
            return self.plot_step(**kwds)
        else:
            raise NotImplementedError, "plot_type was unrecognized"




    def plot_step(self, xmin=0, xmax=1, eps=None, res=None,
             connect=True, **kwds):
        """
        INPUT:
            xmin -- (default: 0) start x position to start plotting
            xmax -- (default: 1) stop x position to stop plotting
            eps -- (default: determined by xmax) we view this vector
                   as defining a function at the points xmin, xmin +
                   eps, xmin + 2*eps, ...,
            res -- (default: all points) total number of points to include
                   in the graph
            connect -- (default: True) if True draws a line; otherwise draw
                       a list of points.

        EXAMPLES:
            sage: eps=0.1
            sage: v = vector(RDF, [sin(n*eps) for n in range(100)])
            sage: v.plot_step(eps=eps, xmax=5, hue=0)
        """
        if res is None:
            res = self.degree()
        if eps is None:
            eps = float(xmax - xmin)/res
        v = []
        x = xmin
        for i in range(0, self.degree(), int(math.ceil(self.degree()/res))):
            y = float(self[i])
            if x > xmax:
                break
            v.append((x,y))
            x += eps
            v.append((x,y))
        from sage.plot.all import line, points
        if connect:
            return line(v, **kwds)
        else:
            return points(v, **kwds)

    def dot_product(self, right):
        """
        Return the dot product of self and right, which is the sum
        of the product of the corresponding entries.

        INPUT:
            right -- vector of the same degree as self.  it need not
                     be in the same vector space as self, as long as
                     the coefficients can be multiplied.

        EXAMPLES:
            sage: V = FreeModule(ZZ, 3)
            sage: v = V([1,2,3])
            sage: w = V([4,5,6])
            sage: v.dot_product(w)
            32

            sage: W = VectorSpace(GF(3),3)
            sage: w = W([0,1,2])
            sage: w.dot_product(v)
            2
            sage: w.dot_product(v).parent()
            Finite Field of size 3

        Implicit coercion is well defined (irregardless of order), so
        we get 2 even if we do the dot product in the other order.

            sage: v.dot_product(w)
            2
        """
        if not PY_TYPE_CHECK(right, FreeModuleElement):
            raise TypeError, "right must be a free module element"
        r = right.list(copy=False)
        l = self.list(copy=False)
        if len(r) != len(l):
            raise ArithmeticError, "degrees (%s and %s) must be the same"%(len(l),len(r))
        if len(r) == 0:
            return self._parent.base_ring()(0)
        sum = l[0] * r[0]
        cdef Py_ssize_t i
        for i from 1 <= i < len(l):
            sum += l[i] * r[i]
        return sum

    def cross_product(self, right):
        """
        Return the cross product of self and right, which is only
        defined for vectors of length 3.

        This product is performed under the assumption that the basis
        vectors are orthonormal.

        EXAMPLES:
            sage: v = vector([1,2,3]); w = vector([0,5,-9])
            sage: v.cross_product(v)
            (0, 0, 0)
            sage: u = v.cross_product(w); u
            (-33, 9, 5)
            sage: u.dot_product(v)
            0
            sage: u.dot_product(w)
            0
        """
        if not PY_TYPE_CHECK(right, FreeModuleElement):
            raise TypeError, "right must be a free module element"
        r = right.list(copy=False)
        l = self.list(copy=False)
        if len(r) != 3 or len(l) != 3:
            raise ArithmeticError, "Cross product only defined for vectors of length three, not (%s and %s)"%(len(l),len(r))
        return vector([l[1]*r[2] - l[2]*r[1],
                       l[2]*r[0] - l[0]*r[2],
                       l[0]*r[1] - l[1]*r[0]])


    def pairwise_product(self, right):
        """
        Return the dot product of self and right, which is a vector of
        of the product of the corresponding entries.

        INPUT:
            right -- vector of the same degree as self.  it need not
                     be in the same vector space as self, as long as
                     the coefficients can be multiplied.

        EXAMPLES:
            sage: V = FreeModule(ZZ, 3)
            sage: v = V([1,2,3])
            sage: w = V([4,5,6])
            sage: v.pairwise_product(w)
            (4, 10, 18)
            sage: sum(v.pairwise_product(w)) == v.dot_product(w)
            True

            sage: W = VectorSpace(GF(3),3)
            sage: w = W([0,1,2])
            sage: w.pairwise_product(v)
            (0, 2, 0)
            sage: w.pairwise_product(v).parent()
            Vector space of dimension 3 over Finite Field of size 3

        Implicit coercion is well defined (irregardless of order), so
        we get 2 even if we do the dot product in the other order.

            sage: v.pairwise_product(w).parent()
            Vector space of dimension 3 over Finite Field of size 3

        TESTS:
            sage: x, y = var('x, y')

            sage: parent(vector(ZZ,[1,2]).pairwise_product(vector(ZZ,[1,2])))
            Ambient free module of rank 2 over the principal ideal domain Integer Ring
            sage: parent(vector(ZZ,[1,2]).pairwise_product(vector(QQ,[1,2])))
            Vector space of dimension 2 over Rational Field
            sage: parent(vector(QQ,[1,2]).pairwise_product(vector(ZZ,[1,2])))
            Vector space of dimension 2 over Rational Field
            sage: parent(vector(QQ,[1,2]).pairwise_product(vector(QQ,[1,2])))
            Vector space of dimension 2 over Rational Field

            sage: parent(vector(QQ,[1,2,3,4]).pairwise_product(vector(ZZ[x],[1,2,3,4])))
            Ambient free module of rank 4 over the principal ideal domain Univariate Polynomial Ring in x over Rational Field
            sage: parent(vector(ZZ[x],[1,2,3,4]).pairwise_product(vector(QQ,[1,2,3,4])))
            Ambient free module of rank 4 over the principal ideal domain Univariate Polynomial Ring in x over Rational Field

            sage: parent(vector(QQ,[1,2,3,4]).pairwise_product(vector(ZZ[x][y],[1,2,3,4])))
            Ambient free module of rank 4 over the integral domain Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
            sage: parent(vector(ZZ[x][y],[1,2,3,4]).pairwise_product(vector(QQ,[1,2,3,4])))
            Ambient free module of rank 4 over the integral domain Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field

            sage: parent(vector(QQ[x],[1,2,3,4]).pairwise_product(vector(ZZ[x][y],[1,2,3,4])))
            Ambient free module of rank 4 over the integral domain Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
            sage: parent(vector(ZZ[x][y],[1,2,3,4]).pairwise_product(vector(QQ[x],[1,2,3,4])))
            Ambient free module of rank 4 over the integral domain Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field

            sage: parent(vector(QQ[y],[1,2,3,4]).pairwise_product(vector(ZZ[x][y],[1,2,3,4])))
            Ambient free module of rank 4 over the integral domain Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field
            sage: parent(vector(ZZ[x][y],[1,2,3,4]).pairwise_product(vector(QQ[y],[1,2,3,4])))
            Ambient free module of rank 4 over the integral domain Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field

            sage: parent(vector(ZZ[x],[1,2,3,4]).pairwise_product(vector(ZZ[y],[1,2,3,4])))
            Traceback (most recent call last):
            ...
            TypeError: no common canonical parent for objects with parents: 'Ambient free module of rank 4 over the integral domain Univariate Polynomial Ring in x over Integer Ring' and 'Ambient free module of rank 4 over the integral domain Univariate Polynomial Ring in y over Integer Ring'
            sage: parent(vector(ZZ[x],[1,2,3,4]).pairwise_product(vector(QQ[y],[1,2,3,4])))
            Traceback (most recent call last):
            ...
            TypeError: no common canonical parent for objects with parents: 'Ambient free module of rank 4 over the integral domain Univariate Polynomial Ring in x over Integer Ring' and 'Ambient free module of rank 4 over the principal ideal domain Univariate Polynomial Ring in y over Rational Field'
            sage: parent(vector(QQ[x],[1,2,3,4]).pairwise_product(vector(ZZ[y],[1,2,3,4])))
            Traceback (most recent call last):
            ...
            TypeError: no common canonical parent for objects with parents: 'Ambient free module of rank 4 over the principal ideal domain Univariate Polynomial Ring in x over Rational Field' and 'Ambient free module of rank 4 over the integral domain Univariate Polynomial Ring in y over Integer Ring'
            sage: parent(vector(QQ[x],[1,2,3,4]).pairwise_product(vector(QQ[y],[1,2,3,4])))
            Traceback (most recent call last):
            ...
            TypeError: no common canonical parent for objects with parents: 'Ambient free module of rank 4 over the principal ideal domain Univariate Polynomial Ring in x over Rational Field' and 'Ambient free module of rank 4 over the principal ideal domain Univariate Polynomial Ring in y over Rational Field'
        """
        if not PY_TYPE_CHECK(right, FreeModuleElement):
            raise TypeError, "right must be a free module element"
        if self._parent is not (<FreeModuleElement>right)._parent:
            self, right = canonical_coercion(self, right)
        return self._pairwise_product_(right)

    def element(self):
        return self

    def get(self, i):
        """
        get is meant to be more efficient
        than getitem, because it does not do
        any error checking.
        """
        return self[i]

    def set(self, i, x):
        """
        set is meant to be more efficient
        than setitem, because it does not do
        any error checking or coercion.  Use with care.
        """
        self[i] = x


    def normalize(self):
        """
        Return this vector divided through by the first nonzero entry of this vector.

        EXAMPLES:
            sage: v = vector(QQ,[0,4/3,5,1,2])
            sage: v.normalize()
            (0, 1, 15/4, 3/4, 3/2)
        """
        cdef Py_ssize_t i
        for i from 0 <= i < self._degree:
            if self[i] != 0:
                return (~self[i]) * self
        return self

    def inner_product(self, right):
        """
        Returns the inner product of self and other, with respect to
        the inner product defined on the parent of self.

        EXAMPLES:

            sage: I = matrix(ZZ,3,[2,0,-1,0,2,0,-1,0,6])
            sage: M = FreeModule(ZZ, 3, inner_product_matrix = I)
            sage: (M.0).inner_product(M.0)
            2
            sage: K = M.span_of_basis([[0/2,-1/2,-1/2], [0,1/2,-1/2],[2,0,0]])
            sage: (K.0).inner_product(K.0)
            2

        """
        if self.parent().is_ambient() and self.parent()._inner_product_is_dot_product():
            return self.dot_product(right)
        if not isinstance(right, FreeModuleElement):
            raise TypeError, "right must be a free module element"
        M = self.parent()
        if M.is_ambient() or M.uses_ambient_inner_product():
            A = M.ambient_module().inner_product_matrix()
            return A.linear_combination_of_rows(self).dot_product(right)
        else:
            A = M.inner_product_matrix()
            v = M.coordinate_vector(self)
            w = M.coordinate_vector(right)
            return A.linear_combination_of_rows(v).dot_product(w)

    def is_dense(self):
        return self.is_dense_c()

    cdef bint is_dense_c(self):
        return self.parent().is_dense()

    def is_sparse(self):
        return self.is_sparse_c()

    cdef bint is_sparse_c(self):
        return self.parent().is_sparse()

    def is_vector(self):
        return True

    def _mathematica_init_(self):
        """
        Returns string representation of this vector as a Mathematica list.

        EXAMPLES:
            sage: vector((1,2,3), QQ)._mathematica_init_()
            '{1/1, 2/1, 3/1}'
            sage: mathematica(vector((1,2,3), QQ))  #optional -- requires mathematica
            {1, 2, 3}
            sage: a = vector(SR, 5, [1, x, x^2, sin(x), pi]); a
            (1, x, x^2, sin(x), pi)
            sage: a._mathematica_init_()
            '{1, x, (x) ^ (2), Sin[x], Pi}'
        """
        return '{' + ', '.join([x._mathematica_init_() for x in self.list()]) + '}'

##     def zero_out_positions(self, P):
##         """
##         Set the positions of self in the list P equal to 0.
##         """
##         z = self.base_ring()(0)
##         d = self.degree()
##         for n in P:
##             self[n] = z

    def nonzero_positions(self):
        """
        Return the sorted list of integers i such that self[i] != 0.
        """
        z = self.base_ring()(0)
        v = self.list()
        cdef Py_ssize_t i
        return [i for i from 0 <= i < self.degree() if v[i] != z]

    def support(self):   # do not override.
        """
        Return the integers i such that self[i] != 0.
        This is the same as the \code{nonzero_positions} function.
        """
        return self.nonzero_positions()

    def _latex_(self):
        """
        Return a latex representation of self.  For example, if self is
        the free module element (1,2,3,4), then following latex is
        generated: "(1,2,3,4)"  (without the quotes).
        """
        s = '\\left('
        for a in self.list():
            s = s + sage.misc.latex.latex(a) + ','
        if len(self.list()) > 0:
            s = s[:-1]  # get rid of last comma
        return s + '\\right)'

    def dense_vector(self):
        if self.is_dense():
            return self
        else:
            return self.parent().ambient_module().dense_module()(self.list())

    def sparse_vector(self):
        if self.is_sparse():
            return self
        else:
            return self.parent().ambient_module().sparse_module()(self.list())


    def apply_map(self, phi, R=None, sparse=None):
        """
        Apply the given map phi (an arbitrary Python function or
        callable object) to this free module element.  If R is not
        given, automatically determine the base ring of the resulting
        element.

        INPUT:
            phi -- arbitrary Python function or callable object
            R -- (optional) ring
            sparse -- True or False will control whether the result
              is sparse.  By default, the result is sparse iff self
              is sparse.

        OUTPUT:
            a free module element over R

        EXAMPLES:
            sage: m = vector([1,x,sin(x+1)])
            sage: m.apply_map(x^2)
            (1, x^2, sin(x + 1)^2)
            sage: m.apply_map(sin)
            (sin(1), sin(x), sin(sin(x + 1)))

            sage: m = vector(ZZ, 9, range(9))
            sage: k.<a> = GF(9)
            sage: m.apply_map(k)
            (0, 1, 2, 0, 1, 2, 0, 1, 2)

        In this example, we explicitly specify the codomain.

            sage: s = GF(3)
            sage: f = lambda x: s(x)
            sage: n = m.apply_map(f, k); n
            (0, 1, 2, 0, 1, 2, 0, 1, 2)
            sage: n.parent()
            Vector space of dimension 9 over Finite Field in a of size 3^2

        If your map sends 0 to a non-zero value, then your resulting
        vector is not mathematically sparse:

            sage: v = vector([0] * 6 + [1], sparse=True); v
            (0, 0, 0, 0, 0, 0, 1)
            sage: v2 = v.apply_map(lambda x: x+1); v2
            (1, 1, 1, 1, 1, 1, 2)

        but it's still represented with a sparse data type:

            sage: parent(v2)
            Ambient sparse free module of rank 7 over the principal ideal domain Integer Ring

        This data type is inefficient for dense vectors, so you may
        want to specify sparse=False:

            sage: v2 = v.apply_map(lambda x: x+1, sparse=False); v2
            (1, 1, 1, 1, 1, 1, 2)
            sage: parent(v2)
            Ambient free module of rank 7 over the principal ideal domain Integer Ring

        Or if you have a map that will result in mostly zeroes, you may
        want to specify sparse=True:

            sage: v = vector(srange(10))
            sage: v2 = v.apply_map(lambda x: 0 if x else 1, sparse=True); v2
            (1, 0, 0, 0, 0, 0, 0, 0, 0, 0)
            sage: parent(v2)
            Ambient sparse free module of rank 10 over the principal ideal domain Integer Ring

        TESTS:
            sage: m = vector(SR,[])
            sage: m.apply_map(lambda x: x*x) == m
            True

        Check that we don't unnecessarily apply phi to 0 in the sparse case:
            sage: m = vector(ZZ, range(1, 4), sparse=True)
            sage: m.apply_map(lambda x: 1/x)
            (1, 1/2, 1/3)

            sage: parent(vector(RDF, (), sparse=True).apply_map(lambda x: x, sparse=True))
            Sparse vector space of dimension 0 over Real Double Field
            sage: parent(vector(RDF, (), sparse=True).apply_map(lambda x: x, sparse=False))
            Vector space of dimension 0 over Real Double Field
            sage: parent(vector(RDF, (), sparse=False).apply_map(lambda x: x, sparse=True))
            Sparse vector space of dimension 0 over Real Double Field
            sage: parent(vector(RDF, (), sparse=False).apply_map(lambda x: x, sparse=False))
            Vector space of dimension 0 over Real Double Field
        """
        if sparse is None:
            sparse = self.is_sparse()

        if self._degree == 0:
            if sparse == self.is_sparse():
                return self.copy()
            elif sparse:
                return self.sparse_vector()
            else:
                return self.dense_vector()

        v = None

        if self.is_sparse():
            if len(self.dict(copy=False)) < self._degree:
                # OK, we have some zero entries.
                zero_res = phi(self.base_ring()(0))
                if not zero_res.is_zero():
                    # And phi maps 0 to a non-zero value.
                    v = [zero_res] * self._degree
                    for i,z in self.dict(copy=False).items():
                        v[i] = phi(z)

            if v is None:
                # phi maps 0 to 0 (or else we don't have any zeroes at all)
                v = dict([(i,phi(z)) for i,z in self.dict(copy=False).items()])
        else:
            v = [phi(z) for z in self.list()]

        if R is None:
            return vector(v, sparse=sparse)
        else:
            return vector(R, v, sparse=sparse)


    def _derivative(self, var=None):
        """
        Differentiate with respect to var by differentiating each
        element with respect to var.

        SEE ALSO:
            self.derivative()

        EXAMPLES:
            sage: v = vector([1,x,x^2])
            sage: v._derivative(x)
            (0, 1, 2*x)
            sage: type(v._derivative()) == type(v)
            True
            sage: v = vector([1,x,x^2], sparse=True)
            sage: v._derivative(x)
            (0, 1, 2*x)
            sage: type(v._derivative()) == type(v)
            True
        """
        # We would just use apply_map, except that Cython doesn't
        # allow lambda functions
        if self._degree == 0:
            return self.copy()

        if self.is_sparse():
            v = dict([(i,z.derivative(var)) for i,z in self.dict().items()])
        else:
            v = [z.derivative(var) for z in self.list()]

        return self.parent().ambient_module()(v)

    def derivative(self, *args):
        """
        Derivative with respect to variables supplied in args.

        Multiple variables and iteration counts may be supplied; see
        documentation for the global derivative() function for more details.

        EXAMPLES:
            sage: v = vector([1,x,x^2])
            sage: v.derivative(x)
            (0, 1, 2*x)
            sage: type(v.derivative()) == type(v)
            True
            sage: v = vector([1,x,x^2], sparse=True)
            sage: v.derivative(x)
            (0, 1, 2*x)
            sage: type(v.derivative()) == type(v)
            True
            sage: v.derivative(x,x)
            (0, 0, 2)
        """
        return multi_derivative(self, args)

#############################################
# Generic dense element
#############################################
def make_FreeModuleElement_generic_dense(parent, entries, degree):
    # If you think you want to change this function, don't.
    # Instead make a new version with a name like
    #    make_FreeModuleElement_generic_dense_v1
    # and changed the reduce method below.
    cdef FreeModuleElement_generic_dense v
    v = FreeModuleElement_generic_dense.__new__(FreeModuleElement_generic_dense)
    v._entries = entries
    v._parent = parent
    v._degree = degree
    return v

def make_FreeModuleElement_generic_dense_v1(parent, entries, degree, is_mutable):
    # If you think you want to change this function, don't.
    # Instead make a new version with a name like
    #    make_FreeModuleElement_generic_dense_v2
    # and changed the reduce method below.
    cdef FreeModuleElement_generic_dense v
    v = FreeModuleElement_generic_dense.__new__(FreeModuleElement_generic_dense)
    v._entries = entries
    v._parent = parent
    v._degree = degree
    v._is_mutable = is_mutable
    return v

cdef class FreeModuleElement_generic_dense(FreeModuleElement):
    """
    A generic dense element of a free module.
    """
    ## these work fine on the command line but fail in doctests :-(
##         TESTS:
##             sage: V = ZZ^3
##             sage: loads(dumps(V)) == V
##             True
##             sage: v = V.0
##             sage: loads(dumps(v)) == v
##             True
##             sage: v = (QQ['x']^3).0
##             sage: loads(dumps(v)) == v
##             True
    cdef _new_c(self, object v):
        # Create a new dense free module element with minimal overhead and
        # no type checking.
        cdef FreeModuleElement_generic_dense x
        x = PY_NEW(FreeModuleElement_generic_dense)
        x._is_mutable = 1
        x._parent = self._parent
        x._entries = v
        x._degree = self._degree
        return x

    cdef bint is_dense_c(self):
        return 1

    cdef bint is_sparse_c(self):
        return 0

    def _hash(self):
        return hash(tuple(list(self)))

    def __copy__(self):
        return self._new_c(list(self._entries))

    def __init__(self, parent, entries, coerce=True, copy=True):
        FreeModuleElement.__init__(self, parent)
        R = self.parent().base_ring()
        if entries == 0:
            entries = [R(0)]*self.degree()
        else:
            if not isinstance(entries, (list, tuple)):

                raise TypeError, "entries (=%s) must be a list"%(entries, )

            if len(entries) != self.degree():
                raise TypeError, "entries must be a list of length %s"%\
                            self.degree()
            if coerce:
                try:
                    entries = [R(x) for x in entries]
                except TypeError:
                    raise TypeError, "Unable to coerce entries (=%s) to %s"%(entries, R)
            elif copy:
                # Make a copy
                entries = list(entries)
        self._entries = entries

    cpdef ModuleElement _add_(left, ModuleElement right):
        """
        Add left and right.
        """
        cdef Py_ssize_t i, n
        n = PyList_Size(left._entries)
        v = [None]*n
        for i from 0 <= i < n:
            v[i] = (<RingElement>left._entries[i])._add_(<RingElement>
                                            ((<FreeModuleElement_generic_dense>right)._entries[i]))
        return left._new_c(v)

    cpdef ModuleElement _sub_(left, ModuleElement right):
        """
        Subtract right from left.

        EXAMPLES:
            sage: V = QQ^5
            sage: W = V.span([V.1, V.2])
            sage: W.0 - V.0
            (-1, 1, 0, 0, 0)
            sage: V.0 - W.0
            (1, -1, 0, 0, 0)
        """
        cdef Py_ssize_t i, n
        n = PyList_Size(left._entries)
        v = [None]*n
        for i from 0 <= i < n:
            v[i] = (<RingElement>left._entries[i])._sub_(<RingElement>
                                            ((<FreeModuleElement_generic_dense>right)._entries[i]))
        return left._new_c(v)

    cpdef ModuleElement _rmul_(self, RingElement left):
        """
        EXAMPLES:
            sage: V = ZZ['x']^5
            sage: 5 * V.0
            (5, 0, 0, 0, 0)
        """
        if left._parent is self._parent._base:
            v = [left._mul_(<RingElement>x) for x in self._entries]
        else:
            v = [left * x for x in self._entries]
        return self._new_c(v)

    cpdef ModuleElement _lmul_(self, RingElement right):
        if right._parent is self._parent._base:
            v = [(<RingElement>x)._mul_(right) for x in self._entries]
        else:
            v = [x * right for x in self._entries]
        return self._new_c(v)

    cpdef Element _dot_product_(left, element_Vector right):
        """
        Return the dot product of left and right.

        EXAMPLES:
            sage: R.<x> = QQ[]
            sage: v = vector([x,x^2,3*x]); w = vector([2*x,x,3+x])
            sage: v*w
            x^3 + 5*x^2 + 9*x
            sage: (x*2*x) + (x^2*x) + (3*x*(3+x))
            x^3 + 5*x^2 + 9*x
            sage: w*v
            x^3 + 5*x^2 + 9*x
        """
        return left.dot_product(right)

    cpdef element_Vector _pairwise_product_(left, element_Vector right):
        """
        EXAMPLES:
            sage: R.<x> = QQ[]
            sage: v = vector([x,x^2,3*x]); w = vector([2*x,x,3+x])
            sage: v.pairwise_product(w)
            (2*x^2, x^3, 3*x^2 + 9*x)
            sage: w.pairwise_product(v)
            (2*x^2, x^3, 3*x^2 + 9*x)
        """
        if not right.parent() == left.parent():
            right = left.parent().ambient_module()(right)
        # Component wise vector * vector multiplication.
        cdef Py_ssize_t i, n
        n = PyList_Size(left._entries)
        v = [None]*n
        for i from 0 <= i < n:
            v[i] = (<RingElement>left._entries[i])._mul_((<FreeModuleElement_generic_dense>right)._entries[i])
        return left._new_c(v)

    def __reduce__(self):
        return (make_FreeModuleElement_generic_dense_v1, (self._parent, self._entries, self._degree, self._is_mutable))

    def __getitem__(self, Py_ssize_t i):
        """
        EXAMPLES:
            sage: v = vector([RR(1), RR(2)]); v
            (1.00000000000000, 2.00000000000000)
            sage: v[0]
            1.00000000000000
            sage: v[-1]
            2.00000000000000
            sage: v[4]
            Traceback (most recent call last):
            ...
            IndexError: index must be between -2 and 1
            sage: v[-4]
            Traceback (most recent call last):
            ...
            IndexError: index must be between -2 and 1

        """
        if isinstance(i, slice):
            return list(self)[i]
        degree = self.degree()
        i = int(i)
        if i < 0:
            i += degree
        if i < 0 or i >= self.degree():
            raise IndexError, "index must be between -%s and %s"%(degree, degree-1)
        return self._entries[i]

    def __setitem__(self, Py_ssize_t i, value):
        """
        Set entry i of self to value.
        """
        if not self._is_mutable:
            raise ValueError, "vector is immutable; please change a copy instead (use self.copy())"
        i = int(i)
        #if not isinstance(i, int):
        #    raise TypeError, "index must an integer"
        if i < 0 or i >= self.degree():
            raise IndexError, "index (i=%s) must be between 0 and %s"%(i,
                            self.degree()-1)
        self._entries[i] = self.base_ring()(value)

    def __setslice__(self, Py_ssize_t i, Py_ssize_t j, value):
        """
        EXAMPLES:
            sage: v = vector(QQ['x,y'], [1,2, 'x*y'])
            sage: v
            (1, 2, x*y)
            sage: v[1:]
            (2, x*y)
            sage: v[1:] = [4,5]; v
            (1, 4, 5)
            sage: v[:2] = [5,(6,2)]; v
            (5, 3, 5)
            sage: v[:2]
            (5, 3)
        """
        if not self._is_mutable:
            raise ValueError, "vector is immutable; please change a copy instead (use self.copy())"
        cdef Py_ssize_t k, n, d
        d = self.degree()
        R = self.base_ring()
        n = 0
        for k from i <= k < j:
            if k >= d:
                return
            if k >= 0:
                self._entries[k] = R(value[n])
                n = n + 1

    def list(self, copy=True):
        if copy:
            return list(self._entries)
        else:
            return self._entries

    cdef int _cmp_c_impl(left, Element right) except -2:
        """
        Compare two free module elements with identical parents.

        Free module elements are compared in lexicographic order on
        the underlying list of coefficients.  A dense a sparse free
        module element are equal if their coefficients are the same.
        """
        return cmp(left._entries, (<FreeModuleElement_generic_dense>right)._entries)

    def __call__(self, *args, **kwargs):
        """
        Calling a free module element returns the result of calling each component.

        EXAMPLES:
            sage: x, y = var('x,y')
            sage: f = x^2 + y^2
            sage: g = f.gradient()
            sage: g
            (2*x, 2*y)
            sage: type(g)
            <type 'sage.modules.free_module_element.FreeModuleElement_generic_dense'>
            sage: g(y=2, x=3)
            (6, 4)
            sage: f(x,y) = x^2 + y^2
            sage: g = f.gradient()
            sage: g(3,2)
            (6, 4)
            sage: g(x=3, y=2)
            (6, 4)
        """
        return vector([e(*args, **kwargs) for e in self])

    def n(self, *args, **kwargs):
        """
        Returns a numerical approximation of self by calling the n()
        method on all of its entries.

        EXAMPLES:
            sage: v = vector(RealField(212), [1,2,3])
            sage: v.n()
            (1.00000000000000, 2.00000000000000, 3.00000000000000)
            sage: _.parent()
            Vector space of dimension 3 over Real Field with 53 bits of precision
            sage: v.n(prec=75)
            (1.000000000000000000000, 2.000000000000000000000, 3.000000000000000000000)
            sage: _.parent()
            Vector space of dimension 3 over Real Field with 75 bits of precision
        """
        return vector([e.n(*args, **kwargs) for e in self])

#############################################
# Generic sparse element
#############################################
def _sparse_dot_product(v, w):
    """
    v and w are dictionaries with integer keys.
    """
    x = set(v.keys()).intersection(set(w.keys()))
    return sum([v[k]*w[k] for k in x])

def make_FreeModuleElement_generic_sparse(parent, entries, degree):
    cdef FreeModuleElement_generic_sparse v
    v = FreeModuleElement_generic_sparse.__new__(FreeModuleElement_generic_sparse)
    v._entries = entries
    v._parent = parent
    v._degree = degree
    return v

def make_FreeModuleElement_generic_sparse_v1(parent, entries, degree, is_mutable):
    cdef FreeModuleElement_generic_sparse v
    v = FreeModuleElement_generic_sparse.__new__(FreeModuleElement_generic_sparse)
    v._entries = entries
    v._parent = parent
    v._degree = degree
    v._is_mutable = is_mutable
    return v

cdef class FreeModuleElement_generic_sparse(FreeModuleElement):
    """
    A generic sparse free module element is a dictionary with keys
    ints i and entries in the base ring.

    EXAMPLES:

    Pickling works:
        sage: v = FreeModule(ZZ, 3, sparse=True).0
        sage: loads(dumps(v)) == v
        True
        sage: v = FreeModule(Integers(8)['x,y'], 5, sparse=True).1
        sage: loads(dumps(v)) - v
        (0, 0, 0, 0, 0)

        sage: a = vector([-1,0,1/1],sparse=True); b = vector([-1/1,0,0],sparse=True)
        sage: a.parent()
        Sparse vector space of dimension 3 over Rational Field
        sage: b - a
        (0, 0, -1)
    """
    cdef _new_c(self, object v):
        # Create a new sparse free module element with minimal overhead and
        # no type checking.
        cdef FreeModuleElement_generic_sparse x
        x = PY_NEW(FreeModuleElement_generic_sparse)
        x._is_mutable = 1
        x._parent = self._parent
        x._entries = v
        x._degree = self._degree
        return x

    cdef bint is_dense_c(self):
        return 0

    cdef bint is_sparse_c(self):
        return 1

    def __copy__(self):
        return self._new_c(dict(self._entries))

    def __init__(self, parent,
                 entries=0,
                 coerce=True,
                 copy=True):
        #WARNING: In creation, we do not check that the i pairs satisfy
        #     0 <= i < degree.
        FreeModuleElement.__init__(self, parent)
        R = self.base_ring()
        if entries == 0:
            entries = {}
        else:
            if isinstance(entries, list):
                if len(entries) != self.degree():
                    raise TypeError, "entries has the wrong length"
                x = entries
                entries = {}
                for i in xrange(self.degree()):
                    if x[i] != 0:
                        entries[i] = x[i]
                copy = False
            if not isinstance(entries, dict):
                raise TypeError, "entries must be a dict"
            if coerce:
                try:
                    for k, x in entries.iteritems():
                        entries[k] = R(x)
                except TypeError:
                    raise TypeError, "Unable to coerce values of entries dict (=%s) to %s"%(entries, R)
            elif copy:
                # Make a copy
                entries = dict(entries)
        self._entries = entries

    cpdef ModuleElement _add_(left, ModuleElement right):
        """
        Add left and right.
        """
        cdef object v, e
        e = dict((<FreeModuleElement_generic_sparse>right)._entries)
        for i, a in left._entries.iteritems():
            if e.has_key(i):
                e[i] = (<RingElement>a)._add_(<RingElement> e[i])
            else:
                e[i] = a
        return left._new_c(e)

    cpdef ModuleElement _sub_(left, ModuleElement right):
        cdef object v, e
        e = dict(left._entries)   # dict to make a copy
        for i, a in (<FreeModuleElement_generic_sparse>right)._entries.iteritems():
            if e.has_key(i):
                e[i] = (<RingElement> e[i])._sub_(<RingElement>a)
            else:
                e[i] = -a
        return left._new_c(e)


    cpdef ModuleElement _lmul_(self, RingElement right):
        cdef object v
        v = PyDict_New()
        for i, a in self._entries.iteritems():
            v[i] = (<RingElement>a)._mul_(right)
        return self._new_c(v)

    cpdef ModuleElement _rmul_(self, RingElement left):
        cdef object v
        v = PyDict_New()
        for i, a in self._entries.iteritems():
            v[i] = left._mul_(a)
        return self._new_c(v)

    cpdef Element _dot_product_(left, element_Vector right):
        """
        Return the dot product of left and right.

        EXAMPLES:
            sage: v = vector([1,2,0], sparse=True); w = vector([0,5,-9], sparse=True)
            sage: v * w
            10
            sage: w * v
            10
        """
        cdef object v, e, z
        e = dict((<FreeModuleElement_generic_sparse>right)._entries)
        z = left.base_ring()(0)
        for i, a in left._entries.iteritems():
            if e.has_key(i):
                z += (<RingElement>a)._mul_(<RingElement> e[i])
        return z

    cpdef element_Vector _pairwise_product_(left, element_Vector right):
        # Component wise vector * vector multiplication.
        cdef object v, e
        e = dict((<FreeModuleElement_generic_sparse>right)._entries)
        for i, a in left._entries.iteritems():
            if e.has_key(i):
                e[i] = (<RingElement>a)._mul_(<RingElement> e[i])
            else:
                e[i] = a
        return left._new_c(e)

    cdef int _cmp_c_impl(left, Element right) except -2:
        """
        Compare two sparse free module elements.

        Free module elements are compared in lexicographic order on
        the underlying list of coefficients.  A dense a sparse free
        module element are equal if their coefficients are the same.
        """
        a = left._entries.items()
        a.sort()
        b = (<FreeModuleElement_generic_dense>right)._entries.items()
        b.sort()
        return cmp(a,b)

    def iteritems(self):
        return self._entries.iteritems()

    def __reduce__(self):
        return (make_FreeModuleElement_generic_sparse_v1, (self._parent, self._entries, self._degree, self._is_mutable))

    def __getitem__(self, i):
        """
        EXAMPLES:
            sage: v = vector([RR(1), RR(2)], sparse=True); v
            (1.00000000000000, 2.00000000000000)
            sage: v[0]
            1.00000000000000
            sage: v[-1]
            2.00000000000000
            sage: v[5]
            Traceback (most recent call last):
            ...
            IndexError: index must be between -2 and 1
            sage: v[-3]
            Traceback (most recent call last):
            ...
            IndexError: index must be between -2 and 1
        """
        i = int(i)
        degree = self.degree()
        if i < 0:
            i += degree
        if i < 0 or i >= degree:
            raise IndexError, "index must be between %s and %s"%(-degree,
                            degree-1)
        if self._entries.has_key(i):
            return self._entries[i]
        return self.base_ring()(0)  # optimize this somehow

    def get(self, i):
        """
        Like __getitem__ but with no type or bounds checking.
        Returns 0 if access is out of bounds.
        """
        i = int(i)
        if self._entries.has_key(i):
            return self._entries[i]
        return self.base_ring()(0)  # optimize this somehow


    def set(self, i, x):
        """
        Like __setitem__ but with no type or bounds checking.
        """
        if not self._is_mutable:
            raise ValueError, "vector is immutable; please change a copy instead (use self.copy())"
        i = int(i)
        if x == 0:
            if self._entries.has_key(i):
                del self._entries[i]
            return
        self._entries[i] = x

    def __setitem__(self, i, value):
        """
        Set the ith entry of self to value.

        EXAMPLES:
            sage: V = VectorSpace(GF(17), 10000000, sparse=True)
            sage: w = V(0)
            sage: w[39893] = 20
            sage: w[39893]
            3
            sage: parent(w[39893])
            Finite Field of size 17
            sage: w[39893] = sqrt(2)
            Traceback (most recent call last):
            ...
            TypeError: unable to convert x (=sqrt(2)) to an integer
        """
        if not self._is_mutable:
            raise ValueError, "vector is immutable; please change a copy instead (use self.copy())"
        i = int(i)
        #if not isinstance(i, int):
        #    raise TypeError, "index must an integer"
        if i < 0 or i >= self.degree():
            raise IndexError, "index (i=%s) must be between 0 and %s"%(i,
                            self.degree()-1)
        self.set(i, self._parent.base_ring()(value))

    def denominator(self):
        R = self.base_ring()
        x = self.entries()
        if len(x) == 0:
            return 1
        Z = x.iteritems()
        d = Z.next()[1].denominator()
        for _, y in Z:
            d = d.lcm(y.denominator())
        return d

    def dict(self, copy=True):
        if copy:
            return dict(self._entries)
        else:
            return self._entries

    def list(self, copy=True):
        cdef Py_ssize_t n
        n = self._parent.degree()
        z = self._parent.base_ring()(0)
        v = [z]*n
        for i, a in self._entries.iteritems():
            v[i] = a
        return v

    def nonzero_positions(self):
        """
        Returns the set of pairs (i,j) such that self[i,j] != 0.
        """
        K = self._entries.keys()
        K.sort()
        return K


    def n(self, *args, **kwargs):
        """
        Returns a numerical approximation of self by calling the n()
        method on all of its entries.

        EXAMPLES:
            sage: v = vector(RealField(200), [1,2,3], sparse=True)
            sage: v.n()
            (1.00000000000000, 2.00000000000000, 3.00000000000000)
            sage: _.parent()
            Sparse vector space of dimension 3 over Real Field with 53 bits of precision
            sage: v.n(prec=75)
            (1.000000000000000000000, 2.000000000000000000000, 3.000000000000000000000)
            sage: _.parent()
            Sparse vector space of dimension 3 over Real Field with 75 bits of precision
        """
        return vector(dict([(e[0],e[1].n(*args, **kwargs)) for e in self._entries.iteritems()]), sparse=True)

