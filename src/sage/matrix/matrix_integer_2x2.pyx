"""
Two by two matrices over the integers.
"""

include "sage/ext/interrupt.pxi"
include "sage/ext/stdsage.pxi"
from cpython.list cimport *
from cpython.number cimport *
from cpython.ref cimport *

from sage.rings.all import polygen, QQ,ZZ
from sage.rings.integer cimport Integer
from sage.structure.element cimport ModuleElement, Element
from sage.structure.mutability cimport Mutability

cimport matrix_dense
import matrix_dense

cimport matrix

from matrix_space import MatrixSpace
from sage.misc.lazy_attribute import lazy_attribute

class MatrixSpace_ZZ_2x2_class(MatrixSpace):
    """
    Return a space of 2x2 matrices over ZZ, whose elements are of
    type sage.matrix.matrix_integer_2x2.Matrix_integer_2x2 instead of
    sage.matrix.matrix_integer_dense.Matrix_integer_dense.

    NOTE: This class exists only for quickly creating and testing
    elements of this type. Once these become the default for 2x2
    matrices over ZZ, this class should be removed.

    EXAMPLES::

        sage: sage.matrix.matrix_integer_2x2.MatrixSpace_ZZ_2x2()
        Space of 2x2 integer matrices

    By trac ticket #12290, it is a unique parent::

        sage: from sage.matrix.matrix_integer_2x2 import MatrixSpace_ZZ_2x2
        sage: MatrixSpace_ZZ_2x2() is MatrixSpace_ZZ_2x2()
        True
        sage: M = MatrixSpace_ZZ_2x2()
        sage: loads(dumps(M)) is M
        True

    """
    @staticmethod
    def __classcall__(cls):
        """
        EXAMPLES::

            sage: from sage.matrix.matrix_integer_2x2 import MatrixSpace_ZZ_2x2
            sage: MatrixSpace_ZZ_2x2() is MatrixSpace_ZZ_2x2()  # indirect doctest
            True

        """
        return super(MatrixSpace,cls).__classcall__(cls)

    def __init__(self):
        """
        EXAMPLES::

            sage: sage.matrix.matrix_integer_2x2.MatrixSpace_ZZ_2x2()
            Space of 2x2 integer matrices
        """
        MatrixSpace.__init__(self, ZZ, 2, 2, False)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sage.matrix.matrix_integer_2x2.MatrixSpace_ZZ_2x2()
            Space of 2x2 integer matrices
        """
        return "Space of 2x2 integer matrices"

    def _get_matrix_class(self):
        """
        EXAMPLES::

            sage: A = sage.matrix.matrix_integer_2x2.MatrixSpace_ZZ_2x2()
            sage: A._get_matrix_class()
            <type 'sage.matrix.matrix_integer_2x2.Matrix_integer_2x2'>
        """
        return Matrix_integer_2x2

ZZ_2x2_parent = MatrixSpace_ZZ_2x2_class()
def MatrixSpace_ZZ_2x2():
    """
    Return the space of 2x2 integer matrices. (This function
    exists to maintain uniqueness of parents.)

    EXAMPLES::

        sage: M = sage.matrix.matrix_integer_2x2.MatrixSpace_ZZ_2x2()
        sage: M
        Space of 2x2 integer matrices
        sage: M is sage.matrix.matrix_integer_2x2.MatrixSpace_ZZ_2x2()
        True
    """
    return ZZ_2x2_parent


cdef class Matrix_integer_2x2(matrix_dense.Matrix_dense):
    r"""
    The \class{Matrix_integer_2x2} class derives from
    \class{Matrix_dense}, and defines fast functionality for 2x2
    matrices over the integers.  Matrices are represented by four gmp
    integer fields: a, b, c, and d.

    EXAMPLES::

        sage: MS = sage.matrix.matrix_integer_2x2.MatrixSpace_ZZ_2x2()
        sage: m = MS([1,2,3,4]) ; m
        [1 2]
        [3 4]
        sage: TestSuite(m).run()

    We check that #13949 is fixed::

        sage: x = MS([1,2,3,4])
        sage: y = MS([4,5,6,7])
        sage: z = x * y
        sage: z.set_immutable()
    """
    ########################################################################
    # LEVEL 1 functionality
    # x * __cinit__
    # x * __init__
    # x * __dealloc__
    # x * set_unsafe
    # x * get_unsafe
    # x * def _pickle
    # x * def _unpickle
    ########################################################################

    def __cinit__(self, parent, entries,copy, coerce):
        mpz_init(self.a)
        mpz_init(self.b)
        mpz_init(self.c)
        mpz_init(self.d)
        self._entries = &self.a
        self._parent = parent
        self._base_ring = ZZ
        self._nrows = 2
        self._ncols = 2

    def __init__(self, parent, entries, copy, coerce):
        """
        EXAMPLES::

            sage: MS = sage.matrix.matrix_integer_2x2.MatrixSpace_ZZ_2x2()
            sage: MS()
            [0 0]
            [0 0]
            sage: MS(5)
            [5 0]
            [0 5]
            sage: MS([3,4,5,6])
            [3 4]
            [5 6]
            sage: MS([11,3])
            Traceback (most recent call last):
            ...
            TypeError: cannot construct an element of
            Space of 2x2 integer matrices from [11, 3]!
        """
        cdef Py_ssize_t i, n

        if entries is None:
            entries = 0

        if not isinstance(entries, list):
            try:
                x = ZZ(entries)
                is_list = 0
            except TypeError:
                if hasattr(entries, "list"):
                    entries = entries.list()
                    is_list = 1
                else:
                    try:
                        entries = list(entries)
                        is_list = 1
                    except TypeError:
                        raise TypeError("entries must be coercible to a list or the base ring")

        else:
            is_list = 1

        if is_list:

            if len(entries) != self._nrows * self._ncols:
                raise TypeError("entries has the wrong length")

            if coerce:
                mpz_set(self.a, (<Integer>ZZ(entries[0])).value)
                mpz_set(self.b, (<Integer>ZZ(entries[1])).value)
                mpz_set(self.c, (<Integer>ZZ(entries[2])).value)
                mpz_set(self.d, (<Integer>ZZ(entries[3])).value)
            else:
                mpz_set(self.a, (<Integer>entries[0]).value)
                mpz_set(self.b, (<Integer>entries[1]).value)
                mpz_set(self.c, (<Integer>entries[2]).value)
                mpz_set(self.d, (<Integer>entries[3]).value)

        else:

            mpz_set(self.a, (<Integer>x).value)
            mpz_set_si(self.b, 0)
            mpz_set_si(self.c, 0)
            mpz_set(self.d, (<Integer>x).value)

    def  __dealloc__(self):
        mpz_clear(self.a)
        mpz_clear(self.b)
        mpz_clear(self.c)
        mpz_clear(self.d)

    def testit(self):
        print "testing",self._base_ring

    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, value):
        mpz_set(self._entries[(i << 1) | j], (<Integer>value).value)

    cdef get_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        cdef Integer z
        z = Integer.__new__(Integer)
        mpz_set(z.value, self._entries[(i << 1) | j])
        return z

    def _pickle(self):
        """
        EXAMPLES:
            sage: M = sage.matrix.matrix_integer_2x2.MatrixSpace_ZZ_2x2()
            sage: m = M([9,8,7,6])
            sage: m == loads(dumps(m))
            True
        """
        return self.list(), 0

    def _unpickle(self, data, int version):
        """
        EXAMPLES::

            sage: M = sage.matrix.matrix_integer_2x2.MatrixSpace_ZZ_2x2()
            sage: m = M(5)
            sage: m == loads(dumps(m))
            True
        """
        if version == 0:
            mpz_set(self.a, (<Integer>ZZ(data[0])).value)
            mpz_set(self.b, (<Integer>ZZ(data[1])).value)
            mpz_set(self.c, (<Integer>ZZ(data[2])).value)
            mpz_set(self.d, (<Integer>ZZ(data[3])).value)
        else:
            raise RuntimeError("unknown matrix version")

    def __richcmp__(matrix.Matrix self, right, int op):
        """
        EXAMPLES::

            sage: M = sage.matrix.matrix_integer_2x2.MatrixSpace_ZZ_2x2()
            sage: m = M([1,2,3,4])
            sage: n = M([-1,2,3,4])
            sage: m < n
            False
            sage: m > n
            True
            sage: m == m
            True
        """
        return self._richcmp(right, op)

    def __hash__(self):
        """
        Return a hash of self.

        EXAMPLES::

            sage: M = sage.matrix.matrix_integer_2x2.MatrixSpace_ZZ_2x2()
            sage: m = M([1,2,3,4])
            sage: m.set_immutable()
            sage: m.__hash__()
            8
        """
        return self._hash()

    def __iter__(self):
        """
        Return an iterator over the entries of self.

        EXAMPLES::

            sage: M = sage.matrix.matrix_integer_2x2.MatrixSpace_ZZ_2x2()
            sage: m = M([1,2,3,4])
            sage: m.__iter__()
            <listiterator object at ...>
        """
        return iter(self.list())

    cdef Matrix_integer_2x2 _new_c(self):
        cdef Matrix_integer_2x2 x
        x = <Matrix_integer_2x2> Matrix_integer_2x2.__new__(Matrix_integer_2x2, self._parent, None, False, False)
        x._parent = self._parent
        x._nrows = 2
        x._ncols = 2
        return x


    ########################################################################
    # LEVEL 2 functionality
    # x  * cdef _add_
    #    * cdef _mul_
    #    * cdef _cmp_c_impl
    # x  * __neg__
    # x  * __invert__
    # x  * __copy__
    # x  * _multiply_classical
    # x  * _list -- copy of the list of underlying elements
    #    * _dict -- copy of the sparse dictionary of underlying elements
    ########################################################################

    def __copy__(self):
        """
        Return a copy of self.

        EXAMPLES::

            sage: M = sage.matrix.matrix_integer_2x2.MatrixSpace_ZZ_2x2()
            sage: m = M([1,3,4,5])
            sage: n = copy(m)
            sage: n == m
            True
            sage: n is m
            False
            sage: m.is_mutable()
            True
            sage: n.is_mutable()
            True

        The following test against a bug fixed in trac ticket #12290::

            sage: n.base_ring()
            Integer Ring

        """
        cdef Matrix_integer_2x2 x
        x = self._new_c()
        mpz_set(x.a, self.a)
        mpz_set(x.b, self.b)
        mpz_set(x.c ,self.c)
        mpz_set(x.d, self.d)
        x._is_immutable = False
        x._base_ring = self._base_ring
        if self._subdivisions is not None:
            x.subdivide(*self.subdivisions())
        return x

    cpdef ModuleElement _add_(left, ModuleElement right):
        """
        EXAMPLES::

            sage: M = sage.matrix.matrix_integer_2x2.MatrixSpace_ZZ_2x2()
            sage: M([8,7,6,5]) + M([4,5,6,7])
            [12 12]
            [12 12]
        """
        cdef Matrix_integer_2x2 A
        A = left._new_c()
        mpz_add(A.a, left.a, (<Matrix_integer_2x2>right).a)
        mpz_add(A.b, left.b, (<Matrix_integer_2x2>right).b)
        mpz_add(A.c, left.c, (<Matrix_integer_2x2>right).c)
        mpz_add(A.d, left.d, (<Matrix_integer_2x2>right).d)
        return A

    cdef int _cmp_c_impl(left, Element right) except -2:
        """
        EXAMPLES::

            sage: M = sage.matrix.matrix_integer_2x2.MatrixSpace_ZZ_2x2()
            sage: m = M([4,3,7,9])
            sage: n = M([3,3,7,9])
            sage: m < n
            False
            sage: m > n
            True
            sage: m == m
            True
            sage: m == n
            False

        TEST:

        Check that :trac:`14688` is fixed::

            sage: from sage.matrix.matrix_integer_2x2 import MatrixSpace_ZZ_2x2
            sage: M2ZSpace = MatrixSpace_ZZ_2x2()
            sage: A = M2ZSpace([-5, -3, 7, 4])
            sage: B = M2ZSpace([1,0,-2,1])
            sage: A < B
            True
        """
        cdef int c = mpz_cmp(left.a, (<Matrix_integer_2x2>right).a) or \
                     mpz_cmp(left.b, (<Matrix_integer_2x2>right).b) or \
                     mpz_cmp(left.c, (<Matrix_integer_2x2>right).c) or \
                     mpz_cmp(left.d, (<Matrix_integer_2x2>right).d)
        if c < 0: return -1
        if c > 0: return 1
        return 0

    def __neg__(self):
        """
        Return the additive inverse of self.

        EXAMPLES::

            sage: M = sage.matrix.matrix_integer_2x2.MatrixSpace_ZZ_2x2()
            sage: m = M([1,-1,-1,1])
            sage: -m
            [-1  1]
            [ 1 -1]
        """
        cdef Matrix_integer_2x2 A
        A = self._new_c()
        mpz_neg(A.a, self.a)
        mpz_neg(A.b, self.b)
        mpz_neg(A.c, self.c)
        mpz_neg(A.d, self.d)
        return A

    def __invert__(self):
        """
        Return the inverse of self, as a 2x2 matrix over QQ.

        EXAMPLES::

            sage: M = sage.matrix.matrix_integer_2x2.MatrixSpace_ZZ_2x2()
            sage: m = M([2,0,0,2])
            sage: m^-1
            [1/2   0]
            [  0 1/2]
            sage: type(m^-1)
            <type 'sage.matrix.matrix_rational_dense.Matrix_rational_dense'>
        """
        MS = self._matrix_(QQ).parent()
        D = self.determinant()
        return MS([self.get_unsafe(1,1)/D, -self.get_unsafe(0,1)/D, -self.get_unsafe(1,0)/D, self.get_unsafe(0,0)/D], coerce=False, copy=False)

    def __invert__unit(self):
        """
        If self is a unit, i.e. if its determinant is 1 or -1, return
        the inverse of self as a 2x2 matrix over ZZ. If not, raise a
        ZeroDivisionError.

        EXAMPLES::

            sage: M = sage.matrix.matrix_integer_2x2.MatrixSpace_ZZ_2x2()
            sage: m = M([1,8,3,25])
            sage: m.__invert__unit()
            [25 -8]
            [-3  1]
            sage: type(m.__invert__unit())
            <type 'sage.matrix.matrix_integer_2x2.Matrix_integer_2x2'>
            sage: m^-1
            [25 -8]
            [-3  1]
            sage: type(m^-1)
            <type 'sage.matrix.matrix_rational_dense.Matrix_rational_dense'>
            sage: m.__invert__unit() == m^-1
            True
            sage: M([2,0,0,2]).__invert__unit()
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Not a unit!
        """
        cdef Matrix_integer_2x2 A
        cdef Integer D
        D = self.determinant()
        if D.is_one():
            A = self._new_c()
            mpz_set(A.a, self.d)
            mpz_neg(A.b, self.b)
            mpz_neg(A.c, self.c)
            mpz_set(A.d, self.a)
            return A

        elif D.is_unit():
            A = self._new_c()
            mpz_neg(A.a, self.d)
            mpz_set(A.b, self.b)
            mpz_set(A.c, self.c)
            mpz_neg(A.d, self.a)
            return A

        else:
            raise ZeroDivisionError("Not a unit!")
    _invert_unit = __invert__unit

    def _multiply_classical(left, matrix.Matrix _right):
        """
        Multiply the matrices left and right using the classical $O(n^3)$
        algorithm.

        EXAMPLES::

            sage: M = sage.matrix.matrix_integer_2x2.MatrixSpace_ZZ_2x2()
            sage: m = M([9,7,6,4])
            sage: n = M([7,1,3,0])
            sage: m*n
            [84  9]
            [54  6]
            sage: m._multiply_classical(n)
            [84  9]
            [54  6]
        """
        cdef Matrix_integer_2x2 A, right
        cdef mpz_t tmp
        mpz_init(tmp)

        right  = _right
        A = left._new_c()

        mpz_mul(A.a, left.a, right.a)
        mpz_mul(tmp, left.b, right.c)
        mpz_add(A.a, A.a, tmp)

        mpz_mul(A.b, left.a, right.b)
        mpz_mul(tmp, left.b, right.d)
        mpz_add(A.b, A.b, tmp)

        mpz_mul(A.c, left.c, right.a)
        mpz_mul(tmp, left.d, right.c)
        mpz_add(A.c, A.c, tmp)

        mpz_mul(A.d, left.c, right.b)
        mpz_mul(tmp, left.d, right.d)
        mpz_add(A.d, A.d, tmp)

        mpz_clear(tmp)
        return A

    def _list(self):
        """
        EXAMPLES::

            sage: M = sage.matrix.matrix_integer_2x2.MatrixSpace_ZZ_2x2()
            sage: M([8,6,0,8])._list()
            [8, 6, 0, 8]
        """
        return [self.get_unsafe(0,0), self.get_unsafe(0,1),
                self.get_unsafe(1,0), self.get_unsafe(1,1)]

    ########################################################################
    # LEVEL 3 functionality (Optional)
    # x  * cdef _sub_
    # x  * __deepcopy__
    #    * Matrix windows -- only if you need strassen for that base
    #    * Other functions (list them here):
    # x  * determinant
    # x  * charpoly
    ########################################################################

    cpdef ModuleElement _sub_(left, ModuleElement right):
        """
        EXAMPLES::

            sage: M = sage.matrix.matrix_integer_2x2.MatrixSpace_ZZ_2x2()
            sage: M([3,4,5,6]) - M([1,2,3,4])
            [2 2]
            [2 2]
        """
        cdef Matrix_integer_2x2 A
        A = left._new_c()
        mpz_sub(A.a, left.a, (<Matrix_integer_2x2>right).a)
        mpz_sub(A.b, left.b, (<Matrix_integer_2x2>right).b)
        mpz_sub(A.c, left.c, (<Matrix_integer_2x2>right).c)
        mpz_sub(A.d, left.d, (<Matrix_integer_2x2>right).d)
        return A

    def determinant(self):
        """
        Return the determinant of self, which is just ad-bc.

        EXAMPLES::

            sage: M = sage.matrix.matrix_integer_2x2.MatrixSpace_ZZ_2x2()
            sage: M([9,0,8,7]).determinant()
            63
        """
        cdef Integer z
        cdef mpz_t tmp
        mpz_init(tmp)
        z = Integer.__new__(Integer)
        mpz_mul(z.value, self.a, self.d)
        mpz_mul(tmp, self.b, self.c)
        mpz_sub(z.value, z.value, tmp)
        mpz_clear(tmp)
        return z

    def trace(self):
        """
        Return the trace of self, which is just a+d.

        EXAMPLES::

            sage: M = sage.matrix.matrix_integer_2x2.MatrixSpace_ZZ_2x2()
            sage: M([9,0,8,7]).trace()
            16
        """
        cdef Integer z = PY_NEW(Integer)
        mpz_add(z.value, self._entries[0], self._entries[3])
        return z

    def charpoly(self, var='x'):
        """
        Return the charpoly of self as a polynomial in var. Since self
        is 2x2, this is just var^2 - self.trace() * var +
        self.determinant().

        EXAMPLES::

            sage: M = sage.matrix.matrix_integer_2x2.MatrixSpace_ZZ_2x2()
            sage: M([9,0,8,7]).charpoly()
            x^2 - 16*x + 63
            sage: M(range(4)).charpoly()
            x^2 - 3*x - 2
            sage: M(range(4)).charpoly('t')
            t^2 - 3*t - 2
        """
        t = polygen(ZZ,name=var)
        return t*t - t*self.trace() + self.determinant()

    def __deepcopy__(self, *args, **kwds):
        """
        EXAMPLES::

            sage: M = sage.matrix.matrix_integer_2x2.MatrixSpace_ZZ_2x2()
            sage: m = M([5,5,5,5])
            sage: n = deepcopy(m)
            sage: n == m
            True
            sage: n is m
            False
            sage: m[0,0] == n[0,0]
            True
            sage: m[0,0] is n[0,0]
            False
        """
        return self.__copy__()


#######################################################################
