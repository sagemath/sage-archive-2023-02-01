"""
Two by two matrices over the integers.
"""

include "../ext/interrupt.pxi"
include "../ext/stdsage.pxi"
include "../ext/python_list.pxi"
include "../ext/python_number.pxi"
include "../ext/python_ref.pxi"

from sage.rings.integer cimport Integer
from sage.rings.rational_field import QQ
from sage.rings.integer_ring import ZZ
from sage.structure.element cimport ModuleElement, Element

cimport matrix_dense
import matrix_dense

cimport matrix


cdef class Matrix_integer_2x2(matrix_dense.Matrix_dense):
    r"""
    The \class{Matrix_generic_dense} class derives from \class{Matrix}, and
    defines fast functionality 2x2 matrices over the integers.
    Matrices are represented by four gmp integer fields: a, b, c, and d.
    """
    ########################################################################
    # LEVEL 1 functionality
    # x * __new__
    # x * __init__
    # x * __dealloc__
    # x * set_unsafe
    # x * get_unsafe
    # x * def _pickle
    # x * def _unpickle
    ########################################################################

    def __new__(self, parent, entries, coerce, copy):
        mpz_init(self.a)
        mpz_init(self.b)
        mpz_init(self.c)
        mpz_init(self.d)
        self._entries = &self.a

    def __init__(self, parent, entries, copy, coerce):
        matrix.Matrix.__init__(self, parent)

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
                        raise TypeError, "entries must be coercible to a list or the basering"

        else:
            is_list = 1

        if is_list:

            if len(entries) != self._nrows * self._ncols:
                raise TypeError, "entries has the wrong length"

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

    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, value):
        mpz_set(self._entries[(i << 1) | j], (<Integer>value).value)

    cdef get_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        cdef Integer z
        z = Integer.__new__(Integer)
        mpz_set(z.value, self._entries[(i << 1) | j])
        return z

    def _pickle(self):
        return self.list(), 0

    def _unpickle(self, data, int version):
        if version == 0:
            mpz_set(self.a, (<Integer>ZZ(data[0])).value)
            mpz_set(self.b, (<Integer>ZZ(data[1])).value)
            mpz_set(self.c, (<Integer>ZZ(data[2])).value)
            mpz_set(self.d, (<Integer>ZZ(data[3])).value)
        else:
            raise RuntimeError, "unknown matrix version"

    def __richcmp__(matrix.Matrix self, right, int op):  # always need for mysterious reasons.
        return self._richcmp(right, op)
    def __hash__(self):
        return self._hash()

    def __iter__(self):
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
        cdef Matrix_integer_2x2 x
        x = self._new_c()
        mpz_set(x.a, self.a)
        mpz_set(x.b, self.b)
        mpz_set(x.c ,self.c)
        mpz_set(x.d, self.d)
        if self.subdivisions is not None:
            x.subdivide(*self.get_subdivisions())
        return x

    cpdef ModuleElement _add_(left, ModuleElement right):
        cdef Matrix_integer_2x2 A
        A = left._new_c()
        mpz_add(A.a, left.a, (<Matrix_integer_2x2>right).a)
        mpz_add(A.b, left.b, (<Matrix_integer_2x2>right).b)
        mpz_add(A.c, left.c, (<Matrix_integer_2x2>right).c)
        mpz_add(A.d, left.d, (<Matrix_integer_2x2>right).d)
        return A

    cdef int _cmp_c_impl(left, Element right) except -2:
        return mpz_cmp(left.a, (<Matrix_integer_2x2>right).a) or mpz_cmp(left.b, (<Matrix_integer_2x2>right).b) or mpz_cmp(left.c, (<Matrix_integer_2x2>right).c) or mpz_cmp(left.d, (<Matrix_integer_2x2>right).d)

    def __neg__(self):
        cdef Matrix_integer_2x2 A
        A = self._new_c()
        mpz_neg(A.a, self.a)
        mpz_neg(A.b, self.b)
        mpz_neg(A.c, self.c)
        mpz_neg(A.d, self.d)
        return A

    def __invert__(self):
        MS = self._matrix_(QQ).parent()
        D = self.determinant()
        return MS([self.get_unsafe(1,1)/D, -self.get_unsafe(0,1)/D, -self.get_unsafe(1,0)/D, self.get_unsafe(0,0)/D], coerce=False, copy=False)

    def __invert__unit(self):
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
            raise ZeroDivisionError, "Not a unit!"
    _invert_unit = __invert__unit

    def _multiply_classical(left, matrix.Matrix _right):
        """
        Multiply the matrices left and right using the classical $O(n^3)$
        algorithm.
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
        return [self.get_unsafe(0,0), self.get_unsafe(0,1), self.get_unsafe(1,0), self.get_unsafe(1,1)]

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
        cdef Matrix_integer_2x2 A
        A = left._new_c()
        mpz_sub(A.a, left.a, (<Matrix_integer_2x2>right).a)
        mpz_sub(A.b, left.b, (<Matrix_integer_2x2>right).b)
        mpz_sub(A.c, left.c, (<Matrix_integer_2x2>right).c)
        mpz_sub(A.d, left.d, (<Matrix_integer_2x2>right).d)
        return A

    def determinant(self):
        cdef Integer z
        cdef mpz_t tmp
        mpz_init(tmp)
        z = Integer.__new__(Integer)
        mpz_mul(z.value, self.a, self.d)
        mpz_mul(tmp, self.b, self.c)
        mpz_sub(z.value, z.value, tmp)
        mpz_clear(tmp)
        return z

    def charpoly(self, var):
        R = ZZ[var]
        t = R.gen(0)
        cdef Integer z
        z = Integer.__new__(Integer)
        mpz_mul(z.value, self.b, self.c)
        return t*t + t*z + self.determinant()

    def __deepcopy__(self):
        return self.__copy__()


#######################################################################
