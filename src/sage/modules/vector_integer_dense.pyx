"""
Vectors with integer entries

AUTHOR:
    -- William Stein (2007)

EXAMPLES:
    sage: v = vector(ZZ,[1,2,3,4,5])
    sage: v
    (1, 2, 3, 4, 5)
    sage: 3*v
    (3, 6, 9, 12, 15)
    sage: v*7
    (7, 14, 21, 28, 35)
    sage: -v
    (-1, -2, -3, -4, -5)
    sage: v - v
    (0, 0, 0, 0, 0)
    sage: v + v
    (2, 4, 6, 8, 10)
    sage: v * v   # dot product.
    55

We make a large zero vector:
    sage: k = ZZ^100000; k
    Ambient free module of rank 100000 over the principal ideal domain Integer Ring
    sage: v = k(0)
    sage: v[:10]
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

TESTS:
    sage: v = vector(ZZ, [1,2,3,4])
    sage: loads(dumps(v)) == v
    True
"""

###############################################################################
#   SAGE: System for Algebra and Geometry Experimentation
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###############################################################################

include '../ext/interrupt.pxi'
include '../ext/stdsage.pxi'

from sage.structure.element cimport Element, ModuleElement, RingElement, Vector

from sage.rings.integer cimport Integer

cimport free_module_element

from free_module_element import vector

cdef class Vector_integer_dense(free_module_element.FreeModuleElement):
    cdef _new_c(self):
        cdef Vector_integer_dense y
        y = PY_NEW(Vector_integer_dense)
        y._init(self._degree, self._parent)
        return y

    cdef bint is_dense_c(self):
        return 1
    cdef bint is_sparse_c(self):
        return 0

    def __copy__(self):
        cdef Vector_integer_dense y
        y = self._new_c()
        cdef Py_ssize_t i
        for i from 0 <= i < self._degree:
            mpz_init_set(y._entries[i], self._entries[i])
        return y

    cdef _init(self, Py_ssize_t degree, parent):
        self._degree = degree
        self._parent = parent
        self._is_mutable = 1
        self._entries = <mpz_t *> sage_malloc(sizeof(mpz_t) * degree)
        if self._entries == NULL:
            raise MemoryError

    def __new__(self, parent=None, x=None, coerce=True,copy=True):
        self._entries = NULL
        self._is_mutable = 1
        if not parent is None:
            self._init(parent.degree(), parent)

    def __init__(self, parent, x, coerce=True, copy=True):
        cdef Py_ssize_t i
        cdef Integer z
        # we have to do this to avoid a garbage collection error in dealloc
        for i from 0 <= i < self._degree:
            mpz_init(self._entries[i])
        if isinstance(x, (list, tuple)):
            if len(x) != self._degree:
                raise TypeError, "x must be a list of the right length"
            for i from 0 <= i < self._degree:
                z = Integer(x[i])
                mpz_set(self._entries[i], z.value)
            return
        if x != 0:
            raise TypeError, "can't initialize vector from nonzero non-list"

    def __dealloc__(self):
        cdef Py_ssize_t i
        if self._entries:
            for i from 0 <= i < self._degree:
                mpz_clear(self._entries[i])
            sage_free(self._entries)

    cdef int _cmp_c_impl(left, Element right) except -2:
        """
        EXAMPLES:
            sage: v = vector(ZZ, [0,0,0,0])
            sage: v == 0
            True
            sage: v == 1
            False
            sage: v == v
            True
            sage: w = vector(ZZ, [-1,0,0,0])
            sage: w < v
            True
            sage: w > v
            False
        """
        cdef Py_ssize_t i
        cdef int c
        for i from 0 <= i < left.degree():
            c = mpz_cmp(left._entries[i], (<Vector_integer_dense>right)._entries[i])
            if c < 0:
                return -1
            elif c > 0:
                return 1
        return 0

    def __len__(self):
        return self._degree

    def __setitem__(self, Py_ssize_t i, x):
        """
        EXAMPLES:

        """
        if not self._is_mutable:
            raise ValueError, "vector is immutable; please change a copy instead (use self.copy())"
        cdef Integer z
        if i < 0 or i >= self._degree:
            raise IndexError
        else:
            z = Integer(x)
            mpz_set(self._entries[i], z.value)

    def __getitem__(self, Py_ssize_t i):
        """
        Return the ith entry of self.
        """
        cdef Integer z
        z = PY_NEW(Integer)
        if i < 0 or i >= self._degree:
            raise IndexError, 'index out of range'
        else:
            mpz_set(z.value, self._entries[i])
            return z

    def __reduce__(self):
        return (unpickle_v0, (self._parent, self.list(), self._degree))

    cdef ModuleElement _add_c_impl(self, ModuleElement right):
        cdef Vector_integer_dense z, r
        r = right
        z = self._new_c()
        cdef Py_ssize_t i
        for i from 0 <= i < self._degree:
            mpz_init(z._entries[i])
            mpz_add(z._entries[i], self._entries[i], r._entries[i])
        return z


    cdef ModuleElement _sub_c_impl(self, ModuleElement right):
        cdef Vector_integer_dense z, r
        r = right
        z = self._new_c()
        cdef Py_ssize_t i
        for i from 0 <= i < self._degree:
            mpz_init(z._entries[i])
            mpz_sub(z._entries[i], self._entries[i], r._entries[i])
        return z

    cdef Element _dot_product_c_impl(self, Vector right):
        """
        Dot product of dense vectors over the integers.

        EXAMPLES:
            sage: v = vector(ZZ, [1,2,-3]); w = vector(ZZ,[4,3,2])
            sage: v*w
            4
            sage: w*v
            4
        """
        cdef Vector_integer_dense r = right
        cdef Integer z
        z = PY_NEW(Integer)
        cdef mpz_t t
        mpz_init(t)
        mpz_set_si(z.value, 0)
        cdef Py_ssize_t i
        for i from 0 <= i < self._degree:
            mpz_mul(t, self._entries[i], r._entries[i])
            mpz_add(z.value, z.value, t)
        mpz_clear(t)
        return z

    cdef Vector _pairwise_product_c_impl(self, Vector right):
        """
        EXAMPLES:
            sage: v = vector(ZZ, [1,2,-3]); w = vector(ZZ,[4,3,2])
            sage: v.pairwise_product(w)
            (4, 6, -6)
        """
        cdef Vector_integer_dense z, r
        r = right
        z = self._new_c()
        cdef Py_ssize_t i
        for i from 0 <= i < self._degree:
            mpz_init(z._entries[i])
            mpz_mul(z._entries[i], self._entries[i], r._entries[i])
        return z

    cdef ModuleElement _rmul_c_impl(self, RingElement left):
        cdef Vector_integer_dense z
        cdef Integer a
        a = left
        z = self._new_c()
        cdef Py_ssize_t i
        for i from 0 <= i < self._degree:
            mpz_init(z._entries[i])
            mpz_mul(z._entries[i], self._entries[i], a.value)
        return z

    cdef ModuleElement _lmul_c_impl(self, RingElement right):
        cdef Vector_integer_dense z
        cdef Integer a
        a = right
        z = self._new_c()
        cdef Py_ssize_t i
        for i from 0 <= i < self._degree:
            mpz_init(z._entries[i])
            mpz_mul(z._entries[i], self._entries[i], a.value)
        return z

    cdef ModuleElement _neg_c_impl(self):
        cdef Vector_integer_dense z
        z = self._new_c()
        cdef Py_ssize_t i
        for i from 0 <= i < self._degree:
            mpz_init(z._entries[i])
            mpz_neg(z._entries[i], self._entries[i])
        return z

    def n(self, *args, **kwargs):
        """
        Returns a numerical approximation of self by calling the n()
        method on all of its entries.

        EXAMPLES:
            sage: v = vector(ZZ, [1,2,3])
            sage: v.n()
            (1.00000000000000, 2.00000000000000, 3.00000000000000)
            sage: _.parent()
            Vector space of dimension 3 over Real Field with 53 bits of precision
            sage: v.n(prec=75)
            (1.000000000000000000000, 2.000000000000000000000, 3.000000000000000000000)
            sage: _.parent()
            Vector space of dimension 3 over Real Field with 75 bits of precision
        """
        return vector( [e.n(*args, **kwargs) for e in self] )



def unpickle_v0(parent, entries, degree):
    # If you think you want to change this function, don't.
    # Instead make a new version with a name like
    #    make_FreeModuleElement_generic_dense_v1
    # and changed the reduce method below.
    cdef Vector_integer_dense v
    v = PY_NEW(Vector_integer_dense)
    v._init(degree, parent)
    cdef Integer z
    for i from 0 <= i < degree:
        z = Integer(entries[i])
        mpz_init_set(v._entries[i], z.value)
    return v
