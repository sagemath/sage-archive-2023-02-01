"""
Vectors with integer mod n entries, with n small.

AUTHOR:
    -- William Stein (2007)

EXAMPLES:
    sage: v = vector(Integers(8),[1,2,3,4,5])
    sage: v
    (1, 2, 3, 4, 5)
    sage: 3*v
    (3, 6, 1, 4, 7)
    sage: v*7
    (7, 6, 5, 4, 3)
    sage: -v
    (7, 6, 5, 4, 3)
    sage: v - v
    (0, 0, 0, 0, 0)
    sage: v + v
    (2, 4, 6, 0, 2)
    sage: v * v
    (1, 4, 1, 0, 1)

We make a large zero vector:
    sage: k = Integers(8)^100000; k
    Ambient free module of rank 100000 over Ring of integers modulo 8
    sage: v = k(0)
    sage: v[:10]
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

We multiply a vector by a matrix:
    sage: a = (GF(97)^5)(range(5))
    sage: m = matrix(GF(97),5,range(25))
    sage: a*m
    (53, 63, 73, 83, 93)

"""

###############################################################################
#   SAGE: System for Algebra and Geometry Experimentation
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###############################################################################

include '../ext/interrupt.pxi'
include '../ext/stdsage.pxi'

from sage.rings.integer_mod cimport IntegerMod_int, IntegerMod_abstract

from sage.structure.element cimport Element, ModuleElement, RingElement, Vector

cimport free_module_element

cdef class Vector_modn_dense(free_module_element.FreeModuleElement):
    cdef _new_c(self):
        cdef Vector_modn_dense y
        y = PY_NEW(Vector_modn_dense)
        y._init(self._degree, self._parent, self._p)
        return y

    def __copy__(self):
        cdef Vector_modn_dense y
        y = self._new_c()
        cdef Py_ssize_t i
        for i from 0 <= i < self._degree:
            y._entries[i] = self._entries[i]
        return y

    cdef _init(self, Py_ssize_t degree, parent, mod_int p):
        self._degree = degree
        self._parent = parent
        self._p = p
        self._entries = <mod_int *> sage_malloc(sizeof(mod_int) * degree)
        if self._entries == NULL:
            raise MemoryError

    def __new__(self, parent=None, x=None, coerce=True, copy=True):
        self._entries = NULL
        if not parent is None:
            self._init(parent.degree(), parent, parent.base_ring().order())

    def __init__(self, parent, x, coerce=True, copy=True):
        cdef Py_ssize_t i
        cdef mod_int a, p
        if isinstance(x, (list, tuple)):
            if len(x) != self._degree:
                raise TypeError, "x must be a list of the right length"
            if coerce:
                R = parent.base_ring()
                p = R.order()
                for i from 0 <= i < self._degree:
                    a = int(R(x[i]))
                    self._entries[i] = a%p
            else:
                for i from 0 <= i < self._degree:
                    self._entries[i] = x[i]
            return
        if x != 0:
            raise TypeError, "can't initialize vector from nonzero non-list"
        else:
            for i from 0 <= i < self._degree:
                self._entries[i] = 0

    def __dealloc__(self):
        cdef Py_ssize_t i
        if self._entries:
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
        cdef mod_int l, r
        for i from 0 <= i < left.degree():
            l = left._entries[i]
            r = right._entries[i]
            if l < r:
                return -1
            elif l > r:
                return 1
        return 0

    def __len__(self):
        return self._degree

    def __setitem__(self, Py_ssize_t i, x):
        cdef IntegerMod_int n
        n = self.base_ring()(x)
        if i < 0 or i >= self._degree:
            raise IndexError
        else:
            self._entries[i] = n.ivalue

    def __getitem__(self, Py_ssize_t i):
        """
        Return the ith entry of self.
        """
        cdef IntegerMod_int n

        if i < 0 or i >= self._degree:
            raise IndexError, 'index out of range'
        else:
            n =  IntegerMod_int.__new__(IntegerMod_int)
            IntegerMod_abstract.__init__(n, self.base_ring())
            n.ivalue = self._entries[i]
            return n

    def __reduce__(self):
        return (unpickle_v0, (self._parent, self.list(), self._degree, self._p))

    cdef ModuleElement _add_c_impl(self, ModuleElement right):
        cdef Vector_modn_dense z, r
        r = right
        z = self._new_c()
        cdef Py_ssize_t i
        for i from 0 <= i < self._degree:
            z._entries[i] = (self._entries[i] + r._entries[i]) % self._p
        return z


    cdef ModuleElement _sub_c_impl(self, ModuleElement right):
        cdef Vector_modn_dense z, r
        r = right
        z = self._new_c()
        cdef Py_ssize_t i
        for i from 0 <= i < self._degree:
            z._entries[i] = (self._entries[i] - r._entries[i]) % self._p
        return z

    cdef Vector _vector_times_vector_c_impl(self, Vector right):
        cdef Vector_modn_dense z, r
        r = right
        z = self._new_c()
        cdef Py_ssize_t i
        for i from 0 <= i < self._degree:
            z._entries[i] = (self._entries[i] * r._entries[i]) % self._p
        return z

    cdef ModuleElement _rmul_c_impl(self, RingElement left):
        cdef IntegerMod_int a
        cdef Vector_modn_dense z

        a = left
        z = self._new_c()
        cdef Py_ssize_t i
        for i from 0 <= i < self._degree:
            z._entries[i] = (self._entries[i] * a.ivalue) % self._p
        return z

    cdef ModuleElement _lmul_c_impl(self, RingElement right):
        return self._rmul_c_impl(right)

    cdef ModuleElement _neg_c_impl(self):
        cdef Vector_modn_dense z
        z = self._new_c()
        cdef Py_ssize_t i
        for i from 0 <= i < self._degree:
            if self._entries[i] > 0:
                z._entries[i] = self._p - self._entries[i]
            else:
                z._entries[i] = 0
        return z



def unpickle_v0(parent, entries, degree, p):
    # If you think you want to change this function, don't.
    # Instead make a new version with a name like
    #    make_FreeModuleElement_generic_dense_v1
    # and changed the reduce method below.
    cdef Vector_modn_dense v
    v = PY_NEW(Vector_modn_dense)
    v._init(degree, parent, p)
    for i from 0 <= i < degree:
        v._entries[i] = entries[i]
    return v
