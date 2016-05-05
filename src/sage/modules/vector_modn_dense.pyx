"""
Vectors with integer mod n entries, with n small.

EXAMPLES::

    sage: v = vector(Integers(8),[1,2,3,4,5])
    sage: type(v)
    <type 'sage.modules.vector_modn_dense.Vector_modn_dense'>
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
    7

    sage: v = vector(Integers(8),[1,2,3,4,5])
    sage: u = vector(Integers(8),[1,2,3,4,4])
    sage: v - u
    (0, 0, 0, 0, 1)
    sage: u - v
    (0, 0, 0, 0, 7)

    sage: v = vector((Integers(5)(1),2,3,4,4))
    sage: u = vector((Integers(5)(1),2,3,4,3))
    sage: v - u
    (0, 0, 0, 0, 1)
    sage: u - v
    (0, 0, 0, 0, 4)

We make a large zero vector::

    sage: k = Integers(8)^100000; k
    Ambient free module of rank 100000 over Ring of integers modulo 8
    sage: v = k(0)
    sage: v[:10]
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

We multiply a vector by a matrix::

    sage: a = (GF(97)^5)(range(5))
    sage: m = matrix(GF(97),5,range(25))
    sage: a*m
    (53, 63, 73, 83, 93)

TESTS::

    sage: v = vector(Integers(8), [1,2,3,4,5])
    sage: loads(dumps(v)) == v
    True
    sage: v = vector(Integers(389), [1,2,3,4,5])
    sage: loads(dumps(v)) == v
    True
    sage: v = vector(Integers(next_prime(10^20)), [1,2,3,4,5])
    sage: loads(dumps(v)) == v
    True

    sage: K = GF(previous_prime(2^31))
    sage: v = vector(K, [42]);  type(v[0])
    <type 'sage.rings.finite_rings.integer_mod.IntegerMod_int64'>
    sage: ~v[0]
    2096353084

    sage: K = GF(next_prime(2^31))
    sage: v = vector(K, [42]);  type(v[0])
    <type 'sage.rings.finite_rings.integer_mod.IntegerMod_gmp'>
    sage: ~v[0]
    1482786336

    sage: w = vector(GF(11), [-1,0,0,0])
    sage: w.set_immutable()
    sage: isinstance(hash(w), int)
    True

AUTHOR:

- William Stein (2007)
"""

#*****************************************************************************
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include "cysignals/memory.pxi"

from sage.rings.finite_rings.stdint cimport INTEGER_MOD_INT64_LIMIT

MAX_MODULUS = INTEGER_MOD_INT64_LIMIT

from sage.rings.finite_rings.integer_mod cimport (
    IntegerMod_int, IntegerMod_int64,
    IntegerMod_abstract, use_32bit_type)

cdef mod_int ivalue(IntegerMod_abstract x) except -1:
    if type(x) is IntegerMod_int:
        return (<IntegerMod_int>x).ivalue
    elif type(x) is IntegerMod_int64:
        return (<IntegerMod_int64>x).ivalue
    else:
        raise TypeError, "non-fixed size integer"

from sage.structure.element cimport Element, ModuleElement, RingElement, Vector

cimport free_module_element
from free_module_element import vector

cdef class Vector_modn_dense(free_module_element.FreeModuleElement):
    cdef _new_c(self):
        cdef Vector_modn_dense y
        y = Vector_modn_dense.__new__(Vector_modn_dense)
        y._init(self._degree, self._parent, self._p)
        return y

    cdef bint is_dense_c(self):
        return 1

    cdef bint is_sparse_c(self):
        return 0

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
        self._entries = <mod_int *>check_allocarray(degree, sizeof(mod_int))

    def __cinit__(self, parent=None, x=None, coerce=True, copy=True):
        self._entries = NULL
        self._is_mutable = 1
        if not parent is None:
            self._init(parent.degree(), parent, parent.base_ring().order())

    def __init__(self, parent, x, coerce=True, copy=True):
        cdef Py_ssize_t i
        cdef mod_int a, p
        if isinstance(x, (list, tuple)):
            if len(x) != self._degree:
                raise TypeError("x must be a list of the right length")
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
            raise TypeError("can't initialize vector from nonzero non-list")
        else:
            for i from 0 <= i < self._degree:
                self._entries[i] = 0

    def __dealloc__(self):
        sig_free(self._entries)

    cpdef int _cmp_(left, Element right) except -2:
        """
        EXAMPLES::

            sage: v = vector(GF(5), [0,0,0,0])
            sage: v == 0
            True
            sage: v == 1
            False
            sage: v == v
            True
            sage: w = vector(GF(11), [-1,0,0,0])
            sage: w == w
            True
        """
        cdef Py_ssize_t i
        cdef mod_int l, r
        for i from 0 <= i < left.degree():
            l = left._entries[i]
            r = (<Vector_modn_dense>right)._entries[i]
            if l < r:
                return -1
            elif l > r:
                return 1
        return 0

    cdef get_unsafe(self, Py_ssize_t i):
        """
        EXAMPLES::

            sage: R = Integers(7)
            sage: v = vector(R, [1,2,3]); v
            (1, 2, 3)
            sage: v[0]
            1
            sage: v[2]
            3
            sage: v[-2]
            2
            sage: v[0:2]
            (1, 2)
            sage: v[5]
            Traceback (most recent call last):
            ...
            IndexError: vector index out of range
            sage: v[-5]
            Traceback (most recent call last):
            ...
            IndexError: vector index out of range
        """
        cdef IntegerMod_int n
        cdef IntegerMod_int64 m

        if use_32bit_type(self._p):
            n = IntegerMod_int.__new__(IntegerMod_int)
            IntegerMod_abstract.__init__(n, self.base_ring())
            n.ivalue = self._entries[i]
            return n
        else:
            m = IntegerMod_int64.__new__(IntegerMod_int64)
            IntegerMod_abstract.__init__(m, self.base_ring())
            m.ivalue = self._entries[i]
            return m

    cdef int set_unsafe(self, Py_ssize_t i, value) except -1:
        """
        EXAMPLES::

            sage: R = Integers(7)
            sage: v = vector(R, [1,2,3]); v
            (1, 2, 3)
            sage: v[0] = 7^7
            sage: v
            (0, 2, 3)
        """
        self._entries[i] = ivalue(<IntegerMod_abstract>value)


    def __reduce__(self):
        return unpickle_v1, (self._parent, self.list(), self._degree, self._p, self._is_mutable)

    cpdef ModuleElement _add_(self, ModuleElement right):
        cdef Vector_modn_dense z, r
        r = right
        z = self._new_c()
        cdef Py_ssize_t i
        for i from 0 <= i < self._degree:
            z._entries[i] = (self._entries[i] + r._entries[i]) % self._p
        return z


    cpdef ModuleElement _sub_(self, ModuleElement right):
        cdef Vector_modn_dense z, r
        r = right
        z = self._new_c()
        cdef Py_ssize_t i
        for i from 0 <= i < self._degree:
            z._entries[i] = (self._p + self._entries[i] - r._entries[i]) % self._p
        return z

    cpdef Element _dot_product_(self, Vector right):
        cdef Py_ssize_t i
        cdef IntegerMod_int n
        cdef Vector_modn_dense r = right
        n =  IntegerMod_int.__new__(IntegerMod_int)
        IntegerMod_abstract.__init__(n, self.base_ring())
        n.ivalue = 0

        for i from 0 <= i < self._degree:
            n.ivalue = (n.ivalue + self._entries[i] * r._entries[i]) % self._p

        return n

    cpdef Vector _pairwise_product_(self, Vector right):
        """
        EXAMPLES:
           sage: v = vector(Integers(8), [2,3]); w = vector(Integers(8), [2,5])
           sage: v * w
           3
           sage: w * v
           3
        """
        cdef Vector_modn_dense z, r
        r = right
        z = self._new_c()
        cdef Py_ssize_t i
        for i from 0 <= i < self._degree:
            z._entries[i] = (self._entries[i] * r._entries[i]) % self._p
        return z

    cpdef ModuleElement _rmul_(self, RingElement left):
        cdef Vector_modn_dense z

        cdef mod_int a = ivalue(left)
        z = self._new_c()
        cdef Py_ssize_t i

        for i from 0 <= i < self._degree:
            z._entries[i] = (self._entries[i] * a) % self._p
        return z

    cpdef ModuleElement _lmul_(self, RingElement right):
        return self._rmul_(right)

    cpdef ModuleElement _neg_(self):
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
    v = Vector_modn_dense.__new__(Vector_modn_dense)
    v._init(degree, parent, p)
    for i from 0 <= i < degree:
        v._entries[i] = entries[i]
    return v

def unpickle_v1(parent, entries, degree, p, is_mutable):
    cdef Vector_modn_dense v
    v = Vector_modn_dense.__new__(Vector_modn_dense)
    v._init(degree, parent, p)
    for i from 0 <= i < degree:
        v._entries[i] = entries[i]
    v._is_mutable = is_mutable
    return v
