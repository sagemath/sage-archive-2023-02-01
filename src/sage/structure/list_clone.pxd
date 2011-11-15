#*****************************************************************************
#  Copyright (C) 2009-2010 Florent Hivert <Florent.Hivert@univ-rouen.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.element cimport Element


cdef class ClonableElement(Element):
    cdef bint _is_immutable
    cdef bint _needs_check
    cdef long int  _hash

    cpdef inline bint _require_mutable(self) except -2
    cpdef inline bint is_mutable(self)
    cpdef inline bint is_immutable(self)
    cpdef inline set_immutable(self)

    cpdef inline _set_mutable(self)

    cpdef inline ClonableElement clone(self, bint check=?)
    cpdef inline ClonableElement __enter__(self)
    cpdef inline bint __exit__(self, typ, value, tracback) except -2

cdef class ClonableArray(ClonableElement):
    cdef list _list

    cpdef inline list _get_list(self)
    cpdef inline _set_list(self, list lst)
    cpdef ClonableArray __copy__(self)
    cpdef check(self)
    cpdef inline object _getitem(self, int key)
    cpdef inline _setitem(self, int key, value)
    cpdef int index(self, key, start=*, stop=*) except -1
    cpdef int count(self, key) except -1
    cpdef long int _hash_(self) except 0

cdef class ClonableList(ClonableArray):
    cpdef append(self, el)
    cpdef extend(self, it)
    cpdef insert(self, int index, el)
    cpdef pop(self, int index=*)
    cpdef remove(self, el)

cdef class ClonableIntArray(ClonableElement):
    cdef int _len
    cdef int* _list

    cpdef _alloc_(self, int size)
    cpdef ClonableIntArray __copy__(self)
    cpdef check(self)
    cpdef inline object _getitem(self, int key)
    cpdef inline _setitem(self, int item, value)
    cpdef int index(self, int item) except -1
    cpdef long int _hash_(self) except 0
    cpdef inline list list(self)
