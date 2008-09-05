from sage.structure.sage_object cimport SageObject

cdef class Generators(SageObject):

    cdef readonly _obj
    cdef readonly _index_set
    cdef readonly _category

    cpdef get_from_index(self, i)
    cpdef index_set(self)
    cpdef category(self)
    cpdef obj(self)
    cpdef count(self)
    cpdef list(self)

cdef class Generators_finite(Generators):
    cdef Py_ssize_t _n

cdef class Generators_list(Generators_finite):
    cdef _List

cdef class Generators_naturals(Generators):
    pass

cdef class Generators_none(Generators):
    pass

cdef class Generators_lazy_all(Generators):
    cdef _f
    cdef int _compute_gens(self) except -1

