
from sage.structure.element cimport Element

cdef class ElementWrapper(Element):
    cdef public object value

    cpdef bint _lt_by_value(self, other)
    cpdef int _cmp_by_value(self, other)

cdef class ElementWrapperCheckWrappedClass(ElementWrapper):
    pass

