from sage.structure.element cimport Element, ModuleElement, RingElement


cdef class LinearTensor(ModuleElement):
    cpdef dict _f
    cpdef ModuleElement _add_(self, ModuleElement b)
    cpdef ModuleElement _sub_(self, ModuleElement b)
    cpdef ModuleElement _neg_(self)
    cpdef ModuleElement _lmul_(self, RingElement b)
    cpdef ModuleElement _rmul_(self, RingElement b)
    cdef _richcmp(left, right, int op)
    cdef int _cmp_c_impl(left, Element right) except -2

    
