include 'decl.pxi'

cimport sage.structure.element
cimport cython

@cython.final
cdef class gen(sage.structure.element.RingElement):
    cdef GEN g
    cdef pari_sp b
    cdef dict refers_to

cdef gen objtogen(object s)
