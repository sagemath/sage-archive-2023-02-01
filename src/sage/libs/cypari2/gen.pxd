from .types cimport *
cimport cython


cdef class Gen_auto:
    cdef GEN g
    cdef pari_sp b
    cdef dict refers_to

@cython.final
cdef class Gen(Gen_auto):
    pass

cdef Gen new_ref(GEN g, Gen parent)
cpdef Gen objtogen(s)
