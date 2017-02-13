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
cdef Gen list_of_Gens_to_Gen(list s)
cpdef Gen objtogen(s)
