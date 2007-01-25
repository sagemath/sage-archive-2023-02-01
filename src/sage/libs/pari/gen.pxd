include 'decl.pxi'

cdef class gen:
    cdef GEN g
    cdef object refers_to
    cdef pari_sp b
    cdef void init(self, GEN g, pari_sp b)
    cdef GEN _gen(self)
    cdef gen new_gen(self, GEN x)
    cdef gen new_gen_noclear(self, GEN x)
    cdef gen pari(self, object x)
    cdef GEN _deepcopy_to_python_heap(self, GEN x, pari_sp* address)
    cdef int get_var(self, v)

cdef class PariInstance:
    cdef gen ZERO, ONE, TWO
    cdef gen new_gen(self, GEN x)
    cdef gen new_gen_noclear(self, GEN x)
    cdef GEN deepcopy_to_python_heap(self, GEN x, pari_sp* address)
    cdef gen new_ref(self, GEN x, g)
    cdef gen adapt(self, s)
    cdef int get_var(self, v)
    cdef object GEN_to_str(self, GEN g)
    cdef GEN toGEN(self, x) except NULL




