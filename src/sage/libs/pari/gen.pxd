include 'decl.pxi'

cimport sage.structure.element
cimport cython

@cython.final
cdef class gen(sage.structure.element.RingElement):
    cdef GEN g
    cdef object _refers_to
    cdef pari_sp b
    cdef void init(self, GEN g, pari_sp b)
    cdef GEN _gen(self)
    cdef gen new_gen(self, GEN x)
    cdef gen new_gen_noclear(self, GEN x)
    cdef gen pari(self, object x)
    cdef GEN _deepcopy_to_python_heap(self, GEN x, pari_sp* address)
    cdef long get_var(self, v)

cdef gen objtogen(object s)
