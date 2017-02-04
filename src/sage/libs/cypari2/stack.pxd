from .types cimport GEN, pari_sp
from .gen cimport Gen

cdef void clear_stack()
cdef GEN deepcopy_to_python_heap(GEN x, pari_sp* address)
cdef Gen new_gen(GEN x)
cdef Gen new_gen_noclear(GEN x)
