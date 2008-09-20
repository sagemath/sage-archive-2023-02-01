include 'decl.pxi'

cimport sage.structure.element

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
    cdef int get_var(self, v)

cimport sage.structure.parent_base

cdef class PariInstance(sage.structure.parent_base.ParentWithBase):
    cdef gen ZERO, ONE, TWO
    cdef gen new_gen(self, GEN x)
    cdef gen new_gen_noclear(self, GEN x)
    cdef gen new_gen_from_mpz_t(self, mpz_t value)
    cdef GEN new_GEN_from_mpz_t(self, mpz_t value)
    cdef gen new_gen_from_int(self, int value)
    cdef gen new_t_POL_from_int_star(self, int *vals, int length, long varnum)
    cdef gen new_gen_from_padic(self, long ordp, long relprec, mpz_t prime, mpz_t p_pow, mpz_t unit)
    cdef void clear_stack(self)
    cdef gen double_to_gen_c(self, double)
    cdef GEN double_to_GEN(self, double)
    cdef GEN deepcopy_to_python_heap(self, GEN x, pari_sp* address)
    cdef gen new_ref(self, GEN g, gen parent)
    cdef gen _empty_vector(self, long n)
    cdef int get_var(self, v)
    cdef object GEN_to_str(self, GEN g)
    cdef GEN toGEN(self, x, int i) except NULL




