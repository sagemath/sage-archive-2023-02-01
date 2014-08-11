include 'decl.pxi'


#clib flint ntl

cimport sage.structure.parent_base
cimport cython

from sage.libs.pari.gen cimport gen

cpdef long prec_bits_to_words(unsigned long prec_in_bits)


cdef extern from "flint/fmpz.h":
    ctypedef void* fmpz_t
    void fmpz_init_set(fmpz_t f, const fmpz_t g)
    void fmpz_init_set_ui(fmpz_t f, const unsigned long g)
    size_t fmpz_size ( const fmpz_t f )
    void fmpz_init(fmpz_t f)
    void fmpz_abs ( fmpz_t f1 , const fmpz_t f2)
    void fmpz_set_ui(fmpz_t f, unsigned long val)
    void fmpz_set_mpz(fmpz_t f, const mpz_t val)
    void fmpz_set(fmpz_t f, const fmpz_t val)
    int fmpz_sgn ( const fmpz_t f )

    void fmpz_get_mpz(mpz_t x, const fmpz_t f)
    void fmpz_fdiv_r( fmpz_t f , const fmpz_t g , const fmpz_t h )
    int fmpz_print(fmpz_t x)
    size_t fmpz_sizeinbase (const fmpz_t f , int b)
    char * fmpz_get_str ( char * str , int b , const fmpz_t f)
    int fmpz_set_str( fmpz_t f , const char * str , int b)
    int fmpz_cmp ( const fmpz_t f , const fmpz_t g)
    void fmpz_set_si ( fmpz_t f , long val )
    void fmpz_addmul ( fmpz_t f , const fmpz_t g , const fmpz_t h )
    void fmpz_clear ( fmpz_t f)

cdef extern from "flint/fmpz_mat.h":
    ctypedef void* fmpz_mat_t
    fmpz_t fmpz_mat_entry(fmpz_mat_t mat ,long i ,long j)
    
    
@cython.final
cdef class PariInstance(sage.structure.parent_base.ParentWithBase):
    cdef long _real_precision
    cdef gen PARI_ZERO, PARI_ONE, PARI_TWO
    cdef inline gen new_gen(self, GEN x)
    cdef inline gen new_gen_noclear(self, GEN x)
    cdef gen new_gen_from_mpz_t(self, mpz_t value)
    cdef inline GEN _new_GEN_from_mpz_t(self, mpz_t value)
    # cdef gen new_gen_from_fmpz_t(self, fmpz_t value)
    # cdef inline GEN _new_GEN_from_fmpz_t(self, fmpz_t value)
    cdef gen new_gen_from_mpq_t(self, mpq_t value)
    cdef inline GEN _new_GEN_from_mpq_t(self, mpq_t value)
    cdef gen new_gen_from_int(self, int value)
    cdef gen new_t_POL_from_int_star(self, int *vals, int length, long varnum)
    cdef gen new_gen_from_padic(self, long ordp, long relprec, mpz_t prime, mpz_t p_pow, mpz_t unit)
    cdef inline void clear_stack(self)
    cdef gen double_to_gen_c(self, double)
    cdef GEN double_to_GEN(self, double)
    cdef GEN deepcopy_to_python_heap(self, GEN x, pari_sp* address)
    cdef gen new_ref(self, GEN g, gen parent)
    cdef gen _empty_vector(self, long n)
    cdef long get_var(self, v)
    cdef GEN _new_GEN_from_mpz_t_matrix(self, mpz_t** B, Py_ssize_t nr, Py_ssize_t nc)
    cdef GEN _new_GEN_from_mpz_t_matrix_rotate90(self, mpz_t** B, Py_ssize_t nr, Py_ssize_t nc)
    cdef gen integer_matrix(self, mpz_t** B, Py_ssize_t nr, Py_ssize_t nc, bint permute_for_hnf)
    cdef GEN _new_GEN_from_fmpz_mat_t(self, fmpz_mat_t B, Py_ssize_t nr, Py_ssize_t nc)
    cdef GEN _new_GEN_from_fmpz_mat_t_rotate90(self, fmpz_mat_t B, Py_ssize_t nr, Py_ssize_t nc)
    cdef gen integer_matrix_from_flint(self, fmpz_mat_t B, Py_ssize_t nr, Py_ssize_t nc, bint permute_for_hnf)
    
    cdef GEN _new_GEN_from_mpq_t_matrix(self, mpq_t** B, Py_ssize_t nr, Py_ssize_t nc)
    cdef gen rational_matrix(self, mpq_t** B, Py_ssize_t nr, Py_ssize_t nc)
