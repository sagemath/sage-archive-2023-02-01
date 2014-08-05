include "sage/ext/cdefs.pxi"
include "sage/libs/ntl/decl.pxi"

cimport matrix_dense
cimport matrix_integer_dense
cimport sage.rings.integer
from sage.rings.integer cimport Integer


#clib flint ntl

ctypedef long* GEN

cdef extern from "flint/fmpz.h":
    ctypedef void* fmpz_t
    void fmpz_init_set(fmpz_t f, const fmpz_t g)
    void fmpz_init_set_ui(fmpz_t f, const unsigned long g)

    void fmpz_init(fmpz_t f)
    void fmpz_abs ( fmpz_t f1 , const fmpz_t f2)
    void fmpz_set_ui(fmpz_t f, unsigned long val)
    void fmpz_set_mpz(fmpz_t f, const mpz_t val)
    void fmpz_set(fmpz_t f, const fmpz_t val)

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

cdef extern from "flint/fmpz_poly.h":
    ctypedef void* fmpz_poly_t
    void fmpz_poly_get_coeff_fmpz ( fmpz_t x , const fmpz_poly_t
poly , long n )
    long fmpz_poly_degree ( const fmpz_poly_t poly )
    void fmpz_poly_init ( fmpz_poly_t poly)    
cdef extern from "flint/fmpz_mat.h":
    ctypedef void* fmpz_mat_t

    void fmpz_mat_init(fmpz_mat_t mat,unsigned long rows,unsigned long cols)
    void fmpz_mat_init_set ( fmpz_mat_t mat , const fmpz_mat_t src )
    int fmpz_mat_print_pretty( const fmpz_mat_t mat )
    fmpz_t fmpz_mat_entry(fmpz_mat_t mat ,long i ,long j)
    void fmpz_mat_zero( fmpz_mat_t mat )
    void fmpz_mat_one( fmpz_mat_t mat )
    void fmpz_mat_scalar_mul_si( fmpz_mat_t B , const fmpz_mat_t A , long c )
    void fmpz_mat_scalar_mul_fmpz( fmpz_mat_t B , const fmpz_mat_t A , const fmpz_t c )
    void fmpz_mat_mul( fmpz_mat_t C , const fmpz_mat_t A , const  fmpz_mat_t B )
    void fmpz_mat_sqr( fmpz_mat_t B , const  fmpz_mat_t A )
    void fmpz_mat_add( fmpz_mat_t C , const fmpz_mat_t A , const  fmpz_mat_t B )
    void fmpz_mat_sub( fmpz_mat_t C , const fmpz_mat_t A , const  fmpz_mat_t B )
    void fmpz_mat_pow( fmpz_mat_t C , const fmpz_mat_t A , unsigned long n )
    void fmpz_mat_clear(fmpz_mat_t mat)
    void fmpz_mat_set(fmpz_mat_t result, fmpz_mat_t mat)
    void fmpz_mat_zero(fmpz_mat_t mat)
    int fmpz_mat_is_zero ( const fmpz_mat_t mat )
    void fmpz_mat_charpoly ( fmpz_poly_t cp , const fmpz_mat_t mat )
    
cdef class Matrix_integer_dense_flint(matrix_dense.Matrix_dense):
    cdef char _initialized
    cdef object _pivots
    cdef fmpz_mat_t _matrix_flint

    cdef int mpz_height(self, mpz_t height) except -1
#     cpdef double _log_avg_sq1(self) except -1.0
#     cdef _mod_int_c(self, mod_int modulus)
#    cdef _mod_two(self)

    cdef _zero_out_matrix(self)
    cdef _new_unitialized_matrix(self, Py_ssize_t nrows, Py_ssize_t ncols)
    cdef _pickle_version0(self)
    cdef _unpickle_version0(self, data)
    cpdef _export_as_string(self, int base=?)

#     cdef void reduce_entry_unsafe(self, Py_ssize_t i, Py_ssize_t j, Integer modulus)

#    cdef _new_uninitialized_matrix(self, Py_ssize_t nrows, Py_ssize_t ncols)

