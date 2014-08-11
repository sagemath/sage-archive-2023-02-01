include "sage/ext/cdefs.pxi"
include "sage/libs/ntl/decl.pxi"

cimport matrix_dense
cimport matrix_integer_dense
cimport sage.rings.integer
from sage.rings.integer cimport Integer
from sage.ext.mod_int cimport *

cdef extern from "flint/fmpz.h":
    ctypedef void* fmpz_t
    void fmpz_init_set(fmpz_t f, const fmpz_t g)
    void fmpz_init_set_ui(fmpz_t f, const unsigned long g)

    void fmpz_init(fmpz_t f)
    void fmpz_abs ( fmpz_t f1 , const fmpz_t f2)
    void fmpz_set_ui(fmpz_t f, unsigned long val)
    void fmpz_set_mpz(fmpz_t f, const mpz_t val)
    void fmpz_set(fmpz_t f, const fmpz_t val)
    void fmpz_add ( fmpz_t f , const fmpz_t g , const fmpz_t h)

    void fmpz_sub ( fmpz_t f , const fmpz_t g , const fmpz_t h)
        
    void fmpz_addmul ( fmpz_t f , const fmpz_t g , const fmpz_t h )
    void fmpz_mul ( fmpz_t f , const fmpz_t g , const fmpz_t h)
    int fmpz_sgn ( const fmpz_t f )
    void fmpz_get_mpz(mpz_t x, const fmpz_t f)
    void fmpz_fdiv_r( fmpz_t f , const fmpz_t g , const fmpz_t h )
    unsigned long fmpz_fdiv_ui ( const fmpz_t g , unsigned long x)    
    void fmpz_divexact ( fmpz_t f , const fmpz_t g , const fmpz_t h)
    void fmpz_gcd ( fmpz_t f , const fmpz_t g , const fmpz_t h)
    
    int fmpz_print(fmpz_t x)
    size_t fmpz_sizeinbase (const fmpz_t f , int b)
    char * fmpz_get_str ( char * str , int b , const fmpz_t f)
    int fmpz_set_str( fmpz_t f , const char * str , int b)
    int fmpz_cmp ( const fmpz_t f , const fmpz_t g)
    int fmpz_cmp_si ( const fmpz_t f , int g)
    int fmpz_cmp_ui ( const fmpz_t f , unsigned int g)

    
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
    ctypedef long slong
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
    long fmpz_mat_rref(fmpz_mat_t R, fmpz_t den, const fmpz_mat_t A)
    # void fmpz_mat_multi_CRT_ui ( fmpz_mat_t mat , nmod_mat_t *const residues , slong nres , int sign )
    void fmpz_mat_det ( fmpz_t det , const fmpz_mat_t A )
    void fmpz_mat_transpose ( fmpz_mat_t B , const fmpz_mat_t A)
    long fmpz_mat_rank ( const fmpz_mat_t A)
    int fmpz_mat_solve ( fmpz_mat_t X , fmpz_t den , const fmpz_mat_t A , const fmpz_mat_t B )

cdef class Matrix_rational_dense(matrix_dense.Matrix_dense):

    cdef mpq_t tmp
    cdef mpq_t *_entries
    cdef mpq_t ** _matrix
    cdef object __pivots

    cdef int mpz_denom(self, mpz_t d) except -1
    cdef int mpz_height(self, mpz_t height) except -1
    cdef int _rescale(self, mpq_t a) except -1

    cdef _pickle_version0(self)
    cdef _unpickle_version0(self, data)
    cpdef _export_as_string(self, int base=?)

    cdef _add_ui_unsafe_assuming_int(self, Py_ssize_t i, Py_ssize_t j, unsigned long int n)
    cdef _sub_ui_unsafe_assuming_int(self, Py_ssize_t i, Py_ssize_t j, unsigned long int n)

cdef class MatrixWindow:
    cdef Matrix_rational_dense _matrix
    cdef int _row, _col, _nrows, _ncols

################################################################
# fast conversion to pari on the stack
################################################################
ctypedef long* GEN
cdef inline GEN pari_GEN(Matrix_rational_dense B)
