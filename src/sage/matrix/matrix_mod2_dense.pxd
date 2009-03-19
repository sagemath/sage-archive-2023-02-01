# choose: dense or sparse

cdef extern from "m4ri/m4ri.h":
    ctypedef unsigned long long word "word"

    ctypedef struct mzd_t:
        int nrows
        int ncols
        int width
        word **rows

    ctypedef struct mzp_t:
        int *values
        int size

    ctypedef int BIT

    cdef int RADIX

    ##############
    # Maintainance
    ##############

    # builds all gray codes up to a certain size
    cdef void m4ri_build_all_codes()
    cdef void m4ri_destroy_all_codes()

    ##############
    # Constructors
    ##############

    # create empty matrix
    cdef mzd_t *mzd_init(int,int)

    # create the identity permutation
    cdef mzp_t *mzp_init(long)

    # free memory for the matrix
    cdef void mzd_free(mzd_t *)

    # free memory for the permutation
    cdef void mzp_free(mzp_t *)

    # filled uniformly random
    cdef void mzd_randomize(mzd_t *)

    # identity matrix if i%2
    cdef void mzd_set_ui(mzd_t *, unsigned int i)

    # [A],[B] -> [AB]
    cdef mzd_t *mzd_concat(mzd_t *C, mzd_t *A, mzd_t *B)

    # [A],[B] -> | A |
    #            | B |
    cdef mzd_t *mzd_stack(mzd_t *C, mzd_t *A, mzd_t *B)

    # returns a submatrix from a
    cdef mzd_t *mzd_submatrix(mzd_t *S, mzd_t *A, int lowr, int lowc, int highr, int highc)

    # return a matrix window to A
    cdef mzd_t *mzd_init_window( mzd_t *A, int lowr, int lowc, int highr, int highc)

    cdef void mzd_free_window(mzd_t *)

    # deep copy
    cdef mzd_t *mzd_copy(mzd_t *, mzd_t *)

    ##############
    # Bit Level IO
    ##############

    # set BIT
    cdef void mzd_write_bit( mzd_t *m, int row, int col, BIT value)

    # get BIT
    cdef BIT mzd_read_bit( mzd_t *m, int row, int col )

    # get BITs (n<=64)
    cdef word mzd_read_bits( mzd_t *m, int row, int col, int n)

    #####################
    # Row/Column Based IO
    #####################

    cdef void mzd_row_swap(mzd_t *, int, int)

    cdef void mzd_col_swap(mzd_t *, int, int)

    cdef void mzd_row_clear_offset(mzd_t *m, int, int)

    cdef void mzd_row_add_offset(mzd_t *m, int, int, int)

    ############
    # Arithmetic
    ############

    # matrix addition
    cdef mzd_t *mzd_add(mzd_t *, mzd_t *, mzd_t *)

    # naive cubic matrix multiply
    cdef mzd_t *mzd_mul_naive(mzd_t *, mzd_t *, mzd_t *)

    # matrix multiply using Gray codes
    cdef mzd_t *mzd_mul_m4rm(mzd_t *, mzd_t *, mzd_t *, int k)

    # matrix multiply using Gray codes (transposed)
    cdef mzd_t *mzd_mul_m4rm_t(mzd_t *, mzd_t *, mzd_t *, int k)

    # matrix multiply and addition using Gray codes: C = C + AB
    cdef mzd_t *mzd_addmul_m4rm(mzd_t *, mzd_t *, mzd_t *, int k)

    # matrix multiplication via Strassen's formula
    cdef mzd_t *mzd_mul(mzd_t *, mzd_t *, mzd_t *, int cutoff)

    # C = C + AB via Strassen's formula
    cdef mzd_t *mzd_addmul(mzd_t *, mzd_t *, mzd_t *, int cutoff)

    # equality testing
    cdef int mzd_equal(mzd_t *, mzd_t *)

    # returns -1,0,1
    cdef int mzd_cmp(mzd_t *, mzd_t *)

    # transpose
    cdef mzd_t *mzd_transpose(mzd_t *, mzd_t *)

    # density with given resolution
    cdef double mzd_density(mzd_t *, int resolution)

    ########################
    # LAPACK Level Functions
    ########################

    # cubic Gaussian elimination
    cdef int mzd_echelonize_naive(mzd_t *, int full)

    # row echelon form using Gray codes
    cdef int mzd_echelonize_m4ri(mzd_t *m, int full, int k)

    # reduced row echelon form from upper triangular form
    cdef void mzd_top_echelonize_m4ri(mzd_t *m, int k)

    # matrix inversion using Gray codes
    cdef mzd_t *mzd_invert_m4ri(mzd_t *m, mzd_t *identity, int k)

    # assymptotically fast PLUQ factorization
    cdef long mzd_pluq(mzd_t *A, mzp_t *P, mzp_t *Q, int cutoff)

    # PLUQ factorization using Gray codes
    cdef long _mzd_pluq_mmpf(mzd_t *A, mzp_t *P, mzp_t *Q, int k)

    # cubic PLUQ factorization
    cdef long _mzd_pluq_naive(mzd_t *A, mzp_t *P, mzp_t *Q)

    # reduced row echelon form using PLUQ factorization
    cdef long mzd_echelonize_pluq(mzd_t *A, int full)

    ##################################
    # Internal Functions for Debugging
    ##################################

    cdef void mzd_write_zeroed_bits(mzd_t *m, int x, int y, int n, word values)

    cdef void mzd_clear_bits(mzd_t *m, int x, int y, int n)


cimport matrix_dense

cdef class Matrix_mod2_dense(matrix_dense.Matrix_dense):
    cdef mzd_t *_entries
    cdef object _one
    cdef object _zero

    cdef set_unsafe_int(self, Py_ssize_t i, Py_ssize_t j, int value)

    cpdef Matrix_mod2_dense _multiply_m4rm(Matrix_mod2_dense self, Matrix_mod2_dense right, int k)
    cpdef Matrix_mod2_dense _multiply_strassen(Matrix_mod2_dense self, Matrix_mod2_dense right, int cutoff)

    # For conversion to other systems
    cpdef _export_as_string(self)
