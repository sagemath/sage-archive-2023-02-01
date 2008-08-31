# choose: dense or sparse

cdef extern from "m4ri/m4ri.h":
    ctypedef unsigned long long word "word"

    ctypedef struct packedmatrix:
        word *values
        int nrows
        int ncols
        int width
        int *rowswap


    ctypedef int BIT

    cdef int RADIX

    ##############
    # Maintainance
    ##############
    cdef void m4ri_build_all_codes() # builds all gray codes up to a certain size
    cdef void m4ri_destroy_all_codes()

    ##############
    # Constructors
    ##############

    # create empty matrix
    cdef packedmatrix *mzd_init(int,int)

    # free memory
    cdef void mzd_free(packedmatrix *)

    # filled uniformly random
    cdef void mzd_randomize(packedmatrix *)

    # identity matrix if i%2
    cdef void mzd_set_ui(packedmatrix *, unsigned int i)

    # [A],[B] -> [AB]
    cdef packedmatrix *mzd_concat(packedmatrix *C, packedmatrix *A, packedmatrix *B)

    # [A],[B] -> | A |
    #            | B |
    cdef packedmatrix *mzd_stack(packedmatrix *C, packedmatrix *A, packedmatrix *B)


    # returns a submatrix from a
    cdef packedmatrix *mzd_submatrix(packedmatrix *S, packedmatrix *A, int lowr, int lowc, int highr, int highc)

    # return a matrix window to A
    cdef packedmatrix *mzd_init_window( packedmatrix *A, int lowr, int lowc, int highr, int highc)

    cdef void mzd_free_window(packedmatrix *)

    # deep copy
    cdef packedmatrix *mzd_copy(packedmatrix *, packedmatrix *)

    ####
    # IO
    ####

    # set BIT
    cdef void mzd_write_bit( packedmatrix *m, int row, int col, BIT value)

    # get BIT
    cdef BIT mzd_read_bit( packedmatrix *m, int row, int col )

    ######################
    # Low-Level Arithmetic
    ######################

    cdef void mzd_row_swap(packedmatrix *, int, int)

    cdef void mzd_col_swap(packedmatrix *, int, int)

    cdef void mzd_row_clear_offset(packedmatrix *m, int, int)

    cdef void mzd_row_add_offset(packedmatrix *m, int, int, int)

    ############
    # Arithmetic
    ############

    # matrix addition
    cdef packedmatrix *mzd_add(packedmatrix *, packedmatrix *, packedmatrix *)

    # naiv cubic matrix multiply
    cdef packedmatrix *mzd_mul_naiv(packedmatrix *, packedmatrix *, packedmatrix *)

    # matrix multiply using Gray codes
    cdef packedmatrix *mzd_mul_m4rm(packedmatrix *, packedmatrix *, packedmatrix *, int k)

    # matrix multiply using Gray codes (transposed)
    cdef packedmatrix *mzd_mul_m4rm_t(packedmatrix *, packedmatrix *, packedmatrix *, int k)

    # matrix multiply and addition using Gray codes: C = C + AB
    cdef packedmatrix *mzd_addmul_m4rm(packedmatrix *, packedmatrix *, packedmatrix *, int k)

    # matrix multiplication via Strassen's formula
    cdef packedmatrix *mzd_mul_strassen(packedmatrix *, packedmatrix *, packedmatrix *, int cutoff)

    # C = C + AB via Strassen's formula
    cdef packedmatrix *mzd_addmul_strassen(packedmatrix *, packedmatrix *, packedmatrix *, int cutoff)

    # equality testing
    cdef int mzd_equal(packedmatrix *, packedmatrix *)

    # returns -1,0,1
    cdef int mzd_cmp(packedmatrix *, packedmatrix *)

    # transpose
    cdef packedmatrix *mzd_transpose(packedmatrix *, packedmatrix *)

    # cubic Gaussian elimination
    cdef int mzd_reduce_naiv(packedmatrix *, int full)

    # matrix reduction using Gray codes
    cdef int mzd_reduce_m4ri(packedmatrix *m, int full, int k, packedmatrix *tablepacked, int *lookuppacked)

    # compute reduced row echelon form from upper triangular form
    cdef void mzd_top_reduce_m4ri(packedmatrix *m, int k, packedmatrix *tablepacked, int *lookuppacked)

    # matrix inversion using Gray codes
    cdef packedmatrix *mzd_invert_m4ri(packedmatrix *m, packedmatrix *identity, int k)


cimport matrix_dense

cdef class Matrix_mod2_dense(matrix_dense.Matrix_dense):
    cdef packedmatrix *_entries
    cdef object _one
    cdef object _zero


    cdef Matrix_mod2_dense _multiply_m4rm_c(Matrix_mod2_dense self, Matrix_mod2_dense right, int k)
