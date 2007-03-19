# choose: dense or sparse

cdef extern from "../libs/m4ri/m4ri.h":
    ctypedef struct packedmatrix:
        void *values
        int rows
        int cols
        int width # rounded as floor(rows/RADIX)+1
        int *rowswap

    ctypedef unsigned long mod_int

    ctypedef int bit

    cdef int RADIX

    ##############
    # Maintainance
    ##############

    cdef void setupPackingMasks()

    cdef void buildAllCodes() # builds all gray codes up to a certain size

    ##############
    # Constructors
    ##############

    # create empty matrix
    cdef packedmatrix *createPackedMatrix(int,int)

    # free memory
    cdef void destroyPackedMatrix(packedmatrix *)

    # filled uniformly random
    cdef void fillRandomlyPacked(packedmatrix *)

    # identity matrix
    cdef void makeIdentityPacked(packedmatrix *)

    # [A],[B] -> [AB]
    cdef packedmatrix *concatPacked(packedmatrix *A, packedmatrix *B)

    # returns a submatrix from a
    cdef packedmatrix *copySubMatrixPacked( packedmatrix *A, int lowr, int lowc, int highr, int highc)

    # deep copy
    cdef packedmatrix *clonePacked(packedmatrix *)

    ####
    # IO
    ####

    # print
    cdef void lazyPrint(packedmatrix *)

    # set bit
    cdef void writePackedCell( packedmatrix *m, int row, int col, bit value)

    # get bit
    cdef bit readPackedCell( packedmatrix *m, int row, int col )

    ######################
    # Low-Level Arithmetic
    ######################

    cdef void rowSwapPacked(packedmatrix *, int, int)

    cdef void rowClearPackedOffset(packedmatrix *m, int, int)

    cdef void rowAddPackedOffset(packedmatrix *m, int, int, int)

    ############
    # Arithmetic
    ############

    # naiv matrix multiply
    cdef packedmatrix *matrixTimesMatrixPacked(packedmatrix *, packedmatrix *)

    # M4RM matrix multiply
    cdef packedmatrix *m4rmPacked(packedmatrix *, packedmatrix *, int k)


    # matrix addition
    cdef packedmatrix *addPacked(packedmatrix *, packedmatrix *)

    # equality testing
    cdef int equalPackedMatrix(packedmatrix *, packedmatrix *)

    # returns -1,0,1
    cdef int comparePackedMatrix(packedmatrix *, packedmatrix *)

    # transpose
    cdef packedmatrix *transposePacked(packedmatrix *)

    # Gaussian Elimination
    cdef int gaussianPacked(packedmatrix *, int full)

    # M4RI
    cdef int simpleFourRussiansPackedFlex(packedmatrix *m, int full, int k)

    # Matrix Inversion
    cdef packedmatrix *invertPackedFlexRussian(packedmatrix *m, packedmatrix *identity, int k)


cimport matrix_dense

cdef class Matrix_mod2_dense(matrix_dense.Matrix_dense):
    cdef packedmatrix *_entries
    cdef object _one
    cdef object _zero


    cdef Matrix_mod2_dense _multiply_m4rm_c(Matrix_mod2_dense self, Matrix_mod2_dense right, int k)
