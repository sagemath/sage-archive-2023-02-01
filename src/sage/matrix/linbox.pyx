
## NOTE: The _sig_on/_sig_off stuff can't go in here -- it has to be in the
## code that calls these functions.  Otherwise strangely objects get left
## in an incorrect state.

include "../ext/interrupt.pxi"

cdef extern from "matrix_modn_dense_linbox.h":
    void linbox_modn_dense_delete_array(mod_int *f)

    int linbox_modn_dense_echelonize(unsigned long modulus,
                              mod_int **matrix, size_t nrows, size_t ncols)
    void linbox_modn_dense_minpoly(unsigned long modulus, mod_int **mp, size_t* degree, size_t n,
                                   mod_int **matrix, int do_minpoly)

    int  linbox_modn_dense_matrix_matrix_multiply(unsigned long modulus, mod_int **ans,
                                                  mod_int **A, mod_int **B,
                                                  size_t A_nr, size_t A_nc,
                                                  size_t B_nr, size_t B_nc)

    int linbox_modn_dense_rank(unsigned long modulus,
                               mod_int** matrix, size_t nrows, size_t ncols)


##########################################################################
## Dense matrices modulo p
##########################################################################
cdef class Linbox_modn_dense:
    cdef set(self, mod_int n, mod_int** matrix,
             size_t nrows, size_t ncols):
        self.n = n
        self.nrows = nrows
        self.ncols = ncols
        self.matrix = matrix

    cdef int echelonize(self):
        cdef int r
        r = linbox_modn_dense_echelonize(self.n, self.matrix,
                                         self.nrows, self.ncols)
        return r

    cdef poly(self, minpoly):
        """
        INPUT:
            as given

        OUTPUT:
            coefficients of minpoly as a Python list
        """
        cdef mod_int *f
        cdef size_t degree
        linbox_modn_dense_minpoly(self.n, &f, &degree,
                                  self.nrows, self.matrix,
                                  minpoly)
        v = []
        cdef Py_ssize_t i
        for i from 0 <= i <= degree:
            v.append(f[i])
        linbox_modn_dense_delete_array(f)
        return v

    cdef matrix_matrix_multiply(self,
                                mod_int **ans,
                                mod_int **B,
                                size_t B_nr, size_t B_nc):
        cdef int e
        e = linbox_modn_dense_matrix_matrix_multiply(self.n, ans,
                                                     self.matrix,  B,
                                                     self.nrows, self.ncols,
                                                     B_nr, B_nc)
        if e:
            raise RuntimError, "error doing matrix matrix multiply modn using linbox"


    cdef unsigned long rank(self) except -1:
        cdef unsigned long r
        r = linbox_modn_dense_rank(self.n,   self.matrix, self.nrows, self.ncols)
        return r


##########################################################################
## Sparse matices modulo p.
##########################################################################
cdef class Linbox_modn_sparse:
    pass


##########################################################################
## Sparse matrices over ZZ
##########################################################################


##########################################################################
## Dense matrices over ZZ
##########################################################################

