
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


cdef class Linbox:
    cdef int modn_dense_echelonize(self, mod_int modulus, mod_int** matrix,
                                   size_t nrows, size_t ncols):
        cdef int r
        _sig_on
        r = linbox_modn_dense_echelonize(modulus, matrix, nrows, ncols)
        _sig_off
        return r

    cdef modn_dense_poly(self, unsigned long modulus, size_t n,
                         mod_int **matrix, minpoly):
        """
        INPUT:
            as given

        OUTPUT:
            coefficients of minpoly as a Python list
        """
        cdef mod_int *f
        cdef size_t degree
        _sig_on
        linbox_modn_dense_minpoly(modulus, &f,
                                  &degree,
                                  n, matrix, minpoly)
        _sig_off
        v = []
        cdef Py_ssize_t i
        for i from 0 <= i <= degree:
            v.append(f[i])
        linbox_modn_dense_delete_array(f)
        return v


##     int  linbox_modn_dense_matrix_matrix_multiply(unsigned long modulus, mod_int **ans,
##                                                   mod_int **A, mod_int **B,
##                                                   size_t A_nr, size_t A_nc,
##                                                   size_t B_nr, size_t B_nc)


##     int linbox_modn_dense_rank(unsigned long modulus,
##                                mod_int** matrix, size_t nrows, size_t ncols)

