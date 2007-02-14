## NOTE: The _sig_on/_sig_off stuff can't go in here -- it has to be in the
## code that calls these functions.  Otherwise strangely objects get left
## in an incorrect state.

from sage.rings.integer cimport Integer
from sage.misc.misc import verbose, get_verbose, cputime, UNAME

##########################################################################
## Dense matrices modulo p
##########################################################################

# LinBox bugs to address:
#  * echelon form over GF(2) -> crash, worked around by using native 'gauss' in that case
#  * charpoly and minpoly don't work randomly

cdef extern from "linbox_wrap.h":
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

    def minpoly(self):
        return self._poly(True)

    def charpoly(self):
        return self._poly(False)

    def _poly(self, minpoly):
        """
        INPUT:
            as given

        OUTPUT:
            coefficients of charpoly or minpoly as a Python list
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

cdef extern from "linbox_wrap.h":
    void linbox_integer_dense_minpoly_hacked(mpz_t* *minpoly, size_t* degree,
                                      size_t n, mpz_t** matrix, int do_minpoly)

    void linbox_integer_dense_minpoly(mpz_t* *minpoly, size_t* degree,
                                      size_t n, mpz_t** matrix)

    void linbox_integer_dense_charpoly(mpz_t* *charpoly, size_t* degree,
                                       size_t n, mpz_t** matrix)

    void linbox_integer_dense_delete_array(mpz_t* f)

    int linbox_integer_dense_matrix_matrix_multiply(mpz_t** ans, mpz_t **A, mpz_t **B,
                                      size_t A_nr, size_t A_nc, size_t B_nr, size_t B_nc)

    unsigned long linbox_integer_dense_rank(mpz_t** matrix, size_t nrows,
                                            size_t ncols)

    unsigned long linbox_integer_dense_det(mpz_t** matrix, size_t nrows,
                                           size_t ncols)

cdef class Linbox_integer_dense:
    cdef set(self, mpz_t** matrix, size_t nrows, size_t ncols):
        self.nrows = nrows
        self.ncols = ncols
        self.matrix = matrix

    def minpoly(self):
        return self._poly(True)

    def charpoly(self):
        return self._poly(False)

    def _poly(self, do_minpoly):
        """
        INPUT:
            as given

        OUTPUT:
            coefficients of charpoly or minpoly as a Python list
        """
        cdef mpz_t* poly
        cdef size_t degree
        if self.nrows % 4 == 0 and UNAME == "Darwin":
            verbose("using hack to get around bug in linbox on OS X since n is divisible by 4")
            if do_minpoly:
                linbox_integer_dense_minpoly_hacked(&poly, &degree, self.nrows, self.matrix, 1)
            else:
                linbox_integer_dense_minpoly_hacked(&poly, &degree, self.nrows, self.matrix, 0)
        else:
            verbose("using linbox poly comp")
            if do_minpoly:
                linbox_integer_dense_minpoly(&poly, &degree, self.nrows, self.matrix)
            else:
                linbox_integer_dense_charpoly(&poly, &degree, self.nrows, self.matrix)
            verbose("computed poly -- now converting back to SAGE")

        v = []
        cdef Integer k
        cdef size_t n
        for n from 0 <= n <= degree:
            k = Integer()
            mpz_set(k.value, poly[n])
            mpz_clear(poly[n])
            v.append(k)
        linbox_integer_dense_delete_array(poly)
        return v

    cdef matrix_matrix_multiply(self,
                                mpz_t **ans,
                                mpz_t **B,
                                size_t B_nr, size_t B_nc):
        cdef int e
        e = linbox_integer_dense_matrix_matrix_multiply(ans,
                                                        self.matrix,  B,
                                                        self.nrows, self.ncols,
                                                        B_nr, B_nc)
        if e:
            raise RuntimError, "error doing matrix matrix multiply over ZZ using linbox"


    cdef unsigned long rank(self) except -1:
        return linbox_integer_dense_rank(self.matrix, self.nrows, self.ncols)

    def det(self):
        linbox_integer_dense_det(self.matrix, self.nrows, self.ncols)
