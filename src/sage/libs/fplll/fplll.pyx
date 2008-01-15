"""
Wrapper for the fpLLL library by Damien Stehle and David Cade found at

   http://perso.ens-lyon.fr/damien.stehle/english.html

This wrapper provides access to all fpLLL LLL implementations except
for integer matrices over C ints as opposed to integer matrices over
multi-precision integers. If that feature is required, please let the
authors of this know.

AUTHORS:
    -- 2007-10 Martin Albrecht <malb@informatik.uni-bremen.de>
       initial release
"""

#*****************************************************************************
#
#   SAGE: System for Algebra and Geometry Experimentation
#
#       Copyright (C) 2007 Martin Albrecht <malb@informatik.uni-bremen.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include "../../ext/interrupt.pxi"

from sage.matrix.matrix_integer_dense cimport Matrix_integer_dense
from sage.rings.integer_ring import ZZ
from sage.matrix.constructor import matrix

cdef class FP_LLL:
    """
    A basic wrapper class to support conversion to/from SAGE integer
    matrices and executing the LLL computation.

    NOTE: Usually you don't want to create this object yourself but
    use the LLL method of the integer matrices.
    """
    def __new__(self, Matrix_integer_dense A):
        """
        Construct a new fpLLL wrapper for the matrix A.

        INPUT:
            A -- a matrix over the integers

        EXAMPLE:
            sage: from sage.libs.fplll.fplll import FP_LLL
            sage: A = random_matrix(ZZ,10,10)
            sage: FP_LLL(A)
            fpLLL wrapper acting on a 10 x 10 matrix

        """
        cdef int i,j
        self._lattice = ZZ_mat_new(A._nrows,A._ncols)

        cdef Z_NR *t

        for i from 0 <= i < A._nrows:
            for j from 0 <= j < A._ncols:
                t = Z_NR_new()
                t.set_mpz_t(A._matrix[i][j])
                self._lattice.Set(i,j,t[0])
                Z_NR_delete(t)

    def __dealloc__(self):
        """
        Destroy internal data.
        """
        ZZ_mat_delete(self._lattice)

    def __repr__(self):
        return "fpLLL wrapper acting on a %d x %d matrix"%(self._lattice.GetNumRows(),self._lattice.GetNumCols())

    def _sage_(self):
        """
        Return a SAGE representation of self's matrix.

        EXAMPLE:
            sage: from sage.libs.fplll.fplll import FP_LLL
            sage: A = random_matrix(ZZ,10,10)
            sage: fpLLL = FP_LLL(A)
            sage: fpLLL._sage_() == A
            True
        """
        return to_sage(self._lattice)

    cdef int _check_precision(self, int precision):
        """
        Check whether the provided precision is within valid
        bounds. If not raise a TypeError.

        INPUT:
            precision -- an integer
        """
        if precision < 0:
            raise TypeError, "precision must be >= 0"

    cdef int _check_eta(self, float eta):
        """
        Check whether the provided parameter $\eta$ is within valid
        bounds. If not raise a TypeError.

        INPUT:
            eta -- a floating point number
        """
        if eta < 0.5:
            raise TypeError, "eta must be >= 0.5"

    cdef int _check_delta(self, float delta):
        """
        Check whether the provided parameter $\delta$ is within valid
        bounds. If not raise a TypeError.

        INPUT:
            delta -- a floating point number
        """
        if delta <= 0.25:
            raise TypeError, "delta must be > 0.25"
        if delta > 1.0:
            raise TypeError, "delta must be <= 1.0"

    def wrapper(self, int precision=0, float eta=0.51, float delta=0.99):
        """
        Perform LLL reduction using fpLLL's \code{wrapper}
        implementation. This implementation invokes a sequence of
        floating point LLL computations such that
          * the computation is reasonably fast (based on an heuristic modell)
          * the result is proven to be LLL reduced.

        INPUT:
            precision -- internal precision (default: auto)
            eta -- LLL parameter eta with 1/2 <= eta < sqrt(delta) (default: 0.51)
            delta -- LLL parameter delta with 1/4 < delta <= 1 (default: 0.99)

        OUTPUT:
            nothing is returned but the internal state modified.

        EXAMPLE:
            sage: from sage.libs.fplll.fplll import FP_LLL
            sage: W = random_matrix(ZZ,10,10)
            sage: W # result random
            [ -3  -1  -1  -2  -1   7  -1   1  -1  -2]
            [  1   4  -1  -1  -1  -1  -1   3   1   1]
            [ -4  -1   1   1  -3  -2  -5   1   4  -1]
            [  1   1  -1  -1  -1  -2   3  -2  -1   1]
            [  4  -1 -12   1  -1   2   1  -2  -2   1]
            [ -2   2  -1  -2   6  -5   5   1  -1  -1]
            [ -1  -1  -1  -2  -2  -4   9   1   3   1]
            [  2  -2   1  -1  -1  -1   1  -1 -22  -1]
            [ -1   1   1  -1  -1 -49  20   2 -16   2]
            [  1  30   1  -1   1  -1  -8  -3  -1  -1]
            sage: F = FP_LLL(W)
            sage: F.wrapper()
            sage: W = F._sage_()
            sage: W # result depending on random input
            [ 1  1 -1 -1 -1 -2  3 -2 -1  1]
            [ 1  4 -1 -1 -1 -1 -1  3  1  1]
            [-3  0  0  0 -4 -4 -2 -1  3  0]
            [ 4 -2 -2  2 -4 -4 -2 -1  0 -2]
            [-2  0 -2 -3 -2  5  2 -1 -2 -1]
            [-1 -1 -4 -2  1 -2  2  1 -2 -5]
            [-2 -2  0 -1 -1 -2  6  3  4  0]
            [ 0  2 -4 -2  6  3  1 -4  2 -3]
            [ 1 -5 -3  5  0 -1  1  1 -3  6]
            [-5  5 -1  6 -2  2  2 -3  1 -3]
        """
        self._check_precision(precision)
        self._check_eta(eta)
        self._check_delta(delta)
        cdef int ret = 0

        cdef wrapper *w = wrapper_new(self._lattice, 0, eta, delta)
        _sig_on
        ret = w.LLL()
        _sig_off
        wrapper_delete(w)
        if ret < 0:
            raise RuntimeError, "fpLLL returned %d < 0"%ret


    def proved(self, int precision=0, float eta=0.51, float delta=0.99, implementation=None):
        """
        Perform LLL reduction using fpLLL's \code{proved}
        implementation. This implementation is the only provable
        correct floating point implementation in the free software
        world. Provability is only guaranteed if the 'mpfr'
        implementation is chosen.

        INPUT:
            precision -- internal precision (default: auto)
            eta -- LLL parameter eta with 1/2 <= eta < sqrt(delta) (default: 0.51)
            delta -- LLL parameter delta with 1/4 < delta <= 1 (default: 0.99)
            implementation -- floating point implementation in ('double', 'dpe', 'mpfr')
                              (default: 'mpfr')

        OUTPUT:
            nothing is returned but the internal state modified.

        EXAMPLE:
            sage: from sage.libs.fplll.fplll import FP_LLL
            sage: P = random_matrix(ZZ,10,10)
            sage: P # result random
            [   10     1     1    22     1    -1     3    -1 -3974    -2]
            [   -3     3    -3    -1     4    -1     4    -1     1     1]
            [    1    10     2     2    -1    -3    -2    -2   -12     2]
            [   -7    -1   -14    -3    17    -5    -9    -1   -82     1]
            [  -11     1     5     1    -1    -1     1     3    -4     1]
            [    8    -1    -2    -2     1    -6     3    -1     2    -1]
            [    1    -1    -1    -3    22     3     1     1     1    -2]
            [    1    -2    -1    32    -1    -1     2     1     1    -1]
            [    2    -5    -2    -1     2    -1    -7    -2   -10    -1]
            [   -4    -3    -2 -1188     2     1    -1    50     1     1]
            sage: F = FP_LLL(P)
            sage: F.proved()
            sage: P = F._sage_()
            sage: P # result depending on random input
            [  -3    3   -3   -1    4   -1    4   -1    1    1]
            [  -3    0    3   -1    0   -7    4    2   -2    0]
            [  -8    1    2    2   -1    6   -3    1   -2    1]
            [  -1   -2   -5   -2    6   -2   -3   -3   -9    0]
            [   2   12    7    4   -7   -1    1    1   -3    2]
            [   2   -2   10   -1   12   -1    4    7    7   -3]
            [   1   -6   -4    0   -2   10   17   13   -7   -2]
            [  -2    1   -4   31    3   -2    6    0    2    0]
            [   0   17  -15  -11    0   -4  -14   37    6  -13]
            [ -18   31    5   -4  -14    9   12  -33   -3 -152]



        """
        self._check_precision(precision)
        self._check_eta(eta)
        self._check_delta(delta)

        cdef proved_double *pdouble
        cdef proved_mpfr *pmpfr
        cdef proved_dpe *pdpe
        cdef int ret = 0

        if implementation is None:
            implementation = "mpfr"

        if implementation == "double":
           pdouble = proved_double_new(self._lattice, precision, eta, delta)
           _sig_on
           ret = pdouble.LLL()
           _sig_off
           proved_double_delete(pdouble)
        elif implementation == "dpe":
           pdpe = proved_dpe_new(self._lattice, precision, <double>eta, <double>delta)
           _sig_on
           ret = pdpe.LLL()
           _sig_off
           proved_dpe_delete(pdpe)
        elif implementation == "mpfr":
           pmpfr = proved_mpfr_new(self._lattice, precision, eta, delta)
           _sig_on
           ret = pmpfr.LLL()
           _sig_off
           proved_mpfr_delete(pmpfr)

        if ret < 0:
            raise RuntimeError, "fpLLL returned %d < 0"%ret

    def fast(self, int precision=0, float eta=0.51, float delta=0.99, implementation=None):
        """
        Perform LLL reduction using fpLLL's \code{fast}
        implementation. This implementation is the fastest floating
        point implementation available in the free software world.

        INPUT:
            precision -- internal precision (default: auto)
            eta -- LLL parameter eta with 1/2 <= eta < sqrt(delta) (default: 0.51)
            delta -- LLL parameter delta with 1/4 < delta <= 1 (default: 0.99)
            implementation -- ignored

        OUTPUT:
            nothing is returned but the internal state modified.

        EXAMPLE:
            sage: from sage.libs.fplll.fplll import FP_LLL
            sage: F = random_matrix(ZZ,10,10)
            sage: F # result random
            [ -1 -20   1  12   2   2   4   1   1  -1]
            [  1  -2  -9   1   2   3  -1  -4 -60  -1]
            [ -4 -15  -8   1   3   2 -60  -1  -2   6]
            [ -5  -7  -1  -2   1   1   1   6   3   1]
            [  5  -6  -2  -1  -2   1   2  -4  -3  -1]
            [  3   1   3  -1  10  -2   1  -1  -2  -1]
            [  4   2 -15   3   1  -9 -12   1  -1  -1]
            [  1   4  16  -1   3   1  -1  -2  -1 -75]
            [-56   1  -1  -1   1   1  -4   3   1   2]
            [  2  -9  -4  -1  -1  -3   4   7  -2   3]
            sage: f = FP_LLL(F)
            sage: f.fast()
            sage: F = f._sage_()
            sage: F # result depending on random input
            [  2   4  -1   2   0  -5   1   5  -2   3]
            [  5  -6  -2  -1  -2   1   2  -4  -3  -1]
            [ -5  -7  -1  -2   1   1   1   6   3   1]
            [  3   1   3  -1  10  -2   1  -1  -2  -1]
            [ -1  -7   4  15   3   0   1  -1   1  -1]
            [  2  -2 -14   1   1  -4 -13  -4   1  -4]
            [-24  -7  -7   0   2  -8   6 -10 -20   2]
            [  3   1   0  -2  -5  22 -13  11 -25  -5]
            [  7  -6  12 -12  -2 -16 -18 -12   9  20]
            [  8  -7  17  -9  -9 -30 -12  15   5 -41]
        """
        self._check_precision(precision)
        self._check_eta(eta)
        self._check_delta(delta)

        cdef fast_double *pdouble
        cdef int ret = 0

        pdouble = fast_double_new(self._lattice, precision, eta, delta)
        _sig_on
        ret = pdouble.LLL()
        _sig_off
        fast_double_delete(pdouble)

        if ret < 0:
            raise RuntimeError, "fpLLL returned %d < 0"%ret

    def fast_early_red(self, int precision=0, float eta=0.51, float delta=0.99, implementation=None):
        """
        Perform LLL reduction using fpLLL's \code{fast}
        implementation with early reduction.

        This implementation inserts some early reduction steps inside
        the execution of the 'fast' LLL algorithm. This sometimes makes the
        entries of the basis smaller very quickly. It occurs in
        particular for lattice bases built from minimal polynomial or
        integer relation detection problems.

        INPUT:
            precision -- internal precision (default: auto)
            eta -- LLL parameter eta with 1/2 <= eta < sqrt(delta) (default: 0.51)
            delta -- LLL parameter delta with 1/4 < delta <= 1 (default: 0.99)
            implementation -- ignored

        OUTPUT:
            nothing is returned but the internal state modified.

        EXAMPLE:
            sage: from sage.libs.fplll.fplll import FP_LLL
            sage: FE = random_matrix(ZZ,10,10)
            sage: FE # result random
            [  -4    3    2   -1    1   -1    6    2   -4  -15]
            [   1   -1    6    4   -1    4   -6   -1    2    2]
            [   1   -1   -1   -7    1    1   -3   -2    1    5]
            [ 232   -2   -1   -1   -4   -1   -6    1    1    1]
            [  -6    4    1   -1   -1   -1    4    3   -2   -2]
            [  -2   -1   -2   -2    6   -2   -2    1    4   -1]
            [  20   -1  -61    1   -2    1    1    2    1    3]
            [   9   -2   -2   -9   -2  -37   11    1   -1    1]
            [  -1   -8    1 -108    1    8    1   -1  -35   -5]
            [  -1   -1   48    2    1    1    8    4    7   -4]

            sage: F = FP_LLL(FE)
            sage: F.fast_early_red()
            sage: FE = F._sage_()
            sage: FE # result depending on random input
            [ -2  -1  -2  -2   6  -2  -2   1   4  -1]
            [  1  -1  -1  -7   1   1  -3  -2   1   5]
            [ -6   4   1  -1  -1  -1   4   3  -2  -2]
            [ -5   3   7   3  -2   3  -2   2   0   0]
            [  3  -2   0  -7   3   1  -1  -3  -1  -8]
            [  2   6  -1  -1   1   6   4  12  11   1]
            [  2  12  -9  -4 -12  -8 -13   8   1  -5]
            [  6 -17   4  -6   0  -2  -3  16 -13   1]
            [ -6   1 -17  -2   6  26   3  -1  -3  -5]
            [ 32  45  11  10  37  -8 -10  17 -35  12]

         """
        self._check_precision(precision)
        self._check_eta(eta)
        self._check_delta(delta)
        cdef int ret = 0

        cdef fast_early_red_double *pdouble

        pdouble = fast_early_red_double_new(self._lattice, precision, eta, delta)
        _sig_on
        ret = pdouble.LLL()
        _sig_off
        fast_early_red_double_delete(pdouble)

        if ret < 0:
            raise RuntimeError, "fpLLL returned %d < 0"%ret

    def heuristic(self, int precision=0, float eta=0.51, float delta=0.99, implementation=None):
        """
        Perform LLL reduction using fpLLL's \code{heuristic}
        implementation.

        INPUT:
            precision -- internal precision (default: auto)
            eta -- LLL parameter eta with 1/2 <= eta < sqrt(delta) (default: 0.51)
            delta -- LLL parameter delta with 1/4 < delta <= 1 (default: 0.99)
            implementation -- floating point implementation in ('double', 'dpe', 'mpfr')
                              (default: 'mpfr')

        OUTPUT:
            nothing is returned but the internal state modified.

        EXAMPLE:
            sage: from sage.libs.fplll.fplll import FP_LLL
            sage: H = random_matrix(ZZ,10,10)
            sage: H # result random
            [  2   1  -1  -6  -3  -8  -4  -4   1   3]
            [ -2   3  -1  -1  -2   1   1  -1   5 -11]
            [ -1  -1   1  -1  -2  -3   1   2  -7  -1]
            [-22  -1  -1  13   5  -3   1   1   2   1]
            [  4  -1  -6   3  -1   2  59   1   1  -1]
            [ -2  -1  -2  -1  -1  20   6   1  -2  -1]
            [ -1   2  -1  -1  -2   1  -3  -2   2  -1]
            [ -1  -3  -1 -10   1  -1  -1 -35   1  -1]
            [  5  -2  -1   1 -44  -1  -4  -1  -1   3]
            [ -1  -1  -4  -1  -1   3  -1 -15   1   1]
            sage: F = FP_LLL(H)
            sage: F.heuristic()
            sage: H = F._sage_()
            sage: H # result depending on random input
            [ -1   2  -1  -1  -2   1  -3  -2   2  -1]
            [ -2   1   0  -2  -4  -2  -2   0  -5  -2]
            [  3  -1   0  -5  -1  -9  -1  -2  -1   4]
            [ -1   1   0   0   0   0   4   1   3 -10]
            [ -1  -1   7  -3   4   2  -2  -4  -3   3]
            [  0  -3  -3   0   1   2   2 -13  -1   2]
            [  1  -2  -2  -6  -2  11   5  -1  -3   3]
            [-15  -3  -4  -2   3  -7   7  -5   7   3]
            [  6  19  13   8  -8   2   9   1  -3   2]
            [  4  -6  17  12 -21   3   2  -1   6  -2]
        """
        self._check_precision(precision)
        self._check_eta(eta)
        self._check_delta(delta)

        cdef heuristic_double *pdouble
        cdef heuristic_mpfr *pmpfr
        cdef heuristic_dpe *pdpe
        cdef int ret = 0

        if implementation is None:
            implementation = "mpfr"

        if implementation == "double":
            pdouble = heuristic_double_new(self._lattice, precision, eta, delta)
            _sig_on
            ret = pdouble.LLL()
            _sig_off
            heuristic_double_delete(pdouble)
        elif implementation == "dpe":
            pdpe = heuristic_dpe_new(self._lattice, precision, <double>eta, <double>delta)
            _sig_on
            ret = pdpe.LLL()
            _sig_off
            heuristic_dpe_delete(pdpe)
        elif implementation == "mpfr":
            pmpfr = heuristic_mpfr_new(self._lattice, precision, eta, delta)
            _sig_on
            ret = pmpfr.LLL()
            _sig_off
            heuristic_mpfr_delete(pmpfr)

        if ret < 0:
            raise RuntimeError, "fpLLL returned %d < 0"%ret


    def heuristic_early_red(self, int precision=0, float eta=0.51, float delta=0.99, implementation=None):
        """
        Perform LLL reduction using fpLLL's \code{heuristic}
        implementation with early reduction.

        This implementation inserts some early reduction steps inside
        the execution of the 'fast' LLL algorithm. This sometimes
        makes the entries of the basis smaller very quickly. It occurs
        in particular for lattice bases built from minimal polynomial
        or integer relation detection problems.

        INPUT:
            precision -- internal precision (default: auto)
            eta -- LLL parameter eta with 1/2 <= eta < sqrt(delta) (default: 0.51)
            delta -- LLL parameter delta with 1/4 < delta <= 1 (default: 0.99)
            implementation -- floating point implementation in ('double', 'dpe', 'mpfr')
                              (default: 'mpfr')

        OUTPUT:
            nothing is returned but the internal state modified.

        EXAMPLE:
            sage: from sage.libs.fplll.fplll import FP_LLL
            sage: HE = random_matrix(ZZ,10,10)
            sage: HE # result random
            [ 11   1  -1   3  -4  -1   1   7   8   2]
            [  5  -1  -1   3  -1  37   1  -4   6   1]
            [ -1   3  -1  -4   1  -1  -1   1   3   1]
            [  3  -1   1   1   1  -2   1   6   1  -5]
            [ -2 -12  10  -3  -1  10   1   1   1  -6]
            [  1   1  -1   1  -3   1  -1  -5   1   2]
            [ -3  11   1  -2  -2   3   6  -1  -1  -9]
            [ 11  10   1  -1   1  -1   3   2   1   7]
            [ -1   1  30  -1   2   1  -1   2  -2  -1]
            [ -5   3   3  -2   1   1   1 -77  -1  -1]

            sage: F = FP_LLL(HE)
            sage: F.heuristic_early_red()
            sage: HE = F._sage_()
            sage: HE # result depending on random input
            [  4   0   0   2  -2  -1   0   1   2  -3]
            [ -1   3  -1  -4   1  -1  -1   1   3   1]
            [ -1  -1   1  -1   3  -1   1   5  -1  -2]
            [  4  -2   0   3  -1   2   2   4   1   7]
            [  4   9   2  -2   3  -1   2  -4  -5   2]
            [ -7   0  -1   1  -4   7   6   6   3  -1]
            [  5   0  11  -7   1   7   2  -1   1  -6]
            [ -2   5  17   5  -3  -7  -5   0   3   5]
            [  1   6   8   6   8   9 -12   4   6  -2]
            [  0  -2  -4   3  12  -8  16 -11  13   0]
        """
        self._check_precision(precision)
        self._check_eta(eta)
        self._check_delta(delta)

        cdef heuristic_early_red_double *pdouble
        cdef heuristic_early_red_mpfr *pmpfr
        cdef heuristic_early_red_dpe *pdpe
        cdef int ret = 0

        if implementation is None:
            implementation = "mpfr"

        if implementation == "double":
           pdouble = heuristic_early_red_double_new(self._lattice, precision, eta, delta)
           _sig_on
           ret = pdouble.LLL()
           _sig_off
           heuristic_early_red_double_delete(pdouble)
        elif implementation == "dpe":
           pdpe = heuristic_early_red_dpe_new(self._lattice, precision, <double>eta, <double>delta)
           _sig_on
           ret = pdpe.LLL()
           _sig_off
           heuristic_early_red_dpe_delete(pdpe)
        elif implementation == "mpfr":
           pmpfr = heuristic_early_red_mpfr_new(self._lattice, precision, eta, delta)
           _sig_on
           ret = pmpfr.LLL()
           _sig_off
           heuristic_early_red_mpfr_delete(pmpfr)

        if ret < 0:
            raise RuntimeError, "fpLLL returned %d < 0"%ret

def gen_intrel(int d, int b):
    """
    Return a (d+1 x d)-dimensional knapsack-type random lattice,
    where the x_i's are random b bits integers.

    INPUT:
        d -- dimension
        b -- bitsize of entries
    """
    cdef ZZ_mat *A = ZZ_mat_new(d,d+1)
    A.gen_intrel(b)

    B = to_sage(A)
    ZZ_mat_delete(A)
    return B

def gen_simdioph(int d, int b, int b2):
    """
    Return a d-dimensional simultaneous diophantine approximation
    random lattice, where the d $x_i$'s are random b bits integers.

    INPUT:
        d -- dimension
        b -- bitsize of entries
        b2 -- bitsize of entries
    """
    cdef ZZ_mat *A = ZZ_mat_new(d,d)
    A.gen_simdioph(b, b2)

    B = to_sage(A)
    ZZ_mat_delete(A)
    return B

def gen_uniform(int nr, int nc, int b):
    """
    Return a (nr x nc) matrix where the entries are random b bits
    integers.

    INPUT:
        nr -- row dimension
        nc -- column dimension
        b-- bitsize of entries
    """
    cdef ZZ_mat *A = ZZ_mat_new(nr,nc)
    A.gen_uniform(b)

    B = to_sage(A)
    ZZ_mat_delete(A)
    return B

def gen_ntrulike(int d, int b, int q):
    """
    Generate a ntru-like lattice of dimension (2*d x 2*d), with the
    coefficients $h_i$ chosen as random b bits integers and parameter
    q:

    \begin{verbatim}
      [[ 1 0 ... 0 h0      h1 ... h_{d-1} ]
       [ 0 1 ... 0 h1      h2 ... h0      ]
       [ ................................ ]
       [ 0 0 ... 1 h_{d-1} h0 ... h_{d-1} ]
       [ 0 0 ... 0 q       0  ...  0      ]
       [ 0 0 ... 0 0       q  ...  0      ]
       [ ................................ ]
       [ 0 0 ... 0 0       0  ...  q      ]]
    \end{verbatim}

    INPUT:
        d -- dimension
        b -- bitsize of entries
        q -- see above
    """
    cdef ZZ_mat *A = ZZ_mat_new(2*d,2*d)
    A.gen_ntrulike(b, q)

    B = to_sage(A)
    ZZ_mat_delete(A)
    return B

def gen_ntrulike2(int d, int b, int q):
    """
    Like gen_ntrulike but with the q vectors coming first.

    INPUT:
        d -- dimension
        b -- bitsize of entries
        q -- see gen_ntrulike
    """
    cdef ZZ_mat *A = ZZ_mat_new(2*d,2*d)
    A.gen_ntrulike2(b,q)

    B = to_sage(A)
    ZZ_mat_delete(A)
    return B

def gen_ajtai(int d, float alpha):
    """
    Return Ajtai-like (d x d)-matrix of fp parameter alpha.  The
    matrix is lower-triangular, Bii is ~2^((d-i+1)^alpha) and Bij is
    ~Bjj/2 for j<i.

    INPUT:
        d -- dimension
        alpha -- see above
    """
    cdef ZZ_mat *A = ZZ_mat_new(d,d)
    A.gen_ajtai(alpha)

    B = to_sage(A)
    ZZ_mat_delete(A)
    return B


cdef to_sage(ZZ_mat *A):
    """
    Return a SAGE integer matrix for A. A is not destroyed.

    INPUT:
        A -- ZZ_mat
    """
    cdef int i,j
    cdef mpz_t *t

    cdef Matrix_integer_dense B = <Matrix_integer_dense>matrix(ZZ,
                                                               A.GetNumRows(),
                                                               A.GetNumCols())

    for i from 0 <= i < A.GetNumRows():
        for j from 0 <= j < A.GetNumCols():
            t = &B._matrix[i][j]
            mpz_set(t[0], A.Get(i,j).GetData())
    return B
