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

        """
        cdef int i,j
        self._lattice = ZZ_mat_new(A._nrows,A._ncols)

        cdef Z_NR *t

        for i from 0 <= i < A._nrows:
            for j from 0 <= j < A._ncols:
                t = Z_NR_new()
                t.set_mpz_t(A._matrix[i][j])
                self._lattice.Set(i,j,t[0])

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
            sage: A = random_matrix(ZZ,10,10); A # result random
            sage: f = FP_LLL(A)
            sage: f.wrapper()
            sage: B = A._sage_(); B # result random
        """
        self._check_precision(precision)
        self._check_eta(eta)
        self._check_delta(delta)

        cdef wrapper *w = wrapper_new(self._lattice, 0, eta, delta)
        _sig_on
        w.LLL()
        _sig_off
        wrapper_delete(w)

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
            sage: A = random_matrix(ZZ,10,10); A # result random
            sage: f = FP_LLL(A)
            sage: f.proved()
            sage: B = A._sage_(); B # result random
        """
        self._check_precision(precision)
        self._check_eta(eta)
        self._check_delta(delta)

        cdef proved_double *pdouble
        cdef proved_mpfr *pmpfr
        cdef proved_dpe *pdpe

        if implementation is None:
            implementation = "mpfr"

        if implementation == "double":
           pdouble = proved_double_new(self._lattice, precision, eta, delta)
           _sig_on
           pdouble.LLL()
           _sig_off
           proved_double_delete(pdouble)
        elif implementation == "dpe":
           pdpe = proved_dpe_new(self._lattice, precision, <double>eta, <double>delta)
           _sig_on
           pdpe.LLL()
           _sig_off
           proved_dpe_delete(pdpe)
        elif implementation == "mpfr":
           pmpfr = proved_mpfr_new(self._lattice, precision, eta, delta)
           _sig_on
           pmpfr.LLL()
           _sig_off
           proved_mpfr_delete(pmpfr)

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
            sage: A = random_matrix(ZZ,10,10); A # result random
            sage: f = FP_LLL(A)
            sage: f.fast()
            sage: B = A._sage_(); B # result random
        """
        self._check_precision(precision)
        self._check_eta(eta)
        self._check_delta(delta)

        cdef fast_double *pdouble

        pdouble = fast_double_new(self._lattice, precision, eta, delta)
        _sig_on
        pdouble.LLL()
        _sig_off
        fast_double_delete(pdouble)

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
            sage: A = random_matrix(ZZ,10,10); A # result random
            sage: f = FP_LLL(A)
            sage: f.fast_early_red()
            sage: B = A._sage_(); B # result random
        """
        self._check_precision(precision)
        self._check_eta(eta)
        self._check_delta(delta)

        cdef fast_early_red_double *pdouble

        pdouble = fast_early_red_double_new(self._lattice, precision, eta, delta)
        _sig_on
        pdouble.LLL()
        _sig_off
        fast_early_red_double_delete(pdouble)

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
            sage: A = random_matrix(ZZ,10,10); A # result random
            sage: f = FP_LLL(A)
            sage: f.heuristic()
            sage: B = A._sage_(); B # result random
        """
        self._check_precision(precision)
        self._check_eta(eta)
        self._check_delta(delta)

        cdef heuristic_double *pdouble
        cdef heuristic_mpfr *pmpfr
        cdef heuristic_dpe *pdpe

        if implementation is None:
            implementation = "mpfr"

        if implementation == "double":
           pdouble = heuristic_double_new(self._lattice, precision, eta, delta)
           _sig_on
           pdouble.LLL()
           _sig_off
           heuristic_double_delete(pdouble)
        elif implementation == "dpe":
           pdpe = heuristic_dpe_new(self._lattice, precision, <double>eta, <double>delta)
           _sig_on
           pdpe.LLL()
           _sig_off
           heuristic_dpe_delete(pdpe)
        elif implementation == "mpfr":
           pmpfr = heuristic_mpfr_new(self._lattice, precision, eta, delta)
           _sig_on
           pmpfr.LLL()
           _sig_off
           heuristic_mpfr_delete(pmpfr)

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
            sage: A = random_matrix(ZZ,10,10); A # result random
            sage: f = FP_LLL(A)
            sage: f.heuristic_early_red()
            sage: B = A._sage_(); B # result random
        """
        self._check_precision(precision)
        self._check_eta(eta)
        self._check_delta(delta)

        cdef heuristic_early_red_double *pdouble
        cdef heuristic_early_red_mpfr *pmpfr
        cdef heuristic_early_red_dpe *pdpe

        if implementation is None:
            implementation = "mpfr"

        if implementation == "double":
           pdouble = heuristic_early_red_double_new(self._lattice, precision, eta, delta)
           _sig_on
           pdouble.LLL()
           _sig_off
           heuristic_early_red_double_delete(pdouble)
        elif implementation == "dpe":
           pdpe = heuristic_early_red_dpe_new(self._lattice, precision, <double>eta, <double>delta)
           _sig_on
           pdpe.LLL()
           _sig_off
           heuristic_early_red_dpe_delete(pdpe)
        elif implementation == "mpfr":
           pmpfr = heuristic_early_red_mpfr_new(self._lattice, precision, eta, delta)
           _sig_on
           pmpfr.LLL()
           _sig_off
           heuristic_early_red_mpfr_delete(pmpfr)


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
