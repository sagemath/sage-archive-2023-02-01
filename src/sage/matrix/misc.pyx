"""
Misc matrix algorithms

Code goes here mainly when it needs access to the internal structure
of several classes, and we want to avoid circular cimports.
"""

include "../ext/interrupt.pxi"
include "../ext/cdefs.pxi"
include '../ext/stdsage.pxi'
include "../ext/gmp.pxi"   # rational reconstruction

# mod_int isn't defined in stdio.h -- this is to fool SageX.
cdef extern from "stdio.h":
    ctypedef int mod_int

from matrix0 cimport Matrix
from matrix_modn_dense cimport Matrix_modn_dense
from matrix_integer_dense cimport Matrix_integer_dense
from matrix_rational_dense cimport Matrix_rational_dense
import matrix_modn_dense

from sage.rings.integer_ring   import ZZ
from sage.rings.rational_field import QQ
from sage.rings.integer cimport Integer
from sage.rings.arith import previous_prime, CRT_basis


from sage.misc.misc import verbose, get_verbose


def matrix_modn_dense_lift(Matrix_modn_dense A):
    cdef Py_ssize_t i, j
    cdef Matrix_integer_dense L
    L = Matrix_integer_dense.__new__(Matrix_integer_dense,
                                     A.parent().change_ring(ZZ),
                                     0, 0, 0)
    cdef mpz_t* L_row
    cdef mod_int* A_row
    for i from 0 <= i < A._nrows:
        L_row = L._matrix[i]
        A_row = A._matrix[i]
        for j from 0 <= j < A._ncols:
            mpz_init_set_si(L_row[j], A_row[j])
    L._initialized = 1
    return L



def matrix_integer_dense_rational_reconstruction(Matrix_integer_dense A, Integer N):
    cdef Matrix_rational_dense R
    R = Matrix_rational_dense.__new__(Matrix_rational_dense,
                                      A.parent().change_ring(QQ), 0,0,0)

    cdef mpz_t denom   # lcm of denoms so far
    for i from 0 <= i < A._nrows:
        for j from 0 <= j < A._ncols:
            mpq_rational_reconstruction(R._matrix[i][j], A._matrix[i][j], N.value)
    return R

def matrix_rational_echelon_form_multimodular(Matrix self, height_guess=None, proof=True):
    """
    Returns reduced row-echelon form using a multi-modular
    algorithm.  Does not change self.

    REFERENCE: Chapter 7 of Stein's "Explicitly Computing Modular Forms".

    INPUT:
        height_guess -- integer or None
        proof -- boolean (default: True)

    ALGORITHM:
    The following is a modular algorithm for computing the echelon
    form.  Define the height of a matrix to be the max of the
    absolute values of the entries.

    Given Matrix A with n columns (self).

     0. Rescale input matrix A to have integer entries.  This does
        not change echelon form and makes reduction modulo lots of
        primes significantly easier if there were denominators.
        Henceforth we assume A has integer entries.

     1. Let c be a guess for the height of the echelon form.  E.g.,
        c=1000, e.g., if matrix is very sparse and application is to
        computing modular symbols.

     2. Let M = n * c * H(A) + 1,
        where n is the number of columns of A.

     3. List primes p_1, p_2, ..., such that the product of
        the p_i is at least M.

     4. Try to compute the rational reconstruction CRT echelon form
        of A mod the product of the p_i.  If rational
        reconstruction fails, compute 1 more echelon forms mod the
        next prime, and attempt again.  Make sure to keep the
        result of CRT on the primes from before, so we don't have
        to do that computation again.  Let E be this matrix.

     5. Compute the denominator d of E.
        Attempt to prove that result is correct by checking that

              H(d*E)*ncols(A)*H(A) < (prod of reduction primes)

        where H denotes the height.   If this fails, do step 4 with
        a few more primes.
    """

    verbose("Multimodular echelon algorithm on %s x %s matrix"%(self._nrows, self._ncols))
    cdef Matrix_rational_dense E
    if self._nrows == 0 or self._ncols == 0:
        self.cache('in_echelon_form', True)
        self.cache('echelon_form', self)
        self.cache('pivots', [])
        return self

    B, _ = self._clear_denom()

    height = self.height()
    if height_guess is None:
        height_guess = 10000000*(height+100)
    tm = verbose("height_guess = %s"%height_guess, level=2)

    if proof:
        M = self._ncols * height_guess * height  +  1
    else:
        M = height_guess + 1

    p = matrix_modn_dense.MAX_MODULUS + 1
    X = []
    best_pivots = []
    prod = 1
    problem = 0
    while True:
        p = previous_prime(p)
        while prod < M:
            problem = problem + 1
            if problem > 50:
                verbose("sparse_matrix multi-modular reduce not converging?")
            t = verbose("echelon modulo p=%s (%.2f%% done)"%(
                       p, 100*float(len(str(prod))) / len(str(M))), level=2)

            # We use denoms=False, since we made self integral by calling clear_denom above.
            A = B._mod_int(p)
            t = verbose("time to reduce matrix mod p:",t, level=2)
            A.echelonize()
            t = verbose("time to put reduced matrix in echelon form:",t, level=2)

            # a very worthwhile check (we will extend the algorithm to
            # switch to a system solving method in the nonsquare full
            # rank case at some point).
            if self._nrows == self._ncols and len(A.pivots()) == self._nrows:
                verbose("done: the echelon form mod p is the identity matrix")
                E = self.parent().identity_matrix()
                E.cache('pivots', range(self._nrows))
                E.cache('in_echelon_form', True)
                self.cache('in_echelon_form', True)
                self.cache('echelon_form', E)
                self.cache('pivots', range(self._nrows))
                return E

            c = cmp_pivots(best_pivots, A.pivots())
            if c <= 0:
                best_pivots = A.pivots()
                X.append(A)
                prod = prod * p
            else:
                # do not save A since it is bad.
                if LEVEL > 1:
                    verbose("Excluding this prime (bad pivots).")
            t = verbose("time for pivot compare", t, level=2)
            p = previous_prime(p)
        # Find set of best matrices.
        Y = []
        # recompute product, since may drop bad matrices
        prod = 1
        for i in range(len(X)):
            if cmp_pivots(best_pivots, X[i].pivots()) <= 0:
                Y.append((matrix_modn_dense_lift(X[i]), X[i].base_ring().order()))
                prod = prod * X[i].base_ring().order()
        try:
            t = verbose("start crt linear combination", level=2)
            a = CRT_basis([w[1] for w in Y])
            # take the linear combination of the lifts of the elements
            # of Y times coefficients in a
            L = a[0]*(Y[0][0])
            for j in range(1,len(Y)):
                L += a[j]*(Y[j][0])
            t = verbose("crt time is",t, level=2)
            E = L.rational_reconstruction(prod)
            L =0  # free memory
            verbose('rational reconstruction time is', t, level=2)
        except ValueError, msg:
            verbose("Redoing with several more primes", level=2)
            M = prod * p*p*p
            continue

        if not proof:
            verbose("Not checking validity of result (since proof=False).", level=2)
            break
        d   = E.denom()
        hdE = long(E.height())
        if True or hdE * self.ncols() * height < prod:
            break
        M = prod * p*p*p
    #end while
    verbose("total time",tm, level=2)
    self.cache('pivots', best_pivots)
    E.cache('pivots', best_pivots)
    return E


###########################

def cmp_pivots(x,y):
    """
    Compare two sequences of pivot columns.
    If x is short than y, return -1, i.e., x < y, "not as good".
    If x is longer than y, x > y, "better"
    If the length is the same then x is better, i.e., x > y
        if the entries of x are correspondingly >= those of y with
        one being greater.
    I
    """
    if len(x) < len(y):
        return -1
    if len(x) > len(y):
        return 1
    if x < y:
        return 1
    elif x == y:
        return 0
    else:
        return -1

