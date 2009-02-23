

#*****************************************************************************
#       Copyright (C) 2007 William Stein and Jonathan Hanke
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


###########################################################################
## TO DO: Add routines for hasse invariants at all places, anisotropic
## places, is_semi_definite, and support for number fields.
###########################################################################


import copy

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.real_rqdf import RR
from sage.rings.arith import prime_divisors, valuation, hilbert_symbol
from sage.quadratic_forms.extras import IsPadicSquare, sgn
from sage.rings.fraction_field import FractionField
from sage.matrix.matrix_space import MatrixSpace
from sage.rings.arith import GCD


## Routines to compute local (p-adic) invariants of a quadratic form Q:
## (Note: Here Q is the matrix so that Q(x) = x^t * Q * x.)
## --------------------------------------------------------------------

def rational_diagonal_form(self, return_matrix=False):
    """
    Returns a diagonal form equivalent to Q over the fraction field of
    its defining ring.  If the return_matrix is True, then we return
    the transformation matrix performing the diagonalization as the
    second argument.

    OUTPUT:
        Q -- the diagonalized form of this quadratic form
        (optional) T -- matrix which diagonalizes Q (over it's fraction field)

    EXAMPLES:
        sage: Q = DiagonalQuadraticForm(ZZ, [1,3,5,7])
        sage: Q.rational_diagonal_form(return_matrix=True)
        (Quadratic form in 4 variables over Rational Field with coefficients:
        [ 1 0 0 0 ]
        [ * 3 0 0 ]
        [ * * 5 0 ]
        [ * * * 7 ]
        ,
         [1 0 0 0]
        [0 1 0 0]
        [0 0 1 0]
        [0 0 0 1])

        sage: Q1 = QuadraticForm(ZZ, 4, [1, 1, 0, 0, 1, 0, 0, 1, 0, 18])
        sage: Q1
        Quadratic form in 4 variables over Integer Ring with coefficients:
        [ 1 1 0 0 ]
        [ * 1 0 0 ]
        [ * * 1 0 ]
        [ * * * 18 ]
        sage: Q1.rational_diagonal_form(return_matrix=True)
        (Quadratic form in 4 variables over Rational Field with coefficients:
        [ 1 0 0 0 ]
        [ * 3/4 0 0 ]
        [ * * 1 0 ]
        [ * * * 18 ]
        ,
         [   1 -1/2    0    0]
        [   0    1    0    0]
        [   0    0    1    0]
        [   0    0    0    1])
    """
    n = self.dim()
    Q = copy.deepcopy(self)
    Q.__init__(FractionField(self.base_ring()), self.dim(), self.coefficients())
    MS = MatrixSpace(Q.base_ring(), n, n)
    T = MS(1)

    ## Construct an integral change of basis matrix T
    ## so that T^t * Q * T is diagonal.

    for i in range(n):
        temp = MS(1)
        for j in range(i+1, n):
            temp[j,j] = 1
            temp[i,j] = -Q[i,j] / ZZ(2) * Q[i,i]

        Q = Q(temp)
        T = temp * T

    ## Return the appropriate output
    if return_matrix:
        return Q, T
    else:
        return Q


def signature(self):
    """
    Returns the signature of the quadratic form, defined as:

       number of positive eigenvalues - number of negative eigenvalues

    of the matrix of the quadratic form.

    OUTPUT:
        an integer

    EXAMPLES:

    """
    diag = self.rational_diagonal_form()
    return sum([sgn(diag[i,i])  for i in range(diag.dim())])



def local_diagonal(self, p):
    """
    Finds a diagonal form equivalent to Q over the p-adic numbers Q_p.

    INPUT:
    p -- prime number in ZZ

    OUTPUT:
    A quadratic form whose only non-zero entries are on the diagonal.

    """
    ## TODO: Check that  is a prime number

    if (p != 2):
        return self.local_normal_form(p)  ## This is a diagonal matrix for Q over Z_p. =)
    else:
        Q2 = self.local_normal_form(p)
        for i in range(Q2.dim() - 1):
            if (Q2[i, i+1] != 0):
                ## This corresponds to [0 1]  =  [0  1/2]   ===>>>  [1  0]
                ##                     [0 0]     [1/2  0]           [0 -1]
                if (Q2[i,i] == 0):
                    Q2[i, i] = Q2[i, i+1]
                    Q2[i+1, i+1] = -Q2[i, i+1]
                    Q2[i, i+1] = 0

                ## This corresponds to [1 1]  =  [1  1/2]   ===>>>  [1 0]
                ##                     [0 1]     [1/2  1]           [0 3]
                else:
                    Q2[i, i+1] = 0
                    Q2[i+1, i+1] = Q2[i+1, i+1] * 3

    return Q2



def hasse_invariant(self, p):
    """
    Computes the Hasse invariant at a prime p, as given on p55 of
    Cassels book.  If Q is diagonal with coefficeints a_i, then the
    (Cassels) Hasse invariant is given by

        c_p = \prod_{i < j} (a_i, a_j)_p

    where (a,b)_p is the Hilbert symbol at p.

    WARNING: This is different from the O'Meara Hasse invariant, which
    allows i <= j in the product.  That is given by the method
    hasse_invariant__OMeara(p).

    NOTE: We should really rename this hasse_invariant__Cassels(), and
    set hasse_invariant() as a front-end to it.
    """
    ## TO DO: Need to deal with the case n=1 separately somewhere!

    Diag = self.local_diagonal(p);

    ## DIAGNOSTIC
    #print "\n Q = " + str(self)
    #print "\n Q diagonalized at p = " + str(p) + " gives " + str(Diag)

    hasse_temp = 1
    n = Diag.dim()

    for j in range(n-1):
        for k in range(j+1, n):
            hasse_temp = hasse_temp * hilbert_symbol(Diag[j,j], Diag[k,k], p)

    return hasse_temp


def hasse_invariant__OMeara(self, p):
    """
    Computes the O'Meara Hasse invariant at a prime p, as given on
    p167 of O'Meara's book.  If Q is diagonal with coefficeints a_i,
    then the (Cassels) Hasse invariant is given by

        c_p = \prod_{i <= j} (a_i, a_j)_p

    where (a,b)_p is the Hilbert symbol at p.

    WARNING: This is different from the (Cassels) Hasse invariant, which
    only allows i < j in the product.  That is given by the method
    hasse_invariant_OMeara(p)
    """
    ## TO DO: Need to deal with the case n=1 separately somewhere!

    Diag = self.local_diagonal(p);

    ## DIAGNOSTIC
    #print "\n Q = " + str(self)
    #print "\n Q diagonalized at p = " + str(p) + " gives " + str(Diag)

    hasse_temp = 1
    n = Diag.dim()

    for j in range(n-1):
        for k in range(j, n):
            hasse_temp = hasse_temp * hilbert_symbol(Diag[j,j], Diag[k,k], p)

    return hasse_temp




def is_hyperbolic(self, p):
    """
    Checks if the quadratic form is a sum of hyperbolic planes over
    the p-adic numbers Q_p.

    REFERENCES:
        This criteria follows from Cassels's "Rational Quadratic Forms":
            - local invariants for hyperbolic plane (Lemma 2.4, p58)
            - direct sum formulas (Lemma 2.3 on p58)

    EXAMPLES:
        sage: Q = DiagonalQuadraticForm(ZZ, [1,1])

        sage: Q.is_hyperbolic("infinity")
        False

        sage: Q.is_hyperbolic(2)
        False

        sage: Q.is_hyperbolic(3)
        False

        sage: Q.is_hyperbolic(5)     ## Here -1 is a square, so it's true.
        True

        sage: Q.is_hyperbolic(7)
        False

        sage: Q.is_hyperbolic(13)    ## Here -1 is a square, so it's true.
        True

    """
    ## False for odd-dim'l forms
    if self.dim() % 2 != 0:
        return False

    ## True for the zero form
    if self.dim == 0:
        return True

    ## Compare local invariants
    ## (Note: since the dimension is even, the extra powers of 2 in
    ##        self.det() := Det(2*Q) don't affect the answer!)
    m = ZZ(self.dim() / 2)
    if p == "infinity":
        return (self.signature() == 0)

    elif p == 2:
        return IsPadicSquare(self.det() * (-1)**m, p) and (self.hasse_invariant(p) == (-1)**m)    ## Actually, this -1 is the Hilbert symbol (-1,-1)_p

    else:
        return IsPadicSquare(self.det() * (-1)**m, p) and (self.hasse_invariant(p) == 1)



def is_anisotropic(self, p):
    """
    Checks if the quadratic form is anisotropic over the p-adic numbers Q_p.
    """
    n = self.dim()
    D = self.det()

    ## TO DO: Should check that p is prime

    if (n >= 5):
        return False;

    if (n == 4):
        return (IsPadicSquare(D, p) and (self.hasse_invariant(p) == - hilbert_symbol(-1,-1,p)) )

    if (n == 3):
        return (self.hasse_invariant(p) != hilbert_symbol(-1, -D, p))

    if (n == 2):
        return (not IsPadicSquare(-D, p))

    if (n == 1):
        return (self[0,0] != 0)

    raise NotImplementedError, "Oops!  We haven't established a convention for 0-dim'l quadratic forms... =("


def is_isotropic(self, p):
    """
    Checks if Q is isotropic over the p-adic numbers Q_p.
    """
    return not self.is_anisotropic(p)


def anisotropic_primes(self):
    """
    Returns a list with all of the anisotropic primes of the quadratic form.


    INPUT:
        None

    OUTPUT:
        Returns a list of prime numbers >0.

    EXAMPLES:
        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1])
        sage: Q.anisotropic_primes()
        [2]

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1,1])
        sage: Q.anisotropic_primes()
        [2]

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1,1,1])
        sage: Q.anisotropic_primes()
        []
    """

    ## Look at all prime divisors of 2 * Det(Q) to find the anisotropic primes...
    possible_primes = prime_divisors(2 * self.det())
    AnisoPrimes = []

    ## DIAGNSOTIC
    #print " Possible anisotropic primes are: " + str(possible_primes)

    for p in possible_primes:
        if (self.is_anisotropic(p)):
            AnisoPrimes += [p]

    ## DIAGNSOTIC
    #print " leaving anisotropic_primes..."

    return AnisoPrimes


def compute_definiteness(self):
    """
    Computes whether the given quadratic form is positive-definite,
    negative-definite, indefinite, degenerate, or the zero form.

    This caches one of the following strings in self.__definiteness_string:
    "pos_def", "neg_def", "indef", "zero", "degenerate".  It is called
    from all routines like:

        is_positive_definite(), is_negative_definite(), is_indefinite(), etc.

    Note:  A degenerate form is considered neither definite nor indefinite.
    Note:  The zero-dim'l form is considered both positive definite and negative definite.

    INPUT:
        None

    OUTPUT:
        None

    EXAMPLES:
        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1,1,1])
        sage: Q.is_positive_definite()
        True
        sage: Q.is_negative_definite()
        False
        sage: Q.is_indefinite()
        False
        sage: Q.is_definite()
        True

        sage: Q = DiagonalQuadraticForm(ZZ, [])
        sage: Q.is_positive_definite()
        True
        sage: Q.is_negative_definite()
        True
        sage: Q.is_indefinite()
        False
        sage: Q.is_definite()
        True

        sage: Q = DiagonalQuadraticForm(ZZ, [1,0,-1])
        sage: Q.is_positive_definite()
        False
        sage: Q.is_negative_definite()
        False
        sage: Q.is_indefinite()
        False
        sage: Q.is_definite()
        False

    """
    ## Sanity Check
    if not ((self.base_ring() == ZZ) or (self.base_ring() == QQ) or (self.base_ring() == RR)):
        raise NotImplementedError, "Oops!  We can only check definiteness over ZZ, QQ, and RR for now."

    ## Some useful variables
    n = self.dim()
    M = self.matrix()


    ## Deal with the zero-diml form
    if n == 0:
        self.__definiteness_string = "zero"
        return

    ## Deal with degenerate forms
    if self.det() == 0:
        self.__definiteness_string = "degenerate"
        return


    ## Check the sign of the ratios of consecutive determinants of the upper triangular r x r submatrices
    first_coeff = self[0,0]
    for r in range(1,n+1):
        I = range(r)
        new_det = M.matrix_from_rows_and_columns(I, I).det()

        ## Check for a (non-degenerate) zero -- so it's indefinite
        if new_det == 0:
            self.__definiteness_string = "indefinite"
            return

        ## Check for a change of signs in the upper r x r submatrix -- so it's indefinite
        if sgn(first_coeff)**r != sgn(new_det):
            self.__definiteness_string = "indefinite"
            return

    ## Here all ratios of determinants have the correct sign, so the matrix is (pos or neg) definite.
    if first_coeff > 0:
        self.__definiteness_string = "pos_def"
    else:
        self.__definiteness_string = "neg_def"




def is_positive_definite(self):
    """
    Determines if the given quadratic form is positive-definite.

    Note:  A degenerate form is considered neither definite nor indefinite.
    Note:  The zero-dim'l form is considered both positive definite and negative definite.

    INPUT:
        None

    OUTPUT:
        boolean -- True or False

    EXAMPLES:
        sage: Q = DiagonalQuadraticForm(ZZ, [1,3,5])
        sage: Q.is_positive_definite()
        True

        sage: Q = DiagonalQuadraticForm(ZZ, [1,-3,5])
        sage: Q.is_positive_definite()
        False

    """
    ## Try to use the cached value
    try:
        def_str = self.__definiteness_string
    except:
        self.compute_definiteness()
        def_str = self.__definiteness_string

    ## Return the answer
    return (def_str == "pos_def") or (def_str == "zero")




def is_negative_definite(self):
    """
    Determines if the given quadratic form is negative-definite.

    Note:  A degenerate form is considered neither definite nor indefinite.
    Note:  The zero-dim'l form is considered both positive definite and negative definite.

    INPUT:
        None

    OUTPUT:
        boolean -- True or False

    EXAMPLES:
        sage: Q = DiagonalQuadraticForm(ZZ, [-1,-3,-5])
        sage: Q.is_negative_definite()
        True

        sage: Q = DiagonalQuadraticForm(ZZ, [1,-3,5])
        sage: Q.is_negative_definite()
        False

    """
    ## Try to use the cached value
    try:
        def_str = self.__definiteness_string
    except:
        self.compute_definiteness()
        def_str = self.__definiteness_string

    ## Return the answer
    return (def_str == "neg_def") or (def_str == "zero")



def is_indefinite(self):
    """
    Determines if the given quadratic form is indefinite.

    Note:  A degenerate form is considered neither definite nor indefinite.
    Note:  The zero-dim'l form is not considered indefinite.

    INPUT:
        None

    OUTPUT:
        boolean -- True or False

    EXAMPLES:
        sage: Q = DiagonalQuadraticForm(ZZ, [-1,-3,-5])
        sage: Q.is_indefinite()
        False

        sage: Q = DiagonalQuadraticForm(ZZ, [1,-3,5])
        sage: Q.is_indefinite()
        True

    """
    ## Try to use the cached value
    try:
        def_str = self.__definiteness_string
    except:
        self.compute_definiteness()
        def_str = self.__definiteness_string

    ## Return the answer
    return def_str == "indefinite"


def is_definite(self):
    """
    Determines if the given quadratic form is (positive or negative) definite.

    Note:  A degenerate form is considered neither definite nor indefinite.
    Note:  The zero-dim'l form is considered indefinite.

    INPUT:
        None

    OUTPUT:
        boolean -- True or False

    EXAMPLES:
        sage: Q = DiagonalQuadraticForm(ZZ, [-1,-3,-5])
        sage: Q.is_indefinite()
        False

        sage: Q = DiagonalQuadraticForm(ZZ, [1,-3,5])
        sage: Q.is_indefinite()
        True

    """
    ## Try to use the cached value
    try:
        def_str = self.__definiteness_string
    except:
        self.compute_definiteness()
        def_str = self.__definiteness_string

    ## Return the answer
    return (def_str == "pos_def") or (def_str == "neg_def") or (def_str == "zero")


