
from random import random
from sage.calculus.calculus import floor
from sage.matrix.constructor import matrix
from sage.matrix.matrix import is_Matrix
from sage.rings.arith import valuation, kronecker_symbol, legendre_symbol, hilbert_symbol, is_prime
from sage.rings.rational_field import QQ
from sage.rings.integer_ring import ZZ
from sage.rings.infinity import infinity
from sage.calculus.calculus import sqrt
from sage.misc.functional import squarefree_part



def hilbert_symbol_rational(a, b, p):
    """
    Extend the usual Hilbert symbol to allow rational entries a and b.

    TO DO: This should really be incorporated into the hilbert_symbol()
    routine.

    EXAMPLES:
        sage: hilbert_symbol_rational(-1, -1, 2) == -1
        True
        sage: hilbert_symbol_rational(QQ(-1)/QQ(4), -1, 2) == -1
        True
        sage: hilbert_symbol_rational(QQ(-1)/QQ(4), -1, 3) == 1
        True
    """
    return hilbert_symbol(squarefree_part(a), squarefree_part(b), p)



def sgn(x):
    """
    Returns the sign of x. defined as:

                   /  1  if  x > 0,
        sgn(x) =   |  0  if  x = 0,
                   \ -1  if  x < 0.

    INPUT:
        a real number

    OUTPUT:
        1, 0, or -1.

    EXAMPLES:
        sage: sgn(pi) == 1
        True
        sage: sgn(5/6) == 1
        True
        sage: sgn(0) == 0
        True
        sage: sgn(-3) == -1
        True
    """
    if x > 0:
        return ZZ(1)
    elif x == 0:
        return ZZ(0)
    elif x < 0:
        return ZZ(-1)
    else:
        raise RuntimeError, "Oops!  We should not be here! =o"



def is_triangular_number(n):
    """
    Determines if the integer n is a triangular number.
    (I.e. determine if n = a*(a+1)/2 for some natural number a.)
    If so, return the number a, otherwise return False.

    Note: As a convetion, n=0 is consigered triangular for the
    number a=0 only (and not for a=-1).

    WARNING: Any non-zero value will return True, so this will test as
    True iff n is truangular and not zero.  If n is zero, then this
    will return the integer zero, which tests as False, so one must test

        if is_triangular_number(n) != False:

    instead of

        if is_triangular_number(n):

    to get zero to appear triangular.


    INPUT:
        an integer

    OUTPUT:
        either False or a non-negative integer

    EXAMPLES:
        sage: is_triangular_number(3)
        2
        sage: is_triangular_number(1)
        1
        sage: is_triangular_number(2)
        False
        sage: is_triangular_number(0)
        0
        sage: is_triangular_number(-1)
        False
        sage: is_triangular_number(-11)
        False
        sage: is_triangular_number(-1000)
        False
        sage: is_triangular_number(-0)
        0
        sage: is_triangular_number(10^6 * (10^6 +1)/2)
        1000000
    """
    if n < 0:
        return False
    elif n == 0:
        return ZZ(0)
    else:
        ## Try to solve for the integer a
        try:
            disc_sqrt = ZZ(sqrt(1+8*n))
            a = ZZ( (ZZ(-1) + disc_sqrt) / ZZ(2) )
            return a
        except:
            return False



def extend_to_primitive(A_input):
    """
    Given a matrix (resp. list of vectors), extend it to a square
    matrix (resp. list of vectors), such that its determinant is the
    gcd of its minors (i.e. extend the basis of a lattice to a
    "maximal" one in Z^n).

    Author(s): Gonzalo Tornaria and Jonathan Hanke.

    INPUT:
        a matrix, or a list of length n vectors (in the same space)

    OUTPUT:
        a square matrix, or a list of n vectors (resp.)

    EXAMPLES:
        sage: A = Matrix(ZZ, 3, 2, range(6))
        sage: extend_to_primitive(A)
        [ 0  1  0]
        [ 2  3  0]
        [ 4  5 -1]

        sage: extend_to_primitive([vector([1,2,3])])
        [(1, 2, 3), (0, 1, 0), (0, 0, 1)]

    """
    ## Deal with a list of vectors
    if not is_Matrix(A_input):
        A = matrix(A_input)      ## Make a matrix A with the given rows.
        vec_output_flag = True
    else:
        A = A_input
        vec_output_flag = False


    ## Arrange for A  to have more columns than rows.
    if A.is_square():
        return A
    if A.nrows() > A.ncols():
        return extend_to_primitive(A.transpose()).transpose()

    ## Setup
    k = A.nrows()
    n = A.ncols()
    R = A.base_ring()

    # Smith normal form transformation, assuming more columns than rows
    V = A.smith_form()[2]

    ## Extend the matrix in new coordinates, then switch back.
    B = A * V
    B_new = matrix(R, n-k, n)
    for i in range(n-k):
        B_new[i, n-i-1] = 1
    C = B.stack(B_new)
    D = C * V**(-1)

    ## DIAGNOSTIC
    #print "A = ", A, "\n"
    #print "B = ", B, "\n"
    #print "C = ", C, "\n"
    #print "D = ", D, "\n"

    # Normalize for a positive determinant
    if D.det() < 0:
        D.rescale_row(n-1, -1)

    ## Return the corrent information
    if  vec_output_flag:
        return D.rows()
    else:
        return D




def random_int_upto(n):
    """
    Returns a random integer x satisfying 0 <= x < n.

    EXAMPLES:
        sage: x = random_int_upto(10)
        sage: x >= 0
        True
        sage: x < 10
        True
    """
    return floor(n * random())



def quadratic_nonresidue(p):
    """
    Returns the smalest positive integer quadratic non-residue in Z/pZ for primes p>2.

    EXAMPLES:
        sage: quadratic_nonresidue(5)
        2
    """
    p1 = abs(p)

    ## Deal with the prime p = 2 and |p| <= 1.
    if p1 == 2:
        raise TypeError, "Oops!  There are no quadratic non-residues in Z/2Z."
    if p1 < 2:
        raise TypeError, "Oops!  p must be a prime number > 2."

    ## TO DO: Test that p is prime???

    ## Find the smallest non-residue mod p
    for r in range(2,p):
        if legendre_symbol(r, p) == -1:
            return r


def IsPadicSquare(m, p):
    """
    Determines whether the (non-zero) rational number m is a square in Q_p.
    When p = infinity this returns the answer for the real numbers.

    EXAMPLES:
        sage: IsPadicSquare(2, 7)
        True
        sage: IsPadicSquare(98, 7)
        True
        sage: IsPadicSquare(2, 5)
        False
    """
    ## Make sure m is non-zero
    if m == 0:
        raise TypeError, "Oops!  We're not allowed to ask about zero right now..."

    ## Deal with p = infinity (i.e. the real numbers)
    if p == infinity:
        return (m > 0)

    ## Check that p is prime
    try:
        is_prime(p)
    except:
        raise TypeError, 'Oops!  p must be "infinity" or a positive prime number.'


    ## Deal with finite primes
    v1 = valuation(QQ(m).numer(), p)
    v2 = valuation(QQ(m).denom(), p)
    v = v1 - v2
    other = m / (p**v)

    #print "p =", p
    #print "m =", m
    #print "val(m,p) =", v

    if ((v % 2) == 0):

        if ((p == 2) and ((other % 8) == 1)):
            return True

        if ((p > 2) and (kronecker_symbol(other, p) == 1)):
            return True

    return False
