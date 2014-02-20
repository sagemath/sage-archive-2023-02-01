
from sage.matrix.constructor import matrix
from sage.matrix.matrix import is_Matrix
from sage.rings.arith import legendre_symbol
from sage.rings.integer_ring import ZZ

def is_triangular_number(n):
    """
    Determines if the integer n is a triangular number.
    (I.e. determine if n = a*(a+1)/2 for some natural number a.)
    If so, return the number a, otherwise return False.

    Note: As a convention, n=0 is considered triangular for the
    number a=0 only (and not for a=-1).

    WARNING: Any non-zero value will return True, so this will test as
    True iff n is triangular and not zero.  If n is zero, then this
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
        from sage.functions.all import sqrt
        ## Try to solve for the integer a
        try:
            disc_sqrt = ZZ(sqrt(1+8*n))
            a = ZZ( (ZZ(-1) + disc_sqrt) / ZZ(2) )
            return a
        except Exception:
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

    ## Return the current information
    if  vec_output_flag:
        return D.rows()
    else:
        return D

def least_quadratic_nonresidue(p):
    """
    Returns the smallest positive integer quadratic non-residue in Z/pZ for primes p>2.

    EXAMPLES::

        sage: least_quadratic_nonresidue(5)
        2
        sage: [least_quadratic_nonresidue(p) for p in prime_range(3,100)]
        [2, 2, 3, 2, 2, 3, 2, 5, 2, 3, 2, 3, 2, 5, 2, 2, 2, 2, 7, 5, 3, 2, 3, 5]

    TESTS:

    Raises an error if input is a positive composite integer.

    ::

        sage: least_quadratic_nonresidue(20)
        Traceback (most recent call last):
        ...
        ValueError: Oops!  p must be a prime number > 2.


    Raises an error if input is 2. This is because every integer is a
    quadratic residue modulo 2.

    ::

        sage: least_quadratic_nonresidue(2)
        Traceback (most recent call last):
        ...
        ValueError: Oops!  There are no quadratic non-residues in Z/2Z.
    """
    from sage.functions.all import floor
    p1 = abs(p)

    ## Deal with the prime p = 2 and |p| <= 1.
    if p1 == 2:
        raise ValueError, "Oops!  There are no quadratic non-residues in Z/2Z."
    if p1 < 2:
        raise ValueError, "Oops!  p must be a prime number > 2."

    ## Find the smallest non-residue mod p
    ## For 7/8 of primes the answer is 2, 3 or 5:
    if p%8 in (3,5):
        return ZZ(2)
    if p%12 in (5,7):
        return ZZ(3)
    if p%5 in (2,3):
        return ZZ(5)
    ## default case (first needed for p=71):
    if not p.is_prime():
        raise ValueError, "Oops!  p must be a prime number > 2."
    from sage.misc.misc import xsrange
    for r in xsrange(7,p):
        if legendre_symbol(r, p) == -1:
            return ZZ(r)

