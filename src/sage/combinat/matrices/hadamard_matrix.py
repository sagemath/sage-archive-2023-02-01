r"""
Hadamard matrices

A Hadamard matrix is an `n\times n` matrix `H` whose entries are either `+1` or `-1`
and whose rows are mutually orthogonal. For example, the matrix `H_2`
defined by

.. MATH::

    \left(\begin{array}{rr}
    1 & 1 \\
    1 & -1
    \end{array}\right)

is a Hadamard matrix. An `n\times n` matrix `H` whose entries are either `+1` or
`-1` is a Hadamard matrix if and only if:

(a) `|det(H)|=n^{n/2}` or

(b)  `H*H^t = n\cdot I_n`, where `I_n` is the identity matrix.

In general, the tensor product of an `m\times m` Hadamard matrix and an
`n\times n` Hadamard matrix is an `(mn)\times (mn)` matrix. In
particular, if there is an `n\times n` Hadamard matrix then there is a
`(2n)\times (2n)` Hadamard matrix (since one may tensor with `H_2`).
This particular case is sometimes called the Sylvester construction.

The Hadamard conjecture (possibly due to Paley) states that a Hadamard
matrix of order `n` exists if and only if `n= 1, 2` or `n` is a multiple
of `4`.

The module below implements the Paley constructions (see for example
[Hora]_) and the Sylvester construction. It also allows you to pull a
Hadamard matrix from the database at [HadaSloa]_.

AUTHORS:

- David Joyner (2009-05-17): initial version

REFERENCES:

.. [HadaSloa] N.J.A. Sloane's Library of Hadamard Matrices, at
   http://neilsloane.com/hadamard/
.. [HadaWiki] Hadamard matrices on Wikipedia, :wikipedia:`Hadamard_matrix`
.. [Hora] K. J. Horadam, Hadamard Matrices and Their Applications,
   Princeton University Press, 2006.
"""
from sage.rings.arith import kronecker_symbol
from sage.rings.integer_ring import ZZ
from sage.rings.integer import Integer
from sage.matrix.constructor import matrix, block_matrix
from urllib import urlopen
from sage.misc.functional import is_even
from sage.rings.arith import is_prime


def H1(i, j, p):
    """
    Returns the i,j-th entry of the Paley matrix, type I case.
    The Paley type I case corresponds to the case `p \cong 3 \mod{4}`
    for a prime `p`.

    .. TODO::

        This construction holds more generally for prime powers `q`
        congruent to `3 \mod{4}`. We should implement these but we
        first need to implement Quadratic character for `GF(q)`.

    EXAMPLES::

        sage: sage.combinat.matrices.hadamard_matrix.H1(1,2,3)
        1
    """
    if i == 0 or j == 0:
        return 1
    # what follows will not be executed for (i, j) = (0, 0).
    if i == j:
        return -1
    return -kronecker_symbol(i - j, p)


def H2(i, j, p):
    """
    Returns the i,j-th entry of the Paley matrix, type II case.

    The Paley type II case corresponds to the case `p \cong 1 \mod{4}`
    for a prime `p` (see [Hora]_).

    .. TODO::

        This construction holds more generally for prime powers `q`
        congruent to `1 \mod{4}`. We should implement these but we
        first need to implement Quadratic character for `GF(q)`.

    EXAMPLES::

        sage: sage.combinat.matrices.hadamard_matrix.H2(1,2,5)
        1
    """
    if i == 0 and j == 0:
        return 0
    if i == 0 or j == 0:
        return 1
    if i == j:
        return 0
    return kronecker_symbol(i - j, p)

def normalise_hadamard(H):
    """
    Return the normalised Hadamard matrix corresponding to ``H``.

    The normalised Hadamard matrix corresponding to a Hadamard matrix `H` is a
    matrix whose every entry in the first row and column is +1.

    EXAMPLES::

        sage: H = sage.combinat.matrices.hadamard_matrix.normalise_hadamard(hadamard_matrix(4))
        sage: H == hadamard_matrix(4)
        True
    """
    Hc1 = H.column(0)
    Hr1 = H.row(0)
    for i in range(H.ncols()):
        if Hc1[i] < 0:
            H.rescale_row(i, -1)
    for i in range(H.nrows()):
        if Hr1[i] < 0:
            H.rescale_col(i, -1)
    return H

def hadamard_matrix_paleyI(n):
    """
    Implements the Paley type I construction.

    The Paley type I case corresponds to the case `p \cong 3 \mod{4}` for a
    prime `p` (see [Hora]_).

    EXAMPLES:

    We note that this method returns a normalised Hadamard matrix ::

        sage: sage.combinat.matrices.hadamard_matrix.hadamard_matrix_paleyI(4)
        [ 1  1  1  1]
        [ 1 -1 -1  1]
        [ 1  1 -1 -1]
        [ 1 -1  1 -1]
    """
    p = n - 1
    if not(is_prime(p) and (p % 4 == 3)):
        raise ValueError("The order %s is not covered by the Paley type I construction." % n)
    H = matrix(ZZ, [[H1(i, j, p) for i in range(n)] for j in range(n)])
    # normalising H so that first row and column have only +1 entries.
    return normalise_hadamard(H)

def hadamard_matrix_paleyII(n):
    """
    Implements the Paley type II construction.

    The Paley type II case corresponds to the case `p \cong 1 \mod{4}` for a
    prime `p` (see [Hora]_).

    EXAMPLES::

        sage: sage.combinat.matrices.hadamard_matrix.hadamard_matrix_paleyII(12).det()
        2985984
        sage: 12^6
        2985984

    We note that the method returns a normalised Hadamard matrix ::

        sage: sage.combinat.matrices.hadamard_matrix.hadamard_matrix_paleyII(12)
        [ 1  1  1  1  1  1| 1  1  1  1  1  1]
        [ 1  1  1 -1 -1  1|-1 -1  1 -1 -1  1]
        [ 1  1  1  1 -1 -1|-1  1 -1  1 -1 -1]
        [ 1 -1  1  1  1 -1|-1 -1  1 -1  1 -1]
        [ 1 -1 -1  1  1  1|-1 -1 -1  1 -1  1]
        [ 1  1 -1 -1  1  1|-1  1 -1 -1  1 -1]
        [-----------------+-----------------]
        [ 1 -1 -1 -1 -1 -1|-1  1  1  1  1  1]
        [ 1 -1  1 -1 -1  1| 1 -1 -1  1  1 -1]
        [ 1  1 -1  1 -1 -1| 1 -1 -1 -1  1  1]
        [ 1 -1  1 -1  1 -1| 1  1 -1 -1 -1  1]
        [ 1 -1 -1  1 -1  1| 1  1  1 -1 -1 -1]
        [ 1  1 -1 -1  1 -1| 1 -1  1  1 -1 -1]
    """
    N = Integer(n/2)
    p = N - 1
    if not(is_prime(p) and (p % 4 == 1)):
        raise ValueError("The order %s is not covered by the Paley type II construction." % n)
    S = matrix(ZZ, [[H2(i, j, p) for i in range(N)] for j in range(N)])
    H = block_matrix([[S + 1, S - 1], [1 - S, S + 1]])
    # normalising H so that first row and column have only +1 entries.
    return normalise_hadamard(H)

def hadamard_matrix(n):
    """
    Tries to construct a Hadamard matrix using a combination of Paley
    and Sylvester constructions.

    EXAMPLES::

        sage: hadamard_matrix(12).det()
        2985984
        sage: 12^6
        2985984
        sage: hadamard_matrix(1)
        [1]
        sage: hadamard_matrix(2)
        [ 1  1]
        [ 1 -1]
        sage: hadamard_matrix(8)
        [ 1  1  1  1  1  1  1  1]
        [ 1 -1  1 -1  1 -1  1 -1]
        [ 1  1 -1 -1  1  1 -1 -1]
        [ 1 -1 -1  1  1 -1 -1  1]
        [ 1  1  1  1 -1 -1 -1 -1]
        [ 1 -1  1 -1 -1  1 -1  1]
        [ 1  1 -1 -1 -1 -1  1  1]
        [ 1 -1 -1  1 -1  1  1 -1]
        sage: hadamard_matrix(8).det() == 8^4
        True

    We note that the method `hadamard_matrix()` returns a normalised Hadamard matrix
    (the entries in the first row and column are all +1) ::

        sage: hadamard_matrix(12)
        [ 1  1  1  1  1  1| 1  1  1  1  1  1]
        [ 1  1  1 -1 -1  1|-1 -1  1 -1 -1  1]
        [ 1  1  1  1 -1 -1|-1  1 -1  1 -1 -1]
        [ 1 -1  1  1  1 -1|-1 -1  1 -1  1 -1]
        [ 1 -1 -1  1  1  1|-1 -1 -1  1 -1  1]
        [ 1  1 -1 -1  1  1|-1  1 -1 -1  1 -1]
        [-----------------+-----------------]
        [ 1 -1 -1 -1 -1 -1|-1  1  1  1  1  1]
        [ 1 -1  1 -1 -1  1| 1 -1 -1  1  1 -1]
        [ 1  1 -1  1 -1 -1| 1 -1 -1 -1  1  1]
        [ 1 -1  1 -1  1 -1| 1  1 -1 -1 -1  1]
        [ 1 -1 -1  1 -1  1| 1  1  1 -1 -1 -1]
        [ 1  1 -1 -1  1 -1| 1 -1  1  1 -1 -1]
    """
    if not(n % 4 == 0) and (n > 2):
        raise ValueError("The Hadamard matrix of order %s does not exist" % n)
    if n == 2:
        return matrix([[1, 1], [1, -1]])
    if is_even(n):
        N = Integer(n / 2)
    elif n == 1:
        return matrix([1])
    if is_prime(N - 1) and (N - 1) % 4 == 1:
        return hadamard_matrix_paleyII(n)
    elif n == 4 or n % 8 == 0:
        had = hadamard_matrix(Integer(n / 2))
        chad1 = matrix([list(r) + list(r) for r in had.rows()])
        mhad = (-1) * had
        R = len(had.rows())
        chad2 = matrix([list(had.rows()[i]) + list(mhad.rows()[i])
                       for i in range(R)])
        return chad1.stack(chad2)
    elif is_prime(N - 1) and (N - 1) % 4 == 3:
        return hadamard_matrix_paleyI(n)
    else:
        raise ValueError("The Hadamard matrix of order %s is not yet implemented." % n)


def hadamard_matrix_www(url_file, comments=False):
    """
    Pulls file from Sloane's database and returns the corresponding Hadamard
    matrix as a Sage matrix.

    You must input a filename of the form "had.n.xxx.txt" as described
    on the webpage http://neilsloane.com/hadamard/, where
    "xxx" could be empty or a number of some characters.

    If comments=True then the "Automorphism..." line of the had.n.xxx.txt
    file is printed if it exists. Otherwise nothing is done.

    EXAMPLES::

        sage: hadamard_matrix_www("had.4.txt")      # optional - internet
        [ 1  1  1  1]
        [ 1 -1  1 -1]
        [ 1  1 -1 -1]
        [ 1 -1 -1  1]
        sage: hadamard_matrix_www("had.16.2.txt",comments=True)   # optional - internet
        Automorphism group has order = 49152 = 2^14 * 3
        [ 1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1]
        [ 1 -1  1 -1  1 -1  1 -1  1 -1  1 -1  1 -1  1 -1]
        [ 1  1 -1 -1  1  1 -1 -1  1  1 -1 -1  1  1 -1 -1]
        [ 1 -1 -1  1  1 -1 -1  1  1 -1 -1  1  1 -1 -1  1]
        [ 1  1  1  1 -1 -1 -1 -1  1  1  1  1 -1 -1 -1 -1]
        [ 1 -1  1 -1 -1  1 -1  1  1 -1  1 -1 -1  1 -1  1]
        [ 1  1 -1 -1 -1 -1  1  1  1  1 -1 -1 -1 -1  1  1]
        [ 1 -1 -1  1 -1  1  1 -1  1 -1 -1  1 -1  1  1 -1]
        [ 1  1  1  1  1  1  1  1 -1 -1 -1 -1 -1 -1 -1 -1]
        [ 1  1  1  1 -1 -1 -1 -1 -1 -1 -1 -1  1  1  1  1]
        [ 1  1 -1 -1  1 -1  1 -1 -1 -1  1  1 -1  1 -1  1]
        [ 1  1 -1 -1 -1  1 -1  1 -1 -1  1  1  1 -1  1 -1]
        [ 1 -1  1 -1  1 -1 -1  1 -1  1 -1  1 -1  1  1 -1]
        [ 1 -1  1 -1 -1  1  1 -1 -1  1 -1  1  1 -1 -1  1]
        [ 1 -1 -1  1  1  1 -1 -1 -1  1  1 -1 -1 -1  1  1]
        [ 1 -1 -1  1 -1 -1  1  1 -1  1  1 -1  1  1 -1 -1]
    """
    n = eval(url_file.split(".")[1])
    rws = []
    url = "http://neilsloane.com/hadamard/" + url_file
    f = urlopen(url)
    s = f.readlines()
    for i in range(n):
        r = []
        for j in range(n):
            if s[i][j] == "+":
                r.append(1)
            else:
                r.append(-1)
        rws.append(r)
    f.close()
    if comments:
        lastline = s[-1]
        if lastline[0] == "A":
            print lastline
    return matrix(rws)
