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

.. [HadaSloa] \N.J.A. Sloane's Library of Hadamard Matrices, at
   http://neilsloane.com/hadamard/
.. [HadaWiki] Hadamard matrices on Wikipedia, :wikipedia:`Hadamard_matrix`
.. [Hora] \K. J. Horadam, Hadamard Matrices and Their Applications,
   Princeton University Press, 2006.
"""

#*****************************************************************************
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.integer_ring import ZZ
from sage.rings.integer import Integer
from sage.matrix.constructor import matrix, block_matrix, block_diagonal_matrix, diagonal_matrix
from urllib import urlopen
from sage.arith.all import is_prime, is_square, is_prime_power, divisors, kronecker_symbol
from math import sqrt
from sage.matrix.constructor import identity_matrix as I
from sage.matrix.constructor import ones_matrix     as J
from sage.misc.unknown import Unknown

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
    for i in range(H.ncols()):
        if H[0,i] < 0:
            H.rescale_col(i, -1)
    for i in range(H.nrows()):
        if H[i,0] < 0:
            H.rescale_row(i, -1)
    return H

def hadamard_matrix_paleyI(n, normalize=True):
    """
    Implements the Paley type I construction.

    The Paley type I case corresponds to the case `p \cong 3 \mod{4}` for a
    prime `p` (see [Hora]_).

    INPUT:

    - ``n`` -- the matrix size 

    - ``normalize`` (boolean) -- whether to normalize the result.

    EXAMPLES:

    We note that this method by default returns a normalised Hadamard matrix ::

        sage: from sage.combinat.matrices.hadamard_matrix import hadamard_matrix_paleyI
        sage: hadamard_matrix_paleyI(4)
        [ 1  1  1  1]
        [ 1 -1  1 -1]
        [ 1 -1 -1  1]
        [ 1  1 -1 -1]

    Otherwise, it returns a skew Hadamard matrix `H`, i.e. `H=S+I`, with
    `S=-S^\top`  ::

        sage: M=hadamard_matrix_paleyI(4, normalize=False); M
        [ 1  1  1  1]
        [-1  1  1 -1]
        [-1 -1  1  1]
        [-1  1 -1  1]
        sage: S=M-identity_matrix(4); -S==S.T
        True

    TESTS::

        sage: from sage.combinat.matrices.hadamard_matrix import is_hadamard_matrix
        sage: test_cases = [x+1 for x in range(100) if is_prime_power(x) and x%4==3]
        sage: all(is_hadamard_matrix(hadamard_matrix_paleyI(n),normalized=True,verbose=True)
        ....:     for n in test_cases)
        True
        sage: all(is_hadamard_matrix(hadamard_matrix_paleyI(n,normalize=False),verbose=True)
        ....:     for n in test_cases)
        True
    """
    p = n - 1
    if not(is_prime_power(p) and (p % 4 == 3)):
        raise ValueError("The order %s is not covered by the Paley type I construction." % n)

    from sage.rings.finite_rings.finite_field_constructor import FiniteField
    K = FiniteField(p,'x')
    K_list = list(K)
    K_list.insert(0,K.zero())
    H = matrix(ZZ, [[(1 if (x-y).is_square() else -1)
                     for x in K_list]
                    for y in K_list])
    for i in range(n):
        H[i,0] = -1
        H[0,i] =  1
    if normalize:
        for i in range(n):
            H[i,i] = -1
        H = normalise_hadamard(H)
    return H

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
        [ 1  1| 1  1| 1  1| 1  1| 1  1| 1  1]
        [ 1 -1|-1  1|-1  1|-1  1|-1  1|-1  1]
        [-----+-----+-----+-----+-----+-----]
        [ 1 -1| 1 -1| 1  1|-1 -1|-1 -1| 1  1]
        [ 1  1|-1 -1| 1 -1|-1  1|-1  1| 1 -1]
        [-----+-----+-----+-----+-----+-----]
        [ 1 -1| 1  1| 1 -1| 1  1|-1 -1|-1 -1]
        [ 1  1| 1 -1|-1 -1| 1 -1|-1  1|-1  1]
        [-----+-----+-----+-----+-----+-----]
        [ 1 -1|-1 -1| 1  1| 1 -1| 1  1|-1 -1]
        [ 1  1|-1  1| 1 -1|-1 -1| 1 -1|-1  1]
        [-----+-----+-----+-----+-----+-----]
        [ 1 -1|-1 -1|-1 -1| 1  1| 1 -1| 1  1]
        [ 1  1|-1  1|-1  1| 1 -1|-1 -1| 1 -1]
        [-----+-----+-----+-----+-----+-----]
        [ 1 -1| 1  1|-1 -1|-1 -1| 1  1| 1 -1]
        [ 1  1| 1 -1|-1  1|-1  1| 1 -1|-1 -1]

    TESTS::

        sage: from sage.combinat.matrices.hadamard_matrix import (hadamard_matrix_paleyII, is_hadamard_matrix)
        sage: test_cases = [2*(x+1) for x in range(50) if is_prime_power(x) and x%4==1]
        sage: all(is_hadamard_matrix(hadamard_matrix_paleyII(n),normalized=True,verbose=True)
        ....:     for n in test_cases)
        True
    """
    q = n//2 - 1
    if not(n%2==0 and is_prime_power(q) and (q % 4 == 1)):
        raise ValueError("The order %s is not covered by the Paley type II construction." % n)

    from sage.rings.finite_rings.finite_field_constructor import FiniteField
    K = FiniteField(q,'x')
    K_list = list(K)
    K_list.insert(0,K.zero())
    H = matrix(ZZ, [[(1 if (x-y).is_square() else -1)
                     for x in K_list]
                    for y in K_list])
    for i in range(q+1):
        H[0,i] = 1
        H[i,0] = 1
        H[i,i] = 0

    tr = { 0: matrix(2,2,[ 1,-1,-1,-1]),
           1: matrix(2,2,[ 1, 1, 1,-1]),
          -1: matrix(2,2,[-1,-1,-1, 1])}

    H = block_matrix(q+1,q+1,[tr[v] for r in H for v in r])

    return normalise_hadamard(H)

def is_hadamard_matrix(M, normalized=False, skew=False, verbose=False):
    r"""
    Test if `M` is a hadamard matrix.

    INPUT:

    - ``M`` -- a matrix

    - ``normalized`` (boolean) -- whether to test if ``M`` is a normalized
      Hadamard matrix, i.e. has its first row/column filled with +1.

    - ``skew`` (boolean) -- whether to test if ``M`` is a skew 
      Hadamard matrix, i.e. `M=S+I` for `-S=S^\top`, and `I` the identity matrix. 

    - ``verbose`` (boolean) -- whether to be verbose when the matrix is not
      Hadamard.

    EXAMPLE::

        sage: from sage.combinat.matrices.hadamard_matrix import is_hadamard_matrix
        sage: h = matrix.hadamard(12)
        sage: is_hadamard_matrix(h)
        True
        sage: from sage.combinat.matrices.hadamard_matrix import skew_hadamard_matrix
        sage: h=skew_hadamard_matrix(12)
        sage: is_hadamard_matrix(h, skew=True)
        True
        sage: h = matrix.hadamard(12)
        sage: h[0,0] = 2
        sage: is_hadamard_matrix(h,verbose=True)
        The matrix does not only contain +1 and -1 entries, e.g. 2
        False
        sage: h = matrix.hadamard(12)
        sage: for i in range(12):
        ....:     h[i,2] = -h[i,2]
        sage: is_hadamard_matrix(h,verbose=True,normalized=True)
        The matrix is not normalized
        False

    TESTS::

        sage: h = matrix.hadamard(12)
        sage: is_hadamard_matrix(h, skew=True)
        False
        sage: is_hadamard_matrix(h, skew=True, verbose=True)
        The matrix is not skew
        False
        sage: h=skew_hadamard_matrix(12)
        sage: is_hadamard_matrix(h, skew=True, verbose=True)
        True
        sage: is_hadamard_matrix(h, skew=False, verbose=True)
        True
        sage: h=-h
        sage: is_hadamard_matrix(h, skew=True, verbose=True)
        The matrix is not skew - diagonal entries must be all 1
        False
        sage: is_hadamard_matrix(h, skew=False, verbose=True)
        True
    """
    n = M.ncols()
    if n != M.nrows():
        if verbose:
            print "The matrix is not square ({}x{})".format(M.nrows(),n)
        return False

    if n == 0:
        return True

    for r in M:
        for v in r:
            if v*v != 1:
                if verbose:
                    print "The matrix does not only contain +1 and -1 entries, e.g. "+str(v)
                return False

    prod = (M*M.transpose()).dict()
    if (len(prod) != n or
        set(prod.itervalues()) != {n} or
        any( (i,i) not in prod for i in range(n) )):
        if verbose:
            print "The product M*M.transpose() is not equal to nI"
        return False

    if normalized:
        if (set(M.row(0)   ) != {1} or
            set(M.column(0)) != {1}):
            if verbose:
                print "The matrix is not normalized"
            return False

    if skew:
        for i in xrange(n-1):
            for j in xrange(i+1, n):
                if M[i,j] != -M[j,i]:
                    if verbose:
                        print "The matrix is not skew"
                    return False
        for i in xrange(n):
            if M[i,i] != 1:
                if verbose:
                    print "The matrix is not skew - diagonal entries must be all 1"
                return False
    return True

from sage.matrix.constructor import matrix_method
@matrix_method
def hadamard_matrix(n,existence=False, check=True):
    r"""
    Tries to construct a Hadamard matrix using a combination of Paley
    and Sylvester constructions.

    INPUT:

    - ``n`` (integer) -- dimension of the matrix

    - ``existence`` (boolean) -- whether to build the matrix or merely query if
      a construction is available in Sage. When set to ``True``, the function
      returns:

        - ``True`` -- meaning that Sage knows how to build the matrix

        - ``Unknown`` -- meaning that Sage does not know how to build the
          matrix, although the matrix may exist (see :mod:`sage.misc.unknown`).

        - ``False`` -- meaning that the matrix does not exist.

    - ``check`` (boolean) -- whether to check that output is correct before
      returning it. As this is expected to be useless (but we are cautious
      guys), you may want to disable it whenever you want speed. Set to ``True``
      by default.

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
        sage: hadamard_matrix(8) # random
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

        sage: hadamard_matrix(12) # random
        [ 1  1| 1  1| 1  1| 1  1| 1  1| 1  1]
        [ 1 -1|-1  1|-1  1|-1  1|-1  1|-1  1]
        [-----+-----+-----+-----+-----+-----]
        [ 1 -1| 1 -1| 1  1|-1 -1|-1 -1| 1  1]
        [ 1  1|-1 -1| 1 -1|-1  1|-1  1| 1 -1]
        [-----+-----+-----+-----+-----+-----]
        [ 1 -1| 1  1| 1 -1| 1  1|-1 -1|-1 -1]
        [ 1  1| 1 -1|-1 -1| 1 -1|-1  1|-1  1]
        [-----+-----+-----+-----+-----+-----]
        [ 1 -1|-1 -1| 1  1| 1 -1| 1  1|-1 -1]
        [ 1  1|-1  1| 1 -1|-1 -1| 1 -1|-1  1]
        [-----+-----+-----+-----+-----+-----]
        [ 1 -1|-1 -1|-1 -1| 1  1| 1 -1| 1  1]
        [ 1  1|-1  1|-1  1| 1 -1|-1 -1| 1 -1]
        [-----+-----+-----+-----+-----+-----]
        [ 1 -1| 1  1|-1 -1|-1 -1| 1  1| 1 -1]
        [ 1  1| 1 -1|-1  1|-1  1| 1 -1|-1 -1]

    TESTS::

        sage: matrix.hadamard(10,existence=True)
        False
        sage: matrix.hadamard(12,existence=True)
        True
        sage: matrix.hadamard(92,existence=True)
        Unknown
        sage: matrix.hadamard(10)
        Traceback (most recent call last):
        ...
        ValueError: The Hadamard matrix of order 10 does not exist
    """
    if not(n % 4 == 0) and (n > 2):
        if existence:
            return False
        raise ValueError("The Hadamard matrix of order %s does not exist" % n)
    if n == 2:
        if existence:
            return True
        M = matrix([[1, 1], [1, -1]])
    elif n == 1:
        if existence:
            return True
        M = matrix([1])
    elif is_prime_power(n//2 - 1) and (n//2 - 1) % 4 == 1:
        if existence:
            return True
        M = hadamard_matrix_paleyII(n)
    elif n == 4 or n % 8 == 0:
        if existence:
            return hadamard_matrix(n//2,existence=True)
        had = hadamard_matrix(n//2,check=False)
        chad1 = matrix([list(r) + list(r) for r in had.rows()])
        mhad = (-1) * had
        R = len(had.rows())
        chad2 = matrix([list(had.rows()[i]) + list(mhad.rows()[i])
                       for i in range(R)])
        M = chad1.stack(chad2)
    elif is_prime_power(n - 1) and (n - 1) % 4 == 3:
        if existence:
            return True
        M = hadamard_matrix_paleyI(n)
    else:
        if existence:
            return Unknown
        raise ValueError("The Hadamard matrix of order %s is not yet implemented." % n)

    if check:
        assert is_hadamard_matrix(M, normalized=True)

    return M

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

_rshcd_cache = {}

def regular_symmetric_hadamard_matrix_with_constant_diagonal(n,e,existence=False):
    r"""
    Return a Regular Symmetric Hadamard Matrix with Constant Diagonal.

    A Hadamard matrix is said to be *regular* if its rows all sum to the same
    value.

    For `\epsilon\in\{-1,+1\}`, we say that `M` is a `(n,\epsilon)-RSHCD` if
    `M` is a regular symmetric Hadamard matrix with constant diagonal
    `\delta\in\{-1,+1\}` and row sums all equal to `\delta \epsilon
    \sqrt(n)`. For more information, see [HX10]_ or 10.5.1 in
    [BH12]_. For the case `n=324`, see :func:`RSHCD_324` and [CP16]_.

    INPUT:

    - ``n`` (integer) -- side of the matrix

    - ``e`` -- one of `-1` or `+1`, equal to the value of `\epsilon`

    EXAMPLES::

        sage: from sage.combinat.matrices.hadamard_matrix import regular_symmetric_hadamard_matrix_with_constant_diagonal
        sage: regular_symmetric_hadamard_matrix_with_constant_diagonal(4,1)
        [ 1  1  1 -1]
        [ 1  1 -1  1]
        [ 1 -1  1  1]
        [-1  1  1  1]
        sage: regular_symmetric_hadamard_matrix_with_constant_diagonal(4,-1)
        [ 1 -1 -1 -1]
        [-1  1 -1 -1]
        [-1 -1  1 -1]
        [-1 -1 -1  1]

    Other hardcoded values::

        sage: for n,e in [(36,1),(36,-1),(100,1),(100,-1),(196, 1)]:
        ....:     print regular_symmetric_hadamard_matrix_with_constant_diagonal(n,e)
        36 x 36 dense matrix over Integer Ring
        36 x 36 dense matrix over Integer Ring
        100 x 100 dense matrix over Integer Ring
        100 x 100 dense matrix over Integer Ring
        196 x 196 dense matrix over Integer Ring

        sage: for n,e in [(324,1),(324,-1)]: # not tested - long time, tested in RSHCD_324
        ....:     print regular_symmetric_hadamard_matrix_with_constant_diagonal(n,e) # not tested - long time
        324 x 324 dense matrix over Integer Ring
        324 x 324 dense matrix over Integer Ring

    From two close prime powers::

        sage: print regular_symmetric_hadamard_matrix_with_constant_diagonal(64,-1)
        64 x 64 dense matrix over Integer Ring

    Recursive construction::

        sage: print regular_symmetric_hadamard_matrix_with_constant_diagonal(144,-1)
        144 x 144 dense matrix over Integer Ring

    REFERENCE:

    .. [BH12] \A. Brouwer and W. Haemers,
      Spectra of graphs,
      Springer, 2012,
      http://homepages.cwi.nl/~aeb/math/ipm/ipm.pdf

    .. [HX10] \W. Haemers and Q. Xiang,
      Strongly regular graphs with parameters `(4m^4,2m^4+m^2,m^4+m^2,m^4+m^2)` exist for all `m>1`,
      European Journal of Combinatorics,
      Volume 31, Issue 6, August 2010, Pages 1553-1559,
      http://dx.doi.org/10.1016/j.ejc.2009.07.009.
    """
    if existence and (n,e) in _rshcd_cache:
        return _rshcd_cache[n,e]

    from sage.graphs.strongly_regular_db import strongly_regular_graph

    def true():
        _rshcd_cache[n,e] = True
        return True

    M = None
    if abs(e) != 1:
        raise ValueError
    if n<0:
        if existence:
            return False
        raise ValueError
    elif n == 4:
        if existence:
            return true()
        if e == 1:
            M = J(4)-2*matrix(4,[[int(i+j == 3) for i in range(4)] for j in range(4)])
        else:
            M = -J(4)+2*I(4)
    elif n ==  36:
        if existence:
            return true()
        if e == 1:
            M = strongly_regular_graph(36, 15, 6, 6).adjacency_matrix()
            M = J(36) - 2*M
        else:
            M = strongly_regular_graph(36,14,4,6).adjacency_matrix()
            M =  -J(36) + 2*M + 2*I(36)
    elif n == 100:
        if existence:
            return true()
        if e == -1:
            M = strongly_regular_graph(100,44,18,20).adjacency_matrix()
            M = 2*M - J(100) + 2*I(100)
        else:
            M = strongly_regular_graph(100,45,20,20).adjacency_matrix()
            M = J(100) - 2*M
    elif n == 196 and e == 1:
        if existence:
            return true()
        M = strongly_regular_graph(196,91,42,42).adjacency_matrix()
        M = J(196) - 2*M
    elif n == 324:
        if existence:
            return true()
        M = RSHCD_324(e)
    elif (  e  == 1                 and
          n%16 == 0                 and
          is_square(n)              and
          is_prime_power(sqrt(n)-1) and
          is_prime_power(sqrt(n)+1)):
        if existence:
            return true()
        M = -rshcd_from_close_prime_powers(int(sqrt(n)))

    # Recursive construction: the kronecker product of two RSHCD is a RSHCD
    else:
        from itertools import product
        for n1,e1 in product(divisors(n)[1:-1],[-1,1]):
            e2 = e1*e
            n2 = n//n1
            if (regular_symmetric_hadamard_matrix_with_constant_diagonal(n1,e1,existence=True) and
                regular_symmetric_hadamard_matrix_with_constant_diagonal(n2,e2,existence=True)):
                if existence:
                    return true()
                M1 = regular_symmetric_hadamard_matrix_with_constant_diagonal(n1,e1)
                M2 = regular_symmetric_hadamard_matrix_with_constant_diagonal(n2,e2)
                M  = M1.tensor_product(M2)
                break

    if M is None:
        from sage.misc.unknown import Unknown
        _rshcd_cache[n,e] = Unknown
        if existence:
            return Unknown
        raise ValueError("I do not know how to build a {}-RSHCD".format((n,e)))

    assert M*M.transpose() == n*I(n)
    assert set(map(sum,M)) == {e*sqrt(n)}

    return M

def RSHCD_324(e):
    r"""
    Return a size 324x324 Regular Symmetric Hadamard Matrix with Constant Diagonal.

    We build the matrix `M` for the case `n=324`, `\epsilon=1` directly from
    :meth:`JankoKharaghaniTonchevGraph
    <sage.graphs.graph_generators.GraphGenerators.JankoKharaghaniTonchevGraph>`
    and for the case `\epsilon=-1` from the "twist" `M'` of `M`, using Lemma 11
    in [HX10]_. Namely, it turns out that the matrix

    .. math::

        M'=\begin{pmatrix} M_{12} & M_{11}\\ M_{11}^\top & M_{21} \end{pmatrix},
        \quad\text{where}\quad
        M=\begin{pmatrix} M_{11} & M_{12}\\ M_{21} & M_{22} \end{pmatrix},

    and the `M_{ij}` are 162x162-blocks, also RSHCD, its diagonal blocks having zero row
    sums, as needed by [loc.cit.]. Interestingly, the corresponding
    `(324,152,70,72)`-strongly regular graph
    has a vertex-transitive automorphism group of order 2592, twice the order of the
    (intransitive) automorphism group of the graph corresponding to `M`. Cf. [CP16]_.

    INPUT:

    - ``e`` -- one of `-1` or `+1`, equal to the value of `\epsilon`

    TESTS::

        sage: from sage.combinat.matrices.hadamard_matrix import RSHCD_324, is_hadamard_matrix
        sage: for e in [1,-1]: # long time
        ....:     M = RSHCD_324(e) # long time
        ....:     print M==M.T,is_hadamard_matrix(M),all([M[i,i]==1 for i in xrange(324)]) # long time
        ....:     print set(map(sum,M)) # long time
        True True True
        set([18])
        True True True
        set([-18])

    REFERENCE:

    .. [CP16] \N. Cohen, D. Pasechnik,
       Implementing Brouwer's database of strongly regular graphs,
       http://arxiv.org/abs/1601.00181
    """

    from sage.graphs.generators.smallgraphs import JankoKharaghaniTonchevGraph as JKTG
    M = JKTG().adjacency_matrix()
    M = J(324) - 2*M
    if e==-1:
        M1=M[:162].T
        M2=M[162:].T
        M11=M1[:162]
        M12=M1[162:].T
        M21=M2[:162].T
        M=block_matrix([[M12,-M11],[-M11.T,M21]])
    return M

def _helper_payley_matrix(n, zero_position=True):
    r"""
    Return the marix constructed in Lemma 1.19 page 291 of [SWW72]_.

    This function return a `n^2` matrix `M` whose rows/columns are indexed by
    the element of a finite field on `n` elements `x_1,...,x_n`. The value
    `M_{i,j}` is equal to `\chi(x_i-x_j)`.

    The elements `x_1,...,x_n` are ordered in such a way that the matrix
    (respectively, its submatrix obtained by removing first row and first column in the case
    ``zero_position=False``) is symmetric with respect to its second diagonal.
    The matrix is symmetric if `n=4k+1`, and skew-symmetric otherwise.

    INPUT:

    - ``n`` -- an odd prime power.

    - ``zero_position`` -- if it is true (default), place 0 of ``F_n`` in the middle,
      otherwise place it first.

    .. SEEALSO::

        :func:`rshcd_from_close_prime_powers`

    EXAMPLE::

        sage: from sage.combinat.matrices.hadamard_matrix import _helper_payley_matrix
        sage: _helper_payley_matrix(5)
        [ 0  1 -1 -1  1]
        [ 1  0  1 -1 -1]
        [-1  1  0  1 -1]
        [-1 -1  1  0  1]
        [ 1 -1 -1  1  0]

    TESTS::

        sage: _helper_payley_matrix(11,zero_position=True)
        [ 0 -1  1 -1 -1 -1  1  1  1 -1  1]
        [ 1  0 -1 -1  1 -1  1 -1  1  1 -1]
        [-1  1  0  1 -1 -1 -1 -1  1  1  1]
        [ 1  1 -1  0  1 -1 -1  1 -1 -1  1]
        [ 1 -1  1 -1  0  1 -1 -1 -1  1  1]
        [ 1  1  1  1 -1  0  1 -1 -1 -1 -1]
        [-1 -1  1  1  1 -1  0  1 -1  1 -1]
        [-1  1  1 -1  1  1 -1  0  1 -1 -1]
        [-1 -1 -1  1  1  1  1 -1  0 -1  1]
        [ 1 -1 -1  1 -1  1 -1  1  1  0 -1]
        [-1  1 -1 -1 -1  1  1  1 -1  1  0]
        sage: _helper_payley_matrix(11,zero_position=False)
        [ 0  1  1  1  1 -1  1 -1 -1 -1 -1]
        [-1  0 -1  1 -1 -1  1  1  1 -1  1]
        [-1  1  0 -1 -1  1  1 -1  1  1 -1]
        [-1 -1  1  0  1 -1 -1 -1  1  1  1]
        [-1  1  1 -1  0  1 -1  1 -1 -1  1]
        [ 1  1 -1  1 -1  0 -1 -1 -1  1  1]
        [-1 -1 -1  1  1  1  0  1 -1  1 -1]
        [ 1 -1  1  1 -1  1 -1  0  1 -1 -1]
        [ 1 -1 -1 -1  1  1  1 -1  0 -1  1]
        [ 1  1 -1 -1  1 -1 -1  1  1  0 -1]
        [ 1 -1  1 -1 -1 -1  1  1 -1  1  0]
    """
    from sage.rings.finite_rings.finite_field_constructor import FiniteField as GF
    K = GF(n,prefix='x')

    # Order the elements of K in K_list
    # so that K_list[i] = -K_list[n-i-1]
    K_pairs = set(frozenset([x,-x]) for x in K)
    K_pairs.discard(frozenset([0]))
    K_list = [None]*n
    if zero_position:
        zero_position=n//2
        shift=0
    else:
        shift=1

    for i,(x,y) in enumerate(K_pairs):
        K_list[i+shift]   = x
        K_list[-i-1] = y
    K_list[zero_position] = K(0)
    M = matrix(n,[[2*((x-y).is_square())-1
                   for x in K_list]
                  for y in K_list])
    M = M-I(n)
    assert (M*J(n)).is_zero()
    assert (M*M.transpose()) == n*I(n)-J(n)
    return M

def rshcd_from_close_prime_powers(n):
    r"""
    Return a `(n^2,1)`-RSHCD when `n-1` and `n+1` are odd prime powers and `n=0\pmod{4}`.

    The construction implemented here appears in Theorem 4.3 from [GS70]_.

    Note that the authors of [SWW72]_ claim in Corollary 5.12 (page 342) to have
    proved the same result without the `n=0\pmod{4}` restriction with a *very*
    similar construction. So far, however, I (Nathann Cohen) have not been able
    to make it work.

    INPUT:

    - ``n`` -- an integer congruent to `0\pmod{4}`

    .. SEEALSO::

        :func:`regular_symmetric_hadamard_matrix_with_constant_diagonal`

    EXAMPLES::

        sage: from sage.combinat.matrices.hadamard_matrix import rshcd_from_close_prime_powers
        sage: rshcd_from_close_prime_powers(4)
        [-1 -1  1 -1  1 -1 -1  1 -1  1 -1 -1  1 -1  1 -1]
        [-1 -1  1  1 -1 -1 -1 -1 -1  1  1 -1 -1  1 -1  1]
        [ 1  1 -1  1  1 -1 -1 -1 -1 -1  1 -1 -1 -1  1 -1]
        [-1  1  1 -1  1  1 -1 -1 -1 -1 -1  1 -1 -1 -1  1]
        [ 1 -1  1  1 -1  1  1 -1 -1 -1 -1 -1  1 -1 -1 -1]
        [-1 -1 -1  1  1 -1  1  1 -1 -1 -1  1 -1  1 -1 -1]
        [-1 -1 -1 -1  1  1 -1 -1  1 -1  1 -1  1  1 -1 -1]
        [ 1 -1 -1 -1 -1  1 -1 -1 -1  1 -1  1 -1  1  1 -1]
        [-1 -1 -1 -1 -1 -1  1 -1 -1 -1  1  1  1 -1  1  1]
        [ 1  1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1  1  1 -1  1]
        [-1  1  1 -1 -1 -1  1 -1  1 -1 -1 -1 -1  1  1 -1]
        [-1 -1 -1  1 -1  1 -1  1  1 -1 -1 -1 -1 -1  1  1]
        [ 1 -1 -1 -1  1 -1  1 -1  1  1 -1 -1 -1 -1 -1  1]
        [-1  1 -1 -1 -1  1  1  1 -1  1  1 -1 -1 -1 -1 -1]
        [ 1 -1  1 -1 -1 -1 -1  1  1 -1  1  1 -1 -1 -1 -1]
        [-1  1 -1  1 -1 -1 -1 -1  1  1 -1  1  1 -1 -1 -1]

    REFERENCE:

    .. [SWW72] A Street, W. Wallis, J. Wallis,
      Combinatorics: Room squares, sum-free sets, Hadamard matrices.
      Lecture notes in Mathematics 292 (1972).
    """
    if n%4:
        raise ValueError("n(={}) must be congruent to 0 mod 4")

    a,b = sorted([n-1,n+1],key=lambda x:-x%4)
    Sa  = _helper_payley_matrix(a)
    Sb  = _helper_payley_matrix(b)
    U   = matrix(a,[[int(i+j == a-1) for i in range(a)] for j in range(a)])

    K = (U*Sa).tensor_product(Sb) + U.tensor_product(J(b)-I(b)) - J(a).tensor_product(I(b))

    F = lambda x:diagonal_matrix([-(-1)**i for i in range(x)])
    G = block_diagonal_matrix([J(1),I(a).tensor_product(F(b))])
    e = matrix(a*b,[1]*(a*b))
    H = block_matrix(2,[-J(1),e.transpose(),e,K])

    HH = G*H*G
    assert len(set(map(sum,HH))) == 1
    assert HH**2 == n**2*I(n**2)
    return HH

def williamson_goethals_seidel_skew_hadamard_matrix(a, b, c, d, check=True):
    r"""
    Williamson-Goethals-Seidel construction of a skew Hadamard matrix

    Given `n\times n` (anti)circulant matrices `A`, `B`, `C`, `D` with 1,-1 entries, 
    and satisfying `A+A^\top = 2I`, `AA^\top + BB^\top + CC^\top + DD^\top = 4nI`,
    one can construct a skew Hadamard matrix of order `4n`, cf. [GS70s]_.

    INPUT:

    - ``a`` -- 1,-1 list specifying the 1st row of `A`

    - ``b`` -- 1,-1 list specifying the 1st row of `B`

    - ``d`` -- 1,-1 list specifying the 1st row of `C`

    - ``c`` -- 1,-1 list specifying the 1st row of `D`

    EXAMPLES::

        sage: from sage.combinat.matrices.hadamard_matrix import williamson_goethals_seidel_skew_hadamard_matrix as WGS
        sage: a=[ 1,  1, 1, -1,  1, -1,  1, -1, -1]
        sage: b=[ 1, -1, 1,  1, -1, -1,  1,  1, -1]
        sage: c=[-1, -1]+[1]*6+[-1]
        sage: d=[ 1,  1, 1, -1,  1,  1, -1,  1,  1]
        sage: M=WGS(a,b,c,d,check=True)

    REFERENCES:

    .. [GS70s] \J.M. Goethals and J. J. Seidel,
      A skew Hadamard matrix of order 36,
      J. Aust. Math. Soc. 11(1970), 343-344
    .. [Wall71] \J. Wallis,
      A skew-Hadamard matrix of order 92,
      Bull. Aust. Math. Soc. 5(1971), 203-204
    .. [KoSt08] \C. Koukouvinos, S. Stylianou
      On skew-Hadamard matrices,
      Discrete Math. 308(2008) 2723-2731


    """
    n = len(a)
    R = matrix(ZZ, n, n, lambda i,j: 1 if i+j==n-1 else 0)
    A,B,C,D=map(matrix.circulant, [a,b,c,d])
    if check:
        assert A*A.T+B*B.T+C*C.T+D*D.T==4*n*I(n)
        assert A+A.T==2*I(n)

    M = block_matrix([[   A,    B*R,    C*R,    D*R],
                      [-B*R,      A, -D.T*R,  C.T*R],
                      [-C*R,  D.T*R,      A, -B.T*R],
                      [-D*R, -C.T*R,  B.T*R,      A]])
    if check:
        assert is_hadamard_matrix(M, normalized=False, skew=True)
    return M

def GS_skew_hadamard_smallcases(n, existence=False, check=True):
    r"""
    Data for Williamson-Goethals-Seidel construction of skew Hadamard matrices

    Here we keep the data for this construction.
    Namely, it needs 4 circulant matrices with extra properties, as described in
    :func:`sage.combinat.matrices.hadamard_matrix.williamson_goethals_seidel_skew_hadamard_matrix`
    Matrices for `n=36` and `52` are given in [GS70s]_. Matrices for `n=92` are given
    in [Wall71]_.

    INPUT:

    - ``n`` -- the order of the matrix

    - ``existence`` -- if true (default), only check that we can do the construction

    - ``check`` -- if true (default), check the result.

    TESTS::

        sage: from sage.combinat.matrices.hadamard_matrix import GS_skew_hadamard_smallcases
        sage: GS_skew_hadamard_smallcases(36)
        36 x 36 dense matrix over Integer Ring...
        sage: GS_skew_hadamard_smallcases(52)
        52 x 52 dense matrix over Integer Ring...
        sage: GS_skew_hadamard_smallcases(92)
        92 x 92 dense matrix over Integer Ring...
        sage: GS_skew_hadamard_smallcases(100)
    """
    from sage.combinat.matrices.hadamard_matrix import\
         williamson_goethals_seidel_skew_hadamard_matrix as WGS
    def pmtoZ(s):
       return map(lambda x: 1 if x=='+' else -1, s)

    if existence:
        return n in [36, 52, 92]

    if n==36:
        a=[ 1,  1, 1, -1,  1, -1,  1, -1, -1]
        b=[ 1, -1, 1,  1, -1, -1,  1,  1, -1]
        c=[-1, -1]+[1]*6+[-1]
        d=[ 1,  1, 1, -1,  1,  1, -1,  1,  1]
        return WGS(a,b,c,d, check=check)
    if n==52:
        a=pmtoZ('++++-++--+---')
        b=pmtoZ('-+-++----++-+')
        c=pmtoZ('--+-+++++-+++')
        return WGS(a,b,c,c, check=check)
    if n==92:
        a = [1,-1,-1,-1,-1,-1,-1,-1, 1, 1,-1, 1,-1, 1,-1,-1, 1, 1, 1, 1, 1, 1, 1]
        b = [1, 1,-1,-1, 1,-1,-1, 1, 1, 1, 1,-1,-1, 1, 1, 1, 1,-1,-1, 1,-1,-1, 1]
        c = [1, 1,-1,-1,-1, 1,-1, 1,-1, 1,-1, 1, 1,-1, 1,-1, 1,-1, 1,-1,-1,-1, 1]
        d = [1,-1,-1,-1,-1, 1,-1,-1, 1,-1,-1, 1, 1,-1,-1, 1,-1,-1, 1,-1,-1,-1,-1]
        return WGS(a,b,c,d, check=check)
    return None

_skew_had_cache={}

def skew_hadamard_matrix(n,existence=False, skew_normalize=True, check=True):
    r"""
    Tries to construct a skew Hadamard matrix

    A Hadamard matrix `H` is called skew if `H=S-I`, for `I` the identity matrix
    and `-S=S^\top`. Currently constructions from Section 14.1 of [Ha83]_ and few
    more exotic ones are implemented.

    INPUT:

    - ``n`` (integer) -- dimension of the matrix

    - ``existence`` (boolean) -- whether to build the matrix or merely query if
      a construction is available in Sage. When set to ``True``, the function
      returns:

        - ``True`` -- meaning that Sage knows how to build the matrix

        - ``Unknown`` -- meaning that Sage does not know how to build the
          matrix, but that the design may exist (see :mod:`sage.misc.unknown`).

        - ``False`` -- meaning that the matrix does not exist.

    - ``skew_normalize`` (boolean) -- whether to make the 1st row all-one, and
      adjust the 1st column accordingly. Set to ``True`` by default.

    - ``check`` (boolean) -- whether to check that output is correct before
      returning it. As this is expected to be useless (but we are cautious
      guys), you may want to disable it whenever you want speed. Set to ``True``
      by default.

    EXAMPLES::

        sage: from sage.combinat.matrices.hadamard_matrix import skew_hadamard_matrix
        sage: skew_hadamard_matrix(12).det()
        2985984
        sage: 12^6
        2985984
        sage: skew_hadamard_matrix(1)
        [1]
        sage: skew_hadamard_matrix(2)
        [ 1  1]
        [-1  1]

    TESTS::

        sage: skew_hadamard_matrix(10,existence=True)
        False
        sage: skew_hadamard_matrix(12,existence=True)
        True
        sage: skew_hadamard_matrix(784,existence=True)
        True
        sage: skew_hadamard_matrix(10)
        Traceback (most recent call last):
        ...
        ValueError: A skew Hadamard matrix of order 10 does not exist
        sage: skew_hadamard_matrix(36)
        36 x 36 dense matrix over Integer Ring...
        sage: skew_hadamard_matrix(36)==skew_hadamard_matrix(36,skew_normalize=False)
        False
        sage: skew_hadamard_matrix(52)
        52 x 52 dense matrix over Integer Ring...
        sage: skew_hadamard_matrix(92)
        92 x 92 dense matrix over Integer Ring...
        sage: skew_hadamard_matrix(816)     # long time
        816 x 816 dense matrix over Integer Ring...
        sage: skew_hadamard_matrix(100)
        Traceback (most recent call last):
        ...
        ValueError: A skew Hadamard matrix of order 100 is not yet implemented.
        sage: skew_hadamard_matrix(100,existence=True)
        Unknown

    REFERENCES:

    .. [Ha83] \M. Hall,
      Combinatorial Theory,
      2nd edition,
      Wiley, 1983
    """
    def true():
        _skew_had_cache[n]=True
        return True
    M = None
    if existence and n in _skew_had_cache:
        return True
    if not(n % 4 == 0) and (n > 2):
        if existence:
            return False
        raise ValueError("A skew Hadamard matrix of order %s does not exist" % n)
    if n == 2:
        if existence:
            return true()
        M = matrix([[1, 1], [-1, 1]])
    elif n == 1:
        if existence:
            return true()
        M = matrix([1])
    elif is_prime_power(n - 1) and ((n - 1) % 4 == 3):
        if existence:
            return true()
        M = hadamard_matrix_paleyI(n, normalize=False)

    elif n % 8 == 0:
        if skew_hadamard_matrix(n//2,existence=True): # (Lemma 14.1.6 in [Ha83]_)
            if existence:
                return true()
            H = skew_hadamard_matrix(n//2,check=False)
            M = block_matrix([[H,H], [-H.T,H.T]])

        else: # try Williamson construction (Lemma 14.1.5 in [Ha83]_)
            for d in divisors(n)[2:-2]: # skip 1, 2, n/2, and n
                n1 = n//d
                if is_prime_power(d - 1) and (d % 4 == 0) and (n1 % 4 == 0)\
                    and skew_hadamard_matrix(n1,existence=True):
                    if existence:
                        return true()
                    H = skew_hadamard_matrix(n1, check=False)-I(n1)
                    U = matrix(ZZ, d, lambda i, j: -1 if i==j==0 else\
                                        1 if i==j==1 or (i>1 and j-1==d-i)\
                                          else 0)
                    A = block_matrix([[matrix([0]), matrix(ZZ,1,d-1,[1]*(d-1))],
                                      [ matrix(ZZ,d-1,1,[-1]*(d-1)),
                                        _helper_payley_matrix(d-1,zero_position=0)]])+I(d)
                    M = A.tensor_product(I(n1))+(U*A).tensor_product(H)
                    break
    if M is None: # try Williamson-Goethals-Seidel construction
        if GS_skew_hadamard_smallcases(n, existence=True):
            if existence:
                return true()
            M = GS_skew_hadamard_smallcases(n)

        else:
            if existence:
                return Unknown
            raise ValueError("A skew Hadamard matrix of order %s is not yet implemented." % n)
    if skew_normalize:
        dd = diagonal_matrix(M[0])
        M = dd*M*dd
    if check:
        assert is_hadamard_matrix(M, normalized=False, skew=True)
        if skew_normalize:
            from sage.modules.free_module_element import vector
            assert M[0]==vector([1]*n)
    _skew_had_cache[n]=True
    return M
