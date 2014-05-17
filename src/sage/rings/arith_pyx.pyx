r"""
Other miscellaneous arithmetic implemented in C for speed.

AUTHORS:

- Vincent Delecroix (2014)
"""
#*****************************************************************************
#       Copyright (C) 2014 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 3 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

cdef extern from "math.h":
    double sqrt(double)

cimport integer
import integer

zero = integer.smallInteger(0)

cdef int two_squares_c(unsigned int n, unsigned int *p, unsigned int *q):
    r"""
    Return ``1`` if ``n`` is a sum of two squares and ``0`` otherwise.

    If ``1`` is returned, the the value of ``p`` and ``q`` are set to the
    lexicographically smallest solution of `p^2 + q^2 = n`.
    """
    cdef unsigned int fac,i,ii,j,jj,nn

    if n == 0:
        p[0] = q[0] = 0
        return 1

    # if n = 0 mod 4 then i and j must be even
    # hence, we first remove the maximum power of 4 from n and will then
    # multiply by the corresponding power of 2 the solution
    fac = 0
    while n%4 == 0:
        n >>= 2
        fac += 1

    # now, n is congruent to 1,2 or 3 mod 4.
    # As a square is congruent to 0,1 mod 4, a sum of square is congruent to
    # 0,1,2 mod 4.
    if n%4 == 3:
        return 0

    # if n=1 mod 4 then exactly one of i or j must be even
    # if n=2 mod 4 then i and j must be odd
    if n%4 == 1:
        i = ii = 0
        j = <unsigned int> sqrt(n) + 1
        jj = j*j
        while ii <= n/2:
            nn = n - ii
            while jj > nn:
                j -= 1
                # strangely enough, the 1-by-1 decreasing above is much faster
                # than integer Newton iteration:
                # j = (j+nn/j)/2
                jj = j*j
            if jj == nn:
                p[0] = i<<fac; q[0] = j<<fac
                return 1
            i += 1
            ii = i*i
    else: # n%4 = 2
        i = ii = 1
        j = <unsigned int> sqrt(n)
        j += (1-j%2)
        jj = j*j
        while ii <= n/2:
            nn = n - ii
            while jj > nn:
                j -= 2
                # strangely enough, the 2-by-2 decreasing above is much faster
                # than integer Newton iteration:
                # j = (j+nn/j)/2
                jj = j*j
            if jj == nn:
                p[0] = i<<fac; q[0] = j<<fac
                return 1
            i += 2
            ii = i*i

    return 0


cdef int three_squares_c(unsigned int n, unsigned int *p, unsigned int *q,  unsigned int *r):
    r"""
    Return ``1`` if ``n`` is a sum of three squares and ``0`` otherwise.

    If ``1`` is returned, the the value of ``p``, ``q`` and ``r`` are set to the
    lexicographically smallest solution of `p^2 + q^2 + r^2 = n`.
    """
    cdef unsigned int i,j,k,fac

    if n == 0:
        p[0] = q[0] = r[0] = 0
        return 1

    # if n == 0 mod 4 then i,j,k must be even
    # hence we remove from n the maximum power of 4 and at the very end we
    # multiply each term of the solution by the appropriate power of 2
    fac = 0
    while n%4 == 0:
        n >>= 2
        fac += 1

    # Legendre's three-square theorem: n is a sum of three squares if and only
    # if it is not of the form 4^a(8b + 7)
    if n%8 == 7:
        return 0

    i = 0
    while not two_squares_c(n-i*i, &j, &k):
        i += 1

    p[0] = i<<fac; q[0] = j<<fac; r[0] = k<<fac
    return 1

def two_squares_pyx(unsigned int n):
    r"""
    Return the lexicographically smallest pair of non-negative integers
    ``(i,j)`` such that `i^2 + j^2 = n`.

    If ``n`` is not a sum of two squares, a ``ValueError`` is raised.

    .. NOTE::

        The algorithm used here is relatively naive and only has interest for
        small values of ``n``. For that reason, the input must fit into an
        ``unsigned int`` (whose limit might be  `2^{32}-1=4294967295` or
        `2^{64}-1=18446744073709551615` depending on your computer and operating
        system).

    .. SEEALSO::

        :func:`~sage.arith.two_squares` is much more suited for large inputs

    EXAMPLES::

        sage: from sage.rings.arith_pyx import two_squares_pyx
        sage: two_squares_pyx(0)
        (0, 0)
        sage: two_squares_pyx(1)
        (0, 1)
        sage: two_squares_pyx(2)
        (1, 1)
        sage: two_squares_pyx(3)
        Traceback (most recent call last):
        ...
        ValueError: 3 is not a sum of 2 squares
        sage: two_squares_pyx(106)
        (5, 9)
    """
    cdef unsigned int i,j

    if two_squares_c(n, &i, &j):
        return (integer.smallInteger(i), integer.smallInteger(j))

    raise ValueError("%d is not a sum of 2 squares"%n)

def three_squares_pyx(unsigned int n):
    r"""
    If ``n`` is a sum of three squares return a 3-tuple ``(i,j,k)`` of Sage integers
    so that `i^2 + j^2 + k^2 = n`. Otherwise raise a ``ValueError``.

    .. NOTE::

        The algorithm used is relatively naive and only has interest for small
        values of ``n``. For that reason, the input must fit into an ``unsigned
        int`` (whose limit might be  `2^{32}-1=4294967295` or
        `2^{64}-1=18446744073709551615` depending on your plateform).

    .. SEEALSO::

        :func:`~sage.arith.three_squares` is much more suited for large input

    EXAMPLES::

        sage: from sage.rings.arith_pyx import three_squares_pyx
        sage: three_squares(0)
        (0, 0, 0)
        sage: three_squares(1)
        (0, 0, 1)
        sage: three_squares(2)
        (0, 1, 1)
        sage: three_squares(3)
        (1, 1, 1)
        sage: three_squares(4)
        (0, 0, 2)
        sage: three_squares(5)
        (0, 1, 2)
        sage: three_squares(6)
        (1, 1, 2)
        sage: three_squares(7)
        Traceback (most recent call last):
        ...
        ValueError: 7 is not a sum of 3 squares
        sage: three_squares(107)
        (1, 5, 9)
    """
    cdef unsigned int i,j,k

    if three_squares_c(n, &i, &j, &k):
        return (integer.smallInteger(i), integer.smallInteger(j), integer.smallInteger(k))

    raise ValueError("%d is not a sum of 3 squares"%n)

def four_squares_pyx(unsigned int n):
    r"""
    Return the lexicographically smallest 4-tuple of non-negative integers
    ``(i,j,k,l)`` such that `i^2 + j^2 + k^2 + l^2 = n`.

    .. NOTE::

        The algorithm used here is relatively naive and only has interest for
        small values of ``n``. For that reason, the input must fit into an
        ``unsigned int`` (whose limit might be  `2^{32}-1=4294967295` or
        `2^{64}-1=18446744073709551615` depending on your plateform).

    .. SEEALSO::

        :func:`~sage.arith.four_squares` is much more suited for large input

    EXAMPLES::

        sage: from sage.rings.arith_pyx import four_squares_pyx
        sage: four_squares_pyx(15447)
        (1, 1, 39, 118)
        sage: 1^2 + 1^2 + 39^2 + 118^2
        15447
        sage: four_squares_pyx(523439)
        (1, 3, 175, 702)
        sage: 1^2 + 3^2 + 273^2 + 670^2
        523439

    Because the output is the lexicographically smallest solution, the function
    can be used to see how many squares are needed to sum up to ``n``::

        sage: four_squares_pyx(169)
        (0, 0, 0, 13)
        sage: four_squares_pyx(113)
        (0, 0, 7, 8)
        sage: four_squares_pyx(114)
        (0, 1, 7, 8)
        sage: four_squares_pyx(119)
        (1, 1, 6, 9)

    TESTS::

        sage: all(sum(i**2 for i in four_squares_pyx(n)) == n for n in xrange(500,1000))
        True
        sage: four_squares_pyx(0)
        (0, 0, 0, 0)
        sage: four_squares_pyx(4)
        (0, 0, 0, 2)
        sage: four_squares_pyx(679)
        (1, 1, 1, 26)
        sage: four_squares_pyx(892)
        (1, 1, 7, 29)
    """
    cdef unsigned int j,k,l,fac

    if n == 0:
        return (zero, zero, zero, zero)

    # Because of Legendre three squares theorem, if we divide n by 4^a and
    # multiply the solution by 2^a we stil obtain the lexicographically smallest
    # solution
    fac = 0
    while n%4 == 0:
        n >>= 2
        fac += 1

    if n%8 == 7:
        three_squares_c(n-1, &j, &k, &l)
        return (integer.smallInteger(1<<fac), integer.smallInteger(j<<fac),
                integer.smallInteger(k<<fac), integer.smallInteger(l<<fac))
    else:
        three_squares_c(n, &j, &k, &l)
        return (zero, integer.smallInteger(j<<fac), integer.smallInteger(k<<fac), integer.smallInteger(l<<fac))
