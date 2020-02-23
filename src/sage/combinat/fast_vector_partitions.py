# -*- coding: utf-8 -*-
r"""
Brent Yorgey's fast algorithm for integer vector (multiset) partitions.

ALGORITHM:

Brent Yorgey, Generating Multiset Partitions, The Monad Reader, Issue 8,
September 2007, p. 5.

https://wiki.haskell.org/The_Monad.Reader/Previous_issues

AUTHORS:

D.K. Sunko (2020-02-19): initial version
F. Chapoton (2020-02-22): conversion to iterators and shorter doctests and doc tweaks
"""
################################################################################
#          Copyright (C) 2020 Denis Sunko <dks@phy.hr>                         #
#          Copyright (C) 2020 Frédéric Chapoton <chapoton@math.unistra.fr>     #
#                                                                              #
# This program is free software: you can redistribute it and/or modify         #
# it under the terms of the GNU General Public License as published by         #
# the Free Software Foundation, either version 2 of the License, or            #
# (at your option) any later version. The text of the license is found at      #
#                  https://www.gnu.org/licenses/                               #
################################################################################
#
# To understand the code below, consult the ALGORITHM.
#
# Use at own risk.

def vector_halve(v):
    r"""
    Internal part of the current implementation of fast_vector_partitions().

    INPUT:

    - ``v`` -- list of non-negative integers, understood as a vector

    OUTPUT:

    A list, understood as the integer vector halfway down the list of
    lexicographically ordered vectors between between ``v`` and zero.

    EXAMPLES::

        sage: from sage.combinat.fast_vector_partitions import vector_halve
        sage: vector_halve([1, 2, 3, 4, 5, 6, 7, 8, 9])
        [0, 2, 3, 4, 5, 6, 7, 8, 9]
        sage: vector_halve([2, 4, 6, 8, 5, 6, 7, 8, 9])
        [1, 2, 3, 4, 2, 6, 7, 8, 9]

    .. NOTE::

        For vectors, ``v=a+b`` implies ``v=b+a``, which means that a
        downward search for such splittings, starting with ``v=v+0``, need
        only look as far as some "v/2", given precise meaning here.

        A similar logic is to stop the search for divisors of ``N`` at
        ``sqrt(N)``, halving the exponents in the prime decomposition.
        However, here "v/2" does not mean halving each coordinate.
    """
    result = []
    for i, vv in enumerate(v):
        result += [vv // 2]
        if vv % 2:
            return result + v[i+1:] # the less significant part is just copied
    return result


def recursive_within_from_to(m, s, e, useS, useE):
    r"""
    Internal part of the current implementation of fast_vector_partitions().

    INPUT:

    - ``m`` -- list of non-negative integers, understood as a vector
    - ``s`` -- list of non-negative integers, understood as a vector
    - ``e`` -- list of non-negative integers, understood as a vector
    - ``useS``  -- boolean
    - ``useE``  -- boolean

    OUTPUT:

    Lexicographically ordered list of lists ``v`` satisfying
    ``e <= v <= s`` and ``v <|= m`` as vectors.

    EXAMPLES::

        sage: from sage.combinat.fast_vector_partitions import recursive_within_from_to
        sage: list(recursive_within_from_to([1, 2, 3],[1, 2, 2],[1, 1, 1],True,True))
        [[1, 2, 2], [1, 2, 1], [1, 2, 0], [1, 1, 3], [1, 1, 2], [1, 1, 1]]

    .. NOTE::

        The flags ``useS`` and ``useE`` are used to implement the condition
        efficiently. Because testing it loops over the vector, re-testing
        for each step as the vector grows is inefficient: all but the last
        comparison have been done cumulatively already. This code tests
        only for the last one, using the flags to accumulate information
        from previous calls.

    .. WARNING::

        Expects to be called with ``s <|= m``.

        Expects to be called first with ``useS==useE==True``.
    """
    if useS:
        start = s[0]
    else:
        start = m[0]

    if useE:
        end = e[0]
    else:
        end = 0

    for x in range(start, end - 1, -1):
        useSS = useS and x == s[0]
        useEE = useE and x == e[0]
        if len(m) > 1:
            out = recursive_within_from_to(m[1:], s[1:], e[1:], useSS, useEE)
            for o in out:
                yield [x] + o
        else:
            yield [x] # we know the answer for singletons


def within_from_to(m, s, e):
    r"""
    Internal part of the current implementation of fast_vector_partitions().

    INPUT:

    - ``m`` -- list of non-negative integers, understood as a vector
    - ``s`` -- list of non-negative integers, understood as a vector
    - ``e`` -- list of non-negative integers, understood as a vector

    OUTPUT:

    Lexicographically ordered list of lists ``v`` satisfying
    ``e <= v <= s`` and ``v <|= m`` as vectors.
    
    EXAMPLES::

        sage: from sage.combinat.fast_vector_partitions import within_from_to
        sage: list(within_from_to([1, 2, 3],[1, 2, 2],[1, 1, 1]))
        [[1, 2, 2], [1, 2, 1], [1, 2, 0], [1, 1, 3], [1, 1, 2], [1, 1, 1]]
    
    .. NOTE::

        The input ``s`` will be "clipped" internally if it does not satisfy
        the condition ``s <|= m``. 

        To understand the input check, some line art is helpful. Assume
        that ``(a,b)`` are the two least significant coordinates of some
        vector. Say

        ``e=(2,3), s=(7,6), m=(9,8)``.

        In the figure, these values are denoted by E, S, and M, while the
        letter X stands for all other allowed values of ``v=(a,b)``::

            b ^
              |
            8 --------X---X---X---X---X-----------M
              |                                   |
            7 -       X   X   X   X   X           |
              |                                   |
            6 -       X   X   X   X   X   S       |
              |                                   |
            5 -       X   X   X   X   X   X       |
              |                                   |
            4 -       X   X   X   X   X   X       |
              |                                   |
            3 -       E   X   X   X   X   X       |
              |                                   |
            2 -           X   X   X   X   X       |
              |                                   |
            1 -           X   X   X   X   X       |
              |                                   |
            0 ----|---|---X---X---X---X---X---|---|--->
              0   1   2   3   4   5   6   7   8   9   a

        If S moves horizontally, the full-height columns fill the box in
        until S reaches M, at which point it remains the limit in the
        b-direction as it moves out of the box, while M takes over as the
        limit in the a-direction, so the M-column remains filled only up to
        S, no matter how much S moves further to the right.

        If S moves vertically, its column will be filled to the top of the
        box, but it remains the relevant limit in the a-direction, while M
        takes over in the b-direction as S goes out of the box upwards.

        Both behaviors are captured by using the smaller coordinate of S
        and M, whenever S is outside the box defined by M. The input will
        be "clipped" accordingly in that case.

    .. WARNING::

        The "clipping" behavior is transparent to the user, but may be puzzling
        when comparing outputs with the function recursive_within_from_to(),
        which has no input protection.
    """
    ss = s
    # if s is not in the box defined by m, we must clip:
    if not all(x <= y for x, y in zip(s, m)): # slightly slower without the if
        ss = [min(x,y) for x,y in zip(s, m)]  # rebuilding the list is costly
    if e > ss:
        return
    yield from recursive_within_from_to(m, ss, e, True, True)


def recursive_vector_partitions(v, vL):
    r"""
    Internal part of the current implementation of fast_vector_partitions().

    INPUT:

    - ``v`` -- list of non-negative integers, understood as a vector
    - ``vL`` -- list of non-negative integers, understood as a vector

    OUTPUT:

    Lexicographically ordered list of lists, each list representing
    a vector partition of ``v``, such that no part of any partition is
    lexicographically smaller than ``vL``.

    EXAMPLES::

        sage: from sage.combinat.fast_vector_partitions import recursive_vector_partitions
        sage: list(recursive_vector_partitions([2, 2, 2],[1, 1, 1]))
        [[[2, 2, 2]], [[1, 1, 1], [1, 1, 1]]]
        sage: list(recursive_vector_partitions([2, 2, 2],[1, 1, 0]))
        [[[2, 2, 2]], [[1, 1, 1], [1, 1, 1]], [[1, 1, 0], [1, 1, 2]]]
        sage: list(recursive_vector_partitions([2, 2, 2],[1, 0, 1]))
        [[[2, 2, 2]],
         [[1, 1, 1], [1, 1, 1]],
         [[1, 1, 0], [1, 1, 2]],
         [[1, 0, 2], [1, 2, 0]],
         [[1, 0, 1], [1, 2, 1]]]
    """
    yield [v]
    for vv in within_from_to(v, vector_halve(v), vL):
        v_minus_vv = [x - y for x, y in zip(v, vv)]
        for pp in recursive_vector_partitions(v_minus_vv, vv):
            yield [vv] + pp


def fast_vector_partitions(v, min=None):
    r"""
    Brent Yorgey's fast algorithm for integer vector (multiset) partitions.

    INPUT:

    - ``v``   -- list of non-negative integers, understood as the vector
                 to be partitioned

    - ``min`` -- optional list of non-negative integers, of same length
                 as ``v``

    OUTPUT:

    A list of lists, each representing a vector partition of ``v``.

    If ``min`` is given, only partitions with parts ``p>=min`` in the
    lexicographic ordering will appear.

    If ``min`` is given and ``len(min)!=len(v)``, ``None`` is returned.

    EXAMPLES:

    The older the computer, the more impressive the comparison::

        sage: from sage.combinat.fast_vector_partitions import fast_vector_partitions
        sage: fastvparts = list(fast_vector_partitions([3, 3, 3]))
        sage: vparts = list(VectorPartitions([3, 3, 3]))
        sage: vparts == fastvparts[::-1]
        True
        sage: len(fastvparts)
        686
        sage: list(fast_vector_partitions([1, 2, 3], min=[0, 1, 1]))
        [[[1, 2, 3]],
         [[0, 2, 3], [1, 0, 0]],
         [[0, 2, 2], [1, 0, 1]],
         [[0, 2, 1], [1, 0, 2]],
         [[0, 2, 0], [1, 0, 3]],
         [[0, 1, 3], [1, 1, 0]],
         [[0, 1, 2], [1, 1, 1]],
         [[0, 1, 1], [1, 1, 2]],
         [[0, 1, 1], [0, 1, 2], [1, 0, 0]],
         [[0, 1, 1], [0, 1, 1], [1, 0, 1]]]
        sage: list(fast_vector_partitions([5, 7, 6], min=[1, 3, 2])) == list(VectorPartitions([5, 7, 6], min = [1, 3, 2]))[::-1]
        True

    .. NOTE::

        The partitions are returned as an iterator.

        In this documentation, ``a <|= b`` means ``a[i] <= b[i]`` for all ``i``
        (notation following B. Yorgey's paper). It is the monomial partial
        ordering in Dickson's lemma: ``a <|= b`` iff ``x^a`` divides ``x^b`` as
        monomials.

    .. WARNING::

        The ordering of the partitions is reversed with respect to the output of
        Sage class VectorPartitions().
    """
    if min is None:
        min = (len(v) - 1) * [0] + [1]  # lexicographically smallest vector > 0
        return recursive_vector_partitions(v, min)
    elif len(v) == len(min):
            return recursive_vector_partitions(v, min)
    else:
        return
