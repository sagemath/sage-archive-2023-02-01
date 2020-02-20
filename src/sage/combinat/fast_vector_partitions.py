r"""
ALGORITHM:
        Brent Yorgey, Generating Multiset Partitions,
        The Monad Reader, Issue 8, September 2007, p. 5.

        https://wiki.haskell.org/The_Monad.Reader/Previous_issues

AUTHORS:
        Denis K. Sunko (2020-02-19): initial version
"""
################################################################################
#            Copyright (C) 2020 Denis Sunko <dks@phy.hr>                       #
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

# Uncomment if using as stand-alone Python file:
#
#import operator

# Uncomment if using as stand-alone Python 3 file:
#
#from functools import reduce

def dickson_le(v, w):
    r"""
    Internal part of the current implementation of fast_vector_partitions().

    INPUT:
            - ``v``, ``w`` -- A pair of lists, understood as vectors.

    OUTPUT:
            True    if ``v`` <= ``w`` coordinatewise.

            False   otherwise.

    NOTES:
            B. Yorgey denotes this operator ``<|=``.

            It is the partial monomial ordering appearing in Dickson's lemma.

    EXAMPLES:
            ``v <|= w`` if and only if ``x^v`` divides ``x^w`` as monomials.

            All vectors ``v`` which make up a vector partition of ``w`` satisfy
            ``v <|= w``.
    """
    return reduce(operator.and_, map(operator.le, v, w))


def vector_unit(v):
    return (len(v) - 1)*[0] + [1]  # lexicographically smallest vector > 0


def vector_minus(a, b):
    return list(map(operator.sub, a, b))  # subtract coordinatewise


def vector_clip(a, b):
    return list(map(min, a, b))  # see line art below


def vector_halve(v):
    r"""
    Internal part of the current implementation of fast_vector_partitions().

    INPUT:
            - ``v`` -- A list of non-negative integers, understood as a vector.

    OUTPUT:
            A list, understood as the integer vector halfway down the list of
            lexicographically ordered vectors between between ``v`` and zero.

    NOTE:
            For vectors, ``v=a+b`` implies ``v=b+a``, which means that a
            downward search for such splittings, starting with ``v=v+0``, need
            only look as far as some "v/2", given precise meaning here.

            A similar logic is to stop the search for divisors of ``N`` at
            ``sqrt(N)``, halving the exponents in the prime decompostion.
            However, here "v/2" does not mean halving each coordinate.
    """
    i = 0
    result = []
    for vv in v:
        result += [vv // 2]
        i += 1
        if vv % 2 != 0:
            return result + v[i:] # the less significant part is just copied
    return result


def within(v):
    r"""
    Internal part of the current implementation of fast_vector_partitions().

    INPUT:
            - ``v`` -- A list of non-negative integers, understood as a vector.

    OUTPUT:
            Lexicographically ordered list of lists ``c`` such that
            ``c <|= v`` as vectors.
    
    NOTE:
            Used only for debugging.
    """
    v0 = [[n] for n in range(v[0], -1, -1)]
    vv = [range(n, -1, -1) for n in v[1:]]
    return reduce( lambda a, b: [aa + [bb] for aa in a for bb in b]
                 , vv, v0)


def debug_within_from_to(m, s, e):
    r"""
    Internal part of the current implementation of fast_vector_partitions().

    INPUT:
            - ``m``, ``s``, ``e``-- Three lists, understood as vectors.

    OUTPUT:
            Lexicographically ordered list of vectors ``v`` satisfying
    
            ``e <= v <= s`` and ``v <|= m``
    
    NOTE:
            Inefficient but obviously correct implementation.

            Used as a drop-in replacement for debugging within_from_to().
    """
    result = []
    for v in within(m):
        if e <= v <= s:
            result += [v]
    return result


def recursive_within_from_to(m, s, e, useS, useE):
    r"""
    Internal part of the current implementation of fast_vector_partitions().

    INPUT:
            - ``m``, ``s``, ``e`` -- Three lists, understood as vectors.
            - ``useS``, ``useE``  -- Two logical flags.

    OUTPUT:
            Lexicographically ordered list of lists ``v`` satisfying

            ``e <= v <= s`` and ``v <|= m`` as vectors.

    WARNING:
            The first call must be with ``useS==useE==True`` and ``s <|= m``.

    NOTE:
            The flags ``useS`` and ``useE`` are used to implement the condition
            efficiently. Because testing it loops over the vector, re-testing
            for each step as the vector grows is inefficient: all but the last
            comparison have been done cumulatively already. This code tests
            only for the last one, using the flags to accumulate information
            from previous calls.
    """
    if useS:
        start = s[0]
    else:
        start = m[0]

    if useE:
        end = e[0]
    else:
        end = 0
    result = []

    for x in range(start, end - 1, -1):
        useSS = useS and x == s[0]
        useEE = useE and x == e[0]
        if len(m) > 1:
            out = recursive_within_from_to(m[1:], s[1:], e[1:], useSS, useEE)
            for o in out:
                result += [[x] + o]
        else:
            result += [[x]] # we know the answer for singletons

    return result


def within_from_to(m, s, e):
    r"""
    Internal part of the current implementation of fast_vector_partitions().

    INPUT:
            - ``m``, ``s``, ``e``-- Three lists, understood as vectors.

    OUTPUT:
            Lexicographically ordered list of lists ``v`` satisfying

            ``e <= v <= s`` and ``v <|= m`` as vectors.

    WARNING:
            The input ``s`` will be "clipped" internally if it does not
            satisfy the condition ``s <|= m``. This behavior is transparent to
            the user, but may be puzzling when comparing outputs with the
            function recursive_within_from_to(), which has no such protection.

    NOTE:
            To understand the input check, some line art is helpful. Assume
            that ``(a,b)`` are the two least significant coordinates of some
            vector. Say

            ``e=(2,3), s=(7,6), m=(9,8)``.

            In the figure, these values are denoted by E, S, and M, while the
            letter X stands for all other allowed values of ``v=(a,b)``:

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

            If S moves horizontally, the full-height columns fill in the box
            until S reaches M, at which point it remains the limit in the
            b-direction as it moves out of the box, while M takes over as the
            limit in the a-direction, so the M-column remains filled only up to
            S, no matter how much S moves further to the right.

            If S moves vertically, its column will be filled to the top of the
            box, but it remains the relevant limit in the a-direction, while M
            takes over in the b-direction as S goes out of the box upwards.

            Both behaviors are captured by using the smaller coordinate of S
            and M, whenever S is outside the box defined by M. This is the
            "clipping" referred to in the warning above.
    """
    ss = s
    if not dickson_le(s, m):
        ss = vector_clip(s, m)
    if e > ss:
        return []
    return recursive_within_from_to(m, ss, e, True, True)


def recursive_vector_partitions(v, vL):
    r"""
    Internal part of the current implementation of fast_vector_partitions().

    INPUT:
            - ``v``, ``vL`` -- Two lists, understood as vectors.

    OUTPUT:
            Lexicographically ordered list of lists, each element representing
            a vector partition of ``v``, such that ``vL`` is lexicographically
            the smallest component used.
    """
    result = [[v]]
    #
    #   changes to within_from_to() may be tested by replacing
    #   within_from_to() with debug_within_from_to() in this line:
    #
    vspan = within_from_to(v, vector_halve(v), vL)
    for vv in vspan:
        for pp in recursive_vector_partitions(vector_minus(v, vv), vv):
            result += [[vv] + pp]
    return result


def fast_vector_partitions(v, min=None):
    r"""
    Brent Yorgey's fast algorithm for integer vector (multiset) partitions.

    INPUT:
            - ``v`` -- A list of non-negative integers, understood as the
                       vector to be partitioned.
            - ``min`` -- An optional list of non-negative integers, of same
                         length as ``v``.

    OUTPUT:
            A list of lists, each representing a vector partition of ``v``.
            
            If ``min`` is given, only partitions with parts ``p>=min`` in the
            lexicographic ordering will appear.
            
            If ``min`` is given and ``len(min)!=len(v)``, ``None`` is returned.

    NOTE:
            The partitions are returned as a list allocated in memory.

    WARNING:        
            The ordering of the partitions is reversed with respect to the
            output of Sage class VectorPartitions().

    EXAMPLES:
            # the older the computer, the more impressive the example:
            sage: %time fastvparts = fast_vector_partitions([6,6,6])
            CPU times: user 17 s, sys: 136 ms, total: 17.1 s
            Wall time: 17.1 s
            sage: %time vparts = list(VectorPartitions([6,6,6]))
            CPU times: user 4min 16s, sys: 319 ms, total: 4min 16s
            Wall time: 4min 16s
            sage: vparts == fastvparts[::-1]
            True
            sage: fast_vector_partitions([1,2,3], min = [0,1,1])
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
            sage: fast_vector_partitions([5,7,6], min = [1,3,2])==
            ....: list(VectorPartitions([5,7,6], min = [1,3,2]))[::-1]
            True
    """
    if min is None:
        return recursive_vector_partitions(v, vector_unit(v))
    else:
        if len(v)==len(min):
            return recursive_vector_partitions(v, min)
        else:
            return None
