r"""
Enumerated set of lists of integers with constraints, in inverse lexicographic order

- :class:`IntegerListsLex`

HISTORY:

This generic tool was originally written by Hivert and Thiery in
MuPAD-Combinat in 2000 and ported over to Sage by Mike Hansen in
2007. It was then completely rewritten in 2015 by Gillespie,
Schilling, and Thiery, with the help of many, to catter for
limitations and lack of robustness w.r.t. input. The old
implementation is still available in
:mod:`sage.combinat.integer_list.old` for benchmarking purposes.
"""
#*****************************************************************************
#       Copyright (C) 2015 Bryan Gillespie <Brg008@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import collections
from warnings import warn
from sage.misc.classcall_metaclass import ClasscallMetaclass, typecall
from sage.misc.constant_function import ConstantFunction
from sage.misc.cachefunc import cached_method
from sage.categories.enumerated_sets import EnumeratedSets
from sage.structure.list_clone import ClonableArray
from sage.structure.parent import Parent
from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
from sage.sets.family import Family
from sage.rings.integer_ring import ZZ

Infinity = float('+inf')

class IntegerListsLex(Parent):
    r"""
    Lists of non negative integers with constraints, in inverse lexicographic order.

    An *integer list* is a list `l` of nonnegative integers, its
    *parts*. The *length* ``len(l)`` of `l` is the number of its
    parts. The *sum* `|l|` of `l` is the sum of its parts. The slope
    (at position `i`) is the difference ``l[i+1]-l[i]`` between two
    consecutive parts.

    This class allows to construct the set of all integer lists `l`
    satisfying specified bounds on the sum, the length, the slope, and
    the individual parts, enumerated in *inverse* lexicographic order,
    that is from largest to smallest in lexicographic order.

    The main purpose is to provide a generic iteration engine for all
    the enumerated sets like :class:`Partitions`,
    :class:`Compositions`, :class:`IntegerVectors`. It can also be
    used to generate many other combinatorial objects like Dyck paths,
    Motzkin paths, etc.

    Mathematically speaking, this is a special case of sets of
    integral points of a polytope (or union thereof, when the length
    is not fixed). The set of allowable constraints has been
    specifically designed to enable iteration with a good time and
    memory complexity in most practical use cases, and in inverse
    lexicographic order (see below).

    INPUT:

    - ``min_sum`` -- a nonnegative integer (default: 0):
      a lower bound on `|l|`.

    - ``max_sum`` -- a nonnegative integer or `\infty` (default: `\infty`):
      an upper bound on `|l|`.

    - ``n`` -- a nonnegative integer (optional): if specified, this
      overrides ``min_sum`` and ``max_sum``. Alternatively a list or
      iterable of nonnegative integers can be specified.

    - ``min_length`` -- a nonnegative integer (default: `0`): a lower
      bound on ``len(l)``.

    - ``max_length`` -- a nonnegative integer or `\infty` (default:
      `\infty`): an upper bound on ``len(l)``.

    - ``length`` -- an integer (optional); overrides ``min_length``
      and ``max_length`` if specified;

    - ``min_part`` -- a nonnegative integer: a lower bounds on the
       parts: ``min_part <= l[i]`` for ``0 <= i < len(l)``.

    - ``floor`` -- a list of nonnegative integers or a function: lower
      bounds on the individual parts `l[i]`.

      If ``floor`` is a list of integers, then ``floor<=l[i]`` for ``0
      <= i < min(len(l), len(floor)``. Similarly, if ``floor`` is a
      function, then ``floor(i) <= l[i]`` for ``0 <= i < len(l)``.

    - ``max_part`` -- a nonnegative integer or `\infty`: an upper
      bound on the parts: ``l[i] <= max_part`` for ``0 <= i < len(l)``.

    - ``ceiling`` -- upper bounds on the individual parts ``l[i]``;
      this takes the same type of input as ``floor``, except that
      `\infty` is allowed in addition to integers, and the default
      value is `\infty`.

    - ``min_slope`` -- an integer or `-\infty` (default: `-\infty`):
      an lower bound on the slope between consecutive parts:
      ``min_slope <= l[i+1]-l[i]`` for ``0 <= i < len(l)-1``

    - ``max_slope`` -- an integer or `+\infty` (defaults: `+\infty`)
      an upper bound on the slope between consecutive parts:
      `` l[i+1]-l[i] <= max_slope`` for ``0 <= i < len(l)-1``

    - ``category`` -- a category (default: :class:`FiniteEnumeratedSets`)

    - ``check`` -- boolean (default: True): whether to display the
      warnings raised when functions are given as input to ``floor``
      or ``ceiling`` and the errors raised when there is no proper
      enumeration.

    .. NOTE::

        When several lists satisfying the constraints differ only by
        trailing zeroes, only the shortest one is enumerated (and
        therefore counted). The others are still considered valid.
        See the examples below.

        This feature is questionable. It is recommended not to rely on
        it, as it may eventually be discontinued.

    EXAMPLES:

    We create the enumerated set of all lists of non negative integers
    of length `3` and sum `2`::

        sage: C = IntegerListsLex(2, length=3)
        sage: C
        Integer lists of sum 2 satisfying certain constraints
        sage: C.cardinality()
        6
        sage: [p for p in C]
        [[2, 0, 0], [1, 1, 0], [1, 0, 1], [0, 2, 0], [0, 1, 1], [0, 0, 2]]

        sage: [2, 0, 0] in C
        True
        sage: [2, 0, 1] in C
        False
        sage: "a" in C
        False
        sage: ["a"] in C
        False
        sage: C.first()
        [2, 0, 0]

    One can specify lower and upper bounds on each part::

        sage: list(IntegerListsLex(5, length=3, floor=[1,2,0], ceiling=[3,2,3]))
        [[3, 2, 0], [2, 2, 1], [1, 2, 2]]

    When the length is fixed as above, one can also use
    :class:`IntegerVectors`::

        sage: IntegerVectors(2,3).list()
        [[2, 0, 0], [1, 1, 0], [1, 0, 1], [0, 2, 0], [0, 1, 1], [0, 0, 2]]

    Using the slope condition, one can generate integer partitions
    (but see :class:`Partitions`)::

        sage: list(IntegerListsLex(4, max_slope=0))
        [[4], [3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1]]

    The following is the list of all partitions of `7` with parts at least `2`::

        sage: list(IntegerListsLex(7, max_slope=0, min_part=2))
        [[7], [5, 2], [4, 3], [3, 2, 2]]


    .. RUBRIC:: floor and ceiling conditions

    Next we list all partitions of `5` of length at most `3` which are
    bounded below by ``[2,1,1]``::

        sage: list(IntegerListsLex(5, max_slope=0, max_length=3, floor=[2,1,1]))
        [[5], [4, 1], [3, 2], [3, 1, 1], [2, 2, 1]]

    Note that ``[5]`` is considered valid, because the floor
    constraints only apply to existing positions in the list. To
    obtain instead the partitions containing ``[2,1,1]``, one needs to
    use ``min_length`` or ``length``::

        sage: list(IntegerListsLex(5, max_slope=0, length=3, floor=[2,1,1]))
        [[3, 1, 1], [2, 2, 1]]

    Here is the list of all partitions of `5` which are contained in
    ``[3,2,2]``::

        sage: list(IntegerListsLex(5, max_slope=0, max_length=3, ceiling=[3,2,2]))
        [[3, 2], [3, 1, 1], [2, 2, 1]]

    This is the list of all compositions of `4` (but see :class:`Compositions`)::

        sage: list(IntegerListsLex(4, min_part=1))
        [[4], [3, 1], [2, 2], [2, 1, 1], [1, 3], [1, 2, 1], [1, 1, 2], [1, 1, 1, 1]]

    This is the list of all integer vectors of sum `4` and length `3`::

        sage: list(IntegerListsLex(4, length=3))
        [[4, 0, 0], [3, 1, 0], [3, 0, 1], [2, 2, 0], [2, 1, 1],
         [2, 0, 2], [1, 3, 0], [1, 2, 1], [1, 1, 2], [1, 0, 3],
         [0, 4, 0], [0, 3, 1], [0, 2, 2], [0, 1, 3], [0, 0, 4]]

    For whatever it's worth, the ``floor`` and ``min_part``
    constraints can be combined::

        sage: L = IntegerListsLex(5, floor=[2,0,2], min_part=1)
        sage: L.list()
        [[5], [4, 1], [3, 2], [2, 3], [2, 1, 2]]

    This is achieved by updating the floor upon constructing ``L``::

        sage: [L._floor(i) for i in range(5)]
        [2, 1, 2, 1, 1]

    Similarly, the ``ceiling`` and ``max_part`` constraints can be
    combined::

        sage: L = IntegerListsLex(4, ceiling=[2,3,1], max_part=2, length=3)
        sage: L.list()
        [[2, 2, 0], [2, 1, 1], [1, 2, 1]]
        sage: [L._ceiling(i) for i in range(5)]
        [2, 2, 1, 2, 2]


    This can be used to generate Motzkin words (see
    :wikipedia:`Motzkin_number`)::

        sage: def motzkin_words(n):
        ....:     return IntegerListsLex(length=n+1, min_slope=-1, max_slope=1,
        ....:                ceiling=[0]+[+oo for i in range(n-1)]+[0])
        sage: motzkin_words(4).list()
        [[0, 1, 2, 1, 0],
         [0, 1, 1, 1, 0],
         [0, 1, 1, 0, 0],
         [0, 1, 0, 1, 0],
         [0, 1, 0, 0, 0],
         [0, 0, 1, 1, 0],
         [0, 0, 1, 0, 0],
         [0, 0, 0, 1, 0],
         [0, 0, 0, 0, 0]]
        sage: [motzkin_words(n).cardinality() for n in range(8)]
        [1, 1, 2, 4, 9, 21, 51, 127]
        sage: oeis(_)                  # optional -- internet
        0: A001006: Motzkin numbers: number of ways of drawing any number
        of nonintersecting chords joining n (labeled) points on a circle.

    or dyck words (see also :class:`DyckWords`), through the bijection
    with paths from (0,0) to (n,n) with left and up steps that remain
    below the diagonal::

        sage: def dyck_words(n):
        ....:     return IntegerListsLex(length=n, ceiling=range(n+1), min_slope=0)
        sage: [dyck_words(n).cardinality() for n in range(8)]
        [1, 1, 2, 5, 14, 42, 132, 429]
        sage: dyck_words(3).list()
        [[0, 1, 2], [0, 1, 1], [0, 0, 2], [0, 0, 1], [0, 0, 0]]



    .. RUBRIC:: Situations with improper lexicographic enumeration

    The set of all lists of integers cannot be enumerated in inverse
    lexicographic order, since there is no largest list (take `[n]`
    for `n` as large as desired)::

        sage: IntegerListsLex().first()
        Traceback (most recent call last):
        ...
        ValueError: The specified parameters do not allow for an
        inverse lexicographic iterator

    Here is a variant which could be enumerated in inverse lexicographically
    increasing order but not in inverse lexicographically decreasing order::

        sage: L = IntegerListsLex(length=2, ceiling=[Infinity, 0], floor=[0,1])
        sage: for l in L: print l
        Traceback (most recent call last):
        ...
        ValueError: infinite upper bound for values of m

    Even when the sum is specified, it is not necessarily possible to
    enumerate all elements inverse lexicographically. In the following
    example, the list `[1, 1, 1]` will never appear in the enumeration::

        sage: IntegerListsLex(3).first()
        Traceback (most recent call last):
        ...
        ValueError: The specified parameters do not allow for an
        inverse lexicographic iterator

    If one wants to proceed anyway, one can sign a waiver by setting
    ``check=False``::

        sage: L = IntegerListsLex(3, check=False)
        sage: it = iter(L)
        sage: [it.next() for i in range(6)]
        [[3], [2, 1], [2, 0, 1], [2, 0, 0, 1], [2, 0, 0, 0, 1], [2, 0, 0, 0, 0, 1]]

    .. RUBRIC:: On trailing zeroes, and their caveats

    As mentioned above, when several lists satisfying the constraints
    differ only by trailing zeroes, only the shortest one is listed::

        sage: L = IntegerListsLex(max_length=4, max_part=1)
        sage: L.list()
        [[1, 1, 1, 1],
         [1, 1, 1],
         [1, 1, 0, 1],
         [1, 1],
         [1, 0, 1, 1],
         [1, 0, 1],
         [1, 0, 0, 1],
         [1],
         [0, 1, 1, 1],
         [0, 1, 1],
         [0, 1, 0, 1],
         [0, 1],
         [0, 0, 1, 1],
         [0, 0, 1],
         [0, 0, 0, 1],
         []]

    and counted::

        sage: L.cardinality()
        16

    Still, the others are considered as elements of `L`::

        sage: L = IntegerListsLex(4,min_length=3,max_length=4)
        sage: L.list()
        [..., [2, 2, 0], ...]

        sage: [2, 2, 0] in L       # in L.list()
        True
        sage: [2, 2, 0, 0] in L    # not in L.list() !
        True
        sage: [2, 2, 0, 0, 0] in L
        False

    .. RUBRIC:: Specifying functions as input for the floor or ceiling

    We construct all lists of sum `4` and length `4` such that ``l[i] <= i``::

        sage: list(IntegerListsLex(4, length=4, ceiling=lambda i: i, check=False))
        [[0, 1, 2, 1], [0, 1, 1, 2], [0, 1, 0, 3], [0, 0, 2, 2], [0, 0, 1, 3]]

    .. WARNING::

        When passing a function as ``floor`` or ``ceiling``, it may
        become undecidable to detect improper inverse lexicographic
        enumeration. For example, the following example has a finite
        enumeration::

            sage: L = IntegerListsLex(3, floor=lambda i: 1 if i>=2 else 0, check=False)
            sage: L.list()
            [[3],
             [2, 1],
             [2, 0, 1],
             [1, 2],
             [1, 1, 1],
             [1, 0, 2],
             [1, 0, 1, 1],
             [0, 3],
             [0, 2, 1],
             [0, 1, 2],
             [0, 1, 1, 1],
             [0, 0, 3],
             [0, 0, 2, 1],
             [0, 0, 1, 2],
             [0, 0, 1, 1, 1]]

        but one cannot decide whether the following has an improper
        inverse lexicographic enumeration without computing the floor
        all the way to ``Infinity``::

            sage: L = IntegerListsLex(3, floor=lambda i: 0, check=False)
            sage: it = iter(L)
            sage: [it.next() for i in range(6)]
            [[3], [2, 1], [2, 0, 1], [2, 0, 0, 1], [2, 0, 0, 0, 1], [2, 0, 0, 0, 0, 1]]

        Hence a warning is raised when a function is specified as
        input, unless the waiver is signed by setting ``check=False``::

            sage: L = IntegerListsLex(3, floor=lambda i: 1 if i>=2 else 0)
            doctest:...
            A function has been given as input of the floor=[...] or ceiling=[...]
            arguments of IntegerListsLex. Please see the documentation for the caveats.
            If you know what you are doing, you can set check=False to skip this warning.

        Similarly, the algorithm may need to search forever for a
        solution when the ceiling is ultimately zero::

            sage: L = IntegerListsLex(2,ceiling=lambda i:0, check=False)
            sage: L.first()           # not tested: will hang forever
            sage: L = IntegerListsLex(2,ceiling=lambda i:0 if i<20 else 1, check=False)
            sage: it = iter(L)
            sage: it.next()
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1]
            sage: it.next()
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1]
            sage: it.next()
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1]


    .. RUBRIC:: Specifying how to construct elements

    This is the list of all monomials of degree `4` which divide the
    monomial `x^3y^1z^2` (a monomial being identified with its
    exponent vector)::

        sage: R.<x,y,z> = QQ[]
        sage: m = [3,1,2]
        sage: def term(exponents):
        ....:     return x^exponents[0] * y^exponents[1] * z^exponents[2]
        sage: list( IntegerListsLex(4, length=len(m), ceiling=m, element_constructor=term) )
        [x^3*y, x^3*z, x^2*y*z, x^2*z^2, x*y*z^2]

    Note the use of the ``element_constructor`` option to specify how
    to construct elements from a plain list.


    .. RUBRIC:: Input list or iterable for the sum

    One may pass a list or iterable `L` as input for the sum. In this
    case, the elements will be generated inverse lexicographically,
    for each sum in `L` in turn::

        sage: C = IntegerListsLex([0,1,2], length=2)
        sage: C.list()
        [[0, 0],  [1, 0], [0, 1],  [2, 0], [1, 1], [0, 2]]

    This is in fact just a short hand for using
    :class:`DisjointUnionEnumeratedSets`::

        sage: C
        Disjoint union of Finite family
        {0: Integer lists of sum 0 satisfying certain constraints,
         1: Integer lists of sum 1 satisfying certain constraints,
         2: Integer lists of sum 2 satisfying certain constraints}

    This feature is mostly here for backward compatibility, as using
    :class:`DisjointEnumeratedSets` is more general and flexible::

        sage: C = DisjointUnionEnumeratedSets(Family([0,1,2],
        ....:         lambda n: IntegerListsLex(n, length=2)))
        sage: C.list()
        [[0, 0],  [1, 0], [0, 1],  [2, 0], [1, 1], [0, 2]]

    ALGORITHM:

    The iteration algorithm uses a depth first search through the
    prefix tree of the list of integers (see also
    :ref:`section-generic-integerlistlex`). While doing so, it does
    some lookahead heuristics to attempt to cut dead branches.

    In most practical use cases, most dead branches are cut. Then,
    roughly speaking, the time needed to iterate through all the
    elements of `S` is proportional to the number of elements, where
    the proportion factor is controlled by the length `l` of the
    longest element of `S`. In addition, the memory usage is also
    controlled by `l`, which is to say negligible in practice.

    Still, there remains much room for efficiency improvements; see
    :trac:`18055`, :trac:`18056`.

    .. NOTE::

        The generation algorithm could in principle be extended to
        deal with non-constant slope constraints and with negative
        parts.

    TESTS:

    This example from the combinatorics tutorial used to fail before
    :trac:`17979` because the floor conditions did not satisfy the
    slope conditions::

        sage: I = IntegerListsLex(16, min_length=2, max_slope=-1, floor=[5,3,3])
        sage: I.list()
        [[13, 3], [12, 4], [11, 5], [10, 6], [9, 7], [9, 4, 3], [8, 5, 3], [8, 4, 3, 1],
         [7, 6, 3], [7, 5, 4], [7, 5, 3, 1], [7, 4, 3, 2], [6, 5, 4, 1], [6, 5, 3, 2],
         [6, 4, 3, 2, 1]]

    ::

        sage: Partitions(2, max_slope=-1, length=2).list()
        []
        sage: list(IntegerListsLex(0, floor=ConstantFunction(1), min_slope=0))
        [[]]
        sage: list(IntegerListsLex(0, floor=ConstantFunction(1), min_slope=0, max_slope=0))
        [[]]
        sage: list(IntegerListsLex(0, max_length=0, floor=ConstantFunction(1), min_slope=0, max_slope=0))
        [[]]
        sage: list(IntegerListsLex(0, max_length=0, floor=ConstantFunction(0), min_slope=0, max_slope=0))
        [[]]
        sage: list(IntegerListsLex(0, min_part=1, min_slope=0))
        [[]]
        sage: list(IntegerListsLex(1, min_part=1, min_slope=0))
        [[1]]
        sage: list(IntegerListsLex(0, min_length=1, min_part=1, min_slope=0))
        []
        sage: list(IntegerListsLex(0, min_length=1, min_slope=0))
        [[0]]
        sage: list(IntegerListsLex(3, max_length=2, ))
        [[3], [2, 1], [1, 2], [0, 3]]
        sage: partitions = {"min_part": 1, "max_slope": 0}
        sage: partitions_min_2 = {"floor": ConstantFunction(2), "max_slope": 0}
        sage: compositions = {"min_part": 1}
        sage: integer_vectors = lambda l: {"length": l}
        sage: lower_monomials = lambda c: {"length": c, "floor": lambda i: c[i]}
        sage: upper_monomials = lambda c: {"length": c, "ceiling": lambda i: c[i]}
        sage: constraints = { "min_part":1, "min_slope": -1, "max_slope": 0}
        sage: list(IntegerListsLex(6, **partitions))
        [[6],
         [5, 1],
         [4, 2],
         [4, 1, 1],
         [3, 3],
         [3, 2, 1],
         [3, 1, 1, 1],
         [2, 2, 2],
         [2, 2, 1, 1],
         [2, 1, 1, 1, 1],
         [1, 1, 1, 1, 1, 1]]
        sage: list(IntegerListsLex(6, **constraints))
        [[6],
         [3, 3],
         [3, 2, 1],
         [2, 2, 2],
         [2, 2, 1, 1],
         [2, 1, 1, 1, 1],
         [1, 1, 1, 1, 1, 1]]
        sage: list(IntegerListsLex(1, **partitions_min_2))
        []
        sage: list(IntegerListsLex(2, **partitions_min_2))
        [[2]]
        sage: list(IntegerListsLex(3, **partitions_min_2))
        [[3]]
        sage: list(IntegerListsLex(4, **partitions_min_2))
        [[4], [2, 2]]
        sage: list(IntegerListsLex(5, **partitions_min_2))
        [[5], [3, 2]]
        sage: list(IntegerListsLex(6, **partitions_min_2))
        [[6], [4, 2], [3, 3], [2, 2, 2]]
        sage: list(IntegerListsLex(7, **partitions_min_2))
        [[7], [5, 2], [4, 3], [3, 2, 2]]
        sage: list(IntegerListsLex(9, **partitions_min_2))
        [[9], [7, 2], [6, 3], [5, 4], [5, 2, 2], [4, 3, 2], [3, 3, 3], [3, 2, 2, 2]]
        sage: list(IntegerListsLex(10, **partitions_min_2))
        [[10],
         [8, 2],
         [7, 3],
         [6, 4],
         [6, 2, 2],
         [5, 5],
         [5, 3, 2],
         [4, 4, 2],
         [4, 3, 3],
         [4, 2, 2, 2],
         [3, 3, 2, 2],
         [2, 2, 2, 2, 2]]
        sage: list(IntegerListsLex(4, **compositions))
        [[4], [3, 1], [2, 2], [2, 1, 1], [1, 3], [1, 2, 1], [1, 1, 2], [1, 1, 1, 1]]
        sage: list(IntegerListsLex(6, min_length=1, floor=[7]))
        []
        sage: L = IntegerListsLex(10**100,length=1)
        sage: L.list()
        [[10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000]]

    Noted on :trac:`17898`::

        sage: list(IntegerListsLex(4, min_part=1, length=3, min_slope=1))
        []
        sage: IntegerListsLex(6, ceiling=[4,2], floor=[3,3]).list()
        []
        sage: IntegerListsLex(6, min_part=1, max_part=3, max_slope=-4).list()
        []

    Noted in :trac:`17548`, which are now fixed::

        sage: IntegerListsLex(10, min_part=2, max_slope=-1).list()
        [[10], [8, 2], [7, 3], [6, 4], [5, 3, 2]]
        sage: IntegerListsLex(5, min_slope=1, floor=[2,1,1], max_part=2).list()
        []
        sage: IntegerListsLex(4, min_slope=0, max_slope=0).list()
        [[4], [2, 2], [1, 1, 1, 1]]
        sage: IntegerListsLex(6, min_slope=-1, max_slope=-1).list()
        [[6], [3, 2, 1]]
        sage: IntegerListsLex(6, min_length=3, max_length=2, min_part=1).list()
        []
        sage: I = IntegerListsLex(3, max_length=2, min_part=1)
        sage: I.list()
        [[3], [2, 1], [1, 2]]
        sage: [1,1,1] in I
        False
        sage: I=IntegerListsLex(10, ceiling=[4], max_length=1, min_part=1)
        sage: I.list()
        []
        sage: [4,6] in I
        False
        sage: I = IntegerListsLex(4, min_slope=1, min_part=1, max_part=2)
        sage: I.list()
        []
        sage: I = IntegerListsLex(7, min_slope=1, min_part=1, max_part=4)
        sage: I.list()
        [[3, 4], [1, 2, 4]]
        sage: I = IntegerListsLex(4, floor=[2,1], ceiling=[2,2], max_length=2, min_slope=0)
        sage: I.list()
        [[2, 2]]
        sage: I = IntegerListsLex(10, min_part=1, max_slope=-1)
        sage: I.list()
        [[10], [9, 1], [8, 2], [7, 3], [7, 2, 1], [6, 4], [6, 3, 1], [5, 4, 1],
         [5, 3, 2], [4, 3, 2, 1]]


    :trac:`17979`, comment 191::

        sage: list(IntegerListsLex(1, min_length=2, min_slope=0, max_slope=0))
        []

    Internally, the iterator works on a single list that is mutated
    along the way. The following test makes sure that we actually make a copy of
    this list before passing it to ``element_constructor`` in order to
    avoid reference effects::

        sage: from sage.misc.c3_controlled import identity
        sage: P = IntegerListsLex(n=3, max_slope=0, min_part=1, element_constructor=identity)
        sage: list(P)
        [[3], [2, 1], [1, 1, 1]]

    Same, step by step::

        sage: it = iter(P)
        sage: a = it.next(); a
        [3]
        sage: b = it.next(); b
        [2, 1]
        sage: a
        [3]
        sage: a is b
        False

    Tests from `MuPAD-Combinat <http://mupad-combinat.svn.sourceforge.net/viewvc/mupad-combinat/trunk/MuPAD-Combinat/lib/COMBINAT/TEST/MachineIntegerListsLex.tst>`_::

        sage: IntegerListsLex(7, min_length=2, max_length=6, floor=[0,0,2,0,0,1], ceiling=[3,2,3,2,1,2]).cardinality()
        83
        sage: IntegerListsLex(7, min_length=2, max_length=6, floor=[0,0,2,0,1,1], ceiling=[3,2,3,2,1,2]).cardinality()
        53
        sage: IntegerListsLex(5, min_length=2, max_length=6, floor=[0,0,2,0,0,0], ceiling=[2,2,2,2,2,2]).cardinality()
        30
        sage: IntegerListsLex(5, min_length=2, max_length=6, floor=[0,0,1,1,0,0], ceiling=[2,2,2,2,2,2]).cardinality()
        43

        sage: IntegerListsLex(0, min_length=0, max_length=7, floor=[1,1,0,0,1,0], ceiling=[4,3,2,3,2,2,1]).first()
        []

        sage: IntegerListsLex(0, min_length=1, max_length=7, floor=[0,1,0,0,1,0], ceiling=[4,3,2,3,2,2,1]).first()
        [0]
        sage: IntegerListsLex(0, min_length=1, max_length=7, floor=[1,1,0,0,1,0], ceiling=[4,3,2,3,2,2,1]).cardinality()
        0

        sage: IntegerListsLex(2, min_length=0, max_length=7, floor=[1,1,0,0,0,0], ceiling=[4,3,2,3,2,2,1]).first()  # Was [1,1], due to slightly different specs
        [2]
        sage: IntegerListsLex(1, min_length=1, max_length=7, floor=[1,1,0,0,0,0], ceiling=[4,3,2,3,2,2,1]).first()
        [1]
        sage: IntegerListsLex(1, min_length=2, max_length=7, floor=[1,1,0,0,0,0], ceiling=[4,3,2,3,2,2,1]).cardinality()
        0
        sage: IntegerListsLex(2, min_length=5, max_length=7, floor=[1,1,0,0,0,0], ceiling=[4,3,2,3,2,2,1]).first()
        [1, 1, 0, 0, 0]
        sage: IntegerListsLex(2, min_length=5, max_length=7, floor=[1,1,0,0,0,1], ceiling=[4,3,2,3,2,2,1]).first()
        [1, 1, 0, 0, 0]
        sage: IntegerListsLex(2, min_length=5, max_length=7, floor=[1,1,0,0,1,0], ceiling=[4,3,2,3,2,2,1]).cardinality()
        0

        sage: IntegerListsLex(4, min_length=3, max_length=6, floor=[2, 1, 2, 1, 1, 1], ceiling=[3, 1, 2, 3, 2, 2]).cardinality()
        0
        sage: IntegerListsLex(5, min_length=3, max_length=6, floor=[2, 1, 2, 1, 1, 1], ceiling=[3, 1, 2, 3, 2, 2]).first()
        [2, 1, 2]
        sage: IntegerListsLex(6, min_length=3, max_length=6, floor=[2, 1, 2, 1, 1, 1], ceiling=[3, 1, 2, 3, 2, 2]).first()
        [3, 1, 2]
        sage: IntegerListsLex(12, min_length=3, max_length=6, floor=[2, 1, 2, 1, 1, 1], ceiling=[3, 1, 2, 3, 2, 2]).first()
        [3, 1, 2, 3, 2, 1]
        sage: IntegerListsLex(13, min_length=3, max_length=6, floor=[2, 1, 2, 1, 1, 1], ceiling=[3, 1, 2, 3, 2, 2]).first()
        [3, 1, 2, 3, 2, 2]
        sage: IntegerListsLex(14, min_length=3, max_length=6, floor=[2, 1, 2, 1, 1, 1], ceiling=[3, 1, 2, 3, 2, 2]).cardinality()
        0

    Indirect tests about specifying an iterable for `n`::

        sage: P = Partitions(NonNegativeIntegers(), max_part = 3)
        sage: P.first()
        []
        sage: P = Partitions(NonNegativeIntegers())
        sage: P.first()
        []
    """
    __metaclass__ = ClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, n=None, **kwargs):
        r"""
        Return a disjoint union if ``n`` is a list or iterable.

        EXAMPLES::

            sage: IntegerListsLex(NN, max_length=3)
            Disjoint union of Lazy family (<lambda>(i))_{i in Non negative integer semiring}
        """
        if isinstance(n, (list,collections.Iterable)):
            return DisjointUnionEnumeratedSets(Family(n, lambda i: IntegerListsLex(i, **kwargs)))
        else:
            return typecall(cls, n=n, **kwargs)

    def __init__(self,
                 n=None,
                 length=None, min_length=0, max_length=Infinity,
                 floor=None, ceiling=None,
                 min_part=0, max_part=Infinity,
                 min_slope=-Infinity, max_slope=Infinity,
                 min_sum=0, max_sum=Infinity,
                 name=None,
                 category=None,
                 element_constructor=None, element_class=None,
                 global_options=None,
                 check=True):
        """
        Initialize ``self``.

        TESTS::

            sage: C = IntegerListsLex(2, length=3)
            sage: C == loads(dumps(C))
            True
            sage: C == loads(dumps(C)) # this did fail at some point, really!
            True
            sage: C is loads(dumps(C)) # todo: not implemented
            True
            sage: C.cardinality().parent() is ZZ
            True
            sage: TestSuite(C).run()
        """
        ## TODO handle finiteness conditions Re: category, when possible, warn when (check parameter)

        if category is None:
            category = EnumeratedSets().Finite()

        # self._warning will be set to ``True`` if a function is given
        # as input for floor or ceiling; in this case a warning will
        # be emitted, unless the user signs the waiver. See the
        # documentation.
        self._warning = False # warning for dangerous (but possibly valid) usage
        self._check = check

        if n is not None:
            n = ZZ(n)
            self._min_sum = n
            self._max_sum = n
        else:
            self._min_sum = min_sum
            self._max_sum = max_sum

        if length is not None:
            length = ZZ(length)
            self._min_length = length
            self._max_length = length
        else:
            min_length = ZZ(min_length)
            if min_length < 0:
                min_length = 0
            self._min_length = min_length
            if max_length != Infinity:
                max_length = ZZ(max_length)
            self._max_length = max_length

        self._min_slope = min_slope
        self._max_slope = max_slope

        min_part = ZZ(min_part)
        if min_part < 0:
            raise NotImplementedError("strictly negative min_part")

        if max_part != Infinity:
            max_part = ZZ(max_part)

        if floor is None:
            self._floor = ConstantFunction(min_part)
            self._floor_type = "constant"
            self._floor_limit = min_part
            self._floor_limit_start = 0
        elif isinstance(floor, (list, tuple)):
            if not all(i in ZZ for i in floor):
                raise TypeError("the parts of floor={} should be non negative integers".format(floor))
            if not all(i >= 0 for i in floor):
                raise NotImplementedError("negative parts in floor={}".format(floor))
            if min_part > 0:
                floor = map(lambda i: max(i, min_part), floor)
            floor = lower_envelope_from_list(floor, min_slope, max_slope)
            self._floor = list_function(floor, min_part)
            self._floor_type = "list"
            self._floor_limit = min_part
            self._floor_limit_start = len(floor)
        elif callable(floor):
            self._warning = True
            if min_part > 0:
                self._floor = lambda i: max(min_part, floor(i))
            else:
                self._floor = floor
            self._floor_type = "function"
            self._floor_limit = None
            self._floor_limit_start = Infinity
        else:
            raise TypeError("floor should be a list, tuple, or function")

        if ceiling is None:
            self._ceiling = ConstantFunction(max_part)
            self._ceiling_type = "constant"
            self._ceiling_limit = max_part
            self._ceiling_limit_start = 0
        elif isinstance(ceiling, (list, tuple)):
            if not all(i==Infinity or i in ZZ for i in ceiling):
                raise TypeError("the parts of ceiling={} should be non negative integers".format(ceiling))
            if not all(i >= 0 for i in ceiling):
                raise NotImplementedError("negative parts in floor={}".format(ceiling))
            if max_part < Infinity:
                ceiling = map(lambda i: min(i, max_part), ceiling)
            ceiling = upper_envelope_from_list(ceiling, min_slope, max_slope)
            self._ceiling = list_function(ceiling, max_part)
            self._ceiling_type = "list"
            self._ceiling_limit = max_part
            self._ceiling_limit_start = len(ceiling)
        elif callable(ceiling):
            self._warning = True
            if max_part < Infinity:
                self._ceiling = lambda i: min(max_part, ceiling(i))
            else:
                self._ceiling = ceiling
            self._ceiling = ceiling
            self._ceiling_type = "function"
            self._ceiling_limit = None
            self._ceiling_limit_start = Infinity
        else:
            raise ValueError("Unable to parse value of parameter ceiling")

        if name is not None:
            self.rename(name)

        if self._warning and self._check:
            warn("""
A function has been given as input of the floor=[...] or ceiling=[...]
arguments of IntegerListsLex. Please see the documentation for the caveats.
If you know what you are doing, you can set check=False to skip this warning.""")

        # In case we want output to be of a different type,
        if element_constructor is not None:
            self._element_constructor_ = element_constructor
        if element_class is not None:
            self.Element = element_class
        if global_options is not None:
            self.global_options = global_options

        Parent.__init__(self, category=category)

    @cached_method
    def _check_lexicographic_iterable(self):
        """
        Checks whether the parameters give a proper inverse lexicographic iterator.

        EXAMPLES::

            sage: IntegerListsLex(4).list()
            Traceback (most recent call last):
            ...
            ValueError: The specified parameters do not allow for an
            inverse lexicographic iterator

            sage: it = iter(IntegerListsLex(4, check=False))
            sage: for _ in range(20): print next(it)
            [4]
            [3, 1]
            [3, 0, 1]
            [3, 0, 0, 1]
            [3, 0, 0, 0, 1]
            [3, 0, 0, 0, 0, 1]
            [3, 0, 0, 0, 0, 0, 1]
            [3, 0, 0, 0, 0, 0, 0, 1]
            [3, 0, 0, 0, 0, 0, 0, 0, 1]
            [3, 0, 0, 0, 0, 0, 0, 0, 0, 1]
            [3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]
            [3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]
            [3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]
            [3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]
            [3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]
            [3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]
            [3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]
            [3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]
            [3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]
            [3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]

            sage: L = IntegerListsLex(ceiling=[0], min_slope=1, max_slope=2)
            sage: L.list()
            Traceback (most recent call last):
            ...
            ValueError: The specified parameters do not allow for an
            inverse lexicographic iterator

            sage: L = IntegerListsLex(ceiling=[0], min_slope=1, max_slope=1)
            sage: L.list()
            Traceback (most recent call last):
            ...
            ValueError: The specified parameters do not allow for an
            inverse lexicographic iterator

        The next example shows a case that is finite since we remove trailing zeroes::

            sage: list(IntegerListsLex(ceiling=[0], max_slope=0))
            [[]]
            sage: L = IntegerListsLex(ceiling=[1], min_slope=1, max_slope=1)
            sage: L.list()
            Traceback (most recent call last):
            ...
            ValueError: The specified parameters do not allow for an
            inverse lexicographic iterator

        In the next examples, there is either no solution, or the region
        is bounded::

            sage: IntegerListsLex(min_sum=10, max_sum=5).list()
            []
            sage: IntegerListsLex(max_part=1, min_slope=10).list()
            [[1], []]
            sage: IntegerListsLex(max_part=100, min_slope=10).first()
            [100]
            sage: I = IntegerListsLex(ceiling=[1,Infinity], max_part=2, min_slope=1)
            sage: I.list()
            [[1, 2], [1], [0, 2], [0, 1, 2], [0, 1], []]
        """
        if self._warning or not self._check:
            return
        message = "The specified parameters do not allow for an inverse lexicographic iterator"
        s = sum(self._floor(i) for i in range(self._floor_limit_start))
        if self._max_sum < Infinity and self._max_length == Infinity and self._floor_limit == 0:
            if self._min_slope < 0 and self._max_slope > 0 and s < self._min_sum and self._min_sum <= self._max_sum:
                raise ValueError(message)
            if self._min_slope == 0 and s==0 and self._max_slope > 0:
                if self._max_sum>0: # this is assuming that we remove trailing zeroes
                    raise ValueError(message)
        elif self._max_sum == Infinity and self._max_length == Infinity:
            if self._max_slope == 0 and min(self._ceiling(i) for i in range(self._ceiling_limit_start+1)) == 0:
                return
            elif self._min_slope > 0 and self._ceiling_limit < Infinity:
                return
            elif self._max_slope >= 0 and self._ceiling_limit > 0:
                raise ValueError(message)

    def __cmp__(self, x):
        """
        Compares two different :class:`IntegerListsLex`.

        TODO: Fix this to make it more robust!

        For now, the comparison is done just on their repr's which is
        not robust!

        EXAMPLES::

            sage: C = IntegerListsLex(2, length=3)
            sage: D = IntegerListsLex(4, length=3)
            sage: repr(C) == repr(D)
            False
            sage: C == D
            False
        """
        return cmp(repr(self), repr(x))

    def _repr_(self):
        """
        Return the name of this combinatorial class.

        EXAMPLES::

            sage: C = IntegerListsLex(2, length=3)
            sage: C # indirect doctest
            Integer lists of sum 2 satisfying certain constraints

            sage: C = IntegerListsLex(2, length=3, name="A given name")
            sage: C
            A given name
        """
        if self._min_sum == self._max_sum:
            return "Integer lists of sum {} satisfying certain constraints".format(self._min_sum)
        elif self._max_sum == Infinity:
            if self._min_sum == 0:
                return "Integer lists with arbitrary sum satisfying certain constraints"
            else:
                return "Integer lists of sum at least {} satisfying certain constraints".format(self._min_sum)
        else:
            return "Integer lists of sum between {} and {} satisfying certain constraints".format(self._min_sum,self._max_sum)

    def __contains__(self, comp):
        """
        Return ``True`` if ``comp`` meets the constraints imposed by the arguments.

        EXAMPLES::

            sage: C = IntegerListsLex(n=2, max_length=3, min_slope=0)
            sage: all([l in C for l in C])
            True
        """
        if len(comp) < self._min_length or len(comp) > self._max_length:
            return False
        n = sum(comp)
        if n < self._min_sum or n > self._max_sum:
            return False
        for i in range(len(comp)):
            if comp[i] < self._floor(i):
                return False
            if comp[i] > self._ceiling(i):
                return False
        for i in range(len(comp)-1):
            slope = comp[i+1] - comp[i]
            if slope < self._min_slope or slope > self._max_slope:
                return False
        return True


    def __iter__(self):
        """
        Return an iterator for the elements of ``self``.

        EXAMPLES::

            sage: C = IntegerListsLex(2, length=3)
            sage: list(C)     # indirect doctest
            [[2, 0, 0], [1, 1, 0], [1, 0, 1], [0, 2, 0], [0, 1, 1], [0, 0, 2]]
        """
        return self._Iter(self)

    class _Iter:
        """
        Iterator class for IntegerListsLex

        The iterator is based on a depth-first search forest. If `I`
        is the iterator, the current state in the forest is given by
        ``I._current_list``. The range for the next entry (which
        corresponds to the next depth in the forest) is stored in
        ``I._search_ranges``. ``I._j`` stores the index of the last
        element of ``I._current_list``, whereas ``I._current_sum`` is
        the sum over all element of ``I._current_list``.
        """
        def __init__(self, parent):
            """
            TESTS::

                sage: C = IntegerListsLex(2, length=3)
                sage: I = IntegerListsLex._Iter(C)
                sage: I._search_ranges
                [(0, 2)]
                sage: I._current_list
                [3]
                sage: I._j
                0
                sage: I.finished
                False
            """
            self.parent = parent

            parent._check_lexicographic_iterable()

            self._search_ranges = []
            self._current_list = []
            self._j = -1     # index of last element of _current_list
            self._current_sum = 0     # sum of parts in _current_list

            self.finished = False

            # initialize for beginning of iteration
            self._push_search()

        def _push_search(self):
            """
            Push search forward. Resetting attributes.

            EXAMPLES::

                sage: C = IntegerListsLex(2, length=3)
                sage: I = C.__iter__()
                sage: I._j
                0
                sage: I._search_ranges
                [(0, 2)]
                sage: I._current_list
                [3]
                sage: I._current_sum
                3
                sage: I._push_search()
                sage: I._j
                1
                sage: I._search_ranges
                [(0, 2), (0, -1)]
                sage: I._current_list
                [3, 0]
                sage: I._current_sum
                3
            """
            if self._j >= 0:
                prev = self._current_list[self._j]
            else:
                prev = None
            self._j += 1
            interval = self._m_interval(self._j, self.parent._max_sum - self._current_sum, prev)
            val = interval[1] + 1 # iterator decrements before acting
            self._search_ranges.append(interval)
            self._current_list.append(val)
            self._current_sum += val

        def _pop_search(self):
            """
            Go back in search tree. Resetting attributes.

            EXAMPLES::

                sage: C = IntegerListsLex(2, length=3)
                sage: I = C.__iter__()
                sage: I._j
                0
                sage: I._search_ranges
                [(0, 2)]
                sage: I._current_sum
                3
                sage: I._current_list
                [3]
                sage: I._pop_search()
                sage: I._j
                -1
                sage: I._search_ranges
                []
                sage: I._current_sum
                0
                sage: I._current_list
                []
            """
            if self._j >= 0:
                self._j -= 1
                self._search_ranges.pop()
                self._current_sum -= self._current_list[-1]
                self._current_list.pop()

        def next(self):
            """
            Return the next element in the iteration.

            EXAMPLES::

                sage: C = IntegerListsLex(2, length=3)
                sage: I = IntegerListsLex._Iter(C)
                sage: I.next()
                [2, 0, 0]
                sage: I.next()
                [1, 1, 0]
            """
            if self.finished:
                raise StopIteration()

            rho = self._search_ranges
            mu = self._current_list
            p = self.parent
            min_sum = p._min_sum
            max_sum = p._max_sum
            min_length = p._min_length
            max_length = p._max_length

            while self._j >= 0: # j = -1 means that we have finished the bottom iteration

                # choose a new value for m to test

                mu[self._j] -= 1
                # m = mu[self._j]
                self._current_sum -= 1

                # check if the new value is valid, if not, pop to prefix, and now check if this is a solution

                if mu[self._j] < rho[self._j][0] or (self._j == max_length-1 and self._current_sum < min_sum):
                    self._pop_search()
                    if self._internal_list_valid():
                        return p._element_constructor(list(mu))
                    else:
                        continue

                # m = mu[self._j] is new and in range:
                # If we're at a leaf node, check if it is a solution to return, pop the layer from the stack.
                # Otherwise, check if any solutions are possible with this value of m.
                #
                # Possible cases to detect leaf nodes:
                # 1. List is of maximal length
                # 2. List is of maximal sum, and of at least minimal length (allow padding zeros)

                if (self._current_sum == max_sum and self._j >= min_length - 1) or self._j == max_length - 1:
                    if self._internal_list_valid():
                        return p._element_constructor(list(mu))
                elif self._possible_m(mu[self._j], self._j,
                                      min_sum - (self._current_sum-mu[self._j]),
                                      max_sum - (self._current_sum-mu[self._j])):
                    self._push_search()

            self.finished = True
            raise StopIteration()

        def _internal_list_valid(self):
            """
            Return whether the current list in the iteration variable ``self._current_list`` is a valid list.

            This method checks whether the sum of the parts in ``self._current_list``
            is in the right range, whether its length is in the
            required range, and whether there are trailing zeroes.  It does not check all of the
            necessary conditions to verify that an arbitrary list satisfies the
            constraints from the corresponding ``IntegerListsLex`` object, and should
            not be used except internally in the iterator class.

            EXAMPLES::

                sage: C = IntegerListsLex(2, length=3)
                sage: I = IntegerListsLex._Iter(C)
                sage: I._current_list
                [3]
                sage: I._internal_list_valid()
                False
                sage: I.next()
                [2, 0, 0]
                sage: I._current_list
                [2, 0, 0]
                sage: I._internal_list_valid()
                True
            """
            p = self.parent
            mu = self._current_list
            nu = self._current_sum
            l = self._j + 1
            good_sum = (nu >= p._min_sum and nu <= p._max_sum)
            good_length = (l >= p._min_length and l <= p._max_length)
            no_trailing_zeros = (l <= max(p._min_length,0) or mu[-1] != 0)
            return good_sum and good_length and no_trailing_zeros

        def _upper_envelope(self, m, j):
            """
            Return the upper envelope starting with value ``m`` at position ``j``.

            INPUT:

            - ``m`` -- a nonnegative integer (starting value)

            - ``j`` -- a nonnegative integer (position)

            This method returns a function of ``i`` which computes the
            upper envelope if the starting value is ``m`` at position
            ``j``. The upper envelope is the minimum of the ceiling
            function and the value restriction given by the slope
            conditions.

            EXAMPLES::

                sage: C = IntegerListsLex(2, length=3)
                sage: I = IntegerListsLex._Iter(C)
                sage: f = I._upper_envelope(1,1)
                sage: type(f)
                <type 'sage.misc.constant_function.ConstantFunction'>
                sage: f(1)
                inf
                sage: f(2)
                inf
                sage: C = IntegerListsLex(6, max_slope=1, max_part=3, max_length=6)
                sage: I = IntegerListsLex._Iter(C)
                sage: f = I._upper_envelope(1,1)
                sage: f(1)
                1
                sage: f(2)
                2
                sage: f(3)
                3
                sage: f(4)
                3
            """
            if self.parent._max_slope == Infinity:
                return self.parent._ceiling
            m = m - j*self.parent._max_slope
            return lambda i: min(m + i*self.parent._max_slope, self.parent._ceiling(i) )

        def _lower_envelope(self, m, j):
            """
            Return the lower envelope starting with value ``m`` at position ``j``.

            INPUT:

            - ``m`` -- a nonnegative integer (starting value)

            - ``j`` -- a nonnegative integer (position)

            This returns a function of ``i`` which compute the lower
            envelope if the starting value is ``m`` at position
            ``j``. The lower envelope is the maximum of the floor
            function and the value restriction given by the slope
            conditions.

            EXAMPLES::

                sage: C = IntegerListsLex(2, length=3)
                sage: I = IntegerListsLex._Iter(C)
                sage: f = I._lower_envelope(1,1)
                sage: type(f)
                <type 'sage.misc.constant_function.ConstantFunction'>
                sage: f(1)
                0
                sage: f(2)
                0
                sage: C = IntegerListsLex(6, min_slope=-1, min_part=1)
                sage: I = IntegerListsLex._Iter(C)
                sage: f = I._lower_envelope(3,1)
                sage: f(1)
                3
                sage: f(2)
                2
                sage: f(3)
                1
                sage: f(4)
                1
            """
            if self.parent._min_slope == -Infinity:
                return self.parent._floor
            m = m-j*self.parent._min_slope
            return lambda i: max( m + i*self.parent._min_slope, self.parent._floor(i) )

        def _m_interval(self, i, max_sum, prev=None):
            r"""
            Return coarse lower and upper bounds for the part ``m`` at position ``i``.

            INPUT:

            - ``i`` -- a nonnegative integer (position)

            - ``max_sum`` -- a nonnegative integer or ``+oo``

            - ``prev`` -- a nonnegative integer or ``None``

            Return coarse lower and upper bounds for the value ``m``
            of the part at position ``i`` so that there could exists
            some list suffix `v_i,\ldots,v_k` of sum bounded by
            ``max_sum`` and satisfying the floor and upper bound
            constraints. If ``prev`` is specified, then the slope
            conditions between ``v[i-1]=prev`` and ``v[i]=m`` should
            also be satisfied.

            This also raises an error if it can be detected that some
            part is neither directly nor indirectly bounded above,
            which implies that the constraints do not allow for an
            inverse lexicographic iterator.

            OUTPUT:

            A tuple of two integers ``(lower_bound, upper_bound)``.

            EXAMPLES::

                sage: C = IntegerListsLex(2, length=3)
                sage: I = IntegerListsLex._Iter(C)
                sage: I._m_interval(1,2)
                (0, 2)

            The second part is not bounded above, hence we can not
            iterate lexicographically through all the elements::

                sage: IntegerListsLex(ceiling=[2,infinity,3], max_length=3).first()
                Traceback (most recent call last):
                ...
                ValueError: infinite upper bound for values of m

            In the following examples, all parts are indirectly
            bounded above::

                sage: IntegerListsLex(ceiling=[2,infinity,2], max_length=3, min_slope=2).cardinality()
                1
                sage: IntegerListsLex(ceiling=[2,infinity,2], max_length=3, max_slope=3).cardinality()
                45

                sage: IntegerListsLex(max_part=2, max_length=3).cardinality()
                27
                sage: IntegerListsLex(3, max_length=3).cardinality()      # parts bounded by n
                10
                sage: IntegerListsLex(max_length=0, min_length=1).list()  # no part!
                []
                sage: IntegerListsLex(length=0).list()                    # no part!
                [[]]
            """
            p = self.parent

            lower_bound = max(0, p._floor(i))
            upper_bound = min(max_sum, p._ceiling(i))
            if prev != None:
                lower_bound = max(lower_bound, prev + p._min_slope)
                upper_bound = min(upper_bound, prev + p._max_slope)

            ## check for infinite upper bound, in case max_sum is infinite
            if p._check and upper_bound == Infinity:
                raise ValueError("infinite upper bound for values of m")

            return (lower_bound, upper_bound)

        def _possible_m(self, m, j, min_sum, max_sum):
            r"""
            Look ahead whether `m` is a possible choice for `v_j` at this stage.

            INPUT:

            - ``m`` -- a nonnegative integer (value)

            - ``j`` -- a nonnegative integer (position)

            - ``min_sum`` -- a nonnegative integer

            - ``max_sum`` -- a nonnegative integer or ``+oo``

            Tries to predict the existence of a list suffix
            `(v_j,...,v_k)` which:

            - satisfies the slope, ceiling, and length
              (``min_length<=j<=max_length``) constraints of ``self``;

            - starts with `v_j=m``;

            - has sum `v_j+\cdots+v_k` between ``min_sum`` and ``max_sum``.

            OUTPUT: ``False`` if it is guaranteed that no such list
            exists and ``True`` otherwise.

            EXAMPLES::

                sage: C = IntegerListsLex(2, length=3)
                sage: I = IntegerListsLex._Iter(C)
                sage: I._possible_m(1,1,2,2)
                True
                sage: I._possible_m(1,1,3,2)
                False

            ALGORITHM::

            The current algorithm computes, for `k=j,j+1,\ldots`, a
            lower bound `l_k` and an upper bound `u_k` for
            `v_j+\dots+v_k`, and stops if none of the invervals `[l_k,
            u_k]` intersect ``[min_sum, max_sum]``.

            The lower bound `l_k` is given by the area below the lower
            envelope between `i` and `k` and starting at `m`. The
            upper bound `u_k` is given by the area below the upper
            envelope between `i` and `k` and starting at `m`.

            The complexity of this algorithm is bounded above by
            ``O(max_length)``. When ``max_length=oo``, the algorithm
            is guaranteed to terminate, unless ``floor`` is a function
            which is eventually constant with value `0`, or which
            reaches the value `0` while ``max_slope=0``.

            Indeed, the lower bound `l_k` is increasing with `k`; in
            fact it is strictly increasing, unless the local lower bound
            at `k` is `0`. Furthermore as soon as ``l_k >= min_sum``,
            we can conclude; we can also conclude if we know that the
            floor is eventually constant with value `0`, or there is a
            local lower bound at `k` is `0` and ``max_slope=0``.

            .. RUBRIC:: Room for improvement

            Improved prediction: the lower bound `l_k` does not take
            into account the slope conditions, except for that imposed
            by the value `m` at `j`. Similarly for `u_k`.

            Improved speed: given that `l_k` is increasing with `k`,
            possibly some dichotomy could be used to search for `k`,
            with appropriate caching / fast calculation of the partial
            sums. Also, some of the information gained at depth `j`
            could be reused at depth `j+1`.

            TESTS::

                sage: C = IntegerListsLex(1, min_length=2, min_slope=0, max_slope=0)
                sage: I = IntegerListsLex._Iter(C)
                sage: I._possible_m(0,0,1,1)
                False
            """
            # Check code for various termination conditions.  Possible cases:
            # 0. interval [lower, upper] intersects interval [min_sum, max_sum] -- terminate True
            # 1. lower sum surpasses max_sum -- terminate False
            # 2. iteration surpasses max_length -- terminate False
            # 3. upper envelope is smaller than lower envelope -- terminate False
            # 4. max_slope <= 0 -- terminate False after upper passes 0
            # 5. ceiling_limit == 0 -- terminate False after reaching larger limit point
            # 6. (uncomputable) ceiling function == 0 for all but finitely many input values -- terminate False after reaching (unknown) limit point -- currently hangs

            if min_sum > max_sum:
                return False

            p = self.parent

            lower_envelope = self._lower_envelope(m,j)
            upper_envelope = self._upper_envelope(m,j)
            lower = 0    # The lower bound `l_k`
            upper = 0    # The upper bound `u_k`

            assert j >= 0
            # get to smallest valid number of parts
            for k in range(j, p._min_length-1):
                # We are looking at lists `v_j,...,v_k`
                lo = lower_envelope(k)
                up = upper_envelope(k)
                if lo > up:
                    return False
                lower += lo
                upper += up

            k = max(p._min_length-1,j)
            # Check if any of the intervals intersect the target interval
            while k < p._max_length:
                lo = lower_envelope(k)
                up = upper_envelope(k)
                if lo > up:
                    # There exists no valid list of length >= k
                    return False
                lower += lo
                upper += up
                assert lower <= upper

                if lower > max_sum:
                    # There cannot exist a valid list `v_j,\dots,v_l` with l>=k
                    return False

                if (p._max_slope <= 0 and up <= 0) or \
                   (p._ceiling_limit == 0 and k > p._ceiling_limit_start):
                    # This implies v_l=0 for l>=k: that is we would be generating
                    # a list with trailing zeroes
                    return False

                if min_sum <= upper and lower <= max_sum:
                    # There could exist a valid list `v_j,\dots,v_k`
                    return True

                k += 1

            return False

    class Element(ClonableArray):
        """
        Element class for :class:`IntegerListsLex`.
        """
        def check(self):
            """
            Check to make sure this is a valid element in its
            :class:`IntegerListsLex` parent.

            EXAMPLES::

                sage: C = IntegerListsLex(4)
                sage: C([4]).check()
                True
                sage: C([5]).check()
                False
            """
            return self.parent().__contains__(self)

##############################################################################
# Utilities to massage list inputs for ceiling and floor
##############################################################################

def list_function(l, default):
    r"""
    Generate a function on the nonnegative integers from a list of values

    INPUT:

    - `l` -- a list

    - `default` -- a default value to use for indices outside of the list

    OUTPUT:

    A function on the nonnegative integers which maps `i` to ``l[i]``
    if ``i<len(l)`` and to ``default`` otherwise.

    EXAMPLES::

        sage: from sage.combinat.integer_list import list_function
        sage: C = IntegerListsLex(2, length=3)
        sage: f = list_function([1,2], Infinity)
        sage: f(1)
        2
        sage: f(3)
        +Infinity
    """
    return lambda i: l[i] if i < len(l) else default

def lower_envelope_from_list(list, min_slope, max_slope):
    """
    Return the lowest list above ``list`` that satisfies the slope constraints.

    EXAMPLES::

        sage: import sage.combinat.integer_list as integer_list
        sage: integer_list.lower_envelope_from_list([4,2,6],-1,1)
        [4, 5, 6]
        sage: integer_list.lower_envelope_from_list([4,2,6],-2, 1)
        [4, 5, 6]
        sage: integer_list.lower_envelope_from_list([4,2,6,3,7],-2, 1)
        [4, 5, 6, 6, 7]
        sage: integer_list.lower_envelope_from_list([4,2,6,1], -2, 1)
        [4, 5, 6, 4]
    """

    new_list = list[:]
    for i in range(1, len(new_list)):
        new_list[i] = max(new_list[i], new_list[i-1] + min_slope)

    for i in reversed(range(len(new_list)-1)):
        new_list[i] = max( new_list[i], new_list[i+1] - max_slope)

    return new_list

def upper_envelope_from_list(list, min_slope, max_slope):
    """
    Return the higest list below ``list`` that satisfies the slope constraints.

    EXAMPLES::

        sage: import sage.combinat.integer_list as integer_list
        sage: integer_list.upper_envelope_from_list([4,2,6], -1, 1)
        [3, 2, 3]
        sage: integer_list.upper_envelope_from_list([4,2,6], -1, infinity)
        [3, 2, 6]
        sage: integer_list.upper_envelope_from_list([1,4,2], -1, 1)
        [1, 2, 2]
        sage: integer_list.upper_envelope_from_list([4,2,6,3,7], -2, 1)
        [4, 2, 3, 3, 4]
        sage: integer_list.upper_envelope_from_list([4,2,infinity,3,7], -2, 1)
        [4, 2, 3, 3, 4]
        sage: integer_list.upper_envelope_from_list([1, infinity, 2], -1, 1)
        [1, 2, 2]
        sage: integer_list.upper_envelope_from_list([infinity, 4, 2], -1, 1)
        [4, 3, 2]
    """

    new_list = list[:]
    for i in range(1, len(new_list)):
        new_list[i] = min(new_list[i], new_list[i-1] + max_slope)

    for i in reversed(range(len(new_list)-1)):
        new_list[i] = min( new_list[i], new_list[i+1] - min_slope)

    return new_list
