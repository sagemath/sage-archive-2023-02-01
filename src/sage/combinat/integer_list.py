r"""
Tools for generating lists of integers in lexicographic order
"""
#*****************************************************************************
#       Copyright (C) 2015 Bryan Gillespie <Brg008@gmail.com>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
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
    Lists of non negative integers with constraints, in lexicographic order.

    An *integer list* is a list `l` of nonnegative integers, its
    *parts*. The *length* ``len(l)`` of `l` is the number of its
    parts. The *sum* `|l|` of `l` is the sum of its parts. The slope
    (at position `i`) is the difference ``l[i+1]-l[i]`` between two
    consecutive entries.

    This class allows to construct the set of all integer lists `l`
    satisfying specified bounds on the sum, the length, the slope, and
    the individual parts, enumerated lexicographically. The main
    purpose is to provide a generic iteration engine for all the
    enumerated sets like :class:`Partitions`, :class:`Compositions`,
    :class:`IntegerVectors`. It can also be used to generate many
    other combinatorial objects like Dyck paths, Motzkin paths, etc.

    Mathematically speaking, this is a special case of sets of
    integral points of a polytope (or union thereof, when the length
    is not fixed). The set of allowable has been specifically designed
    to enable iteration with a good time and memory complexity, in
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

    - ``waiver`` -- boolean (default: False): whether to suppress the
      warning raised when functions are given as input to ``floor`` or
      ``ceiling``

    .. NOTE::

       Two valid integer lists are considered equivalent if they only
       differ by trailing zeroes. In this case, only the list with the
       least number of trailing zeroes will be produced.

       .. TODO:: Do we really want to keep this "feature"?

    EXAMPLES:

    We create the enumerated set of all lists of length `3` and sum
    `2`::

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

        sage: [L.floor(i) for i in range(5)]
        [2, 1, 2, 1, 1]

    Similarly, the ``ceiling`` and ``max_part`` constraints can be
    combined::

        sage: L = IntegerListsLex(4, ceiling=[2,3,1], max_part=2, length=3)
        sage: L.list()
        [[2, 2, 0], [2, 1, 1], [1, 2, 1]]
        sage: [L.ceiling(i) for i in range(5)]
        [2, 2, 1, 2, 2]


    .. RUBRIC:: Situations with improper lexicographic enumeration

    The set of all lists of integers cannot be enumerated
    lexicographically, since there is no largest list (take `[n]` for
    `n` as large as desired)::

        sage: IntegerListsLex().first()
        Traceback (most recent call last):
        ...
        ValueError: infinite upper bound for values of m

    Here is a variant which could be enumerated in lexicographically
    increasing order but not in lexicographically decreasing order::

        sage: L = IntegerListsLex(length=2, ceiling=[Infinity, 0], floor=[0,1])
        sage: for l in L: print l
        Traceback (most recent call last):
        ...
        ValueError: infinite upper bound for values of m

    Even when the sum is specified, it is not necessarily possible to
    enumerate all elements lexicographically. In the following
    example, the list `[1, 1, 1]` will never appear in the enumeration::

        sage: IntegerListsLex(3).first()
        Traceback (most recent call last):
        ...
        ValueError: The specified parameters do not allow for a lexicographic iterator!

    .. TODO:: Maybe this should be ``check=False`` in this case?

    If one wants to proceed anyway, one can sign a waiver by setting
    ``waiver=True``::

        sage: L = IntegerListsLex(3, waiver=True)
        sage: it = iter(L)
        sage: [it.next() for i in range(6)]
        [[3], [2, 1], [2, 0, 1], [2, 0, 0, 1], [2, 0, 0, 0, 1], [2, 0, 0, 0, 0, 1]]


    .. RUBRIC:: Specifying functions as input for the floor or ceiling

    We construct all lists of sum `4` and length `4` such that ``l[i] <= i``::

        sage: list(IntegerListsLex(4, length=4, ceiling=lambda i: i, waiver=True))
        [[0, 1, 2, 1], [0, 1, 1, 2], [0, 1, 0, 3], [0, 0, 2, 2], [0, 0, 1, 3]]

    .. WARNING::

        When passing a function as ``floor`` or ``ceiling``, it may
        become undecidable to detect improper lexicographic
        enumeration. For example, the following example as a finite
        enumeration::

            sage: L = IntegerListsLex(3, floor=lambda i: 1 if i>=2 else 0, waiver=True)
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
        lexicographic enumeration without computing the floor all the
        way to Infinity::

            sage: L = IntegerListsLex(3, floor=lambda i: 0, waiver=True)
            sage: it = iter(L)
            sage: [it.next() for i in range(6)]
            [[3], [2, 1], [2, 0, 1], [2, 0, 0, 1], [2, 0, 0, 0, 1], [2, 0, 0, 0, 0, 1]]

        Hence a warning is raised when a function is specified as
        input, unless the waiver is signed::

            sage: L = IntegerListsLex(3, floor=lambda i: 1 if i>=2 else 0)
            doctest:...
            A function has been given as input of the floor=[...] or ceiling=[...]
            arguments of IntegerListsLex. Please see the documentation for the caveats.
            If you know what you are doing, you can set waiver=True to skip this warning.

        Similarly, the algorithm may need to search forever for a
        solution when the ceiling is ultimately zero::

            sage: L = IntegerListsLex(2,ceiling=lambda i:0, waiver=True)
            sage: L.first()           # not tested: will hang forever
            sage: L = IntegerListsLex(2,ceiling=lambda i:0 if i<20 else 1, waiver=True)
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


    .. RUBRIC:: list or iterable as input for the sum

    One may pass a list or iterable `L` as input for the sum. In this
    case, the elements will be generated lexicographically, for each
    sum in `L` in turn::

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
    :ref:`_section-generic-integerlistlex`). While doing so, it does
    some lookahead which allows for cutting most dead branches.

    The complexity of the algorithm has not been formally proven, but
    the average runtime for producing each list `l` is suspected to be
    bounded by a low-degree polynomial in ``lmax``, where ``lmax`` is
    the length of the longest list. Similarly, the space complexity of
    the algorithm is bounded by a low-degree polynomial in ``lmax``.

    .. NOTE::

        The generation algorithm could in principle be extended to
        deal with non-constant slope constraints and with negative
        parts.

    .. TODO::

        - Integrate all remaining tests from
          http://mupad-combinat.svn.sourceforge.net/viewvc/mupad-combinat/trunk/MuPAD-Combinat/lib/COMBINAT/TEST/MachineIntegerListsLex.tst

        - Integrate all tests from the ticket, sage-devel, ...

        - Improve the lookahead to get a proper complexity in the
          following example; even just for n=1000 this is taking a
          long time::

              sage: n = 100
              sage: IntegerListsLex(binomial(n,2)-1, length=n, min_slope=1).list()
              []

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
                 n=None, min_sum=0, max_sum=Infinity,
                 length=None, min_length=0, max_length=Infinity,
                 floor=None, ceiling=None,
                 min_part=0, max_part=Infinity,
                 min_slope=-Infinity, max_slope=Infinity,
                 name=None,
                 category=None,
                 element_constructor=None, element_class=None,
                 global_options=None,
                 waiver=False):
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
        ## TODO handle finiteness conditions Re: category, when possible, warn when not (waiver parameter)

        if category is None:
            category = EnumeratedSets().Finite()

        # self._warning will be set to ``True`` if a function is given
        # as input for floor or ceiling; in this case a warning will
        # be emitted, unless the user signs the waiver. See the
        # documentation.
        self._warning = False # warning for dangerous (but possibly valid) usage
        self._waiver = waiver

        if n is not None:
            if n in ZZ:
                self.min_sum = n
                self.max_sum = n
            else:
                raise TypeError("invalid value for parameter n")
        else:
            self.min_sum = min_sum
            self.max_sum = max_sum

        if length is not None:
            self.min_length = length
            self.max_length = length
        else:
            self.min_length = min_length
            self.max_length = max_length

        self.min_slope = min_slope
        self.max_slope = max_slope

        if min_part not in ZZ:
            raise TypeError("min_part (={}) should be an integer".format(min_part))
        elif min_part <0:
            raise NotImplementedError("strictly negative min_part")

        if max_part != Infinity and max_part not in ZZ:
            raise TypeError("max_part (={}) should be an integer or +oo".format(max_part))

        if floor is None:
            self.floor = ConstantFunction(min_part)
            self.floor_type = "constant"
            self.floor_limit = min_part
            self.floor_limit_start = 0
        elif isinstance(floor, (list, tuple)):
            if not all(i in ZZ for i in floor):
                raise TypeError("the parts of floor={} should be non negative integers".format(floor))
            if not all(i >= 0 for i in floor):
                raise NotImplementedError("negative parts in floor={}".format(floor))
            if min_part > 0:
                floor = map(lambda i: max(i, min_part), floor)
            self.floor = IntegerListsLex._list_function(floor, min_part)
            self.floor_type = "list"
            self.floor_limit = min_part
            self.floor_limit_start = len(floor)
        elif callable(floor):
            self._warning = True
            if min_part > 0:
                self.floor = lambda i: max(min_part, floor(i))
            else:
                self.floor = floor
            self.floor_type = "function"
            self.floor_limit = None
            self.floor_limit_start = Infinity
        else:
            raise TypeError("floor should be a list, tuple, or function")

        if ceiling is None:
            self.ceiling = ConstantFunction(max_part)
            self.ceiling_type = "constant"
            self.ceiling_limit = max_part
            self.ceiling_limit_start = 0
        elif isinstance(ceiling, (list, tuple)):
            if not all(i==Infinity or i in ZZ for i in ceiling):
                raise TypeError("the parts of ceiling={} should be non negative integers".format(ceiling))
            if not all(i >= 0 for i in ceiling):
                raise NotImplementedError("negative parts in floor={}".format(ceiling))
            if max_part < Infinity:
                ceiling = map(lambda i: min(i, max_part), ceiling)
            self.ceiling = IntegerListsLex._list_function(ceiling, max_part)
            self.ceiling_type = "list"
            self.ceiling_limit = max_part
            self.ceiling_limit_start = len(ceiling)
        elif callable(ceiling):
            self._warning = True
            if max_part < Infinity:
                self.ceiling = lambda i: min(max_part, ceiling(i))
            else:
                self.ceiling = ceiling
            self.ceiling = ceiling
            self.ceiling_type = "function"
            self.ceiling_limit = None
            self.ceiling_limit_start = Infinity
        else:
            raise ValueError("unable to parse value of parameter ceiling")

        if name is not None:
            self.rename(name)

        if self._warning and not self._waiver:
            warn("""
A function has been given as input of the floor=[...] or ceiling=[...]
arguments of IntegerListsLex. Please see the documentation for the caveats.
If you know what you are doing, you can set waiver=True to skip this warning.""")

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
        Checks whether the parameters give a proper lexicographic iterator.

        EXAMPLES::

            sage: IntegerListsLex(4).list()
            Traceback (most recent call last):
            ...
            ValueError: The specified parameters do not allow for a lexicographic iterator!

            sage: it = iter(IntegerListsLex(4, waiver=True))
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
        """
        if self._warning or self._waiver:
            return
        s = sum(self.floor(i) for i in range(self.floor_limit_start))
        if self.max_sum < Infinity and self.max_length == Infinity and self.floor_limit == 0:
            if self.min_slope<0 and self.max_slope>0 and s<self.min_sum:
                raise ValueError("The specified parameters do not allow for a lexicographic iterator!")
            if self.min_slope == 0 and s==0 and self.max_slope>0:
                if self.max_sum>0: # this is assuming that we remove trailing zeroes
                    raise ValueError("The specified parameters do not allow for a lexicographic iterator!")

    @staticmethod
    def _list_function(l, default):
        r"""
        Generate a function on the nonnegative integers from input.

        This method generates a function on the nonnegative integers
        whose values are taken from ``l`` when the input is a valid index
        in the list ``l``, and has a default value ``default`` otherwise.

        INPUT:

        - `l` -- a list to use as a source of values

        - `default` -- a default value to use for indices outside of the list

        OUTPUT:

        A function on the nonnegative integers.

        EXAMPLES::

            sage: C = IntegerListsLex(2, length=3)
            sage: f = C._list_function([1,2], Infinity)
            sage: f(1)
            2
            sage: f(3)
            +Infinity
        """
        return lambda i: l[i] if i < len(l) else default

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
        if self.min_sum == self.max_sum:
            return "Integer lists of sum {} satisfying certain constraints".format(self.min_sum)
        elif self.max_sum == Infinity:
            if self.min_sum == 0:
                return "Integer lists with arbitrary sum satisfying certain constraints"
            else:
                return "Integer lists of sum at least {} satisfying certain constraints".format(self.min_sum)
        else:
            return "Integer lists of sum between {} and {} satisfying certain constraints".format(self.min_sum,self.max_sum)

    def __contains__(self, comp):
        """
        Return ``True`` if ``comp`` meets the constraints imposed by the arguments.

        EXAMPLES::

            sage: C = IntegerListsLex(n=2, max_length=3, min_slope=0)
            sage: all([l in C for l in C])
            True
        """
        if len(comp) < self.min_length or len(comp) > self.max_length:
            return False
        n = sum(comp)
        if n < self.min_sum or n > self.max_sum:
            return False
        for i in range(len(comp)):
            if comp[i] < self.floor(i):
                return False
            if comp[i] > self.ceiling(i):
                return False
        for i in range(len(comp)-1):
            slope = comp[i+1] - comp[i]
            if slope < self.min_slope or slope > self.max_slope:
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
        return self._IntegerListsLexIter(self)

    class _IntegerListsLexIter:
        """
        Iterator class for IntegerListsLex
        """
        def __init__(self, parent):
            """
            TESTS::

                sage: C = IntegerListsLex(2, length=3)
                sage: I = IntegerListsLex._IntegerListsLexIter(C)
                sage: I.rho
                [(0, 2)]
                sage: I.mu
                [3]
                sage: I.j
                0
                sage: I.finished
                False
            """
            self.parent = parent

            parent._check_lexicographic_iterable()

            self.rho = []   # list of current search ranges
            self.mu = []    # list of integers
            self.j = -1     # index of last element of mu
            self.nu = 0     # sum of values in mu

            self.finished = False

            # initialize for beginning of iteration
            self._push_search()

        def _push_search(self):
            """
            Push search forward. Resetting attributes.

            EXAMPLES::

                sage: C = IntegerListsLex(2, length=3)
                sage: I = C.__iter__()
                sage: I.j
                0
                sage: I.rho
                [(0, 2)]
                sage: I.mu
                [3]
                sage: I.nu
                3
                sage: I._push_search()
                sage: I.j
                1
                sage: I.rho
                [(0, 2), (0, -1)]
                sage: I.mu
                [3, 0]
                sage: I.nu
                3
            """
            if self.j >= 0:
                prev = self.mu[self.j]
            else:
                prev = None
            self.j += 1
            interval = self._m_interval(self.j, self.parent.max_sum - self.nu, prev)
            val = interval[1] + 1 # iterator decrements before acting
            self.rho.append(interval)
            self.mu.append(val)
            self.nu += val

        def _pop_search(self):
            """
            Go back in search tree. Resetting attributes.

            EXAMPLES::

                sage: C = IntegerListsLex(2, length=3)
                sage: I = C.__iter__()
                sage: I.j
                0
                sage: I.rho
                [(0, 2)]
                sage: I.nu
                3
                sage: I.mu
                [3]
                sage: I._pop_search()
                sage: I.j
                -1
                sage: I.rho
                []
                sage: I.nu
                0
                sage: I.mu
                []
            """
            if self.j >= 0:
                self.j -= 1
                self.rho.pop()
                self.nu -= self.mu[-1]
                self.mu.pop()

        def next(self):
            """
            Return the next element in the iteration.

            EXAMPLES::

                sage: C = IntegerListsLex(2, length=3)
                sage: I = IntegerListsLex._IntegerListsLexIter(C)
                sage: I.next()
                [2, 0, 0]
                sage: I.next()
                [1, 1, 0]
            """
            if self.finished:
                raise StopIteration()

            rho = self.rho
            mu = self.mu
            p = self.parent
            min_sum = p.min_sum
            max_sum = p.max_sum
            min_length = p.min_length
            max_length = p.max_length

            while self.j >= 0: # j = -1 means that we have finished the bottom iteration

                # choose a new value of m to test

                mu[self.j] -= 1
                # m = mu[self.j]
                self.nu -= 1

                # check if the new value is valid, if not, pop to prefix, and now check if this is a solution

                if mu[self.j] < rho[self.j][0] or (self.j == max_length-1 and self.nu < min_sum):
                    self._pop_search()
                    if self._internal_list_valid():
                        return p._element_constructor(list(mu))
                    else:
                        continue

                # m is new and in range:
                # If we're at a leaf node, check if its a solution to return, pop the layer from the stack.
                # Otherwise, check if any solutions are possible with this value of m.
                #
                # Possible cases to detect leaf nodes:
                # 1. List is of maximal length
                # 2. List is of maximal sum, and of at least minimal length (allow padding zeros)

                if (self.nu == max_sum and self.j >= min_length - 1) or self.j == max_length - 1:
                    if self._internal_list_valid():
                        return p._element_constructor(list(mu))
                elif self._possible_m(mu[self.j], self.j,
                                      min_sum - (self.nu-mu[self.j]),
                                      max_sum - (self.nu-mu[self.j])):
                    self._push_search()

            self.finished = True
            raise StopIteration()

        def _internal_list_valid(self):
            """
            Return whether the current list in the iteration variable ``self.mu`` is a valid list.

            This method checks whether the sum of the parts in ``self.mu``
            is in the right range, whether its length is in the
            required range, and whether there are trailing zeroes.  It does not check all of the
            necessary conditions to verify that an arbitrary list satisfies the
            constraints from the corresponding ``IntegerListsLex`` object, and should
            not be used except internally in the iterator class.

            EXAMPLES::

                sage: C = IntegerListsLex(2, length=3)
                sage: I = IntegerListsLex._IntegerListsLexIter(C)
                sage: I.mu
                [3]
                sage: I._internal_list_valid()
                False
                sage: I.next()
                [2, 0, 0]
                sage: I.mu
                [2, 0, 0]
                sage: I._internal_list_valid()
                True
            """
            p = self.parent
            mu = self.mu
            nu = self.nu
            l = self.j + 1
            good_sum = (nu >= p.min_sum and nu <= p.max_sum)
            good_length = (l >= p.min_length and l <= p.max_length)
            no_trailing_zeros = (l <= max(p.min_length,0) or mu[-1] != 0)
            return good_sum and good_length and no_trailing_zeros

        def _upper_envelope(self, m, j):
            """
            Return function of the upper envelope starting with value ``m`` at position ``j``.

            INPUT:

            - ``m`` -- a nonnegative integer (starting value)

            - ``j`` -- a nonnegative integer (position)

            This method returns a function of ``i`` which returns the upper envelope if the starting value
            is ``m`` at position ``j``. The upper envelope is the minimum of the ceiling function and the
            value restriction given by the slope conditions.

            EXAMPLES::

                sage: C = IntegerListsLex(2, length=3)
                sage: I = IntegerListsLex._IntegerListsLexIter(C)
                sage: f = I._upper_envelope(1,1)
                sage: type(f)
                <type 'sage.misc.constant_function.ConstantFunction'>
                sage: f(1)
                inf
                sage: f(2)
                inf
                sage: C = IntegerListsLex(6, max_slope=1, max_part=3, max_length=6)
                sage: I = IntegerListsLex._IntegerListsLexIter(C)
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
            if self.parent.max_slope == Infinity:
                return self.parent.ceiling
            return lambda i: min(m + (i-j)*self.parent.max_slope, self.parent.ceiling(i) )

        def _lower_envelope(self, m, j):
            """
            Return function of the lower envelope starting with value ``m`` at position ``j``.

            INPUT:

            - ``m`` -- a nonnegative integer (starting value)

            - ``j`` -- a nonnegative integer (position)

            This method returns a function of ``i`` which returns the lower envelope if the starting value
            is ``m`` at position ``j``. The lower envelope is the maximum of the floor function and the
            value restriction given by the slope conditions.

            EXAMPLES::

                sage: C = IntegerListsLex(2, length=3)
                sage: I = IntegerListsLex._IntegerListsLexIter(C)
                sage: f = I._lower_envelope(1,1)
                sage: type(f)
                <type 'sage.misc.constant_function.ConstantFunction'>
                sage: f(1)
                0
                sage: f(2)
                0
                sage: C = IntegerListsLex(6, min_slope=-1, min_part=1)
                sage: I = IntegerListsLex._IntegerListsLexIter(C)
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
            if self.parent.min_slope == -Infinity:
                return self.parent.floor
            return lambda i: max( m + (i-j)*self.parent.min_slope, self.parent.floor(i) )

        def _m_interval(self, i, target_max, prev=None):
            r"""
            Return an interval for possible values of ``m`` (entry) at position ``i``.

            INPUT:

            - ``i`` -- a nonnegative integer (position)

            - ``target_max`` -- a nonnegative integer or +oo, the largest valid sum of a list tail.
              If +oo, the ``max_slope`` or ``ceiling`` restrictions must give a finite bound for the current part

            - ``prev`` -- a nonnegative integer, the last entry in the integer sequence prior to the desired tail,
              if the sequence is non-empty

            OUTPUT:

            A tuple of two integers which bounds an interval containing (in general this is proper containment)
            all possible suffix integers `m` which lead to a valid integer tail.

            EXAMPLES::

                sage: C = IntegerListsLex(2, length=3)
                sage: I = IntegerListsLex._IntegerListsLexIter(C)
                sage: I._m_interval(1,2)
                (0, 2)
            """
            p = self.parent

            lower_bounds = [0, p.floor(i)]
            upper_bounds = [target_max, p.ceiling(i)]
            if prev != None:
                lower_bounds.append(prev + p.min_slope)
                upper_bounds.append(prev + p.max_slope)
            lower_bound = max(lower_bounds)
            upper_bound = min(upper_bounds)

            ## check for infinite upper bound, in case target_max is infinite
            if upper_bound == Infinity:
                raise ValueError("infinite upper bound for values of m")

            return (lower_bound, upper_bound)

        def _possible_m(self, m, j, target_min, target_max):
            r"""
            INPUT:

            - ``m`` -- a nonnegative integer (value)

            - ``j`` -- a nonnegative integer (position)

            - ``target_min`` -- a nonnegative integer

            - ``target_max`` -- a nonnegative integer or +oo

            OUTPUT:

            Whether there exists a vector suffix `(v_j,...)` satisfying
            `v_j` equals ``m`` and the other conditions of ``self``, and
            with sum between ``target_min`` and ``target_max``.

            EXAMPLES::

                sage: C = IntegerListsLex(2, length=3)
                sage: I = IntegerListsLex._IntegerListsLexIter(C)
                sage: I._possible_m(1,1,2,2)
                True
                sage: I._possible_m(1,1,3,2)
                False
            """
            # Check code for various termination conditions.  Possible cases:
            # 0. interval [lower, upper] intersects interval [target_min, target_max] -- terminate True
            # 1. lower sum surpasses target_max -- terminate False
            # 2. iteration surpasses max_length -- terminate False
            # 3. upper envelope is smaller than lower envelope -- terminate False
            # 4. max_slope <= 0 -- terminate False after upper passes 0
            # 5. ceiling_limit == 0 -- terminate False after reaching larger limit point
            # 6. (uncomputable) ceiling function == 0 for all but finitely many input values -- terminate False after reaching (unknown) limit point -- currently hangs

            if target_min > target_max:
                return False

            p = self.parent

            loEnv = self._lower_envelope(m,j)
            upEnv = self._upper_envelope(m,j)
            lower = 0
            upper = 0

            # get to smallest valid number of parts
            for lo,up in [(loEnv(i) if i!=j else m, upEnv(i) if i!=j else m) for i in range(j, p.min_length-1)]:
                if lo > up:
                    return False
                lower += lo
                upper += up

            i = max(p.min_length-1,j)
            # check if any of the intervals intersect the target interval
            while i <= p.max_length - 1:
                lo = loEnv(i) if i!=j else m
                up = upEnv(i) if i!=j else m
                if lo > up:
                    break
                elif p.max_slope <= 0 and up <= 0 and self.j >= p.min_length:
                    break
                elif lower > target_max:
                    break
                elif p.ceiling_limit == 0 and i > p.ceiling_limit_start:
                    break
                else:
                    lower += lo
                    upper += up
                    i += 1
                if lower <= target_max and upper >= target_min and lower <= upper:
                    return True
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
