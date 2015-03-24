r"""
Tools for generating lists of integers in lexicographic order
"""
#*****************************************************************************
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

from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.misc.constant_function import ConstantFunction
from sage.categories.enumerated_sets import EnumeratedSets
from sage.structure.list_clone import ClonableArray
from sage.structure.parent import Parent
from sage.structure.list_clone import ClonableArray
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from warnings import warn

infinity = float('+inf')

class IntegerListsLex(Parent):
    r"""
    A combinatorial class `C` for integer lists satisfying certain
    sum, length, upper/lower bound and regularity constraints. The
    purpose of this tool is mostly to provide a Constant Amortized
    Time iterator through these lists, in lexicographic order.

    INPUT:

    - ``min_n`` -- a nonnegative integer specifying the minimum number to which
      the elements in the list sum; defaults to ``0``
    - ``max_n`` -- a nonnegative integer or `\infty` specifying the maximum number to which
      the elements in the list sum; defaults to ``0``
    - ``n`` -- a nonnegative integer or list of nonnegative integers; overrides min_n and max_n if specified
    - ``min_length`` -- a nonnegative integer specifying the minimal length of the vectors;
      defaults to ``0``
    - ``max_length`` -- a nonnegative integer or `\infty` specifying the maximal length of the vectors;
      defaults to `\infty`
    - ``length`` -- an integer; overrides min_length and max_length if specified
    - ``floor`` -- an integer, a list of integers, or a function `f`;
      defaults to the constant zero function
    - ``ceiling`` -- an integer, a list of integers, or a function `f` (or list);
      defaults to the constant `\infty` function
    - ``min_part`` -- a nonnegative integer (default: None, in which case internally it defaults to 0)
    - ``max_part`` -- a nonnegative integer (default: None, in which case internally it defaults to `\infty`)
    - ``min_slope`` -- an integer or `-\infty`; defaults to `-\infty`
    - ``max_slope`` -- an integer or `+\infty`; defaults to `+\infty`
    - ``category`` -- a category (default: FiniteEnumeratedSets)
    - ``waiver`` -- boolean (default: False)

    An *integer list* is a list `l` of nonnegative integers, its
    *parts*. The *length* of `l` is the number of its parts;
    the *sum* of `l` is the sum of its parts.

    .. NOTE::

       Two valid integer lists are considered equivalent if they only
       differ by trailing zeroes. In this case, only the list with the
       least number of trailing zeroes will be produced.

    The constraints on the lists are as follows:

    - Sum: `min_n \le sum(l) \le max_n` (with ``n = min_n = max_n`` if ``n`` is specified)

    - Length: ``min_length <= len(l) <= max_length`` (with ``length = min_length = max_length`` if
      ``length`` is specified)

    - Lower and upper bounds: ``max(floor(i), min_part) <= l[i] <= min(ceiling(i), max_part)``, for
      ``i`` from 0 to ``len(l)``. If ``floor`` (resp. ``ceiling``) is an
      integer `k`, then ``floor`` (resp. ``ceiling``) is considered to be
      the constant function `k`.
      If ``floor`` is a list of integers `v`, then ``floor(i) = v(i)``
      for ``i`` from ``0`` to ``len(v)-1``. For ``i`` from ``len(v)`` to ``max_length-1``,
      the ``floor`` function is considered to be zero.
      Similarly, if ``ceiling`` is a list of integers `v`, then ``ceiling(i) = v(i)``
      for ``i`` from ``0`` to ``len(v)-1``. For ``i`` from ``len(v)`` to ``max_length-1``,
      the ``ceiling`` function is considered to be `\infty`.

    - Regularity condition: ``minSlope <= l[i+1]-l[i] <= maxSlope``,
      for ``i`` from 0 to ``len(l)-1``

    This is a generic low level tool. The interface has been designed
    with efficiency in mind. It is subject to incompatible changes in
    the future. More user friendly interfaces are provided by high
    level tools like :class:`Partitions` or :class:`Compositions`.

    EXAMPLES:

    We create the combinatorial class of lists of length 3 and sum 2::

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

        sage: list(IntegerListsLex(5, length = 3, floor = [1,2,0], ceiling = [3,2,3]))
        [[3, 2, 0], [2, 2, 1], [1, 2, 2]]

    Using the slope condition, one can generate integer partitions
    (but see :mod:`sage.combinat.partition.Partitions`)::

        sage: list(IntegerListsLex(4, max_slope=0))
        [[4], [3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1]]

    The following is the list of all partitions of `7` with parts at least `2`::

        sage: list(IntegerListsLex(7, max_slope = 0, floor = 2))
        [[7], [5, 2], [4, 3], [3, 2, 2]]

    Next we list all partitions of `5` and length at most 3
    which are bounded below by [2,1,1]::

        sage: list(IntegerListsLex(5, max_slope = 0, max_length = 3, floor = [2,1,1]))
        [[5], [4, 1], [3, 2], [3, 1, 1], [2, 2, 1]]

    Note that ``[5]`` is considered valid, because the lower bound
    constraint only apply to existing positions in the list. To
    obtain instead the partitions containing ``[2,1,1]``, one need to
    use ``min_length``::

        sage: list(IntegerListsLex(5, max_slope = 0, min_length = 3, max_length = 3, floor = [2,1,1]))
        [[3, 1, 1], [2, 2, 1]]

    This is the list of all partitions of `5` which are contained in
    ``[3,2,2]``::

        sage: list(IntegerListsLex(5, max_slope = 0, max_length = 3, ceiling = [3,2,2]))
        [[3, 2], [3, 1, 1], [2, 2, 1]]

    This is the list of all compositions of `4` (but see Compositions)::

        sage: list(IntegerListsLex(4, floor = 1))
        [[4], [3, 1], [2, 2], [2, 1, 1], [1, 3], [1, 2, 1], [1, 1, 2], [1, 1, 1, 1]]

    This is the list of all integer vectors of sum `4` and length `3`::

        sage: list(IntegerListsLex(4, length = 3))
        [[4, 0, 0], [3, 1, 0], [3, 0, 1], [2, 2, 0], [2, 1, 1],
         [2, 0, 2], [1, 3, 0], [1, 2, 1], [1, 1, 2], [1, 0, 3],
         [0, 4, 0], [0, 3, 1], [0, 2, 2], [0, 1, 3], [0, 0, 4]]

    Next we obtain all lists of sum 4 and length 4 such that l[i] <= i::

        sage: list(IntegerListsLex(4, length = 4, ceiling = lambda i: i, waiver=True))
        [[0, 1, 2, 1], [0, 1, 1, 2], [0, 1, 0, 3], [0, 0, 2, 2], [0, 0, 1, 3]]

    Note that when passing a function as the ceiling, the existence of a vector
    with the specified conditions might be undecidable. For example, when ``ceiling`` is
    a function and is zero for a long time, there might not be a solution up to an arbitrarily
    high position. But the user specified function could in principle then increase and allow
    for a solution::

        sage: it = IntegerListsLex(2,ceiling=lambda i:0 if i<20 else 1).__iter__()
        doctest:...
        UserWarning: When the user specifies a method, then (s)he is responsible that the algorithm
        will not hang. Also note that the specified function should start at 0 rather than 1.
        Before trac#17979 the indexing was ambiguous and sometimes started at 1.
        sage: it.next()
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1]
        sage: it.next()
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1]
        sage: it.next()
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1]

    In the above example, there are infinitely many solutions, but they are returned.
    In the next example, there are no solutions, but since the ceiling is given as a method,
    the computer cannot decide whether there will be a solution or not. The code hangs::

        sage: IntegerListsLex(2,ceiling=lambda i:0).list() # not tested

    For this reason, the user needs to sign a waiver when providing a method rather than a
    list. The same example does work when the computer is told that the function is a constant
    zero:::

        sage: IntegerListsLex(2,ceiling=0).list()
        []

    This is the list of all monomials of degree `4` which divide the
    monomial `x^3y^1z^2` (a monomial being identified with its
    exponent vector)::

        sage: R.<x,y,z> = QQ[]
        sage: m = [3,1,2]
        sage: def term(exponents):
        ...       return x^exponents[0] * y^exponents[1] * z^exponents[2]
        ...
        sage: list( IntegerListsLex(4, length = len(m), ceiling = m, element_constructor = term) )
        [x^3*y, x^3*z, x^2*y*z, x^2*z^2, x*y*z^2]

    Note the use of the element_constructor feature.

    In general, the complexity of the iteration algorithm is constant
    time amortized for each integer list produced.  There is one
    degenerate case though where the algorithm may run forever without
    producing anything. If max_length is `+\infty` and max_slope `>0`,
    testing whether there exists a valid integer list of sum `n` may
    be only semi-decidable. In the following example, the algorithm
    will enter an infinite loop, because it needs to decide whether
    `ceiling(i)` is nonzero for some `i`::

        sage: list( IntegerListsLex(1, ceiling = lambda i: 0) ) # todo: not implemented

    .. NOTE::

       Caveat: counting is done by brute force generation. In some
       special cases, it would be possible to do better by counting
       techniques for integral point in polytopes.

    In the following example, the floor conditions do not satisfy the
    slope conditions since the floor for the third part is also 3. The algorithm
    will nonetheless give the correct result::

        sage: I = IntegerListsLex(16, min_length=2, max_slope=-1, floor=[5,3,3])
        sage: I.list()
        [[13, 3], [12, 4], [11, 5], [10, 6], [9, 7], [9, 4, 3], [8, 5, 3], [8, 4, 3, 1],
         [7, 6, 3], [7, 5, 4], [7, 5, 3, 1], [7, 4, 3, 2], [6, 5, 4, 1], [6, 5, 3, 2],
         [6, 4, 3, 2, 1]]

    .. NOTE::

        The generation algorithm could in principle be extended to deal with non-constant
        slope constraints and with negative parts.

    .. TODO:

        Integrate all remaining tests from
        http://mupad-combinat.svn.sourceforge.net/viewvc/mupad-combinat/trunk/MuPAD-Combinat/lib/COMBINAT/TEST/MachineIntegerListsLex.tst

    TESTS::

        sage: g = lambda x: lambda i: x
        sage: list(IntegerListsLex(0, floor = g(1), min_slope = 0))
        [[]]
        sage: list(IntegerListsLex(0, floor = g(1), min_slope = 0, max_slope = 0))
        [[]]
        sage: list(IntegerListsLex(0, max_length=0, floor = g(1), min_slope = 0, max_slope = 0))
        [[]]
        sage: list(IntegerListsLex(0, max_length=0, floor = g(0), min_slope = 0, max_slope = 0))
        [[]]
        sage: list(IntegerListsLex(0, floor = 1, min_slope = 0))
        [[]]
        sage: list(IntegerListsLex(1, floor = 1, min_slope = 0))
        [[1]]
        sage: list(IntegerListsLex(0, min_length = 1, floor = 1, min_slope = 0))
        []
        sage: list(IntegerListsLex(0, min_length = 1, min_slope = 0))
        [[0]]
        sage: list(IntegerListsLex(3, max_length=2, ))
        [[3], [2, 1], [1, 2], [0, 3]]
        sage: partitions = {"floor": 1, "max_slope": 0}
        sage: partitions_min_2 = {"floor": g(2), "max_slope": 0}
        sage: compositions = {"floor": 1}
        sage: integer_vectors = lambda l: {"length": l}
        sage: lower_monomials = lambda c: {"length": c, "floor": lambda i: c[i]}
        sage: upper_monomials = lambda c: {"length": c, "ceiling": lambda i: c[i]}
        sage: constraints = { "floor":1, "min_slope": -1, "max_slope": 0}
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

    Noted on :trac:`17898`::

        sage: list(IntegerListsLex(4, floor=1, length=3, min_slope=1))
        []
        sage: IntegerListsLex(6, ceiling=[4,2], floor=[3,3]).list()
        []
        sage: IntegerListsLex(6, floor=1, ceiling=3, max_slope=-4).list()
        []

    Noted in :trac:`17548`, which are now fixed::

        sage: IntegerListsLex(10, floor=2, max_slope=-1).list()
        [[10], [8, 2], [7, 3], [6, 4], [5, 3, 2]]
        sage: IntegerListsLex(5, ceiling=2, min_slope=1, floor=[2,1,1]).list()
        []
        sage: IntegerListsLex(4, min_slope=0, max_slope=0).list()
        [[4], [2, 2], [1, 1, 1, 1]]
        sage: IntegerListsLex(6, min_slope=-1, max_slope=-1).list()
        [[6], [3, 2, 1]]
        sage: IntegerListsLex(6, min_length=3, max_length=2, floor=1).list()
        []
        sage: I = IntegerListsLex(3, max_length=2, floor=1)
        sage: I.list()
        [[3], [2, 1], [1, 2]]
        sage: [1,1,1] in I
        False
        sage: I = IntegerListsLex(10, ceiling=[4], max_length=1, floor = 1)
        sage: I.list()
        []
        sage: [4,6] in I
        False
        sage: I = IntegerListsLex(4, ceiling=2, min_slope=1, floor=1)
        sage: I.list()
        []
        sage: I = IntegerListsLex(7, ceiling=4, min_slope=1, floor=1)
        sage: I.list()
        [[3, 4], [1, 2, 4]]
        sage: I = IntegerListsLex(4, floor=[2,1], ceiling=[2,2], max_length=2, min_slope=0)
        sage: I.list()
        [[2, 2]]
        sage: I = IntegerListsLex(10, floor=1, max_slope=-1)
        sage: I.list()
        [[10], [9, 1], [8, 2], [7, 3], [7, 2, 1], [6, 4], [6, 3, 1], [5, 4, 1],
         [5, 3, 2], [4, 3, 2, 1]]

    .. RUBRIC:: Some comparative timings with the former implementation

        ::

            sage: from sage.combinat.integer_list_old import IntegerListsLex as IntegerListsLexOld

            sage: P = IntegerListsLex(n=20, max_slope=0, floor=1)
            sage: %time x = list(P)  # not tested
            CPU times: user 264 ms, sys: 0 ns, total: 264 ms
            Wall time: 276 ms
            sage: P = IntegerListsLexOld(n=20, max_slope=0, min_part=1)
            sage: %time x = list(P)  # not tested
            CPU times: user 249 ms, sys: 15.8 ms, total: 265 ms
            Wall time: 401 ms
            sage: P.cardinality()
            627

            sage: P = IntegerListsLex(n=30, max_slope=0, floor=1)
            sage: %time x = list(P)  # not tested
            CPU times: user 2.98 s, sys: 6.52 ms, total: 2.99 s
            Wall time: 3.81 s
            sage: P = IntegerListsLexOld(n=30, max_slope=0, min_part=1)
            sage: %time x = list(P)  # not tested
            CPU times: user 2.55 s, sys: 7.05 ms, total: 2.55 s
            Wall time: 3.49 s
            sage: P.cardinality()
            5604

    TESTS:

    Internally, the iterator works on a single list that is mutated
    along the way. This tests makes sure that we do make a copy of
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

    def __init__(self,
                 n=None, min_n=0, max_n=0,
                 length=None, min_length=0, max_length=infinity,
                 floor=None, ceiling=None,
                 min_part=None, max_part=None,
                 min_slope=-infinity, max_slope=infinity,
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

        self.waiver = waiver
        self.warning = False # warning for dangerous (but possibly valid) usage

        if n is not None:
            if n in ZZ:
                #if n < 0:
                #    print n
                #    raise ValueError("value of n can't be less than 0")
                self.n_list = [n]
                self.min_n = n
                self.max_n = n
            elif isinstance(n, list):
                self.n_list = n
                self.min_n = min(self.n_list)
                #if self.min_n < 0:
                #    raise ValueError("can't have negative values of n")
                self.max_n = max(self.n_list)
            else:
                raise ValueError("invalid value for parameter n")
        else:
            self.n_list = None
            self.min_n = min_n
            self.max_n = max_n

        if length is not None:
            self.min_length = length
            self.max_length = length
        else:
            self.min_length = min_length
            self.max_length = max_length

        self.min_slope = min_slope
        self.max_slope = max_slope

        if min_part is None:
            self.min_part = 0
        elif min_part in ZZ:
            if min_part < 0:
                raise(ValueError("min_part can't be negative"))
            else:
                ## TODO warn user of deprecation
                self.min_part = min_part
        else:
            if min_part == infinity:
                raise(ValueError("min_part can't be infinite"))
            elif min_part == -infinity:
                raise(ValueError("min_part can't be negative or infinite"))
            else:
                raise(ValueError("unable to parse value of min_part"))

        if max_part is None:
            self.max_part = infinity
        elif max_part in ZZ:
            if max_part < 0:
                raise(ValueError("max_part can't be negative"))
            else:
                ## TODO warn user of deprecation
                self.max_part = max_part
        else:
            if max_part == -infinity:
                raise(ValueError("max_part can't be negative"))
            else:
                raise(ValueError("unable to parse value of max_part"))

        ## TODO: Warn user about using arbitrary functions?
        # indexing for floor and ceiling functions is base 0 internally
        if floor is None:
            self.floor_type = "none"
            self.floor = ConstantFunction(self.min_part)
            self.floor_limit = 0
            self.floor_limit_start = 0
        elif floor in ZZ:
            if floor < 0:
                raise(ValueError("floor value can't be negative"))
            self.floor_type = "constant"
            self.floor = ConstantFunction(max(floor,self.min_part))
            self.floor_limit = floor
            self.floor_limit_start = 0
        elif isinstance(floor, list) or isinstance(floor, tuple):
            self.floor_type = "list"
            if isinstance(floor, tuple):
                (floor,default) = floor
            else:
                (floor,default) = (floor, 0)
            floor_fn = IntegerListsLex._list_function(floor, default)
            self.floor_limit = 0
            self.floor_limit_start = len(floor)
            if (self.min_part > 0):
                self.floor = lambda i: max(self.min_part, floor_fn(i))
            else:
                self.floor = floor_fn
        elif callable(floor):
            # indexing is base 0
            self.warning = True
            self.floor_type = "function"
            self.floor_limit = None
            self.floor_limit_start = float('+inf')
            if (self.min_part > 0):
                self.floor = lambda i: max(self.min_part, floor(i))
            else:
                self.floor = floor
        else:
            raise(ValueError("unable to parse value of parameter floor"))

        if ceiling is None:
            self.ceiling_type = "none"
            self.ceiling = ConstantFunction(min(infinity,self.max_part))
            self.ceiling_limit = infinity
            self.ceiling_limit_start = 0
        elif ceiling in ZZ or ceiling == infinity:
            if ceiling < 0:
                raise(ValueError("ceiling value can't be negative"))
            self.ceiling_type = "constant"
            self.ceiling = ConstantFunction(min(ceiling,self.max_part))
            self.ceiling_limit = ceiling
            self.ceiling_limit_start = 0
        elif isinstance(ceiling, list) or isinstance(ceiling, tuple):
            self.ceiling_type = "list"
            if isinstance(ceiling, tuple):
                (ceiling, default) = ceiling
            else:
                (ceiling, default) = (ceiling, infinity)
            ceiling_fn = IntegerListsLex._list_function(ceiling, default)
            self.ceiling_limit = infinity
            self.ceiling_limit_start = len(ceiling)
            if (self.max_part < infinity):
                self.ceiling = lambda i: min(self.max_part, ceiling_fn(i))
            else:
                self.ceiling = ceiling_fn
        elif callable(ceiling):
            self.warning = True
            self.ceiling_type = "function"
            self.ceiling_limit = None
            self.ceiling_limit_start = infinity
            if (self.max_part < infinity):
                self.ceiling = lambda i: min(self.max_part, ceiling(i))
            else:
                self.ceiling = ceiling
        else:
            raise(ValueError("unable to parse value of parameter ceiling"))

        if name is not None:
            self.rename(name)

        if self.warning and not self.waiver:
            warn("""When the user specifies a method, then (s)he is responsible that the algorithm
            will not hang. Also note that the specified function should start at 0 rather than 1.
            Before trac#17979 the indexing was ambiguous and sometimes started at 1.""")

        # In case we want output to be of a different type,
        if element_constructor is not None:
            self._element_constructor_ = element_constructor
        if element_class is not None:
            self.Element = element_class
        if global_options is not None:
            self.global_options = global_options

        Parent.__init__(self, category=category)

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
            sage: f = C._list_function([1,2],Infinity)
            sage: f(1)
            2
            sage: f(3)
            +Infinity
        """
        return lambda i: l[i] if (i >= 0 and i < len(l) and i in ZZ) else default

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

            sage: C = IntegerListsLex([1,2,4], length=3)
            sage: C # indirect doctest
            Integer lists of sum in [1, 2, 4] satisfying certain constraints

            sage: C = IntegerListsLex([1,2,4], length=3, name="A given name")
            sage: C
            A given name
        """
        if self.min_n == self.max_n:
            return "Integer lists of sum %s satisfying certain constraints" % self.min_n
        elif self.max_n == infinity:
            if self.min_n == 0:
                return "Integer lists with arbitrary sum satisfying certain constraints"
            else:
                return "Integer lists of sum at least %s satisfying certain constraints" % self.min_n
        else:
            if self.n_list is not None:
                return "Integer lists of sum in %s satisfying certain constraints" % self.n_list
            else:
                return "Integer lists of sum between %s and %s satisfying certain constraints" % (self.min_n,self.max_n)

    def __contains__(self, comp):
        """
        Return ``True`` if ``comp`` meets the constraints imposed by the arguments.

        EXAMPLES::

            sage: C = IntegerListsLex(n=2,max_length=3,min_slope=0)
            sage: all([l in C for l in C])
            True
        """
        if len(comp) < self.min_length or len(comp) > self.max_length:
            return False
        n = sum(comp)
        if self.n_list is not None:
            if not n in self.n_list:
                return False
        else:
            if n < self.min_n or n > self.max_n:
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
            sage: list(C) #indirect doctest
            [[2, 0, 0], [1, 1, 0], [1, 0, 1], [0, 2, 0], [0, 1, 1], [0, 0, 2]]
        """
        return IntegerListsLex._IntegerListsLexIter(self)

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
                []
                sage: I.mu
                []
                sage: I.j
                0
                sage: I.finished
                False
            """
            self.parent = parent

            ## TODO either use SearchForest or use ranges / indices rather than endpoints / values (more modular)
            self.rho = []   # vector of current search ranges
            self.mu = []    # vector of elements
            self.j = 0      # length of current list

            self.finished = False

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
            min_n = p.min_n
            max_n = p.max_n
            min_length = p.min_length
            max_length = p.max_length

            while self.j >= 0: # j = -1 means that we have finished the bottom iteration
                # choose a new value of m to test

                if len(rho) > self.j:   # decreasing the value of a given m
                    mu[self.j] -= 1
                else:                   # add new range of m-values for index j
                    new_interval = None
                    if self.j > 0:
                        new_interval = self._m_interval(self.j, max_n - sum(mu), mu[self.j-1])
                    else:
                        new_interval = self._m_interval(self.j, max_n - sum(mu))
                    rho.append( new_interval )
                    mu.append( rho[self.j][1] ) # start at largest value

                m = mu[self.j]
                nu = sum(mu)

                # check if the new value is valid, if not, pop to prefix, and now check if this is a solution

                if m < rho[self.j][0]:
                    mu.pop()
                    rho.pop()
                    self.j -= 1
                    if self._internal_list_valid(mu):
                        return p._element_constructor(list(mu))
                    else:
                        continue

                # m is new and in range
                # if we're at a solution that's obviously next in order, return it
                # otherwise, check if any solutions are possible with this value of m
                if (nu == max_n and self.j >= min_length - 1) or self.j == max_length - 1:
                    if self._internal_list_valid(mu):
                        return p._element_constructor(list(mu))
                elif self._possible_m(m, self.j, min_n - (nu-m), max_n - (nu-m)):
                    self.j += 1

            self.finished = True
            raise StopIteration()

        def _internal_list_valid(self, mu):
            """
            Return whether ``mu`` is a valid list.

            INPUT:

            - ``mu`` -- a list of integers

            This method checks whether the sum of the entries in ``mu``
            is in the right range, whether its length is in the required
            range, and whether there are trailing zeroes.

            EXAMPLES::

                sage: C = IntegerListsLex(2, length=3)
                sage: I = IntegerListsLex._IntegerListsLexIter(C)
                sage: I._internal_list_valid([2,0,0])
                True
                sage: I._internal_list_valid([1,0,0])
                False
            """
            p = self.parent
            nu = sum(mu)
            j = len(mu)
            if p.n_list is not None:
                good_sum = (nu in p.n_list)
            else:
                good_sum = (nu >= p.min_n and nu <= p.max_n)
            good_length = (j >= p.min_length and j <= p.max_length)
            no_trailing_zeros = (j <= max(p.min_length,0) or mu[-1] != 0)
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
                sage: C = IntegerListsLex(6, max_slope=1, ceiling=3)
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
            if self.parent.max_slope == infinity:
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
                sage: C = IntegerListsLex(6, min_slope=-1, floor=1)
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
            if self.parent.min_slope == -infinity:
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
            if upper_bound == infinity:
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
