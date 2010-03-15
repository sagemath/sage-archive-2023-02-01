"""
Integer Range

AUTHORS:

 - Nicolas Borie  (2010-03): First release.
 - Florent Hivert (2010-03): Added a class factory + cardinality method.
"""
#*****************************************************************************
#  Copyright (C) 2010 Nicolas Borie <nicolas.borie@math.u-psud.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.parent import Parent
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.sets.finite_enumerated_set import FiniteEnumeratedSet
from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.integer import Integer
from sage.rings.infinity import Infinity, MinusInfinity, PlusInfinity

class IntegerRange(UniqueRepresentation, Parent):
    r"""
    The class of :class:`Integer <sage.rings.integer.Integer>` ranges

    Returns an enumerated set containing an arithmetic progression of integers.

    INPUT:

    - ``begin``        -- an integer, Infinity or -Infinity
    - ``end``          -- an integer, Infinity or -Infinity
    - ``step``         -- a non zero integer (default to 1)
    - ``middle_point`` -- an integer inside the set (default to ``None``)

    OUTPUT:

    A parent in the category :class:`FiniteEnumeratedSets()
    <sage.categories.finite_enumerated_sets.FiniteEnumeratedSets>` or
    :class:`InfiniteEnumeratedSets()
    <sage.categories.infinite_enumerated_sets.InfiniteEnumeratedSets>`
    depending on the arguments defining ``self``.

    ``IntegerRange(i, j)`` returns the set of `\{i, i+1, i+2, \dots , j-1\}`.
    ``start`` (!) defaults to 0. When ``step`` is given, it specifies the
    increment. The default increment is `1`. IntegerRange allows ``begin`` and
    ``end`` to be infinite.

    ``IntegerRange`` is designed to have similar interface Python
    range. However, whereas ``range`` accept and returns Python ``int``,
    ``IntegerRange`` deals with :class:`Integer <sage.rings.integer.Integer>`.

    If ``middle_point`` is given, then the elements are generated starting
    from it, in a alternating way: `\{m, m+1, m-2, m+2, m-2 \dots \}`.

    EXAMPLES::

        sage: list(IntegerRange(5))
        [0, 1, 2, 3, 4]
        sage: list(IntegerRange(2,5))
        [2, 3, 4]
        sage: I = IntegerRange(2,100,5); I
        Integer progression from 2 to 100 with increment 5
        sage: list(I)
        [2, 7, 12, 17, 22, 27, 32, 37, 42, 47, 52, 57, 62, 67, 72, 77, 82, 87, 92, 97]
        sage: I.category()
        Category of finite enumerated sets
        sage: I[1].parent()
        Integer Ring

    When ``begin`` and ``end`` are both finite, ``IntegerRange(begin, end,
    step)`` is the set whose list of elements is equivalent to the python
    construction ``range(begin, end, step)``::

        sage: list(IntegerRange(4,105,3)) == range(4,105,3)
        True
        sage: list(IntegerRange(-54,13,12)) == range(-54,13,12)
        True

    Except for the type of the numbers::

        sage: type(IntegerRange(-54,13,12)[0]), type(range(-54,13,12)[0])
        (<type 'sage.rings.integer.Integer'>, <type 'int'>)

    When ``begin`` is finite and ``end`` is +Infinity, ``self`` is the infinite
    arithmetic progression starting from the ``begin`` by step ``step``::

        sage: I = IntegerRange(54,Infinity,3); I
        Integer progression from 54 to +Infinity with increment 3
        sage: I.category()
        Category of infinite enumerated sets
        sage: p = iter(I)
        sage: (p.next(), p.next(), p.next(), p.next(), p.next(), p.next())
        (54, 57, 60, 63, 66, 69)

        sage: I = IntegerRange(54,-Infinity,-3); I
        Integer progression from 54 to -Infinity with increment -3
        sage: I.category()
        Category of infinite enumerated sets
        sage: p = iter(I)
        sage: (p.next(), p.next(), p.next(), p.next(), p.next(), p.next())
        (54, 51, 48, 45, 42, 39)

    When ``begin`` and ``end`` are both infinite, you will have to specify the
    extra argument ``middle_point``. ``self`` is then defined by a point
    and a progression/regression setting by ``step``. The enumeration
    is done this way: (let us call `m` the ``middle_point``)
    `\{m, m+step, m-step, m+2step, m-2step, m+3step, \dots \}`::

        sage: I = IntegerRange(-Infinity,Infinity,37,-12); I
        Integer progression containing -12 with increment 37 and bounded with -Infinity and +Infinity
        sage: I.category()
        Category of infinite enumerated sets
        sage: -12 in I
        True
        sage: -15 in I
        False
        sage: p = iter(I)
        sage: (p.next(), p.next(), p.next(), p.next(), p.next(), p.next(), p.next(), p.next())
        (-12, 25, -49, 62, -86, 99, -123, 136)

    It is also possible to use the argument ``middle_point`` for other cases, finite
    or infinite. The set will be the same as if you didn't give this extra argument
    but the enumeration will begin with this ``middle_point``::

        sage: I = IntegerRange(123,-12,-14); I
        Integer progression from 123 to -12 with increment -14
        sage: list(I)
        [123, 109, 95, 81, 67, 53, 39, 25, 11, -3]
        sage: J = IntegerRange(123,-12,-14,25); J
        Integer progression containing 25 with increment -14 and bounded with 123 and -12
        sage: list(J)
        [25, 11, 39, -3, 53, 67, 81, 95, 109, 123]

    Remember that, like for range, if you define a non empty set, ``begin`` is
    supposed to be included and ``end`` is supposed to be excluded. In the same
    way, when you define a set with a ``middle_point``, the ``begin`` bound will
    be supposed to be included and the ``end`` bound supposed to be excluded::

        sage: I = IntegerRange(-100,100,10,0)
        sage: J = range(-100,100,10)
        sage: 100 in I
        False
        sage: 100 in J
        False
        sage: -100 in I
        True
        sage: -100 in J
        True
        sage: list(I)
        [0, 10, -10, 20, -20, 30, -30, 40, -40, 50, -50, 60, -60, 70, -70, 80, -80, 90, -90, -100]

    TESTS::

        sage: # Some category automatic tests
        sage: TestSuite(IntegerRange(2,100,3)).run()
        sage: TestSuite(IntegerRange(564,-12,-46)).run()
        sage: TestSuite(IntegerRange(2,Infinity,3)).run()
        sage: TestSuite(IntegerRange(732,-Infinity,-13)).run()
        sage: TestSuite(IntegerRange(-Infinity,Infinity,3,2)).run()
        sage: TestSuite(IntegerRange(56,Infinity,12,80)).run()
        sage: TestSuite(IntegerRange(732,-12,-2743,732)).run()
        sage: # 20 random tests: range and IntegerRange give the same set for finite cases
        sage: for i in range(20):
        ...       begin = Integer(randint(-300,300))
        ...       end = Integer(randint(-300,300))
        ...       step = Integer(randint(-20,20))
        ...       if step == 0:
        ...           step = Integer(1)
        ...       assert list(IntegerRange(begin, end, step)) == range(begin, end, step)
        sage: # 20 random tests: range and IntegerRange with middle point for finite cases
        sage: for i in range(20):
        ...       begin = Integer(randint(-300,300))
        ...       end = Integer(randint(-300,300))
        ...       step = Integer(randint(-15,15))
        ...       if step == 0:
        ...           step = Integer(-3)
        ...       I = IntegerRange(begin, end, step)
        ...       if I.cardinality() == 0:
        ...           assert len(range(begin, end, step)) == 0
        ...       else:
        ...           TestSuite(I).run()
        ...           L1 = list(IntegerRange(begin, end, step, I.an_element()))
        ...           L2 = range(begin, end, step)
        ...           L1.sort()
        ...           L2.sort()
        ...           assert L1 == L2

    TODO: Fix the following test as soon as categories test allows for it
    see #8543::

        sage: TestSuite(IntegerRange(0, 5, -1)).run()
        Failure in _test_an_element:
        Traceback (most recent call last):
        ...
        AssertionError: FIXME (#8543): categories doesn't allow empty sets
        ------------------------------------------------------------
        The following tests failed: _test_an_element, _test_elements
    """

    @staticmethod
    def __classcall_private__(cls, begin, end=None, step=Integer(1), middle_point=None):
        """
        TESTS::

            sage: IntegerRange(2,5,0)
            Traceback (most recent call last):
            ...
            ValueError: IntegerRange() step argument must not be zero
            sage: IntegerRange(2) is IntegerRange(0, 2)
            True
        """
        if end is None:
            end = begin
            begin = Integer(0)
        # check of the arguments
        assert isinstance(begin, (Integer, MinusInfinity, PlusInfinity))
        assert isinstance(end, (Integer, MinusInfinity, PlusInfinity))
        assert isinstance(step, Integer)
        if step.is_zero():
            raise ValueError("IntegerRange() step argument must not be zero")

        # If begin and end are infinite, middle_point and step will defined the set.
        if begin == -Infinity and end == Infinity:
            if middle_point == None:
                raise ValueError("Can't iterate over this set, please provide middle_point")

        # If we have a middle point, we go on the special enumeration way...
        if middle_point != None:
            return IntegerRangeFromMiddle(begin, end, step, middle_point)

        if (begin == -Infinity) or (begin == Infinity):
            raise ValueError("Can't iterate over this set: It is impossible to begin an enumeration with plus/minus Infinity")

        if end != Infinity and end != -Infinity:
            return IntegerRangeFinite(begin, end, step)
        else:
            return IntegerRangeInfinite(begin, end, step)

    def __repr__(self):
        r"""
        TESTS::

            sage: from sage.sets.integer_range import IntegerRangeFinite
            sage: IntegerRangeFinite(123,12,-4)
            Integer progression from 123 to 12 with increment -4
            sage: IntegerRangeFinite(-57,1,3)
            Integer progression from -57 to 1 with increment 3

            sage: from sage.sets.integer_range import IntegerRangeInfinite
            sage: IntegerRangeInfinite(-57,Infinity,8)
            Integer progression from -57 to +Infinity with increment 8
            sage: IntegerRangeInfinite(-112,-Infinity,-13)
            Integer progression from -112 to -Infinity with increment -13
        """
        return "Integer progression from %s to %s with increment %s"%(self._begin,self._end,self._step)

    element_class = Integer


class IntegerRangeFinite(IntegerRange):
    r"""
    The class of finite enumerated sets of integers defined by finite
    arithmetic progressions

    See :class:`IntegerRange` for more details.
    """
    def __init__(self, begin, end, step=Integer(1)):
        r"""
        TESTS::

            sage: from sage.sets.integer_range import IntegerRangeFinite
            sage: I = IntegerRangeFinite(123,12,-4)
            sage: I.category()
            Category of finite enumerated sets
            sage: TestSuite(I).run()
        """
        self._begin = begin
        self._end = end
        self._step = step
        Parent.__init__(self, category = FiniteEnumeratedSets())

    def __contains__(self, elt):
        r"""
        Returns True if ``elt`` is in ``self``.

        EXAMPLES::

            sage: from sage.sets.integer_range import IntegerRangeFinite
            sage: I = IntegerRangeFinite(123,12,-4)
            sage: 123 in I
            True
            sage: 127 in I
            False
            sage: 12 in I
            False
            sage: 13 in I
            False
            sage: 14 in I
            False
            sage: 15 in I
            True
            sage: 11 in I
            False
        """
        if isinstance(elt, Integer):
            if (self._step.__abs__()).divides(Integer(elt)-self._begin):
                return (self._begin <= elt < self._end and self._step > 0) or \
                       (self._begin >= elt > self._end and self._step < 0)
            else:
                return False
        else:
            return False

    def cardinality(self):
        """
        Return the cardinality of ``self``

        EXAMPLES::

            sage: from sage.sets.integer_range import IntegerRangeFinite
            sage: IntegerRangeFinite(123,12,-4).cardinality()
            28
            sage: IntegerRangeFinite(-57,12,8).cardinality()
            9
            sage: IntegerRangeFinite(123,12,4).cardinality()
            0
        """
        if (self._begin > self._end and self._step > 0) or \
           (self._begin < self._end and self._step < 0):
            return 0
        return (abs((self._end+self._step-self._begin))-1) // abs(self._step)


    def __iter__(self):
        r"""
        Returns an iterator over the elements of ``self``

        EXAMPLES::

            sage: from sage.sets.integer_range import IntegerRangeFinite
            sage: I = IntegerRangeFinite(123,12,-4)
            sage: p = iter(I)
            sage: (p.next(), p.next(), p.next(), p.next(), p.next(), p.next(), p.next(), p.next())
            (123, 119, 115, 111, 107, 103, 99, 95)
            sage: I = IntegerRangeFinite(-57,12,8)
            sage: p = iter(I)
            sage: (p.next(), p.next(), p.next(), p.next(), p.next(), p.next(), p.next(), p.next())
            (-57, -49, -41, -33, -25, -17, -9, -1)
        """
        n = self._begin
        if self._step > 0:
            while n < self._end:
                yield n
                n += self._step
        else:
            while n > self._end:
                yield n
                n += self._step

    def _an_element_(self):
        r"""
        Returns an element of ``self``.

        EXAMPLES::

            sage: from sage.sets.integer_range import IntegerRangeFinite
            sage: I = IntegerRangeFinite(123,12,-4)
            sage: I.an_element()
            115
            sage: I = IntegerRangeFinite(-57,12,8)
            sage: I.an_element()
            -41
        """
        assert self.cardinality() <> 0, "FIXME (#8543): categories doesn't allow empty sets"
        p = (self._begin + 2*self._step)
        if p in self:
            return p
        else:
            return self._begin

class IntegerRangeInfinite(IntegerRange):
    r""" The class of infinite enumerated sets of integers defined by infinite
    arithmetic progressions.

    .. warning::

         To be well defined, the Integer ``begin`` has to be finite and
         ``end`` infinite.

    See :class:`IntegerRange` for more details.
    """
    def __init__(self, begin, end, step=Integer(1)):
        r"""
        TESTS::

            sage: from sage.sets.integer_range import IntegerRangeInfinite
            sage: I = IntegerRangeInfinite(-57,Infinity,8)
            sage: I.category()
            Category of infinite enumerated sets
            sage: TestSuite(I).run()
        """
        assert isinstance(begin, Integer)
        assert isinstance(end, (MinusInfinity, PlusInfinity))
        self._begin = begin
        self._end = end
        self._step = step
        Parent.__init__(self, category = InfiniteEnumeratedSets())

    def __contains__(self, elt):
        r"""
        Returns True if ``elt`` is in ``self``.

        EXAMPLES::

            sage: from sage.sets.integer_range import IntegerRangeInfinite
            sage: I = IntegerRangeInfinite(-57,Infinity,8)
            sage: -57 in I
            True
            sage: -65 in I
            False
            sage: -49 in I
            True
            sage: 743 in I
            True
        """
        if isinstance(elt, Integer):
            if (self._step.__abs__()).divides(Integer(elt)-self._begin):
                return (self._step > 0 and elt >= self._begin) or \
                       (self._step < 0 and elt <= self._begin)
            else:
                return False
        else:
            return False

    def __iter__(self):
        r"""
        Returns an iterator over the elements of ``self``.

        EXAMPLES::

            sage: from sage.sets.integer_range import IntegerRangeInfinite
            sage: I = IntegerRangeInfinite(-57,Infinity,8)
            sage: p = iter(I)
            sage: (p.next(), p.next(), p.next(), p.next(), p.next(), p.next(), p.next(), p.next())
            (-57, -49, -41, -33, -25, -17, -9, -1)
            sage: I = IntegerRangeInfinite(-112,-Infinity,-13)
            sage: p = iter(I)
            sage: (p.next(), p.next(), p.next(), p.next(), p.next(), p.next(), p.next(), p.next())
            (-112, -125, -138, -151, -164, -177, -190, -203)
        """
        n = self._begin
        while True:
            yield n
            n += self._step

    def _an_element_(self):
        r"""
        Returns an element of ``self``.

        EXAMPLES::

            sage: from sage.sets.integer_range import IntegerRangeInfinite
            sage: I = IntegerRangeInfinite(-57,Infinity,8)
            sage: I.an_element()
            191
            sage: I = IntegerRangeInfinite(-112,-Infinity,-13)
            sage: I.an_element()
            -515
        """
        return self._begin+31*self._step


class IntegerRangeFromMiddle(IntegerRange):
    r"""
    The class of finite or infinite enumerated sets defined with
    an inside point, a progression and two limits.

    See :class:`IntegerRange` for more details.
    """
    def __init__(self, begin, end, step=Integer(1), middle_point=Integer(1)):
        r"""
        TESTS::

            sage: from sage.sets.integer_range import IntegerRangeFromMiddle
            sage: I = IntegerRangeFromMiddle(-100,100,10,0)
            sage: I.category()
            Category of finite enumerated sets
            sage: TestSuite(I).run()
            sage: I = IntegerRangeFromMiddle(Infinity,-Infinity,-37,0)
            sage: I.category()
            Category of infinite enumerated sets
            sage: TestSuite(I).run()

            sage: IntegerRange(0, 5, 1, -3)
            Traceback (most recent call last):
            ...
            AssertionError: middle_point is not in the interval
        """
        self._begin = begin
        self._end = end
        self._step = step
        self._middle_point = middle_point
        assert middle_point in self, "middle_point is not in the interval"

        if (begin != Infinity and begin != -Infinity) and \
             (end != Infinity and end != -Infinity):
            Parent.__init__(self, category = FiniteEnumeratedSets())
        else:
            Parent.__init__(self, category = InfiniteEnumeratedSets())

    def __repr__(self):
        r"""
        TESTS::

            sage: from sage.sets.integer_range import IntegerRangeFromMiddle
            sage: IntegerRangeFromMiddle(Infinity,-Infinity,-37,0)
            Integer progression containing 0 with increment -37 and bounded with +Infinity and -Infinity
            sage: IntegerRangeFromMiddle(-100,100,10,0)
            Integer progression containing 0 with increment 10 and bounded with -100 and 100
        """
        return "Integer progression containing %s with increment %s and bounded with %s and %s"%(self._middle_point,self._step,self._begin,self._end)

    def __contains__(self, elt):
        r"""
        Returns True if ``elt`` is in ``self``.

        EXAMPLES::

            sage: from sage.sets.integer_range import IntegerRangeFromMiddle
            sage: I = IntegerRangeFromMiddle(-100,100,10,0)
            sage: -110 in I
            False
            sage: -100 in I
            True
            sage: 30 in I
            True
            sage: 90 in I
            True
            sage: 100 in I
            False
        """
        if isinstance(elt, Integer):
            if (self._step.__abs__()).divides(Integer(elt)-self._middle_point):
                return (self._begin <= elt and elt < self._end) or \
                       (self._begin >= elt and elt > self._end)
            else:
                return False
        else:
            return False

    def next(self, elt):
        r"""
        Return the next element of ``elt`` in ``self``.

        EXAMPLES::

            sage: from sage.sets.integer_range import IntegerRangeFromMiddle
            sage: I = IntegerRangeFromMiddle(-100,100,10,0)
            sage: (I.next(0), I.next(10), I.next(-10), I.next(20), I.next(-100))
            (10, -10, 20, -20, None)
            sage: I = IntegerRangeFromMiddle(-Infinity,Infinity,10,0)
            sage: (I.next(0), I.next(10), I.next(-10), I.next(20), I.next(-100))
            (10, -10, 20, -20, 110)
        """
        assert (elt in self)
        n = self._middle_point
        if (elt <= n and self._step > 0) or (elt >= n and self._step < 0):
            right = 2*n-elt+self._step
            if right in self:
                return right
            else:
                left = elt-self._step
                if left in self:
                    return left
        else:
            left = 2*n-elt
            if left in self:
                return left
            else:
                right = elt+self._step
                if right in self:
                    return right

    def __iter__(self):
        r"""
        Returns an iterator over the elements of ``self``.

        EXAMPLES::

            sage: from sage.sets.integer_range import IntegerRangeFromMiddle
            sage: I = IntegerRangeFromMiddle(Infinity,-Infinity,-37,0)
            sage: p = iter(I)
            sage: (p.next(), p.next(), p.next(), p.next(), p.next(), p.next(), p.next(), p.next())
            (0, -37, 37, -74, 74, -111, 111, -148)
            sage: I = IntegerRangeFromMiddle(-12,214,10,0)
            sage: p = iter(I)
            sage: (p.next(), p.next(), p.next(), p.next(), p.next(), p.next(), p.next(), p.next())
            (0, 10, -10, 20, 30, 40, 50, 60)
        """
        n = self._middle_point
        while n != None:
            yield n
            n = self.next(n)

    def _an_element_(self):
        r"""
        Returns an element of ``self``.

        EXAMPLES::

           sage: from sage.sets.integer_range import IntegerRangeFromMiddle
           sage: I = IntegerRangeFromMiddle(Infinity,-Infinity,-37,0)
           sage: I.an_element()
           0
           sage: I = IntegerRangeFromMiddle(-12,214,10,0)
           sage: I.an_element()
           0
        """
        return self._middle_point
