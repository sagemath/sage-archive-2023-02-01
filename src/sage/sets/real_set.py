# -*- coding: utf-8 -*-
"""
Subsets of the Real Line

This module contains subsets of the real line that can be constructed
as the union of a finite set of open and closed intervals.

EXAMPLES::

    sage: RealSet(0,1)
    (0, 1)
    sage: RealSet((0,1), [2,3])
    (0, 1) + [2, 3]
    sage: RealSet(-oo, oo)
    (-Infinity, +Infinity)

Brackets must be balanced in Python, so the naive notation for
half-open intervals does not work::

    sage: RealSet([0,1))
    Traceback (most recent call last):
    ...
    SyntaxError: invalid syntax
    
Instead, you can use the following construction functions::

    sage: RealSet.open_closed(0,1)
    (0, 1]
    sage: RealSet.closed_open(0,1)
    [0, 1)
    sage: RealSet.point(1/2)
    {1/2}
    sage: RealSet.unbounded_below_open(0)
    (-Infinity, 0)
    sage: RealSet.unbounded_below_closed(0)
    (-Infinity, 0]
    sage: RealSet.unbounded_above_open(1)
    (1, +Infinity)
    sage: RealSet.unbounded_above_closed(1)
    [1, +Infinity)
"""

########################################################################
#       Copyright (C) 2013 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#   as published by the Free Software Foundation; either version 2 of
#   the License, or (at your option) any later version.
#                   http://www.gnu.org/licenses/
########################################################################

from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.real_lazy import LazyFieldElement, RLF
from sage.rings.infinity import infinity, minus_infinity


class RealInterval(UniqueRepresentation, Parent):

    def __init__(self, lower, lower_closed, upper, upper_closed, check=True):
        """
        A real interval

        You are not supposed to create :class:`RealInterval` objects
        yourself. Always use :class:`RealSet` instead.

        EXAMPLES::

            sage: RealSet([0, oo])
            Traceback (most recent call last):
            ...
            ValueError: interval cannot be closed at +oo
        """
        self._lower = lower
        self._upper = upper
        self._lower_closed = lower_closed
        self._upper_closed = upper_closed
        if check:
            if not isinstance(lower, LazyFieldElement):
                raise ValueError('lower bound must be in RLF')
            if not isinstance(upper, LazyFieldElement):
                raise ValueError('upper bound must be in RLF')
            if not isinstance(lower_closed, bool):
                raise ValueError('lower_closed must be boolean')
            if not isinstance(upper_closed, bool):
                raise ValueError('upper_closed must be boolean')
            if lower == infinity or upper == minus_infinity:
                raise ValueError('only lower/upper bound can be minus/plus infinity')
            if lower > upper:
                raise ValueError('lower/upper bounds are not sorted')
            if (lower_closed and lower == minus_infinity):
                raise ValueError('interval cannot be closed at -oo')
            if (upper_closed and upper == infinity):
                raise ValueError('interval cannot be closed at +oo')
            
    def is_empty(self):
        """
        Return whether the interval is empty
        
        The normalized form of :class:`RealSet` has all intervals
        non-empty, so this method usually returns ``False``.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: I = RealSet(0, 1)[0]
            sage: I.is_empty()
            False
        """
        return (self._lower == self._upper) and not (self._lower_closed or self._upper_closed)

    def is_point(self):
        """
        Return whether the interval consists of a single point

        OUTPUT:
        
        Boolean.

        EXAMPLES::

            sage: I = RealSet(0, 1)[0]
            sage: I.is_point()
            False
        """
        return (self._lower == self._upper) and self._lower_closed and self._upper_closed

    def lower(self):
        """
        Return the lower bound
        
        OUTPUT:

        The lower bound as it was originally specified.
        
        EXAMPLES::

            sage: I = RealSet(0, 1)[0]
            sage: I.lower()
            0
            sage: I.upper()
            1
        """
        return self._lower._value

    def upper(self):
        """
        Return the upper bound

        OUTPUT:

        The upper bound as it was originally specified.
        
        EXAMPLES::

            sage: I = RealSet(0, 1)[0]
            sage: I.lower()
            0
            sage: I.upper()
            1
        """
        return self._upper._value
                
    def lower_closed(self):
        """
        Return whether the interval is open at the lower bound
        
        OUTPUT:

        Boolean.
        
        EXAMPLES::

            sage: I = RealSet.open_closed(0, 1)[0];  I
            (0, 1]
            sage: I.lower_closed()
            False
            sage: I.lower_open()
            True
            sage: I.upper_closed()
            True
            sage: I.upper_open()
            False
        """
        return self._lower_closed

    def upper_closed(self):
        """
        Return whether the interval is closed at the lower bound

        OUTPUT:

        Boolean.
        
        EXAMPLES::

            sage: I = RealSet.open_closed(0, 1)[0];  I
            (0, 1]
            sage: I.lower_closed()
            False
            sage: I.lower_open()
            True
            sage: I.upper_closed()
            True
            sage: I.upper_open()
            False
        """
        return self._upper_closed
                
    def lower_open(self):
        """
        Return whether the interval is closed at the upper bound

        OUTPUT:

        Boolean.
        
        EXAMPLES::

            sage: I = RealSet.open_closed(0, 1)[0];  I
            (0, 1]
            sage: I.lower_closed()
            False
            sage: I.lower_open()
            True
            sage: I.upper_closed()
            True
            sage: I.upper_open()
            False
        """
        return not self._lower_closed

    def upper_open(self):
        """
        Return whether the interval is closed at the upper bound

        OUTPUT:

        Boolean.
        
        EXAMPLES::

            sage: I = RealSet.open_closed(0, 1)[0];  I
            (0, 1]
            sage: I.lower_closed()
            False
            sage: I.lower_open()
            True
            sage: I.upper_closed()
            True
            sage: I.upper_open()
            False
        """
        return not self._upper_closed
                
    def __cmp__(self, other):
        """
        Intervals are sorted by lower bound, then upper bound

        OUTPUT:

        `-1`, `0`, or `+1` depending on how the intervals compare.
        
        EXAMPLES::

             sage: I1 = RealSet.open_closed(1, 3)[0];  I1
             (1, 3]
             sage: I2 = RealSet.open_closed(0, 5)[0];  I2
             (0, 5]
             sage: cmp(I1, I2)
             1
             sage: sorted([I1, I2])
             [(0, 5], (1, 3]]
        """
        return cmp([self._lower, self._lower_closed, self._upper, self._upper_closed],
                   [other._lower, other._lower_closed, other._upper, other._upper_closed])
        
    element_class = LazyFieldElement

    def _repr_(self):
        """
        Return a string representation

        OUTPUT:

        String.
        
        EXAMPLES::

            sage: RealSet.open_closed(0, 1)
            (0, 1]
            sage: RealSet.point(0)
            {0}
        """
        if self.is_point():
            return '{' + str(self._lower._value) + '}'
        s =  '[' if self._lower_closed else '('
        s += str(self._lower._value)
        s += ', '
        s += str(self._upper._value)
        s +=  ']' if self._upper_closed else ')'
        return s

    def closure(self):
        """
        Return the closure

        OUTPUT:

        The closure as a new :class:`RealInterval`

        EXAMPLES:
        
            sage: RealSet.open(0,1)[0].closure()
            [0, 1]
            sage: RealSet.open(-oo,1)[0].closure()
            (-Infinity, 1]
            sage: RealSet.open(0, oo)[0].closure()
            [0, +Infinity)
        """
        lower_closed = (self._lower != minus_infinity)
        upper_closed = (self._upper != infinity)
        return RealInterval(self._lower, lower_closed, self._upper, upper_closed)
        
    def interior(self):
        """
        Return the interior

        OUTPUT:

        The interior as a new :class:`RealInterval`

        EXAMPLES:
        
            sage: RealSet.closed(0, 1)[0].interior()
            (0, 1)
            sage: RealSet.open_closed(-oo, 1)[0].interior()
            (-Infinity, 1)
            sage: RealSet.closed_open(0, oo)[0].interior()
            (0, +Infinity)
        """
        return RealInterval(self._lower, False, self._upper, False)
        
    def is_connected(self, other):
        """
        Test whether two intervals are connected
        
        OUTPUT:

        Boolean. Whether the set-theoretic union of the two intervals
        has a single connected component.

        EXAMPLES::
        
            sage: I1 = RealSet.open(0, 1)[0];  I1
            (0, 1)
            sage: I2 = RealSet.closed(1, 2)[0];  I2
            [1, 2]
            sage: I1.is_connected(I2)
            True
            sage: I1.is_connected(I2.interior())
            False
            sage: I1.closure().is_connected(I2.interior())
            True
            sage: I2.is_connected(I1)
            True
            sage: I2.interior().is_connected(I1)
            False
            sage: I2.closure().is_connected(I1.interior())
            True
            sage: I3 = RealSet.closed(1/2, 3/2)[0]; I3
            [1/2, 3/2]
            sage: I1.is_connected(I3)
            True
            sage: I3.is_connected(I1)
            True
        """
        cmp_lu = cmp(self._lower, other._upper)
        cmp_ul = cmp(self._upper, other._lower)
        # self is seperated and below other
        if cmp_ul == -1:
            return False
        # self is adjacent and below other 
        if cmp_ul == 0:
            return self._upper_closed or other._lower_closed
        # self is seperated and above other
        if cmp_lu == +1:
            return False
        # self is adjacent and above other 
        if cmp_lu == 0:
            return self._lower_closed or other._upper_closed
        # They are not separated
        return True
        
    def convex_hull(self, other):
        """
        Return the convex hull of the two intervals

        OUTPUT:

        The convex hull as a new :class:`RealInterval`.

        EXAMPLES::

            sage: I1 = RealSet.open(0, 1)[0];  I1
            (0, 1)
            sage: I2 = RealSet.closed(1, 2)[0];  I2
            [1, 2]
            sage: I1.convex_hull(I2)
            (0, 2]
            sage: I2.convex_hull(I1)
            (0, 2]
            sage: I1.convex_hull(I2.interior())
            (0, 2)
            sage: I1.closure().convex_hull(I2.interior())
            [0, 2)
            sage: I1.closure().convex_hull(I2)
            [0, 2]
            sage: I3 = RealSet.closed(1/2, 3/2)[0]; I3
            [1/2, 3/2]
            sage: I1.convex_hull(I3)
            (0, 3/2]
        """
        cmp_ll = cmp(self._lower, other._lower)
        cmp_uu = cmp(self._upper, other._upper)
        lower = upper = None
        lower_closed = upper_closed = None
        if cmp_ll == -1:
            lower = self._lower
            lower_closed = self._lower_closed
        elif cmp_ll == +1:
            lower = other._lower
            lower_closed = other._lower_closed
        else:
            assert(cmp_ll == 0)
            lower = self._lower
            lower_closed = self._lower_closed or other._lower_closed
        if cmp_uu == +1:
            upper = self._upper
            upper_closed = self._upper_closed
        elif cmp_uu == -1:
            upper = other._upper
            upper_closed = other._upper_closed
        else:
            assert(cmp_uu == 0)
            upper = self._upper
            upper_closed = self._upper_closed or other._upper_closed
        return RealInterval(lower, lower_closed, upper, upper_closed)

    def intersection(self, other):
        """
        Return the intersection of the two intervals
        
        INPUT:
        
        - ``other`` -- a :class:`RealInterval`

        OUTPUT:

        The intersection as a new :class:`RealInterval`

        EXAMPLES::

            sage: I1 = RealSet.open(0, 2)[0];  I1
            (0, 2)
            sage: I2 = RealSet.closed(1, 3)[0];  I2
            [1, 3]
            sage: I1.intersection(I2)
            [1, 2)
            sage: I2.intersection(I1)
            [1, 2)
            sage: I1.closure().intersection(I2.interior())
            (1, 2]
            sage: I2.interior().intersection(I1.closure())
            (1, 2]

            sage: I3 = RealSet.closed(10, 11)[0];  I3
            [10, 11]
            sage: I1.intersection(I3)
            (0, 0)
            sage: I3.intersection(I1)
            (0, 0)
        """
        cmp_ll = cmp(self._lower, other._lower)
        cmp_uu = cmp(self._upper, other._upper)
        lower = upper = None
        lower_closed = upper_closed = None
        if cmp_ll == -1:
            lower = other._lower
            lower_closed = other._lower_closed
        elif cmp_ll == +1:
            lower = self._lower
            lower_closed = self._lower_closed
        else:
            assert(cmp_ll == 0)
            lower = self._lower
            lower_closed = self._lower_closed and other._lower_closed
        if cmp_uu == +1:
            upper = other._upper
            upper_closed = other._upper_closed
        elif cmp_uu == -1:
            upper = self._upper
            upper_closed = self._upper_closed
        else:
            assert(cmp_uu == 0)
            upper = self._upper
            upper_closed = self._upper_closed and other._upper_closed
        if lower > upper:
            lower = upper = RLF(0)
            lower_closed = upper_closed = False
        return RealInterval(lower, lower_closed, upper, upper_closed)
        


class RealSet(UniqueRepresentation, Parent):

    @staticmethod
    def __classcall__(cls, *args):
        """
        Normalize the input.

        INPUT:
        
        See :class:`RealSet`.
          
        OUTPUT:

        A :class:`RealSet`.

        EXAMPLES::

            sage: RealSet(RealSet.open_closed(0,1), RealSet.closed_open(2,3))
            (0, 1] + [2, 3)
        """
        intervals = []
        if len(args) == 2:
            # allow RealSet(0,1) interval constructor
            try:
                lower, upper = args
                lower.n()
                upper.n()
                args = (RealSet._prep(lower, upper), )
            except (AttributeError, ValueError, TypeError):
                pass
        for arg in args:
            if isinstance(arg, tuple):
                lower, upper = RealSet._prep(*arg)
                intervals.append(RealInterval(lower, False, upper, False))
            elif isinstance(arg, list):
                lower, upper = RealSet._prep(*arg)
                intervals.append(RealInterval(lower, True, upper, True))
            elif isinstance(arg, RealInterval):
                intervals.append(arg)
            elif isinstance(arg, RealSet):
                intervals.extend(arg._intervals)
            else:
                raise ValueError(str(arg) + ' does not determine real interval')
        intervals = RealSet.normalize(intervals)
        return UniqueRepresentation.__classcall__(cls, intervals)
                
    def __init__(self, intervals):
        """
        A subset of the real line

        INPUT:

        Arguments defining a real set. Possibilities are either two
        real numbers to construct an open set or a list/tuple/iterable
        of intervals. The individual intervals can be specified by
        either a :class:`RealInterval`, a tuple of two real numbers
        (constructing an open interval), or a list of two number
        (constructing a closed interval).

        EXAMPLES::

            sage: RealSet(0,1)    # open set from two numbers
            (0, 1)
            sage: i = RealSet(0,1)[0]
            sage: RealSet(i)      # interval
            (0, 1)
            sage: RealSet(i, (3,4))    # tuple of two numbers = open set
            (0, 1) + (3, 4)
            sage: RealSet(i, [3,4])    # list of two numbers = closed set
            (0, 1) + [3, 4]
        """
        self._intervals = intervals
    
    def __iter__(self):
        """
        Iterate over the component intervals is ascending order
        
        OUTPUT:

        An iterator over the intervals.

        EXAMPLES::

            sage: s = RealSet(RealSet.open_closed(0,1), RealSet.closed_open(2,3))
            sage: i = iter(s)
            sage: i.next()
            (0, 1]
            sage: i.next()
            [2, 3)
        """
        return iter(self._intervals)

    def n_components(self):
        """
        Return the number of connected components

        See also :meth:`get_interval`

        EXAMPLES::

            sage: s = RealSet(RealSet.open_closed(0,1), RealSet.closed_open(2,3))
            sage: s.n_components()
            2
            sage: len(s)   # shorthand
            2
        """
        return len(self._intervals)

    __len__ = n_components

    def get_interval(self, i):
        """
        Return the ``i``-th connected component.

        Note that the intervals representing the real set are always
        normalized, see :meth:`normalize`.

        INPUT:
        
        - ``i`` -- integer.

        OUTPUT:

        The $i$-th connected component as a :class:`RealInterval`.

        EXAMPLES::
        
            sage: s = RealSet(RealSet.open_closed(0,1), RealSet.closed_open(2,3))
            sage: s.get_interval(0)
            (0, 1]
            sage: s[0]    # shorthand
            (0, 1]
            sage: s.get_interval(1)
            [2, 3)
            sage: s[0] == s.get_interval(0)
            True
        """
        return self._intervals[i]

    __getitem__ = get_interval

    @staticmethod
    def normalize(intervals):
        """
        Bring a collection of intervals into canonical form

        INPUT:

        - ``intervals`` -- a list/tuple/iterable of intervals.

        OUTPUT:
        
        A tuple of intervals such that

        * they are sorted in ascending order (by lower bound)
        
        * there is a gap between each interval
        
        * all intervals are non-empty

        EXAMPLES::

            sage: i1 = RealSet((0, 1))[0]
            sage: i2 = RealSet([1, 2])[0]
            sage: i3 = RealSet((2, 3))[0]
            sage: RealSet.normalize([i1, i2, i3])
            ((0, 3),)

            sage: RealSet((0, 1), [1, 2], (2, 3))
            (0, 3)
            sage: RealSet((0, 1), (1, 2), (2, 3))
            (0, 1) + (1, 2) + (2, 3)
            sage: RealSet([0, 1], [2, 3])
            [0, 1] + [2, 3]
            sage: RealSet((0, 2), (1, 3))
            (0, 3)
            sage: RealSet(0,0)
            {}
        """
        # sort by lower bound
        intervals = sorted(intervals)
        if len(intervals) == 0:
            return tuple()
        merged = []
        curr = intervals.pop(0)
        while len(intervals) != 0:
            next = intervals.pop(0)
            cmp_ul = cmp(curr._upper, next._lower)
            if cmp_ul == +1 or (
                cmp_ul == 0 and (curr._upper_closed or next._lower_closed)):
                curr = curr.convex_hull(next)
            else:
                if not curr.is_empty():
                    merged.append(curr)
                curr = next
        if not curr.is_empty():
            merged.append(curr)
        return tuple(merged)

    def _repr_(self):
        """
        Return a string representation
        
        OUTPUT:

        A string representation.

        EXAMPLES::

            sage: RealSet(0, 1)._repr_()
            '(0, 1)'
        """
        if self.n_components() == 0:
            return '{}'
        else:
            # Switch to u'\u222A' (cup sign) with Python 3
            return ' + '.join(map(repr, self._intervals))

    @staticmethod
    def _prep(lower, upper):
        """
        Helper to prepare the lower and upper bound

        EXAMPLES::

            sage: RealSet._prep(1, 0)
            (0, 1)
        """
        lower = RLF(lower)
        upper = RLF(upper)
        if upper < lower:
            return upper, lower
        else:
            return lower, upper

    @staticmethod
    def open(lower, upper):
        """
        Construct an open interval

        INPUT:

        - ``lower``, ``upper`` -- two real numbers or infinity. They
          will be sorted if necessary.

        OUTPUT:

        A new :class:`RealSet`.

        EXAMPLES::

            sage: RealSet.open(1, 0)
            (0, 1)
        """
        lower, upper = RealSet._prep(lower, upper)
        return RealSet(RealInterval(lower, False, upper, False))

    @staticmethod
    def closed(lower, upper):
        """
        Construct a closed interval

        INPUT:

        - ``lower``, ``upper`` -- two real numbers or infinity. They
          will be sorted if necessary.

        OUTPUT:

        A new :class:`RealSet`.

        EXAMPLES::

            sage: RealSet.closed(1, 0)
            [0, 1]
        """
        lower, upper = RealSet._prep(lower, upper)
        return RealSet(RealInterval(lower, True, upper, True))

    @staticmethod
    def point(p):
        """
        Construct an interval containing a single point

        INPUT:

        - ``p``, ``upper`` -- a real number.

        OUTPUT:

        A new :class:`RealSet`.

        EXAMPLES::

            sage: RealSet.open(1, 0)
            (0, 1)
        """
        p = RLF(p)
        return RealSet(RealInterval(p, True, p, True))
    
    @staticmethod
    def open_closed(lower, upper):
        """
        Construct a half-open interval

        INPUT:

        - ``lower``, ``upper`` -- two real numbers or infinity. They
          will be sorted if necessary.

        OUTPUT:

        A new :class:`RealSet` that is open at the lower bound and
        closed at the upper bound.

        EXAMPLES::

            sage: RealSet.open_closed(1, 0)
            (0, 1]
        """
        lower, upper = RealSet._prep(lower, upper)
        return RealSet(RealInterval(lower, False, upper, True))

    @staticmethod
    def closed_open(lower, upper):
        """
        Construct an half-open interval

        INPUT:

        - ``lower``, ``upper`` -- two real numbers or infinity. They
          will be sorted if necessary.

        OUTPUT:

        A new :class:`RealSet` that is closed at the lower bound and
        open an the upper bound.

        EXAMPLES::

            sage: RealSet.closed_open(1, 0)
            [0, 1)
        """
        lower, upper = RealSet._prep(lower, upper)
        return RealSet(RealInterval(lower, True, upper, False))

    @staticmethod
    def unbounded_below_closed(bound):
        """
        Construct a semi-infinite interval

        INPUT:

        - ``bound`` -- a real number.

        OUTPUT:

        A new :class:`RealSet` from minus infinity to the bound (including).

        EXAMPLES::

            sage: RealSet.unbounded_below_closed(1)
            (-Infinity, 1]
        """
        return RealSet(RealInterval(RLF(minus_infinity), False, RLF(bound), True))

    @staticmethod
    def unbounded_below_open(bound):
        """
        Construct a semi-infinite interval

        INPUT:

        - ``bound`` -- a real number.

        OUTPUT:

        A new :class:`RealSet` from minus infinity to the bound (excluding).

        EXAMPLES::

            sage: RealSet.unbounded_below_open(1)
            (-Infinity, 1)
        """
        return RealSet(RealInterval(RLF(minus_infinity), False, RLF(bound), False))

    @staticmethod
    def unbounded_above_closed(bound):
        """
        Construct a semi-infinite interval

        INPUT:

        - ``bound`` -- a real number.

        OUTPUT:

        A new :class:`RealSet` from the bound (including) to plus
        infinity.

        EXAMPLES::

            sage: RealSet.unbounded_above_closed(1)
            [1, +Infinity)
        """
        return RealSet(RealInterval(RLF(bound), True, RLF(infinity), False))

    @staticmethod
    def unbounded_above_open(bound):
        """
        Construct a semi-infinite interval

        INPUT:

        - ``bound`` -- a real number.

        OUTPUT:

        A new :class:`RealSet` from the bound (excluding) to plus
        infinity.

        EXAMPLES::

            sage: RealSet.unbounded_above_open(1)
            (1, +Infinity)
        """
        return RealSet(RealInterval(RLF(bound), False, RLF(infinity), False))

    def union(self, other):
        """
        Return the union of the two sets

        INPUT:
        
        - ``other`` -- a :class:`RealSet`.

        OUTPUT:
        
        The set-theoretic union as a new :class:`RealSet`.

        EXAMPLES::

            sage: s1 = RealSet(0,2)
            sage: s2 = RealSet(1,3)
            sage: s1.union(s2)
            (0, 3)
            sage: s1 | s2    # syntactic sugar
            (0, 3)
            sage: s1 + s2    # syntactic sugar
            (0, 3)
        """
        intervals = self._intervals + other._intervals
        return RealSet(*intervals)
    
    __or__ = union
    __add__ = union

    def intersection(self, other):
        """
        Return the intersection of the two sets

        INPUT:
        
        - ``other`` -- a :class:`RealSet`.

        OUTPUT:
        
        The set-theoretic intersection as a new :class:`RealSet`.

        EXAMPLES::

            sage: s1 = RealSet(0,2) + RealSet.unbounded_above_closed(10);  s1
            (0, 2) + [10, +Infinity)
            sage: s2 = RealSet(1,3) + RealSet.unbounded_below_closed(-10);  s2
            (-Infinity, -10] + (1, 3)
            sage: s1.intersection(s2)
            (1, 2)
            sage: s1 & s2    # syntactic sugar
            (1, 2)
        """
        # TODO: this can be done in linear time since the intervals are already sorted
        intervals = []
        for i1 in self._intervals:
            for i2 in other._intervals:
                intervals.append(i1.intersection(i2))
        return RealSet(*intervals)

    __and__ = intersection

    def inf(self):
        """
        Return the infimum

        OUTPUT:

        A real number or infinity.

        EXAMPLES::

            sage: s1 = RealSet(0,2) + RealSet.unbounded_above_closed(10);  s1
            (0, 2) + [10, +Infinity)
            sage: s1.inf()
            0

            sage: s2 = RealSet(1,3) + RealSet.unbounded_below_closed(-10);  s2
            (-Infinity, -10] + (1, 3)
            sage: s2.inf()
            -Infinity
        """
        if self.n_components() == 0:
            return infinity
        return self._intervals[0].lower()

    def sup(self):
        """
        Return the supremum

        OUTPUT:

        A real number or infinity.

        EXAMPLES::

            sage: s1 = RealSet(0,2) + RealSet.unbounded_above_closed(10);  s1
            (0, 2) + [10, +Infinity)
            sage: s1.sup()
            +Infinity

            sage: s2 = RealSet(1,3) + RealSet.unbounded_below_closed(-10);  s2
            (-Infinity, -10] + (1, 3)
            sage: s2.sup()
            3
        """
        if self.n_components() == 0:
            return minus_infinity
        return self._intervals[-1].upper()

    def complement(self):
        """
        Return the complement

        OUTPUT:

        The set-theoretic complement as a new :class:`RealSet`.
        
        EXAMPLES::

            sage: RealSet(0,1).complement()
            (-Infinity, 0] + [1, +Infinity)
       
            sage: s1 = RealSet(0,2) + RealSet.unbounded_above_closed(10);  s1
            (0, 2) + [10, +Infinity)
            sage: s1.complement()
            (-Infinity, 0] + [2, 10)

            sage: s2 = RealSet(1,3) + RealSet.unbounded_below_closed(-10);  s2
            (-Infinity, -10] + (1, 3)
            sage: s2.complement()
            (-10, 1] + [3, +Infinity)
        """
        n = self.n_components()
        if n == 0:
            return RealSet(minus_infinity, infinity)
        intervals = []
        if self.inf() != minus_infinity:
            first = self._intervals[0]
            intervals.append(RealInterval(RLF(minus_infinity), False, 
                                          first._lower, first.lower_open()))
        if self.sup() != infinity:
            last = self._intervals[-1]
            intervals.append(RealInterval(last._upper, last.upper_open(), 
                                          RLF(infinity), False))
        for i in range(1,n):
            prev = self._intervals[i-1]
            next = self._intervals[i]
            i = RealInterval(prev._upper, prev.upper_open(),
                             next._lower, next.lower_open())
            intervals.append(i)
        return RealSet(*intervals)
                             
    def difference(self, other):
        """
        Return ``self`` with ``other`` subtracted

        INPUT:
        
        - ``other`` -- a :class:`RealSet`.

        OUTPUT:
        
        The set-theoretic difference of ``self`` with ``other``
        removed as a new :class:`RealSet`.

        OUTPUT:

        A new :class:`RealSet`

            sage: s1 = RealSet(0,2) + RealSet.unbounded_above_closed(10);  s1
            (0, 2) + [10, +Infinity)
            sage: s2 = RealSet(1,3) + RealSet.unbounded_below_closed(-10);  s2
            (-Infinity, -10] + (1, 3)
            sage: s1.difference(s2)
            (0, 1] + [10, +Infinity)
            sage: s1 - s2    # syntactic sugar
            (0, 1] + [10, +Infinity)
            sage: s2.difference(s1)
            (-Infinity, -10] + [2, 3)
            sage: s2 - s1    # syntactic sugar
            (-Infinity, -10] + [2, 3)
        """
        return self.intersection(other.complement())

    __sub__ = difference
