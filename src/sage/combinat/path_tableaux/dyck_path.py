r"""
Dyck Paths

This is an implementation of the abstract base class
:class:`sage.combinat.pathtableau.pathtableaux`.
This is the simplest implementation of PathTableaux and is included to
provide a convenient test case and for pedagogical purposes.

In this implementation we have sequences of nonnegative integers. These
are required to be the heights Dyck words (except that we do not require
the sequence to start or end at height zero). These are in bijection
with skew standard tableaux with at most two rows. Sequences which start
and end at height zero are in bijection with noncrossing perfect matchings.

AUTHORS:

- Bruce Westbury (2018): initial version

Here we illustrate the slogan that promotion = rotation.

EXAMPLES::

    sage: t = DyckPath([0,1,2,3,2,1,0])
    sage: t.to_perfect_matching()
    [(0, 5), (1, 4), (2, 3)]

    sage: t = t.promotion()
    sage: t.to_perfect_matching()
    [(0, 3), (1, 2), (4, 5)]

    sage: t = t.promotion()
    sage: t.to_perfect_matching()
    [(0, 1), (2, 5), (3, 4)]

    sage: t = t.promotion()
    sage: t.to_perfect_matching()
    [(0, 5), (1, 4), (2, 3)]

    sage: TestSuite(t).run()
"""

#*****************************************************************************
#       Copyright (C) 2018 Bruce Westbury <bruce.westbury@gmail.com>,
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.list_clone import ClonableArray
from sage.combinat.path_tableaux.path_tableau import PathTableau, PathTableaux
from sage.combinat.combinatorial_map import combinatorial_map
from sage.combinat.dyck_word import DyckWord
from sage.combinat.perfect_matching import PerfectMatching
from sage.combinat.skew_tableau import SkewTableau
from sage.combinat.tableau import Tableau, StandardTableau
from sage.rings.integer import Integer

###############################################################################

class DyckPath(PathTableau):
    r"""
    An instance is the sequence of nonnegative
    integers given by the heights of a Dyck word.

    INPUT:

    * a sequence of nonnegative integers
    * a two row standard skew tableau
    * a Dyck word
    * a noncrossing perfect matching

    EXAMPLES::

        sage: DyckPath([0,1,2,1,0])
        [0, 1, 2, 1, 0]

        sage: w = DyckWord([1,1,0,0])
        sage: DyckPath(w)
        [0, 1, 2, 1, 0]

        sage: p = PerfectMatching([(1,2),(3,4)])
        sage: DyckPath(p)
        [0, 1, 0, 1, 0]

        sage: t = Tableau([[1,2],[3,4]])
        sage: DyckPath(t)
        [0, 1, 2, 1, 0]
    """

    @staticmethod
    def __classcall_private__(cls, ot):
        r"""
        This ensures that a tableau is only ever constructed as an
        ``element_class`` call of an appropriate parent.

        TESTS::

            sage: t = DyckPath([0,1,2,1,0])

            sage: t.parent()
            <sage.combinat.path_tableaux.catalan.DyckPaths_with_category object at ...>
            sage: t.category()
            Category of elements of <sage.combinat.path_tableaux.catalan.DyckPaths_with_category object at ...>
            sage: type(t)
            <class 'sage.combinat.path_tableaux.catalan.DyckPaths_with_category.element_class'>
        """
        return DyckPaths()(ot)

    def __init__(self, parent, ot, check=True):
        r"""
        Initialize a Catalan tableau.

        INPUT:

        Can be any of:

        * word of nonnegative integers with successive differences '\pm 1'
        * a DyckWord
        * a noncrossing perfect matching
        * a two row standard tableau
        * a two row standard skew tab;eau

        TESTS::

            sage: DyckPath([0,1,2,1,0])
            [0, 1, 2, 1, 0]
            sage: DyckPath(DyckWord([1, 0, 1, 0]))
            [0, 1, 0, 1, 0]
            sage: DyckPath(PerfectMatching([(1, 4), (2, 3), (5, 6)]))
            [0, 1, 2, 1, 0, 1, 0]
            sage: DyckPath(Tableau([[1,2,4],[3,5,6]]))
            [0, 1, 2, 1, 2, 1, 0]
            sage: DyckPath(SkewTableau([[None, 1,4],[2,3]]))
            [1, 2, 1, 0, 1]
            sage: DyckPath(PerfectMatching([(1, 3), (2, 4), (5, 6)]))
            Traceback (most recent call last):
            ...
            ValueError: the perfect matching must be non crossing
            sage: DyckPath(Tableau([[1,2,5],[3,5,6]]))
            Traceback (most recent call last):
            ...
            ValueError: the tableau must be standard
            sage: DyckPath(Tableau([[1,2,4],[3,5,6],[7]]))
            Traceback (most recent call last):
            ...
            ValueError: the tableau must have at most two rows
            sage: DyckPath(SkewTableau([[None, 1,4],[2,3],[5]]))
            Traceback (most recent call last):
            ...
            ValueError: the skew tableau must have at most two rows
            sage: DyckPath([0,1,2.5,1,0])
            Traceback (most recent call last):
            ...
            ValueError: [0, 1, 2.50000000000000, 1, 0] is not a sequence of integers
            sage: DyckPath(Partition([3,2,1]))
            Traceback (most recent call last):
            ...
            ValueError: invalid input [3, 2, 1]
        """
        w = None

        if isinstance(ot, DyckWord):
            w = ot.heights()

        elif isinstance(ot, PerfectMatching):
            if ot.is_noncrossing():
                u = [1]*ot.size()
                for a in ot.arcs():
                    u[a[1]-1] = 0
                w = DyckWord(u).heights()
            else:
                raise ValueError("the perfect matching must be non crossing")

        elif isinstance(ot, Tableau):
            if len(ot) > 2:
                raise ValueError("the tableau must have at most two rows")
            if ot.is_standard():
                u = [1] * ot.size()
                for i in ot[1]:
                    u[i-1] = 0
                w = DyckWord(u).heights()
            else:
                raise ValueError("the tableau must be standard")


        elif isinstance(ot, SkewTableau):
            if len(ot) > 2:
                raise ValueError("the skew tableau must have at most two rows")
            # The check that ot is standard is not implemented
            c = ot.to_chain()
            w = [0]*len(c)
            for i,a in enumerate(c):
                if len(a) == 1:
                    w[i] = a[0]
                else:
                    w[i] = a[0]-a[1]

        elif isinstance(ot, (list,tuple)):
            try:
                w = tuple([Integer(a) for a in ot])
            except TypeError:
                raise ValueError("%s is not a sequence of integers" % ot)

        if w is None:
            raise ValueError("invalid input %s" % ot)

        ClonableArray.__init__(self, parent, w, check=check)

    def check(self):
        """ Checks that ``self`` is a valid path.

        TESTS::

            sage: DyckPath([0,1,0,-1,0]) # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: [0, 1, 0, -1, 0] has a negative entry

            sage: DyckPath([0,1,3,1,0]) # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: [0, 1, 3, 1, 0] is not a Dyck path
        """
        n = len(self)
        if any(a < 0 for a in self):
           raise ValueError( "%s has a negative entry" % (str(self)) )
        for i in range(n-1):
            if abs(self[i+1]-self[i]) != 1:
                raise ValueError( "%s is not a Dyck path" % str(self) )

    def _local_rule(self,i):
        """
        This has input a list of objects. This method first takes
        the list of objects of length three consisting of the `(i-1)`-st,
        `i`-th and `(i+1)`-term and applies the rule. It then replaces
        the `i`-th object  by the object returned by the rule.

        EXAMPLES::

            sage: t = DyckPath([0,1,2,3,2,1,0])
            sage: t._local_rule(3)
            [0, 1, 2, 1, 2, 1, 0]

        TESTS::

            sage: t = DyckPath([0,1,2,3,2,1,0])
            sage: t._local_rule(0)
            Traceback (most recent call last):
            ...
            ValueError: 0 is not a valid integer
            sage: t._local_rule(5)
            [0, 1, 2, 3, 2, 1, 0]
            sage: t._local_rule(6)
            Traceback (most recent call last):
            ...
            ValueError: 6 is not a valid integer
        """

        def _rule(x):
            """
            This is the rule on a sequence of three letters.
            """
            return abs(x[0]-x[1]+x[2])

        if not (i > 0 and i < len(self)-1):
            raise ValueError("%d is not a valid integer" % i)

        with self.clone() as result:
            result[i] = _rule(self[i-1:i+2])

        return result

    def is_skew(self):
        r"""
        Return ``True`` if ``self`` is skew and ``False`` if not.

        EXAMPLES::

            sage: DyckPath([0,1,2,1]).is_skew()
            False

            sage: DyckPath([1,0,1,2,1]).is_skew()
            True
        """
        return self[0] != 0

    @combinatorial_map(name='to Dyck word')
    def to_DyckWord(self):
        r"""
        Converts ``self`` to a Dyck word.

        EXAMPLES::

            sage: c = DyckPath([0,1,2,1,0])
            sage: c.to_DyckWord()
            [1, 1, 0, 0]
        """
        return DyckWord(heights_sequence = list(self))

    def descents(self):
        r"""
        Return the descent set of ``self``.

        EXAMPLES::

            sage: DyckPath([0,1,2,1,2,1,0,1,0]).descents()
            {3, 6}
        """
        result = set()

        for i in range(1,len(self)-1):
            if self[i] < self[i-1] and self[i] < self[i+1]:
                result.add(i)

        return result

    def to_word(self):
        r"""
        Return the word in the alphabet `\{0,1\}` associated to ``self``.

        EXAMPLES::

            sage: DyckPath([1,0,1,2,1]).to_word()
            [0, 1, 1, 0]
        """
        return [ (self[i+1]-self[i]+1)/2 for i in range(self.size()-1) ]

    def to_perfect_matching(self):
        r"""
        Return the perfect matching associated to ``self``.

        EXAMPLES::

            sage: DyckPath([0,1,2,1,2,1,0,1,0]).to_perfect_matching()
            [(0, 5), (1, 2), (3, 4), (6, 7)]

        TESTS::

            sage: DyckPath([1,2,1,2,1,0,1]).to_perfect_matching()
            Traceback (most recent call last):
            ...
            ValueError: [1, 2, 1, 2, 1, 0, 1] does not start at 0
        """
        if self.is_skew():
            raise ValueError( "%s does not start at 0" % (str(self)) )
        w = self.to_word()
        y = DyckWord(w)
        pairs = set()
        for i, a in enumerate(y):
            c = y.associated_parenthesis(i)
            if i < c:
                pairs.add((i,c))
        return PerfectMatching(pairs)

    def to_tableau(self):
        r"""
        Return the skew tableau associated to ``self``.

        EXAMPLES::

            sage: T = DyckPath([0,1,2,3,2,3])
            sage: T.to_tableau()
            [[1, 2, 3, 5], [4]]

            sage: U = DyckPath([2,3,2,3])
            sage: U.to_tableau()
            [[None, None, 1, 3], [2]]
        """
        w = self.to_word()
        top = [ i+1 for i, a in enumerate(w) if a == 1 ]
        bot = [ i+1 for i, a in enumerate(w) if a == 0 ]
        if self.is_skew():
            return SkewTableau([[None]*self[0]+top,bot])
        else:
            return StandardTableau([top,bot])

class DyckPaths(PathTableaux):
    """
    The parent class for DyckPath.
    """
    Element = DyckPath

