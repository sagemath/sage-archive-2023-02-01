r"""
Semistandard Tableaux
=====================

This is an implementation of the abstract base class
:class:`sage.combinat.path_tableaux.path_tableau.PathTableau`.

This implementation is for semistandard tableaux, represented as a chain of partitions
(essentially, the Gelfand-Tsetlin pattern).
This generalises the jeu-de-taquin operations of rectification, promotion, evacuation from
standard tableaux to semistandard tableaux. The local rule is the Bender-Knuth involution.

EXAMPLES::

    sage: pt = path_tableaux.SemistandardPath([[],[3],[3,2],[3,3,1],[3,3,2,1],[4,3,3,1,0]])
    sage: pt.promotion()
    [(), (2,), (3, 1), (3, 2, 1), (4, 3, 1, 0), (4, 3, 3, 1, 0)]
    sage: pt.evacuation()
    [(), (2,), (4, 0), (4, 2, 0), (4, 3, 1, 0), (4, 3, 3, 1, 0)]

    sage: pt = path_tableaux.SemistandardPath([[],[3],[3,2],[3,3,1],[3,3,2,1],[9/2,3,3,1,0]])
    sage: pt.promotion()
    [(), (2,), (3, 1), (3, 2, 1), (9/2, 3, 1, 0), (9/2, 3, 3, 1, 0)]
    sage: pt.evacuation()
    [(), (5/2,), (9/2, 0), (9/2, 2, 0), (9/2, 3, 1, 0), (9/2, 3, 3, 1, 0)]

    sage: pt = path_tableaux.SemistandardPath([[],[3],[4,2],[5,4,1]])
    sage: path_tableaux.CylindricalDiagram(pt)
    [       (),      (3,),    (4, 2), (5, 4, 1)]
    [         ,        (),      (3,),    (5, 2), (5, 4, 1)]
    [         ,          ,        (),      (4,),    (4, 3), (5, 4, 1)]
    [         ,          ,          ,        (),      (3,),    (5, 1), (5, 4, 1)]

    sage: pt2 = path_tableaux.SemistandardPath([[3,2],[3,3,1],[3,3,2,1],[4,3,3,1,0]])
    sage: pt1 = path_tableaux.SemistandardPath([[],[3],[3,2]])
    sage: pt1.commutor(pt2)
    ([(), (2,), (2, 2), (4, 2, 0)], [(4, 2, 0), (4, 3, 2, 0), (4, 3, 3, 1, 0)])
    sage: pt1.commutor(pt2,verbose=True)
    [(3, 2), (3, 3, 1), (3, 3, 2, 1), (4, 3, 3, 1, 0)]
    [(3,), (3, 2), (3, 2, 2), (4, 3, 2, 0)]
    [(), (2,), (2, 2), (4, 2, 0)]
    ([(), (2,), (2, 2), (4, 2, 0)], [(4, 2, 0), (4, 3, 2, 0), (4, 3, 3, 1, 0)])

    sage: st = SkewTableau([[None,None,None,4,4,5,6,7],[None,2,4,6,7,7,7],[None,4,5,8,8,9],[None,6,7,10],[None,8,8,11],[None],[4]])
    sage: pt = path_tableaux.SemistandardPath(st)
    sage: bk = [SkewTableau(st.bender_knuth_involution(i+1)) for i in range(10)]
    sage: lr = [pt.local_rule(i+1) for i in range(10)]
    sage: all(r.to_tableau() == s for r,s in zip(lr,bk))
    True

TESTS::

    sage: pt = path_tableaux.SemistandardPath([[],[3],[3,2],[3,3,1],[3,3,2,1],[4,3,3,1,0]])
    sage: TestSuite(pt).run()

    sage: pt = path_tableaux.SemistandardPath([[],[3],[3,2],[3,3,1],[7/2,3,2,1],[4,3,3,1,0]])
    sage: TestSuite(pt).run()

    sage: pt = path_tableaux.SemistandardPath([[3,2],[3,3,1],[3,3,2,1],[4,3,3,1,0]])
    sage: pt.promotion()
    [(3, 2), (3, 2, 2), (4, 3, 2, 0), (4, 3, 3, 1, 0)]

AUTHORS:

- Bruce Westbury (2020): initial version
"""

#*****************************************************************************
#       Copyright (C) 2020 Bruce Westbury <bruce.westbury@gmail.com>,
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.combinat.path_tableaux.path_tableau import PathTableau, PathTableaux
from sage.combinat.combinatorial_map import combinatorial_map
from sage.combinat.skew_tableau import SkewTableau, SkewTableaux
from sage.combinat.tableau import Tableau
from sage.combinat.gelfand_tsetlin_patterns import GelfandTsetlinPattern
from sage.combinat.partition import Partition

###############################################################################

class SemistandardPath(PathTableau):
    r"""
    An instance is the sequence of partitions which is the
    chain of partitions of a skew semistandard tableau.

    The :class:`SemistandardSkewTableau is not implemented.

    INPUT:

    * a sequence of partitions
    * a sequence of lists/tuples
    * a semistandard tableau
    * a Gelfand-Tsetlin pattern

    EXAMPLES::

        sage: path_tableaux.SemistandardPath([[],[2],[2,1]])
        [(), (2,), (2, 1)]

        sage: gt = GelfandTsetlinPattern([[2,1],[2]])
        sage: path_tableaux.SemistandardPath(gt)
        [(), (2,), (2, 1)]

        sage: st = SemistandardTableau([[1,1],[2]])
        sage: path_tableaux.SemistandardPath(st)
        [(), (2,), (2, 1)]

        sage: st = SkewTableau([[1,1],[2]])
        sage: path_tableaux.SemistandardPath(st)
        [(), (2,), (2, 1)]

        sage: st = SkewTableau([[None,1,1],[2]])
        sage: path_tableaux.SemistandardPath(st)
        [(1,), (3, 0), (3, 1, 0)]

        sage: path_tableaux.SemistandardPath([[],[5/2],[7/2,2]])
        [(), (5/2,), (7/2, 2)]

        sage: path_tableaux.SemistandardPath([[],[2.5],[3.5,2]])
        [(), (2.50000000000000,), (3.50000000000000, 2)]
    """

    @staticmethod
    def __classcall_private__(cls, st, check=True):
        r"""
        Ensure that a tableau is only ever constructed as an
        ``element_class`` call of an appropriate parent.

        EXAMPLES::

            sage: t = path_tableaux.SemistandardPath([[],[2]])
            sage: t.parent()
            <sage.combinat.path_tableaux.semistandard.SemistandardPaths_with_category object at ...>
        """
        return SemistandardPaths()(st, check=check)

    def __init__(self, parent, st, check=True):
        r"""
        Initialize a semistandard tableau.

        Don't require entries to be positive or to be integral.
        The operations we require are addition, subtraction and comparisons.

        TESTS::

            sage: path_tableaux.SemistandardPath([[],[3/2],[2,5/2]])
            Traceback (most recent call last):
            ...
            ValueError: [(), (3/2,), (2, 5/2)] does not satisfy the required inequalities in row 1

            sage: path_tableaux.SemistandardPath([(), 3, (3, 2)])
            Traceback (most recent call last):
            ...
            ValueError: [(), 3, (3, 2)] is not a sequence of lists
        """
        w = None

        if isinstance(st, SemistandardPath):
            w = list(st)

        elif isinstance(st, GelfandTsetlinPattern):
            w = list(st)
            w.reverse()
            w = [()] + w

        elif isinstance(st, (Tableau,SkewTableau)):
            w = st.to_chain()

        elif isinstance(st, (list,tuple)):
            if any(not isinstance(a,(list,tuple)) for a in st):
                raise ValueError(f"{st} is not a sequence of lists")
            w = st

        else:
            raise ValueError(f"invalid input {st} is of type {type(st)}")

        # Pad with zeroes, if necessary
        m = max(len(a)-i for i,a in enumerate(w))
        w = [list(a)+[0]*(m+i-len(a)) for i,a in enumerate(w)]
        # Convert to immutable
        w = tuple([tuple(a) for a in w])

        PathTableau.__init__(self, parent, w, check=check)

    def check(self):
        """
        Check that ``self`` is a valid path.

        TESTS::

            sage: path_tableaux.SemistandardPath([[],[2],[1,2]])
            Traceback (most recent call last):
            ...
            ValueError: [(), (2,), (1, 2)] does not satisfy the required inequalities in row 1

            sage: path_tableaux.SemistandardPath([[],[2],[1,2]],check=False)
            [(), (2,), (1, 2)]
        """
        for i in range(1,len(self)-1):
            if not all(r>=s for r,s in zip(self[i+1],self[i])):
                raise ValueError(f"{self} does not satisfy the required inequalities in row {i}")
            if not all(r>=s for r,s in zip(self[i],self[i+1][1:])):
                raise ValueError(f"{self} does not satisfy the required inequalities in row {i}")

    def size(self):
        r"""
        Return the size or length of ``self``.

        EXAMPLES::

            sage: path_tableaux.SemistandardPath([[],[3],[3,2],[3,3,1],[3,3,2,1]]).size()
            5
        """
        return len(self)

    def is_skew(self):
        """
        Return ``True`` if ``self`` is skew.

        EXAMPLES::

            sage: path_tableaux.SemistandardPath([[],[2]]).is_skew()
            False
            sage: path_tableaux.SemistandardPath([[2,1]]).is_skew()
            True
        """
        return self[0] != ()

    def local_rule(self,i):
        r"""
        This has input a list of objects. This method first takes
        the list of objects of length three consisting of the `(i-1)`-st,
        `i`-th and `(i+1)`-term and applies the rule. It then replaces
        the `i`-th object  by the object returned by the rule.

        This is the Bender-Knuth involution.

        EXAMPLES::

            sage: pt = path_tableaux.SemistandardPath([[],[3],[3,2],[3,3,1],[3,3,2,1]])
            sage: pt.local_rule(1)
            [(), (2,), (3, 2), (3, 3, 1), (3, 3, 2, 1)]
            sage: pt.local_rule(2)
            [(), (3,), (3, 2), (3, 3, 1), (3, 3, 2, 1)]
            sage: pt.local_rule(3)
            [(), (3,), (3, 2), (3, 2, 2), (3, 3, 2, 1)]

        TESTS::

            sage: pt = path_tableaux.SemistandardPath([[],[3],[3,2],[3,3,1],[3,3,2,1]])
            sage: pt.local_rule(0)
            Traceback (most recent call last):
            ...
            ValueError: 0 is not defined on [(), (3,), (3, 2), (3, 3, 1), (3, 3, 2, 1)]
            sage: pt.local_rule(4)
            Traceback (most recent call last):
            ...
            ValueError: 4 is not defined on [(), (3,), (3, 2), (3, 3, 1), (3, 3, 2, 1)]
        """
        def toggle(i,j):
            """
            Return the toggle of entry 'self[i][j]'.
            """

            if j == 0:
                left = self[i+1][0]
            else:
                left = min(self[i+1][j], self[i-1][j-1])
            if j == len(self[i])-1:
                right = self[i+1][j+1]
            else:
                right = max(self[i+1][j+1], self[i-1][j])

            return left + right - self[i][j]

        if not 0 < i < self.size()-1:
            raise ValueError(f"{i} is not defined on {self}")

        with self.clone() as result:
            result[i] = tuple([toggle(i,k) for k in range(len(self[i]))])

        return result

    def rectify(self,inner=None,verbose=False):
        """
        Rectify ``self``.

        EXAMPLES::

            sage: st = SkewTableau([[None, None, None, 4],[None,None,1,6],[None,None,5],[2,3]])
            sage: path_tableaux.SemistandardPath(st).rectify()
            [(), (1,), (1, 1), (2, 1, 0), (3, 1, 0, 0), (3, 2, 0, 0, 0), (4, 2, 0, 0, 0, 0)]
            sage: path_tableaux.SemistandardPath(st.rectify())
            [(), (1,), (1, 1), (2, 1, 0), (3, 1, 0, 0), (3, 2, 0, 0, 0), (4, 2, 0, 0, 0, 0)]
            sage: path_tableaux.SemistandardPath(st).rectify(verbose=True)
            [[(3, 2, 2), (3, 3, 2, 0), (3, 3, 2, 1, 0), (3, 3, 2, 2, 0, 0), (4, 3, 2, 2, 0, 0, 0), (4, 3, 3, 2, 0, 0, 0, 0), (4, 4, 3, 2, 0, 0, 0, 0, 0)],
            [(3, 2), (3, 3, 0), (3, 3, 1, 0), (3, 3, 2, 0, 0), (4, 3, 2, 0, 0, 0), (4, 3, 3, 0, 0, 0, 0), (4, 4, 3, 0, 0, 0, 0, 0)],
            [(3,), (3, 1), (3, 1, 1), (3, 2, 1, 0), (4, 2, 1, 0, 0), (4, 3, 1, 0, 0, 0), (4, 4, 1, 0, 0, 0, 0)],
            [(), (1,), (1, 1), (2, 1, 0), (3, 1, 0, 0), (3, 2, 0, 0, 0), (4, 2, 0, 0, 0, 0)]]

        TESTS::

            sage: S = SemistandardSkewTableaux([[5,3,3],[3,1]],[3,2,2])
            sage: LHS = [path_tableaux.SemistandardPath(st.rectify()) for st in S]
            sage: RHS = [path_tableaux.SemistandardPath(st).rectify() for st in S]
            sage: all(r==s for r,s in zip(LHS,RHS) if r != s)
            True

            sage: st = SkewTableau([[None, None, None, 4],[None,None,1,6],[None,None,5],[2,3]])
            sage: pt = path_tableaux.SemistandardPath(st)
            sage: SP = [path_tableaux.SemistandardPath(it) for it in StandardTableaux([3,2,2])]
            sage: set(pt.rectify(inner=ip) for ip in SP)
            {[(), (1,), (1, 1), (2, 1, 0), (3, 1, 0, 0), (3, 2, 0, 0, 0), (4, 2, 0, 0, 0, 0)]}
        """
        if not self.is_skew():
            return self

        n = len(self)
        pp = self[0]

        if inner == None:
            initial = [pp[:r] for r in range(len(pp))]
        elif Partition(inner[-1]) == Partition(pp):
            initial = list(inner)[:-1]
        else:
            raise ValueError(f"The final shape{inner[-1]} must agree with the initial shape {pp}")

        r = len(initial)
        path = SemistandardPath(initial + list(self))
        if verbose:
            rect = [self]

        for i in range(r):
            for j in range(n-1):
                path = path.local_rule(r+j-i)
            if verbose:
                rect.append(SemistandardPath(list(path)[r-i-1:r+n-i-1]))

        if verbose:
            return rect
        else:
            return SemistandardPath(list(path)[:n])

    @combinatorial_map(name='to semistandard tableau')
    def to_tableau(self):
        r"""
        Convert ``self`` to a :class:`SemistandardTableau`.

        The :class:`SemistandardSkewTableau` is not implemented.

        EXAMPLES::

            sage: pt = path_tableaux.SemistandardPath([[],[3],[3,2],[3,3,1],[3,3,2,1],[4,3,3,1,0]])
            sage: pt.to_tableau()
            [[1, 1, 1, 5], [2, 2, 3], [3, 4, 5], [4]]

        TESTS::

            sage: SST = SemistandardTableaux(shape=[5,5,3],eval=[2,2,3,4,2])
            sage: all(st == path_tableaux.SemistandardPath(st).to_tableau() for st in SST)
            True
        """
        from sage.combinat.tableau import from_chain

        parts = [Partition(a) for a in self]
        if self.is_skew():
            return SkewTableaux().from_chain(parts)
        else:
            return from_chain(parts)


    @combinatorial_map(name='to Gelfand-Tsetlin pattern')
    def to_pattern(self):
        r"""
        Convert ``self`` to a Gelfand-Tsetlin pattern.

        EXAMPLES::

            sage: pt = path_tableaux.SemistandardPath([[],[3],[3,2],[3,3,1],[3,3,2,1],[4,3,3,1]])
            sage: pt.to_pattern()
            [[4, 3, 3, 1, 0], [3, 3, 2, 1], [3, 3, 1], [3, 2], [3]]

        TESTS::

            sage: GT = GelfandTsetlinPatterns(top_row=[5,5,3])
            sage: all(gt == path_tableaux.SemistandardPath(gt).to_pattern() for gt in GT)
            True

            sage: GT = GelfandTsetlinPatterns(top_row=[5,5,3])
            sage: all(gt.to_tableau() == path_tableaux.SemistandardPath(gt).to_tableau() for gt in GT)
            True
        """
        lt = list(self)
        if lt[0] == ():
            lt = lt[1:]
        lt.reverse()
        return GelfandTsetlinPattern([list(a) for a in lt])

    def _test_jdt_promotion(self, **options):
        """
        Check that promotion agrees with :meth:`Tableau.promotion_inverse`
        constructed using jeu de taquin.

        TESTS::

            sage: pt = path_tableaux.SemistandardPath([(),(1,),(2,1),(4,2),(4,3,1),(4,3,3)])
            sage: pt._test_jdt_promotion()

            sage: pt = path_tableaux.SemistandardPath([(),(1,),(2,1),(4,2),(4,3,1),(9/2,3,3)])
            sage: pt._test_jdt_promotion()
        """
        tester = self._tester(**options)
        try:
            LHS = self.promotion().to_tableau()
            RHS = self.to_tableau().promotion_inverse(len(self)-2)
            tester.assertEqual(LHS,RHS)
        except ValueError:
            pass

class SemistandardPaths(PathTableaux):
    """
    The parent class for :class:`SemistandardTableau`.
    """

    def _an_element_(self):
        """
        Return an element of ``self``.

        EXAMPLES::

            sage: path_tableaux.SemistandardPaths()._an_element_()
            [(), (2,), (2, 1)]
        """
        return SemistandardPath([[],[2],[2,1]])

    Element = SemistandardPath
