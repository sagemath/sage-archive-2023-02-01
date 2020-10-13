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

TESTS::

    sage: pt = path_tableaux.SemistandardPath([[],[3],[3,2],[3,3,1],[3,3,2,1],[4,3,3,1,0]])
    sage: TestSuite(pt).run()

    sage: pt = path_tableaux.SemistandardPath([[],[3],[3,2],[3,3,1],[7/2,3,2,1],[4,3,3,1,0]])
    sage: TestSuite(pt).run()

    sage: pt = path_tableaux.SemistandardPath([[3,2],[3,3,1],[3,3,2,1],[4,3,3,1,0]])
    sage: pt.promotion()
    Traceback (most recent call last):
    ...
    IndexError: tuple index out of range

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
from sage.combinat.skew_tableau import SkewTableau
from sage.combinat.tableau import Tableau
from sage.combinat.gelfand_tsetlin_patterns import GelfandTsetlinPattern

###############################################################################

class SemistandardPath(PathTableau):
    r"""
    An instance is the sequence of partitions which is the
    chain of partitions of a skew semistandard tableau.

    The :class:'SemistandardSkewTableau` is not implemented.

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
        [(2,), (2, 1)]

        sage: path_tableaux.SemistandardPath([[],[5/2],[7/2,2]])
        [(), (5/2,), (7/2, 2)]

        sage: path_tableaux.SemistandardPath([[],[2.5],[3.5,2]])
        [(), (2.50000000000000,), (3.50000000000000, 2)]
    """

    @staticmethod
    def __classcall_private__(cls, st, check=True):
        r"""
        This ensures that a tableau is only ever constructed as an
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
            ValueError: [(), (3/2,), (2, 5/2)] does not satisfy the required inequalities

            sage: path_tableaux.SemistandardPath([(), 3, (3, 2)])
            Traceback (most recent call last):
            ...
            ValueError: [(), 3, (3, 2)] is not a sequence of lists
        """
        w = None

        if isinstance(st, GelfandTsetlinPattern):
            w = [tuple(a) for a in st]
            w.reverse()
            w = [()] + w
            w = tuple(w)

        elif isinstance(st, Tableau):
            w = [tuple(p) for p in st.to_chain()]
            w = tuple(w)

        elif isinstance(st, SkewTableau):
            w = [tuple(p) for p in st.to_chain()]
            w = w[1:]
            w = tuple(w)

        elif isinstance(st, (list,tuple)):
            try:
                w = tuple([tuple(a) for a in st])
            except TypeError:
                raise ValueError(f"{st} is not a sequence of lists")

        if w is None:
            raise ValueError(f"invalid input {st}")

        PathTableau.__init__(self, parent, w, check=check)

    def check(self):
        """
        Checks that ``self`` is a valid path.

        TESTS::

            sage: path_tableaux.SemistandardPath([[],[2],[1,2]])
            Traceback (most recent call last):
            ...
            ValueError: [(), (2,), (1, 2)] does not satisfy the required inequalities

            sage: path_tableaux.SemistandardPath([[],[2],[1,2]],check=False)
            [(), (2,), (1, 2)]
        """
        if not all( self[i+1][j] >= self[i][j] >= self[i+1][j+1]
                    for i in range(len(self)-1) for j in range(len(self[i+1])-1) ):
            raise ValueError(f"{self} does not satisfy the required inequalities")

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
        Return `True` if `self` is skew.

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
        if not 0 < i < self.size()-1:
            raise ValueError(f"{i} is not defined on {self}")

        #st = self.to_tableau()
        #result = st.bender_knuth_involution(i)
        gt = self.to_pattern()
        result = gt.bender_knuth_involution(i)
        return SemistandardPath(result)


    @combinatorial_map(name='to semistandard tableau')
    def to_tableau(self):
        r"""
        Converts ``self`` to a :class:`SemistandardTableau`.

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
        from sage.combinat.partition import Partition

        parts = [Partition(a) for a in self]
        return from_chain(parts)

    @combinatorial_map(name='to Gelfand-Tsetlin pattern')
    def to_pattern(self):
        r"""
        Converts ``self`` to a Gelfand-Tsetlin pattern.

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

        m = len(lt[0])
        pt = [list(a)+[0]*(m+i-len(a)) for i,a in enumerate(lt)]
        pt.reverse()
        return GelfandTsetlinPattern(pt)

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
