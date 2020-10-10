r"""
Semistandard Tableaux
=====================

This is an implementation of the abstract base class
:class:`sage.combinat.path_tableaux.path_tableau.PathTableau`.

This implementation is for semistandard tableaux, represented as a chain of partitions
(essentially, the Gelfand-Tsetlin pattern).
This generalises the jeu-de-taquin operations of rectification, promotion, evacuation from
standard tableaux to semistandard tableaux. The local rule is the Bender-Knuth involution.

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
from sage.combinat.partition import Partition
from sage.combinat.skew_tableau import SkewTableau, SkewTableaux, SemistandardSkewTableaux
from sage.combinat.tableau import Tableau, StandardTableau
from sage.combinat.gelfand_tsetlin_patterns import GelfandTsetlinPattern
from itertools import zip_longest

###############################################################################

class SemistandardPath(PathTableau):
    r"""
    An instance is the sequence of partitions which is the
    chain of partitions of a skew semistandard tableau.

    The :class:'SemistandardSkewTableau` is not implemented.

    INPUT:

    * a sequence of partitions
    * a semistandard tableau
    * a skew semistandard tableau    
    * a Gelfand-Tsetlin pattern

    """

    @staticmethod
    def __classcall_private__(cls, st):
        r"""
        This ensures that a tableau is only ever constructed as an
        ``element_class`` call of an appropriate parent.

        TESTS::

            sage: t = path_tableaux.SemistandardPath([[],[2]])
            sage: t.parent()
            <sage.combinat.path_tableaux.semistandard.SemistandardPaths_with_category object at ...>
        """
        return SemistandardPaths()(st)

    def __init__(self, parent, st, check=True):
        r"""
        Initialize a semistandard tableau.

        TESTS::

            sage: path_tableaux.SemistandardPath([[],[2],[2,1]])
            [[], [2], [2, 1]]

            sage: gt = GelfandTsetlinPattern([[2,1],[2],[]])
            sage: path_tableaux.SemistandardPath(gt)
            [[], [2], [2, 1]]

            sage: st = SemistandardTableau([[1,1],[2]])
            sage: path_tableaux.SemistandardPath(st)
            [[], [2], [2, 1]]

            sage: st = SkewTableau([[1,1],[2]])
            sage: path_tableaux.SemistandardPath(st)
            [[], [2], [2, 1]]

            sage: path_tableaux.SemistandardPath([[],[2],[1,2]])
            Traceback (most recent call last):
            ...
            ValueError: [1, 2] is not an element of Partitions
        """
        w = None

        if isinstance(st, GelfandTsetlinPattern):
            w = [Partition(p) for p in st]
            w.reverse()
            w = tuple(w)

        elif isinstance(st, Tableau):
            w = [Partition(p) for p in st.to_chain()]
            w = tuple(w)

        elif isinstance(st, SkewTableau):
            w = [Partition(p) for p in st.to_chain()]
            w = tuple(w)

        elif isinstance(st, (list,tuple)):
            try:
                w = tuple([Partition(a) for a in st])
            except TypeError:
                raise ValueError(f"{st} is not a sequence of partitions")

        if w is None:
            raise ValueError(f"invalid input {st}")

        PathTableau.__init__(self, parent, w, check=check)

    def check(self):
        """
        Checks that ``self`` is a valid path.

        TESTS::

        """        
        for x, y in zip(self,self[1:]):
            if not y.contains(x):
                raise ValueError(f"{y} does not contain {x}")
            for u, v in zip_longest(x,y[1:],fillvalue=0):
                if v > u:
                    print(u,v)
                    raise ValueError(f"the skew partition {y}\{x} is not a horizontal strip")

    def is_skew(self):
        """
        Return `True` if `self` is skew.

        EXAMPLES::

            sage: path_tableaux.SemistandardPath([[],[2]]).is_skew()
            False
            sage: path_tableaux.SemistandardPath([[2],[2,1]]).is_skew()
            True
        """
        return self[0] != Partition([])

    def local_rule(self,i):
        r"""
        This has input a list of objects. This method first takes
        the list of objects of length three consisting of the `(i-1)`-st,
        `i`-th and `(i+1)`-term and applies the rule. It then replaces
        the `i`-th object  by the object returned by the rule.

        This is the Bender-Knuth involution.

        EXAMPLES::

        TESTS::

        """
        gt = self.to_pattern()
        result = gt.bender_knuth_involution(i)
        return SemistandardPath(result)


    @combinatorial_map(name='to semistandard tableau')
    def to_tableau(self):
        r"""
        Converts ``self`` to a :class:`SemistandardTableau`.

        The :class:`SemistandardSkewTableau` is not implemented.

        EXAMPLES::

            sage: pt = path_tableaux.SemistandardPath([[],[3],[3,2],[3,3,1],[3,3,2,1],[4,3,3,1]])
            sage: pt.to_tableau()
            [[1, 1, 1, 5], [2, 2, 3], [3, 4, 5], [4]]
        """
        from sage.combinat.tableau import from_chain
        if self.is_skew():
            return SkewTableaux().from_chain(self)
        else:
            return from_chain(self)

    @combinatorial_map(name='to Gelfand-Tsetlin pattern')
    def to_pattern(self):
        r"""
        Converts ``self`` to a Gelfand-Tsetlin pattern.

        EXAMPLES::

            sage: pt = path_tableaux.SemistandardPath([[],[3],[3,2],[3,3,1],[3,3,2,1],[4,3,3,1]])
            sage: pt.to_pattern()
            [[4, 3, 3, 1, 0], [3, 3, 2, 1], [3, 3, 1], [3, 2], [3], []]
        """
        m = len(self[0])
        pt = [list(a)+[0]*(m+i-len(a)) for i,a in enumerate(self)]
        pt.reverse()
        return GelfandTsetlinPattern(pt)

class SemistandardPaths(PathTableaux):
    """
    The parent class for :class:`SemistandardTableau`.
    """

    def _an_element_(self):
        """
        Return an element of ``self``.

        EXAMPLES::

            sage: path_tableaux.SemistandardPaths()._an_element_()
            [[], [2], [2, 1]]

        """
        return SemistandardPath([[],[2],[2,1]])

    Element = SemistandardPath
