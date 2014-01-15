"""
Interval-posets

This module implements the combinatorial object interval-poset which 
represents an interval of the Tamari order.

**AUTHORS:**

- Viviane Pons 2014: initial implementation
"""
#*****************************************************************************
#       Copyright (C) 2010 Florent Hivert <viviane.pons@univie.ac.at>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.categories.finite_posets import FinitePosets
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.misc.cachefunc import cached_method
from sage.rings.integer import Integer
from sage.sets.non_negative_integers import NonNegativeIntegers
from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
from sage.sets.family import Family
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation

# Abstract class to serve as a Factory no instance are created.
class IntervalPosets(UniqueRepresentation, Parent):
    """
    Factory for interval-posets.

    INPUT:

    - ``size`` -- (optional) an integer

    OUPUT:

    - the set of all interval-posets (of the given ``size`` if specified)

    EXAMPLES::

        sage: IntervalPosets()
        Interval-posets

        sage: IntervalPosets(2)
        Interval-posets of size 2

    .. NOTE:: this in a factory class whose constructor returns instances of
              subclasses.
    """
    @staticmethod
    def __classcall_private__(cls, n=None):
        """
        TESTS::

            sage: from sage.combinat.interval_posets import IntervalPosets_all, IntervalPosets_size
            sage: isinstance(IntervalPosets(2), IntervalPosets)
            True
            sage: isinstance(IntervalPosets(), IntervalPosets)
            True
            sage: IntervalPosets(2) is IntervalPosets_size(2)
            True
            sage: IntervalPosets() is IntervalPosets_all()
            True
        """
        if n is None:
            return IntervalPosets_all()
        else:
            if not (isinstance(n, (Integer, int)) and n >= 0):
                raise ValueError("n must be a non negative integer")
            return IntervalPosets_size(Integer(n))
            
#################################################################
# Enumerated set of all Interval-posets
#################################################################
class IntervalPosets_all(DisjointUnionEnumeratedSets, IntervalPosets):

    def __init__(self):
        """
        TESTS::

            sage: from sage.combinat.binary_tree import BinaryTrees_all
            sage: B = BinaryTrees_all()
            sage: B.cardinality()
            +Infinity

            sage: it = iter(B)
            sage: (it.next(), it.next(), it.next(), it.next(), it.next())
            (., [., .], [., [., .]], [[., .], .], [., [., [., .]]])
            sage: it.next().parent()
            Binary trees
            sage: B([])
            [., .]

            sage: B is BinaryTrees_all()
            True
            sage: TestSuite(B).run()
            """
        DisjointUnionEnumeratedSets.__init__(
            self, Family(NonNegativeIntegers(), IntervalPosets_size),
            facade=True, keepkey = False)

    def _repr_(self):
        """
        TEST::

            sage: BinaryTrees()   # indirect doctest
            Binary trees
        """
        return "Interval-posets"
            
#################################################################
# Enumerated set of binary trees of a given size
#################################################################
class IntervalPosets_size(IntervalPosets):
    """
    The enumerated sets of interval-posets of a given size

    TESTS::

        sage: from sage.combinat.binary_tree import BinaryTrees_size
        sage: for i in range(6): TestSuite(BinaryTrees_size(i)).run()
    """
    def __init__(self, size):
        """
        TESTS::

            sage: S = BinaryTrees(3)
            sage: S == loads(dumps(S))
            True

            sage: S is BinaryTrees(3)
            True
        """
        # there is a natural order on interval-posets throught inclusions
        # that is why we use the FinitePosets category
        super(IntervalPosets_size, self).__init__(category = FinitePosets())

        self._size = size
        
    def _repr_(self):
        """
        TESTS::

            sage: BinaryTrees(3)   # indirect doctest
            Binary trees of size 3
        """
        return "Interval-posets of size %s"%(self._size)
        
    def cardinality(self):
        """
        The cardinality of ``self``

        The formula was given in [CHA]_ `\frac{2(4n+1)!}{(n+1)!(3n+2)!}

        REFERENCES:

        ..[CHA] Sur le nombre d'intervalles dans les treillis de Tamari, F. Chapoton, 2008

        EXAMPLES::

            sage: IntervalPosets(2).cardinality()
            3
            sage: IntervalPosets(3).cardinality()
            13
            sage: IntervalPosets(4).cardinality()
            68
            sage: IntervalPosets(5).cardinality()
            399
        """
        from sage.functions.other import factorial
        n = self._size
        return 2*factorial(4*n+1)/(factorial(n+1)*factorial(3*n+2))


