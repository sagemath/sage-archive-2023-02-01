"""
Tamari Interval-posets

This module implements the combinatorial object Tamari interval-poset which 
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
from sage.combinat.posets.posets import Poset 
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute, lazy_class_attribute
from sage.rings.integer import Integer
from sage.sets.non_negative_integers import NonNegativeIntegers
from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
from sage.sets.family import Family
from sage.structure.element import Element
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation

class TamariIntervalPoset(Element):
    """
    The class of Tamari Interval-posets.
    
    An interval-poset is a labelled poset of size n, with labelled `1,\dots,n`
    satisfying the following conditions:

    - if a<c (as a number) and a preceds c in the poset, then, 
      for all b such that a<b<c, b preceds c,
    
    - if a<c (as a number) and c preceds a in the poset, then,
      for all b such that a<b<c, b preceds a.

    They are in bijection with intervals of the Tamari lattice.

    INPUT:

    - ``size``, an integer, the size of the interval-posets (number of 
      vertices)

    - ``relations``, an iterable of couples (a,b) (list, tuple or iterable) 
      representing a relation 'a preceds b' in the poset. 
      
    - ``check`` (default: True) whether to check the interval-poset 
      condition or not.

    EXAMPLES::

        sage: TamariIntervalPoset(0,[])
        The tamari interval of size 0 induced by relations []
        sage: TamariIntervalPoset(3,[])
        The tamari interval of size 3 induced by relations []
        sage: TamariIntervalPoset(3,[(1,2)])
        The tamari interval of size 3 induced by relations [(1, 2)]
        sage: TamariIntervalPoset(3,[(1,2),(2,3)])
        The tamari interval of size 3 induced by relations [(1, 2), (2, 3)]
        sage: TamariIntervalPoset(3,[(1,2),(2,3),(1,3)])
        The tamari interval of size 3 induced by relations [(1, 2), (2, 3)]
        sage: TamariIntervalPoset(3,[(1,2),(3,2)])
        The tamari interval of size 3 induced by relations [(1, 2), (3, 2)]

        sage: TamariIntervalPoset(3,[(3,4)])
        Traceback (most recent call last):
        ...
        ValueError: The relations do not correspond to the size of the poset.

        sage: TamariIntervalPoset(2,[(2,1),(1,2)])
        Traceback (most recent call last):
        ...
        ValueError: Hasse diagram contains cycles.

    """

    def __init__(self, size, relations, check=True):
        """
        TESTS::

            sage: TamariIntervalPoset(3,[(1,2),(3,2)]).parent()
            Interval-posets

        """
        parent = TamariIntervalPosets()
        self._size = size
        self._poset = Poset( ([i for i in xrange(1,size+1)], relations) )
        if(self._poset.cardinality()!=size):
            raise ValueError, "The relations do not correspond to the size of the poset."%()

        Element.__init__(self, parent)

        # latex parameters initialization
        self.latex_color_decreasing = "red"
        self.latex_color_increasing = "blue"
        self.latex_hspace = 1
        self.latex_vspace = 1

    @cached_method
    def increasing_cover_relations(self):
        """
        Return the cover relations of the increasing poset of ``self`` (the poset
        formed by keeping only relations a preceded b with a<b)

        EXAMPLES::

            sage: TamariIntervalPoset(4,[(1,2),(3,2),(2,4),(3,4)]).increasing_cover_relations()
            [(1, 2), (2, 4), (3, 4)]
            sage: TamariIntervalPoset(3,[(1,2),(1,3),(2,3)]).increasing_cover_relations()
            [(1, 2), (2, 3)]

        """
        relations = []
        for i in xrange(1,self.size()):
            for j in xrange(i+1, self.size()+1):
                if self.le(i,j):
                    relations.append((i,j))
                    break
        return relations
    
    @cached_method
    def decreasing_cover_relations(self):
        """
        Return the cover relations of the increasing poset of ``self`` (the poset
        formed by keeping only relations a preceded b with a<b)

        EXAMPLES::
        
            sage: TamariIntervalPoset(4,[(2,1),(3,2),(3,4),(4,2)]).decreasing_cover_relations()
            [(4, 2), (3, 2), (2, 1)]
            sage: TamariIntervalPoset(4,[(2,1),(4,3),(2,3)]).decreasing_cover_relations()
            [(4, 3), (2, 1)]
            sage: TamariIntervalPoset(3,[(2,1),(3,1),(3,2)]).decreasing_cover_relations()
            [(3, 2), (2, 1)]
        """
        relations = []
        for i in xrange(self.size(),1,-1):
            for j in xrange(i-1,0,-1):
                if self.le(i,j):
                    relations.append((i,j))
                    break
        return relations

    def le(self, e1, e2):
        """
        Return whether ``e1`` precedes or equals ``e2`` in ``self``

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(1,2),(2,3)])
            sage: ip.le(1,2)
            True
            sage: ip.le(1,3)
            True
            sage: ip.le(2,3)
            True
            sage: ip.le(3,4)
            False
            sage: ip.le(1,1)
            True

        """
        return self._poset.le(e1,e2)
        
    def lt(self, e1, e2):
        """
        Return whether ``e1`` strictly precedes ``e2`` in ``self``

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(1,2),(2,3)])
            sage: ip.lt(1,2)
            True
            sage: ip.lt(1,3)
            True
            sage: ip.lt(2,3)
            True
            sage: ip.lt(3,4)
            False
            sage: ip.lt(1,1)
            False

        """
        return self._poset.lt(e1,e2)
        
    def ge(self, e1, e2):
        """
        Return whether ``e2`` precedes or equals ``e1`` in ``self``

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(1,2),(2,3)])
            sage: ip.ge(2,1)
            True
            sage: ip.ge(3,1)
            True
            sage: ip.ge(3,2)
            True
            sage: ip.ge(4,3)
            False
            sage: ip.ge(1,1)
            True
        """
        return self._poset.ge(e1, e2)
        
    def gt(self, e1, e2):
        """
        Return whether ``e2`` strictly precedes ``e1`` in ``self``

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(1,2),(2,3)])
            sage: ip.gt(2,1)
            True
            sage: ip.gt(3,1)
            True
            sage: ip.gt(3,2)
            True
            sage: ip.gt(4,3)
            False
            sage: ip.gt(1,1)
            False
        """
        return self._poset.gt(e1,e2)
        
    def size(self):
        """
        Return the size (number of vertices) of the interval-poset.

        EXAMPLES::

            sage: TamariIntervalPoset(3,[(2,1),(3,1)]).size()
            3
        """
        return self._size
    
    def _repr_(self):
        """
        TESTS::

            sage: TamariIntervalPoset(3,[(2,1),(3,1)])
            The tamari interval of size 3 induced by relations [(3, 1), (2, 1)]
            sage: TamariIntervalPoset(3,[(3,1),(2,1)])
            The tamari interval of size 3 induced by relations [(3, 1), (2, 1)]
            sage: TamariIntervalPoset(3,[(2,3),(2,1)])
            The tamari interval of size 3 induced by relations [(2, 3), (2, 1)]
        """
        return "The tamari interval of size %s induced by relations %s"%(self.size(),str(self.increasing_cover_relations() + self.decreasing_cover_relations()))

# Abstract class to serve as a Factory no instance are created.
class TamariIntervalPosets(UniqueRepresentation, Parent):
    """
    Factory for interval-posets.

    INPUT:

    - ``size`` -- (optional) an integer

    OUPUT:

    - the set of all interval-posets (of the given ``size`` if specified)

    EXAMPLES::

        sage: TamariIntervalPosets()
        Interval-posets

        sage: TamariIntervalPosets(2)
        Interval-posets of size 2

    .. NOTE:: this in a factory class whose constructor returns instances of
              subclasses.
    """
    @staticmethod
    def __classcall_private__(cls, n=None):
        """
        TESTS::

            sage: from sage.combinat.interval_posets import TamariIntervalPosets_all, TamariIntervalPosets_size
            sage: isinstance(TamariIntervalPosets(2), TamariIntervalPosets)
            True
            sage: isinstance(TamariIntervalPosets(), TamariIntervalPosets)
            True
            sage: TamariIntervalPosets(2) is TamariIntervalPosets_size(2)
            True
            sage: TamariIntervalPosets() is TamariIntervalPosets_all()
            True
        """
        if n is None:
            return TamariIntervalPosets_all()
        else:
            if not (isinstance(n, (Integer, int)) and n >= 0):
                raise ValueError("n must be a non negative integer")
            return TamariIntervalPosets_size(Integer(n))
            
#################################################################
# Enumerated set of all Tamari Interval-posets
#################################################################
class TamariIntervalPosets_all(DisjointUnionEnumeratedSets, TamariIntervalPosets):

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
            self, Family(NonNegativeIntegers(), TamariIntervalPosets_size),
            facade=True, keepkey = False)

    def _repr_(self):
        """
        TEST::

            sage: BinaryTrees()   # indirect doctest
            Binary trees
        """
        return "Interval-posets"
        
    def _element_constructor_(self, size, relations):
        """
            EXAMPLES::

                sage: TIP = TamariIntervalPosets()
                sage: TIP(3,[(1,2)]) # indirect doctest
                The tamari interval of size 3 induced by relations [(1, 2)]
        """
        return self.element_class(size,relations) 

    Element = TamariIntervalPoset

#################################################################
# Enumerated set of Tamari interval-posets of a given size
#################################################################
class TamariIntervalPosets_size(TamariIntervalPosets):
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
        super(TamariIntervalPosets_size, self).__init__(category = FinitePosets())

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

            sage: TamariIntervalPosets(2).cardinality()
            3
            sage: TamariIntervalPosets(3).cardinality()
            13
            sage: TamariIntervalPosets(4).cardinality()
            68
            sage: TamariIntervalPosets(5).cardinality()
            399
        """
        from sage.functions.other import factorial
        n = self._size
        return 2*factorial(4*n+1)/(factorial(n+1)*factorial(3*n+2))
        
    @lazy_attribute
    def _parent_for(self):
        """
        The parent of the element generated by ``self``

        TESTS::

            sage: TIP3 = TamariIntervalPosets(3)
            sage: TIP3._parent_for
            Interval-posets
        """
        return TamariIntervalPosets_all()

    @lazy_attribute
    def element_class(self):
        """
        TESTS::

            sage: S = BinaryTrees(3)
            sage: S.element_class
            <class 'sage.combinat.binary_tree.BinaryTrees_all_with_category.element_class'>
            sage: S.first().__class__ == BinaryTrees().first().__class__
            True
        """
        return self._parent_for.element_class

    def _element_constructor_(self, relations):
        """
        EXAMPLES::

            sage: TIP3 = TamariIntervalPosets(3)
            sage: TIP3([(1,2)]) # indirect doctest
            The tamari interval of size 3 induced by relations [(1, 2)]
            sage: TIP3([(3,4)])
            Traceback (most recent call last):
            ...
            ValueError: The relations do not correspond to the size of the poset.

        """
        return self.element_class(self._size,relations) 

