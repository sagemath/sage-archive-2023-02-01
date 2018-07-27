r"""
Path Tableaux

This is an abstract base class for using local rules to construct rectification
and the action of the cactus group. This is an effective version
of the Henriques-Kamnitzer construction of the action of the cactus
group on tensor powers of a crystal. This is a generalisation of
the Fomin growth rules which are an effective version of the operations
on standard tableaux which were previously constructed using jeu-de-taquin.

The basic operations are rectification, evacuation and promotion.
Rectification of standard skew tableaux agrees with the rectification
by jeu-de-taquin as does evacuation. Promotion agrees with promotion
by jeu-de-taquin on rectangular tableaux but in general they are different.

AUTHORS:

- Bruce Westbury (2018): initial version
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

from six import add_metaclass
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.misc.abstract_method import abstract_method
from sage.structure.list_clone import ClonableList
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
#from sage.combinat.partition import Partition, _Partitions

@add_metaclass(InheritComparisonClasscallMetaclass)
class PathTableau(ClonableList):
    @abstract_method(optional=False)
    def _rule_(self,p):
        """
        This is an abstract method. It must be overwritten.
        This rule provides the functionality. It is called in
        :method:`_local_rule`.

        The key property is that the following operation on lists
        of length three is an involution: apply the rule to a list
        and replace the middle term with the output.
        """

################################# Book Keeping ################################

    def size(self):
        """
        Return the size or length of ``self``.

        EXAMPLES::

            sage: t = CatalanTableau([0,1,2,3,2,1,0])
            sage: t.size()
            7
        """
        return len(self)

    def initial_shape(self):
        """
        Return the initial shape of ``self``.

        EXAMPLES::

            sage: t = CatalanTableau([0,1,2,3,2,1,0])
            sage: t.initial_shape()
            0
        """
        return self[0]

    def final_shape(self):
        """
        Return the final shape of ``self``.

        EXAMPLES::

            sage: t = CatalanTableau([0,1,2,3,2,1,0])
            sage: t.final_shape()
            0
        """
        return self[-1]

############################# Jeu de taquin ###################################

    def _local_rule(self,i):
        """
        This is the local that is used for the remaining constructions.
        This has input a list of objects. This method first takes
        the list of objects of length three consisting of the `(i-1)`-st,
        `i`-th and `(i+1)`-term and applies the rule. It then replaces
        the `i`-th object  by the object returned by the rule.

        EXAMPLES::

            sage: t = CatalanTableau([0,1,2,3,2,1,0])
            sage: t._local_rule(3)
            [0, 1, 2, 1, 2, 1, 0]
        """
        if not (i > 0 and i < len(self) ):
            raise ValueError("%d is not a valid integer" % i)

        result = list(self)
        result[i] = self._rule_(self[i-1:i+2])

        return self.parent()(result)

    def promotion(self):
        """
        Return the promotion operator applied to ``self``.

        EXAMPLES::

            sage: t = CatalanTableau([0,1,2,3,2,1,0])
            sage: t.promotion()
            [0, 1, 2, 1, 0, 1, 0]
        """
        result = list(self)
        for i in range(1,len(result)-1):
            result[i] = self._rule_(result[i-1:i+2])
        return self.parent()(result)

    def evacuation(self):
        """
        Return the evacuation operator applied to ``self``.

        EXAMPLES::

            sage: t = CatalanTableau([0,1,2,3,2,1,0])
            sage: t.evacuation()
            [0, 1, 2, 3, 2, 1, 0]
        """
        if self.size() < 3:
            return self

        T = self
        L = list(T)
        result = []
        for i in range(len(self)):
            T = self.parent()(L).promotion()
            L = list(T)
            result.append( L.pop() )
        result.reverse()
        return self.parent()(result)

    def commutor(self,other,display=False):
        """
        Return the commutor of ``self`` with ``other``

        If ``display=True`` then the function will print
        the rectangle.

        EXAMPLES::

            sage: t1 = CatalanTableau([0,1,2,3,2,1,0])
            sage: t2 = CatalanTableau([0,1,2,1,0])
            sage: t1.commutor(t2)
            ([0, 1, 2, 1, 0], [0, 1, 2, 3, 2, 1, 0])
            sage: t1.commutor(t2,display=True)
            [0, 1, 2, 1, 0]
            [1, 2, 3, 2, 1]
            [2, 3, 4, 3, 2]
            [3, 4, 5, 4, 3]
            [2, 3, 4, 3, 2]
            [1, 2, 3, 2, 1]
            [0, 1, 2, 1, 0]
            ([0, 1, 2, 1, 0], [0, 1, 2, 3, 2, 1, 0])
        """
        n = len(self)
        m = len(other)
        if n == 0 or m == 0:
            raise ValueError("this requires nonempty lists")
        if n == 1 or m == 1:
            return (other,self)

        row = list(other)
        col = list(self)
        if col[-1] != row[0]:
            raise ValueError("%s,%s is not a composable pair" % (self,other))

        path = self.parent()(col + row[1:])

        for i in range(1,n):
            if display:
                print path[n-i:n+m-i]
            for j in range(m-1):
                path = path._local_rule(n+j-i)
        if display:
            print path[:m]


        return (self.parent()(path[:m]),self.parent()(path[m-1:]))

    def cactus(self,i,j):
        """
        Return the action of the generators of the cactus group on ``self``.
        These generators are involutions and are usually denoted by
        `s_{i,j}`.

        INPUT:

        - ``i`` -- a positive integer

        - ``j`` -- a positive integer strictly greater than ``i``


        EXAMPLES::

            sage: t = CatalanTableau([0,1,2,3,2,1,0])
            sage: t.cactus(1,5)
            [0, 1, 0, 1, 2, 1, 0]

            sage: t.cactus(1,6)
            [0, 1, 2, 1, 0, 1, 0]

            sage: t.cactus(1,7) == t.evacuation()
            True

            sage: t.cactus(1,7).cactus(1,6) == t.promotion()
            True

        """
        if not 0 < i < j <= self.size():
            raise ValueError("integers out of bounds")

        if i == j:
            return self

        if i == 1:
            h = list(self)[:j]
            t = list(self)[j:]
            T = self.parent()(h)
            L = list(T.evacuation()) + t
            return self.parent()(L)

        return self.cactus(1,j).cactus(1,j-i+1).cactus(1,j)

########################### Visualisation and checking ########################

    def cylindrical_diagram(self):
        """
        Return the cylindrical growth diagram associated to ``self``.

        This provides a visual summary of several operations simultaneously.
        The operations which can be read off directly from this diagram
        are the powers of the promotion operator (which form the rows)
        and the cactus group operators `s_{1,j}` (which form the
        first half of the columns).

        EXAMPLES::

            sage: t = CatalanTableau([0,1,2,3,2,1,0])
            sage: SkewTableau(t.cylindrical_diagram()).pp()
            0  1  2  3  2  1  0
            .  0  1  2  1  0  1  0
            .  .  0  1  0  1  2  1  0
            .  .  .  0  1  2  3  2  1  0
            .  .  .  .  0  1  2  1  0  1  0
            .  .  .  .  .  0  1  0  1  2  1  0
            .  .  .  .  .  .  0  1  2  3  2  1  0
        """
        n = len(self)
        result = [[None]*(2*n-1)]*n
        T = self
        for i in range(n):
            result[i] = [None]*i + list(T)
            T = T.promotion()

        return result

    def _test_involution_rule(self, **options):
        """
        Check that the local rule gives an involution.

        TESTS::

            sage: t = CatalanTableau([0,1,2,3,2,1,0])
            sage: t._test_involution_rule()

        """
        tester = self._tester(**options)
        tester.assertTrue(all( self._local_rule(i+1)._local_rule(i+1) == self
                               for i in range(self.size()-2) ))

    def _test_involution_cactus(self, **options):
        """
        Check that the cactus group generators are involutions.

        TESTS::

            sage: t = CatalanTableau([0,1,2,3,2,1,0])
            sage: t._test_involution_cactus()
        """
        tester = self._tester(**options)
        tester.assertTrue(all( self.cactus(1,i).cactus(1,i) == self
                               for i in range(2,self.size()+1) ))

    def _test_promotion(self, **options):
        """
        Check that promotion can be expressed in terms of the cactus generators.

        TESTS::

            sage: t = CatalanTableau([0,1,2,3,2,1,0])
            sage: t._test_promotion()
        """
        tester = self._tester(**options)
        n = self.size()-1
        tester.assertTrue( self.cactus(1,n-1).cactus(1,n).promotion() == self )

    def _test_commutation(self, **options):
        """
        Check the commutation relations in the presentation of the cactus group.

        TESTS::

            sage: t = CatalanTableau([0,1,2,3,2,1,0])
            sage: t._test_commutation()
        """
        from itertools import combinations
        tester = self._tester(**options)

        n = self.size()
        if n < 5:
            return True
        for i,j,r,s in combinations(range(1,n+1),4):
            lhs = self.cactus(i,j).cactus(r,s)
            rhs = self.cactus(r,s).cactus(i,j)
            tester.assertTrue(lhs == rhs)

    def _test_coboundary(self, **options):
        """
        Check the coboundary relations in the presentation of the cactus group.

        TESTS::

            sage: t = CatalanTableau([0,1,2,3,2,1,0])
            sage: t._test_coboundary()
        """
        from itertools import combinations
        tester = self._tester(**options)

        n = self.size()
        if n < 4:
            return True
        for i,j,r,s in combinations(range(1,n+3),4):
            lhs = self.cactus(i,s-2).cactus(j-1,r-1)
            rhs = self.cactus(i+s-r-1,i+s-j-1).cactus(i,s-2)
            tester.assertTrue(lhs == rhs)

    def orbit(self):
        """
        Return the orbit of ``self`` under the action of the cactus group.

        EXAMPLES::

            sage: t = CatalanTableau([0,1,2,3,2,1,0])
            sage: t.orbit()
            {[0, 1, 0, 1, 0, 1, 0],
             [0, 1, 0, 1, 2, 1, 0],
             [0, 1, 2, 1, 0, 1, 0],
             [0, 1, 2, 1, 2, 1, 0],
             [0, 1, 2, 3, 2, 1, 0]}
        """

        orb = set([])
        rec = set([self])
        new = set([])
        while rec != set([]):
            for a in rec:
                for i in range(2,self.size()):
                    b = a.cactus(1,i)
                    if (b not in orb) and (b not in rec):
                        new.add(b)
            orb = orb.union(rec)
            rec = new.copy()
            new = set([])

        return orb

    def dual_equivalence_graph(self):
        """
        Return the graph with vertices the orbit of ``self``
        and edges given by the action of the cactus group generators.

        In most implementations the generators `s_{i,i+1}` will act
        as the identity operators. The usual dual equivalence graphs
        are given by replacing the label `i,i+2` by `i` and removing
        edges with other labels.

        EXAMPLES::

            sage: s = CatalanTableau([0,1,2,3,2,3,2,1,0])
            sage: s.dual_equivalence_graph().adjacency_matrix()
            [0 1 1 1 0 1 0 1 1 0 0 0 0 0]
            [1 0 1 1 1 1 1 0 1 0 0 1 1 0]
            [1 1 0 1 1 1 0 1 0 1 1 1 0 0]
            [1 1 1 0 1 0 1 1 1 1 0 1 1 0]
            [0 1 1 1 0 0 1 0 0 1 1 0 1 1]
            [1 1 1 0 0 0 1 1 1 1 1 0 1 0]
            [0 1 0 1 1 1 0 1 0 1 1 1 0 1]
            [1 0 1 1 0 1 1 0 1 1 1 1 1 0]
            [1 1 0 1 0 1 0 1 0 1 0 1 1 0]
            [0 0 1 1 1 1 1 1 1 0 0 1 1 1]
            [0 0 1 0 1 1 1 1 0 0 0 1 1 1]
            [0 1 1 1 0 0 1 1 1 1 1 0 1 1]
            [0 1 0 1 1 1 0 1 1 1 1 1 0 1]
            [0 0 0 0 1 0 1 0 0 1 1 1 1 0]

        """
        from sage.graphs.graph import Graph
        from itertools import combinations

        G = Graph()
        orb = self.orbit()

        for a in orb:
            for i,j in combinations(range(1,self.size()+1),2):
                b = a.cactus(i,j)
                if a != b:
                    G.add_edge(a,b,"%d,%d" % (i,j))
        return G

class PathTableaux(UniqueRepresentation,Parent):
#
#    def __init__(self):
#        Parent.__init__(self, category = Sets())
#
    def _element_constructor_(self, *args, **keywords):
        return self.element_class(self, *args, **keywords)
#
#    Element = PathTableau
