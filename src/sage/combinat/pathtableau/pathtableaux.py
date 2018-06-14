r"""
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
from sage.combinat.partition import Partition
#from sage.modules.free_module_element import vector


@add_metaclass(InheritComparisonClasscallMetaclass)
class PathTableau(ClonableList):

    @staticmethod
    def __classcall_private__(cls, t):

        if isinstance(t, cls):
            return t

        raise NotImplementedError("This needs to be overwritten.")

    @abstract_method(optional=False)
    def check(self):
        """
        This is an abstract method. It must be overwritten
        Typically an instance of
        an Element class is a sequence of partitions with conditions
        on adjacent partitions in the sequence. This function checks
        that these conditions are met.
        """

    @abstract_method(optional=False)
    def _rule(self,p):
        """
        This is an abstract method. It must be overwritten.
        This rule provides the functionality. It is called in local_rule.

        The key property is that the following operation on lists
        of length three is an involution: apply the rule to a list
        and replace the middle term with the output.
        """

################################# Book Keeping ################################

    def size(self):
        """
        Returns the size or length.
        """
        return len(self)

    def initial_shape(self):
        """
        Returns the initial shape.
        """
        return self[0]

    def final_shape(self):
        """
        Returns the final shape.
        """
        return self[-1]

############################# Jeu de taquin ###################################

    def local_rule(self,i):
        """
        This is the local that is used for the remaining constructions.
        This has input a list of objects. This method first takes
        the list of objects of length three consisting of the $(i-1)$-st,
        $i$-th and $(i+1)$-term and applies the rule. It then replaces
        the $i$-th object  by the object returned by the rule.
        """
        if not (i > 0 and i < len(self) ):
            raise ValueError("%d is not a valid integer." % i)

        result = list(self)
        result[i] = self._rule(self[i-1:i+2])

        return self.parent()(result)

    def promotion(self):
        """
        The promotion operator. This is given by a two row diagram.
        """
        result = list(self)
        for i in range(1,len(result)-1):
            result[i] = self._rule(result[i-1:i+2])
        return self.parent()(result)

    def evacuation(self):
        """
        The evacuation operator. This is given by a triangular diagram.
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

    def path_rule(self,other,display=False):
        """
        This constructs the commutor of a pair of tableau.
        This is given by a rectangular diagram.

        If display=True then the function will print
        the rectangle.
        """
        n = len(self)
        m = len(other)
        if n == 0 or m == 0:
            raise ValueError("This requires nonempty lists.")
        if n == 1 or m == 1:
            return (other,self)

        row = list(other)
        col = list(self)
        if col[-1] != row[0]:
            raise ValueError("%s,%s is not a composable pair." % (self,other))

        path = self.parent()(col + row[1:])

        for i in range(1,n):
            if display:
                print path[n-i:n+m-i]
            for j in range(m-1):
                path = path.local_rule(n+j-i)
        if display:
            print path[:m]


        return (self.parent()(path[:m]),self.parent()(path[m-1:]))

    def cactus(self,i,j):
        """
        This constructs the action of the generators of the cactus group.
        These generators are involutions and are usually denoted by
        $s_{i,\,j$}$.
        """
        if not 0 < i < j <= self.size():
            raise ValueError("Integers out of bounds.")

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
        This constructs the cylindrical growth diagram. This provides
        a visual summary of several operations simultaneously. The
        operations which can be read off directly from this diagram
        are the powers of the promotion operator (which form the rows)
        and the cactus group operators $s_{1,\,j$}$ (which form the
        first half of the columns).
        """
        n = len(self)
        result = [[None]*(2*n-1)]*n
        T = self
        for i in range(n):
            result[i] = [None]*i + list(T)
            T = T.promotion()

        return result

    def check_involution_rule(self):
        """
        This is to check that the local rule gives an involution.
        This is crucial.
        """
        for i in range(self.size()-2):
            if self.local_rule(i+1).local_rule(i+1) != self:
                return False
        return True

    def check_involution_cactus(self):
        """
        This is to check that the cactus group generators are
        involutions..
        """
        return all([ self.cactus(1,i).cactus(1,i) == self for i in range(2,self.size()+1 ) ])

    def check_promotion(self):
        """
        Promotion can be expressed in terms of the cactus generators.
        Here we check this relation.
        """
        n = self.size()-1
        return self == self.cactus(1,n-1).cactus(1,n).promotion()

    def check_commutation(self):
        """
        This is to check the commutation relations in the presentation
        of the cactus group.
        """
        from itertools import combinations

        n = self.size()
        if n < 5:
            return True
        for i,j,r,s in combinations(range(1,n+1),4):
            lhs = self.cactus(i,j).cactus(r,s)
            rhs = self.cactus(r,s).cactus(i,j)
            if lhs != rhs:
                return False
        return True

    def check_coboundary(self):
        """
        This is to check the coboundary relations in the presentation
        of the cactus group.
        """
        from itertools import combinations

        n = self.size()
        if n < 4:
            return True
        for i,j,r,s in combinations(range(1,n+3),4):
            lhs = self.cactus(i,s-2).cactus(j-1,r-1)
            rhs = self.cactus(i+s-r-1,i+s-j-1).cactus(i,s-2)
            if lhs != rhs:
                return False
        return True

    def orbit(self):
        """
        Constructs the orbit under the action of the cactus group.
        """
        n = self.size()
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
        This constructs the graph with vertices the orbit of self
        and edges given by the action of the cactus group generators.

        In most implementations the generators $s_{i,\,i+1}$ will act
        as the identity operators. The usual dual equivalence graphs
        are given by replacing the label $i,i+2$ by $i$ and removing
        edges with other labels.

        PLOT::

            sage: t = SkewTableau([[None,1,1],[2,2]])
            sage: s = DualSemistandardTableau(t)
            sage: s.dual_equivalence_graph().show()
            Launched png viewer for Graphics object consisting of 4 graphics primitives

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

    def csp(self):
        import sage.combinat.cyclic_sieving_phenomenon

#### These functions don't belong here but I don't have a home for them. ####



class PathTableaux(UniqueRepresentation,Parent):
#
#    def __init__(self):
#        Parent.__init__(self, category = Sets())
#
    def _element_constructor_(self, *args, **keywords):
        return self.element_class(self, *args, **keywords)
#
#    Element = PathTableau

class PathTableau_partitions(PathTableau):
    """
    This is an abstract base class. This class assumes that we have
    a sequence of partitions. The main examples are the minuscule
    representations of classical groups.
    """

    @staticmethod
    def _rule(x):
        y = map(list,x)
        m = max([ len(u) for u in y ])
        z = map( lambda u: u + [0]*(m-len(u)), y )
        result = [ abs(a-b+c) for a,b,c in zip(z[0],z[1],z[2]) ]
        result.sort(reverse=True)
        return Partition(result)

    def _plotL(self):
        """
        This draws a plot of the sequence of partitions.
        This plot assumes we do not have a chain of partitions
        and plots the partitions in a line.

        PLOT::

            sage: t = SkewTableau([[None,1,1],[2,2]])
            sage: s = DualSemistandardTableau(t)
            sage: s._plotL()
            Launched png viewer for Graphics object consisting of 11 graphics primitives

        """
        from sage.plot.graphics import Graphics
        from sage.plot.line import line
        from copy import copy

        global gap
        gap = 1

        def draw_partition(p,origin):

            global gap

            if p == Partition([]):
                return point(origin,axes=False,size=60)

            r = origin[0]
            s = origin[1]

            u = p.to_dyck_word()
            u = u[u.index(0):]
            u.reverse()
            u = u[u.index(1):]
            u.reverse()
            x = u.count(0)
            y = u.count(1)

            gap = max(x,gap)
            n = len(u)

            edge = []
            edge.append([r,-y+s])
            for i in range(n):
                v = copy(edge[i])
                if u[i] == 1:
                    v[1] += 1
                else:
                    v[0] += 1
                edge.append(v)

            G = Graphics()
            G += line([(r,-y+s),(r,s),(r+x,s)],axes=False,thickness=2)
            G += line(edge,color='red',axes=False,thickness=3)

            for i, a in enumerate(p[1:]):
                G += line([(r,s-i-1),(r+a,s-i-1)],color='green')

            for i, a in enumerate(p.conjugate()[1:]):
                G += line([(r+i+1,s),(r+i+1,s-a)],color='green')

            return G

        G = Graphics()

        for i, x in enumerate(self):
            G += draw_partition(x, (i*gap+1.5*i,0))

        G.set_aspect_ratio(1)

        return G

    def _plotC(self):
        """
        This draws a plot of the sequence of partitions.
        This plot assumes the sequence is not a chain and so
        plots the sequence.

        PLOT::

            sage: t = SkewTableau([[None,1,1],[2,2]])
            sage: s = DualSemistandardTableau(t)
            sage: s._plotC()
            Launched png viewer for Graphics object consisting of 10 graphics primitives

        """
        from sage.plot.graphics import Graphics
        from sage.plot.line import line
        from copy import copy

        def draw_partition(p):

            if p == Partition([]):
                return point((0,0),axes=False,size=60)

            u = p.to_dyck_word()
            u = u[u.index(0):]
            u.reverse()
            u = u[u.index(1):]
            u.reverse()
            x = u.count(0)
            y = u.count(1)

            n = len(u)

            edge = []
            edge.append([0,-y])
            for i in range(n):
                v = copy(edge[i])
                if u[i] == 1:
                    v[1] += 1
                else:
                    v[0] += 1
                edge.append(v)

            return line(edge,color='red',axes=False,thickness=3)

        p = self.final_shape()

        G = line([(0,-len(p)),(0,0),(p[0],0)],axes=False)

        for i, a in enumerate(p[1:]):
            G += line([(0,-i-1),(a,-i-1)],color='green')

        for i, a in enumerate(p.conjugate()[1:]):
            G += line([(i+1,0),(i+1,-a)],color='green')

        for i, x in enumerate(self):
            G += draw_partition(x)

        for p in self:
            G += draw_partition(p)

        G.set_aspect_ratio(1)

        return G