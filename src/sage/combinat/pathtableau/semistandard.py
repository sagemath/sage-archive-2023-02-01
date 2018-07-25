r"""
In this implementation we have sequences of partitions. These are in
bijection with dual semistandard tableaux. This gives an effective
version of operations on tableaux constructed using jeu-de-taquin.
In the standard constructions of these operations one usually assumes
the tableau is standard.

For rectification and evacuation the operations here
agree with the standard construction. For promotion the construction
here agrees with the standard construction on rectangular standard
tableaux, but, in general, they are different.

The operations here also give the Bender-Knuth involutions and
dual equivalence graphs.

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
from sage.combinat.pathtableau.pathtableaux import PathTableau_partitions, PathTableaux
from sage.combinat.skew_tableau import SkewTableau
from sage.combinat.tableau import SemistandardTableau
from sage.combinat.partition import Partition

"""
This implementation is on dual semistandard tableaux. This is the standard
context for jeu-de-taquin operations. Here we show that the constructions
here agree with the jeu-de-taquin constructions. Even in this standard context
our operations extend the standard definitions in the sense that the
constructions here are naturally defined for skew semistandard tableaux.
The only caveat is that the two constructions of promotion agree on rectangular
tableaux but are, in general, different.

-promotion
-evacuation
-rectify
-multiply
-Bender-Knuth

EXAMPLES::

    sage: T = DualSemistandardTableau([[],[1],[2],[2,1]])
    sage: T.evacuation()
    [[], [1], [1, 1], [2, 1]]

    sage: Tableau([[1,2],[3]]).evacuation()
    [[1, 3], [2]]

    sage: ST = SemistandardTableaux([5,3,3,2,1],[2,1,4,2,2,2,1])
    sage: ST.cardinality()
    84

    sage: t = ST.an_element()
    sage: s = DualSemistandardTableau(t.conjugate().to_chain())
    sage: v = Tableau(list(SkewTableau(chain=s.evacuation())))
    sage: v.conjugate() == t.evacuation()
    True

    sage: ST = SemistandardTableaux([5,3,3,2,1],[2,1,4,2,2,2,1])
    sage: s = DualSemistandardTableau(ST.an_element())
    sage: s.check_involution_cactus()
    True

    sage: s.check_commutation()
    True

    sage: s.check_coboundary()
    True

    sage: ST = StandardTableaux([3,3,3])
    sage: ST.cardinality()
    42

    sage: t = ST.an_element()
    sage: t.promotion()
    [[1, 2, 5], [3, 6, 8], [4, 7, 9]]

    sage: ST = StandardTableaux([3,3,3])
    sage: t = ST.an_element()
    sage: s = DualSemistandardTableau(t.to_chain())
    sage: u = StandardTableau(list(SkewTableau(chain=s.promotion())))
    sage: u.promotion() == t
    True

    sage: ST.cardinality()
    42
    sage: t = ST.an_element()
    sage: s = DualSemistandardTableau(t.to_chain())
    sage: len(s.orbit())
    42

"""

@add_metaclass(InheritComparisonClasscallMetaclass)
class DualSemistandardTableau(PathTableau_partitions):
    """
       An instance is the sequence of partitions correspond to the
       chain of partitions of a dual semistandard skew tableau.

    The acceptable inputs are:
        - a sequence such that each term defines a partition
        - a semistandard skew tableau

    EXAMPLES::

        sage: DualSemistandardTableau([[],[1],[2],[2,1]])
        [[], [1], [2], [2, 1]]

        sage: t = SkewTableau([[None,None,None,4,4,5,6,7],[None,2,4,6,7,7,7],[None,4,5,8,8,9],[None,6,7,10],[None,8,8,11],[None],[4]])
        sage: DualSemistandardTableau(t)
        [[6, 1, 1], [6, 1, 1], [6, 2, 1], [6, 2, 1], [7, 3, 2, 1, 1], [7, 3, 3, 1, 1, 1], [7, 4, 3, 2, 1, 1, 1], [7, 4, 4, 2, 2, 2, 2, 1], [7, 5, 5, 3, 3, 2, 2, 1], [7, 5, 5, 3, 3, 3, 2, 1], [7, 5, 5, 4, 3, 3, 2, 1], [7, 5, 5, 5, 3, 3, 2, 1]]

    """
    @staticmethod
    def __classcall_private__(self, ot):

        w = None

        if isinstance(ot,(SkewTableau,SemistandardTableau)):
            w = ot.conjugate().to_chain()

        if isinstance(ot,(list,tuple)):
            try:
                w = tuple([ Partition(a) for a in ot ])
            except TypeError:
                raise ValueError("%s is not a sequence of partitions." % str(ot) )

        if w == None:
            raise ValueError( "Sorry, not sorry; I don't know what to do with %s." % str(ot) )

        return DualSemistandardTableaux()(w)

    def _hash_(self):
        return hash(tuple(map(tuple, self)))

    def check(self):
        n = len(self)
        for i in range(n-1):
            h = self[i]
            t = self[i+1]
            if not t.contains(h):
                raise ValueError( "%s must contain %s" % (str(t),str(h)) )
            for r, s in zip(h,t):
                if s > r+1:
                    raise ValueError( "%s / %s is not a vertical strip" % (str(t),str(h)) )
            for a in t[len(h):]:
                if a > 1:
                    raise ValueError( "%s / %s is not a vertical strip" % (str(t),str(h)) )

    def evaluation(self):
        z = [ p.size() for p in self ]
        return [ z[i+1] - z[i] for i in range(len(self)-1) ]

    def to_tableau(self):
        """
        Returns the conjugate skew tableau. This will be semistandard.
        """
        ch = [ p.conjugate() for p in self]
        s = SkewTableau(chain=ch)
        if self.is_skew():
            return s
        else:
            return SemistandardTableau(list(s))

    def is_skew(self):
        """
        Returns True if Tableau is skew and False if not.

        EXAMPLE::
            sage: t = DualSemistandardTableau([[],[1],[2],[2,1]])
            sage: t.is_skew()
            False
        """
        return self[0] != Partition([])

    def rectify(self,display=False):
        """
        This is the same function as skew_tableau.rectify

        EXAMPLE::

            sage: t = SkewTableau([[None,None,2,4],[1,3,5]])
            sage: s = DualSemistandardTableau(t.to_chain())
            sage: s.rectify(display=True)
            [[2], [2, 1], [3, 1], [3, 2], [4, 2], [4, 3]]
            [[1], [1, 1], [2, 1], [2, 2], [3, 2], [3, 3]]
            [[], [1], [2], [2, 1], [3, 1], [3, 2]]
            [[], [1], [2], [2, 1], [3, 1], [3, 2]]

        """
        p = self[0].conjugate()
        path = [[]]
        for i in range(len(p)):
            path += [Partition(p[:i+1]).conjugate()]

        return DualSemistandardTableau(path).path_rule(self,display=display)[0]

    def multiply(self,other):
        """
        This is the same function as tableau.slide_multiply and tableau.bump_multiply.

        EXAMPLE::

            sage: t = DualSemistandardTableau([[2],[3,1],[4,1,1]])
            sage: t.multiply(t)
            [[], [1, 1, 1, 1], [2, 2, 2, 1, 1]]

        """

        left = list(self)
        right = list(other)

        m = max([len(a) for a in right])
        n = max([ a[0] for a in left])

        right = [a+[0]*(m-len(a)) for a in right]

        p = max(len(left),len(right))
        left = left + [left[-1]]*(p-len(left))
        right = right + [right[-1]]*(p-len(right))

        result = [Partition([a+n for a in y]+list(x)) for x,y in zip(left,right)]

        return DualSemistandardTableau(result).rectify()

    def check_bender_knuth(self,i):
        """
        Check that the i-th Bender-Knuth move on the conjugate
        tableau is the i-th local rule.

        EXAMPLE::

            sage: t = SemistandardTableaux(8).random_element()
            sage: s = DualSemistandardTableau(t)
            sage: s.check_bender_knuth(5)
            True
            sage: s.check_bender_knuth(4)
            True

        """

        lhs = self.local_rule(i).to_tableau()
        rhs = self.to_tableau().bender_knuth_involution(i)
        return lhs == rhs

    def check_rectify(self):
        """
        Check that jdt-rectification on the conjugate tableaux
        is the rectification defined here.

        EXAMPLE::

            sage: t = SkewTableau([[None,None,1,3,3],[None,1,2,4],[2,2]])
            sage: s = DualSemistandardTableau(t)
            sage: s.check_rectify()
            True

        """

        lhs = self.rectify().to_tableau()
        rhs = self.to_tableau().rectify()
        return lhs == rhs

    def check_evacuation(self):
        """
        Check that jdt-evacuation on the conjugate tableaux
        is the evacuation defined here.

        EXAMPLE::

            sage: t = SemistandardTableaux(6).random_element()
            sage: s = DualSemistandardTableau(t)
            sage: s.check_evacuation()
            True

        """
        lhs = self.evacuation().to_tableau()
        rhs = self.to_tableau().evacuation()
        return lhs == rhs

    def check_promotion(self):
        """
        Check that jdt-promotion on the conjugate tableaux
        is the promotion defined here.

        EXAMPLE::

            sage: t = SemistandardTableaux(shape=[4,4,4],eval=[1]*12).an_element()
            sage: s = DualSemistandardTableau(t)
            sage: s.check_promotion()
            True

        """
        lhs = self.promotion().to_tableau()
        rhs = self.to_tableau().promotion_inverse(11)
        return lhs == rhs


    def plot(self):
        """
        This provides a plot of the dual semistandard tableau.
        PLOT::

            sage: t = SkewTableau([[None,None,2,2],[3,4,4],[4]])
            sage: DualSemistandardTableau(t).plot()
            Graphics object consisting of 16 graphics primitives

        """
        return self._plotC()

class DualSemistandardTableaux(PathTableaux):

    Element = DualSemistandardTableau

#### These functions don't belong here but I don't have a home for them. ####

#############################################################################



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

    TESTS::

        sage: F = Foo()
        sage: TestSuite(F).run()
    """

    @staticmethod
    def _rule_(x):
        y = map(list, x)
        m = max(len(u) for u in y) # FIXME?: This will fail if y is empty
        z = [u + [0]*(m-len(u)) for u in y]
        result = [ abs(a-b+c) for a,b,c in zip(z[0],z[1],z[2]) ]
        result.sort(reverse=True)
        return _Partitions(result)

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

            if p == _Partitions([]):
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

    def _plotC_(self):
        """
        This draws a plot of the sequence of partitions.
        This plot assumes the sequence is not a chain and so
        plots the sequence.

        PLOT::

            sage: t = SkewTableau([[None,1,1],[2,2]])
            sage: s = DualSemistandardTableau(t)
            sage: s._plotC_()
            Launched png viewer for Graphics object consisting of 10 graphics primitives

        """
        from sage.plot.graphics import Graphics
        from sage.plot.line import line
        from copy import copy

        def draw_partition(p):

            if p == _Partitions([]):
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