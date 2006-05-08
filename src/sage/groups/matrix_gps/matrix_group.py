"""
Matrix Groups

AUTHORS:
   William Stein -- initial version
   David Joyner  -- degree, base_ring, _contains_, list, random, order
                    methods; examples (2006-03-15)

Generally, the instances of the class defined in this file is designed
to represent matrix groups given by a (small) finite set of generating
matrices.
"""

##############################################################################
#       Copyright (C) 2006 David Joyner and William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
##############################################################################

from matrix_group_element import MatrixGroupElement
from sage.groups.group import Group
from sage.rings.all import IntegerRing, is_Ring
from sage.misc.functional import is_field
from sage.rings.finite_field import is_FiniteField
from sage.interfaces.gap import gap, GapElement
from sage.matrix.matrix_space import MatrixSpace, is_MatrixSpace
import sage.rings.integer as integer
from sage.misc.latex import latex
from sage.structure.sequence import Sequence


#################################################################

class MatrixGroup(Group):
    """
    EXAMPLES:

    A ValueError is raised if one of the generators is not invertible.

        sage: F = GF(5); MS = MatrixSpace(F,2,2)
        sage: G = MatrixGroup([MS.0])
        Traceback (most recent call last):
        ...
        ValueError: each generator must be an invertible matrix but one is not:
        [1 0]
        [0 0]
    """
    def __init__(self, gensG):
        # TODO: gensG should be private
        v = Sequence(gensG, immutable=True)
        M = v.universe()
        if not is_MatrixSpace(M):
            raise TypeError, "universe of sequence (=%s) of generators must be a matrix space"%M
        if M.nrows() != M.ncols():
            raise ValueError, "matrices must be square."
        for x in v:
            if not x.is_invertible():
                raise ValueError, "each generator must be an invertible matrix but one is not:\n%s"%x
        self.__ngens = len(v)

        #gensG is a list of matrices generating the group
        self.__n = M.nrows()
        self.__R = M.base_ring()
        self.__matrix_space = M

        self.gensG = v

    def __call__(self, x):
        """
        EXAMPLES:
            sage: F = GF(5); MS = MatrixSpace(F,2,2)
            sage: G = MatrixGroup([MS(1), MS([1,2,3,4])])
            sage: G.matrix_space()
            Full MatrixSpace of 2 by 2 dense matrices over Finite Field of size 5
            sage: G(1)
            [1 0]
            [0 1]
        """
        if isinstance(x, MatrixGroupElement) and x.parent() is self:
            return x
        return MatrixGroupElement(self.__matrix_space(x), self)

    def hom(self, x):
        v = Sequence(x)
        U = v.universe()
        if not isinstance(U, MatrixGroup):
            raise TypeError, "u (=%s) must have universe a matrix group."%U
        return self.Hom(U)(x)

    def _Hom_(self, G, cat=None):
        if not (cat is None or (cat is G.category() and cat is self.category())):
            raise NotImplementedError
        if not isinstance(G, MatrixGroup):
            raise TypeError, "G (=%s) must be a matrix group."%G
        import homset
        return homset.MatrixGroupHomset(self, G)

    def matrix_space(self):
        """
        EXAMPLES:
            sage: F = GF(5); MS = MatrixSpace(F,2,2)
            sage: G = MatrixGroup([MS(1), MS([1,2,3,4])])
            sage: G.matrix_space()
            Full MatrixSpace of 2 by 2 dense matrices over Finite Field of size 5
        """
        return self.__matrix_space

    def degree(self):
        return self.__n

    def base_ring(self):
        return self.__R

    def ngens(self):
        return self.__ngens

    def gens(self):
        """
        EXAMPLES:
            sage: F = GF(3); MS = MatrixSpace(F,2,2)
            sage: gens = [MS([[1,0],[0,1]]),MS([[1,1],[0,1]])]
	    sage: G = MatrixGroup(gens)
            sage: gens[0] in G
            False
            sage: gens = G.gens()
            sage: gens[0] in G
            True
            sage: gens = [MS([[1,0],[0,1]]),MS([[1,1],[0,1]])]


            sage: F = GF(5); MS = MatrixSpace(F,2,2)
            sage: G = MatrixGroup([MS(1), MS([1,2,3,4])])
            sage: G
            Matrix group over Finite Field of size 5 with 2 generators
            sage: G.gens()
            [[1 0]
            [0 1], [1 2]
            [3 4]]
        """
        try:
            return self.__gens
        except AttributeError:
            t = Sequence([MatrixGroupElement(x, self) for x in self.gensG],
                         immutable=True, universe=self)
            self.__gens = t
            return t

    def gen(self, n):
        return self.gens()[n]

    def __contains__(self, x):
        """
        Return True if $x$ is an element of this abelian group.

        EXAMPLES:
            sage: F = GF(5); MS = MatrixSpace(F,2,2)
            sage: gens = [MS([[1,2],[-1,1]]),MS([[1,1],[0,1]])]
            sage: G = MatrixGroup(gens)
            sage: g = G.random()
            sage: g.parent()
            Matrix group over Finite Field of size 5 with 2 generators
            sage: g in G
            True

        """
        from matrix_group_element import MatrixGroupElement
        return isinstance(x, MatrixGroupElement) and x.parent() == self

    def _gap_init_(self):
        """
        Returns the string representation of the corresponding GAP command.

        EXAMPLES:
            sage: F = GF(5); MS = MatrixSpace(F,2,2)
            sage: gens = [MS([[1,2],[-1,1]]),MS([[1,1],[0,1]])]
            sage: G = MatrixGroup(gens)
            sage: G._gap_init_()
            'Group([ [ [ Z(5)^0, Z(5) ], [ Z(5)^2, Z(5)^0 ] ],
                     [ [ Z(5)^0, Z(5)^0 ], [ 0*Z(5), Z(5)^0 ] ] ])'
        """
        gens_gap = [gap(x) for x in self.gensG]
        return gap.eval("Group("+str(gens_gap)+")").replace("\n","")

    def _repr_(self):
        return "Matrix group over %s with %s generators"%(self.base_ring(), self.ngens())

    def _latex_(self):
        gens = ', '.join([latex(x) for x in self.gens()])
        return '\\left\\langle %s \\right\\rangle'%gens

    def is_finite(self):
        """
        EXAMPLES:
            sage: G = GL(2,GF(3))
            sage: G.is_finite()
            True
        """
        return self.base_ring().is_finite()

    def random(self):
        """
        Returns a random matrix from the matrix group.

        EXAMPLES:
            sage: F = GF(5); MS = MatrixSpace(F,2,2)
	    sage: gens = [MS([[1,2],[-1,1]]),MS([[1,1],[0,1]])]
	    sage: G = MatrixGroup(gens)
	    sage: G.random()
            [2 4]
            [0 1]
            sage: G.random()
            [3 1]
            [2 2]
            sage: G.random()
            [3 2]
            [4 2]

        """
        from sage.interfaces.gap import GapElement
        if not(self.is_finite()):
            return Infinity()
        F = self.base_ring()
        deg = self.degree()
        gp_gens = self.gens()
        L = [gap(A) for A in gp_gens]
        sL = ','.join(str(x) for x in L)
        ans = gap("Random(Group(["+sL+"]))")
        return MatrixGroupElement(ans._matrix_(F),self)

    def order(self):
        """
        EXAMPLES:
            sage: F = GF(5); MS = MatrixSpace(F,2,2)
            sage: gens = [MS([[1,2],[-1,1]]),MS([[1,1],[0,1]])]
            sage: G = MatrixGroup(gens)
            sage: G.order()
            480
        """
        #from sage.misc.misc import add
        from sage.interfaces.gap import GapElement
        if not(self.is_finite()):
            return Infinity()
        F = self.base_ring()
        deg = self.degree()
        gp_gens = self.gens()
        L = [gap(A) for A in gp_gens]
        sL = ','.join(str(x) for x in L)
        ans = gap.eval('Size(Group([%s]))'%sL)
        return eval(ans)

    def conjugacy_class_representatives(self):
        """
        Wraps GAP Representative+ConjugactClasses.

        EXAMPLES:
            sage: F = GF(3); MS = MatrixSpace(F,2,2)
            sage: gens = [MS([[1,0],[-1,1]]),MS([[1,1],[0,1]])]
            sage: G = MatrixGroup(gens)
            sage: G.conjugacy_class_representatives()
            [[1 0]
            [0 1], [0 1]
            [2 1], [0 1]
            [2 2], [0 2]
            [1 1], [0 2]
            [1 2], [0 1]
            [2 0], [2 0]
            [0 2]]

        AUTHOR: David Joyner (1-2006)
        """
        #from sage.misc.misc import add
        #from sage.interfaces.gap import GapElement
        F = self.base_ring()
        deg = self.degree()
        gp_gens = self.gens()
        L = [gap(A) for A in gp_gens]
        sL = '[' + ','.join(str(x) for x in L) + ']'
        #print sL
        if is_FiniteField(F):
            q = F.order()
            gap.eval("cl:=ConjugacyClasses(Group("+sL+"))")
            m = eval(gap.eval("Length(cl)"))
            #print m, L
            gap.eval("reps:=List(cl,x->Representative(x))")
            sreps = [gap("reps["+str(i+1)+"]") for i in range(m)]
            reps = [MatrixGroupElement(s._matrix_(F),self) for s in sreps]
            return reps
        raise TypeError, "R (=%s) must be a finite field"%R

    def conjugacy_class_representatives_gap(self):
        """
        Wraps GAP Representative+ConjugactClasses but returns a list of
        strings representing the GAP matrices which form a complete
        set of representatives of the conjugacy classes of the group.

        EXAMPLES:
            sage: F = GF(3); MS = MatrixSpace(F,2,2)
            sage: gens = [MS([[1,0],[-1,1]]),MS([[1,1],[0,1]])]
            sage: G = MatrixGroup(gens)
            sage: G.conjugacy_class_representatives_gap()
	    ['[ [ Z(3)^0, 0*Z(3) ], [ 0*Z(3), Z(3)^0 ] ]',
 	    '[ [ 0*Z(3), Z(3)^0 ], [ Z(3), Z(3)^0 ] ]',
            '[ [ 0*Z(3), Z(3)^0 ], [ Z(3), Z(3) ] ]',
            '[ [ 0*Z(3), Z(3) ], [ Z(3)^0, Z(3)^0 ] ]',
     	    '[ [ 0*Z(3), Z(3) ], [ Z(3)^0, Z(3) ] ]',
     	    '[ [ 0*Z(3), Z(3)^0 ], [ Z(3), 0*Z(3) ] ]',
 	    '[ [ Z(3), 0*Z(3) ], [ 0*Z(3), Z(3) ] ]']

        AUTHOR: David Joyner (1-2006)
        """
        F = self.__R
        deg = self.__n
        gp_gens = self.gens()
        L = [gap(A) for A in gp_gens]
        sL = ','.join(str(x) for x in L)
        if is_FiniteField(F):
            q = F.order()
            gap.eval("cl:=ConjugacyClasses(Group(["+sL+"]))")
            m = eval(gap.eval("Length(cl)"))
            gap.eval("reps:=List(cl,x->Representative(x))")
            sreps = [gap.eval("reps["+str(i+1)+"]") for i in range(m)]
            return sreps
        raise TypeError, "R (=%s) must be a finite field"%R


    def list(self):
        """
        Return list of all elements of this group.

        EXAMPLES:
            sage: F = GF(3); MS = MatrixSpace(F,2,2)
            sage: gens = [MS([[1,0],[0,1]]),MS([[1,1],[0,1]])]
	    sage: G = MatrixGroup(gens)
	    sage: G.order()
	    3
	    sage: G.list()
            [[1 0]
            [0 1], [1 1]
            [0 1], [1 2]
            [0 1]]
            sage: G.list()[0] in G
	    True
        """
        X = self._gap_().Elements()
        n = X.Length()
        F = self.base_ring()
        return [MatrixGroupElement(X[i]._matrix_(F), self)
                            for i in range(1,n+1)]

