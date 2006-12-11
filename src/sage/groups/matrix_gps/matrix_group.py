"""
Matrix Groups

AUTHORS:
   William Stein -- initial version
   David Joyner  -- degree, base_ring, _contains_, list, random, order
                    methods; examples (2006-03-15)
   William Stein (2006-12) -- rewrite

This class is designed for computing with matrix groups defined by a
relatively (small) finite set of generating matrices.

EXAMPLES:
    sage: F = GF(3)
    sage: gens = [matrix(F,2, [1,0, -1,1]), matrix(F, 2, [1,1,0,1])]
    sage: G = MatrixGroup(gens)
    sage: G.conjugacy_class_representatives()
    [
    [1 0]
    [0 1],
    [2 0]
    [0 2]
    ]

Loading and saving work:
    sage: G = GL(2,5); G
    General Linear Group of degree 2 over Finite Field of size 5
    sage: loads(dumps(G)) == G
    True
    sage: g = G.1; g
    [4 1]
    [4 0]
    sage: loads(dumps(g)) == g
    True
"""

##############################################################################
#       Copyright (C) 2006 David Joyner and William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
##############################################################################

from matrix_group_element import MatrixGroupElement
from sage.groups.group import Group
from sage.rings.all import IntegerRing, is_Ring, infinity
from sage.misc.functional import is_field
from sage.rings.finite_field import is_FiniteField
from sage.interfaces.gap import gap, GapElement
from sage.matrix.all import MatrixSpace, is_MatrixSpace, is_Matrix
import sage.rings.integer as integer
from sage.misc.latex import latex
from sage.structure.sequence import Sequence


#################################################################


class MatrixGroup_generic(Group):
    pass

def is_MatrixGroup(x):
    """
    EXAMPLES:
        sage: is_MatrixGroup(MatrixSpace(QQ,3))
        False
        sage: is_MatrixGroup(Mat(QQ,3))
        False
        sage: is_MatrixGroup(GL(2,ZZ))
        True
        sage: is_MatrixGroup(MatrixGroup([matrix(2,[1,1,0,1])]))
        True
    """
    return isinstance(x, MatrixGroup_generic)

def MatrixGroup(gens):
    r"""
    Return the matrix group with given generators.

    INPUT:
         gens -- list of matrices in a matrix space or matrix group

    EXAMPLES:
        sage: F = GF(5)
        sage: gens = [matrix(F,2,[1,2, -1, 1]), matrix(F,2, [1,1, 0,1])]
        sage: G = MatrixGroup(gens); G
        Matrix group over Finite Field of size 5 with 2 generators:
         [[[1, 2], [4, 1]], [[1, 1], [0, 1]]]

    In the second example, the generators are a matrix over $\ZZ$, a
    matrix over a finite field, and the integer $2$.  SAGE determines
    that they both canonically map to matrices over the finite field,
    so creates that matrix group there.
        sage: gens = [matrix(2,[1,2, -1, 1]), matrix(GF(7), 2, [1,1, 0,1]), 2]
        sage: G = MatrixGroup(gens); G
        Matrix group over Finite Field of size 7 with 3 generators:
         [[[1, 2], [6, 1]], [[1, 1], [0, 1]], [[2, 0], [0, 2]]]

    Each generator must be invertible:
        sage: G = MatrixGroup([matrix(ZZ,2,[1,2,3,4])])
        Traceback (most recent call last):
        ...
        ValueError: each generator must be an invertible matrix but one is not:
        [1 2]
        [3 4]
    """
    if len(gens) == 0:
        raise ValueError, "gens must have positive length"
    try:
        R = gens[0].base_ring()
    except AttributeError:
        raise TypeError, "gens must be a list of matrices"
    for i in range(len(gens)):
        if not is_Matrix(gens[i]):
            try:
                gens[i] = gens[i].matrix()
            except AttributeError:
                pass
    if is_FiniteField(R):
        return MatrixGroup_gens_finite_field(gens)
    else:
        return MatrixGroup_gens(gens)

class MatrixGroup_gap(MatrixGroup_generic):
    def __init__(self, n, R, var='a'):
        """
        INPUT:
            n -- the degree
            R -- the base ring
            var -- variable used to define field of definition of
                   actual matrices in this group.
        """
        if not is_Ring(R):
            raise TypeError, "R (=%s) must be a ring"%R
        self._var = var
        self.__n = integer.Integer(n)
        self.__R = R

    def __cmp__(self, H):
        if not isinstance(H, MatrixGroup_gap):
            return cmp(type(self), type(H))
        return cmp(gap(self), gap(H))

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
        M = self.matrix_space()(x)
        g = MatrixGroupElement(M, self)
        if not gap(g) in gap(self):
            raise TypeError, "no way to coerce element to self."
        return g

    def _Hom_(self, G, cat=None):
        if not (cat is None or (cat is G.category() and cat is self.category())):
            raise NotImplementedError
        if not is_MatrixGroup(G):
            raise TypeError, "G (=%s) must be a matrix group."%G
        import homset
        return homset.MatrixGroupHomset(self, G)

    def hom(self, x):
        v = Sequence(x)
        U = v.universe()
        if not is_MatrixGroup(U):
            raise TypeError, "u (=%s) must have universe a matrix group."%U
        return self.Hom(U)(x)

    def matrix_space(self):
        """
        Return the matrix space corresponding to this matrix group.

        This is a matrix space over the field of definition of this
        matrix group.

        EXAMPLES:
            sage: F = GF(5); MS = MatrixSpace(F,2,2)
            sage: G = MatrixGroup([MS(1), MS([1,2,3,4])])
            sage: G.matrix_space()
            Full MatrixSpace of 2 by 2 dense matrices over Finite Field of size 5
        """
        try:
            return self.__matrix_space
        except AttributeError:
            pass
        self.__matrix_space = MatrixSpace(self.field_of_definition(), self.__n)
        return self.__matrix_space

    def degree(self):
        """
        Return the degree of this matrix group.

        EXAMPLES:
            sage: SU(5,5).degree()
            5
        """
        return self.__n

    def field_of_definition(self, var='a'):
        """
        Return a field that contains all the matrices in this matrix group.

        EXAMPLES:
            sage: G = SU(3,GF(5))
            sage: G.base_ring()
	    Finite Field of size 5
	    sage: G.field_of_definition()
            Finite Field in a of size 5^2
            sage: G = GO(4,GF(7),1)
            sage: G.field_of_definition()
            Finite Field of size 7
            sage: G.base_ring()
            Finite Field of size 7
        """
        return self.__R

    def base_ring(self):
        """
        Return the base ring of this matrix group.

        EXAMPLES:
            sage: GL(2,GF(3)).base_ring()
            Finite Field of size 3
            sage: G = SU(3,GF(5))
            sage: G.base_ring()
            Finite Field of size 5
            sage: G.field_of_definition()
            Finite Field in a of size 5^2
        """
        return self.__R

    base_field = base_ring

    def is_finite(self):
        """
        Return True if this matrix group is finite.

        EXAMPLES:
            sage: G = GL(2,GF(3))
            sage: G.is_finite()
            True
            sage: SL(2,ZZ).is_finite()
            False
        """
        if self.base_ring().is_finite():
            return True
        return self._gap_().IsFinite().bool()

    def order(self):
        """
        EXAMPLES:
            sage: G = Sp(4,GF(3))
	    sage: G.order()
            51840
            sage: G = SL(4,GF(3))
            sage: G.order()
            12130560
            sage: F = GF(5); MS = MatrixSpace(F,2,2)
            sage: gens = [MS([[1,2],[-1,1]]),MS([[1,1],[0,1]])]
            sage: G = MatrixGroup(gens)
            sage: G.order()
            480
            sage: G = MatrixGroup([matrix(ZZ,2,[1,1,0,1])])
            sage: G.order()
            Infinity
        """
        g = self._gap_()
        if g.IsFinite().bool():
            return integer.Integer(gap(self).Size())
        return infinity

    def gens(self):
        """
        Return generators for this matrix group.

        EXAMPLES:
            sage: G = GO(3,GF(5))
            sage: G.gens()
            [
            [2 0 0]
            [0 3 0]
            [0 0 1],
            [0 1 0]
            [1 4 4]
            [0 2 1]
            ]
        """
        try:
            return self.__gens
        except AttributeError:
            pass
        F = self.field_of_definition()
        gap_gens = list(gap(self).GeneratorsOfGroup())
        gens = Sequence([MatrixGroupElement(g._matrix_(F), self, check=False) for g in gap_gens],
                        cr=True, universe=self, check=False)
        self.__gens = gens
        return gens

    def ngens(self):
        """
        Return the number of generators of this linear group.

        EXAMPLES:
            sage: G = GO(3,GF(5))
            sage: G.ngens()
            2
        """
        return len(self.gens())


    def gen(self, n):
        """
        Return the n-th generator.

        EXAMPLES:
            sage: G = GU(4,GF(5), var='beta')
            sage: G.gen(0)
            [  beta      0      0      0]
            [     0      1      0      0]
            [     0      0      1      0]
            [     0      0      0 3*beta]
        """
        return self.gens()[n]

    def as_matrix_group(self):
        """
        Return this group, but as a general matrix group, i.e., throw
        away the extra structure of general unitary group.

        EXAMPLES:
            sage: G = SU(4,GF(5))
            sage: G.as_matrix_group()
            Matrix group over Finite Field in a of size 5^2 with 2 generators:
            [[[a, 0, 0, 0], [0, 2*a + 3, 0, 0], [0, 0, 4*a + 1, 0], [0, 0, 0, 3*a]], [[1, 0, 4*a + 3, 0], [1, 0, 0, 0], [0, 2*a + 4, 0, 1], [0, 3*a + 1, 0, 0]]]

            sage: G = GO(3,GF(5))
            sage: G.as_matrix_group()
            Matrix group over Finite Field of size 5 with 2 generators:
            [[[2, 0, 0], [0, 3, 0], [0, 0, 1]], [[0, 1, 0], [1, 4, 4], [0, 2, 1]]]
        """
        from sage.groups.matrix_gps.matrix_group import MatrixGroup
        return MatrixGroup([g.matrix() for g in self.gens()])

    def list(self):
        """
        Return list of all elements of this group.

        EXAMPLES:
            sage: F = GF(3)
            sage: gens = [matrix(F,2, [1,0, -1,1]), matrix(F, 2, [1,1,0,1])]
	    sage: G = MatrixGroup(gens)
	    sage: G.order()
	    24
	    sage: v = G.list()
            sage: len(v)
            24
            sage: v[:2]
            [[0 1]
            [2 0], [0 1]
            [2 1]]
            sage: G.list()[0] in G
	    True

            sage: GL(2,ZZ).list()
            Traceback (most recent call last):
            ...
            ValueError: group must be finite
        """
        if not self.is_finite():
            raise ValueError, "group must be finite"
        F = self.field_of_definition()
        X = list(self._gap_().Elements())
        return [MatrixGroupElement(a._matrix_(F), self) for a in X]

class MatrixGroup_gap_finite_field(MatrixGroup_gap):
    def order(self):
        """
        EXAMPLES:
            sage: G = Sp(4,GF(3))
	    sage: G.order()
            51840
            sage: G = SL(4,GF(3))
            sage: G.order()
            12130560
            sage: F = GF(5); MS = MatrixSpace(F,2,2)
            sage: gens = [MS([[1,2],[-1,1]]),MS([[1,1],[0,1]])]
            sage: G = MatrixGroup(gens)
            sage: G.order()
            480
            sage: G = MatrixGroup([matrix(ZZ,2,[1,1,0,1])])
            sage: G.order()
            Infinity
        """
        return integer.Integer(gap(self).Size())

    def random(self):
        """
        Return a random element of this group.

        EXAMPLES:
            sage: G = Sp(4,GF(3))
            sage: G.random()        ## random output
            [0 2 1 1]
            [1 1 1 2]
            [2 2 2 2]
            [0 2 0 2]

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
        from matrix_group_element import MatrixGroupElement
        F = self.field_of_definition()
        return MatrixGroupElement(gap(self).Random()._matrix_(F), self, check=False)

    def __contains__(self, x):
        """
        Return True if $x$ is an element of this abelian group.

        EXAMPLES:
            sage: F = GF(5); MS = MatrixSpace(F,2,2)
            sage: G = MatrixGroup([MS(1), MS([1,2,3,4])])
            sage: G
            Matrix group over Finite Field of size 5 with 2 generators:
             [[[1, 0], [0, 1]], [[1, 2], [3, 4]]]
            sage: G.order()
            8
            sage: G(1)
            [1 0]
            [0 1]
            sage: G.1 in G
            True
            sage: 1 in G
            True
            sage: [1,2,3,4] in G
            True
            sage: matrix(GF(5),2,[1,2,3,5]) in G
            False
            sage: G(matrix(GF(5),2,[1,2,3,5]))
            Traceback (most recent call last):
            ...
            TypeError: no way to coerce element to self.
        """
        from matrix_group_element import MatrixGroupElement
        if isinstance(x, MatrixGroupElement):
            if x.parent() == self:
                return True
            return gap(x) in gap(self)
        try:
            self(x)
            return True
        except TypeError:
            return False

    def conjugacy_class_representatives(self):
        """
        Return a set of representatives for each of the conjugacy
        classes of the group.

        EXAMPLES:
            sage: G = SU(3,GF(2))
            sage: G.conjugacy_class_representatives()
            [
            [1 0 0]
            [0 1 0]
            [0 0 1],
            [a 0 0]
            [0 a 0]
            [0 0 a],
            [a + 1     0     0]
            [    0 a + 1     0]
            [    0     0 a + 1]
            ]
            sage: GL(2,GF(3)).conjugacy_class_representatives()
            [
            [1 0]
            [0 1],
            [2 0]
            [0 2]
            ]
            sage: GU(2,GF(5)).conjugacy_class_representatives()
            [
            [1 0]
            [0 1],
            [4 0]
            [0 4],
            [2*a + 2       0]
            [      0 2*a + 2],
            [2*a + 1       0]
            [      0 2*a + 1],
            [3*a + 3       0]
            [      0 3*a + 3],
            [3*a + 4       0]
            [      0 3*a + 4]
            ]
        """
        try:
            return self.__reps
        except AttributeError:
            pass
        from matrix_group_element import MatrixGroupElement
        G    = self._gap_().Center().ConjugacyClasses()
        reps = list(gap.List(G, 'x -> Representative(x)'))
        F    = self.field_of_definition()
        self.__reps = Sequence([self(g._matrix_(F)) for g in reps], cr=True, universe=self, check=False)
        return self.__reps

    def center(self):
        """
        Return the center of this linear group as a matrix group.

        EXAMPLES:
            sage: G = SU(3,GF(2))
            sage: G.center()
            Matrix group over Finite Field in a of size 2^2 with 1 generators:
             [[[a, 0, 0], [0, a, 0], [0, 0, a]]]
            sage: GL(2,GF(3)).center()
            Matrix group over Finite Field of size 3 with 1 generators:
             [[[2, 0], [0, 2]]]
            sage: GL(3,GF(3)).center()
            Matrix group over Finite Field of size 3 with 1 generators:
             [[[2, 0, 0], [0, 2, 0], [0, 0, 2]]]
            sage: GU(3,GF(2)).center()
            Matrix group over Finite Field in a of size 2^2 with 1 generators:
             [[[a + 1, 0, 0], [0, a + 1, 0], [0, 0, a + 1]]]
        """
        try:
            return self.__center
        except AttributeError:
            pass
        G = list(self._gap_().Center().GeneratorsOfGroup())
        F = self.field_of_definition()
        from sage.groups.matrix_gps.matrix_group import MatrixGroup
        self.__center = MatrixGroup([g._matrix_(F) for g in G])
        return self.__center


class MatrixGroup_gens(MatrixGroup_gap):
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
        v = Sequence(gensG, immutable=True)
        M = v.universe()
        if not is_MatrixSpace(M):
            raise TypeError, "universe of sequence (=%s) of generators must be a matrix space"%M
        if M.nrows() != M.ncols():
            raise ValueError, "matrices must be square."
        for x in v:
            if not x.is_invertible():
                raise ValueError, "each generator must be an invertible matrix but one is not:\n%s"%x
        self._gensG = v
        MatrixGroup_gap.__init__(self, M.nrows(), M.base_ring())

    def gens(self):
        """
        EXAMPLES:
            sage: F = GF(3); MS = MatrixSpace(F,2,2)
            sage: gens = [MS([[1,0],[0,1]]),MS([[1,1],[0,1]])]
	    sage: G = MatrixGroup(gens)
            sage: gens[0] in G
            True
            sage: gens = G.gens()
            sage: gens[0] in G
            True
            sage: gens = [MS([[1,0],[0,1]]),MS([[1,1],[0,1]])]

            sage: F = GF(5); MS = MatrixSpace(F,2,2)
            sage: G = MatrixGroup([MS(1), MS([1,2,3,4])])
            sage: G
            Matrix group over Finite Field of size 5 with 2 generators:
             [[[1, 0], [0, 1]], [[1, 2], [3, 4]]]
            sage: G.gens()
            [[1 0]
            [0 1], [1 2]
            [3 4]]
        """
        try:
            return self.__gens
        except AttributeError:
            t = Sequence([MatrixGroupElement(x, self) for x in self._gensG],
                         immutable=True, universe=self)
            self.__gens = t
            return t

    def _gap_init_(self):
        """
        Returns the string representation of the corresponding GAP command.

        EXAMPLES:
            sage: F = GF(5); MS = MatrixSpace(F,2,2)
            sage: gens = [MS([[1,2],[-1,1]]),MS([[1,1],[0,1]])]
            sage: G = MatrixGroup(gens)
            sage: G._gap_init_()
            'Group([[Z(5)^0,Z(5)^1],[Z(5)^2,Z(5)^0]],[[Z(5)^0,Z(5)^0],[0*Z(5),Z(5)^0]])'
        """
        gens_gap = ','.join([x._gap_init_() for x in self._gensG])
        return 'Group(%s)'%gens_gap

    def _repr_(self):
        """
        EXAMPLES:
            sage: F = GF(5); MS = MatrixSpace(F,2,2)
            sage: gens = [MS([[1,2],[-1,1]]),MS([[1,1],[0,1]])]
            sage: G = MatrixGroup(gens)
            sage: G
            Matrix group over Finite Field of size 5 with 2 generators:
             [[[1, 2], [4, 1]], [[1, 1], [0, 1]]]
        """
        gns = [x.list() for x in self.gens()]
        return "Matrix group over %s with %s generators: \n %s"%(self.base_ring(), self.ngens(), gns)

    def _latex_(self):
        """
        EXAMPLES:
            sage: F = GF(5); MS = MatrixSpace(F,2,2)
            sage: gens = [MS([[1,2],[-1,1]]),MS([[1,1],[0,1]])]
            sage: G = MatrixGroup(gens)
            sage: latex(G)
            \left\langle \left(\begin{array}{rr}
            1&2\\
            4&1
            \end{array}\right), \left(\begin{array}{rr}
            1&1\\
            0&1
            \end{array}\right) \right\rangle
        """
        gens = ', '.join([latex(x) for x in self.gens()])
        return '\\left\\langle %s \\right\\rangle'%gens

class MatrixGroup_gens_finite_field(MatrixGroup_gens, MatrixGroup_gap_finite_field):
    pass


##     def conjugacy_class_representatives_gap(self):
##         """
##         Wraps GAP Representative+ConjugactClasses but returns a list of
##         strings representing the GAP matrices which form a complete
##         set of representatives of the conjugacy classes of the group.

##         EXAMPLES:
##             sage: F = GF(3); MS = MatrixSpace(F,2,2)
##             sage: gens = [MS([[1,0],[-1,1]]),MS([[1,1],[0,1]])]
##             sage: G = MatrixGroup(gens)
##             sage: G.conjugacy_class_representatives_gap()
## 	    ['[ [ Z(3)^0, 0*Z(3) ], [ 0*Z(3), Z(3)^0 ] ]',
##  	    '[ [ 0*Z(3), Z(3)^0 ], [ Z(3), Z(3)^0 ] ]',
##             '[ [ 0*Z(3), Z(3)^0 ], [ Z(3), Z(3) ] ]',
##             '[ [ 0*Z(3), Z(3) ], [ Z(3)^0, Z(3)^0 ] ]',
##      	    '[ [ 0*Z(3), Z(3) ], [ Z(3)^0, Z(3) ] ]',
##      	    '[ [ 0*Z(3), Z(3)^0 ], [ Z(3), 0*Z(3) ] ]',
##  	    '[ [ Z(3), 0*Z(3) ], [ 0*Z(3), Z(3) ] ]']

##         AUTHOR: David Joyner (1-2006)
##         """
##         F = self.__R
##         deg = self.__n
##         gp_gens = self.gens()
##         L = [gap(A) for A in gp_gens]
##         sL = ','.join(str(x) for x in L)
##         if is_FiniteField(F):
##             q = F.order()
##             gap.eval("cl:=ConjugacyClasses(Group(["+sL+"]))")
##             m = eval(gap.eval("Length(cl)"))
##             gap.eval("reps:=List(cl,x->Representative(x))")
##             sreps = [gap.eval("reps["+str(i+1)+"]") for i in range(m)]
##             return sreps
##         raise TypeError, "R (=%s) must be a finite field"%R



