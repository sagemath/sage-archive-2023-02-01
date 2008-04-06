"""
Matrix Group Elements

AUTHORS:
   David Joyner -- initial version
   David Joyner -- (2006-05) various modifications to address William
                   Stein's TODO's.
   William Stein (2006-12-09): many revisions.

EXAMPLES:

    sage: F = GF(3); MS = MatrixSpace(F,2,2)
    sage: gens = [MS([[1,0],[0,1]]),MS([[1,1],[0,1]])]
    sage: G = MatrixGroup(gens); G
    Matrix group over Finite Field of size 3 with 2 generators:
     [[[1, 0], [0, 1]], [[1, 1], [0, 1]]]
    sage: g = G([[1,1],[0,1]])
    sage: h = G([[1,2],[0,1]])
    sage: g*h
    [1 0]
    [0 1]

You cannot add two matrices, since this is not a group
operation.  You can coerce matrices back to the
matrix space and add them there:
    sage: g + h
    Traceback (most recent call last):
    ...
    TypeError: unsupported operand type(s) for +: 'MatrixGroupElement' and 'MatrixGroupElement'

    sage: g.matrix() + h.matrix()
    [2 0]
    [0 2]

    sage: 2*g
    Traceback (most recent call last):
    ...
    TypeError: unsupported operand parent(s) for '*': 'Integer Ring' and 'Matrix group over Finite Field of size 3 with 2 generators:
     [[[1, 0], [0, 1]], [[1, 1], [0, 1]]]'
    sage: 2*g.matrix()
    [2 2]
    [0 2]
"""

#*****************************************************************************
#       Copyright (C) 2006 David Joyner and William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.groups.group import Group
from sage.rings.all import IntegerRing, is_Ring, Integer
from sage.interfaces.gap import gap
import sage.structure.element as element
from sage.matrix.matrix import Matrix
from sage.structure.factorization import Factorization

def is_MatrixGroupElement(x):
    return isinstance(x, MatrixGroupElement)

class MatrixGroupElement(element.MultiplicativeGroupElement):
    """
    An element of a matrix group.

    EXAMPLES:
	sage: F = GF(3); MS = MatrixSpace(F,2,2)
        sage: gens = [MS([[1,0],[0,1]]),MS([[1,1],[0,1]])]
	sage: G = MatrixGroup(gens)
	sage: g = G.random()
	sage: type(g)
	<class 'sage.groups.matrix_gps.matrix_group_element.MatrixGroupElement'>
    """
    def __init__(self, g, parent, check = True):
        r"""
        Create element of a matrix group.

        INPUT:
            g -- a Matrix
            parent -- defines parent group (g must be in parent or a TypeError is raised).
            check -- bool (default: True), if true does some type checking.
        """
        if check:
            if not isinstance(g, Matrix):
                raise TypeError, "g must be a matrix"
            if not (g.parent() is parent):
                g = parent.matrix_space()(g)
        element.Element.__init__(self, parent)
        self.__mat = g
        self.__gap = gap(g)

    def matrix(self):
        """
        Obtain the usual matrix (as an element of a matrix space)
        associated to this matrix group element.

        One reason to compute the associated matrix is that matrices
        support a huge range of functionality.

        EXAMPLES:
            sage: k = GF(7); G = MatrixGroup([matrix(k,2,[1,1,0,1]), matrix(k,2,[1,0,0,2])])
            sage: g = G.0
            sage: g.matrix()
            [1 1]
            [0 1]
            sage: parent(g.matrix())
            Full MatrixSpace of 2 by 2 dense matrices over Finite Field of size 7

        Matrices have extra functionality that matrix group elements do not have.
            sage: g.matrix().charpoly('t')
            t^2 + 5*t + 1
        """
        return self.__mat.copy()

    def _gap_init_(self):
        """
        Return a string representation of a Gap object corresponding to this matrix group element.

        EXAMPLES:
            sage: k = GF(7); G = MatrixGroup([matrix(k,2,[1,1,0,1]), matrix(k,2,[1,0,0,2])]); g = G.1
            sage: g._gap_init_() # The variable $sage27 belongs to gap(k) and is somehow random
            '[[Z(7)^0,0*Z(7)],[0*Z(7),Z(7)^2]]*One($sage27)'
            sage: gap(g._gap_init_())
            [ [ Z(7)^0, 0*Z(7) ], [ 0*Z(7), Z(7)^2 ] ]

        It may be better to use gap(the matrix), since the result is cached.
            sage: gap(G.1)
            [ [ Z(7)^0, 0*Z(7) ], [ 0*Z(7), Z(7)^2 ] ]
            sage: gap(G.1).IsMatrix()
            true
        """
        return self.__mat._gap_init_()

    def _gap_latex_(self):
        r"""
        Return the GAP latex version of this matrix.

        AUTHOR: S. Kohl wrote the GAP function.

        EXAMPLES:
            sage: F = GF(3); MS = MatrixSpace(F,2,2)
            sage: gens = [MS([[1,0],[0,1]]),MS([[1,1],[0,1]])]
            sage: G = MatrixGroup(gens)
            sage: g = G([[1, 1], [0, 1]])
            sage: print g._gap_latex_()
            \left(\begin{array}{rr}%
            Z(3)^{0}&Z(3)^{0}\\%
            0*Z(3)&Z(3)^{0}\\%
            \end{array}\right)%
            <BLANKLINE>

        Type view(g._latex_()) to see the object in an xdvi window (assuming you
        have latex and xdvi installed).
        """
        s1 = self.__mat._gap_init_()
        s2 = gap.eval("LaTeX("+s1+")")
        return eval(s2)

    def _repr_(self):
        """
        Return string representation of this matrix.

        EXAMPLES:
            sage: F = GF(3); MS = MatrixSpace(F,2,2)
            sage: gens = [MS([[1,0],[0,1]]),MS([[1,1],[0,1]])]
            sage: G = MatrixGroup(gens)
            sage: g = G([[1, 1], [0, 1]])
            sage: g
            [1 1]
            [0 1]
        """
        return str(self.__mat)

    def _latex_(self):
        r"""
        EXAMPLES:
            sage: F = GF(3); MS = MatrixSpace(F,2,2)
            sage: gens = [MS([[1,0],[0,1]]),MS([[1,1],[0,1]])]
            sage: G = MatrixGroup(gens)
            sage: g = G([[1, 1], [0, 1]])
            sage: print g._latex_()
            \left(\begin{array}{rr}
            1 & 1 \\
            0 & 1
            \end{array}\right)

        Type \code{view(g._latex_())} to see the object in an xdvi window (assuming you
        have latex and xdvi installed).
        """
        return self.__mat._latex_()

    def _mul_(self,other):
        """
        Return the product of self and other, which must have
        identical parents.

        EXAMPLES:
            sage: F = GF(3); MS = MatrixSpace(F,2)
            sage: gens = [MS([1,0, 0,1]), MS([1,1, 0,1])]
            sage: G = MatrixGroup(gens)
            sage: g = G([1,1, 0,1])
            sage: h = G([1,1, 0,1])
            sage: g*h
            [1 2]
            [0 1]
        """
        return MatrixGroupElement(self.__mat * other.__mat, self.parent(), check=False)

    def __invert__(self):
        return MatrixGroupElement(~self.__mat, self.parent(), check=False)

    def __cmp__(self, other):
        return cmp(self.__mat, other.__mat)

    def order(self):
        """
        Return the order of this group element, which is the smallest
        positive integer $n$ such that $g^n = 1$.

        EXAMPLES:
            sage: k = GF(7); G = MatrixGroup([matrix(k,2,[1,1,0,1]), matrix(k,2,[1,0,0,2])])
            sage: G
            Matrix group over Finite Field of size 7 with 2 generators:
             [[[1, 1], [0, 1]], [[1, 0], [0, 2]]]
            sage: G.order()
            21
        """
        try:
            return self.__order
        except AttributeError:
            self.__order = Integer(self._gap_().Order())
            return self.__order

    def word_problem(self, words=None):
        r"""
        Right this group element in terms of the elements of the list
        \code{words}.

        If G and H are permutation groups (with G the parent of self),
        H is a subgroup of G generated by a list \code{words} of
        elements of G. If g is in H, return the expression for g as a
        word in the elements of \code{words}.

        ALGORITHM: Use GAP, which has optimized algorithms for solving
        the word problem (the GAP functions EpimorphismFromFreeGroup
        and PreImagesRepresentative).

        INPUT:
            words -- list (default: None) a list of elements
                     of the parent group; if words is empty, uses
                     the gens of the parent group.
        EXAMPLE:
            sage: G = GL(2,5); G
            General Linear Group of degree 2 over Finite Field of size 5
            sage: G.gens()
            [
            [2 0]
            [0 1],
            [4 1]
            [4 0]
            ]
            sage: G(1).word_problem([G.1, G.0])
            1

        Next we construct a more complicated element of the group from
        the generators:
            sage: (G.0).order(), (G.1).order()
            (4, 3)
            sage: g = G.0^3 * G.1^2

        We then ask to solve the word problem, with the generators
        reversed (just to make it trickier):
            sage: s = g.word_problem([G.1, G.0]); s
            ([2 0]
            [0 1])^-1 *
            ([4 1]
            [4 0])^-1

        It worked!
            sage: s.prod() == g
            True

        AUTHORS: David Joyner and William Stein
        """
        G = self.parent()
        gg = gap(G)
        gens = gg.GeneratorsOfGroup()
        H = gens.Group()
        hom = H.EpimorphismFromFreeGroup()
        in_terms_of_gens = str(hom.PreImagesRepresentative(self))
        if 'identity' in in_terms_of_gens:
            return Factorization([])
        v = list(H.GeneratorsOfGroup())
        F = G.field_of_definition()
        w = [G(x._matrix_(F)) for x in v]
        factors = [Factorization([(x,1)]) for x in w]
        d = dict([('x%s'%(i+1),factors[i]) for i in range(len(factors))])

        from sage.misc.sage_eval import sage_eval
        F = sage_eval(str(in_terms_of_gens), d)
        F._set_cr(True)
        return F

    def list(self):
        """
        Return list representation of this matrix.

        EXAMPLES:
            sage: F = GF(3); MS = MatrixSpace(F,2,2)
            sage: gens = [MS([[1,0],[0,1]]),MS([[1,1],[0,1]])]
            sage: G = MatrixGroup(gens)
            sage: g = G.0
            sage: g.list()
            [[1, 0], [0, 1]]
        """
        rws = (self.__mat).rows()
        lrws = [list(x) for x in rws]
        return lrws


