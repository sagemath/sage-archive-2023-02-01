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
    TypeError: unsupported operand parent(s) for '*': 'Integer Ring' and 'MatrixGroup( [[[1, 0], [0, 1]], [[1, 1], [0, 1]]] )'
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
        Return a string that when evaluated in GAP defines this matrix group element.

        EXAMPLES:
            sage: k = GF(7); G = MatrixGroup([matrix(k,2,[1,1,0,1]), matrix(k,2,[1,0,0,2])]); g = G.1
            sage: g._gap_init_()
            '[[Z(7)^0,0*Z(7)],[0*Z(7),Z(7)^2]]'
            sage: gap(g._gap_init_())
            [ [ Z(7)^0, 0*Z(7) ], [ 0*Z(7), Z(7)^2 ] ]

        It may be better to use gap(the matrix), since the result is cached.
            sage: gap(G.1)
            [ [ Z(7)^0, 0*Z(7) ], [ 0*Z(7), Z(7)^2 ] ]
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
            1&1\\
            0&1
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

    def word_problem(g, words, display=True):
        """
        G and H are permutation groups, g in G, H is a subgroup of G generated by a
        list (words) of elements of G. If g is in H, return the expression
        for g as a word in the elements of (words).

        ALGORITHM: This function does not solve the word problem in
        SAGE. Rather it pushes it over to GAP, which has optimized
        algorithms for the word problem. Essentially, this function is
        a wrapper for the GAP functions "EpimorphismFromFreeGroup" and
        "PreImagesRepresentative".

        EXAMPLE:
        sage: ?
        """
        import copy
        from sage.groups.perm_gps.permgroup import PermutationGroup
        from sage.interfaces.all import gap
        G = g.parent()
        gap.eval("l:=One(Rationals)")
        s1 = "gens := GeneratorsOfGroup(%s)"%G._gap_init_()
        gap.eval(s1)
        s2 = "g0:=%s; gensH:=%s"%(gap(g),words)
        gap.eval(s2)
        s3 = 'G:=Group(gens); H:=Group(gensH)'
        gap.eval(s3)
        phi = gap.eval("hom:=EpimorphismFromFreeGroup(G)")
        l1 = gap.eval("ans:=PreImagesRepresentative(hom,g0)")
        l2 = l1
        l4 = []
        l3 = l1.split("*")
        for i in range(1,len(words)+1):
            l2 = l2.replace("x"+str(i),str(words[i-1]))
        if display:
            for i in range(len(l3)):    ## parsing the word for display
                if len(l3[i].split("^"))==2:
                    l4.append([l3[i].split("^")[0],int(l3[i].split("^")[1])])
                if len(l3[i].split("^"))==1:
                    l4.append([l3[i].split("^")[0],1])
            l5 = copy.copy(l4)
            for i in range(len(l4)):
                for j in range(1,len(words)+1):
                    l5[i][0] = l5[i][0].replace("x"+str(j),str(words[j-1]))
            print "         ",l1
            print "         ",l5
        return l2

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


