"""
To contain general base classes for elements of the classical groups
GL, SL, Sp.

AUTHORS:
   David Joyner -- initial version
   David Joyner -- (2006-05) various modifications to address William Stein's
                   TODO's.

TODO: I don't like this at all.
  The matrix group element class should derive from matrix, so that
  all functionality of matrices (e.g., multiplication, charpoly, det)
  is automatically available.  Certain functionality should return an
  element not in the group anymore, e.g., +.
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
from linear import LinearGroup_generic
import sage.structure.element as element
from general_linear import GeneralLinearGroup_generic
from sage.matrix.matrix_space import MatrixSpace, is_MatrixSpace
from sage.matrix.matrix_generic import Matrix

def is_MatrixGroupElement(x):
    return isinstance(x, MatrixGroupElement)

class MatrixGroupElement(element.Element_cmp_,
                              element.MultiplicativeGroupElement):
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
            g -- defines element
            parent (optional) -- defines parent group (g must be in parent if
                                 specified, or a TypeError is raised).
            check -- bool (default: True), if False assumes g is a
                     gap element in parent (if specified).

        EXAMPLES:

        """
        from sage.groups.matrix_gps.matrix_group import MatrixGroup
        if check:
            if not (parent is None or isinstance(parent, MatrixGroup)):
                raise TypeError, 'parent must be a matrix group'
            if not (g.parent() is parent):
                g = parent.matrix_space()(g)
        element.Element.__init__(self, parent)
        self.__mat = g
        self.__gap = gap(g)

    def _gap_init_(self):
        return self.__mat._gap_init_()

    def _gap_latex_(self):
        """
        wraps GAP's LaTeX command (written by S. Kohl)

        EXAMPLES:
            sage: F = GF(3); MS = MatrixSpace(F,2,2)
            sage: gens = [MS([[1,0],[0,1]]),MS([[1,1],[0,1]])]
            sage: G = MatrixGroup(gens)
            sage: g = G([[1, 1], [0, 1]])
            sage: g._gap_latex_()
            '\\left(\\begin{array}{rr}%\nZ(3)^{0}&Z(3)^{0}\\\\%\n0*Z(3)&Z(3)^{0}\\\\%\n\\end{array}\\right)%\n'

        Type view(g._latex_()) to see the object in an xdvi window (assuming you
        have latex and xdvi installed).
        """
        s1 = self.__mat._gap_init_()
        s2 = gap.eval("LaTeX("+s1+")")
        return eval(s2)
        #s2 = s2.replace("\\\\left(\\\\","\left(\ ")
        #s3 = s2.replace("\\\\left(\\","\left( \n ")
        #s4 = s3.replace("\\\\right)","\right) ")
        #s5 = s4.replace("\\\\\\"," ")
        #s6 = s5.replace("\\\\end","\end")

    def __repr__(self):
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

    def __str__(self):
        """
        Return string representation of this matrix (called by print).

        EXAMPLES:
            sage: F = GF(3); MS = MatrixSpace(F,2,2)
            sage: gens = [MS([[1,0],[0,1]]),MS([[1,1],[0,1]])]
            sage: G = MatrixGroup(gens)
            sage: g = G([[1, 1], [0, 1]])
            sage: print g
            [[1, 1], [0, 1]]

        """
        rws = (self.__mat).rows()
        lrws = [list(x) for x in rws]
        return str(lrws)

    def _latex_(self):
        """
        EXAMPLES:
            sage: F = GF(3); MS = MatrixSpace(F,2,2)
            sage: gens = [MS([[1,0],[0,1]]),MS([[1,1],[0,1]])]
            sage: G = MatrixGroup(gens)
            sage: g = G([[1, 1], [0, 1]])
            sage: g._latex_()
            '\\left(\\begin{array}{rr}\n1&1\\\\\n0&1\n\\end{array}\\right)'

        Type view(g._latex_()) to see the object in an xdvi window (assuming you
        have latex and xdvi installed).
        """
        return self.__mat._latex_()

    def __mul__(self,other):
        """

        EXAMPLES:
            sage: F = GF(3); MS = MatrixSpace(F,2,2)
            sage: gens = [MS([[1,0],[0,1]]),MS([[1,1],[0,1]])]
            sage: G = MatrixGroup(gens)
            sage: g = G([[1,1],[0,1]])
            sage: h = G([[1,1],[0,1]])
            sage: g*h
            [1 2]
            [0 1]

        """
        m = self.parent().degree()
        n = other.parent().degree()
        if m != n:
            raise TypeError, "Matrices must have compatible dimensions"
        F = self.parent().base_ring()
        M = MatrixSpace(F,n,n)
        x = M(self.__mat)*M(other.__mat)
        return MatrixGroupElement((self.parent()).matrix_space()(x), self.parent())

    def __add__(self,other):
        """
        You can add and scalar multiply even though it is not a group
        operation.

        EXAMPLES:
            sage: F = GF(3); MS = MatrixSpace(F,2,2)
            sage: gens = [MS([[1,0],[0,1]]),MS([[1,1],[0,1]])]
            sage: G = MatrixGroup(gens)
            sage: g = G([[1,1],[0,1]])
            sage: h = G([[1,2],[0,1]])
            sage: g+h
            [2 0]
            [0 2]
            sage: 2*g
            [2 2]
            [0 2]
            sage: MS = MatrixSpace(QQ,2,2)
            sage: gens = [MS([[1,0],[0,1]]),MS([[1,1],[0,1]])]
            sage: G = MatrixGroup(gens)
            sage: g = G([[1,1],[0,1]])
            sage: h = G([[1,2],[0,1]])
            sage: g+h
            [2 3]
            [0 2]
            sage: 5/3*h
            [ 5/3 10/3]
            [   0  5/3]

        """

        m = self.parent().degree()
        n = other.parent().degree()
        if m != n:
            raise TypeError, "Matrices must have compatible dimensions"
        F = self.parent().base_ring()
        M = MatrixSpace(F,n,n)
        return M(self.__mat)+M(other.__mat)

    def order(self):
        """
        Return the order of this group element, which is the smallest
        positive integer $n$ for which $g^n = 1$.

        EXAMPLES:
        """
        return Integer(self._gap_().Order())

    def word_problem(words, g, display=True):
        """
        G and H are permutation groups, g in G, H is a subgroup of G generated by a
        list (words) of elements of G. If g is in H, return the expression
        for g as a word in the elements of (words).

        This function does not solve the word problem in SAGE. Rather
        it pushes it over to GAP, which has optimized algorithms for
        the word problem. Essentially, this function is a wrapper for the GAP
        functions "EpimorphismFromFreeGroup" and "PreImagesRepresentative".

        EXAMPLE:
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
            sage: g = G.random()
            sage.: g.list()   # random output
            [[1, 1], [0, 1]]

        """
        rws = (self.__mat).rows()
        lrws = [list(x) for x in rws]
        return lrws


