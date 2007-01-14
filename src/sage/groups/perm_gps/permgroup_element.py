"""
Permutation group elements

AUTHORS:
    - David Joyner (2006-02)
    - David Joyner (2006-03), word problem method and reorganization

EXAMPLES:
The Rubik's cube group:
    sage: f= [(17,19,24,22),(18,21,23,20),(6,25,43,16),(7,28,42,13),(8,30,41,11)]
    sage: b=[(33,35,40,38),(34,37,39,36),( 3, 9,46,32),( 2,12,47,29),( 1,14,48,27)]
    sage: l=[( 9,11,16,14),(10,13,15,12),( 1,17,41,40),( 4,20,44,37),( 6,22,46,35)]
    sage: r=[(25,27,32,30),(26,29,31,28),( 3,38,43,19),( 5,36,45,21),( 8,33,48,24)]
    sage: u=[( 1, 3, 8, 6),( 2, 5, 7, 4),( 9,33,25,17),(10,34,26,18),(11,35,27,19)]
    sage: d=[(41,43,48,46),(42,45,47,44),(14,22,30,38),(15,23,31,39),(16,24,32,40)]
    sage: cube = PermutationGroup([f,b,l,r,u,d])
    sage: F=cube.gens()[0]
    sage: B=cube.gens()[1]
    sage: L=cube.gens()[2]
    sage: R=cube.gens()[3]
    sage: U=cube.gens()[4]
    sage: D=cube.gens()[5]
    sage: cube.order()
    43252003274489856000
    sage: F.order()
    4

The interested user may wish to explore the following commands:
move = cube.random() and time word_problem([F,B,L,R,U,D], move, False).
This typically takes about 5 minutes (on a 2 Ghz machine) and outputs
a word ('solving' the cube in the position move) with about 60 terms
or so.

"""

###########################################################################
#  Copyright (C) 2006 William Stein <wstein@gmail.com>
#  Copyright (C) 2006 David Joyner
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################


import random

import sage.structure.element as element
import sage.groups.group as group

from sage.rings.all      import ZZ, Integer, is_MPolynomial, MPolynomialRing, is_Polynomial
from sage.matrix.all     import MatrixSpace
from sage.interfaces.all import gap, is_GapElement, is_ExpectElement

import sage.structure.coerce as coerce

import operator

from sage.rings.integer import Integer
from sage.structure.element import MonoidElement
from sage.rings.arith import *   # todo: get rid of this -- "from blah import *" is evil.

import permgroup

def is_PermutationGroupElement(x):
    return isinstance(x, PermutationGroupElement)


def gap_format(x):
    """
    Put a permutation in Gap format, as a string.
    """
    x = str(x).replace(' ','')
    return x.replace('),(',')(').replace('[','').replace(']','')

class PermutationGroupElement(element.MultiplicativeGroupElement):
    """
    An element of a permutation group.

    EXAMPLES:
        sage: G = PermutationGroup(['(1,2,3)(4,5)'])
        sage: G
        Permutation Group with generators [(1,2,3)(4,5)]
        sage: g = G.random()
        sage: g in G
        True
        sage: g = G.gen(0); g
        (1,2,3)(4,5)
        sage: print g
        (1,2,3)(4,5)
        sage: g*g
        (1,3,2)
        sage: g**(-1)
        (1,3,2)(4,5)
        sage: g**2
        (1,3,2)
        sage: G = PermutationGroup([(1,2,3)])
        sage: g = G.gen(0); g
        (1,2,3)
        sage: g.order()
        3

    This example illustrates how permutations act on multivariate
    polynomials.

        sage: R = MPolynomialRing(RationalField(), 5, ["x","y","z","u","v"])
        sage: x, y, z, u, v = R.gens()
        sage: f = x**2 - y**2 + 3*z**2
        sage: G = PermutationGroup(['(1,2,3)(4,5)', '(1,2,3,4,5)'])
        sage: sigma = G.gen(0)
        sage: f * sigma
        -1*z^2 + y^2 + 3*x^2

    """
    def __init__(self, g, parent = None, check = True):
        r"""
        Create element of a permutation group.

        There are several ways to define a permutation group element:
        \begin{itemize}
        \item  Define a permutation group $G$, then use
               \code{G.gens()} and multiplication * to construct
               elements.
        \item Define a permutation group $G$, then use e.g.,
               \code{G([(1,2),(3,4,5)])} to construct an element of
               the group.  You could also use \code{G('(1,2)(3,4,5)')}
        \item Use e.g., \code{PermutationGroupElement([(1,2),(3,4,5)])}
        or \code{PermutationGroupElement('(1,2)(3,4,5)')}
              to make a permutation group element with parent $S_5$.
        \end{itemize}

        INPUT:
            g -- defines element
            parent (optional) -- defines parent group (g must be in parent if
                                 specified, or a TypeError is raised).
            check -- bool (default: True), if False assumes g is a
                     gap element in parent (if specified).

        EXAMPLES:
        We illustrate construction of permutation using several
        different methods.

        First we construct elements by multiplying together generators
        for a group.

            sage: G = PermutationGroup(['(1,2)(3,4)', '(3,4,5,6)'])
            sage: s = G.gens()
            sage: s[0]
            (1,2)(3,4)
            sage: s[1]
            (3,4,5,6)
            sage: s[0]*s[1]
            (1,2)(3,5,6)
            sage: (s[0]*s[1]).parent()
            Permutation Group with generators [(1,2)(3,4), (3,4,5,6)]

        Next we illustrate creation of a permutation using
        coercion into an already-created group.

            sage: g = G([(1,2),(3,5,6)])
            sage: g
            (1,2)(3,5,6)
            sage: g.parent()
            Permutation Group with generators [(1,2)(3,4), (3,4,5,6)]
            sage: g == s[0]*s[1]
            True

        We can also use a string instead of a list to specify
        the permutation.

            sage: h = G('(1,2)(3,5,6)')
            sage: g == h
            True

        We can also make a permutation group element directly
        using the \code{PermutationGroupElement} command.  Note
        that the parent is then the full symmetric group $S_n$,
        where $n$ is the largest integer that is moved by the
        permutation.

            sage: k = PermutationGroupElement('(1,2)(3,5,6)')
            sage: k
            (1,2)(3,5,6)
            sage: k.parent()
            Symmetric group of order 6! as a permutation group

        Note the comparison of permutations doesn't require that the
        parent groups are the same.

            sage: k == g
            True

        Arithmetic with permutations having different parents is also defined:

            sage: k*g
            (3,6,5)
            sage: (k*g).parent()
            Symmetric group of order 6! as a permutation group

            sage: G = PermutationGroup([[(1,2,3),(4,5)],[(3,4)]])
            sage: loads(dumps(G.0)) == G.0
            True

        EXAMPLES:
            sage: k = PermutationGroupElement('(1,2)(3,5,6)')
            sage: k._gap_()
            (1,2)(3,5,6)
            sage: k._gap_().parent()
            Gap
        """
        from sage.groups.perm_gps.permgroup import PermutationGroup_generic
        if check:
            if not (parent is None or isinstance(parent, PermutationGroup_generic)):
                raise TypeError, 'parent must be a permutation group'
            self.__gap = gap_format(g)
            if not parent is None:
                P = parent._gap_()
                if not self._gap_(P.parent()) in P:
                    raise TypeError, 'permutation %s not in %s'%(self.__gap, parent)
        else:
            self.__gap = str(g)
        if parent is None:
            parent = permgroup.SymmetricGroup(self._gap_().LargestMovedPoint())
        element.Element.__init__(self, parent)

    def _gap_init_(self):
        return self.__gap


    def _repr_(self):
        """
        Return string representation of this permutation.

        EXAMPLES:
        We create the permutation $(1,2,3)(4,5)$ and print it.

            sage: g = PermutationGroupElement([(1,2,3),(4,5)])
            sage: g._repr_()
            '(1,2,3)(4,5)'

        Permutation group elements support renaming them so
        they print however you want, as illustred below:

            sage: g.rename('sigma')
            sage: g
            sigma
            sage: g.rename()
            sage: g
            (1,2,3)(4,5)
        """
        return self.__gap

    def _latex_(self):
        return str(self)

    def __getitem__(self, i):
        """
        Return the ith permutation cycle in the disjoint cycle
        representation of self.

        INPUT:
            i -- integer

        OUTPUT:
            a permutation group element

        EXAMPLE:
            sage: G = PermutationGroup([[(1,2,3),(4,5)]],5)
            sage: g = G.gen(0)
            sage: g[0]
            (1,2,3)
            sage: g[1]
            (4,5)
        """
        S = str(self).split(')(')
        if i >= 0 and i < len(S):
            T = S[i]
            if i > 0:
                T = '(' + T
            if i < len(S)-1:
                T += ')'
        else:
            raise IndexError, "i (=%s) must be between 0 and %s, inclusive"%(i, len(S)-1)
        return PermutationGroupElement(gap(T), check = False)

    def __cmp__(self, right):
        """
        Compare group elements self and right.

        EXAMPLES:
            sage: G = PermutationGroup([[(1,2,3),(4,5)],[(3,4)]])
            sage: G.gen(0) < G.gen(1)
            False
            sage: G.gen(0) > G.gen(1)
            True
        """
        r = right._gap_()
        G = r.parent()
        l = self._gap_(G)
        return cmp(l,r)

    def __call__(self, i):
        """
        Returns the image of the integer i under this permutation.

        EXAMPLE:
            sage: G = PermutationGroup(['(1,2,3)(4,5)'])
            sage: G
            Permutation Group with generators [(1,2,3)(4,5)]
            sage: g = G.gen(0)
            sage: g(5)
            4
        """
        return int(gap.eval('%s^%s'%(i, self._gap_().name())))

    #def __mul__(self, other):
    #    """
    #    This overloaded operator implements multiplication *and*
    #    permutation action on polynomials.
    #    EXAMPLES:
    #        sage: G = PermutationGroup(['(1,2)(3,4)', '(3,4,5,6)'])
    #        sage: g = G.gens()
    #        sage: g[0] * g[1]
    #        (1,2)(3,5,6)
    #    """
    #    if isinstance(other, MPolynomial):
    #        return self.right_action_on_polynomial(other)
    #    else:
    #        return element.MultiplicativeGroupElement.__mul__(self, other)

    def _r_action(self, left):
        """
        Return the right action of self on left.

        For example, if f=left is a polynomial, then this function
        returns f(sigma*x), which is image of f under the right action
        of sigma on the indeterminates.  This is a right action since
        the image of f(sigma*x) under tau is f(sigma*tau*x).

        INPUT:
            left -- element of space on which permutations act from the right

        EXAMPLES:
            sage: G = PermutationGroup(['(1,2,3)(4,5)', '(1,2,3,4,5)'])
            sage: R.<x,y,z,u,v> = MPolynomialRing(QQ,5)
            sage: f = x^2 + y^2 - z^2 + 2*u^2
            sage: sigma, tau = G.gens()
            sage: f*sigma
            2*v^2 + z^2 + y^2 - x^2
            sage: f*tau
            2*v^2 - u^2 + z^2 + y^2
            sage: f*(sigma*tau)
            u^2 + z^2 - y^2 + 2*x^2
            sage: (f*sigma)*tau
            u^2 + z^2 - y^2 + 2*x^2
        """
        if is_Polynomial(left):
            if self != 1:
                raise ValueError, "%s does not act on %s"%(self, left.parent())
            return left
        elif isinstance(left, PermutationGroupElement):
            return PermutationGroupElement(self._gap_()*left._gap_(),
                                           parent = None, check = True)
        elif is_MPolynomial(left):
            F = left.base_ring()
            R = left.parent()
            x = R.gens()
            vars = list(x)
            try:
                sigma_x  = [vars[int(self(i+1)-1)] for i in range(len(x))]
            except IndexError:
                raise ValueError, "%s does not act on %s"%(self, left.parent())
            return left(tuple(sigma_x))
        else:
            raise TypeError, "left (=%s) must be a polynomial."%left


    def _mul_(self, other):
        return PermutationGroupElement(self._gap_()*other._gap_(),
                                       self.parent(), check = True)

    def _div_(self, other):
        """
        Returns self divided by other, i.e., self times the inverse
        of other.

        EXAMPLES:
            sage: g = PermutationGroupElement('(1,2,3)(4,5)')
            sage: h = PermutationGroupElement('(1,2,3)')
            sage: g/h
            (4,5)
        """
        return PermutationGroupElement(self._gap_()/other._gap_(),
                                           self.parent(),
                                           check = True)

    def __invert__(self):
        """
        Return the inverse of this permutation.

        EXAMPLES:
            sage: g = PermutationGroupElement('(1,2,3)(4,5)')
            sage: ~g
            (1,3,2)(4,5)
            sage: (~g) * g
            ()
        """
        return PermutationGroupElement(self._gap_().Inverse(),
                          self.parent(), check=False)

    def list(self):
        """
        Returns list of the images of the integers from 1 to n under
        this permutation as a list of Python ints.

        EXAMPLES:
            sage: G = SymmetricGroup(4)
            sage: x = G([2,1,4,3]); x
            (1,2)(3,4)
            sage: v = x.list(); v
            [2, 1, 4, 3]
            sage: type(v[0])
            <type 'int'>
        """
        return eval(gap.eval('ListPerm(%s)'%self.__gap))

    def order(self):
        """
        Return the order of this group element, which is the smallest
        positive integer $n$ for which $g^n = 1$.

        EXAMPLES:
            sage: s = PermutationGroupElement('(1,2)(3,5,6)')
            sage: s.order()
            6
        """
        return int(self._gap_().Order())

    def sign(self):
        """
        Returns the sign of self, which is $(-1)^{s}$, where $s$ is
        the number of swaps.

        EXAMPLES:
            sage: s = PermutationGroupElement('(1,2)(3,5,6)')
            sage: s.sign()
            -1
        """
        return int(self._gap_().SignPerm())


    def orbit(self, n):
        """
        Returns the orbit of the integer $n$ under this group element,
        as a sorted list of integers.

        EXAMPLES:
            sage: G = PermutationGroup(['(1,2,3)(4,5)'])
            sage: g = G.gen(0)
            sage: g.orbit(4)
            [4, 5]
            sage: g.orbit(3)
            [1, 2, 3]
            sage: g.orbit(10)
            [10]
        """
        n = Integer(n)
        # We use eval to avoid creating intermediate gap objects (so
        # this is slightly faster.
        s = gap.eval('Orbit(Group(%s), %s)'%(self._gap_().name(), Integer(n)))
        v = [Integer(k) for k in eval(s)]
        v.sort()
        return v

    def matrix(self):
        """
        Returns deg x deg permutation matrix associated
        to the permutation self

        EXAMPLES:
            sage: G = PermutationGroup(['(1,2,3)(4,5)'])
            sage: g = G.gen(0)
            sage: g.matrix()
            [0 1 0 0 0]
            [0 0 1 0 0]
            [1 0 0 0 0]
            [0 0 0 0 1]
            [0 0 0 1 0]
        """
        deg = self.parent().degree()
        M = MatrixSpace(ZZ,deg,deg)
        A = M(0)
        for i in range(deg):
            A[i, self(i+1) - 1] = 1
        return A

    def word_problem(g, words, display=True):
        """
        G and H are permutation groups, g in G, H is a subgroup of G
        generated by a list (words) of elements of G. If g is in H,
        return the expression for g as a word in the elements of
        (words).

        This function does not solve the word problem in SAGE. Rather
        it pushes it over to GAP, which has optimized algorithms for
        the word problem. Essentially, this function is a wrapper for the GAP
        functions "EpimorphismFromFreeGroup" and "PreImagesRepresentative".

        EXAMPLE:
            sage: G = PermutationGroup([[(1,2,3),(4,5)],[(3,4)]])
            sage: g1 = G.gens()[0]
            sage: g2 = G.gens()[1]
            sage: h = g1^2*g2*g1
            sage: h.word_problem([g1,g2], False)
            '(1,2,3)(4,5)^2*(3,4)^-1*(1,2,3)(4,5)'
            sage: h.word_problem([g1,g2])
               x1^2*x2^-1*x1
               [['(1,2,3)(4,5)', 2], ['(3,4)', -1], ['(1,2,3)(4,5)', 1]]
            '(1,2,3)(4,5)^2*(3,4)^-1*(1,2,3)(4,5)'
        """
        import copy
        from sage.groups.perm_gps.permgroup import PermutationGroup
        from sage.interfaces.all import gap
        G = g.parent()
        #print G
        gap.eval("l:=One(Rationals)")
        s1 = "gens := GeneratorsOfGroup(%s)"%G._gap_init_()
        #print s1
        gap.eval(s1)
        s2 = "g0:=%s; gensH:=%s"%(gap(g),words)
        gap.eval(s2)
        s3 = 'G:=Group(gens); H:=Group(gensH)'
        #print s3
        gap.eval(s3)
        phi = gap.eval("hom:=EpimorphismFromFreeGroup(G)")
        #print phi
        l1 = gap.eval("ans:=PreImagesRepresentative(hom,g0)")
        l2 = copy.copy(l1)
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
