"""
Abelian group elements

AUTHORS:
    - David Joyner (2006-02); based on free_abelian_monoid_element.py, written by David Kohel.
    - David Joyner (2006-05); bug fix in order
    -              (2006-08); bug fix+new method in pow for negatives+fixed corresponding examples.

EXAMPLES:
Recall an example from abelian groups.
    sage: F = AbelianGroup(5,[4,5,5,7,8],names = list("abcde"))
    sage: (a,b,c,d,e) = F.gens()
    sage: x = a*b^2*e*d^20*e^12
    sage: x
    a*b^2*d^6*e^5
    sage: x = a^10*b^12*c^13*d^20*e^12
    sage: x
    a^2*b^2*c^3*d^6*e^4
    sage: y = a^13*b^19*c^23*d^27*e^72
    sage: y
    a*b^4*c^3*d^6
    sage: x*y
    a^3*b*c*d^5*e^4
    sage: x.list()
    [2, 2, 3, 6, 4]

It is important to note that lists are mutable and the
returned list is not a copy.  As a result, reassignment
of an element of the list changes the object.
    sage: x.list()[0] = 3
    sage: x.list()
    [3, 2, 3, 6, 4]
    sage: x
    a^3*b^2*c^3*d^6*e^4

"""

###########################################################################
#  Copyright (C) 2006 William Stein <wstein@gmail.com>
#  Copyright (C) 2006 David Joyner
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

import operator

from sage.rings.integer import Integer
from sage.structure.element import MultiplicativeGroupElement
from sage.rings.infinity import infinity
from sage.misc.misc import prod
from sage.rings.arith import LCM, GCD


def is_AbelianGroupElement(x):
    """
    Return true if x is an abelian group element, i.e., an element of type
    AbelianGroupElement.

    EXAMPLES:
    Though the integer 3 is in the integers, and the integers have an abelian
    group structure, 3 is not an AbelianGroupElement:
         sage: from sage.groups.abelian_gps.abelian_group_element import is_AbelianGroupElement
         sage: is_AbelianGroupElement(3)
         False
         sage: F = AbelianGroup(5, [3,4,5,8,7], 'abcde')
         sage: is_AbelianGroupElement(F.0)
         True
    """
    return isinstance(x, AbelianGroupElement)

class AbelianGroupElement(MultiplicativeGroupElement):
    def __init__(self, F, x):
        """
        Create the element x of the AbelianGroup F.

        EXAMPLES:
            sage: F = AbelianGroup(5, [3,4,5,8,7], 'abcde')
            sage: a, b, c, d, e = F.gens()
            sage: a^2 * b^3 * a^2 * b^-4
            a*b^3
            sage: b^-11
            b
            sage: a^-11
            a
            sage: a*b in F
            True

        """
        MultiplicativeGroupElement.__init__(self, F)
        self.__repr = None
        n = F.ngens()
        if isinstance(x, (int, Integer)) and x == 1:
            self.__element_vector = [ 0 for i in range(n) ]
        elif isinstance(x, list):
            if len(x) != n:
                raise IndexError, \
                      "Argument length (= %s) must be %s."%(len(x), n)
            self.__element_vector = x
        else:
            raise TypeError, "Argument x (= %s) is of wrong type."%x

    def _repr_(self):
        """
        EXAMPLES:
            sage: AbelianGroupElement(AbelianGroup(1, [1], names='e'),[0])
            1
            sage: AbelianGroupElement(AbelianGroup(3, [2,3,4], names='e'),[1,2,3])
            e0*e1^2*e2^3

        """
        s = ""
        A = self.parent()
        n = A.ngens()
        x = A.variable_names()
        v = self.__element_vector
        for i in range(n):
            if v[i] == 0:
                continue
            elif v[i] == 1:
                if len(s) > 0: s += "*"
                s += "%s"%x[i]
            else:
                if len(s) > 0: s += "*"
                s += "%s^%s"%(x[i],v[i])
        if len(s) == 0: s = str(1)
        return s

    def _div_(self, y):
        """
        TESTS:
            sage: G.<a,b> = AbelianGroup(2)
            sage: a/b
            a*b^-1

        """
        return self*y.inverse()

    def _mul_(self, y):
        #Same as _mul_ in FreeAbelianMonoidElement except that the
        #exponents get reduced mod the invariant.

        M = self.parent()
        n = M.ngens()
        invs = M.invariants()
        z = M(1)
        xelt = self.__element_vector
        yelt = y.__element_vector
        zelt = [ xelt[i]+yelt[i] for i in range(len(xelt)) ]
        if len(invs) >= n:
            L =  []
            for i in range(len(xelt)):
                if invs[i]!=0:
                    L.append(zelt[i]%invs[i])
                if invs[i]==0:
                    L.append(zelt[i])
            z.__element_vector = L
        if len(invs) < n:
            L1 =  []
            for i in range(len(invs)):
                if invs[i]!=0:
                    L1.append(zelt[i]%invs[i])
                if invs[i]==0:
                    L1.append(zelt[i])
            L2 =  [ zelt[i] for i in range(len(invs),len(xelt)) ]
            z.__element_vector = L1+L2
        return z

    def __pow__(self, _n):
        """
        requires that len(invs) = n
        """
        n = int(_n)
        if n != _n:
            raise TypeError, "Argument n (= %s) must be an integer."%n
        M = self.parent()
        N = M.ngens()
        invs = M.invariants()
        z = M(1)
        for i in xrange(len(invs)):
            if invs[i] == 0:
                z.__element_vector[i] = self.__element_vector[i]*n
            else:
                z.__element_vector[i] = (self.__element_vector[i]*n)%invs[i]
        return z

    def inverse(self):
        """
        Returns the inverse element.

        EXAMPLE:
            sage: G.<a,b> = AbelianGroup(2)
            sage: a^-1
            a^-1
            sage: a.list()
            [1, 0]
            sage: a.inverse().list()
            [-1, 0]
            sage: a.inverse()
            a^-1

        """
        new_list = [-1*n for n in self.list()]
        par_inv = self.parent().invariants()
        for i in xrange(len(par_inv)):
            if par_inv[i] != 0 and new_list[i] < 0:
                while new_list[i] < 0:
                    new_list[i] += abs(par_inv[i])
        return AbelianGroupElement(self.parent(), new_list)

    def as_permutation(self):
        r"""
        Return the element of the permutation group G (isomorphic to the
        abelian group A) associated to a in A.

        EXAMPLES:
            sage: G = AbelianGroup(3,[2,3,4],names="abc"); G
            Multiplicative Abelian Group isomorphic to C2 x C3 x C4
            sage: a,b,c=G.gens()
            sage: Gp = G.permutation_group(); Gp
            Permutation Group with generators [(1,13)(2,14)(3,15)(4,16)(5,17)(6,18)(7,19)(8,20)(9,21)(10,22)(11,23)(12,24), (1,5,9)(2,6,10)(3,7,11)(4,8,12)(13,17,21)(14,18,22)(15,19,23)(16,20,24), (1,3,2,4)(5,7,6,8)(9,11,10,12)(13,15,14,16)(17,19,18,20)(21,23,22,24)]
            sage: a.as_permutation()
            (1,13)(2,14)(3,15)(4,16)(5,17)(6,18)(7,19)(8,20)(9,21)(10,22)(11,23)(12,24)
            sage: ap = a.as_permutation(); ap
            (1,13)(2,14)(3,15)(4,16)(5,17)(6,18)(7,19)(8,20)(9,21)(10,22)(11,23)(12,24)
            sage: ap in Gp
            True
        """
        from sage.groups.perm_gps.permgroup import PermutationGroup
        from sage.interfaces.all import gap
        G = self.parent()
        invs = G.invariants()
        s1 = 'A:=AbelianGroup(%s)'%invs
        gap.eval(s1)
        s2 = 'phi:=IsomorphismPermGroup(A)'
        gap.eval(s2)
        s3 = "gens := GeneratorsOfGroup(A)"
        gap.eval(s3)
        L = self.list()
        gap.eval("L1:="+str(L))
        s4 = "L2:=List([1..%s], i->gens[i]^L1[i]);"%len(L)
        gap.eval(s4)
        pg = gap.eval("Image(phi,Product(L2))")
        Gp = G.permutation_group()
        gp = Gp(pg)
        return gp

    def list(self):
        """
        Return (a reference to) the underlying list used to represent
        this element.  If this is a word in an abelian group on $n$
        generators, then this is a list of nonnegative integers of
        length $n$.

        EXAMPLES:
            sage: F = AbelianGroup(5, [3,4,5,8,7], 'abcde')
            sage: (a, b, c, d, e) = F.gens()
            sage: a.list()
            [1, 0, 0, 0, 0]
        """
        return self.__element_vector

    def __cmp__(self,other):
        if (self.list() != other.list()):
                return -1
        return 0

    def order(self):
        """
        Returns the (finite) order of this element or Infinity if this element
        does not have finite order.

        EXAMPLES:
            sage: F = AbelianGroup(3,[7,8,9]); F
            Multiplicative Abelian Group isomorphic to C7 x C8 x C9
            sage: F.gens()[2].order()
            9
            sage: a,b,c = F.gens()
            sage: (b*c).order()
            72
        """
        M = self.parent()
        if self == M(1):
            return Integer(1)
        invs = M.invariants()
        if self in M.gens():
            o = invs[list(M.gens()).index(self)]
            if o == 0:
                return infinity
            return o
        L = list(self.list())
        N = LCM([invs[i]/GCD(invs[i],L[i]) for i in range(len(invs)) if L[i]!=0])   ####### error here????
        if N == 0:
            return infinity
        else:
            return N

    def random_element(self):
        """
        Return a random element of this dual group.
        """
        if not(self.is_finite()):
            raise NotImplementedError, "Only implemented for finite groups"
        gens = self.gens()
        g = gens[0]**0
        for i in range(len(gens)):
            g = g*gens[i]**(random(gens[i].order()))
        return g



    def word_problem(self, words):
        """
        TODO -- this needs a rewrite -- see stuff in the matrix_grp directory.

        G and H are abelian groups, g in G, H is a subgroup of G generated by a
        list (words) of elements of G. If self is in H, return the expression
        for self as a word in the elements of (words).

        This function does not solve the word problem in SAGE. Rather
        it pushes it over to GAP, which has optimized algorithms for
        the word problem.

        WANRING: Don't use E (or other GAP-reserved letters) as a generator
        name.

        EXAMPLE:
            sage: G = AbelianGroup(2,[2,3], names="xy")
            sage: x,y = G.gens()
            sage: x.word_problem([x,y])
            [[x, 1]]
            sage: y.word_problem([x,y])
            [[y, 1]]
            sage: (y*x).word_problem([x,y])
            [[x, 1], [y, 1]]

        """
        from sage.groups.abelian_gps.abelian_group import AbelianGroup, word_problem
        return word_problem(words,self)
