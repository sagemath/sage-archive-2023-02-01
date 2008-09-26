r"""
Basic functionality for dual groups of finite multiplicative Abelian groups

AUTHOR:
    -- David Joyner (2006-08) (based on abelian_groups)
    -- David Joyner (2006-10) modifications suggested by William Stein

TODO:
   * additive abelian groups should also be supported.

The basic idea is very simple. Let G be an abelian group and
$G^*$ its dual (i.e., the group of homomorphisms from G to
$\CC^\times$). Let $g_j$, $j=1,..,n$, denote generators of $G$ - say $g_j$
is of order $m_j>1$. There are generators $X_j$, $j=1,..,n$, of $G^*$ for which
$X_j(g_j)=\exp(2\pi i/m_j)$ and $X_i(g_j)=1$ if $i\not= j$. These
are used to construct $G^*$ in the DualAbelianGroup class below.

SAGE supports multiplicative abelian groups on any prescribed finite
number $n > 0$ of generators.  Use the \code{AbelianGroup} function
to create an abelian group, the \code{DualAbelianGroup} function
to create its dual, and then the \code{gen} and \code{gens}
functions to obtain the corresponding generators.  You can print the
generators as arbitrary strings using the optional \code{names}
argument to the \code{DualAbelianGroup} function.


"""

##########################################################################
#  Copyright (C) 2006 David Joyner <wdjoyner@gmail.com> and William Stein
#
#  Distributed under the terms of the GNU General Public License (GPL):
#
#                  http://www.gnu.org/licenses/
##########################################################################

import weakref
import copy

from sage.rings.integer import Integer

from sage.rings.infinity import infinity
from sage.rings.arith import factor,is_prime_power,LCM
from abelian_group_element import AbelianGroupElement,is_AbelianGroupElement
from sage.misc.misc import add, prod
import sage.groups.group as group
from abelian_group import AbelianGroup
from dual_abelian_group_element import DualAbelianGroupElement,is_DualAbelianGroupElement
from sage.misc.mrange import mrange
from sage.rings.integer_ring import IntegerRing
ZZ = IntegerRing()
from sage.rings.complex_field import ComplexField
CC = ComplexField()
from sage.rings.number_field.number_field import CyclotomicField

def DualAbelianGroup(G, names="X", base_ring=CC):
    r"""
    Create the dual group of the multiplicative abelian group $G$.

    INPUT:
        G     -- a finite abelian group
        names -- (optional) names of generators

    OUTPUT:
        The dual group of G.

    EXAMPLES:
        sage: F = AbelianGroup(5, [2,5,7,8,9], names='abcde')
        sage: (a, b, c, d, e) = F.gens()
        sage: Fd = DualAbelianGroup(F,names='ABCDE')
        sage: A,B,C,D,E = Fd.gens()
        sage: A(a)    ## random
        -1.0000000000000000 + 0.00000000000000013834419720915037*I
        sage: A(b); A(c); A(d); A(e)
        1.00000000000000
        1.00000000000000
        1.00000000000000
        1.00000000000000
    """
    if G.order() is infinity:
        NotImplementedError, "The group must be finite"
    #infac = G.invariants()
    #n = G.ngens()
    #namesG = [G.gen(i) for i in range(n)]
    M = DualAbelianGroup_class(G, names, base_ring)
    return M

def is_DualAbelianGroup(x):
    """
    Return True if $x$ is the dual group of an abelian group.

    EXAMPLES:
        sage: from sage.groups.abelian_gps.dual_abelian_group import is_DualAbelianGroup
        sage: F = AbelianGroup(5,[3,5,7,8,9],names = list("abcde"))
        sage: Fd = DualAbelianGroup(F)
        sage: is_DualAbelianGroup(Fd)
        True
        sage: F = AbelianGroup(3,[1,2,3],names='a')
        sage: Fd = DualAbelianGroup(F)
        sage: Fd.gens()
        (X0, X1)
        sage: F.gens()
        (a0, a1)

    """
    return isinstance(x, DualAbelianGroup_class)


class DualAbelianGroup_class(group.AbelianGroup):
    """
    Dual of abelian group.

    EXAMPLES:
        sage: F = AbelianGroup(5,[3,5,7,8,9],names = list("abcde"))
        sage: DualAbelianGroup(F)
        Dual of Abelian Group isomorphic to Z/3Z x Z/5Z x Z/7Z x Z/8Z x Z/9Z  over Complex Field with 53 bits of precision
        sage: F = AbelianGroup(4,[15,7,8,9],names = list("abcd"))
        sage: DualAbelianGroup(F)
        Dual of Abelian Group isomorphic to Z/3Z x Z/5Z x Z/7Z x Z/8Z x Z/9Z  over Complex Field with 53 bits of precision

    """
    def __init__(self, G, names="X", bse_ring=None):
        """
        If G has invariants invs = [n1,...,nk] then
        the default base_ring is CyclotoomicField(N), where
        N = LCM(n1,...,nk).
        """
        if bse_ring == None:
            base_ring = CyclotomicField(LCM(G.invariants()))
        else:
            base_ring = bse_ring
        self.__group = G
        self._assign_names(names)
        self._base_ring = base_ring

    def group(self):
        return self.__group

    def base_ring(self):
        return self._base_ring

    def __str__(self):
        """
        Print method.

        EXAMPLES:
            sage: F = AbelianGroup(3,[5,64,729],names = list("abc"))
            sage: Fd = DualAbelianGroup(F)
	    sage: print Fd
            DualAbelianGroup( AbelianGroup ( 3, [5, 64, 729] ) )

        """
        s = "DualAbelianGroup( AbelianGroup ( %s, %s ) )"%(self.ngens(), self.invariants())
        return s

    def _repr_(self):
        """
        EXAMPLES:
            sage: F = AbelianGroup(5, [2,5,7,8,9], names='abcde')
            sage: Fd = DualAbelianGroup(F,names='ABCDE',base_ring = CyclotomicField(2*5*7*8*9))
            sage: Fd
            Dual of Abelian Group isomorphic to Z/2Z x Z/5Z x Z/7Z x Z/8Z x Z/9Z  over Cyclotomic Field of order 5040 and degree 1152
            sage: Fd = DualAbelianGroup(F,names='ABCDE')
            sage: Fd
            Dual of Abelian Group isomorphic to Z/2Z x Z/5Z x Z/7Z x Z/8Z x Z/9Z  over Complex Field with 53 bits of precision
        """
        G = self.group()
        eldv = G.elementary_divisors()
        gp = ""
        for x in eldv:
            if x!=0:
                gp = gp + "Z/%sZ x "%x
            if x==0:
                gp = gp + "Z x "
        gp = gp[:-2]
        s = "Dual of Abelian Group isomorphic to "+gp+" over %s"%self.base_ring()
        return s

    def _latex_(self):
        r"""
        Return latex representation of this group.

        EXAMPLES:
            sage: F = AbelianGroup(3, [2]*3)
            sage: Fd = DualAbelianGroup(F)
            sage: Fd._latex_()
            '${\rm DualAbelianGroup}( AbelianGroup ( 3, [2, 2, 2] ) )$'

        """
        s = "${\rm DualAbelianGroup}( AbelianGroup ( %s, %s ) )$"%(self.ngens(), self.invariants())
        return s

    def __call__(self, x):
        """
        Create an element of this abelian group from $x$.

        EXAMPLES:
            sage: F = AbelianGroup(10, [2]*10)
            sage: Fd = DualAbelianGroup(F)
            sage: Fd(Fd.2)
            X2
            sage: Fd(1)
            1
        """
        if isinstance(x, DualAbelianGroupElement) and x.parent() is self:
            return x
        return DualAbelianGroupElement(self, x)

    def random_element(self):
        """
        Return a random element of this dual group.

        EXAMPLES:
            sage: G = AbelianGroup([2,3,9])
            sage: Gd = DualAbelianGroup(G)
            sage: Gd.random_element()
            X0*X1^2*X2
            sage: N = 43^2-1
            sage: G = AbelianGroup([N],names="a")
            sage: Gd = DualAbelianGroup(G,names="A")
            sage: a, = G.gens()
            sage: A, = Gd.gens()
            sage: x = a^(N/4); y = a^(N/3); z = a^(N/14)
            sage: X = Gd.random_element(); X
            A^615
            sage: len([a for a in [x,y,z] if abs(X(a)-1)>10^(-8)])
            2

        """
        from sage.misc.prandom import randint
        gens = self.gens()
        g = gens[0]**0
        for i in range(len(gens)):
            g = g*gens[i]**(randint(1,gens[i].order()))
        return g

    def random(self):
        """
        Deprecated.  Use self.random_element() instead.
        """
        raise NotImplementedError, "Deprecated: use random_element() instead"


    def gen(self, i=0):
        """
        The $i$-th generator of the abelian group.

        EXAMPLES:
            sage: F = AbelianGroup(3,[1,2,3],names='a')
            sage: Fd = DualAbelianGroup(F, names="A")
            sage: Fd.0
            A0
            sage: Fd.1
            A1
            sage: Fd.invariants()
            [2, 3]
        """
        n = self.group().ngens()
        if i < 0 or i >= n:
            raise IndexError, "Argument i (= %s) must be between 0 and %s."%(i, n-1)
        x = [0]*int(n)
        x[int(i)] = 1
        return DualAbelianGroupElement(self, x)

    #def gens(self):
    #    """
    #    Return the generators for this subgroup.
    #
    #    """
    #    n = self.group().ngens()
    #    return tuple(self.gen(i) for i in range(n))

    def ngens(self):
        """
        The number of generators of the dual group.

        EXAMPLES:
            sage: F = AbelianGroup(1000)
            sage: Fd = DualAbelianGroup(F)
            sage: Fd.ngens()
            1000

        This can be slow for 10000 or more generators.
        """
        return self.group().ngens()

    def invariants(self):
        """
        The invariants of the dual group.

        EXAMPLES:
            sage: F = AbelianGroup(1000)
            sage: Fd = DualAbelianGroup(F)
            sage: invs = Fd.invariants(); len(invs)
            1000

        This can be slow for 10000 or more generators.
        """
        return self.group().invariants()

    def __contains__(self,X):
        """
        Implements "in".

        EXAMPLES:
            sage: F = AbelianGroup(5,[2, 3, 5, 7, 8], names="abcde")
            sage: a,b,c,d,e = F.gens()
            sage: Fd = DualAbelianGroup(F, names = "ABCDE")
            sage: A,B,C,D,E = Fd.gens()
            sage: A*B^2*D^7 in Fd
            True

        """
        return X.parent() == self and is_DualAbelianGroupElement(X)

    def order(self):
        """
        Return the order of this group.

        EXAMPLES:
            sage: G = AbelianGroup([2,3,9])
            sage: Gd = DualAbelianGroup(G)
            sage: Gd.order()
            54

        """
        G = self.group()
        return G.order()

    def is_commutative(self):
        """
        Return True since this group is commutative.

        EXAMPLES:
            sage: G = AbelianGroup([2,3,9])
            sage: Gd = DualAbelianGroup(G)
            sage: Gd.is_commutative()
            True
            sage: Gd.is_abelian()
            True

        """
        return True

    def list(self):
        """
        Return list of all elements of this group.

        EXAMPLES:
            sage: G = AbelianGroup([2,3], names = "ab")
            sage: Gd = DualAbelianGroup(G, names = "AB")
            sage: Gd.list()
            [1, B, B^2, A, A*B, A*B^2]

        """
        try:
            return list(self.__list)
        except AttributeError:
            pass
        if not(self.is_finite()):
           raise NotImplementedError, "Group must be finite"
        invs = self.invariants()
        T = mrange(invs)
        n = self.order()
        L = [DualAbelianGroupElement(self, t) for t in T]
        self.__list = L
        return list(self.__list)

    def __iter__(self):
        """
        Return an iterator over the elements of this group.

        EXAMPLES:
            sage: G = AbelianGroup([2,3], names = "ab")
            sage: Gd = DualAbelianGroup(G, names = "AB")
            sage: [X for X in Gd]
            [1, B, B^2, A, A*B, A*B^2]
            sage: N = 43^2-1
            sage: G = AbelianGroup([N],names="a")
            sage: Gd = DualAbelianGroup(G,names="A")
            sage: a, = G.gens()
            sage: A, = Gd.gens()
            sage: x = a^(N/4)
            sage: y = a^(N/3)
            sage: z = a^(N/14)
            sage: len([X for X in Gd if abs(X(x)-1)>0.01 and abs(X(y)-1)>0.01 and abs(X(z)-1)>0.01])
            880

        """
        for g in self.list():
            yield g



