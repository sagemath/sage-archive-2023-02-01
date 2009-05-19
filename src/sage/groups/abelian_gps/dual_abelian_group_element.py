"""
Elements (characters) of the dual group of a finite Abelian group.

AUTHORS:
    - David Joyner (2006-07); based on abelian_group_element.py.
    - David Joyner (2006-10); modifications suggested by William Stein.

EXAMPLES:
    sage: F = AbelianGroup(5,[2, 3, 5, 7, 8], names="abcde")
    sage: a,b,c,d,e = F.gens()
    sage: Fd = DualAbelianGroup(F, names = "ABCDE")
    sage: A,B,C,D,E = Fd.gens()
    sage: A*B^2*D^7
    A*B^2
    sage: A(a)    ## random last few digits
    -1.0000000000000000 + 0.00000000000000013834419720915037*I
    sage: B(b)
    -0.500000000000000 + 0.866025403784439*I
    sage: A(a*b)    ## random last few digits
    -1.0000000000000000 + 0.00000000000000013834419720915037*I
    sage: (A*B*C^2*D^20*E^65).list()
    [1, 1, 2, 6, 1]
    sage: B^(-1)
    B^2

It is important to note that lists are mutable and the
returned list is not a copy.  As a result, reassignment
of an element of the list changes the object.

    sage: X = A*B*C^2*D^2*E^-6
    sage: X.list()
    [1, 1, 2, 2, 2]
    sage: X.list()[1] = -1
    sage: X
    A*B^-1*C^2*D^2*E^2

"""

###########################################################################
#  Copyright (C) 2006 William Stein <wstein@gmail.com>
#  Copyright (C) 2006 David Joyner  <wdjoyner@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

import operator

from sage.rings.integer import Integer
from sage.structure.element import MonoidElement
from sage.rings.infinity import infinity
from sage.rings.arith import *
from sage.misc.misc import prod, add
from abelian_group_element import AbelianGroupElement,is_AbelianGroupElement
from abelian_group import AbelianGroup
from sage.misc.functional import exp
from sage.rings.complex_field import is_ComplexField


def add_strings(x, z=0):
    """
    This was in sage.misc.misc but commented out. Needed to add
    lists of strings in the word_problem method below.

    Return the sum of the elements of x.  If x is empty,
    return z.

    INPUT:
        x -- iterable
        z -- the "0" that will be returned if x is empty.

    OUTPUT:
        object
    """
    if len(x) == 0:
        return z
    if not isinstance(x, list):
        m = x.__iter__()
        y = m.next()
        return reduce(operator.add, m, y)
    else:
        return reduce(operator.add, x[1:], x[0])


def is_DualAbelianGroupElement(x):
    return isinstance(x, DualAbelianGroupElement)

class DualAbelianGroupElement(MonoidElement):
    def __init__(self, F, X):
        """
        Create an element X of the DualAbelianGroup of F.

        EXAMPLES:
            sage: F = AbelianGroup(3,[7,8,9])
            sage: Fd = DualAbelianGroup(F,names="ABC")
            sage: A,B,C = Fd.gens()
            sage: A*B^-1 in Fd
            True

        """
        MonoidElement.__init__(self, F)
        self.__repr = None
        G = F.group()
        n = G.ngens()
        if isinstance(X, (int, Integer)) and X == 1:
            self.__element_vector = [ 0 for i in range(n) ]
        elif isinstance(X, list):
            if len(X) != n:
                raise IndexError, \
                      "Argument length (= %s) must be %s."%(len(X), n)
            self.__element_vector = X
        else:
            raise TypeError, "Argument X (= %s) is of wrong type."%X

    def list(self):
        """
        Return (a reference to) the underlying list used to represent
        this element.  If this is a word in an abelian group on $n$
        generators, then this is a list of nonnegative integers of
        length $n$.

        EXAMPLES:
            sage: F = AbelianGroup(5,[2, 3, 5, 7, 8], names="abcde")
            sage: a,b,c,d,e = F.gens()
            sage: Ad = DualAbelianGroup(F, names = "ABCDE")
            sage: A,B,C,D,E = Ad.gens()
            sage: (A*B*C^2*D^20*E^65).list()
            [1, 1, 2, 6, 1]
            sage: X = A*B*C^2*D^2*E^-6
            sage: X.list()
            [1, 1, 2, 2, 2]
            sage: X.list()[1] = -1
            sage: X
            A*B^-1*C^2*D^2*E^2
        """
        return self.__element_vector

    def _repr_(self):
        s = ""
        A = self.parent()
        n = A.ngens()
        x = A.variable_names()
        v = self.list()
        for i in range(n):
            if v[i] == 0:
                continue
            elif v[i] == 1:
                if len(s) > 0: s += "*"
                s += "%s"%x[i]
            else:
                if len(s) > 0: s += "*"
                s += "%s^%s"%(x[i],v[i])
        if len(s) == 0: s = "1"
        return s

    def __mul__(self, y):
        #Same as _mul_ in AbelianGroupElement

        M = self.parent()
        n = M.ngens()
        invs = M.invariants()
        z = M(1)
        xelt = self.list()
        yelt = y.list()
        zelt = [ xelt[i]+yelt[i] for i in range(len(xelt)) ]
        if len(invs) >= n:
            L =  []
            for i in range(len(xelt)):
                if invs[i]!=0:
                    L.append(zelt[i]%invs[i])
                if invs[i]==0:
                    L.append(zelt[i])
            z.__element_vector = L
            #print z.__element_vector
        if len(invs) < n:
            L1 =  []
            for i in range(len(invs)):
                if invs[i]!=0:
                    L1.append(zelt[i]%invs[i])
                if invs[i]==0:
                    L1.append(zelt[i])
            L2 =  [ zelt[i] for i in range(len(invs),len(xelt)) ]
            z.__element_vector = L1+L2
        return M(z.__element_vector)

    def __pow__(self, n):
        """
        requires that len(invs) = n
        """
        if not isinstance(n, (int, long, Integer)):
            raise TypeError, "Argument n (= %s) must be an integer."%n
        n = int(n)
        M = self.parent()
        N = M.ngens()
        invs = M.invariants()
        if n < 0:
            L =[n*self.list()[i]%M.gen(i).order() for i in range(M.ngens())]
            return prod([M.gen(i)**L[i] for i in range(M.ngens())])
            #m = LCM(invs) ## Not very efficient version
            #pw = (n)%m
            #x = self**pw
            #return x
        elif n == 0:
            return M(1)
        elif n == 1:
            return self
        elif n == 2:
            return self * self
        k = n//2
        return self**k * self**(n-k)

    def __cmp__(self,other):
        if (self.list() != other.list()):
                return -1
        return 0

    def order(self):
        """
        Returns the (finite) order of this element.

        EXAMPLES:
            sage: F = AbelianGroup(3,[7,8,9])
            sage: Fd = DualAbelianGroup(F)
            sage: A,B,C = Fd.gens()
            sage: (B*C).order()
            72
        """
        M = self.parent()
        #print self, M
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

    def __call__(self,g):
        """
        Computes the value of a character self on a group element
        g (g must be an element of self.group())

        EXAMPLES:
            sage: F = AbelianGroup(5, [2,3,5,7,8], names="abcde")
            sage: a,b,c,d,e = F.gens()
            sage: Fd = DualAbelianGroup(F, names="ABCDE")
            sage: A,B,C,D,E = Fd.gens()
            sage: A*B^2*D^7
            A*B^2
            sage: A(a)    ## random last few digits
            -1.0000000000000000 + 0.00000000000000013834419720915037*I
            sage: B(b)    ## random last few digits
            -0.49999999999999983 + 0.86602540378443871*I
            sage: A(a*b)    ## random last few digits
            -1.0000000000000000 + 0.00000000000000013834419720915037*I
        """
        F = self.parent().base_ring()
        expsX = list(self.list())
        expsg = list(g.list())
        invs = self.parent().invariants()
        N = LCM(invs)
        if is_ComplexField(F):
            from sage.symbolic.constants import pi
            I = F.gen()
            PI = F(pi)
            ans = prod([exp(2*PI*I*expsX[i]*expsg[i]/invs[i]) for i in range(len(expsX))])
            return ans
        ans = F(1)  ## assumes F is the cyclotomic field
        zeta = F.gen()
        #print F,zeta
        for i in range(len(expsX)):
            inv_noti = N/invs[i]
            ans = ans*zeta**(expsX[i]*expsg[i]*inv_noti)
        return ans

    def word_problem(self, words, display=True):
        """
        This is a rather hackish method and is included for completeness.

        The word problem for an instance of DualAbelianGroup as it can
        for an AbelianGroup. The reason why is that word problem
        for an instance of AbelianGroup simply calls GAP (which
        has abelian groups implemented) and invokes "EpimorphismFromFreeGroup"
        and "PreImagesRepresentative". GAP does not have duals of
        abelian groups implemented. So, by using the same name
        for the generators, the method below converts the problem for
        the dual group to the corresponding problem on the group
        itself and uses GAP to solve that.

        EXAMPLES:
            sage: G = AbelianGroup(5,[3, 5, 5, 7, 8],names="abcde")
            sage: Gd = DualAbelianGroup(G,names="abcde")
            sage: a,b,c,d,e = Gd.gens()
            sage: u = a^3*b*c*d^2*e^5
            sage: v = a^2*b*c^2*d^3*e^3
            sage: w = a^7*b^3*c^5*d^4*e^4
            sage: x = a^3*b^2*c^2*d^3*e^5
            sage: y = a^2*b^4*c^2*d^4*e^5
            sage: e.word_problem([u,v,w,x,y],display=False)
            [[b^2*c^2*d^3*e^5, 245]]

        The command e.word_problem([u,v,w,x,y],display=True) returns
        the same list but also prints $e = (b^2*c^2*d^3*e^5)^245$.

        """
        ## First convert the problem to one using AbelianGroups
        import copy
        from sage.groups.abelian_gps.abelian_group import AbelianGroup
        from sage.interfaces.all import gap
        M = self.parent()
        G = M.group()
        gens = M.variable_names()
        g = prod([G.gen(i)**(self.list()[i]) for i in range(G.ngens())])
        gap.eval("l:=One(Rationals)")            ## trick needed for LL line below to keep SAGE from parsing
        s1 = "gens := GeneratorsOfGroup(%s)"%G._gap_init_()
        gap.eval(s1)
        for i in range(len(gens)):
           cmd = ("%s := gens["+str(i+1)+"]")%gens[i]
           gap.eval(cmd)
        s2 = "g0:=%s; gensH:=%s"%(str(g),words)
        gap.eval(s2)
        s3 = 'G:=Group(gens); H:=Group(gensH)'
        gap.eval(s3)
        phi = gap.eval("hom:=EpimorphismFromFreeGroup(H)")
        l1 = gap.eval("ans:=PreImagesRepresentative(hom,g0)")
        l2 = copy.copy(l1)
        l4 = []
        l3 = l1.split("*")
        for i in range(1,len(words)+1):
            l2 = l2.replace("x"+str(i),"("+str(words[i-1])+")")
        l3 = eval(gap.eval("L3:=ExtRepOfObj(ans)"))
        nn = eval(gap.eval("n:=Int(Length(L3)/2)"))
        LL1 = eval(gap.eval("L4:=List([l..n],i->L3[2*i])"))         ## note the l not 1
        LL2 = eval(gap.eval("L5:=List([l..n],i->L3[2*i-1])"))       ## note the l not 1
        if display:
            s = str(g)+" = "+add_strings(["("+str(words[LL2[i]-1])+")^"+str(LL1[i])+"*" for i in range(nn)])
            m = len(s)
            print "      ",s[:m-1],"\n"
        return [[words[LL2[i]-1],LL1[i]] for i in range(nn)]




