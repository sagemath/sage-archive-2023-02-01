"""nodoctest
Abelian group elements

TODO:
    * remove nodoctest above to doctest this file

AUTHORS:
    - David Joyner (2006-02); based on free_abelian_monoid_element.py, written by David Kohel.

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
#  Copyright (C) 2006 William Stein <wstein@ucsd.edu>
#  Copyright (C) 2006 David Joyner
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

import operator

from sage.ext.integer import Integer
from sage.ext.element import MonoidElement
from sage.rings.arith import *

def is_AbelianGroupElement(x):
    return isinstance(x, AbelianGroupElement)

class AbelianGroupElement(MonoidElement):
    def __init__(self, F, x):
        """
        Create the element x of the AbelianGroup F.

        EXAMPLES:
            sage: F = AbelianGroup(5, [3,4,5,8,7], 'abcde')
            sage: a, b, c, d, e = F.gens()
            sage: a^2 * b^3 * a^2 * b^-4
            a*b^3
            sage: b^-11
            b^3
            sage: a^-11
            a^2
            sage: a*b in F
            True

        """
        MonoidElement.__init__(self, F)
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
        if len(s) == 0: s = "1"
        return s

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
            z.__element_vector =  [ zelt[i]%invs[i] for i in range(len(xelt)) ]
        if len(invs) < n:
            L1 =  [ zelt[i]%invs[i] for i in range(len(invs)) ]
            L2 =  [ zelt[i] for i in range(len(invs),len(xelt)) ]
            z.__element_vector = L1+L2
        return z

    def __pow__(self, n):
        """
        requires that len(invs) = n
        """
        if not isinstance(n, (int, long, Integer)):
            raise TypeError, "Argument n (= %s) must be an integer."%n
        M = self.parent()
        N = M.ngens()
        invs = M.invariants()
        if n < 0:
            m = LCM(invs)
            pw = (-n)%m
            x = self**pw
            return x
        elif n == 0:
            return self.parent()(1)
        elif n == 1:
            return self
        elif n == 2:
            return self * self
        k = n//2
        return self**k * self**(n-k)

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
            Abelian group on 3 generators (f_0, f_1, f_2) with invariants [7, 8, 9]
            sage: F.gens()[2].order()
            9
            sage: a,b,c = F.gens()
            sage: (b*c).order()
            72
        """
        M = self.parent()
        if self == M(1):
            return 1
        invs = M.invariants()
        if self in M.gens():
            return invs[list(M.gens()).index(self)]
        L = list(self.list())
        N = LCM([invs[i]*L[i] for i in range(len(invs)) if L[i]!=0])
        return N
