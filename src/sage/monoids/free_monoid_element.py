"""
Monoid Elements

AUTHOR: David Kohel <kohel@maths.usyd.edu.au>, 2005/09/29

Elements of free monoids are represented internally as lists of pairs
of integers.
"""

#*****************************************************************************
#  Copyright (C) 2005 David Kohel <kohel@maths.usyd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty
#    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#  See the GNU General Public License for more details; the full text
#  is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import operator
from sage.ext.integer import Integer
from sage.ext.element import MonoidElement

def is_FreeMonoidElement(x):
    return isinstance(x, FreeMonoidElement)

class FreeMonoidElement(MonoidElement):
    """
    Element of a free monoid.
    """
    def __init__(self, F, x, check=True):
        """
        Create the element $x$ of the FreeMonoid $F$.

        This should typically be called by a FreeMonoid.
        """
        MonoidElement.__init__(self, F)
        if isinstance(x, (int, long, Integer)):
            if x == 1:
                self._element_list = []
            else:
                raise TypeError, "Argument x (= %s) is of the wrong type."%x
        elif isinstance(x, list):
            if check:
                x2 = []
                for v in x:
                    if not isinstance(v, tuple) and len(v) == 2:
                        raise TypeError, "x (= %s) must be a list of 2-tuples or 1."%x
                    if not (isinstance(v[0], (int,long,Integer)) and \
                            isinstance(v[1], (int,long,Integer))):
                        raise TypeError, "x (= %s) must be a list of 2-tuples of integers or 1."%x
                    if len(x2) > 0 and v[0] == x2[len(x2)-1][0]:
                        x2[len(x2)-1] = (v[0], v[1]+x2[len(x2)-1][1])
                    else:
                        x2.append(v)
                self._element_list = x2
            else:
                self._element_list = list(x)  # make copy, so user can't accidently change monoid.

        else:
            # TODO: should have some other checks here...
            raise TypeError, "Argument x (= %s) is of the wrong type."%x

    def __repr__(self):
        s = ""
        v = self._element_list
        x = self.parent().variable_names()
        for i in range(len(v)):
            if len(s) > 0: s += "*"
            g = x[int(v[i][0])]
            e = v[i][1]
            if e == 1:
                s += "%s"%g
            else:
                s += "%s^%s"%(g,e)
        if len(s) == 0: s = "1"
        return s

    def _latex_(self):
        """
        Return latex representation of self.

        EXAMPLES:
            sage: F = FreeMonoid(3, 'a')
            sage: z = F([(0,5),(1,2),(0,10),(0,2),(1,2)])
            sage: z._latex_()
            'a0^{5}a1^{2}a0^{12}a1^{2}'
        """
        s = ""
        v = self._element_list
        x = self.parent().variable_names()
        for i in range(len(v)):
            g = x[int(v[i][0])]
            e = v[i][1]
            if e == 1:
                s += "%s"%g
            else:
                s += "%s^{%s}"%(g,e)
        if len(s) == 0: s = "1"
        return s

    def __mul__(self, y):
        """
        Multiply 2 free monoid elements.

        EXAMPLES:
            sage: a = FreeMonoid(5, 'a').gens()
            sage: x = a[0] * a[1] * a[4]**3
            sage: y = a[4] * a[0] * a[1]
            sage: x*y
            a0*a1*a4^4*a0*a1
        """
        if not isinstance(y, FreeMonoidElement):
            raise TypeError, "Argument y (= %s) is of wrong type."%y
        M = self.parent()
        z = M(1)
        x_elt = self._element_list
        y_elt = y._element_list
        if len(x_elt) == 0:
            z._element_list = y_elt
        elif len(y_elt) == 0:
            z._element_list = x_elt
        else:
            k = len(x_elt)-1
            if x_elt[k][0] != y_elt[0][0]:
                z._element_list = x_elt + y_elt
            else:
                m = (y_elt[0][0],x_elt[k][1]+y_elt[0][1])
                z._element_list = x_elt[0:k] + [ m ] + y_elt[1:]
        return z

    def __pow__(self, n):
        """
        Return the $n$-th power of this monoid element.

        EXAMPLES:
            sage: a = FreeMonoid(5, 'a').gens()
            sage: x = a[0]*a[1]*a[4]**3
            sage: x**3
            a0*a1*a4^3*a0*a1*a4^3*a0*a1*a4^3
            sage: x**0
            1

        Note that raising to a negative power is \emph{not} a constructor
        for an element of the corresponding free group (yet).
            sage: x**(-1)
            Traceback (most recent call last):
            ...
            IndexError: Argument n (= -1) must be non-negative.
        """
        if not isinstance(n, (int, long, Integer)):
            raise TypeError, "Argument n (= %s) must be an integer."%n
        if n < 0:
            raise IndexError, "Argument n (= %s) must be non-negative."%n
        elif n == 0:
            return self.parent()(1)
        elif n == 1:
            return self
        elif n == 2:
            return self * self
        k = n//2
        return self**k * self**(n-k)

    def __len__(self):
        """
        Return the number of products that occur in this monoid element.
        For example, the length of the identity is 0, and the length
        of the monoid $x_0^2x_1$ is three.

        EXAMPLES:
            sage: F = FreeMonoid(3, 'a')
            sage: z = F(1)
            sage: len(z)
            0
            sage: a = F.gens()
            sage: len(a[0]**2 * a[1])
            3
        """
        s = 0
        for x in self._element_list:
            s += x[1]
        return s

    def __cmp__(self,y):
##         """
##         The comparison operator, defined via x = self:
##             x < y <=> x.__cmp__(y) == -1
##             x == y <=> x.__cmp__(y) == 0
##             x > y <=> x.__cmp__(y) == 1
##         It is not possible to use __cmp__ to define a
##         non-totally ordered poset.
##         Question: How can the operators <, >, ==, !=,
##         <=, and >= be defined for a general poset?
##         N.B. An equal operator __equal__ may or may not
##         have been introduced to define == and != but can
##         not be used in conjuction with __cmp__.
##        """
        if not isinstance(y,FreeMonoidElement) or y.parent() != self.parent():
            #raise TypeError, "Argument y (= %s) is of the wrong type."%y
            return 1
        n = len(self)
        m = len(y)
        if n < m:
            return -1
        elif m < n:
            return 1
        elif n == 0:
            return 0 # n = m = 0 hence x = y = 1
        x_elt = self._element_list
        y_elt = y._element_list
        for i in range(len(x_elt)):
            k = x_elt[i][0]
            l = y_elt[i][0]
            if k < l:
                return -1
            elif k > l:
                return 1
            e = x_elt[i][1]
            f = y_elt[i][1]
            if e < f:
                # x_elt is longer so compare next index
                if x_elt[i+1][0] < l:
                    return -1
                else:
                    return 1
            elif f < e:
                # y_elt is longer so compare next index
                if k < y_elt[i+1][0]:
                    return -1
                else:
                    return 1
        return 0 # x = self and y are equal


