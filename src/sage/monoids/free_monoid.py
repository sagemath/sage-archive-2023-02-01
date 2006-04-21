r"""
Free Monoids

AUTHOR:
    -- David Kohel <kohel@maths.usyd.edu.au>, 2005/09

SAGE supports free monoids on any prescribed finite number $n\geq 0$
of generators.  Use the \code{FreeMonoid} function to create a free
monoid, and the \code{gen} and \code{gens} functions to obtain the
corresponding generators.  You can print the generators as arbitrary
strings using the optional \code{names} argument to the
\code{FreeMonoid} function.
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

from sage.ext.integer import Integer
from sage.structure.gens import Generators
from free_monoid_element import FreeMonoidElement

from monoid import Monoid_class

import weakref

_cache = {}
def FreeMonoid(n, names=None):
    """
    Returns a free monoid on $n$ generators.

    INPUT:
        n -- integer
        names -- (optional) names of generators

    OUTPUT:
        free abelian monoid

    EXAMPLES:
        sage: FreeMonoid(0)
        Free monoid on 0 generators ()
        sage: F = FreeMonoid(5, names = list("abcde"))
        sage: F
        Free monoid on 5 generators (a, b, c, d, e)
        sage: F(1)
        1
        sage: (a, b, c, d, e) = F.gens()
        sage: mul([ a, b, a, c, b, d, c, d ])
        a*b*a*c*b*d*c*d
    """
    global _cache
    if isinstance(names, list):
        key = (n, tuple(names))
    else:
        key = (n, names)
    if _cache.has_key(key):
        M = _cache[key]()
        if not M is None:
            return M
    M = FreeMonoid_class(n, names)
    _cache[key] = weakref.ref(M)
    return M

def is_FreeMonoid(x):
    """
    Return True if $x$ is a free monoid.

    EXAMPLES:
        sage: is_FreeMonoid(5)
        False
        sage: is_FreeMonoid(FreeMonoid(7))
        True
        sage: is_FreeMonoid(FreeAbelianMonoid(7))
        False
        sage: is_FreeMonoid(FreeAbelianMonoid(0))
        False
    """
    return isinstance(x, FreeMonoid_class)

class FreeMonoid_class(Monoid_class):
    """
    The free monoid on $n$ generators.
    """
    def __init__(self, n, names=None):
        """
        Create free monoid on $n$ generators.

        INPUT:
            n -- integer
            names -- (optional) variable name or list of variable names

        EXAMPLES:
            sage: F = FreeMonoid(3)
            sage: F
            Free monoid on 3 generators (x0, x1, x2)
            sage: x = F.gens()
            sage: x[0]*x[1]**5 * (x[0]*x[2])
            x0*x1^5*x0*x2
            sage: F = FreeMonoid(3, 'a')
            sage: F
            Free monoid on 3 generators (a0, a1, a2)


            sage: M = FreeMonoid(3, names=['a','b','c'])
            sage: loads(M.dumps()) == M
            True
        """
        if not isinstance(n, (int, long, Integer)):
            raise TypeError, "n (=%s) must be an integer."%n
        if n < 0:
            raise ValueError, "n (=%s) must be nonnegative."%n
        self.__ngens = int(n)
        self.assign_names(names)

    def __cmp__(self, other):
        if not isinstance(other, FreeMonoid_class):
            return -1
        c = cmp(self.__ngens, other.__ngens)
        if c: return c
        if self.variable_names() == other.variable_names():
            return 0
        return 1

    def __repr__(self):
        return "Free monoid on %s generators %s"%(self.__ngens,self.gens())

    def __call__(self, x, check=True):
        """
        Return $x$ coerced into this free monoid.

        One can create a free monoid from the integer 1 and from a
        list of 2-tuples of integers $(i,j)$, where $(i,j)$
        corresponds to $x_i^j$, where $x_i$ is the $i$th generator.

        EXAMPLES:
            sage: F = FreeMonoid(3, 'a')
            sage: F(1)
            1
            sage: F(F.gen(0))
            a0
            sage: F(0)
            Traceback (most recent call last):
            ...
            TypeError: Argument x (= 0) is of the wrong type.

            sage: F([(0,5),(1,2),(0,10),(0,2),(1,2)])
            a0^5*a1^2*a0^12*a1^2

        """
        ## There should really some careful type checking here...
        if isinstance(x, FreeMonoidElement) and x.parent() == self:
            return x
        elif isinstance(x, (int, long, Integer)) and x == 1:
            return FreeMonoidElement(self, x, check)
        elif isinstance(x, list):
            return FreeMonoidElement(self, x, check)
        else:
            raise TypeError, "Argument x (= %s) is of the wrong type."%x

    def __contains__(self, x):
        return isinstance(x, FreeMonoidElement) and x.parent() == self

    def gen(self,i=0):
        """
        The $i$-th generator of the monoid.

        INPUT:
            i -- integer (default: 0)

        EXAMPLES:
            sage: F = FreeMonoid(3, 'a')
            sage: F.gen(1)
            a1
            sage: F.gen(2)
            a2
            sage: F.gen(5)
            Traceback (most recent call last):
            ...
            IndexError: Argument i (= 5) must be between 0 and 2.
        """
        n = self.__ngens
        if i < 0 or not i < n:
            raise IndexError, "Argument i (= %s) must be between 0 and %s."%(i, n-1)
        return FreeMonoidElement(self,[(Integer(i),Integer(1))])

    def ngens(self):
        """
        The number of free generators of the monoid.

        EXAMPLES:
            sage: F = FreeMonoid(2005)
            sage: F.ngens()
            2005
        """
        return self.__ngens



