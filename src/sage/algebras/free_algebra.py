"""
Free algebras

AUTHOR: David Kohel, 2005-09
"""

#*****************************************************************************
#  Copyright (C) 2005 David Kohel <kohel@maths.usyd.edu>
#  Copyright (C) 2005 William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.ring import Ring
from sage.ext.integer import Integer
from sage.monoids.free_monoid import FreeMonoid
from sage.monoids.free_monoid_element import FreeMonoidElement
from sage.algebras.algebra import Algebra
from sage.algebras.free_algebra_element import FreeAlgebraElement


def FreeAlgebra(R, n, names = None):
    """
    Return the free algebra over the ring $R$ on $n$ generators with
    given names.

    INPUT:
        R -- ring
        n -- integer
        names -- string or list/tuple of n strings

    OUTPUT:
        a free algebra

    EXAMPLES:
        sage: FreeAlgebra(GF(5),3)
        Free Algebra on 3 generators (x_0, x_1, x_2) over Finite Field of size 5
        sage: FreeAlgebra(GF(5),3, ['xx', 'zba', 'Y'])
        Free Algebra on 3 generators (xx, zba, Y) over Finite Field of size 5
        sage: FreeAlgebra(GF(5),3, 'abc')
        Free Algebra on 3 generators (a, b, c) over Finite Field of size 5
        sage: FreeAlgebra(GF(5),1, 'z')
        Free Algebra on 1 generators (z,) over Finite Field of size 5
        sage: FreeAlgebra(GF(5),1, ['alpha'])
        Free Algebra on 1 generators (alpha,) over Finite Field of size 5
        sage: FreeAlgebra(FreeAlgebra(ZZ,1), 2)
        Free Algebra on 2 generators (x_0, x_1) over Free Algebra on 1 generators (x,) over Integer Ring
    """
    return cache(R,n, cache.format_names(names, n))

def is_FreeAlgebra(x):
    return isinstance(x, FreeAlgebra_generic)


class FreeAlgebra_generic(Algebra):
    def __init__(self, R, n, names = None):
        """
        Returns the free algebra on $n$ generators.

        EXAMPLES:
            sage: F = FreeAlgebra(QQ,ZZ(3),names=("x","y","z"))
            sage: mul([ F.gen(i) for i in range(3) ])
            x*y*z
            sage: mul([ F.gen(i%3) for i in range(12) ])
            x*y*z*x*y*z*x*y*z*x*y*z
            sage: (x,y,z) = F.gens()
            sage: (2 + x*z + x**2)**2 + (x - y)**2
            4 + 3*x^2 - x*y + 2*x*z - y*x + y^2 + x^4 + x^3*z + x*z*x^2 + x*z*x*z
        """
        if not isinstance(R, Ring):
            raise TypeError, "Argument R must be a ring."
        self.__base_ring = R
        self.__monoid = FreeMonoid(n, names=names)
        self.__ngens = n
        self.assign_names(names)

    def __cmp__(self, other):
        if not isinstance(other, FreeAlgebra_generic):
            return -1
        c = cmp(self.__ngens, other.__ngens)
        if c: return c
        if self.variable_names() == other.variable_names():
            return 0
        return 1

    def __repr__(self):
        return "Free Algebra on %s generators %s over %s"%(
            self.__ngens, self.gens(), self.__base_ring)

    def __call__(self, x, canonical=False):
        if isinstance(x, FreeAlgebraElement) and x.parent() == self:
            return x
        F = self.__monoid
        R = self.__base_ring
        if isinstance(x, (int, long, Integer)):
            if x == 0:
                return FreeAlgebraElement(self,{})
            else:
                return FreeAlgebraElement(self,{F(1):R(x)})
        elif isinstance(x, Ring) and x.parent() is R:
            if x == 0:
                return FreeAlgebraElement(self,{})
            else:
                return FreeAlgebraElement(self,{F(1):x})
        elif isinstance(x, FreeMonoidElement) and x.parent() is F:
            return FreeAlgebraElement(self,{x:R(1)})
        else:
            raise TypeError, "Argument x (= %s) if of invalid type."%x

    def __contains__(self, x):
        return isinstance(x, FreeAlgebraElement) and x.parent() == self

    def gen(self,i):
        """
        The i-th generator of the algebra.
        """
        n = self.__ngens
        if i < 0 or not i < n:
            raise IndexError, "Argument i (= %s) must be between 0 and %s."%(i, n-1)
        R = self.__base_ring
        F = self.__monoid
        return FreeAlgebraElement(self,{F.gen(i):R(1)})

    def ngens(self):
        """
        The number of generators of the algebra.
        """
        return self.__ngens

    def assign_names(self,names):
        """
        Assign the printing names for the generators; this will have the unfortunate
        effect of overwriting the names for the covering algebra; this also does not
        overwrite the return value of names() for the Algebra.
        """
        #self._Algebra.assign_names(names) # where is this located?
        self.monoid().assign_names(names)


    def base_ring(self):
        return self.__base_ring

    def monoid(self):
        """
        The free monoid of generators of the algebra.
        """
        return self.__monoid


from sage.misc.cache import Cache
cache = Cache(FreeAlgebra_generic)
