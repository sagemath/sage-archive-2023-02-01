"""
Free algebras

AUTHOR: David Kohel, 2005-09
    William Stein 2006-11-01 -- add all doctests; implemented many things.

EXAMPLES:
    sage: F = FreeAlgebra(ZZ,3,'x,y,z')
    sage: F.base_ring()
    Integer Ring
    sage: G = FreeAlgebra(F, 2, 'm,n'); G
    Free Algebra on 2 generators (m, n) over Free Algebra on 3 generators (x, y, z) over Integer Ring
    sage: G.base_ring()
    Free Algebra on 3 generators (x, y, z) over Integer Ring

TESTS:
    sage: F = FreeAlgebra(GF(5),3,'x')
    sage: F == loads(dumps(F))
    True

    sage: F.<x,y,z> = FreeAlgebra(GF(5),3)
    sage: F == loads(dumps(F))
    True

    sage: F = FreeAlgebra(GF(5),3, ['xx', 'zba', 'Y'])
    sage: F == loads(dumps(F))
    True

    sage: F = FreeAlgebra(GF(5),3, 'abc')
    sage: F == loads(dumps(F))
    True

    sage: F = FreeAlgebra(FreeAlgebra(ZZ,1,'a'), 2, 'x')
    sage: F == loads(dumps(F))
    True

"""

#*****************************************************************************
#  Copyright (C) 2005 David Kohel <kohel@maths.usyd.edu>
#  Copyright (C) 2005,2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.ring import Ring
from sage.rings.integer import Integer

from sage.monoids.free_monoid import FreeMonoid
from sage.monoids.free_monoid_element import FreeMonoidElement

from sage.algebras.algebra import Algebra
from sage.algebras.free_algebra_element import FreeAlgebraElement

import sage.structure.parent_gens


def FreeAlgebra(R, n, names):
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
        sage: FreeAlgebra(GF(5),3,'x')
        Free Algebra on 3 generators (x0, x1, x2) over Finite Field of size 5
        sage: F.<x,y,z> = FreeAlgebra(GF(5),3)
        sage: (x+y+z)^2
        x^2 + x*y + x*z + y*x + y^2 + y*z + z*x + z*y + z^2
        sage: FreeAlgebra(GF(5),3, ['xx', 'zba', 'Y'])
        Free Algebra on 3 generators (xx, zba, Y) over Finite Field of size 5
        sage: FreeAlgebra(GF(5),3, 'abc')
        Free Algebra on 3 generators (a, b, c) over Finite Field of size 5
        sage: FreeAlgebra(GF(5),1, 'z')
        Free Algebra on 1 generators (z,) over Finite Field of size 5
        sage: FreeAlgebra(GF(5),1, ['alpha'])
        Free Algebra on 1 generators (alpha,) over Finite Field of size 5
        sage: FreeAlgebra(FreeAlgebra(ZZ,1,'a'), 2, 'x')
        Free Algebra on 2 generators (x0, x1) over Free Algebra on 1 generators (a,) over Integer Ring

    Free algebras are globally unique:
        sage: F = FreeAlgebra(ZZ,3,'x,y,z')
        sage: G = FreeAlgebra(ZZ,3,'x,y,z')
        sage: F is G
        True

    Free algebras commute with their base ring.
        sage: K.<a,b> = FreeAlgebra(QQ,2)
        sage: K.is_commutative()
        False
        sage: L.<c> = FreeAlgebra(K,1)
        sage: L.is_commutative()
        False
        sage: s = a*b^2 * c^3; s
        a*b^2*c^3
        sage: parent(s)
        Free Algebra on 1 generators (c,) over Free Algebra on 2 generators (a, b) over Rational Field
        sage: c^3 * a * b^2
        a*b^2*c^3
    """
    names = sage.structure.parent_gens.normalize_names(n, names)
    return cache(R, n, names)

def is_FreeAlgebra(x):
    """
    Return True if x is a free algebra; otherwise, return False.

    EXAMPLES:
        sage: is_FreeAlgebra(5)
        False
        sage: is_FreeAlgebra(ZZ)
        False
        sage: is_FreeAlgebra(FreeAlgebra(ZZ,100,'x'))
        True
    """
    return isinstance(x, FreeAlgebra_generic)


class FreeAlgebra_generic(Algebra):
    """
    The free algebra on $n$ generators over a base ring.

    EXAMPLES:
        sage: F.<x,y,z> = FreeAlgebra(QQ, 3); F
        Free Algebra on 3 generators (x, y, z) over Rational Field
        sage: mul(F.gens())
        x*y*z
        sage: mul([ F.gen(i%3) for i in range(12) ])
        x*y*z*x*y*z*x*y*z*x*y*z
        sage: mul([ F.gen(i%3) for i in range(12) ]) + mul([ F.gen(i%2) for i in range(12) ])
        x*y*x*y*x*y*x*y*x*y*x*y + x*y*z*x*y*z*x*y*z*x*y*z

        sage: (2 + x*z + x^2)^2 + (x - y)^2
        4 + 3*x^2 - x*y + 2*x*z - y*x + y^2 + x^4 + x^3*z + x*z*x^2 + x*z*x*z
    """
    def __init__(self, R, n, names):
        """
        INPUT:
            R -- ring
            n -- an integer
            names -- generator names
        """
        if not isinstance(R, Ring):
            raise TypeError, "Argument R must be a ring."
        self.__monoid = FreeMonoid(n, names=names)
        self.__ngens = n
        sage.structure.parent_gens.ParentWithGens.__init__(self, R, names)

    def is_field(self):
        if self.__ngens == 0:
            return self.base_ring().is_field()
        return False

    def is_commutative(self):
        """
        Return True if this free algebra is commutative.

        EXAMPLES:

        """
        return self.__ngens <= 1 and self.base_ring().is_commutative()

    def __cmp__(self, other):
        """
        Two free algebras are considered the same if they have the
        same base ring, number of generators and variable names.

        EXAMPLES:
            sage: F = FreeAlgebra(QQ,3,'x')
            sage: F ==  FreeAlgebra(QQ,3,'x')
            True
            sage: F is  FreeAlgebra(QQ,3,'x')
            True
            sage: F == FreeAlgebra(ZZ,3,'x')
            False
            sage: F == FreeAlgebra(QQ,4,'x')
            False
            sage: F == FreeAlgebra(QQ,3,'y')
            False
        """
        if not isinstance(other, FreeAlgebra_generic):
            return -1
        c = cmp(self.base_ring(), other.base_ring())
        if c: return c
        c = cmp(self.__ngens, other.__ngens)
        if c: return c
        c = cmp(self.variable_names(), other.variable_names())
        if c: return c
        return 0

    def _repr_(self):
        """
        Text representation of this free algebra.

        EXAMPLES:
            sage: F = FreeAlgebra(QQ,3,'x')
            sage: print F
            Free Algebra on 3 generators (x0, x1, x2) over Rational Field
            sage: F.rename('QQ<<x0,x1,x2>>')
            sage: print F
            QQ<<x0,x1,x2>>
        """
        return "Free Algebra on %s generators %s over %s"%(
            self.__ngens, self.gens(), self.base_ring())

    def __call__(self, x):
        """
        Coerce x into self.
        """
        if isinstance(x, FreeAlgebraElement):
            P = x.parent()
            if P is self:
                return x
            if not (P is self.base_ring()):
                return FreeAlgebraElement(self, x)
        # ok, not a free algebra element (or should not be viewed as one).
        F = self.__monoid
        R = self.base_ring()
        # coercion from free monoid
        if isinstance(x, FreeMonoidElement) and x.parent() is F:
            return FreeAlgebraElement(self,{x:R(1)})
        # coercion via base ring
        x = R(x)
        if x == 0:
            return FreeAlgebraElement(self,{})
        else:
            return FreeAlgebraElement(self,{F(1):x})

    def _coerce_impl(self, x):
        """
        Canonical coercion of x into self.

        Here's what canonically coerces to self:
            * this free algebra
            * the underlying monoid
            * anything that coerces to the base ring of this free algebra
            * any free algebra whose base ring coerces to the base ring of this free algebra

        EXAMPLES:
            sage: F.<x,y,z> = FreeAlgebra(GF(7),3); F
            Free Algebra on 3 generators (x, y, z) over Finite Field of size 7

        Elements of the free algebra canonically coerce in.
            sage: F._coerce_(x*y)
            x*y

        Elements of the integers coerce in, since there is a coerce map from ZZ to GF(7).
            sage: F._coerce_(1)
            1

        There is no coerce map from QQ to GF(7).
            sage: F._coerce_(2/3)
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion of element into self

        Elements of the base ring coerce in.
            sage: F._coerce_(GF(7)(5))
            5

        Elements of the correspondining moind (of monomials) coerce in:
            sage: M = F.monoid(); m = M.0*M.1^2; m
            x*y^2
            sage: F._coerce_(m)
            x*y^2

        The free algebra over ZZ on x,y,z coerces in, since ZZ coerces
        to GF(7):
            sage: G = FreeAlgebra(ZZ,3,'x,y,z')
            sage: F._coerce_(G.0^3 * G.1)
            x^3*y

        However, GF(7) doesn't coerce to ZZ, so the free algebra over
        GF(7) doesn't coerce to the one over ZZ:
            sage: G._coerce_(x^3*y)
            Traceback (most recent call last):
            ...
            TypeError: no natural map between bases of free algebras
        """
        try:
            R = x.parent()

            # monoid
            if R is self.__monoid:
                return self(x)

            # polynomial rings in the same variable over any base that coerces in:
            if is_FreeAlgebra(R):
                if R.variable_names() == self.variable_names():
                    if self.has_coerce_map_from(R.base_ring()):
                        return self(x)
                    else:
                        raise TypeError, "no natural map between bases of free algebras"

        except AttributeError:
            pass

        # any ring that coerces to the base ring of this free algebra.
        return self._coerce_try(x, [self.base_ring()])

    def gen(self,i):
        """
        The i-th generator of the algebra.

        EXAMPLES:
            sage: F = FreeAlgebra(ZZ,3,'x,y,z')
            sage: F.gen(0)
            x
        """
        n = self.__ngens
        if i < 0 or not i < n:
            raise IndexError, "Argument i (= %s) must be between 0 and %s."%(i, n-1)
        R = self.base_ring()
        F = self.__monoid
        return FreeAlgebraElement(self,{F.gen(i):R(1)})

    def quotient(self, mons, mats, names):
        """
        Returns a quotient algebra defined via the action of a free algebra A
        on a (finitely generated) free module.  The input for the quotient algebra
        is a list of monomials (in the underlying monoid for A) which form a free
        basis for the module of A, and a list of matrices, which give the action
        of the free generators of A on this monomial basis.
        """
        import free_algebra_quotient
        return free_algebra_quotient.FreeAlgebraQuotient(self, mons, mats, names)
    quo = quotient

    def ngens(self):
        """
        The number of generators of the algebra.

        EXAMPLES:
            sage: F = FreeAlgebra(ZZ,3,'x,y,z')
            sage: F.ngens()
            3
        """
        return self.__ngens

    def monoid(self):
        """
        The free monoid of generators of the algebra.

        EXAMPLES:
            sage: F = FreeAlgebra(ZZ,3,'x,y,z')
            sage: F.monoid()
            Free monoid on 3 generators (x, y, z)
        """
        return self.__monoid


from sage.misc.cache import Cache
cache = Cache(FreeAlgebra_generic)
