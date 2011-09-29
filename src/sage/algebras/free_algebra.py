"""
Free algebras

AUTHORS:

- David Kohel (2005-09)

- William Stein (2006-11-01): add all doctests; implemented many
  things.

- Simon King (2011-04): Put free algebras into the category framework.

EXAMPLES::

    sage: F = FreeAlgebra(ZZ,3,'x,y,z')
    sage: F.base_ring()
    Integer Ring
    sage: G = FreeAlgebra(F, 2, 'm,n'); G
    Free Algebra on 2 generators (m, n) over Free Algebra on 3 generators (x, y, z) over Integer Ring
    sage: G.base_ring()
    Free Algebra on 3 generators (x, y, z) over Integer Ring

TESTS::

    sage: F = FreeAlgebra(GF(5),3,'x')
    sage: TestSuite(F).run()

::

    sage: F.<x,y,z> = FreeAlgebra(GF(5),3)
    sage: TestSuite(F).run()

::

    sage: F = FreeAlgebra(GF(5),3, ['xx', 'zba', 'Y'])
    sage: TestSuite(F).run()

::

    sage: F = FreeAlgebra(GF(5),3, 'abc')
    sage: TestSuite(F).run()

::

    sage: F = FreeAlgebra(FreeAlgebra(ZZ,1,'a'), 2, 'x')
    sage: TestSuite(F).run()

"""

#*****************************************************************************
#  Copyright (C) 2005 David Kohel <kohel@maths.usyd.edu>
#  Copyright (C) 2005,2006 William Stein <wstein@gmail.com>
#  Copyright (C) 2011 Simon King <simon.king@uni-jena.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.ring import Ring

from sage.monoids.free_monoid import FreeMonoid
from sage.monoids.free_monoid_element import FreeMonoidElement

from sage.algebras.algebra import Algebra
from sage.algebras.free_algebra_element import FreeAlgebraElement

import sage.structure.parent_gens

def FreeAlgebra(R, n, names):
    """
    Return the free algebra over the ring `R` on `n`
    generators with given names.

    INPUT:

    -  ``R`` - ring
    -  ``n`` - integer
    -  ``names`` - string or list/tuple of n strings

    OUTPUT:

    A free algebra.

    EXAMPLES::

        sage: FreeAlgebra(GF(5),3,'x')
        Free Algebra on 3 generators (x0, x1, x2) over Finite Field of size 5
        sage: F.<x,y,z> = FreeAlgebra(GF(5),3)
        sage: (x+y+z)^2
        x^2 + x*y + x*z + y*x + y^2 + y*z + z*x + z*y + z^2
        sage: FreeAlgebra(GF(5),3, 'xx, zba, Y')
        Free Algebra on 3 generators (xx, zba, Y) over Finite Field of size 5
        sage: FreeAlgebra(GF(5),3, 'abc')
        Free Algebra on 3 generators (a, b, c) over Finite Field of size 5
        sage: FreeAlgebra(GF(5),1, 'z')
        Free Algebra on 1 generators (z,) over Finite Field of size 5
        sage: FreeAlgebra(GF(5),1, ['alpha'])
        Free Algebra on 1 generators (alpha,) over Finite Field of size 5
        sage: FreeAlgebra(FreeAlgebra(ZZ,1,'a'), 2, 'x')
        Free Algebra on 2 generators (x0, x1) over Free Algebra on 1 generators (a,) over Integer Ring

    Free algebras are globally unique::

        sage: F = FreeAlgebra(ZZ,3,'x,y,z')
        sage: G = FreeAlgebra(ZZ,3,'x,y,z')
        sage: F is G
        True

    Free algebras commute with their base ring.
    ::

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

    EXAMPLES::

        sage: from sage.algebras.free_algebra import is_FreeAlgebra
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
    The free algebra on `n` generators over a base ring.

    EXAMPLES::

        sage: F.<x,y,z> = FreeAlgebra(QQ, 3); F
        Free Algebra on 3 generators (x, y, z) over Rational Field
        sage: mul(F.gens())
        x*y*z
        sage: mul([ F.gen(i%3) for i in range(12) ])
        x*y*z*x*y*z*x*y*z*x*y*z
        sage: mul([ F.gen(i%3) for i in range(12) ]) + mul([ F.gen(i%2) for i in range(12) ])
        x*y*x*y*x*y*x*y*x*y*x*y + x*y*z*x*y*z*x*y*z*x*y*z
        sage: (2 + x*z + x^2)^2 + (x - y)^2
        4 + 5*x^2 - x*y + 4*x*z - y*x + y^2 + x^4 + x^3*z + x*z*x^2 + x*z*x*z
    """
    Element = FreeAlgebraElement
    def __init__(self, R, n, names):
        """
        The free algebra on `n` generators over a base ring.

        INPUT:

        -  ``R`` - ring
        -  ``n`` - an integer
        -  ``names`` - generator names

        EXAMPLES::

            sage: F.<x,y,z> = FreeAlgebra(QQ, 3); F # indirect doctet
            Free Algebra on 3 generators (x, y, z) over Rational Field
        """
        if not isinstance(R, Ring):
            raise TypeError, "Argument R must be a ring."
        self.__monoid = FreeMonoid(n, names=names)
        self.__ngens = n
        #sage.structure.parent_gens.ParentWithGens.__init__(self, R, names)
        Algebra.__init__(self, R,names=names)

    def is_field(self, proof = True):
        """
        Return True if this Free Algebra is a field, which is only if the
        base ring is a field and there are no generators

        EXAMPLES::

            sage: A=FreeAlgebra(QQ,0,'')
            sage: A.is_field()
            True
            sage: A=FreeAlgebra(QQ,1,'x')
            sage: A.is_field()
            False
        """
        if self.__ngens == 0:
            return self.base_ring().is_field(proof)
        return False

    def is_commutative(self):
        """
        Return True if this free algebra is commutative.

        EXAMPLES::

            sage: R.<x> = FreeAlgebra(QQ,1)
            sage: R.is_commutative()
            True
            sage: R.<x,y> = FreeAlgebra(QQ,2)
            sage: R.is_commutative()
            False
        """
        return self.__ngens <= 1 and self.base_ring().is_commutative()

    def __cmp__(self, other):
        """
        Two free algebras are considered the same if they have the same
        base ring, number of generators and variable names.

        EXAMPLES::

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

        EXAMPLES::

            sage: F = FreeAlgebra(QQ,3,'x')
            sage: print F  # indirect doctest
            Free Algebra on 3 generators (x0, x1, x2) over Rational Field
            sage: F.rename('QQ<<x0,x1,x2>>')
            sage: print F #indirect doctest
            QQ<<x0,x1,x2>>
            sage: FreeAlgebra(ZZ,1,['a'])
            Free Algebra on 1 generators (a,) over Integer Ring
        """
        return "Free Algebra on %s generators %s over %s"%(
            self.__ngens, self.gens(), self.base_ring())

    def _element_constructor_(self, x):
        """
       Coerce x into self.

       EXAMPLES::

           sage: R.<x,y> = FreeAlgebra(QQ,2)
           sage: R(3) # indirect doctest
           3
       """
        if isinstance(x, FreeAlgebraElement):
            P = x.parent()
            if P is self:
                return x
            if not (P is self.base_ring()):
                return self.element_class(self, x)
        # ok, not a free algebra element (or should not be viewed as one).
        F = self.__monoid
        R = self.base_ring()
        # coercion from free monoid
        if isinstance(x, FreeMonoidElement) and x.parent() is F:
            return self.element_class(self,{x:R(1)})
        # coercion via base ring
        x = R(x)
        if x == 0:
            return self.element_class(self,{})
        else:
            return self.element_class(self,{F(1):x})

    def _coerce_impl(self, x):
        """
        Canonical coercion of x into self.

        Here's what canonically coerces to self:

        - this free algebra

        - the underlying monoid

        - anything that coerces to the base ring of this free algebra

        - any free algebra whose base ring coerces to the base ring of
          this free algebra

        EXAMPLES::

            sage: F.<x,y,z> = FreeAlgebra(GF(7),3); F
            Free Algebra on 3 generators (x, y, z) over Finite Field of size 7

        Elements of the free algebra canonically coerce in.

        ::

            sage: F._coerce_(x*y) # indirect doctest
            x*y

        Elements of the integers coerce in, since there is a coerce map
        from ZZ to GF(7).

        ::

            sage: F._coerce_(1)       # indirect doctest
            1

        There is no coerce map from QQ to GF(7).

        ::

            sage: F._coerce_(2/3)
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Rational Field to Free Algebra
            on 3 generators (x, y, z) over Finite Field of size 7

        Elements of the base ring coerce in.

        ::

            sage: F._coerce_(GF(7)(5))
            5

        Elements of the corresponding monoid (of monomials) coerce in::

            sage: M = F.monoid(); m = M.0*M.1^2; m
            x*y^2
            sage: F._coerce_(m)
            x*y^2

        The free algebra over ZZ on x,y,z coerces in, since ZZ coerces to
        GF(7)::

            sage: G = FreeAlgebra(ZZ,3,'x,y,z')
            sage: F._coerce_(G.0^3 * G.1)
            x^3*y

        However, GF(7) doesn't coerce to ZZ, so the free algebra over GF(7)
        doesn't coerce to the one over ZZ::

            sage: G._coerce_(x^3*y)
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Free Algebra on 3 generators
            (x, y, z) over Finite Field of size 7 to Free Algebra on 3
            generators (x, y, z) over Integer Ring
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

    def _coerce_map_from_(self, R):
        if self.__monoid.has_coerce_map_from(R):
            return True

        # polynomial rings in the same variable over any base that coerces in:
        if is_FreeAlgebra(R):
            if R.variable_names() == self.variable_names():
                if self.base_ring().has_coerce_map_from(R.base_ring()):
                    return True
                else:
                    return False

        return self.base_ring().has_coerce_map_from(R)

    def gen(self,i):
        """
        The i-th generator of the algebra.

        EXAMPLES::

            sage: F = FreeAlgebra(ZZ,3,'x,y,z')
            sage: F.gen(0)
            x
        """
        n = self.__ngens
        if i < 0 or not i < n:
            raise IndexError, "Argument i (= %s) must be between 0 and %s."%(i, n-1)
        R = self.base_ring()
        F = self.__monoid
        return self.element_class(self,{F.gen(i):R(1)})

    def quotient(self, mons, mats, names):
        """
        Returns a quotient algebra.

        The quotient algebra is defined via the action of a free algebra
        A on a (finitely generated) free module. The input for the quotient
        algebra is a list of monomials (in the underlying monoid for A)
        which form a free basis for the module of A, and a list of
        matrices, which give the action of the free generators of A on this
        monomial basis.

        EXAMPLE:

        Here is the quaternion algebra defined in terms of three generators::

            sage: n = 3
            sage: A = FreeAlgebra(QQ,n,'i')
            sage: F = A.monoid()
            sage: i, j, k = F.gens()
            sage: mons = [ F(1), i, j, k ]
            sage: M = MatrixSpace(QQ,4)
            sage: mats = [M([0,1,0,0, -1,0,0,0, 0,0,0,-1, 0,0,1,0]),  M([0,0,1,0, 0,0,0,1, -1,0,0,0, 0,-1,0,0]),  M([0,0,0,1, 0,0,-1,0, 0,1,0,0, -1,0,0,0]) ]
            sage: H.<i,j,k> = A.quotient(mons, mats); H
            Free algebra quotient on 3 generators ('i', 'j', 'k') and dimension 4 over Rational Field
        """
        import free_algebra_quotient
        return free_algebra_quotient.FreeAlgebraQuotient(self, mons, mats, names)
    quo = quotient

    def ngens(self):
        """
        The number of generators of the algebra.

        EXAMPLES::

            sage: F = FreeAlgebra(ZZ,3,'x,y,z')
            sage: F.ngens()
            3
        """
        return self.__ngens

    def monoid(self):
        """
        The free monoid of generators of the algebra.

        EXAMPLES::

            sage: F = FreeAlgebra(ZZ,3,'x,y,z')
            sage: F.monoid()
            Free monoid on 3 generators (x, y, z)
        """
        return self.__monoid

    def g_algebra(self, relations, names=None, order='degrevlex', check = True):
        """
        The G-Algebra derived from this algebra by relations.
        By default is assumed, that two variables commute.

        TODO:

        - Coercion doesn't work yet, there is some cheating about assumptions
        - The optional argument ``check`` controls checking the degeneracy
          conditions. Furthermore, the default values interfere with
          non-degeneracy conditions.

        EXAMPLES::

            sage: A.<x,y,z>=FreeAlgebra(QQ,3)
            sage: G=A.g_algebra({y*x:-x*y})
            sage: (x,y,z)=G.gens()
            sage: x*y
            x*y
            sage: y*x
            -x*y
            sage: z*x
            x*z
            sage: (x,y,z)=A.gens()
            sage: G=A.g_algebra({y*x:-x*y+1})
            sage: (x,y,z)=G.gens()
            sage: y*x
            -x*y + 1
            sage: (x,y,z)=A.gens()
            sage: G=A.g_algebra({y*x:-x*y+z})
            sage: (x,y,z)=G.gens()
            sage: y*x
            -x*y + z
        """
        from sage.matrix.constructor  import Matrix

        base_ring=self.base_ring()
        n=self.ngens()
        cmat=Matrix(base_ring,n)
        dmat=Matrix(self,n)
        for i in xrange(n):
            for j in xrange(i+1,n):
                cmat[i,j]=1
        for (to_commute,commuted) in relations.iteritems():
            #This is dirty, coercion is broken
            assert isinstance(to_commute,FreeAlgebraElement), to_commute.__class__
            assert isinstance(commuted,FreeAlgebraElement), commuted
            ((v1,e1),(v2,e2))=list(list(to_commute)[0][1])
            assert e1==1
            assert e2==1
            assert v1>v2
            c_coef=None
            d_poly=None
            for (c,m) in commuted:
                if list(m)==[(v2,1),(v1,1)]:
                    c_coef=c
                    #buggy coercion workaround
                    d_poly=commuted-self(c)*self(m)
                    break
            assert not c_coef is None,list(m)
            v2_ind = self.gens().index(v2)
            v1_ind = self.gens().index(v1)
            cmat[v2_ind,v1_ind]=c_coef
            if d_poly:
                dmat[v2_ind,v1_ind]=d_poly
        from sage.rings.polynomial.plural import g_Algebra
        return g_Algebra(base_ring, cmat, dmat, names = names or self.variable_names(),
                         order=order, check=check)


from sage.misc.cache import Cache
cache = Cache(FreeAlgebra_generic)
