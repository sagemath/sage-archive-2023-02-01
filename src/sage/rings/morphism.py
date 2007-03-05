r"""
Homomorphisms of rings

We give a large number of examples of ring homomorphisms.

EXAMPLE: Natural inclusion $\Z \hookrightarrow \Q$.
    sage: H = Hom(ZZ, QQ)
    sage: phi = H([1])
    sage: phi(10)
    10
    sage: phi(3/1)
    3
    sage: phi(2/3)
    Traceback (most recent call last):
    ...
    TypeError: 2/3 must be coercible into Integer Ring

There is no homomorphism in the other direction:
    sage: H = Hom(QQ, ZZ)
    sage: H([1])
    Traceback (most recent call last):
    ...
    TypeError: images do not define a valid homomorphism

EXAMPLE: Reduction to finite field.
    sage: H = Hom(ZZ, GF(9, 'a'))
    sage: phi = H([1])
    sage: phi(5)
    2
    sage: psi = H([4])
    sage: psi(5)
    2

EXAMPLE: Map from single variable polynomial ring.
    sage: R, x = PolynomialRing(ZZ, 'x').objgen()
    sage: phi = R.hom([2], GF(5))
    sage: phi
    Ring morphism:
      From: Univariate Polynomial Ring in x over Integer Ring
      To:   Finite Field of size 5
      Defn: x |--> 2
    sage: phi(x + 12)
    4

EXAMPLE: Identity map on the real numbers.
    sage: f = RR.hom([RR(1)]); f
    Ring endomorphism of Real Field with 53 bits of precision
      Defn: 1.00000000000000 |--> 1.00000000000000
    sage: f(2.5)
    2.50000000000000
    sage: f = RR.hom( [2.0] )
    Traceback (most recent call last):
    ...
    TypeError: images do not define a valid homomorphism

EXAMPLE: Homomorphism from one precision of field to another.

From smaller to bigger doesn't make sense:
    sage: R200 = RealField(200)
    sage: f = RR.hom( R200 )
    Traceback (most recent call last):
    ...
    TypeError: Natural coercion morphism from Real Field with 53 bits of precision to Real Field with 200 bits of precision not defined.

From bigger to small does:
    sage: f = RR.hom( RealField(15) )
    sage: f(2.5)
    2.500
    sage: f(RR.pi())
    3.142

EXAMPLE: Inclusion map from the reals to the complexes:
    sage: i = RR.hom([CC(1)]); i
    Ring morphism:
      From: Real Field with 53 bits of precision
      To:   Complex Field with 53 bits of precision
      Defn: 1.00000000000000 |--> 1.00000000000000
    sage: i(RR('3.1'))
    3.10000000000000

EXAMPLE: A map from a multivariate polynomial ring to itself:
    sage: R.<x,y,z> = PolynomialRing(QQ,3)
    sage: phi = R.hom([y,z,x^2]); phi
    Ring endomorphism of Polynomial Ring in x, y, z over Rational Field
      Defn: x |--> y
            y |--> z
            z |--> x^2
    sage: phi(x+y+z)
    z + y + x^2

EXAMPLE: An endomorphism of a quotient of a multi-variate polynomial ring:
    sage: R.<x,y> = PolynomialRing(QQ)
    sage: S.<a,b> = quo(R, ideal(1 + y^2))
    sage: phi = S.hom([a^2, -b])
    sage: phi
    Ring endomorphism of Quotient of Polynomial Ring in x, y over Rational Field by the ideal (1 + y^2)
      Defn: a |--> a^2
            b |--> -1*b
    sage: phi(b)
    -1*b
    sage: phi(a^2 + b^2)
    -1 + a^4

EXAMPLE: The reduction map from the integers to the integers modulo 8,
viewed as a quotient ring:

    sage: R = ZZ.quo(8*ZZ)
    sage: pi = R.cover()
    sage: pi
    Ring morphism:
      From: Integer Ring
      To:   Ring of integers modulo 8
      Defn: Natural quotient map
    sage: pi.domain()
    Integer Ring
    sage: pi.codomain()
    Ring of integers modulo 8
    sage: pi(10)
    2
    sage: pi.lift()
    Set-theoretic ring morphism:
      From: Ring of integers modulo 8
      To:   Integer Ring
      Defn: Choice of lifting map
    sage: pi.lift(13)
    5


EXAMPLE: Inclusion of GF(2) into GF(4,'a').
    sage: k = GF(2)
    sage: i = k.hom(GF(4, 'a'))
    sage: i
    Coercion morphism:
      From: Finite Field of size 2
      To:   Finite Field in a of size 2^2
    sage: i(0)
    0
    sage: a = i(1); a.parent()
    Finite Field in a of size 2^2

We next compose the inclusion with reduction from the integers to GF(2).
    sage: pi = ZZ.hom(k)
    sage: pi
    Coercion morphism:
      From: Integer Ring
      To:   Finite Field of size 2
    sage: f = i * pi
    sage: f
    Composite morphism:
      From: Integer Ring
      To:   Finite Field in a of size 2^2
      Defn:   Coercion morphism:
              From: Integer Ring
              To:   Finite Field of size 2
            then
              Coercion morphism:
              From: Finite Field of size 2
              To:   Finite Field in a of size 2^2
    sage: a = f(5); a
    1
    sage: a.parent()
    Finite Field in a of size 2^2

EXAMPLE: Inclusion from $\Q$ to the 3-adic field.
    sage: phi = QQ.hom(pAdicField(3))
    sage: phi
    Coercion morphism:
      From: Rational Field
      To:   3-adic Field
    sage: phi.codomain()
    3-adic Field
    sage: phi(394)
    1 + 2*3 + 3^2 + 2*3^3 + 3^4 + 3^5 + ... + O(3^Infinity)

EXAMPLE: An automorphism of a quotient of a univariate polynomial ring.
    sage: R.<x> = PolynomialRing(QQ)
    sage: S.<sqrt2> = R.quo(x^2-2)
    sage: sqrt2^2
    2
    sage: (3+sqrt2)^10
    993054*sqrt2 + 1404491
    sage: c = S.hom([-sqrt2])
    sage: c(1+sqrt2)
    -sqrt2 + 1

Note that \sage verifies that the morphism is valid:
    sage: (1 - sqrt2)^2
    -2*sqrt2 + 3
    sage: c = S.hom([1-sqrt2])    # this is not valid
    Traceback (most recent call last):
    ...
    TypeError: images do not define a valid homomorphism

EXAMPLE: Endomorphism of power series ring.
    sage: R.<t> = PowerSeriesRing(QQ); R
    Power Series Ring in t over Rational Field
    sage: f = R.hom([t^2]); f
    Ring endomorphism of Power Series Ring in t over Rational Field
      Defn: t |--> t^2
    sage: R.set_default_prec(10)
    sage: s = 1/(1 + t); s
    1 - t + t^2 - t^3 + t^4 - t^5 + t^6 - t^7 + t^8 - t^9 + O(t^10)
    sage: f(s)
    1 - t^2 + t^4 - t^6 + t^8 - t^10 + t^12 - t^14 + t^16 - t^18 + O(t^20)

EXAMPLE: Frobenious on a power series ring over a finite field.
    sage: R.<t> = PowerSeriesRing(GF(5))
    sage: f = R.hom([t^5]); f
    Ring endomorphism of Power Series Ring in t over Finite Field of size 5
      Defn: t |--> t^5
    sage: a = 2 + t + 3*t^2 + 4*t^3 + O(t^4)
    sage: b = 1 + t + 2*t^2 + t^3 + O(t^5)
    sage: f(a)
    2 + t^5 + 3*t^10 + 4*t^15 + O(t^20)
    sage: f(b)
    1 + t^5 + 2*t^10 + t^15 + O(t^25)
    sage: f(a*b)
    2 + 3*t^5 + 3*t^10 + t^15 + O(t^20)
    sage: f(a)*f(b)
    2 + 3*t^5 + 3*t^10 + t^15 + O(t^20)

EXAMPLE: Homomorphism of Laurent series ring.
    sage: R.<t> = LaurentSeriesRing(QQ)
    sage: f = R.hom([t^3 + t]); f
    Ring endomorphism of Laurent Series Ring in t over Rational Field
      Defn: t |--> t + t^3
    sage: R.set_default_prec(10)
    sage: s = 2/t^2 + 1/(1 + t); s
    2*t^-2 + 1 - t + t^2 - t^3 + t^4 - t^5 + t^6 - t^7 + t^8 - t^9 + O(t^10)
    sage: f(s)
    2*t^-2 - 3 - t + 7*t^2 - 2*t^3 - 5*t^4 - 4*t^5 + 16*t^6 - 9*t^7 + O(t^8)
    sage: f = R.hom([t^3]); f
    Ring endomorphism of Laurent Series Ring in t over Rational Field
      Defn: t |--> t^3
    sage: f(s)
    2*t^-6 + 1 - t^3 + t^6 - t^9 + t^12 - t^15 + t^18 - t^21 + t^24 - t^27
    sage: s = 2/t^2 + 1/(1 + t); s
    2*t^-2 + 1 - t + t^2 - t^3 + t^4 - t^5 + t^6 - t^7 + t^8 - t^9 + O(t^10)
    sage: f(s)
    2*t^-6 + 1 - t^3 + t^6 - t^9 + t^12 - t^15 + t^18 - t^21 + t^24 - t^27

Note that the homomorphism must result in a converging Laurent series,
so the valuation of the image of the generator must be positive:
    sage: R.hom([1/t])
    Traceback (most recent call last):
    ...
    TypeError: images do not define a valid homomorphism
    sage: R.hom([1])
    Traceback (most recent call last):
    ...
    TypeError: images do not define a valid homomorphism


EXAMPLE: Complex conjugation on cyclotomic fields.
    sage: K.<zeta7> = CyclotomicField(7)
    sage: c = K.hom([1/zeta7]); c
    Ring endomorphism of Cyclotomic Field of order 7 and degree 6
      Defn: zeta7 |--> -zeta7^5 - zeta7^4 - zeta7^3 - zeta7^2 - zeta7 - 1
    sage: a = (1+zeta7)^5; a
    zeta7^5 + 5*zeta7^4 + 10*zeta7^3 + 10*zeta7^2 + 5*zeta7 + 1
    sage: c(a)
    5*zeta7^5 + 5*zeta7^4 - 4*zeta7^2 - 5*zeta7 - 4
    sage: c(zeta7 + 1/zeta7)       # this element is obviously fixed by inversion
    -zeta7^5 - zeta7^4 - zeta7^3 - zeta7^2 - 1
    sage: zeta7 + 1/zeta7
    -zeta7^5 - zeta7^4 - zeta7^3 - zeta7^2 - 1

EXAMPLE: Embedding a number field into the reals.
    sage: R.<x> = PolynomialRing(QQ)
    sage: K.<beta> = NumberField(x^3 - 2)
    sage: alpha = RR(2)^(1/3); alpha
    1.25992104989487
    sage: i = K.hom([alpha],check=False); i
    Ring morphism:
      From: Number Field in beta with defining polynomial x^3 - 2
      To:   Real Field with 53 bits of precision
      Defn: beta |--> 1.25992104989487
    sage: i(beta)
    1.25992104989487
    sage: i(beta^3)
    2.00000000000000
    sage: i(beta^2 + 1)
    2.58740105196820
"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.all import Morphism, is_Homset, Sets
import ideal

import homset

def is_RingHomomorphism(phi):
    return isinstance(phi, RingHomomorphism)

class RingMap(Morphism):
    """
    Set-theoretic map between rings.
    """
    def __init__(self, parent):
        Morphism.__init__(self, parent)

    def _repr_type(self):
        return "Set-theoretic ring"

    def __call__(self, x):
        if ideal.is_Ideal(x):
            R = self.codomain()
            return R.ideal([self(y) for y in x.gens()])
        return Morphism.__call__(self, x)


class RingMap_lift(RingMap):
    r"""
    Given rings $R$ and $S$ such that for any $x \in R$ the function
    \code{x.lift()} is an element that naturally coerces to $S$, this
    returns the set-theoretic ring map $R \to S$ sending $x$ to
    \code{x.lift()}.

    EXAMPLES:
        sage: R, (x,y) = PolynomialRing(QQ, 2, 'xy').objgens()
        sage: S.<xbar,ybar> = R.quo( (x^2 + y^2, y) )
        sage: S.lift()
        Set-theoretic ring morphism:
          From: Quotient of Polynomial Ring in x, y over Rational Field by the ideal (y, y^2 + x^2)
          To:   Polynomial Ring in x, y over Rational Field
          Defn: Choice of lifting map
    """
    def __init__(self, R, S):
        H = R.Hom(S, Sets())
        RingMap.__init__(self, H)
        self.__S = S  # for efficiency
        try:
            S._coerce_(R(0).lift())
        except TypeError:
            raise TypeError, "No natural lift map"

    def _repr_defn(self):
        return "Choice of lifting map"

    def _call_(self, x):
        return self.__S._coerce_(x.lift())

class RingHomomorphism(RingMap):
    """
    Homomorphism of rings.
    """
    def __init__(self, parent):
        if not homset.is_RingHomset(parent):
            raise TypeError, "parent must be a ring homset"
        RingMap.__init__(self, parent)

    def _repr_type(self):
        return "Ring"

    def _set_lift(self, lift):
        if not isinstance(lift, RingMap):
            raise TypeError, "lift must be a RingMap"
        if lift.domain() != self.codomain():
            raise TypeError, "lift must have correct domain"
        if lift.codomain() != self.domain():
            raise TypeError, "lift must have correct codomain"
        self.__lift = lift

    def is_injective(self):
        ## TODO -- actually implement this in some generality (!)
        raise NotImplementedError

    def is_zero(self):
        r"""
        Return True if this is the zero map and False otherwise.

        A *ring* homomorphism is considered to be 0 if and only if
        it sends the 1 element of the domain to the 0 element of the
        codomain.

        EXAMPLES:
        First an example of a map that is obviously nonzero.
            sage: h = Hom(ZZ, QQ)
            sage: f = h.natural_map()
            sage: f.is_zero()
            False

        Next we make the zero ring as $\ZZ/1\ZZ$.
            sage: R = Integers(1)
            sage: R
            Ring of integers modulo 1
            sage: h = Hom(ZZ, R)
            sage: f = h.natural_map()
            sage: f.is_zero()
            True

        Finally we check an example in characteristic 2.
            sage: h = Hom(ZZ, GF(2))
            sage: f = h.natural_map()
            sage: f.is_zero()
            False
        """
        return self(self.domain()(1)) == self.codomain()(0)

    def inverse_image(self, I):
        """
        Return the inverse image of the ideal $I$ under this ring
        homomorphism.
        """
        raise NotImplementedError

    def lift(self, x=None):
        """
        Return a lifting homomorphism associated to this homomorphism,
        if it has been defined.

        If x is not None, return the value of the lift morphism on x.
        """
        if not (x is None):
            return self.lift()(x)
        try:
            return self.__lift
        except AttributeError:
            raise ValueError, "No lift map defined."


from sage.categories.morphism import FormalCoercionMorphism
class RingHomomorphism_coercion(FormalCoercionMorphism, RingHomomorphism):
    pass

import sage.structure.all

class RingHomomorphism_im_gens(RingHomomorphism):
    """
    A ring homomorphism determined by the images of generators.
    """
    def __init__(self, parent, im_gens, check=True):
        RingHomomorphism.__init__(self, parent)
        if not isinstance(im_gens, (tuple, list)):
            im_gens = [im_gens]
        im_gens = sage.structure.all.Sequence(im_gens, parent.codomain())
        if len(im_gens) != parent.domain().ngens():
            raise ValueError, "number of images must equal number of generators"
        if check:
            t = parent.domain()._is_valid_homomorphism_(parent.codomain(), im_gens)
            if not t:
                raise ValueError, "relations do not all (canonically) map to 0 under map determined by images of generators."
        self.__im_gens = im_gens

    def im_gens(self):
        return self.__im_gens

    def _repr_defn(self):
        D = self.domain()
        ig = self.__im_gens
        return '\n'.join(['%s |--> %s'%(D.gen(i), ig[i]) for\
                       i in range(D.ngens())])

    def _call_(self, x):
        return x._im_gens_(self.codomain(), self.im_gens())

class RingHomomorphism_cover(RingHomomorphism):
    r"""
    A homomorphism induced by quotienting a ring out by an ideal.

    EXAMPLES:
        sage: R.<x,y> = PolynomialRing(QQ, 2)
        sage: S.<a,b> = R.quo(x^2 + y^2)
        sage: phi = S.cover(); phi
        Ring morphism:
          From: Polynomial Ring in x, y over Rational Field
          To:   Quotient of Polynomial Ring in x, y over Rational Field by the ideal (y^2 + x^2)
          Defn: Natural quotient map
        sage: phi(x+y)
        b + a
    """
    def __init__(self, ring, quotient_ring):
        RingHomomorphism.__init__(self, ring.Hom(quotient_ring))

    def _call_(self, x):
        return self.codomain()(x)

    def _repr_defn(self):
        return "Natural quotient map"

    def kernel(self):
        return self.codomain().defining_ideal()


class RingHomomorphism_from_quotient(RingHomomorphism):
    r"""
    A ring homomorphism with domain a generic quotient ring.

    INPUT:
        parent -- a ring homset Hom(R,S)
        phi -- a ring homomorphism C --> S, where C is the
               domain of R.cover()
    OUTPUT:
        a ring homomorphism

    The domain $R$ is a quotient object $C \to R$, and
    \code{R.cover()} is the ring homomorphism $\varphi: C \to R$.  The
    condition on the elements \code{im_gens} of $S$ is that they
    define a homomorphism $C \to S$ such that each generator of the
    kernel of $\varphi$ maps to $0$.

    EXAMPLES:
        sage: R.<x, y, z> = PolynomialRing(QQ, 3)
        sage: S.<a, b, c> = R.quo(x^3 + y^3 + z^3)
        sage: phi = S.hom([b, c, a]); phi
        Ring endomorphism of Quotient of Polynomial Ring in x, y, z over Rational Field by the ideal (z^3 + y^3 + x^3)
          Defn: a |--> b
                b |--> c
                c |--> a
        sage: phi(a+b+c)
        c + b + a


     Validity of the homomorphism is determined, when possible, and a
     TypeError is raised if there is no homomorphism sending the
     generators to the given images.
        sage: S.hom([b^2, c^2, a^2])
        Traceback (most recent call last):
        ...
        TypeError: images do not define a valid homomorphism
    """
    def __init__(self, parent, phi):
        RingHomomorphism.__init__(self, parent)
        R = parent.domain()
        pi = R.cover()  # the covering map
        if pi.domain() != phi.domain():
            raise ValueError, "Domain of phi must equal domain of covering."
        for x in pi.kernel().gens():
            if phi(x) != 0:
                raise ValueError, "relations do not all (canonically) map to 0 under map determined by images of generators."
        self.__lift = pi.lift()
        self.__phi = phi

    def morphism_from_cover(self):
        return self.__phi

    def _repr_defn(self):
        D = self.domain()
        ig = self.__phi.im_gens()
        return '\n'.join(['%s |--> %s'%(D.gen(i), ig[i]) for\
                          i in range(D.ngens())])

    def _call_(self, x):
        return self.__phi(self.__lift(x))




