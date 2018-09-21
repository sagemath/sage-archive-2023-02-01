r"""
Derivations

Let `A` be a ring and `B` be an algebra over `A`.
A derivation `d : A \to B` is an additive map that satisfies
the Leibniz rule

.. MATH::

    d(xy) = x d(y) + d(x) y.

If you are given in addition a ring homomorphism `\theta : A \to B`,
a twisted derivation w.r.t. `\theta` (or a `\theta`-derivation) is
an additive map `d : A \to B` such that

.. MATH::

    d(xy) = \theta(x) d(y) + d(x) y.

When `\theta` is the morphism defining the structure of `A`-algebra
on `B`, a `\theta`-derivation is nothing but a derivation.
One easily checks that `\theta - id` is a `\theta`-derivation.

The set of derivations (resp. `\theta`-derivations) is a module 
over `B`.


This file provides support for derivations and twisted derivations
over commutative rings.

Given a ring `A`, the module of derivations over `A` can be created
as follows::

    sage: A.<x,y,z> = QQ[]
    sage: M = A.derivation_module()
    sage: M
    Module of derivations over Multivariate Polynomial Ring in x, y, z over Rational Field

A codomain can be specified::

    sage: B = A.fraction_field()
    sage: A.derivation_module(B)
    Module of derivations from Multivariate Polynomial Ring in x, y, z over Rational Field to Fraction Field of Multivariate Polynomial Ring in x, y, z over Rational Field

The method :meth:`gens` return generators of these modules::

    sage: M.gens()
    (d/dx, d/dy, d/dz)

We can combine them in order to create all derivations::

    sage: d = 2*M.gen(0) + z*M.gen(1) + (x^2 + y^2)*M.gen(2)
    sage: d
    2*d/dx + z*d/dy + (x^2 + y^2)*d/dz

and now play with them::

    sage: d(x + y + z)
    x^2 + y^2 + z + 2
    sage: P = A.random_element()
    sage: Q = A.random_element()
    sage: d(P*Q) == P*d(Q) + d(P)*Q
    True

Alternatively we can use the method :meth:`derivation` of the ring `A`
to create derivations::

    sage: A.derivation(x)
    d/dx
    sage: A.derivation(y)
    d/dy
    sage: A.derivation(z)
    d/dz
    sage: A.derivation([2, z, x^2+y^2])
    2*d/dx + z*d/dy + (x^2 + y^2)*d/dz

::

Twisted derivations and handled similarly::

    sage: theta = B.hom([B(y),B(z),B(x)])
    sage: theta
    Ring endomorphism of Fraction Field of Multivariate Polynomial Ring in x, y, z over Rational Field
      Defn: x |--> y
            y |--> z
            z |--> x

    sage: M = B.derivation_module(twist=theta)
    sage: M
    Module of twisted derivations over Fraction Field of Multivariate Polynomial Ring in x, y, z over Rational Field (twisting morphism: x |--> y, y |--> z, z |--> x)

Over a field, one proves that every `\theta`-derivation is a multiple
of `\theta - id`, so that::

    sage: d = M.gen(); d
    [x |--> y, y |--> z, z |--> x] - id

and then::

    sage: d(x)
    -x + y
    sage: d(y)
    -y + z
    sage: d(z)
    x - z
    sage: d(x + y + z)
    0

AUTHOR:

- Xavier Caruso (2018-09)
"""

#############################################################################
#    Copyright (C) 2018 Xavier Caruso <xavier.caruso@normalesup.org>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#                  http://www.gnu.org/licenses/
#****************************************************************************

from sage.structure.unique_representation import UniqueRepresentation
from sage.modules.module import Module
from sage.structure.element import ModuleElement
from sage.rings.integer_ring import ZZ

from sage.categories.map import Map
from sage.categories.all import Rings, Algebras
from sage.rings.morphism import RingMap, RingHomomorphism


class RingDerivationModule(Module, UniqueRepresentation):
    """
    A class for modules of derivations over a commutative ring.
    """
    def __init__(self, domain, codomain, twist=None, element_class=None):
        """
        Initialize this module of derivation.

        TESTS::

            sage: from sage.rings.derivation import RingDerivationModule
            sage: R5.<x> = GF(5)[]
            sage: R25.<x> = GF(25)[]
            sage: R7.<x> = GF(7)[]

            sage: RingDerivationModule(R5, R25)
            Module of derivations from Univariate Polynomial Ring in x over Finite Field of size 5 to Univariate Polynomial Ring in x over Finite Field in z2 of size 5^2
            sage: RingDerivationModule(R5, R5^2)
            Traceback (most recent call last):
            ...
            TypeError: The codomain must be an algebra over the domain
            sage: RingDerivationModule(R5, R7)
            Traceback (most recent call last):
            ...
            TypeError: The codomain must be an algebra over the domain

            sage: theta = R5.hom([R5.gen()^2])
            sage: RingDerivationModule(R5, R25, twist=theta)
            Module of twisted derivations from Univariate Polynomial Ring in x over Finite Field of size 5 to Univariate Polynomial Ring in x over Finite Field in z2 of size 5^2 (twisting morphism: x |--> x^2)
            sage: RingDerivationModule(R7, R7, twist=theta)
            Traceback (most recent call last):
            ...
            TypeError: The domain of the derivation must coerce to the domain of the twisting homomorphism

        """
        if not domain in Rings().Commutative():
            raise TypeError("The domain must be a commutative ring")
        if not (codomain in Rings().Commutative() and codomain.has_coerce_map_from(domain)):
            raise TypeError("The codomain must be an algebra over the domain")
        if twist is not None:
            if not (isinstance(twist, Map) and twist.category_for().is_subcategory(Rings())):
                raise TypeError("The twisting homorphism must be an homomorphism of rings")
            if twist.domain() is not domain:
                map = twist.domain().coerce_map_from(domain)
                if map is None:
                    raise TypeError("The domain of the derivation must coerce to the domain of the twisting homomorphism")
                twist = twist * map
            if twist.codomain() is not codomain:
                map = codomain.coerce_map_from(twist.codomain())
                if map is None:
                    raise TypeError("The codomain of the twisting homomorphism must coerce to the codomain of the derivation")
                twist = map * twist
            # We check if the twisting morphism is the identity
            try:
                if twist.is_identity():
                    twist = None
                else:
                    for g in domain.gens():
                        if self._twist(g) != g:
                            break
                    else:
                        twist = None
            except (AttributeError, NotImplementedError):
                pass
        self._domain = domain
        self._codomain = codomain
        self._twist = twist
        if element_class is None:
            if twist is None:
                self.element_class = RingDerivationWithoutTwist_im_gens
            else:
                self.element_class = RingDerivationWithTwist_generic
        else:
            self.element_class = element_class
        Module.__init__(self, codomain)

    def __hash__(self):
        """
        Returns a hash of this parent.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: M = R.derivation_module()
            sage: hash(M)  # random
            2727832899085333035

        """
        return hash((self._domain, self._codomain, self._twist))

    def _repr_(self):
        """
        Returns a string representation of this module of derivations.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: R.derivation_module()
            Module of derivations over Multivariate Polynomial Ring in x, y over Integer Ring

        ::

            sage: theta = R.hom([y,x])
            sage: R.derivation_module(twist=theta)
            Module of twisted derivations over Multivariate Polynomial Ring in x, y over Integer Ring (twisting morphism: x |--> y, y |--> x)

        """
        t = ""
        if self._twist is None:
            s = "Module of derivations"
        else:
            s = "Module of twisted derivations"
            try:
                t = " (twisting morphism: %s)" % self._twist._repr_short()
            except AttributeError:
                pass
        if self._domain is self._codomain:
            s += " over %s" % self._domain
        else:
            s += " from %s to %s" % (self._domain, self._codomain)
        return s + t

    def domain(self):
        """
        Returns the domain of the derivations in this module.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: M = R.derivation_module(); M
            Module of derivations over Multivariate Polynomial Ring in x, y over Integer Ring
            sage: M.domain()
            Multivariate Polynomial Ring in x, y over Integer Ring
        """
        return self._domain

    def codomain(self):
        """
        Returns the codomain of the derivations in this module.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: M = R.derivation_module(); M
            Module of derivations over Multivariate Polynomial Ring in x, y over Integer Ring
            sage: M.codomain()
            Multivariate Polynomial Ring in x, y over Integer Ring
        """
        return self._codomain

    def twisting_homomorphism(self):
        """
        Returns the twisting homorphism of the derivations in this module.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: theta = R.hom([y,x])
            sage: M = R.derivation_module(twist=theta); M
            Module of twisted derivations over Multivariate Polynomial Ring in x, y over Integer Ring (twisting morphism: x |--> y, y |--> x)
            sage: M.twisting_homomorphism()
            Ring endomorphism of Multivariate Polynomial Ring in x, y over Integer Ring
              Defn: x |--> y
                    y |--> x
        """
        if self._twist is None:
            return self._codomain.coerce_map_from(self._domain)
        else:
            return self._twist

    def ngens(self):
        """
        Returns the number of generators of this module of derivations.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: M = R.derivation_module(); M
            Module of derivations over Multivariate Polynomial Ring in x, y over Integer Ring
            sage: M.ngens()
            2

        Indeed, generators are::

            sage: M.gens()
            (d/dx, d/dy)

        We check that, For a nontrivial twist over a field, the module of 
        twisted derivation is a vector space of dimension 1 generated by 
        ``twist - id``::

            sage: K = R.fraction_field()
            sage: theta = K.hom([K(y),K(x)])
            sage: M = K.derivation_module(twist=theta); M
            Module of twisted derivations over Fraction Field of Multivariate Polynomial Ring in x, y over Integer Ring (twisting morphism: x |--> y, y |--> x)
            sage: M.ngens()
            1
            sage: M.gen()
            [x |--> y, y |--> x] - id

        """
        domain = self.domain()
        if self._twist is None:
            return domain.ngens()
        if self._codomain.is_field():
            return ZZ(1)
        raise NotImplementedError("Generators are not implemented for twisted derivations over rings")

    def gens(self):
        """
        Returns the generators of this module of derivations.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: M = R.derivation_module(); M
            Module of derivations over Multivariate Polynomial Ring in x, y over Integer Ring
            sage: M.gens()
            (d/dx, d/dy)

        We check that, For a nontrivial twist over a field, the module of 
        twisted derivation is a vector space of dimension 1 generated by 
        ``twist - id``::

            sage: K = R.fraction_field()
            sage: theta = K.hom([K(y),K(x)])
            sage: M = K.derivation_module(twist=theta); M
            Module of twisted derivations over Fraction Field of Multivariate Polynomial Ring in x, y over Integer Ring (twisting morphism: x |--> y, y |--> x)
            sage: M.gens()
            ([x |--> y, y |--> x] - id,)

        """
        domain = self.domain()
        if self._twist is None:
            c = [0] * domain.ngens()
            gens = []
            for i in range(domain.ngens()):
                c[i] = 1
                gens.append(self(c))
                c[i] = 0
            return tuple(gens)
        if self._codomain.is_field():
            return (self(1),)
        raise NotImplementedError("Generators are not implemented for twisted derivations over rings")

    def gen(self, n=0):
        """
        Returns ``n``th generator of this module of derivations.

        INPUT::

        ``n`` - an integer (default: ``0``)

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: M = R.derivation_module(); M
            Module of derivations over Multivariate Polynomial Ring in x, y over Integer Ring
            sage: M.gen()
            d/dx
            sage: M.gen(1)
            d/dy
        """
        domain = self.domain()
        if self._twist is None:
            c = [0] * domain.ngens()
            try:
                c[n] = 1
            except IndexError:
                raise ValueError("Generator not defined")
            return self(c)
        if self._codomain.is_field():
            return self(1)
        raise NotImplementedError("Generators are not implemented for twisted derivations over rings")


# The class RingDerivation does not derive from Map (or RingMap)
# because we don't want to see derivations as morphisms in some
# category since they are not stable by composition.
class RingDerivation(ModuleElement):
    """
    An abstract class for twisted and untwisted derivations over 
    commutative rings.

    TESTS::

        sage: R.<x,y> = ZZ[]
        sage: f = R.derivation(x) + 2*R.derivation(y); f
        d/dx + 2*d/dy
        sage: f(x*y)
        2*x + y

    """
    def __call__(self, x):
        """
        Returns the image of ``x`` under this derivation.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: f = x*R.derivation(x) + y*R.derivation(y)
            sage: f(x^2 + 3*x*y - y^2)
            2*x^2 + 6*x*y - 2*y^2
        """
        arg = self.parent().domain()(x)
        return self._call_(arg)


class RingDerivationWithoutTwist_im_gens(RingDerivation):
    """
    A class for untwisted derivations over rings with generators.
    """
    def __init__(self, parent, images):
        """
        Initialize this derivation.

        TESTS::

            sage: R.<x,y> = ZZ[]
            sage: R.derivation(x)  # indirect doctest
            d/dx
            sage: R.derivation([1,2])  # indirect doctest
            d/dx + 2*d/dy

        """
        domain = parent.domain()
        codomain = parent.codomain()
        self._images = [ codomain(x) for x in images ]
        if len(self._images) != domain.ngens():
            raise ValueError("The number of images is incorrect")
        RingDerivation.__init__(self, parent)

    def __hash__(self):
        """
        Returns a hash of this derivation.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: f = R.derivation(x)
            sage: hash(f)  # random
            -5740389279366308785

        """
        parent = self.parent()
        return hash((parent._domain, parent._codomain, tuple(self._images)))

    def _repr_(self):
        """
        Returns a string representation of this derivation.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: R.derivation(x)  # indirect doctest
            d/dx
            sage: R.derivation(y)  # indirect doctest
            d/dy
        """
        s = ""
        domain = self.parent().domain()
        for i in range(len(self._images)):
            c = self._images[i]
            sc = str(c)
            if sc == "0":
                continue
            ddx = "d/d%s" % domain.gen(i)
            if sc == "1":
                s += " + " + ddx
            elif sc == "-1":
                s += " - " + ddx
            elif c._is_atomic():
                s += " + %s*%s" % (sc, ddx)
            elif (-c)._is_atomic():
                s += " - %s*%s" % (-c, ddx)
            else:
                s += " + (%s)*%s" % (sc, ddx)
        if s[:3] == " + ":
            return s[3:]
        elif s[:3] == " - ":
            return "-" + s[3:]
        elif s == "":
            return "0"
        else:
            return s

    def _add_(self, other):
        """
        Returns the sum of this derivation and ``other``.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: derx = R.derivation(x)
            sage: dery = R.derivation(y)
            sage: derx + dery  # indirect doctest
            d/dx + d/dy
        """
        im = [ self._images[i] + other._images[i] for i in range(self.parent().domain().ngens()) ]
        return self.parent()(im)

    def _sub_(self, other):
        """
        Returns the subtraction of this derivation and ``other``.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: derx = R.derivation(x)
            sage: dery = R.derivation(y)
            sage: derx - dery  # indirect doctest
            d/dx - d/dy
        """
        im = [ self._images[i] - other._images[i] for i in range(self.parent().domain().ngens()) ]
        return self.parent()(im)

    def _rmul_(self, factor):
        """
        Returns the product of this derivation by the scalar ``factor``.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: derx = R.derivation(x)
            sage: 2 * derx
            2*d/dx
            sage: x^2 * derx
            x^2*d/dx
        """
        factor = self.parent().codomain()(factor)
        im = [ factor*x  for x in self._images ]
        return self.parent()(im)

    def _lmul_(self, factor):
        """
        Returns the product of this derivation by the scalar ``factor``.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: derx = R.derivation(x)
            sage: derx * 2
            2*d/dx
            sage: derx * x^2
            x^2*d/dx
        """
        return self._rmul_(factor)

    def _call_(self, x):
        """
        Returns the image of ``x`` under this derivation.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: f = x*R.derivation(x) + y*R.derivation(y)
            sage: f(x^2 + 3*x*y - y^2)
            2*x^2 + 6*x*y - 2*y^2
        """
        res = self.parent().codomain()(0)
        domain = self.parent().domain()
        for i in range(len(self._images)):
            res += x.derivative(domain.gen(i)) * self._images[i]
        return res

    def is_zero(self):
        """
        Return ``True`` if this derivation is zero.

        EXEMPLES::

            sage: R.<x,y> = ZZ[]
            sage: f = R.derivation(); f
            d/dx
            sage: f.is_zero()
            False

            sage: (f-f).is_zero()
            True
        """
        for im in self._images:
            if im != 0: return False
        return True


class RingDerivationWithTwist_generic(RingDerivation):
    r"""
    The class handles `\theta`-derivations of the form
    `\lambda*(\theta - id)` for a scalar `\lambda` varying
    in the codomain of `\theta`.
    """
    def __init__(self, parent, scalar=0):
        """
        Initialize this derivation.

        TESTS::

            sage: R.<x,y> = ZZ[]
            sage: theta = R.hom([y,x])
            sage: R.derivation_twisted(theta)  # indirect doctest
            0
            sage: R.derivation_twisted(theta, 1)  # indirect doctest
            [x |--> y, y |--> x] - id
            sage: R.derivation_twisted(theta, x)  # indirect doctest
            x*([x |--> y, y |--> x] - id)

        """
        codomain = parent.codomain()
        self._scalar = codomain(scalar)
        RingDerivation.__init__(self, parent)

    def __hash__(self):
        """
        Returns a hash of this derivation.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: theta = R.hom([y,x])
            sage: f = R.derivation_twisted(theta, 1)
            sage: hash(f)  # random
            -6511057926760520014
        """
        parent = self.parent()
        return hash((parent._domain, parent._codomain, self._scalar))

    def _repr_(self):
        """
        Returns a string representation of this derivation.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: theta = R.hom([y,x])
            sage: R.derivation_twisted(theta, 1)
            [x |--> y, y |--> x] - id
        """
        scalar = self._scalar
        sc = str(scalar)
        if sc == "0":
            return "0"
        try:
            t = "[%s] - id" % self.parent().twisting_homomorphism()._repr_short();
        except AttributeError:
            t = "twisting_morphism - id"
        if sc == "1":
            return t
        elif sc == "-1":
            s = "-"
        elif scalar._is_atomic():
            s = "%s*" % sc
        elif (-scalar)._is_atomic():
            s = "-%s*" % (-c)
        else:
            s = "(%s)*" % sc
        return "%s(%s)" % (s,t)

    def _add_(self, other):
        """
        Returns the sum of this derivation and ``other``.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: theta = R.hom([y,x])
            sage: der1 = R.derivation_twisted(theta, x); der1
            x*([x |--> y, y |--> x] - id)
            sage: der2 = R.derivation_twisted(theta, y); der2
            y*([x |--> y, y |--> x] - id)
            sage: der1 + der2
            (x + y)*([x |--> y, y |--> x] - id)
        """
        return self.parent()(self._scalar + other._scalar)

    def _sub_(self, other):
        """
        Returns the subtraction of this derivation and ``other``.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: theta = R.hom([y,x])
            sage: der1 = R.derivation_twisted(theta, x); der1
            x*([x |--> y, y |--> x] - id)
            sage: der2 = R.derivation_twisted(theta, y); der2
            y*([x |--> y, y |--> x] - id)
            sage: der1 - der2
            (x - y)*([x |--> y, y |--> x] - id)

        TESTS::

            sage: der1 - der1
            0
            sage: der2 - der2
            0
        """
        return self.parent()(self._scalar - other._scalar)

    def _rmul_(self, factor):
        """
        Returns the product of this derivation by the scalar ``factor``.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: theta = R.hom([y,x])
            sage: der1 = R.derivation_twisted(theta, x); der1
            x*([x |--> y, y |--> x] - id)
            sage: y * der1
            x*y*([x |--> y, y |--> x] - id)
        """
        return self.parent()(factor * self._scalar)

    def _lmul_(self, factor):
        """
        Returns the product of this derivation by the scalar ``factor``.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: theta = R.hom([y,x])
            sage: der1 = R.derivation_twisted(theta, x); der1
            x*([x |--> y, y |--> x] - id)
            sage: der1 * y
            x*y*([x |--> y, y |--> x] - id)
        """
        return self._rmul_(factor)

    def _call_(self, x):
        """
        Returns the image of ``x`` under this derivation.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: theta = R.hom([y,x])
            sage: f = R.derivation_twisted(theta, 1); f
            [x |--> y, y |--> x] - id
            sage: f(x)
            -x + y
        """
        return self._scalar * (self.parent().twisting_homomorphism()(x) - x)
