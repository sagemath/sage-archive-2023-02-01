# -*- coding: utf-8 -*-
"""
Gauss valuations on polynomial rings

This file implements Gauss valuations for polynomial rings, i.e. discrete
valuations which assign to a polynomial the minimal valuation of its
coefficients.

AUTHORS:

- Julian Rüth (15-04-2013): initial version

EXAMPLES::

    sage: from mac_lane import * # optional: standalone
    sage: R.<x> = QQ[]
    sage: v0 = pAdicValuation(QQ, 2)
    sage: v = GaussValuation(R, v0); v
    Gauss valuation induced by 2-adic valuation
    sage: v(2*x + 2)
    1

Gauss valuations can also be defined iteratively based on valuations over
polynomial rings::

    sage: v = v.augmentation(x, 1/4); v
    [ Gauss valuation induced by 2-adic valuation, v(x) = 1/4 ]
    sage: v = v.augmentation(x^4+2*x^3+2*x^2+2*x+2, 4/3); v
    [ Gauss valuation induced by 2-adic valuation, v(x) = 1/4, v(x^4 + 2*x^3 + 2*x^2 + 2*x + 2) = 4/3 ]
    sage: S.<T> = R[]
    sage: w = GaussValuation(S, v); w
    Gauss valuation induced by [ Gauss valuation induced by 2-adic valuation, v(x) = 1/4, v(x^4 + 2*x^3 + 2*x^2 + 2*x + 2) = 4/3 ]
    sage: w(2*T + 1)
    0

"""
#*****************************************************************************
#       Copyright (C) 2013-2016 Julian Rüth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

# Fix doctests so they work in standalone mode (when invoked with sage -t, they run within the mac_lane/ directory)
import sys, os
if hasattr(sys.modules['__main__'], 'DC') and 'standalone' in sys.modules['__main__'].DC.options.optional:
    sys.path.append(os.getcwd())
    sys.path.append(os.path.dirname(os.getcwd()))

from inductive_valuation import FiniteInductiveValuation

from sage.misc.cachefunc import cached_method
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.factory import UniqueFactory

class GaussValuationFactory(UniqueFactory):
    r"""
    Create a Gauss valuation on ``domain``.

    INPUT:

    - ``domain`` -- a univariate polynomial ring

    - ``v`` -- a valuation on the base ring of ``domain``, the underlying
      valuation on the constants of the polynomial ring (if unspecified take
      the natural valuation on the valued ring ``domain``.)

    EXAMPLES:

    The Gauss valuation is the minimum of the valuation of the coefficients::

        sage: from mac_lane import * # optional: standalone
        sage: v = pAdicValuation(QQ, 2)
        sage: R.<x> = QQ[]
        sage: w = GaussValuation(R, v)
        sage: w(2)
        1
        sage: w(x)
        0
        sage: w(x + 2)
        0

    """
    def create_key(self, domain, v = None):
        r"""
        Normalize and check the parameters to create a Gauss valuation.

        TESTS::

            sage: from mac_lane import * # optional: standalone
            sage: v = pAdicValuation(QQ, 2)
            sage: R.<x> = ZZ[]
            sage: GaussValuation.create_key(R, v)
            Traceback (most recent call last):
            ...
            ValueError: the domain of v must be the base ring of domain but 2-adic valuation is not defined over Integer Ring but over Rational Field

        """
        from sage.rings.polynomial.polynomial_ring import is_PolynomialRing
        if not is_PolynomialRing(domain):
            raise TypeError("GaussValuations can only be created over polynomial rings but %r is not a polynomial ring"%(domain,))
        if not domain.ngens() == 1:
            raise NotImplementedError("domain must be univariate but %r is not univariate"%(domain,))

        if v is None:
            v = domain.base_ring().valuation()

        if not v.domain() is domain.base_ring():
            raise ValueError("the domain of v must be the base ring of domain but %r is not defined over %r but over %r"%(v, domain.base_ring(), v.domain()))
        if not v.is_discrete_valuation():
            raise ValueError("v must be a discrete valuation but %r is not"%(v,))

        return (domain, v)

    def create_object(self, version, key, **extra_args):
        r"""
        Create a Gauss valuation from normalized parameters.

        TESTS::

            sage: from mac_lane import * # optional: standalone
            sage: v = pAdicValuation(QQ, 2)
            sage: R.<x> = QQ[]
            sage: GaussValuation.create_object(0, (R, v))
            Gauss valuation induced by 2-adic valuation

        """
        domain, v = key
        from sage.rings.valuation.valuation_space import DiscretePseudoValuationSpace
        parent = DiscretePseudoValuationSpace(domain)
        return parent.__make_element_class__(GaussValuation_generic)(parent, v)

GaussValuation = GaussValuationFactory("GaussValuation")

class GaussValuation_generic(FiniteInductiveValuation):
    """
    A Gauss valuation on a polynomial ring ``domain``.

    INPUT:

    - ``domain`` -- a univariate polynomial ring over a valued ring `R`

    - ``v`` -- a discrete valuation on `R`

    EXAMPLES::

        sage: from mac_lane import * # optional: standalone
        sage: R = Zp(3,5)
        sage: S.<x> = R[]
        sage: v0 = pAdicValuation(R)
        sage: v = GaussValuation(S, v0); v
        Gauss valuation induced by 3-adic valuation

        sage: S.<x> = QQ[]
        sage: v = GaussValuation(S, pAdicValuation(QQ, 5)); v
        Gauss valuation induced by 5-adic valuation

    TESTS::

        sage: TestSuite(v).run()

    """
    def __init__(self, parent, v):
        """
        TESTS::

            sage: from mac_lane import * # optional: standalone
            sage: from mac_lane.gauss_valuation import GaussValuation_generic # optional: standalone
            sage: S.<x> = QQ[]
            sage: v = GaussValuation(S, pAdicValuation(QQ, 5))
            sage: isinstance(v, GaussValuation_generic)
            True

        """
        FiniteInductiveValuation.__init__(self, parent, parent.domain().gen())

        self._base_valuation = v

    def value_group(self):
        """
        Return the value group of this valuation.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: S.<x> = QQ[]
            sage: v = GaussValuation(S, pAdicValuation(QQ, 5))
            sage: v.value_group()
            Additive Abelian Group generated by 1

        """
        return self._base_valuation.value_group()

    def _repr_(self):
        """
        Return a printable representation of this valuation.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: S.<x> = QQ[]
            sage: v = GaussValuation(S, pAdicValuation(QQ, 5))
            sage: v # indirect doctest
            Gauss valuation induced by 5-adic valuation

        """
        return "Gauss valuation induced by %r"%self._base_valuation

    @cached_method
    def uniformizer(self):
        """
        Return a uniformizer of this valuation, i.e., a uniformizer of the
        valuation of the base ring.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: S.<x> = QQ[]
            sage: v = GaussValuation(S, pAdicValuation(QQ, 5))
            sage: v.uniformizer()
            5
            sage: v.uniformizer().parent() is S
            True

        """
        return self.domain()(self._base_valuation.uniformizer())

    def shift(self, f, s):
        """
        Multiply ``f`` by the ``s``th power of the uniformizer of this
        valuation.

        INPUT:

        - ``f`` -- a polynomial in the domain of this valuation

        - ``s`` -- an element of the :meth:`value_group`

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: S.<x> = QQ[]
            sage: v = GaussValuation(S, pAdicValuation(QQ, 5))
            sage: v.shift(x, -2)
            1/25*x

        It is an error to perform a shift if the result is not in the domain of
        the valuation anymore::

            sage: S.<x> = Zp(2,5)[]
            sage: v = GaussValuation(S)
            sage: f = v.shift(x, 2); f
            (2^2 + O(2^7))*x
            sage: f.parent() is S
            True
            sage: f = v.shift(x, -2)
            Traceback (most recent call last):
            ...
            ValueError: since 2-adic Ring with capped relative precision 5 is not a field, -s must not exceed the valuation of f but 2 does exceed 0

        Of course, the above example works over a field::

            sage: S.<x> = Qp(2,5)[]
            sage: v = GaussValuation(S)
            sage: f = v.shift(x, -2); f
            (2^-2 + O(2^3))*x

        """
        f = self.domain().coerce(f)
        s = self.value_group()(s)

        from sage.categories.fields import Fields
        if -s > self(f) and self.domain().base_ring() not in Fields():
            raise ValueError("since %r is not a field, -s must not exceed the valuation of f but %r does exceed %r"%(self.domain().base_ring(), -s, self(f)))

        return f.map_coefficients(lambda c:self._base_valuation.shift(c, s))

    # TODO: declare this upstairs
    def valuations(self, f):
        """
        Return the valuations of the `f_i\phi^i` in the expansion `f=\sum f_i\phi^i`.

        INPUT:

        - ``f`` -- a polynomial in the domain of this valuation

        OUTPUT:

        A list of rational numbers, the valuations of `f_0, f_1\phi, \dots`

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R = Qp(2,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S, pAdicValuation(R))
            sage: f = x^2 + 2*x + 16
            sage: list(v.valuations(f))
            [4, 1, 0]

        """
        f = self.domain().coerce(f)

        for c in self.coefficients(f):
            yield self._base_valuation(self.domain().base_ring()(c))

    @cached_method
    def residue_ring(self):
        """
        Return the residue ring of this valuation, i.e., the elements of
        valuation zero module the elements of positive valuation.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: S.<x> = Qp(2,5)[]
            sage: v = GaussValuation(S)
            sage: v.residue_ring()
            Univariate Polynomial Ring in x over Finite Field of size 2 (using NTL)

        """
        return self.domain().change_ring(self._base_valuation.residue_ring())

    def reduce(self, f):
        """
        Return the reduction of ``f`` modulo this valuation.

        INPUT:

        - ``f`` -- an integral element of the domain of this valuation

        OUTPUT:

        A polynomial in the :meth:`residue_ring` of this valuation.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: S.<x> = Qp(2,5)[]
            sage: v = GaussValuation(S)
            sage: f = x^2 + 2*x + 16
            sage: v.reduce(f)
            x^2
            sage: v.reduce(f).parent() is v.residue_ring()
            True

        The reduction is only defined for integral elements::

            sage: f = x^2/2
            sage: v.reduce(f)
            Traceback (most recent call last):
            ...
            ValueError: reduction not defined for non-integral elements and (2^-1 + O(2^4))*x^2 is not integral over Gauss valuation induced by 2-adic valuation

        .. SEEALSO::

            :meth: `lift`

        """
        f = self.domain().coerce(f)
        if not all([v>=0 for v in self.valuations(f)]):
            raise ValueError("reduction not defined for non-integral elements and %r is not integral over %r"%(f, self))

        return f.map_coefficients(lambda c:self._base_valuation.reduce(c), self.restriction(self.domain().base_ring()).residue_field())

    def lift(self, F):
        """
        Return a lift of ``F``.

        INPUT::

        - ``F`` -- a polynomial over the :meth:`residue_ring` of this valuation

        OUTPUT:

        a (possibly non-monic) polynomial in the domain of this valuation which
        reduces to ``F``

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: S.<x> = Qp(3,5)[]
            sage: v = GaussValuation(S)
            sage: f = x^2 + 2*x + 16
            sage: F = v.reduce(f); F
            x^2 + 2*x + 1
            sage: g = v.lift(F); g
            (1 + O(3^5))*x^2 + (2 + O(3^5))*x + (1 + O(3^5))
            sage: v.is_equivalent(f,g)
            True
            sage: g.parent() is v.domain()
            True

        .. SEEALSO::

            :meth:`reduce`

        """
        F = self.residue_ring().coerce(F)

        return F.map_coefficients(lambda c:self._base_valuation.lift(c), self._base_valuation.domain())

    # TODO: declare this upstairs
    def lift_to_key(self, F):
        """
        Lift the irreducible polynomial ``F`` to a key polynomial.

        INPUT:

        - ``F`` -- an irreducible non-constant monic polynomial in
          :meth:`residue_ring` of this valuation

        OUTPUT:

        A polynomial `f` in the domain of this valuation which is a key
        polynomial for this valuation and which, for a suitable equivalence
        unit `R`, satifies that the reduction of `Rf` is ``F``

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R.<u> = Qq(4,10)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: y = v.residue_ring().gen()
            sage: u0 = v.residue_ring().base_ring().gen()
            sage: f = v.lift_to_key(y^2 + y + u0); f
            (1 + O(2^10))*x^2 + (1 + O(2^10))*x + u + O(2^10)

        """
        F = self.residue_ring().coerce(F)

        if F.is_constant():
            raise ValueError("F must not be constant but %r is constant"%(F,))
        if not F.is_monic():
            raise ValueError("F must be monic but %r is not monic"%(F,))
        if not F.is_irreducible():
            raise ValueError("F must be irreducible but %r factors"%(F,))

        return self.lift(F)

    # TODO: declare this upstairs
    def equivalence_unit(self, s):
        """
        Return an equivalence unit of valuation ``s``.

        INPUT:

        - ``s`` -- an element of the :meth:`value_group`

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: S.<x> = Qp(3,5)[]
            sage: v = GaussValuation(S)
            sage: v.equivalence_unit(2)
            (3^2 + O(3^7))
            sage: v.equivalence_unit(-2)
            (3^-2 + O(3^3))

        """
        one = self._base_valuation.domain().one()
        ret = self._base_valuation.shift(one, s)
        return self.domain()(ret)

    # TODO: eliminate this
    element_with_valuation = equivalence_unit

    # TODO: declare this upstairs
    def is_commensurable_inductive(self):
        """
        Return whether this valuation is a commensurable inductive valuation
        over the discrete valuation of the base ring of the polynomial ring, as
        defined in section 4 of [ML1936].

        OUTPUT:

        ``True`` since a Gauss valuation always is commensurable inductive.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.is_commensurable_inductive()
            True

        REFERENCES:

        .. [ML1936] Mac Lane, S. (1936). A construction for prime ideals as absolute
        values of an algebraic field. Duke Mathematical Journal, 2(3), 492-510.

        """
        return True

    # TODO: declare this uptstairs
    def E(self):
        """
        Return the ramification index of this valuation over its underlying
        Gauss valuation, i.e., 1.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.E()
            1

        """
        from sage.rings.all import ZZ
        return ZZ.one()

    # TODO: declare this upstairs
    def F(self):
        """
        Return the degree of the residue field extension of this valuation
        over the Gauss valuation, i.e., 1.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.F()
            1

        """
        from sage.rings.all import ZZ
        return ZZ.one()

    def change_domain(self, ring):
        r"""
        Return this valuation as a valuation over ``ring``.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: v = pAdicValuation(ZZ, 2)
            sage: R.<x> = ZZ[]
            sage: w = GaussValuation(R, v)
            sage: w.change_domain(QQ['x'])
            Gauss valuation induced by 2-adic valuation

        """
        from sage.rings.polynomial.polynomial_ring import is_PolynomialRing
        if is_PolynomialRing(ring) and ring.ngens() == 1:
            base_valuation = self._base_valuation.change_domain(ring.base_ring())
            return GaussValuation(self.domain().change_ring(ring.base_ring()), base_valuation)
        return super(GaussValuation_generic, self).change_domain(ring)

    def extensions(self, ring):
        r"""
        Return the extensions of this valuation to ``ring``.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: v = pAdicValuation(ZZ, 2)
            sage: R.<x> = ZZ[]
            sage: w = GaussValuation(R, v)
            sage: w.extensions(GaussianIntegers()['x'])
            [Gauss valuation induced by 2-adic valuation]

        """
        from sage.rings.polynomial.polynomial_ring import is_PolynomialRing
        if is_PolynomialRing(ring) and ring.ngens() == 1:
            if self.domain().is_subring(ring):
                return [GaussValuation(ring, w) for w in self._base_valuation.extensions(ring.base_ring())]
        return super(GaussValuation_generic, self).extensions(ring)

    def restriction(self, ring):
        r"""
        Return the restriction of this valuation to ``ring``.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: v = pAdicValuation(ZZ, 2)
            sage: R.<x> = ZZ[]
            sage: w = GaussValuation(R, v)
            sage: w.restriction(ZZ)
            2-adic valuation

        """
        if ring is self.domain().base_ring():
            return self._base_valuation
        return super(GaussValuation_generic, self).restriction(ring)

    # TODO: declare this upstairs
    def is_gauss_valuation(self):
        r"""
        Return whether this valuation is a Gauss valuation.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.is_gauss_valuation()
            True

        """
        return True

    # TODO: declare this upstairs under a better name
    def _augmentations(self):
        r"""
        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v._augmentations()
            [Gauss valuation induced by 2-adic valuation]

        """
        return [self]

    def is_trivial(self):
        r"""
        Return whether this is a trivial valuation (sending everything but zero
        to zero.)

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, TrivialValuation(QQ))
            sage: v.is_trivial()
            True

        """
        return self._base_valuation.is_trivial()

    # TODO: declare this upstairs under a better name
    def _make_monic_integral(self, G):
        r"""
        Return a polynomial ``G`` which defines the self extension of the base
        ring of the domain of this valuation but which is monic and integral.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, pAdicValuation(QQ, 2))
            sage: v._make_monic_integral(5*x^2 + 1/2*x + 1/4)
            x^2 + 1/5*x + 1/5

        """
        if not G.is_monic():
            # this might fail if the base ring is not a field
            G = G / G.leading_coefficient()
        while self(G) < 0:
            u = self._base_valuation.uniformizer()
            x = G.parent().gen()
            # this might fail if the base ring is not a field
            G = G.parent(G(x/u) * (u ** G.degree()))
        assert G.is_monic()
        return G
            
    def _ge_(self, other):
        r"""
        Return whether this valuation is greater than or equal to ``other``
        everywhere.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, pAdicValuation(QQ, 2))
            sage: w = GaussValuation(R, pAdicValuation(QQ, 3))
            sage: v >= w
            False
            sage: w >= v
            False

        """
        if isinstance(other, GaussValuation_generic):
            return self._base_valuation >= other._base_valuation
        from augmented_valuation import AugmentedValuation_base
        if isinstance(other, AugmentedValuation_base):
            return False
        if other.is_trivial():
            return other.is_discrete_valuation()
        return super(GaussValuation_generic, self)._ge_(other)

    def is_discrete_valuation(self):
        r"""
        Return whether this is a discrete valuation.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, TrivialValuation(QQ))
            sage: v.is_discrete_valuation()
            True

        """
        return True
