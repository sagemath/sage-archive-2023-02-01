# -*- coding: utf-8 -*-
r"""
Augmented valuations on polynomial rings

Implements augmentations of valutions as defined in [ML1936].

TESTS::

    sage: from mac_lane import * # optional: standalone
    sage: R.<x> = QQ[]
    sage: v = GaussValuation(R, pAdicValuation(QQ, 2))
    sage: w = v.augmentation(x, 1)
    sage: TestSuite(w).run() # long time

    sage: w = v.augmentation(x, 2)
    sage: TestSuite(w).run() # long time

Run the test suite for a valuation with a residual extension::

    sage: R.<x> = QQ[]
    sage: v = GaussValuation(R, pAdicValuation(QQ, 2))
    sage: w = v.augmentation(x^2 + x + 1, 1)
    sage: TestSuite(w).run() # long time

Run the test suite for an iterated residual extension starting from a
non-prime residue field::

    sage: R.<u> = Qq(4, 40)
    sage: S.<x> = R[]
    sage: v = GaussValuation(S)
    sage: w = v.augmentation(x^2 + x + u, 1/2)
    sage: TestSuite(w).run() # long time

    sage: ww = w.augmentation(x^8 + 4*x^7 + 2*x^6 + 2*x^5 + x^4 + 2*x^3 + 4*(u + 1)*x^2 + 6*(u + 1)*x + 4 + 3*u, 10)
    sage: TestSuite(ww).run() # long time

Run the test suite for an augmentation of a ramified augmentation::

    sage: R.<u> = Qq(4, 5)
    sage: S.<x> = R[]
    sage: v = GaussValuation(S)
    sage: w = v.augmentation(x, 3/4)
    sage: TestSuite(w).run() # long time

    sage: ww = w.augmentation(x^4 + 8, 5)
    sage: TestSuite(ww).run() # long time


Run the test suite for a ramified augmentation of an unramified augmentation::

    sage: R.<x> = QQ[]
    sage: v = GaussValuation(R, pAdicValuation(QQ, 2))
    sage: w = v.augmentation(x^2 + x + 1, 1)
    sage: TestSuite(w).run() # long time

    sage: ww = w.augmentation(x^4 + 2*x^3 + 5*x^2 + 8*x + 3, 16/3)
    sage: TestSuite(ww).run() # long time

Run the test suite for a ramified augmentation of a ramified augmentation::

    sage: R.<u> = Qq(4, 20)
    sage: S.<x> = R[]
    sage: v = GaussValuation(S)
    sage: w = v.augmentation(x^2 + x + u, 1/2)
    sage: TestSuite(w).run() # long time

    sage: ww = w.augmentation((x^2 + x + u)^2 + 2, 5/3)
    sage: TestSuite(ww).run() # long time

Run the test suite for another augmentation with iterated residue field extensions::

    sage: R.<u> = Qq(4, 10)
    sage: S.<x> = R[]
    sage: v = GaussValuation(S)
    sage: w = v.augmentation(x^2 + x + u, 1)
    sage: TestSuite(w).run() # long time

    sage: ww = w.augmentation((x^2 + x + u)^2 + 2*x*(x^2 + x + u) + 4*x, 3)
    sage: TestSuite(ww).run() # long time

Run the test suite for a rather trivial pseudo-valuation::

    sage: R.<u> = Qq(4, 5)
    sage: S.<x> = R[]
    sage: v = GaussValuation(S)
    sage: w = v.augmentation(x, infinity)
    sage: TestSuite(w).run() # long time

Run the test suite for an infinite valuation which extends the residue field::

    sage: R.<u> = Qq(4, 5)
    sage: S.<x> = R[]
    sage: v = GaussValuation(S)
    sage: w = v.augmentation(x^2 + x + u, infinity)
    sage: TestSuite(w).run() # long time

Run the test suite for an infinite valuation which extends a valuation which
extends the residue field::

    sage: R.<u> = Qq(4, 5)
    sage: S.<x> = R[]
    sage: v = GaussValuation(S)
    sage: w = v.augmentation(x^2 + x + u, 1/2)
    sage: TestSuite(w).run() # long time

    sage: ww = w.augmentation((x^2 + x + u)^2 + 2, infinity)
    sage: TestSuite(ww).run() # long time

REFERENCES:

.. [ML1936] Mac Lane, S. (1936). A construction for prime ideals as absolute
values of an algebraic field. Duke Mathematical Journal, 2(3), 492-510.

.. [ML1936'] MacLane, S. (1936). A construction for absolute values in
polynomial rings. Transactions of the American Mathematical Society, 40(3),
363-395.

AUTHORS:

- Julian Rüth (15-04-2013): initial version

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

from inductive_valuation import _lift_to_maximal_precision
from inductive_valuation import NonFinalInductiveValuation, FiniteInductiveValuation, InfiniteInductiveValuation, InductiveValuation
from valuation import InfiniteDiscretePseudoValuation, DiscreteValuation

from sage.misc.cachefunc import cached_method
from sage.rings.all import infinity
from sage.structure.factory import UniqueFactory

class AugmentedValuationFactory(UniqueFactory):
    r"""
    Factory for augmented valuations.

    EXAMPLES:

    This factory is not meant to be called directly. Instead,
    :meth:`augmentation` of a valuation should be called::

        sage: from mac_lane import * # optional: standalone
        sage: R.<x> = QQ[]
        sage: v = GaussValuation(R, pAdicValuation(QQ, 2))
        sage: w = v.augmentation(x, 1) # indirect doctest

    Note that trivial parts of the augmented valuation might be dropped, so you
    should not rely on ``_base_valuation`` to be the valuation you started
    with::

        sage: ww = w.augmentation(x, 2)
        sage: ww._base_valuation is v
        True

    """
    def create_key(self, base_valuation, phi, mu, check=True):
        r"""
        Create a key which uniquely identifies the valuation over
        ``base_valuation`` which sends ``phi`` to ``mu``.

        .. NOTE::

            The uniqueness that this factory provides is not why we chose to
            use a factory.  However, it makes pickling and equality checks much
            easier. At the same time, going through a factory makes it easier
            to enforce that all instances correctly inherit methods from the
            parent Hom space.

        TESTS::

            sage: from mac_lane import * # optional: standalone
            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, pAdicValuation(QQ, 2))
            sage: w = v.augmentation(x, 1) # indirect doctest
            sage: ww = v.augmentation(x, 1)
            sage: w is ww
            True

        """
        if check:
            is_key, reason = base_valuation.is_key(phi, explain=True)
            if not is_key:
                raise ValueError(reason)
            if mu <= base_valuation(phi):
                raise ValueError("the value of the key polynomial must strictly increase but `%s` does not exceed `%s`."%(mu, base_valuation(phi)))
            if not isinstance(base_valuation, InductiveValuation):
                raise TypeError("base_valuation must be inductive")

        phi = base_valuation.domain().coerce(phi)
        from sage.rings.all import QQ, infinity
        if mu is not infinity:
            mu = QQ(mu)

        if isinstance(base_valuation, AugmentedValuation_base):
            if phi.degree() == base_valuation.phi().degree():
                # drop base_valuation and extend base_valuation._base_valuation instead
                return self.create_key(base_valuation._base_valuation, phi, mu, check=check)

        return base_valuation, phi, mu

    def create_object(self, version, key):
        r"""
        Create the augmented valuation represented by ``key``.

        TESTS::

            sage: from mac_lane import * # optional: standalone
            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, pAdicValuation(QQ, 2))
            sage: w = v.augmentation(x^2 + x + 1, 1) # indirect doctest

        """
        base_valuation, phi, mu = key

        from valuation_space import DiscretePseudoValuationSpace
        parent = DiscretePseudoValuationSpace(base_valuation.domain())
        if mu < infinity:
            if base_valuation.is_trivial():
                return parent.__make_element_class__(FiniteAugmentedValuationWithTrivialResidueRing)(parent, base_valuation, phi, mu)
            else:
                return parent.__make_element_class__(FiniteAugmentedValuationWithNonTrivialResidueRing)(parent, base_valuation, phi, mu)
        else:
            return parent.__make_element_class__(InfiniteAugmentedValuation)(parent, base_valuation, phi, mu)

AugmentedValuation = AugmentedValuationFactory("AugmentedValuation")

class AugmentedValuation_base(InductiveValuation):
    """
    An augmented valuation is a discrete valuation on a polynomial ring. It
    extends another discrete valuation `v` by setting the valuation of a
    polynomial `f` to the minumum of `v(f_i)i\mu` when writing `f=\sum_i
    f_i\phi^i`.

    INPUT:

    - ``v`` -- a :class:`InductiveValuation` on a polynomial ring

    - ``phi`` -- a key polynomial over ``v`` (see :meth:`is_key`)

    - ``mu`` -- a rational number such that ``mu > v(phi)`` or ``infinity``

    EXAMPLES::

        sage: from mac_lane import * # optional: standalone
        sage: K.<u> = CyclotomicField(5)
        sage: R.<x> = K[]
        sage: v = GaussValuation(R, pAdicValuation(K, 2))
        sage: w = v.augmentation(x, 1/2); w # indirect doctest
        [ Gauss valuation induced by 2-adic valuation, v(x) = 1/2 ]
        sage: ww = w.augmentation(x^4 + 2*x^2 + 4*u, 3); ww
        [ Gauss valuation induced by 2-adic valuation, v(x) = 1/2, v(x^4 + 2*x^2 + 4*u) = 3 ]

    TESTS::

        sage: TestSuite(w).run() # long time
        sage: TestSuite(ww).run() # long time

    """
    def __init__(self, parent, v, phi, mu):
        """
        TESTS::

            sage: from mac_lane import * # optional: standalone
            sage: K.<u> = Qq(4, 5)
            sage: R.<x> = K[]
            sage: v = GaussValuation(R)
            sage: w = AugmentedValuation(v, x, 1/2)
            sage: isinstance(w, AugmentedValuation_base)
            True

            sage: TestSuite(w).run() # long time

        """
        InductiveValuation.__init__(self, parent, phi)

        self._base_valuation = v
        self._mu = mu

    def _call_(self, f):
        """
        Evaluate this valuation at ``f``.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R = Qp(2, 5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: w = v.augmentation(x^2 + x + 1, 1)

            sage: w(x^2 + 2*x + 3)
            0
            sage: w((x^2 + 2*x + 3) * w.phi()^3 )
            3
            sage: w(0)
            +Infinity

        """
        f = self.domain().coerce(f)

        from sage.rings.all import infinity
        if f.is_zero():
            return infinity

        if f.degree() < self.phi().degree():
            return self._base_valuation(f)

        # We can slightly optimize the approach of DevelopingValuation._call_
        # We know that self(f) >= self._base_valuation(f)
        # as soon as we find a coefficient of f with self._base_valuation(c) ==
        # self._base_valuation(f) we know that this is the valuation of f

        # this optimization does only pay off for polynomials of large degree:
        if f.degree() // self.phi().degree() <= 3:
            return super(AugmentedValuation_base, self)._call_(f)

        ret = infinity

        lower_bound = self._base_valuation(f)

        for v in self.valuations(f):
            ret = min(ret, v)
            if ret == lower_bound:
                break

        return ret

    def equivalence_unit(self, s):
        """
        Return an equivalence unit of minimal degree and valuation ``s``.

        INPUT:

        - ``s`` -- a rational number

        OUTPUT:

        A polynomial in the domain of this valuation which
        :meth:`is_equivalence_unit` for this valuation.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R.<u> = Qq(4, 5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: w = v.augmentation(x^2 + x + u, 1)

            sage: w.equivalence_unit(0)
            1 + O(2^5)
            sage: w.equivalence_unit(-4)
            2^-4 + O(2)

        Since an equivalence unit is of effective degree zero, `\phi` must not
        divide it. Therefore, its valuation is in the value group of the base
        valuation::

            sage: w = v.augmentation(x, 1/2)

            sage: w.equivalence_unit(3/2)
            Traceback (most recent call last):
            ...
            ValueError: 3/2 is not in the value group of 2-adic valuation
            sage: w.equivalence_unit(1)
            2 + O(2^6)

        An equivalence unit might not be integral, even if ``s >= 0``::

            sage: w = v.augmentation(x, 3/4)
            sage: ww = w.augmentation(x^4 + 8, 5)

            sage: ww.equivalence_unit(1/2)
            (2^-1 + O(2^4))*x^2

        """
        from sage.categories.fields import Fields
        if s < 0 and not self.domain().base_ring() in Fields():
            raise NotImplementedError("only implemented for polynomial rings over fields")

        ret = self._base_valuation.element_with_valuation(s)

        assert self.is_equivalence_unit(ret)
        assert self(ret) == s
        return ret

    def element_with_valuation(self, s):
        """
        Create an element of minimal degree and of valuation ``s``.

        INPUT:

        - ``s`` -- a rational number in the value group of this valuation

        OUTPUT:

        An element in the domain of this valuation

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R.<u> = Qq(4, 5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: w = v.augmentation(x^2 + x + u, 1/2)
            sage: w.element_with_valuation(0)
            1 + O(2^5)
            sage: w.element_with_valuation(1/2)
            (1 + O(2^5))*x^2 + (1 + O(2^5))*x + u + O(2^5)
            sage: w.element_with_valuation(1)
            2 + O(2^6)
            sage: c = w.element_with_valuation(-1/2); c
            (2^-1 + O(2^4))*x^2 + (2^-1 + O(2^4))*x + u*2^-1 + O(2^4)
            sage: w(c)
            -1/2
            sage: w.element_with_valuation(1/3)
            Traceback (most recent call last):
            ...
            ValueError: s must be in the value group of the valuation

        """
        if s not in self.value_group():
            raise ValueError("s must be in the value group of the valuation")
        from sage.categories.fields import Fields
        if s < 0 and not self.domain().base_ring() in Fields():
            raise NotImplementedError("only implemented for polynomial rings over fields")

        ret = self.domain().one()
        while s not in self._base_valuation.value_group():
            ret *= self._phi
            s -= self._mu
        return ret * self._base_valuation.element_with_valuation(s)

    def shift(self, x, s):
        r"""
        Return a modified version of ``x`` whose valuation is increased by ``s``.

        The element returned is such that repeated shifts which go back to
        the original valuation produce the same element in reduction.

        EXAMPLES:

            sage: from mac_lane import * # optional: standalone
            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, pAdicValuation(QQ, 2))
            sage: w = v.augmentation(x, 1)
            sage: w.shift(1, 1)
            2

        Whenever there is ramification, a shift with such consistency is not
        possible::

            sage: w = v.augmentation(x, 1/2)

            sage: w.shift(1, -1/2)
            Traceback (most recent call last):
            ...
            NotImplementedError: Shifts with consistent reduction not implemented for this augmented valuation

        Multiplication by an :meth:`element_with_valuation` might sometimes
        produce useful results in such cases::

            sage: 1 * w.element_with_valuation(-1/2)
            1/2*x

        However, this does not preserve the element in reduction::

            sage: 1 * w.element_with_valuation(-1/2) * w.element_with_valuation(1/2)
            1/2*x^2

        In general this is only possible by using an
        :meth:`equivalence_unit` and its :meth:`equialence_reciprocal`.
        These do, however, not exist for all values of ``s``.

        """
        if s not in self.value_group():
            raise ValueError("s must be in the value group of the valuation")

        if self.value_group() == self._base_valuation.value_group():
            return self._base_valuation.shift(x, s)

        if self._base_valuation.value_group().is_trivial():
            # We could implement a consistent shift in this case by multplying
            # and dividing by powers of the key polynomial. Since an element of
            # positive valuation has to be a power of the key polynomial, there
            # can be no ambiguity here
            raise NotImplementedError("Shifts with consistent reduction not implemented for augmented valuations over trivial valuations")

        # Except for very few special cases, it is not possible to implement a
        # consistent shift for augmented valuations
        raise NotImplementedError("Shifts with consistent reduction not implemented for this augmented valuation")

    def _repr_(self):
        """
        Return a printable representation of this valuation.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R.<u> = Qq(4, 5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: w = v.augmentation(x^2 + x + u, 1/2)
            sage: w # indirect doctest
            [ Gauss valuation induced by 2-adic valuation, v((1 + O(2^5))*x^2 + (1 + O(2^5))*x + u + O(2^5)) = 1/2 ]

        """
        vals = self.augmentation_chain()
        vals.reverse()
        vals = [ "v(%s) = %s"%(v._phi, v._mu) if isinstance(v, AugmentedValuation_base) else str(v) for v in vals ]
        return "[ %s ]"%", ".join(vals)

    def augmentation_chain(self):
        r"""
        Return a list with the chain of augmentations down to the underlying
        :class:`GaussValuation`.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, pAdicValuation(QQ, 2))
            sage: w = v.augmentation(x, 1)
            sage: w.augmentation_chain()
            [[ Gauss valuation induced by 2-adic valuation, v(x) = 1 ],
                 Gauss valuation induced by 2-adic valuation]

        For performance reasons, (and to simplify the underlying
        implementation,) trivial augmentations might get dropped. You should
        not rely on :meth:`augmentation_chain` to contain all the steps that
        you specified to create the current valuation::

            sage: ww = w.augmentation(x, 2)
            sage: ww.augmentation_chain()
            [[ Gauss valuation induced by 2-adic valuation, v(x) = 2 ],
                 Gauss valuation induced by 2-adic valuation]

        """
        return [self] + self._base_valuation.augmentation_chain()

    @cached_method
    def psi(self):
        """
        Return the minimal polynomial of the residue field extension of this valuation.

        OUTPUT:

        A polynomial in the residue ring of the base valuation

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R.<u> = Qq(4, 5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)

            sage: w = v.augmentation(x^2 + x + u, 1/2)
            sage: w.psi()
            x^2 + x + u0

            sage: ww = w.augmentation((x^2 + x + u)^2 + 2, 5/3)
            sage: ww.psi()
            x + 1

        """
        R = self._base_valuation.equivalence_unit(-self._base_valuation(self._phi))
        F = self._base_valuation.reduce(self._phi*R)
        assert(F.is_irreducible())
        return F

    def E(self):
        """
        Return the ramification index of this valuation over its underlying
        Gauss valuation.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R.<u> = Qq(4, 5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)

            sage: w = v.augmentation(x^2 + x + u, 1)
            sage: w.E()
            1

            sage: w = v.augmentation(x, 1/2)
            sage: w.E()
            2

        """
        if self.augmentation_chain()[-1].is_trivial():
            raise NotImplementedError("ramification index is not defined over a trivial Gauss valuation")
        return self.value_group().index(self._base_valuation.value_group()) * self._base_valuation.E()

    def F(self):
        """
        Return the degree of the residue field extension of this valuation
        over the underlying Gauss valuation.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R.<u> = Qq(4, 5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)

            sage: w = v.augmentation(x^2 + x + u, 1)
            sage: w.F()
            2

            sage: w = v.augmentation(x, 1/2)
            sage: w.F()
            1

        """
        return self.psi().degree() * self._base_valuation.F()

    def extensions(self, ring):
        r"""
        Return the extensions of this valuation to ``ring``.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, pAdicValuation(QQ, 2))
            sage: w = v.augmentation(x^2 + x + 1, 1)

            sage: w.extensions(GaussianIntegers().fraction_field()['x'])
            [[ Gauss valuation induced by 2-adic valuation, v(x^2 + x + 1) = 1 ]]
            
        """
        if ring is self.domain():
            return [self]

        from sage.rings.polynomial.polynomial_ring import is_PolynomialRing
        if is_PolynomialRing(ring) and ring.ngens() == 1:
            base_valuations = self._base_valuation.extensions(ring)
            phi = self.phi().change_ring(ring.base_ring())

            ret = []
            for v in base_valuations:
                if v.is_key(phi):
                    ret.append(AugmentedValuation(v, phi, self._mu))
                else:
                    F = v.equivalence_decomposition(phi)
                    mu0 = v(phi)
                    for f,e in F:
                        # We construct a valuation with [v, w(phi) = mu] which should be such that
                        # self(phi) = self._mu, i.e., w(phi) = w(unit) + sum e_i * w(f_i) where
                        # the sum runs over all the factors in the equivalence decomposition of phi
                        # Solving for mu gives
                        mu = (self._mu - v(F.unit()) - sum([ee*v(ff) for ff,ee in F if ff != f])) / e
                        ret.append(AugmentedValuation(v, f, mu))
            return ret

        return super(AugmentedValuation_base, self).extensions(ring)

    def restriction(self, ring):
        r"""
        Return the restriction of this valuation to ``ring``.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: K = GaussianIntegers().fraction_field()
            sage: R.<x> = K[]
            sage: v = GaussValuation(R, pAdicValuation(K, 2))
            sage: w = v.augmentation(x^2 + x + 1, 1)

            sage: w.restriction(QQ['x'])
            [ Gauss valuation induced by 2-adic valuation, v(x^2 + x + 1) = 1 ]

        """
        if ring.is_subring(self.domain()):
            base = self._base_valuation.restriction(ring)
            if ring.is_subring(self.domain().base_ring()):
                return base
            from sage.rings.polynomial.polynomial_ring import is_PolynomialRing
            if is_PolynomialRing(ring) and ring.ngens() == 1:
                return base.augmentation(self.phi().change_ring(ring.base()), self._mu)
        return super(AugmentedValuation_base, self).restriction(ring)

    def uniformizer(self):
        r"""
        Return a uniformizing element for this valuation.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, pAdicValuation(QQ, 2))
            sage: w = v.augmentation(x^2 + x + 1, 1)

            sage: w.uniformizer()
            2

        """
        return self.element_with_valuation(self.value_group()._generator)

    def is_gauss_valuation(self):
        r"""
        Return whether this valuation is a Gauss valuation.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, pAdicValuation(QQ, 2))
            sage: w = v.augmentation(x^2 + x + 1, 1)

            sage: w.is_gauss_valuation()
            False

        """
        assert(self._mu > 0)
        return False

    def monic_integral_model(self, G):
        r"""
        Return a monic integral irreducible polynomial which defines the same
        extension of the base ring of the domain as the irreducible polynomial
        ``G``.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, pAdicValuation(QQ, 2))
            sage: w = v.augmentation(x^2 + x + 1, 1)

            sage: w.monic_integral_model(5*x^2 + 1/2*x + 1/4)
            x^2 + 1/5*x + 1/5

        """
        return self._base_valuation.monic_integral_model(G)
            
    def _ge_(self, other):
        r"""
        Return whether this valuation is greater or equal than ``other``
        everywhere.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, pAdicValuation(QQ, 2))
            sage: w = v.augmentation(x^2 + x + 1, 1)
            sage: w >= v
            True
            sage: ww = v.augmentation(x^2 + x + 1, 2)
            sage: ww >= w
            True
            sage: www = w.augmentation(x^4 + 2*x^3 + 5*x^2 + 8*x + 3, 16/3)
            sage: www >= w
            True
            sage: www >= ww
            False

        """
        from gauss_valuation import GaussValuation_generic
        if other.is_trivial():
            return other.is_discrete_valuation()
        if isinstance(other, GaussValuation_generic):
            return self._base_valuation >= other
        if isinstance(other, AugmentedValuation_base):
            if self(other._phi) >= other._mu:
                return self >= other._base_valuation
            else:
                return False

        return super(AugmentedValuation_base, self)._ge_(other)


class AugmentedValuationWithTrivialResidueRing(AugmentedValuation_base, InductiveValuation):
    @cached_method
    def residue_ring(self):
        r"""
        Return the residue ring of this valuation, i.e., the elements of
        non-negative valuation modulo the elements of positive valuation.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, TrivialValuation(QQ))

            sage: w = v.augmentation(x, 1)
            sage: w.residue_ring()
            Rational Field

            sage: w = v.augmentation(x^2 + x + 1, infinity)
            sage: w.residue_ring()
            Number Field in u1 with defining polynomial x^2 + x + 1

        An example with a non-trivial base valuation::

            sage: v = GaussValuation(R, pAdicValuation(QQ, 2))
            sage: w = v.augmentation(x^2 + x + 1, infinity)
            sage: w.residue_ring()
            Finite Field in u1 of size 2^2

        Since trivial extensions of finite fields are not implemented, the
        resulting ring might be identical to the residue ring of the underlying
        valuation::

            sage: w = v.augmentation(x, infinity)
            sage: w.residue_ring()
            Finite Field of size 2

        """
        generator = 'u' + str(len(self.augmentation_chain()) - 1)

        base = self._base_valuation.residue_ring().base()
        if self.psi().degree() > 1:
            # Do not call extension() if self.psi().degree() == 1:
            # In that case the resulting field appears to be the same as the original field,
            # however, it is not == to the original field (for finite fields at
            # least) but a distinct copy (this is a bug in finite field's
            # extension() implementation.)
            return base.extension(self.psi(), names=generator)
        else:
            return base

    def reduce(self, f):
        r"""
        Reduce ``f`` module this valuation.

        OUTPUT:

        an element of the :meth:`residue_ring` of this valuation, the reduction
        modulo the ideal of elements of positive valuation

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, TrivialValuation(QQ))

            sage: w = v.augmentation(x, 1)
            sage: w.reduce(x^2 + x + 1)
            1

            sage: w = v.augmentation(x^2 + x + 1, infinity)
            sage: w.reduce(x)
            u1

        TESTS:

        Cases with non-trivial base valuation::

            sage: R.<u> = Qq(4, 10)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.reduce(x)
            x
            sage: v.reduce(S(u))
            u0

            sage: w = v.augmentation(x^2 + x + u, 1/2)
            sage: w.reduce(S.one())
            1
            sage: w.reduce(S(2))
            0
            sage: w.reduce(S(u))
            u0
            sage: w.reduce(x) # this gives the generator of the residue field extension of w over v
            u1
            sage: f = (x^2 + x + u)^2 / 2
            sage: w.reduce(f)
            x
            sage: w.reduce(f + x + 1)
            x + u1 + 1

            sage: ww = w.augmentation((x^2 + x + u)^2 + 2, 5/3)
            sage: g = ((x^2 + x + u)^2 + 2)^3 / 2^5
            sage: ww.reduce(g)
            x
            sage: ww.reduce(f)
            1
            sage: ww.is_equivalent(f, 1)
            True
            sage: ww.reduce(f * g)
            x
            sage: ww.reduce(f + g)
            x + 1

        """
        f = self.domain().coerce(f)

        if self(f) < 0:
            raise ValueError("f must have non-negative valuation")
        elif self(f) > 0:
            return self.residue_ring().zero()

        constant_term = self.coefficients(f).next()
        constant_term_reduced = self._base_valuation.reduce(constant_term)
        return constant_term_reduced(self._residue_field_generator())

    @cached_method
    def _residue_field_generator(self):
        r"""
        Return a root of :meth:`psi` in :meth:`residue_ring`.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, TrivialValuation(QQ))

            sage: w = v.augmentation(x, 1)
            sage: w._residue_field_generator()
            0

            sage: w = v.augmentation(x^2 + x + 1, infinity)
            sage: w._residue_field_generator()
            u1

        A case with non-trivial base valuation::
            
            sage: R.<u> = Qq(4, 10)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: w = v.augmentation(x^2 + x + u, infinity)
            sage: w._residue_field_generator()
            u1

        """
        if self.psi().degree() == 1:
            ret = self.residue_ring()(-self.psi()[0])
        else:
            ret = self.residue_ring().gen()

        assert self.psi()(ret).is_zero()
        return ret

    def lift(self, F):
        """
        Return a polynomial which :meth:`reduce`s to ``F``.

        INPUT:

        - ``F`` -- an element of the :meth:`residue_ring`

        ALGORITHM:

        We simply undo the steps performed in :meth:`reduce`.

        OUTPUT:

        A polynomial in the domain of the valuation with reduction ``F``.

        """
        F = self.residue_ring().coerce(F)

        if F.is_zero():
            return self.domain().zero()
        if F.is_one():
            return self.domain().one()

        # Write F as a polynomial in self._residue_field_generator()
        # We only have to do that if psi is non-trivial
        if self.psi().degree() > 1:
            from sage.rings.polynomial.polynomial_quotient_ring_element import PolynomialQuotientRingElement
            if isinstance(F, PolynomialQuotientRingElement):
                G = F.lift().change_variable_name(self._base_valuation.residue_ring().variable_name())
            else:
                G = F.polynomial(self._base_valuation.residue_ring().variable_name())
            assert(G(self._residue_field_generator()) == F)
            F = G

        H = self._base_valuation.lift(F)
        return self.domain()(H)


class AugmentedValuationWithNonTrivialResidueRing(AugmentedValuation_base, NonFinalInductiveValuation):
    @cached_method
    def residue_ring(self):
        r"""
        Return the residue ring of this valuation, i.e., the elements of
        non-negative valuation modulo the elements of positive valuation.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, pAdicValuation(QQ, 2))

            sage: w = v.augmentation(x^2 + x + 1, 1)
            sage: w.residue_ring()
            Univariate Polynomial Ring in x over Finite Field in u1 of size 2^2

        Since trivial valuations of finite fields are not implemented, the
        resulting ring might be identical to the residue ring of the underlying
        valuation::

            sage: w = v.augmentation(x, 1)
            sage: w.residue_ring()
            Univariate Polynomial Ring in x over Finite Field of size 2 (using NTL)

        """
        generator = 'u' + str(len(self.augmentation_chain()) - 1)

        base = self._base_valuation.residue_ring().base()
        if self.psi().degree() > 1:
            # Do not call extension() if self.psi().degree() == 1:
            # In that case the resulting field appears to be the same as the original field,
            # however, it is not == to the original field (for finite fields at
            # least) but a distinct copy (this is a bug in finite field's
            # extension() implementation.)
            base = base.extension(self.psi(), names=generator)
        return base[self.domain().variable_name()]

    def reduce(self, f):
        r"""
        Reduce ``f`` module this valuation.

        OUTPUT:

        an element of the :meth:`residue_ring` of this valuation, the reduction
        modulo the ideal of elements of positive valuation

        ALGORITHM:

        We follow the algorithm given in the proof of Theorem 12.1 of [ML1936]:
        If ``f`` has positive valuation, the reduction is simply zero.
        Otherwise, let `f=\sum f_i\phi^i` be the expansion of `f`, as computed
        by :meth:`coefficients`. Since the valuation is zero, the exponents `i`
        must all be multiples of `\tau`, the index the value group of the base
        valuation in the value group of this valuation.
        Hence, there is an :meth:`equivalence_unit` `Q` with the same valuation
        as `\phi^\tau`. Let `Q'` be its :meth:`reciprocal_inverse`.
        Now, rewrite each term `f_i\phi^{i\tau}=(f_iQ^i)(\phi^\tauQ^{-1})^i`;
        it turns out that the second factor in this expression is a lift of the
        generator of the :meth:`residue_field`. The reduction of the first
        factor can be computed recursively.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R.<u> = Qq(4, 10)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.reduce(x)
            x
            sage: v.reduce(S(u))
            u0

            sage: w = v.augmentation(x^2 + x + u, 1/2)
            sage: w.reduce(S.one())
            1
            sage: w.reduce(S(2))
            0
            sage: w.reduce(S(u))
            u0
            sage: w.reduce(x) # this gives the generator of the residue field extension of w over v
            u1
            sage: f = (x^2 + x + u)^2 / 2
            sage: w.reduce(f)
            x
            sage: w.reduce(f + x + 1)
            x + u1 + 1

            sage: ww = w.augmentation((x^2 + x + u)^2 + 2, 5/3)
            sage: g = ((x^2 + x + u)^2 + 2)^3 / 2^5
            sage: ww.reduce(g)
            x
            sage: ww.reduce(f)
            1
            sage: ww.is_equivalent(f, 1)
            True
            sage: ww.reduce(f * g)
            x
            sage: ww.reduce(f + g)
            x + 1

        """
        f = self.domain().coerce(f)

        if self(f) < 0:
            raise ValueError("f must have non-negative valuation")
        elif self(f) > 0:
            return self.residue_ring().zero()

        CV = zip(self.coefficients(f), self.valuations(f))
        # rewrite as sum of f_i phi^{i tau}, i.e., drop the coefficients that
        # can have no influence on the reduction
        tau = self.value_group().index(self._base_valuation.value_group())
        assert not any([v==0 for i,(c,v) in enumerate(CV) if i % tau != 0])
        CV = CV[::tau]

        # replace f_i by f_i Q^{i tau}
        vQ = self._mu * tau
        CV = [(c*self._Q()**i, v - vQ*i) for i,(c,v) in enumerate(CV)]
        assert all([self._base_valuation(c)>=0 for c,v in CV])

        # recursively reduce the f_i Q^{i tau}
        C = [self._base_valuation.reduce(c)(self._residue_field_generator()) for c,v in CV]

        # reduce the Q'^i phi^i
        return self.residue_ring()(C)

    @cached_method
    def _residue_field_generator(self):
        r"""
        Return a root of :meth:`psi` in :meth:`residue_ring`.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R.<u> = Qq(4, 10)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: w = v.augmentation(x^2 + x + u, 1/2)
            sage: w._residue_field_generator()
            u1
            
        """
        if self.residue_ring() == self._base_valuation.residue_ring():
            assert self.psi().degree() == 1
            ret = self.residue_ring().base()(-self.psi()[0])
        else:
            ret = self.residue_ring().base().gen()

        assert ret.parent() is self.residue_ring().base()
        assert self.psi()(ret).is_zero()
        return ret

    def lift(self, F):
        """
        Return a polynomial which :meth:`reduce`s to ``F``.

        INPUT:

        - ``F`` -- an element of the :meth:`residue_ring`

        OUTPUT:

        A polynomial in the domain of the valuation with reduction ``F``, monic
        if ``F`` is monic.

        ALGORITHM:

        Since this is the inverse of :meth:`reduce`, we only have to go backwards
        through the algorithm described there.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R.<u> = Qq(4, 10)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)

            sage: w = v.augmentation(x^2 + x + u, 1/2)
            sage: y = w.residue_ring().gen()
            sage: u1 = w.residue_ring().base().gen()

            sage: w.lift(1)
            1 + O(2^10)
            sage: w.lift(0)
            0
            sage: w.lift(u1)
            (1 + O(2^10))*x
            sage: w.reduce(w.lift(y)) == y
            True
            sage: w.reduce(w.lift(y + u1 + 1)) == y + u1 + 1
            True

            sage: ww = w.augmentation((x^2 + x + u)^2 + 2, 5/3)
            sage: y = ww.residue_ring().gen()
            sage: u2 = ww.residue_ring().base().gen()

            sage: ww.reduce(ww.lift(y)) == y
            True
            sage: ww.reduce(ww.lift(1)) == 1
            True
            sage: ww.reduce(ww.lift(y + 1)) == y +  1
            True

        A more complicated example::

            sage: v = GaussValuation(S)
            sage: w = v.augmentation(x^2 + x + u, 1)
            sage: ww = w.augmentation((x^2 + x + u)^2 + 2*x*(x^2 + x + u) + 4*x, 3)
            sage: u = ww.residue_ring().base().gen()

            sage: F = ww.residue_ring()(u); F
            u2
            sage: f = ww.lift(F); f
            (2^-1 + O(2^9))*x^2 + (2^-1 + O(2^9))*x + u*2^-1 + O(2^9)
            sage: F == ww.reduce(f)
            True

        """
        F = self.residue_ring().coerce(F)

        from sage.categories.fields import Fields
        if not self.domain().base_ring() in Fields():
            raise NotImplementedError("only implemented for polynomial rings over fields")

        if F.is_constant():
            if F.is_zero():
                return self.domain().zero()
            if F.is_one():
                return self.domain().one()

        R0 = self._base_valuation.residue_ring()

        # in the last step of reduce, the f_iQ^i are reduced, and evaluated at
        # the generator of the residue field
        # here, we undo this:
        coeffs = [ R0(c if self.psi().degree()==1 else list(c._vector_() if hasattr(c, '_vector_') else c.list())) for c in F.coefficients(sparse=False) ]
        coeffs = [ self._base_valuation.lift(c) for c in coeffs ]
        # now the coefficients correspond to the expansion with (f_iQ^i)(Q^{-1} phi)^i

        # now we undo the factors of Q^i (the if else is necessary to handle the case when mu is infinity, i.e., when _Q_reciprocal() is undefined)
        coeffs = [ (c if i == 0 else c*self._Q_reciprocal()**i).map_coefficients(lambda d:_lift_to_maximal_precision(d)) for i,c in enumerate(coeffs) ]

        RR = self.domain().change_ring(self.domain())

        tau = self.value_group().index(self._base_valuation.value_group())
        ret = RR(coeffs)(self.phi()**tau)
        ret = ret.map_coefficients(lambda c:_lift_to_maximal_precision(c))
        return ret

    def lift_to_key(self, F):
        """
        Lift the irreducible polynomial ``F`` to a key polynomial.

        INPUT:

        - ``F`` -- an irreducible non-constant polynomial in the
          :meth:`residue_ring` of this valuation

        OUTPUT:

        A polynomial `f` in the domain of this valuation which is a key
        polynomial for this valuation and which, for a suitable equivalence
        unit `R`, satifies that the reduction of `Rf` is ``F``

        ALGORITHM:

        We follow the algorithm described in Theorem 13.1 [ML1936] which, after
        a :meth:`lift` of ``F``, essentially shifts the valuations of all terms
        in the `\phi`-adic expansion up and then kills the leading coefficient.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R.<u> = Qq(4, 10)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)

            sage: w = v.augmentation(x^2 + x + u, 1/2)
            sage: y = w.residue_ring().gen()
            sage: f = w.lift_to_key(y + 1); f
            (1 + O(2^10))*x^4 + (2 + O(2^11))*x^3 + (1 + u*2 + O(2^10))*x^2 + (u*2 + O(2^11))*x + (u + 1) + u*2 + u*2^2 + u*2^3 + u*2^4 + u*2^5 + u*2^6 + u*2^7 + u*2^8 + u*2^9 + O(2^10)
            sage: w.is_key(f)
            True

        A more complicated example::

            sage: v = GaussValuation(S)
            sage: w = v.augmentation(x^2 + x + u, 1)
            sage: ww = w.augmentation((x^2 + x + u)^2 + 2*x*(x^2 + x + u) + 4*x, 3)

            sage: u = ww.residue_ring().base().gen()
            sage: y = ww.residue_ring().gen()
            sage: f = ww.lift_to_key(y^3+y+u)
            sage: f.degree()
            12
            sage: ww.is_key(f)
            True

        """
        F = self.residue_ring().coerce(F)

        from sage.categories.fields import Fields
        if not self.domain().base_ring() in Fields():
            raise NotImplementedError("only implemented for polynomial rings over fields")
        if self._base_valuation.is_gauss_valuation() and self._mu == infinity:
            raise TypeError("there are no keys over this valuation")

        if F.is_constant():
            raise ValueError("F must not be constant")
        if not F.is_monic():
            raise ValueError("F must be monic")
        if not F.is_irreducible():
            raise ValueError("F must be irreducible")
        if F == F.parent().gen():
            return self.phi()

        f = self.lift(F)
        assert self(f) == 0
        assert self.reduce(f) == F

        f *= self._Q()**F.degree()
        CV = zip(self.coefficients(f), self.valuations(f))
        vf = self(f)
        CV = [(c,v) if v==vf else (c.parent().zero(),infinity) for c,v in CV]
        while CV[-1][1] is infinity:
            CV.pop()

        CV[-1] = (CV[-1][0].parent().one(), vf)
        ret = self.domain().change_ring(self.domain())([c for c,v in CV])(self.phi())
        ret = ret.map_coefficients(lambda c:_lift_to_maximal_precision(c))
        assert (ret == self.phi()) == (F == F.parent().gen())
        assert self.is_key(ret)
        return ret

    @cached_method
    def _Q(self):
        r"""
        Return the polynomial `Q` used in the construction to :meth:`reduce` an
        element to the :meth:`residue_ring`.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, pAdicValuation(QQ, 2))
            sage: w = v.augmentation(x^2 + x + 1, 1)

            sage: w._Q()
            2

        """
        tau = self.value_group().index(self._base_valuation.value_group())
        return self.equivalence_unit(self._mu * tau)

    @cached_method
    def _Q_reciprocal(self):
        r"""
        Return the :meth:`equivalence_reciprocal` of :meth:`_Q`.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, pAdicValuation(QQ, 2))
            sage: w = v.augmentation(x^2 + x + 1, 1)

            sage: w._Q_reciprocal()
            1/2

        """
        ret = self.equivalence_reciprocal(self._Q())

        assert self.is_equivalence_unit(ret)
        # esentially this checks that the reduction of Q'*phi^tau is the
        # generator of the residue field
        assert self._base_valuation.reduce(self._Q()*ret)(self._residue_field_generator()).is_one()

        return ret


class FiniteAugmentedValuation(AugmentedValuation_base, FiniteInductiveValuation):
    @cached_method
    def value_group(self):
        """
        Return the value group of this valuation.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R.<u> = Qq(4, 5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)

            sage: w = v.augmentation(x^2 + x + u, 1/2)
            sage: w.value_group()
            Additive Abelian Group generated by 1/2

            sage: ww = w.augmentation((x^2 + x + u)^2 + 2, 5/3)
            sage: ww.value_group()
            Additive Abelian Group generated by 1/6

        """
        base = self._base_valuation.value_group()
        from value_group import DiscreteValueGroup
        return base + DiscreteValueGroup(self._mu)

    def valuations(self, f):
        """
        Return the valuations of the `f_i\phi^i` in the expansion `f=\sum_i
        f_i\phi^i`.

        INPUT:

        - ``f`` -- a polynomial in the domain of this valuation

        OUTPUT:

        An iterator over rational numbers (or infinity) `[v(f_0), v(f_1\phi), \dots]`

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R.<u> = Qq(4, 5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)

            sage: w = v.augmentation(x^2 + x + u, 1/2)
            sage: list(w.valuations( x^2 + 1 ))
            [0, 1/2]

            sage: ww = w.augmentation((x^2 + x + u)^2 + 2, 5/3)
            sage: list(ww.valuations( ((x^2 + x + u)^2 + 2)^3 ))
            [+Infinity, +Infinity, +Infinity, 5]

        """
        f = self.domain().coerce(f)

        for i,c in enumerate(self.coefficients(f)):
            yield self._base_valuation(c) + i*self._mu


class FiniteAugmentedValuationWithTrivialResidueRing(FiniteAugmentedValuation, AugmentedValuationWithTrivialResidueRing):
    pass


class FiniteAugmentedValuationWithNonTrivialResidueRing(FiniteAugmentedValuation, AugmentedValuationWithNonTrivialResidueRing):
    pass


class InfiniteAugmentedValuation(AugmentedValuationWithTrivialResidueRing, InfiniteInductiveValuation):
    @cached_method
    def value_group(self):
        """
        Return the value group of this valuation.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R.<u> = Qq(4, 5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: w = v.augmentation(x, infinity)
            sage: w.value_group()
            Additive Abelian Group generated by 1

        """
        return self._base_valuation.value_group()

    def valuations(self, f):
        """
        Return the valuations of the `f_i\phi^i` in the expansion `f=\sum_i
        f_i\phi^i`.

        INPUT:

        - ``f`` -- a polynomial in the domain of this valuation

        OUTPUT:

        An iterator over rational numbers (or infinity) `[v(f_0), v(f_1\phi), \dots]`

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R.<u> = Qq(4, 5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: w = v.augmentation(x, infinity)
            sage: list(w.valuations(x^2 + 1))
            [0, +Infinity, +Infinity]

        """
        f = self.domain().coerce(f)

        num_infty_coefficients = f.degree() // self.phi().degree()
        yield self._base_valuation(self.coefficients(f).next())
        for i in range(num_infty_coefficients):
            yield infinity
