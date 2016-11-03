# -*- coding: utf-8 -*-
r"""
Inductive valuations on polynomial rings

This module provides functionality for inductive valuations, i.e., finite
chains of :class:`AugmentedValuation`s on top of a :class:`GaussValuation`.

AUTHORS:

- Julian Rüth (01-11-2016): initial version

REFERENCES:

.. [ML1936] Mac Lane, S. (1936). A construction for prime ideals as absolute
values of an algebraic field. Duke Mathematical Journal, 2(3), 492-510.

.. [ML1936'] MacLane, S. (1936). A construction for absolute values in
polynomial rings. Transactions of the American Mathematical Society, 40(3),
363-395.

"""
#*****************************************************************************
#       Copyright (C) 2016 Julian Rüth <julian.rueth@fsfe.org>
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

from valuation import DiscreteValuation, InfiniteDiscretePseudoValuation
from developing_valuation import DevelopingValuation

from sage.misc.cachefunc import cached_method
from sage.misc.abstract_method import abstract_method

class InductiveValuation(DevelopingValuation):
    r"""
    Abstract base class for iterated :class:`AugmentedValuation` on top of a
    :class:`GaussValuation`.

    EXAMPLES::

        sage: from mac_lane import * # optional: standalone
        sage: R.<x> = QQ[]
        sage: v = GaussValuation(R, pAdicValuation(QQ, 5))

    TESTS::

        sage: TestSuite(v).run() # long time

    """
    def is_equivalence_unit(self, f):
        r"""
        Return whether ``f`` is an equivalence unit, i.e., an element of
        :meth:`effective_degree` zero (see [ML1936'] p.497.)

        INPUT:

        - ``f`` -- a polynomial in the domain of this valuation

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R = Zp(2,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.is_equivalence_unit(x)
            False
            sage: v.is_equivalence_unit(S.zero())
            False
            sage: v.is_equivalence_unit(2*x + 1)
            True

        """
        f = self.domain().coerce(f)

        if f.is_zero():
            return False
        return self.effective_degree(f) == 0

    def equivalence_reciprocal(self, f):
        r"""
        Return an equivalence reciprocal of ``f``.

        An equivalence reciprocal of `f` is a polynomial `h` such that `f\cdot
        h` is equivalent to 1 modulo this valuation (see [ML1936'] p.497.)

        INPUT:

        - ``f`` -- a polynomial in the domain of this valuation which is an
          :meth:`equivalence_unit`

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R = Zp(3,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: f = 3*x + 2
            sage: h = v.equivalence_reciprocal(f); h # optional: integrated (needs xgcd for polynomials with p-adic coefficients)
            2 + 3 + 3^2 + 3^3 + 3^4 + O(3^5)
            sage: v.is_equivalent(f*h, 1) # optional: integrated
            True

        In an extended valuation over an extension field::

            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v = v.augmentation(x^2 + x + u, 1)
            sage: f = 2*x + u
            sage: h = v.equivalence_reciprocal(f); h
            (u + 1) + (u + 1)*2 + 2^2 + u*2^3 + 2^4 + O(2^5)
            sage: v.is_equivalent(f*h, 1)
            True

        Extending the valuation once more::

            sage: v = v.augmentation((x^2 + x + u)^2 + 2*x*(x^2 + x + u) + 4*x, 3)
            sage: h = v.equivalence_reciprocal(f); h
            (u + 1) + (u + 1)*2 + (u + 1)*2^2 + (u + 1)*2^3 + u*2^4 + O(2^5)
            sage: v.is_equivalent(f*h, 1)
            True

        TESTS:

        A case that caused problems at some point::

            sage: K = Qp(2, 4)
            sage: R.<x> = K[]
            sage: L.<a> = K.extension(x^4 + 4*x^3 + 6*x^2 + 4*x + 2)
            sage: R.<t> = L[]
            sage: v = GaussValuation(R)
            sage: w = v.augmentation(t + 1, 5/16)
            sage: w = w.augmentation(t^4 + (a^8 + a^12 + a^14 + a^16 + a^17 + a^19 + a^20 + a^23)*t^3 + (a^6 + a^9 + a^13 + a^15 + a^18 + a^19 + a^21)*t^2 + a^10*t + 1 + a^4 + a^5 + a^8 + a^13 + a^14 + a^15, 17/8)
            sage: f = a^-15*t^2 + (a^-11 + a^-9 + a^-6 + a^-5 + a^-3 + a^-2)*t + a^-15
            sage: f_ = w.equivalence_reciprocal(f)
            sage: w.reduce(f*f_)
            1
            sage: f = f.parent()([f[0], f[1].add_bigoh(1), f[2]])
            sage: f_ = w.equivalence_reciprocal(f)
            sage: w.reduce(f*f_)
            1

        """
        f = self.domain().coerce(f)

        if not self.is_equivalence_unit(f):
            raise ValueError("f must be an equivalence unit but %r is not"%(f,))

        e0 = self.coefficients(f).next()
        one,g,h = self.phi().xgcd(e0)
        assert one.is_one()

        # it might be the case that f*h has non-zero valuation because h has
        # insufficient precision, so we must not assert that here but only
        # until we lifted to higher precision

        # We do not actually need g*phi + h*e0 = 1, it is only important that
        # the RHS is 1 in reduction.
        # This allows us to do two things:
        # - we may lift h to arbitrary precision
        # - we can add anything which times e0 has positive valuation, e.g., we
        # may drop coefficients of positive valuation
        h = h.map_coefficients(lambda c:_lift_to_maximal_precision(c))
        h = h.parent()([ c if self(e0*c) <= 0 else c.parent().zero() for c in h.coefficients(sparse=False)])

        return h

    @abstract_method
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

        Note that this might fail for negative ``s`` if the domain is not
        defined over a field::

            sage: v = pAdicValuation(ZZ, 2)
            sage: R.<x> = ZZ[]
            sage: w = GaussValuation(R, v)
            sage: w.equivalence_unit(1)
            2
            sage: w.equivalence_unit(-1)
            Traceback (most recent call last):
            ...
            TypeError: no conversion of this rational to integer

        """
        
    @abstract_method
    def augmentation_chain(self):
        r"""
        Return a list with the chain of augmentations down to the underlying
        :class:`GaussValuation`.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.augmentation_chain()
            [Gauss valuation induced by 2-adic valuation]

        """

    @abstract_method
    def is_gauss_valuation(self):
        r"""
        Return whether this valuation is a Gauss valuation over the domain.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.is_gauss_valuation()
            True

        """

    @abstract_method
    def E(self):
        """
        Return the ramification index of this valuation over its underlying
        Gauss valuation.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.E()
            1

        """

    @abstract_method
    def F(self):
        """
        Return the residual degree of this valuation over its Gauss extension.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.F()
            1

        """

    @abstract_method
    def monic_integral_model(self, G):
        r"""
        Return a monic integral irreducible polynomial which defines the same
        extension of the base ring of the domain as the irreducible polynomial
        ``G``.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, pAdicValuation(QQ, 2))
            sage: v.monic_integral_model(5*x^2 + 1/2*x + 1/4)
            x^2 + 1/5*x + 1/5

        """

    @abstract_method
    def element_with_valuation(self, s):
        r"""
        Return a polynomial of minimal degree with valuation ``s``.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, pAdicValuation(QQ, 2))
            sage: v.element_with_valuation(-2)
            1/4
            
        """

    def _test_element_with_valuation(self, **options):
        r"""
        Test the correctness of :meth:`element_with_valuation`.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, pAdicValuation(QQ, 2))
            sage: v._test_element_with_valuation()

        """
        tester = self._tester(**options)
        chain = self.augmentation_chain()
        for s in tester.some_elements(self.value_group().some_elements()):
            R = self.element_with_valuation(s)
            tester.assertEqual(self(R), s)
            if chain != [self]:
                base = chain[1]
                if s in base.value_group():
                    S = base.element_with_valuation(s)
                    tester.assertEqual(self(S), s)
                    tester.assertGreaterEqual(S.degree(), R.degree())

    def _test_EF(self, **options):
        r"""
        Test the correctness of :meth:`E` and :meth:`F`.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v._test_EF()

        """
        tester = self._tester(**options)
        chain = self.augmentation_chain()
        for w,v in zip(chain, chain[1:]):
            from sage.rings.all import infinity, ZZ
            if w(w.phi()) == infinity:
                tester.assertEqual(w.E(), v.E())
            tester.assertIn(w.E(), ZZ)
            tester.assertIn(w.F(), ZZ)
            tester.assertGreaterEqual(w.E(), v.E())
            tester.assertGreaterEqual(w.F(), v.F())

    def _test_augmentation_chain(self, **options):
        r"""
        Test the correctness of :meth:`augmentation_chain`.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, TrivialValuation(QQ))
            sage: v._test_augmentation_chain()
            
        """
        tester = self._tester(**options)
        chain = self.augmentation_chain()
        tester.assertIs(chain[0], self)
        tester.assertTrue(chain[-1].is_gauss_valuation())
        for w,v in zip(chain, chain[1:]):
            tester.assertGreaterEqual(w, v)

    def _test_equivalence_unit(self, **options):
        r"""
        Test the correctness of :meth:`lift_to_key`.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, TrivialValuation(QQ))
            sage: v._test_equivalence_unit()

        """
        tester = self._tester(**options)
        for s in tester.some_elements(self.value_group().some_elements()):
            try:
                R = self.equivalence_unit(s)
            except ValueError:
                if s >= 0 or self.domain().base_ring() in Fields():
                    raise

            tester.assertIs(R.parent(), self.domain())
            tester.assertEqual(self(R), s)
            tester.assertTrue(self.is_equivalence_unit(R))

    def _test_is_equivalence_unit(self, **options):
        r"""
        Test the correctness of :meth:`is_equivalence_unit`.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, TrivialValuation(QQ))
            sage: v._test_is_equivalence_unit()

        """
        tester = self._tester(**options)
        tester.assertFalse(self.is_equivalence_unit(self.phi()))

    def _test_equivalence_reciprocal(self, **options):
        r"""
        Test the correctness of :meth:`equivalence_reciprocal`.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, TrivialValuation(QQ))
            sage: v._test_equivalence_reciprocal()

        """
        tester = self._tester(**options)
        S = tester.some_elements(self.domain().some_elements())
        for f in S:
            if self.is_equivalence_unit(f):
                g = self.equivalence_reciprocal(f)
                tester.assertEqual(self.reduce(f*g), 1)

    def _test_inductive_valuation_inheritance(self, **options):
        r"""
        Test that every instance that is a :class:`InductiveValuation` is
        either a :class:`FiniteInductiveValuation` or a
        :class:`InfiniteInductiveValuation`.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, TrivialValuation(QQ))
            sage: v._test_inductive_valuation_inheritance()

        """
        tester = self._tester(**options)
        tester.assertTrue(isinstance(self, InfiniteInductiveValuation) != isinstance(self, FiniteInductiveValuation))

class FiniteInductiveValuation(InductiveValuation, DiscreteValuation):
    r"""
    Abstract base class for iterated :class:`AugmentedValuation` on top of a
    :class:`GaussValuation` which is a discrete valuation, i.e., the last key
    polynomial has finite valuation.

    EXAMPLES::

        sage: from mac_lane import * # optional: standalone
        sage: R.<x> = QQ[]
        sage: v = GaussValuation(R, TrivialValuation(QQ))

    TESTS::

        sage: TestSuite(v).run()

    """
    def augmentation(self, phi, mu, check=True):
        r"""
        Return the inductive valuation which extends this valuation by mapping
        ``phi`` to ``mu``.

        INPUT:

        - ``phi`` -- a polynomial in the domain of this valuation; this must be
          a key polynomial, see :meth:`is_key` for properties of key
          polynomials.

        - ``mu`` -- a rational number or infinity, the valuation of ``phi`` in
          the extended valuation

        - ``check`` -- a boolean (default: ``True``), whether or not to check
          the correctness of the parameters

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v = v.augmentation(x^2 + x + u, 1)
            sage: v = v.augmentation((x^2 + x + u)^2 + 2*x*(x^2 + x + u) + 4*x, 3)
            sage: v
            [ Gauss valuation induced by 2-adic valuation,
              v((1 + O(2^5))*x^2 + (1 + O(2^5))*x + u + O(2^5)) = 1,
              v((1 + O(2^5))*x^4 + (2^2 + O(2^6))*x^3 + (1 + (u + 1)*2 + O(2^5))*x^2 + ((u + 1)*2^2 + O(2^6))*x + (u + 1) + (u + 1)*2 + (u + 1)*2^2 + (u + 1)*2^3 + (u + 1)*2^4 + O(2^5)) = 3 ]

        TESTS:

        Make sure that we do not make the assumption that the degrees of the
        key polynomials are strictly increasing::
            
            sage: v_K = pAdicValuation(QQ,3)
            sage: A.<t> = QQ[]
            sage: v0 = GaussValuation(A,v_K)

            sage: v1 = v0.augmentation(t, 1/12)
            sage: v2 = v1.augmentation(t^12 + 3, 7/6)
            sage: v3 = v2.augmentation(t^12 + 3*t^2 + 3, 9/4)
            sage: v4 = v1.augmentation(t^12 + 3*t^2 + 3, 9/4)
            sage: v3 <= v4 and v3 >= v4 
            True

        .. SEEALSO::

            :meth:`AugmentedValuation`

        """
        from augmented_valuation import AugmentedValuation
        return AugmentedValuation(self, phi, mu, check)

    def is_key(self, phi, explain=False):
        r"""
        Return whether ``phi`` is a key polynomial for this valuation, i.e.,
        whether it is monic, whether it :meth:`is_equivalence_irreducible`, and
        whether it is :meth:`is_minimal`.

        INPUT:

        - ``phi`` -- a polynomial in the domain of this valuation

        - ``explain`` -- a boolean (default: ``False``), if ``True``, return a
          string explaining why ``phi`` is not a key polynomial

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R.<u> = Qq(4, 5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.is_key(x)
            True
            sage: v.is_key(2*x, explain = True)
            (False, 'phi must be monic')
            sage: v.is_key(x^2, explain = True)
            (False, 'phi must be equivalence irreducible')

            sage: w = v.augmentation(x, 1)
            sage: w.is_key(x + 1, explain = True)
            (False, 'phi must be minimal')

        """
        phi = self.domain().coerce(phi)

        reason = None

        if not phi.is_monic():
            reason = "phi must be monic"
        elif not self.is_equivalence_irreducible(phi):
            reason = "phi must be equivalence irreducible"
        elif not self.is_minimal(phi):
            reason = "phi must be minimal"

        if explain:
            return reason is None, reason
        else:
            return reason is None

    def is_minimal(self, f):
        r"""
        Return whether the polynomial ``f`` is minimal with respect to this
        valuation, i.e., whether ``f`` is not constant any non-constant
        polynomial `h` has at least the degree of ``f`` or ``f`` is not
        divisible by `h` with respect to this valuation, i.e., there is no `c`
        such that `c h` :meth:`is_equivalent` to `f`.

        ALGORITHM:

        Based on Theorem 9.4 of [ML1936'].

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R.<u> = Qq(4, 5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.is_minimal(x + 1)
            True
            sage: w = v.augmentation(x, 1)
            sage: w.is_minimal(x + 1)
            False

        TESTS::

            sage: K = Qp(2, 10)
            sage: R.<x> = K[]
            sage: vp = pAdicValuation(K)
            sage: v0 = GaussValuation(R, vp)
            sage: v1 = v0.augmentation(x, 1/4)
            sage: v2 = v1.augmentation(x^4 + 2, 5/4)
            sage: v2.is_minimal(x^5 + x^4 + 2) # long time
            False

        Polynomials which are equivalent to the key polynomial are minimal if
        and only if they have the same degree as the key polynomial::

            sage: v2.is_minimal(x^4 + 2)
            True
            sage: v2.is_minimal(x^4 + 4)
            False

        """
        f = self.domain().coerce(f)

        if f.is_constant():
            return False
    
        if not self.is_equivalence_irreducible(f):
            # any factor divides f with respect to this valuation
            return False

        if not f.is_monic():
            # divide out the leading factor, it does not change minimality
            v = self
            if not self.domain().base_ring().is_field():
                domain = self.domain().change_ring(self.domain().base_ring().fraction_field())
                v = self.extension(domain)
                f = domain(f)
            return v.is_minimal(f / f.leading_coefficient())
        
        if self.is_gauss_valuation():
            if self(f) == 0:
                F = self.reduce(f)
                assert not F.is_constant()
                return F.is_irreducible()
            else:
                assert(self(f) <= 0) # f is monic
                # f is not minimal:
                # Let g be f stripped of its leading term, i.e., g = f - x^n.
                # Then g and f are equivalent with respect to this valuation
                # and in particular g divides f with respect to this valuation
                return False
        
        if self.is_equivalent(self.phi(), f):
            assert f.degree() >= self.phi().degree()
            # If an h divides f with respect to this valuation, then it also divides phi:
            # v(f - c*h) > v(f) = v(c*h) => v(phi - c*h) = v((phi - f) + (f - c*h)) > v(phi) = v(c*h)
            # So if f were not minimal then phi would not be minimal but it is.
            return f.degree() == self.phi().degree()

        else:
            # see Theorem 9.4 of [ML1936']
            return list(self.valuations(f))[-1] == self(f) and \
                   list(self.coefficients(f))[-1].is_constant() and \
                   list(self.valuations(f))[0] == self(f) and \
                   self.tau().divides(len(list(self.coefficients(f))) - 1)

    def mac_lane_step(self, G, assume_squarefree=False):
        r"""
        Perform an approximation step towards the squarefree non-constant
        polynomial ``G`` with this valuation.

        This performs the individual steps that are used in
        :meth:`mac_lane_approximants`.

        TESTS::

            sage: from mac_lane import * # optional: standalone
            sage: K.<x> = FunctionField(QQ)
            sage: S.<y> = K[]
            sage: F = y^2 - x^2 - x^3 - 3
            sage: v0 = GaussValuation(K._ring,pAdicValuation(QQ, 3))
            sage: v1 = v0.augmentation(K._ring.gen(), 1/3)
            sage: mu0 = FunctionFieldValuation(K, v1)
            sage: eta0 = GaussValuation(S, mu0)
            sage: eta1 = eta0.mac_lane_step(F)[0]
            sage: eta2 = eta1.mac_lane_step(F)[0]
            sage: eta2
            [ Gauss valuation induced by Valuation on rational function field induced by [ Gauss valuation induced by 3-adic valuation, v(x) = 1/3 ], v(y + x) = 2/3 ]

        """
        G = self.domain().coerce(G)

        if G.is_constant():
            raise ValueError("G must not be constant")

        from sage.misc.misc import verbose
        verbose("Expanding %s towards %s"%(self, G), caller_name = "mac_lane_step")

        if not G.is_monic() or self(G) < 0:
            # G must be monic, there is no fundamental reason for this, but the implementation makes this assumption in some places.
            # G must be integral, otherwise, e.g., the effective degree is too low
            # We try to turn G into a monic integral polynomial that describes the same extension
            # This might fail if the constants of our polynomial ring do not form a field
            return self.mac_lane_step(self.monic_integral_model(G), assume_squarefree=assume_squarefree)

        if not assume_squarefree and not G.is_squarefree():
            raise ValueError("G must be squarefree")

        from sage.rings.all import infinity
        assert self(G) != infinity # this is a valuation and G is non-zero

        if self.is_key(G):
            return [self.augmentation(G, infinity)]

        F = self.equivalence_decomposition(G)
        assert len(F), "%s equivalence-decomposese as an equivalence-unit %s"%(G, F)

        ret = []
        for phi,e in F:
            if G == phi:
                # Something strange happened here:
                # G is not a key (we checked that before) but phi==G is; so phi must have less precision than G
                # this can happen if not all coefficients of G have the same precision
                # if we drop some precision of G then it will be a key (but is
                # that really what we should do?)
                assert not G.base_ring().is_exact()
                prec = min([c.precision_absolute() for c in phi.list()])
                g = G.map_coefficients(lambda c:c.add_bigoh(prec))
                assert self.is_key(g)
                return [self.augmentation(g, infinity)]

            if phi == self.phi():
                # a factor phi in the equivalence decomposition means that we
                # found an actual factor of G, i.e., we can set
                # v(phi)=infinity
                # However, this should already have happened in the last step
                # (when this polynomial had -infinite slope in the Newton
                # polygon.)
                if self.is_gauss_valuation(): # unless in the first step
                    pass
                else:
                    continue

            verbose("Determining the valuation for %s"%phi, level=2, caller_name="mac_lane_step")
            w = self.augmentation(phi, self(phi), check=False)
            NP = w.newton_polygon(G).principal_part()
            verbose("Newton-Polygon for v(phi)=%s : %s"%(self(phi), NP), level=2, caller_name="mac_lane_step")
            slopes = NP.slopes(repetition=False)
            if NP.vertices()[0][0] != 0:
                slopes = [-infinity] + slopes

            if not slopes:
                q,r = G.quo_rem(phi)
                assert not r.is_zero()
                phi = phi.coefficients(sparse=False)
                for i,c in enumerate(r.coefficients(sparse=False)):
                    if not c.is_zero():
                        v = w(c)
                        # for a correct result we need to add O(pi^v) in degree i
                        # we try to find the coefficient of phi where such an
                        # error can be introduced without losing much absolute
                        # precision on phi
                        best = i
                        for j in range(i):
                            if w(q[j]) < w(q[best]):
                                best = j
                        # now add the right O() to phi in degree i - best
                        phi[i-best] = phi[i-best].add_bigoh(w(c)-w(q[best]))

                phi = G.parent()(phi)
                w = self._base_valuation.augmentation(phi, infinity)
                ret.append(w)
                continue

            for i in range(len(slopes)):
                slope = slopes[i]
                verbose("Slope = %s"%slope, level=3, caller_name="mac_lane_step")
                new_mu = self(phi) - slope
                base = self
                if phi.degree() == base.phi().degree():
                    assert new_mu > self(phi)
                    if not base.is_gauss_valuation():
                        base = base._base_valuation
                new_leaf = base.augmentation(phi, new_mu)
                assert slope is -infinity or 0 in new_leaf.newton_polygon(G).slopes(repetition=False)
                ret.append(new_leaf)

        assert ret
        return ret

    @cached_method
    def _equivalence_reduction(self, f):
        r"""
        Helper method for :meth:`is_equivalence_irreducible` and
        :meth:`equivalence_decomposition` which essentially returns the
        reduction of ``f`` after multiplication with an ``R`` which
        :meth:`is_equivalence_unit`.

        This only works when ``f`` is not divisible by :meth:`phi` with respect
        to this valuation. Therefore, we also return the number of times that
        we took out :meth:`phi` of ``f`` before we computed the reduction.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, pAdicValuation(QQ, 2))
            sage: v._equivalence_reduction(2*x^6 + 4*x^5 + 2*x^4 + 8)
            (1/2, 4, x^2 + 1)

        """
        f = self.domain().coerce(f)

        # base change from R[x] to K[x], so divisions work and sufficient
        # elements of negative valuation exist
        if not self.domain().base_ring().is_field():
            domain = self.domain().change_ring(self.domain().base_ring().fraction_field())
            v = self.extension(domain)
            assert self.residue_ring() is v.residue_ring()
            return v._equivalence_reduction(f)

        phi_divides = 0
        while self.valuations(f).next() > self(f):
            # phi is an equivalence-factor of f
            f = f-self.coefficients(f).next()
            assert self.phi().divides(f)
            f,_ = f.quo_rem(self.phi())
            phi_divides += 1

        R = self.equivalence_unit(-self(f))
        return R, phi_divides, self.reduce(f*R)

    @cached_method
    def is_equivalence_irreducible(self, f):
        r"""
        Return whether the polynomial ``f`` is equivalence-irreducible, i.e.,
        whether its :meth:`equivalence_decomposition` is trivial.

        ALGORITHM:

        We use the same algorithm as in :meth:`equivalence_decomposition` we
        just do not lift the result to key polynomials.

        INPUT:

        - ``f`` -- a non-constant polynomial in the domain of this valuation

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.is_equivalence_irreducible(x)
            True
            sage: v.is_equivalence_irreducible(x^2)
            False
            sage: v.is_equivalence_irreducible(x^2 + 2)
            False

        """
        f = self.domain().coerce(f)

        if f.is_constant():
            raise ValueError("f must not be constant")

        _, phi_divides, F = self._equivalence_reduction(f)
        if phi_divides == 0:
            return F.is_constant() or F.is_irreducible()
        if phi_divides == 1:
            return F.is_constant()
        if phi_divides > 1:
            return False

    @cached_method
    def equivalence_decomposition(self, f, partial=False):
        r"""
        Return an equivalence decomposition of ``f``, i.e., a polynomial
        `g(x)=e(x)\prod_i \phi_i(x)` with `e(x)` an equivalence unit (see
        :meth:`is_equivalence_unit()`) and the `\phi_i` key polynomials (see
        :meth:`is_key`) such that ``f`` :meth:`is_equivalent` to `g`.

        INPUT:

        - ``f`` -- a non-zero polynomial in the domain of this valuation

        ALGORITHM:

        We use the algorithm described in Theorem 4.4 of [ML1936']. After
        removing all factors `\phi` from a polynomial `f`, there is an
        equivalence unit `R` such that `Rf` has valuation zero. Now `Rf` can be
        factored as `\prod_i \alpha_i` over the :meth:`residue_field`. Lifting
        all `\alpha_i` to key polynomials `\phi_i` gives `Rf=\prod_i R_i f_i`
        for suitable equivalence units `R_i` (see :meth:`lift_to_key`). Taking
        `R'` an :meth:`equivalence_reciprocal` of `R`, we have `f` equivalent
        to `(R'\prod_i R_i)\prod_i\phi_i`.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R.<u> = Qq(4,10)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.equivalence_decomposition(S.zero())
            Traceback (most recent call last):
            ...
            ValueError: equivalence decomposition of zero is not defined
            sage: v.equivalence_decomposition(S.one())
            1 + O(2^10)
            sage: v.equivalence_decomposition(x^2+2)
            ((1 + O(2^10))*x)^2
            sage: v.equivalence_decomposition(x^2+1)
            ((1 + O(2^10))*x + 1 + O(2^10))^2

        A polynomial that is an equivalence unit, is returned as the unit part
        of a :class:`sage.structure.factorization.Factorization`, leading to a unit
        non-minimal degree::

            sage: w = v.augmentation(x, 1)
            sage: F = w.equivalence_decomposition(x^2+1); F
            (1 + O(2^10))*x^2 + 1 + O(2^10)
            sage: F.unit()
            (1 + O(2^10))*x^2 + 1 + O(2^10)

        However, if the polynomial has a non-unit factor, then the unit might
        be replaced by a factor of lower degree::

            sage: f = x * (x^2 + 1)
            sage: F = w.equivalence_decomposition(f); F
            (1 + O(2^10))*x
            sage: F.unit()
            1 + O(2^10)

        Examples over an iterated unramified extension:

            sage: v = v.augmentation(x^2 + x + u, 1)
            sage: v = v.augmentation((x^2 + x + u)^2 + 2*x*(x^2 + x + u) + 4*x, 3)

            sage: v.equivalence_decomposition(x)
            (1 + O(2^10))*x
            sage: F = v.equivalence_decomposition( v.phi() )
            sage: len(F)
            1
            sage: F = v.equivalence_decomposition( v.phi() * (x^4 + 4*x^3 + (7 + 2*u)*x^2 + (8 + 4*u)*x + 1023 + 3*u) )
            sage: len(F)
            2

        TESTS::

            sage: R.<x> = QQ[]
            sage: K1.<pi>=NumberField(x^3-2)
            sage: K.<alpha>=K1.galois_closure()
            sage: R.<x>=K[]
            sage: vp=pAdicValuation(QQ,2)
            sage: vp=vp.extension(K)
            sage: v0=GaussValuation(R,vp)
            sage: G=x^36 + 36*x^35 + 630*x^34 + 7144*x^33 + 59055*x^32 + 379688*x^31 +1978792*x^30 + 8604440*x^29 + 31895428*x^28 + 102487784*x^27 + 289310720*x^26 + 725361352*x^25 + 1629938380*x^24 + 3307417800*x^23 + 6098786184*x^22+10273444280*x^21 + 15878121214*x^20 + 22596599536*x^19 + 29695703772*x^18 +36117601976*x^17 + 40722105266*x^16 + 42608585080*x^15 + 41395961848*x^14 +37344435656*x^13 + 31267160756*x^12 + 24271543640*x^11 + 17439809008*x^10 + 11571651608*x^9 + 7066815164*x^8 + 3953912472*x^7 + 2013737432*x^6 + 925014888*x^5 + 378067657*x^4 + 134716588*x^3 + 40441790*x^2 + 9532544*x + 1584151
            sage: v1=v0.mac_lane_step(G)[0] # long time
            sage: V=v1.mac_lane_step(G) # long time
            sage: v2=V[0] # long time
            sage: v2.equivalence_decomposition(G) # long time
            (x^4 + 4*x^3 + 6*x^2 + 4*x + alpha^4 + alpha^3 + 1)^3 * (x^4 + 4*x^3 + 6*x^2 + 4*x + 1/2*alpha^4 + alpha^3 - 27*alpha + 1)^3 * (x^4 + 4*x^3 + 6*x^2 + 4*x + 3/2*alpha^4 + alpha^3 - 27*alpha + 1)^3

        REFERENCES:

        .. [ML1936'] MacLane, S. (1936). A construction for absolute values in
        polynomial rings. Transactions of the American Mathematical Society, 40(3),
        363-395.

        """
        f = self.domain().coerce(f)

        if f.is_zero():
            raise ValueError("equivalence decomposition of zero is not defined")

        from sage.structure.factorization import Factorization
        if self.is_equivalence_unit(f):
            return Factorization([],unit=f)

        if not self.domain().base_ring().is_field():
            domain = self.domain().change_ring(self.domain().base_ring().fraction_field())
            v = self.extension(domain)
            ret = v.equivalence_decomposition(v.domain()(f))
            return Factorization([(g.change_ring(self.domain().base_ring()),e) for g,e in ret], unit=ret.unit().change_ring(self.domain().base_ring()))

        R, phi_divides, F = self._equivalence_reduction(f)
        F = F.factor()
        from sage.misc.misc import verbose
        verbose("%s factors as %s = %s in reduction"%(f, F.prod(), F), caller_name="equivalence_decomposition")

        unit = self.lift(self.residue_ring()(F.unit())) * self.equivalence_reciprocal(R)
        F = list(F)

        from sage.misc.all import prod
        unit *= self.lift(self.residue_ring()(prod([ psi.leading_coefficient()**e for psi,e in F ])))
        F = [(self.lift_to_key(psi/psi.leading_coefficient()),e) for psi,e in F]
        unit *= prod([self.equivalence_unit(-self(g))**e for g,e in F])

        if phi_divides:
            for i,(g,e) in enumerate(F):
                if g == self.phi():
                    F[i] = (self.phi(),e+phi_divides)
                    break
            else:
                F.append((self.phi(),phi_divides))

        ret = Factorization(F, unit=unit)

        assert self.is_equivalent(ret.prod(), f) # this might fail because of leading zeros in inexact rings
        assert self.is_equivalence_unit(ret.unit())

        return ret

    def minimal_representative(self, f):
        r"""
        Return a minimal representative for ``f``, i.e., a pair `e, a` such
        that ``f`` :meth:`is_equivalent` to `e a`, `e` is an
        :meth:`equivalence_unit`, and `a` :meth:`is_minimal` and monic.

        INPUT:

        - ``f`` -- a non-zero polynomial which is not an equivalence unit

        OUTPUT:

        A factorization which has `e` as its unit and `a` as its unique factor.

        ALGORITHM:

        We use the algorithm described in the proof of Lemma 4.1 of [ML1936'].
        In the expansion `f=\sum_i f_i\phi^i` take `e=f_i` for the largest `i`
        with `f_i\phi^i` minimal (see :meth:`effective_degree`).
        Let `h` be the :meth:`equivalence_reciprocal` of `e` and take `a` given
        by the terms of minimal valuation in the expansion of `e f`.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R.<u> = Qq(4,10)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.minimal_representative(x + 2)
            (1 + O(2^10))*x

            sage: v = v.augmentation(x, 1)
            sage: v.minimal_representative(x + 2)
            (1 + O(2^10))*x + 2 + O(2^11)
            sage: f = x^3 + 6*x + 4
            sage: F = v.minimal_representative(f); F
            (2 + 2^2 + O(2^11)) * ((1 + O(2^10))*x + 2 + 2^2 + 2^4 + 2^6 + 2^8 + 2^10 + O(2^11))
            sage: v.is_minimal(F[0][0])
            True
            sage: v.is_equivalent(F.prod(), f)
            True

        """
        f = self.domain().coerce(f)

        if f.is_zero():
            raise ValueError("zero has no minimal representative")

        degree = self.effective_degree(f)
        if degree == 0:
            raise ValueError("equivalence units can not have a minimal representative")

        e = list(self.coefficients(f))[degree]
        h = self.equivalence_reciprocal(e).map_coefficients(lambda c:_lift_to_maximal_precision(c))
        g = h*f
        vg = self(g)

        coeffs = [c if v == vg else c.parent().zero() for v,c in zip(self.valuations(g), self.coefficients(g))]
        coeffs[degree] = self.domain().base_ring().one()
        ret = sum([c*self._phi**i for i,c in enumerate(coeffs)])

        assert self.effective_degree(ret) == degree
        assert ret.is_monic()
        assert self.is_minimal(ret)

        from sage.structure.factorization import Factorization
        ret = Factorization([(ret, 1)], unit=e)

        assert self.is_equivalent(ret.prod(), f) # this might fail because of leading zeros
        return ret

    @abstract_method
    def lift_to_key(self, F):
        """
        Lift the irreducible polynomial ``F`` from the ;meth:`residue_ring` to
        a key polynomial over this valuation.

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

    def _test_lift_to_key(self, **options):
        r"""
        Test the correctness of :meth:`lift_to_key`.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, TrivialValuation(QQ))
            sage: v._test_lift_to_key()

        """
        tester = self._tester(**options)
        S = tester.some_elements(self.residue_ring().some_elements())
        for F in S:
            if F.is_monic() and not F.is_constant() and F.is_irreducible():
                f = self.lift_to_key(F)
                tester.assertIs(f.parent(), self.domain())
                tester.assertTrue(self.is_key(f))

                from sage.categories.fields import Fields
                if self.domain().base_ring() in Fields():
                    R = self.equivalence_unit(-self(f))
                    tester.assertEqual(self.reduce(f*R), F)

    def _test_is_equivalence_irreducible(self, **options):
        r"""
        Test the correctness of :meth:`is_equivalence_irreducible`.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, TrivialValuation(QQ))
            sage: v._test_is_equivalence_irreducible()

        """
        tester = self._tester(**options)
        S = tester.some_elements(self.domain().some_elements())
        for f in S:
            if f.is_constant(): continue
            is_equivalence_irreducible = self.is_equivalence_irreducible(f)
            F = self.equivalence_decomposition(f)
            tester.assertEqual(is_equivalence_irreducible, len(F)==0 or (len(F)==1 and F[0][1]==1))
            if self.is_equivalence_unit(f):
                tester.assertTrue(f.is_constant() or self.is_equivalence_irreducible(f))

        tester.assertTrue(self.is_equivalence_irreducible(self.phi()))
        tester.assertTrue(self.is_equivalence_irreducible(-self.phi()))
        tester.assertFalse(self.is_equivalence_irreducible(self.phi() ** 2))


class InfiniteInductiveValuation(InductiveValuation, InfiniteDiscretePseudoValuation):
    r"""
    Abstract base class for iterated :class:`AugmentedValuation` on top of a
    :class:`GaussValuation` which is not discrete valuation, i.e., the last key
    polynomial has infinite valuation.

    EXAMPLES::

        sage: from mac_lane import * # optional: standalone
        sage: R.<x> = QQ[]
        sage: v = GaussValuation(R, pAdicValuation(QQ, 2))
        sage: w = v.augmentation(x^2 + x + 1, infinity)

    TESTS::

        sage: TestSuite(w).run()

    """
    pass


def _lift_to_maximal_precision(c):
    r"""
    Lift ``c`` to maximal precision if the parent is not exact.

    EXAMPLES::

        sage: from mac_lane import * # optional: standalone
        sage: R = Zp(2,5)
        sage: x = R(1,2); x
        1 + O(2^2)
        sage: _lift_to_maximal_precision(x)
        1 + O(2^5)

        sage: x = 1
        sage: _lift_to_maximal_precision(x)
        1

    """
    return c if c.parent().is_exact() else c.lift_to_precision()

