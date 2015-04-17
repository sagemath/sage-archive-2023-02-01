"""
Valuations on polynomial rings based on `\phi`-adic expansions

This file implements a base class for discrete valuations on polynomial rings,
defined by a `\phi`-adic expansion.

AUTHORS:

- Julian Rueth (15-04-2013): initial version

REFERENCES:

.. [ML1936] Mac Lane, S. (1936). A construction for prime ideals as absolute
values of an algebraic field. Duke Mathematical Journal, 2(3), 492-510.

.. [ML1936'] MacLane, S. (1936). A construction for absolute values in
polynomial rings. Transactions of the American Mathematical Society, 40(3),
363-395.

TODO: Check that things work out when v is a pseudo-valuation!

"""
#*****************************************************************************
#       Copyright (C) 2013 Julian Rueth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from discrete_valuation import DiscreteValuation
from sage.misc.abstract_method import abstract_method

from sage.misc.cachefunc import cached_method

def _lift_to_maximal_precision(c):
    """
    Lift ``c`` to maximal precision if the parent is not exact.

    EXAMPLES::

        sage: from sage.rings.padics.developing_valuation import _lift_to_maximal_precision
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

class DevelopingValuation(DiscreteValuation):
    """
    An abstract base class for a discrete valuation of polynomials defined over
    the polynomial ring ``domain`` by the `\phi`-adic development.

    INPUT:

    - ``domain`` -- a polynomial ring

    - ``phi`` -- a monic element of ``domain``

    EXAMPLES::

        sage: from sage.rings.padics.developing_valuation import DevelopingValuation
        sage: R = Zp(2,5)
        sage: S.<x> = R[]
        sage: DevelopingValuation(S, x)
        `(1 + O(2^5))*x`-adic valuation of Univariate Polynomial Ring in x over 2-adic Ring with capped relative precision 5

    """
    def __init__(self, domain, phi):
        """
        Initialization.

        TESTS::

            sage: from sage.rings.padics.developing_valuation import DevelopingValuation
            sage: R = Zp(2,5)
            sage: S.<x> = R[]
            sage: v = DevelopingValuation(S, x)
            sage: type(v)
            <class 'sage.rings.padics.developing_valuation.DevelopingValuation'>

        """
        if phi.parent() is not domain:
            raise ValueError("phi must be in the domain of the valuation")
        if phi.is_constant():
            raise ValueError("phi must not be constant")
        if not phi.leading_coefficient().is_one():
            raise ValueError("phi must be monic")

        DiscreteValuation.__init__(self, domain)

        self._phi = phi

    def phi(self):
        """
        Return the polynomial `\phi`, the key polynomial of this valuation.

        EXAMPLES::

            sage: R = Zp(2,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.phi()
            (1 + O(2^5))*x

        """
        return self._phi

    def effective_degree(self, f):
        """
        Return the effective degree of ``f`` with respect to this valuation.

        The effective degree of `f` is the largest `i` such that the valuation
        of `f` and the valuation of `f_i\phi^i` in the development `f=\sum_j
        f_j\phi^j` coincide.

        INPUT:

        - ``f`` -- a non-zero polynomial in the domain of this valuation

        EXAMPLES::

            sage: R = Zp(2,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.effective_degree(x)
            1
            sage: v.effective_degree(2*x + 1)
            0

        """
        # defined on p.497 of [ML1936']
        if f.parent() is not self.domain():
            raise ValueError("f must be in the domain of the valuation")
        if f.is_zero():
            raise ValueError("the effective degree is only defined for non-zero polynomials")

        v = self(f)
        return [i for i,w in enumerate(self.valuations(f)) if w == v][-1]

    @cached_method
    def is_equivalence_irreducible(self, f):
        """
        Return whether the polynomial ``f`` is equivalence irreducible, i.e.,
        whether its :meth:`equivalence_decomposition` is irreducible.

        INPUT:

        - ``f`` -- a polynomial in the domain of this valuation

        EXAMPLES::

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
        if f.parent() is not self.domain():
            raise ValueError("f must be in the domain of the valuation")
        if f.is_constant():
            raise ValueError("f must not be constant")

        from sage.misc.cachefunc import _cache_key
        key = _cache_key(f)

        if self.equivalence_decomposition.is_in_cache(key):
            F = self.equivalence_decomposition(f)
            return len(F) <= 1 and (len(F) == 0 or F[0][1] == 1)

        if self.is_commensurable_inductive():
            # use the characterization of Theorem 13.1 in [ML1936]
            if not f.is_monic():
                raise NotImplementedError("is_equivalence_irreducible() only implemented for monic polynomials")

            # special case: phi is factor of f
            if self.valuations(f).next() > self(f):
                f = f-self.coefficients(f).next()
                assert self.phi().divides(f)
                f,_ = f.quo_rem(self.phi())
                return f.is_constant()

            R = self.equivalence_unit(-self(f))

            # check irreducibility in reduction
            F = self.reduce(f*R)
            F = F.factor()
            if len(F) > 1 or (len(F) and F[0][1] > 1):
                return False

            return True

        raise NotImplementedError("is_equivalence_irreducible() only implemented for inductive values")

    def is_equivalence_unit(self, f):
        """
        Return whether ``f`` is an equivalence unit, i.e., an element of
        :meth:`effective_degree` zero.

        INPUT:

        - ``f`` -- a polynomial in the domain of this valuation

        EXAMPLES::

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
        # defined on p.497 of [ML1936']
        if f.parent() is not self.domain():
            raise ValueError("f must be in the domain of the valuation")

        if f.is_zero():
            return False
        return self.effective_degree(f) == 0

    def equivalence_reciprocal(self, f):
        """
        Return an equivalence reciprocal of ``f``.

        An equivalence reciprocal of `f` is a polynomial `h` such that `f\cdot
        h` is equivalent to 1 modulo this valuation.

        INPUT:

        - ``f`` -- a polynomial in the domain of this valuation which is an
          equivalence unit

        EXAMPLES::

            sage: R = Zp(3,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: f = 3*x + 2
            sage: h = v.equivalence_reciprocal(f); h
            2 + 3 + 3^2 + 3^3 + 3^4 + O(3^5)
            sage: v.is_equivalent(f*h, 1)
            True

        In an extended valuation over an extension field::

            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v = v.extension(x^2 + x + u, 1)
            sage: f = 2*x + u
            sage: h = v.equivalence_reciprocal(f); h
            (u + 1) + (u + 1)*2 + 2^2 + u*2^3 + 2^4 + O(2^5)
            sage: v.is_equivalent(f*h, 1)
            True

        Extending the valuation once more::

            sage: v = v.extension((x^2 + x + u)^2 + 2*x*(x^2 + x + u) + 4*x, 3)
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
            sage: w = v.extension(t + 1, 5/4)
            sage: w = w.extension(t^4 + (a^8 + a^12 + a^14 + a^16 + a^17 + a^19 + a^20 + a^23)*t^3 + (a^6 + a^9 + a^13 + a^15 + a^18 + a^19 + a^21)*t^2 + a^10*t + 1 + a^4 + a^5 + a^8 + a^13 + a^14 + a^15, 17/2)
            sage: f = a^-15*t^2 + (a^-11 + a^-9 + a^-6 + a^-5 + a^-3 + a^-2)*t + a^-15
            sage: f_ = w.equivalence_reciprocal(f); f_
            (a^10 + a^13 + a^14 + a^17 + a^18 + a^19 + a^20 + a^24 + a^25 + O(a^26))*t^2 + a^10 + a^13 + a^18 + a^19 + a^22 + a^23 + a^24 + a^25 + O(a^26)
            sage: w.reduce(f*f_)
            1
            sage: f = f.parent()([f[0],f[1].add_bigoh(1),f[2]])
            sage: f_ = w.equivalence_reciprocal(f); f_
            (a^10 + a^13 + a^14 + a^17 + a^18 + a^19 + a^20 + a^24 + a^25 + O(a^26))*t^2 + a^10 + a^13 + a^18 + a^19 + a^22 + a^23 + a^24 + a^25 + O(a^26)
            sage: w.reduce(f*f_)
            1

        .. SEEALSO::

            :meth:`is_equivalence_unit`

        """
        # defined on p.497 of [ML1936']
        if f.parent() is not self.domain():
            raise ValueError("f must be in the domain of the valuation")
        if not self.is_equivalence_unit(f):
            raise ValueError("f must be an equivalence unit")

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

        assert self(f*h) == 0
        assert self(f*h - 1) > 0

        return h

    @cached_method
    def extension(self, phi, mu, check=True):
        """
        Return the inductive valuation which extends this valuation by mapping
        ``phi`` to ``mu``.

        INPUT:

        - ``phi`` -- a polynomial in the domain of this valuation; this must be
          a key polynomial, see :meth:`is_key` for properties of key
          polynomials.

        - ``mu`` -- a rational number, the valuation of ``phi`` in the extended
          valuation

        - ``check`` -- whether or not to check the correctness of the
          parameters

        EXAMPLES::

            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v = v.extension(x^2 + x + u, 1)
            sage: v = v.extension((x^2 + x + u)^2 + 2*x*(x^2 + x + u) + 4*x, 3)
            sage: v
            [ Gauss valuation induced by 2-adic valuation, v((1 + O(2^5))*x^2 + (1 + O(2^5))*x + u + O(2^5)) = 1, v((1 + O(2^5))*x^4 + (2^2 + O(2^6))*x^3 + (1 + (u + 1)*2 + O(2^5))*x^2 + ((u + 1)*2^2 + O(2^6))*x + (u + 1) + (u + 1)*2 + (u + 1)*2^2 + (u + 1)*2^3 + (u + 1)*2^4 + O(2^5)) = 3 ]
            sage: v.residue_field()
            Univariate Quotient Polynomial Ring in u2 over Univariate Quotient Polynomial Ring in u1 over Finite Field in u0 of size 2^2 with modulus u1^2 + u1 + u0 with modulus u2^2 + u1*u2 + u1

        .. SEEALSO::

            :class:`AugmentedValuation`

        """
        from augmented_valuation import AugmentedValuation
        return AugmentedValuation(self, phi, mu, check)

    def is_key(self, phi, explain=False):
        """
        Return whether ``phi`` is a key polynomial for this valuation.

        A key polynomial must satisfy the following conditions:

        - it must be monic
        - it must be equivalence-irreducible (see :meth:`is_equivalence_irreducible`)
        - it must be minimal (see :meth:`is_minimal`)

        INPUT:

        - ``phi`` -- a polynomial in the domain of this valuation

        - ``explain`` -- a boolean (default: ``False``), if ``True``, return a
          string explaining why ``phi`` is not a key polynomial

        EXAMPLES::

            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.is_key(x)
            True
            sage: v.is_key(2*x, explain = True)
            (False, 'phi must be monic')
            sage: v.is_key(x^2, explain = True)
            (False, 'phi must be equivalence irreducible')

            sage: w = v.extension(x, 1)
            sage: w.is_key(x + 1, explain = True)
            (False, 'phi must be minimal')

        """
        if phi.parent() is not self.domain():
            raise ValueError("phi must be in the domain of the valuation")

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

    @abstract_method
    def is_commensurable_inductive(self):
        """
        Return whether this valuation is a commensurable inductive valuation
        over the discrete valuation of the base ring of the polynomial ring, as
        defined in section 4 of [ML1936].

        EXAMPLES::

            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.is_commensurable_inductive()
            True
            sage: w = v.extension(x, 1)
            sage: w.is_commensurable_inductive()
            True

        REFERENCES:

        .. [ML1936] Mac Lane, S. (1936). A construction for prime ideals as absolute
        values of an algebraic field. Duke Mathematical Journal, 2(3), 492-510.

        """
        pass

    def is_minimal(self, f):
        """
        Return whether the polynomial ``f`` is minimal with respect to this
        valuation, as defined in definition 4.1 of [ML1936].

        INPUT:

        - ``f`` -- a polynomial in the domain of this valuation

        ALGORITHM:

        When ``f`` :meth:`is_equivalence_irreducible` for this valuation, then
        Theorem 9.4 of [ML1936'] describes what to do. TODO: what if not?

        EXAMPLES::

            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.is_minimal(x + 1)
            True
            sage: w = v.extension(x, 1)
            sage: w.is_minimal(x + 1)
            False

        TODO: An example that failed for Stefan:

            sage: K = Qp(2,10)
            sage: R.<x> = K[]
            sage: vp=pAdicValuation(K)
            sage: v0 = GaussValuation(R,vp)
            sage: f=x^5+x^4+2
            sage: v1 = v0.extension(x,1/4)
            sage: v2 = v1.extension(x^4+2,5/4)
            sage: v2.is_minimal(f)
            False

        """
        if f.parent() is not self.domain():
            raise ValueError("f must be in the domain of the valuation")
        if f.is_constant():
            raise ValueError("f must not be constant")

        if self.is_commensurable_inductive():
            # use the characterization of theorem 9.4 in [ML1936]
            if not f.is_monic():
                raise NotImplementedError("is_minimal() only implemented for monic polynomials")
            if not self.is_equivalence_irreducible(f):
                raise NotImplementedError("is_minimal() only implemented for equivalence-irreducible polynomials")
            from gauss_valuation import GaussValuation
            if isinstance(self,GaussValuation):
                return f.is_monic() and self.reduce(f).is_irreducible()
            return list(self.valuations(f))[-1] == self(f) and list(self.coefficients(f))[-1].is_constant() and list(self.valuations(f))[0] == self(f) and self.tau().divides(len(list(self.coefficients(f)))-1)

        raise NotImplementedError("is_minimal() only implemented for commensurable inductive values")

    @cached_method
    def equivalence_decomposition(self, f):
        """
        Return an equivalence decomposition of ``f``, i.e., a polynomial
        `g(x)=e(x)\prod_i \phi_i(x)` with `e(x)` an equivalence unit (see
        :meth:`is_equivalence_unit()`) and the `\phi_i` key polynomials (see
        :meth:`is_key`) such that ``f`` :meth:`is_equivalent` to `g`.

        INPUT:

        - ``f`` -- a polynomial in the domain of this valuation

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

            sage: w = v.extension(x, 1)
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

            sage: v = v.extension(x^2 + x + u, 1)
            sage: v = v.extension((x^2 + x + u)^2 + 2*x*(x^2 + x + u) + 4*x, 3)

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
            sage: v1=v0.mac_lane_step(G)[0]
            sage: V=v1.mac_lane_step(G)
            sage: v2=V[0]
            sage: v2.equivalence_decomposition(G)
            (x^4 + 4*x^3 + 6*x^2 + 4*x + alpha^4 + alpha^3 + 1)^3 * (x^4 + 4*x^3 + 6*x^2 + 4*x + 1/2*alpha^4 + alpha^3 - 27*alpha + 1)^3 * (x^4 + 4*x^3 + 6*x^2 + 4*x + 3/2*alpha^4 + alpha^3 - 27*alpha + 1)^3

        REFERENCES:

        .. [ML1936'] MacLane, S. (1936). A construction for absolute values in
        polynomial rings. Transactions of the American Mathematical Society, 40(3),
        363-395.

        """
        if f.parent() is not self.domain():
            raise ValueError("f must be in the domain of the valuation")
        if f.is_zero():
            raise ValueError("equivalence decomposition of zero is not defined")
        if any([self.constant_valuation()(c)<0 for c in f.list()]):
            raise ValueError("f must be integral")

        from sage.structure.factorization import Factorization
        if not self.domain().base_ring().is_field():
            v = self.change_ring(self.domain().base_ring().fraction_field())
            ret = v.equivalence_decomposition(v.domain()(f))
            return Factorization([(g.change_ring(self.domain().base_ring()),e) for g,e in ret], unit=ret.unit().change_ring(self.domain().base_ring()))

        if self.is_equivalence_unit(f):
            return Factorization([],unit=f)

        if not self.is_commensurable_inductive():
            raise NotImplementedError("only implemented for inductive valuations")

        f0 = f # used to check correctness of the output

        phi_divides = 0
        while self.valuations(f).next() > self(f):
            f = f-self.coefficients(f).next()
            assert self.phi().divides(f)
            f,_ = f.quo_rem(self.phi())
            phi_divides += 1

        R = self.equivalence_unit(-self(f))
        R_ = self.equivalence_reciprocal(R)

        F = self.reduce(f*R)
        F = F.factor()
        from sage.misc.misc import verbose
        verbose("%s factors as %s = %s in reduction"%(f0,F.prod(),F),caller_name="equivalence_decomposition")
        unit = F.unit()

        F = list(F)
        unit = self.lift( self.residue_ring()(unit) )

        from sage.misc.all import prod
        unit *= self.lift(self.residue_ring()(prod([ psi.leading_coefficient()**e for psi,e in F ])))
        F = [(self.lift_to_key(psi/psi.leading_coefficient()),e) for psi,e in F]

        unit *= R_ * prod([self.equivalence_unit(-self(g))**e for g,e in F])

        if phi_divides:
            for i,(g,e) in enumerate(F):
                if g == self.phi():
                    F[i] = (self.phi(),e+phi_divides)
                    break
            else:
                F.append((self.phi(),phi_divides))

        ret = Factorization(F, unit=unit)
        # assert self.is_equivalent(ret.prod(), f0) -- this might fail because of leading zeros
        assert self((ret.prod() - f0).map_coefficients(lambda c:_lift_to_maximal_precision(c)))
        assert self.is_equivalence_unit(ret.unit())
        return ret

    def minimal_representative(self, f):
        """
        Return a minimal representative for ``f``, i.e., a pair `e, a`
        such that ``f`` :meth:`is_equivalent`` to `e a`, `e` is
        an equivalence unit and `a` is minimal and monic.

        INPUT:

        - ``f`` -- a polynomial in the domain of this valuation

        OUTPUT:

        A factorization which has `e` as its unit and `a` as its unique factor.

        ALGORITHM:

        We use the algorithm described in the proof of Lemma 4.1 of [ML1936'].
        In the expansion `f=\sum_i f_i\phi^i` take `e=f_i` for the largest `i`
        with `f_i\phi^i` minimal (see :meth:`effective_degree`).
        Let `h` be the :meth:`equivalence_reciprocal` of `e` and take `a` given
        by the terms of minimal valuation in the expansion of `e f`.

        EXAMPLES::

            sage: R.<u> = Qq(4,10)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.minimal_representative(x + 2)
            (1 + O(2^10))*x

            sage: v = v.extension(x, 1)
            sage: v.minimal_representative(x + 2)
            (1 + O(2^10))*x + 2 + O(2^11)
            sage: f = x^3 + 6*x + 4
            sage: F = v.minimal_representative(f); F
            (2 + 2^2 + O(2^11)) * ((1 + O(2^10))*x + 2 + 2^2 + 2^4 + 2^6 + 2^8 + 2^10 + O(2^11))
            sage: v.is_minimal(F[0][0])
            True
            sage: v.is_equivalent(F[0][0], f)
            True

        REFERENCES:

        .. [ML1936'] MacLane, S. (1936). A construction for absolute values in
        polynomial rings. Transactions of the American Mathematical Society, 40(3),
        363-395.

        """
        if f.parent() is not self.domain():
            raise ValueError("f must be in the domain of the valuation")
        if f.is_zero():
            raise ValueError("the minimal representative of zero is not defined")

        if not self.is_commensurable_inductive():
            raise NotImplemented("only implemented for inductive valuations")

        f0 = f
        e = list(self.coefficients(f))[self.effective_degree(f)]
        f *= self.equivalence_reciprocal(e).map_coefficients(lambda c:_lift_to_maximal_precision(c))

        coeffs = [c if v == self(f) else c.parent().zero() for v,c in zip(self.valuations(f),self.coefficients(f))]
        coeffs[self.effective_degree(f0)] = self.domain().base_ring().one()
        ret = sum([c*self._phi**i for i,c in enumerate(coeffs)])
        assert self.effective_degree(ret) == self.effective_degree(f0)
        assert ret.is_monic(), coeffs
        assert self.is_minimal(ret)
        from sage.structure.factorization import Factorization
        ret = Factorization([(ret,1)],unit=e)
        # assert self.is_equivalent(ret.prod(), f0) -- this might fail because of leading zeros
        assert self((ret.prod() - f0).map_coefficients(lambda c:_lift_to_maximal_precision(c)))
        return ret

    def _normalize_leading_coefficients(self, f):
        """
        This method removes leading zero coefficients from ``f`` when
        appropriate.

        INPUT:

        - ``f`` -- a polynomial in the domain of this valuation

        OUTPUT:

        ``f`` with leading zero coefficients removed.

        .. NOTE::

            When ``f`` has leading zero coefficients, one could argue that we
            should never strip these but they often arise naturally, e.g., when
            when working with expressions as ``g-g`` or ``(g+c)-g``. We strip
            such coefficients if they are zero to sufficient precision. To be
            precise, if their precision exceeds the valuation of any other
            coefficient.
            It is not clear that this is the right way to do this.

        EXAMPLES::

            sage: R = Qp(2,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: g = x
            sage: list(v.coefficients(g-g)) # indirect doctest
            []
            sage: (g-g).list()
            [0, O(2^5)]
            sage: f = x*R(0,1) + R(1,2); f
            1 + O(2^2)
            sage: list(v.coefficients(f)) # indirect doctest
            [1 + O(2^2)]
            sage: f = x*R(0,1) + R(2,2); f
            2 + O(2^2)
            sage: list(v.coefficients(f)) # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: f must not have leading zero coefficients

        """
        if len(f.list()) > f.degree()+1:
            from sage.rings.all import infinity
            # f has leading zero coefficients
            m = min([self.constant_valuation()(c) for c in f.list()[f.degree()+1:]])
            if f.is_zero():
                f= f.parent().zero()
            elif m is infinity or m > max([self.constant_valuation()(c) for c in f.list()[:f.degree()+1]]):
                f= self.domain()(f.list()[:f.degree()+1])
            else:
                print "WARNING: DROPPING LEADING ZEROS!"
                #raise ValueError("f must not have leading zero coefficients")

        return f

    def coefficients(self, f):
        """
        Return the `\phi`-adic expansion of ``f``.

        INPUT:

        - ``f`` -- a monic polynomial in the domain of this valuation

        OUTPUT:

        An iterator `[f_0,f_1,\dots]` of polynomials in the domain of this
        valuation such that `f=\sum_i f_i\phi^i`

        EXAMPLES::

            sage: R = Qp(2,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: f = x^2 + 2*x + 3
            sage: list(v.coefficients(f)) # note that these constants are in the polynomial ring
            [1 + 2 + O(2^5), 2 + O(2^6), 1 + O(2^5)]
            sage: v = v.extension( x^2 + x + 1, 1)
            sage: list(v.coefficients(f))
            [(1 + O(2^5))*x + 2 + O(2^5), 1 + O(2^5)]

        """
        if f.parent() is not self.domain():
            raise ValueError("f must be in the domain of the valuation")
        f = self._normalize_leading_coefficients(f)

        if self.phi().degree() == 1:
            from itertools import imap
            return imap(f.parent(), f(self.phi().parent().gen() - self.phi()[0]).coefficients(sparse=False))
        else:
            return self.__coefficients(f)

    def __coefficients(self, f):
        """
        Helper method for :meth:`coefficients` to create an iterator if `\phi`
        is not linear.

        INPUT:

        - ``f`` -- a monic polynomial in the domain of this valuation

        OUTPUT:

        An iterator `[f_0,f_1,\dots]` of polynomials in the domain of this
        valuation such that `f=\sum_i f_i\phi^i`

        EXAMPLES::

            sage: R = Qp(2,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v = v.extension( x^2 + x + 1, 1)
            sage: f = x^2 + 2*x + 3
            sage: list(v.coefficients(f)) # indirect doctest
            [(1 + O(2^5))*x + 2 + O(2^5), 1 + O(2^5)]
        """
        while f.degree() >= 0:
            f,r = self.__quo_rem(f)
            yield r

    def __quo_rem(self, f):
        qr = [ self.__quo_rem_monomial(i) for i in range(f.degree()+1) ]
        q = [ f[i]*g for i,(g,_) in enumerate(qr) ]
        r = [ f[i]*h for i,(_,h) in enumerate(qr) ]
        return sum(q), sum(r)

    @cached_method
    def __quo_rem_monomial(self, degree):
        f = self.domain().one() << degree
        return f.quo_rem(self.phi())

    def newton_polygon(self, f):
        """
        Return the newton polygon the `\phi`-adic development of ``f``.

        INPUT::

        - ``f`` -- a polynomial in the domain of this valuation

        EXAMPLES::

            sage: R = Qp(2,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: f = x^2 + 2*x + 3
            sage: v.newton_polygon(f)
            Newton Polygon with vertices [(0, 0), (2, 0)]

            sage: v = v.extension( x^2 + x + 1, 1)
            sage: v.newton_polygon(f)
            Newton Polygon with vertices [(0, 0), (1, 1)]
            sage: v.newton_polygon( f * v.phi()^3 )
            Newton Polygon with vertices [(0, +Infinity), (3, 3), (4, 4)]

        .. SEEALSO::

            :class:`newton_polygon.NewtonPolygon`

        """
        if f.parent() is not self.domain():
            raise ValueError("f must be in the domain of the valuation")

        from newton_polygon import NewtonPolygon
        return NewtonPolygon(self.valuations(f))

    def _call_(self, f):
        """
        Evaluate this valuation at ``f``.

        INPUT::

        - ``f`` -- a polynomial in the domain of this valuation

        EXAMPLES::

            sage: R = Qp(2,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: f = x^2 + 2*x + 3
            sage: v(f)
            0

            sage: v = v.extension( x^2 + x + 1, 1)
            sage: v(f)
            0
            sage: v(f * v.phi()^3 )
            3
            sage: v(S.zero())
            +Infinity

        """
        if f.parent() is not self.domain():
            raise ValueError("f must be in the domain of the valuation %s but is in %s"%(self.domain(),f.parent()))

        if f.is_zero():
            from sage.rings.all import infinity
            return infinity

        return min(self.valuations(f))

    def _repr_(self):
        """
        Return a printable representation of this valuation.

        EXAMPLES::

            sage: R = Qp(2,5)
            sage: S.<x> = R[]
            sage: from sage.rings.padics.developing_valuation import DevelopingValuation
            sage: DevelopingValuation(S, x)
            `(1 + O(2^5))*x`-adic valuation of Univariate Polynomial Ring in x over 2-adic Field with capped relative precision 5

        """
        return "`%s`-adic valuation of %s"%(self._phi, self.domain())

    def residue_ring(self):
        """
        Return the residue ring of this valuation, i.e., a polynomial ring over
        the :meth:`residue_field`

        EXAMPLES::

            sage: R = Qp(2,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.residue_ring()
            Univariate Polynomial Ring in x over Finite Field of size 2 (using NTL)

        """
        return self.domain().change_ring(self.residue_field())

    def mac_lane_step(self, G, assume_squarefree=False):
        r"""

        TESTS::

            sage: K.<x>=FunctionField(QQ)
            sage: S.<y>=K[]
            sage: F=y^2-x^2-x^3-3
            sage: v0=GaussValuation(K._ring,pAdicValuation(QQ,3))
            sage: v1=v0.extension(K._ring.gen(),1/3)
            sage: from sage.rings.padics.function_field_valuation import RationalFunctionFieldValuation
            sage: mu0=RationalFunctionFieldValuation(K,v1)
            sage: eta0=GaussValuation(S,mu0)
            sage: eta1=eta0.mac_lane_step(F)[0]
            sage: eta2=eta1.mac_lane_step(F)[0]

        """
        from sage.misc.misc import verbose
        verbose("Expanding %s towards %s"%(self,G),caller_name="mac_lane_step")
        assert not G.is_constant()
        R = G.parent()
        if R is not self.domain():
            raise ValueError("G must be defined over the domain of this valuation")
        if not assume_squarefree and not G.is_squarefree():
            raise ValueError("G must be squarefree")

        from sage.rings.all import infinity

        if self(G) is infinity:
            raise ValueError("G must not have valuation infinity")

        if self.is_key(G):
            return [self.extension(G, infinity)]

        F = self.equivalence_decomposition(G)
        assert len(F), "%s factored as a unit %s"%(G,F)

        ret = []
        for phi,e in F:
            if G == phi:
                # something strange happened here:
                # G is not a key (we checked that before) but phi==G is; so phi must have less precision than G
                # this can happen if not all coefficients of G have the same precision
                # if we drop some precision of G then it will be a key
                assert not G.base_ring().is_exact()
                prec = min([c.precision_absolute() for c in phi.list()])
                g = G.map_coefficients(lambda c:c.add_bigoh(prec))
                assert self.is_key(g)
                return [self.extension(g, infinity)]

            if phi == self.phi():
                # self.phi() always is a key over self but it will not lead to an extension of this valuation
                from gauss_valuation import GaussValuation
                if isinstance(self,GaussValuation): # unless in the first step
                    pass
                elif len(F)==1: # unless this is the only factor, a terminating case which should give a valuation with v(phi)=infinity
                    pass
                else:
                    continue

            verbose("Determining the valuation for %s"%phi,level=2,caller_name="mac_lane_step")
            w = self.extension(phi, self(phi), check=False)
            NP = w.newton_polygon(G).principal_part()
            verbose("Newton-Polygon for v(phi)=%s : %s"%(self(phi),NP),level=2,caller_name="mac_lane_step")
            # assert len(NP)
            if not NP:
                q,r = G.quo_rem(phi)
                assert not r.is_zero()
                phi = phi.coefficients(sparse=False)
                for i,c in enumerate(r.coefficients(sparse=False)):
                    if not c.is_zero():
                        v = w(c)
                        # for a correct result we need to add O(pi^v) in degree i
                        # we try to find the coefficient of phi where such an error can be introduced without losing much absolute precision on phi
                        best = i
                        for j in range(i):
                            if w(q[j]) < w(q[best]):
                                best = j
                        # now add the right O() to phi in degree i-best
                        phi[i-best] = phi[i-best].add_bigoh(w(c)-w(q[best]))

                phi = G.parent()(phi)
                w = self._base_valuation.extension(phi, infinity)
                ret.append(w)
                continue

            for i in range(len(NP.slopes())):
                slope = NP.slopes()[i]
                verbose("Slope = %s"%slope,level=3,caller_name="mac_lane_step")
                side = NP.sides()[i]
                verbose("Left end is %s"%(list(w.coefficients(G))[side[0][0]]),level=3,caller_name="mac_lane_step")
                new_mu = self(phi) - slope
                base = self
                if phi.degree() == base.phi().degree():
                    assert new_mu > self(phi)
                    from gauss_valuation import GaussValuation
                    if not isinstance(base, GaussValuation):
                        base = base._base_valuation

                new_leaf = base.extension(phi, new_mu)
                assert slope is -infinity or 0 in new_leaf.newton_polygon(G).slopes()
                ret.append(new_leaf)

        assert ret
        return ret
