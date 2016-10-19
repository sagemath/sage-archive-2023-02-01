"""
Gauss valuations on polynomial rings

This file implements Gauss valuations for polynomial rings, i.e. discrete
valuations which assign to a polynomial the minimal valuation of its
coefficients.

AUTHORS:

- Julian Rueth (15-04-2013): initial version

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

    sage: v = v.extension(x, 1/4); v
    [ Gauss valuation induced by 2-adic valuation, v(x) = 1/4 ]
    sage: v = v.extension(x^4+2*x^3+2*x^2+2*x+2, 4/3); v
    [ Gauss valuation induced by 2-adic valuation, v(x) = 1/4, v(x^4 + 2*x^3 + 2*x^2 + 2*x + 2) = 4/3 ]
    sage: S.<T> = R[]
    sage: w = GaussValuation(S, v); w
    Gauss valuation induced by [ Gauss valuation induced by 2-adic valuation, v(x) = 1/4, v(x^4 + 2*x^3 + 2*x^2 + 2*x + 2) = 4/3 ]
    sage: w(2*T + 1)
    0

"""
#*****************************************************************************
#       Copyright (C) 2013 Julian Rueth <julian.rueth@fsfe.org>
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

from developing_valuation import DevelopingValuation

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
            ValueError: the domain of v must be the base ring of domain

        """
        if not domain.ngens() == 1:
            raise NotImplementedError("only implemented for univariate rings")

        if v is None:
            v = domain.base_ring().valuation()

        if not v.domain() is domain.base_ring():
            raise ValueError("the domain of v must be the base ring of domain")

        return (domain, v)

    def create_object(self, version, key, **extra_args):
        r"""
        Create a Gauss valuation from the normalized parameters.

        TESTS::

            sage: from mac_lane import * # optional: standalone
            sage: v = pAdicValuation(QQ, 2)
            sage: R.<x> = QQ[]
            sage: GaussValuation.create_object(0, (R, v))
            Gauss valuation induced by 2-adic valuation

        """
        domain, v = key
        from sage.rings.valuation.valuation_space import DiscreteValuationSpace
        parent = DiscreteValuationSpace(domain)
        return GaussValuation_generic(parent, v)

GaussValuation = GaussValuationFactory("GaussValuation")

class GaussValuation_generic(DevelopingValuation):
    """
    A Gauss valuation on a polynomial ring ``domain``.

    INPUT:

    - ``domain`` -- a polynomial ring over a valued ring `R`

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

        sage: TestSuite(v).run(skip="_test_category")

    """
    def __init__(self, parent, v):
        """
        Initialization.

        EXAMPLES::

            sage: from sage.rings.padics import GaussValuation # optional: integrated
            sage: from sage.rings.padics.gauss_valuation import GaussValuation_generic # optional: integrated
            sage: from mac_lane import * # optional: standalone
            sage: from mac_lane.gauss_valuation import GaussValuation_generic # optional: standalone
            sage: S.<x> = QQ[]
            sage: v = GaussValuation(S, pAdicValuation(QQ, 5))
            sage: isinstance(v, GaussValuation_generic)
            True

        """
        domain = parent.domain()
        DevelopingValuation.__init__(self, parent, domain.gen())

        self._base_valuation = v

    def value_group(self):
        """
        Return the value group of this valuation.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: S.<x> = QQ[]
            sage: v = GaussValuation(S, pAdicValuation(QQ, 5))
            sage: v.value_group()
            DiscreteValueGroup(1)

        """
        return self._base_valuation.value_group()

    def _repr_(self):
        """
        Return a printable representation of this valuation.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: S.<x> = QQ[]
            sage: v = GaussValuation(S, pAdicValuation(QQ, 5))
            sage: v
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
            ValueError: -s must not exceed the valuation of f

        Of course, the above example works over a field::

            sage: S.<x> = Qp(2,5)[]
            sage: v = GaussValuation(S)
            sage: f = v.shift(x, -2); f
            (2^-2 + O(2^3))*x

        """
        if f.parent() is not self.domain():
            raise ValueError("f must be in the domain of this valuation")
        if -s > self(f) and self.domain().base_ring() is not self.domain().base_ring().fraction_field():
            raise ValueError("-s must not exceed the valuation of f")
        if s not in self.value_group():
            raise ValueError("s must be in the value group of this valuation")

        return f.map_coefficients(lambda c:self._base_valuation.shift(c, s))

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

        TESTS:

        The treatment of (inexact) zero values is slightly complicated, see
        :meth:`DevelopingValuation._normalize_leading_coefficients`::

            sage: list(v.valuations(S.zero()))
            []
            sage: list(v.valuations(S([R(0,1),R(0,2)])))
            []
            sage: list(v.valuations(S([R(0,2),R(0,1)])))
            []
            sage: list(v.valuations(S([R(1,1),R(0,1)])))
            [0]
            sage: list(v.valuations(S([R(4,3),R(0,1)])))
            Traceback (most recent call last):
            ...
            ValueError: f must not have leading zero coefficients

        """
        if f.parent() is not self.domain():
            raise ValueError("f must be in the domain of this valuation")

        for c in self.coefficients(f):
            yield self._base_valuation(self.domain().base_ring()(c))

    @cached_method
    def residue_ring(self):
        """
        Return the residue ring of this valuation, i.e., the elements of
        valuation zero module the elements of positive valuation.

        EXAMPLES::

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
            ValueError: reduction not defined for non-integral elements

        .. SEEALSO::

            :meth: `lift`

        """
        if f.parent() is not self.domain():
            raise ValueError("f must be in the domain of this valuation")
        if not all([v>=0 for v in self.valuations(f)]):
            raise ValueError("reduction not defined for non-integral elements")

        return f.map_coefficients(lambda c:self._base_valuation.reduce(c), self.constant_valuation().residue_field())

    def lift(self, reduction):
        """
        Return a lift of ``reduction``.

        INPUT::

        - ``reduction`` -- a polynomial over the :meth:`residue_ring` of this valuation

        OUTPUT:

        a (possibly non-monic) polynomial in the domain of this valuation which
        reduces to ``reduction``

        EXAMPLES::

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
        if reduction.parent() is not self.residue_ring():
            raise ValueError("f must be in the residue ring of this valuation")

        return reduction.map_coefficients(lambda c:self._base_valuation.lift(c), self._base_valuation.domain())

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

            sage: R.<u> = Qq(4,10)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: y = v.residue_ring().gen()
            sage: u0 = v.residue_field().gen()
            sage: f = v.lift_to_key(y^2 + y + u0); f
            (1 + O(2^10))*x^2 + (1 + O(2^10))*x + u + O(2^10)

        """
        if F.parent() is not self.residue_ring():
            raise ValueError("F must be an element of the residue ring of the valuation")
        if F.is_constant():
            raise ValueError("F must not be constant")
        if not F.is_monic():
            raise ValueError("F must be monic")
        if not F.is_irreducible():
            raise ValueError("F must be irreducible")

        return self.lift(F)

    def constant_valuation(self):
        """
        Return the restriction of this valuations to the constants of its
        domain.

        EXAMPLES::

            sage: S.<x> = Qp(3,5)[]
            sage: v = GaussValuation(S)
            sage: v.constant_valuation()
            3-adic valuation

        """
        return self._base_valuation

    def equivalence_unit(self, s):
        """
        Return an equivalence unit of valuation ``s``.

        INPUT:

        - ``s`` -- an element of the :meth:`value_group`

        EXAMPLES::

            sage: S.<x> = Qp(3,5)[]
            sage: v = GaussValuation(S)
            sage: v.equivalence_unit(2)
            (3^2 + O(3^7))
            sage: v.equivalence_unit(-2)
            (3^-2 + O(3^3))

        """
        ret = self._base_valuation.domain().one()
        ret = self._base_valuation.shift(ret, s)
        return self.domain()(ret)

    element_with_valuation = equivalence_unit

    def is_commensurable_inductive(self):
        """
        Return whether this valuation is a commensurable inductive valuation
        over the discrete valuation of the base ring of the polynomial ring, as
        defined in section 4 of [ML1936].

        OUTPUT:

        ``True`` since a Gauss valuation always is commensurable inductive.

        EXAMPLES::

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

    def E(self):
        """
        Return the ramification index of this valuation over its underlying
        Gauss valuation, i.e., 1.

        EXAMPLES::

            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.E()
            1

        """
        from sage.rings.all import ZZ
        return ZZ.one()

    def F(self):
        """
        Return the degree of the residue field extension of this valuation
        over the Gauss valuation, i.e., 1.

        EXAMPLES::

            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.F()
            1

        """
        from sage.rings.all import ZZ
        return ZZ.one()

    def change_ring(self, base_ring):
        r"""
        Change the base ring of this valuation to ``base_ring``.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: v = pAdicValuation(ZZ, 2)
            sage: R.<x> = ZZ[]
            sage: w = GaussValuation(R, v)
            sage: w.change_ring(QQ)
            Gauss valuation induced by 2-adic valuation

        """
        base_valuation = self._base_valuation.change_ring(base_ring)
        return GaussValuation(self.domain().change_ring(base_ring), base_valuation)

    def extension(self, ring):
        from sage.rings.polynomial.polynomial_ring import is_PolynomialRing
        if is_PolynomialRing(ring) and len(ring.gens()) == 1:
            return GaussValuation(ring, self._base_valuation.extension(ring.base()))

    def is_gauss_valuation(self):
        r"""
        Return whether this valuation is a Gauss valuation.

        EXAMPLES::

            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.is_gauss_valuation()
            True

        """
        return True

    def _augmentations(self):
        r"""
        EXAMPLES::

            sage: TODO
        """
        return [self]

    def is_trivial(self):
        return self._base_valuation.is_trivial()
