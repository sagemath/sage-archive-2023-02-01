r"""
Discrete valuations

This file defines abstract base classes for discrete (pseudo-)valuations.

AUTHORS:

- Julian Rueth (2013-03-16): initial version

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

from sage.categories.morphism import Morphism

class DiscretePseudoValuation(Morphism):
    r"""
    Abstract base class for discrete pseudo-valuations, i.e., discrete
    valuations which might send more that just zero to infinity.

    INPUT:

    - ``domain`` -- an integral domain

    EXAMPLES::

        sage: from mac_lane import * # optional: standalone
        sage: v = pAdicValuation(ZZ, 2); v # indirect doctest
        2-adic valuation

    TESTS::

        sage: TestSuite(v).run()

    """
    def __init__(self, parent):
        r"""
        TESTS::

            sage: from mac_lane import * # optional: standalone
            sage: isinstance(pAdicValuation(ZZ, 2), DiscretePseudoValuation)
            True

        """
        Morphism.__init__(self, parent=parent)

    def is_equivalent(self, f, g):
        r"""
        Return whether ``f`` and ``g`` are equivalent.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: v = pAdicValuation(QQ, 2)
            sage: v.is_equivalent(2, 1)
            False
            sage: v.is_equivalent(2, -2)
            True
            sage: v.is_equivalent(2, 0)
            False
            sage: v.is_equivalent(0, 0)
            True

        """
        vf = self(f)
        vg = self(g)
        from sage.rings.all import infinity
        if self(f) == infinity and self(g) == infinity:
            return True
        if self(f) == infinity or self(g) == infinity:
            return False
        return self(f-g) > vf

    def __hash__(self):
        r"""
        The hash value of this valuation.

        We redirect to :meth:`_hash_`, so that subclasses can only override
        :meth:`_hash_` and :meth:`_eq` if they want to provide a different
        notion of equality but they can leave the partial and total operators
        untouched.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: v = pAdicValuation(QQ, 2)
            sage: hash(v) == hash(v) # indirect doctest
            True

        """
        return self._hash_()

    def _hash_(self):
        r"""
        Return a hash value for this valuation.

        We override the strange default provided by
        ``sage.categories.marphism.Morphism`` here and implement equality by
        ``id``. This works fine for objects which use unique representation
        which is the case for most valuations.
        """
        return id(self)

    def _cmp_(self, other):
        r"""
        Compare this element to ``other``.

        Since there is no reasonable total order on valuations, this method
        just throws an exception.

        EXAMPLES::

            sage: from mac_lane import *
            sage: v = pAdicValuation(QQ, 2)
            sage: v < v
            Traceback (most recent call last):
            ...
            NotImplementedError: Operator not implemented for this valuation.

        Note that this does not affect comparison of valuations which do not
        coerce into a common parent. This is by design in Sage, see
        :meth:`sage.structure.element.Element.__cmp__`. When the valuations do
        not coerce into a common parent, a rather random comparison of ``id``
        happens::

            sage: w = TrivialValuation(GF(2))
            sage: w < v # random output
            True
            sage: v < w # random output
            False

        """
        raise NotImplementedError("No total order for these valuations.")

    def _richcmp_(self, other, op):
        r"""
        Compare this element to ``other``.

        We redirect to methods :meth:`_eq_`, :meth:`_lt_`, and :meth:`_gt_` to
        make it easier for subclasses to override only parts of this
        functionality.

        Note that valuations usually implement ``x == y`` as ``x`` and ``y``
        are indistinguishable. Whereas ``x <= y`` and ``x >= y`` are
        implemented with respect to the natural partial order of valuations.
        As a result, ``x <= y and x >= y`` does not imply ``x == y``.

        EXAMPLES::

            sage: from mac_lane import *
            sage: v = pAdicValuation(QQ, 2)
            sage: v == v
            True
            sage: v != v
            False
            sage: w = pAdicValuation(QQ, 3)
            sage: v == w
            False
            sage: v != w
            True

        Note that this does not affect comparison of valuations which do not
        coerce into a common parent. This is by design in Sage, see
        :meth:`sage.structure.element.Element.__richcmp__`. When the valuations
        do not coerce into a common parent, a rather random comparison of
        ``id`` happens::

            sage: w = TrivialValuation(GF(2))
            sage: w <= v # random output
            True
            sage: v <= w # random output
            False

        """
        if op == 1: # <=
            return self._lt_(other)
        if op == 2: # ==
            return self._eq_(other)
        if op == 3: # !=
            return not self == other
        if op == 5: # >=
            return self._gt_(other)
        raise NotImplementedError("Operator not implemented for this valuation.")

    def _eq_(self, other):
        r"""
        Return whether this valuation and ``other`` are indistinguishable.

        We override the strange default provided by
        ``sage.categories.marphism.Morphism`` here and implement equality by
        ``id``. This is the right behaviour in many cases.

        When overriding this method, you can assume that ``other`` is a
        (pseudo-)valuation on the same domain.
        """
        return self is other

    def _lt_(self, other):
        r"""
        Return whether this valuation is less than or equal to ``other``
        pointwise.

        When overriding this method, you can assume that ``other`` is a
        (pseudo-)valuation on the same domain.
        """
        return other >= self

    def _gt_(self, other):
        r"""
        Return whether this valuation is greater than or equal to ``other``
        pointwise.

        When overriding this method, you can assume that ``other`` is a
        (pseudo-)valuation on the same domain.
        """
        if self == other: return True
        raise NotImplementedError("Operator not implemented for this valuation.")

    # Remove the default implementation of Map.__reduce__ that does not play
    # nice with factories (a factory, does not override Map.__reduce__ because
    # it is not the generic reduce of object) and that does not match equality
    # by id.
    __reduce__ = object.__reduce__

class InfiniteDiscretePseudoValuation(DiscretePseudoValuation):
    r"""
    sage: TODO
    """
    def is_discrete_valuation(self):
        r"""
        sage: TODO
        """
        return False

class DiscreteValuation(DiscretePseudoValuation):
    r"""
    sage: TODO
    """
    def is_discrete_valuation(self):
        r"""
        Return whether this valuation is a discrete valuation.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: v = TrivialValuation(ZZ)
            sage: v.is_discrete_valuation()
            True

        """
        return True

    def mac_lane_approximants(self, G, precision_cap=None, assume_squarefree=False):
        r"""
        sage: TODO

        TESTS:

        Some difficult cases provided by Mark van Hoeij::

            sage: from mac_lane import * # optional: standalone
            sage: k = GF(2)
            sage: K.<x> = FunctionField(k)
            sage: R.<y> = K[]
            sage: F = y^21 + x*y^20 + (x^3 + x + 1)*y^18 + (x^3 + 1)*y^17 + (x^4 + x)*y^16 + (x^7 + x^6 + x^3 + x + 1)*y^15 + x^7*y^14 + (x^8 + x^7 + x^6 + x^4 + x^3 + 1)*y^13 + (x^9 + x^8 + x^4 + 1)*y^12 + (x^11 + x^9 + x^8 + x^5 + x^4 + x^3 + x^2)*y^11 + (x^12 + x^9 + x^8 + x^7 + x^5 + x^3 + x + 1)*y^10 + (x^14 + x^13 + x^10 + x^9 + x^8 + x^7 + x^6 + x^3 + x^2 + 1)*y^9 + (x^13 + x^9 + x^8 + x^6 + x^4 + x^3 + x)*y^8 + (x^16 + x^15 + x^13 + x^12 + x^11 + x^7 + x^3 + x)*y^7 + (x^17 + x^16 + x^13 + x^9 + x^8 + x)*y^6 + (x^17 + x^16 + x^12 + x^7 + x^5 + x^2 + x + 1)*y^5 + (x^19 + x^16 + x^15 + x^12 + x^6 + x^5 + x^3 + 1)*y^4 + (x^18 + x^15 + x^12 + x^10 + x^9 + x^7 + x^4 + x)*y^3 + (x^22 + x^21 + x^20 + x^18 + x^13 + x^12 + x^9 + x^8 + x^7 + x^5 + x^4 + x^3)*y^2 + (x^23 + x^22 + x^20 + x^17 + x^15 + x^14 + x^12 + x^9)*y + x^25 + x^23 + x^19 + x^17 + x^15 + x^13 + x^11 + x^5
            sage: x = K._ring.gen()
            sage: v0 = FunctionFieldValuation(K, GaussValuation(K._ring, TrivialValuation(k)).augmentation(x,1))
            sage: v0.mac_lane_approximants(F, assume_squarefree=True) # optional: integrated; assumes squarefree for speed
            [[ Gauss valuation induced by Valuation on rational function field induced by [ Gauss valuation induced by Trivial valuation, v(x) = 1 ], v(y + x + 1) = 3/2 ],
             [ Gauss valuation induced by Valuation on rational function field induced by [ Gauss valuation induced by Trivial valuation, v(x) = 1 ], v(y) = 4/3, v(y^3 + x^4) = 13/3 ],
             [ Gauss valuation induced by Valuation on rational function field induced by [ Gauss valuation induced by Trivial valuation, v(x) = 1 ], v(y + x) = 2 ],
             [ Gauss valuation induced by Valuation on rational function field induced by [ Gauss valuation induced by Trivial valuation, v(x) = 1 ], v(y^15 + y^13 + (x + 1)*y^12 + x*y^11 + (x + 1)*y^10 + y^9 + y^8 + x*y^6 + x*y^5 + y^4 + y^3 + y^2 + (x + 1)*y + x + 1) = 2 ]]
            sage: v0 = FunctionFieldValuation(K, GaussValuation(K._ring, TrivialValuation(k)).augmentation(x+1,1))
            sage: v0.mac_lane_approximants(F, assume_squarefree=True) # optional: integrated; assumes squarefree for speed
            [[ Gauss valuation induced by Valuation on rational function field induced by [ Gauss valuation induced by Trivial valuation, v(x + 1) = 1 ], v(y) = 7/2, v(y^2 + x^7 + x^6 + x^5 + x^4 + x^3 + x^2 + x + 1) = 15/2 ],
             [ Gauss valuation induced by Valuation on rational function field induced by [ Gauss valuation induced by Trivial valuation, v(x + 1) = 1 ], v(y + x^2 + 1) = 7/2 ],
             [ Gauss valuation induced by Valuation on rational function field induced by [ Gauss valuation induced by Trivial valuation, v(x + 1) = 1 ], v(y) = 3/4, v(y^4 + x^3 + x^2 + x + 1) = 15/4 ],
             [ Gauss valuation induced by Valuation on rational function field induced by [ Gauss valuation induced by Trivial valuation, v(x + 1) = 1 ], v(y^13 + x*y^12 + y^10 + (x + 1)*y^9 + (x + 1)*y^8 + x*y^7 + x*y^6 + (x + 1)*y^4 + y^3 + (x + 1)*y^2 + 1) = 2 ]]
            sage: v0 = FunctionFieldValuation(K, GaussValuation(K._ring, TrivialValuation(k)).augmentation(x^3+x^2+1,1))
            sage: v0.mac_lane_approximants(F, assume_squarefree=True) # optional: integrated; assumes squarefree for speed
            [[ Gauss valuation induced by Valuation on rational function field induced by [ Gauss valuation induced by Trivial valuation, v(x^3 + x^2 + 1) = 1 ], v(y + x^3 + x^2 + x) = 2, v(y^2 + (x^6 + x^4 + 1)*y + x^14 + x^10 + x^9 + x^8 + x^5 + x^4 + x^3 + x^2 + x) = 5 ],
             [ Gauss valuation induced by Valuation on rational function field induced by [ Gauss valuation induced by Trivial valuation, v(x^3 + x^2 + 1) = 1 ], v(y^2 + (x^7 + x^5 + x^4 + x^3 + x^2 + x)*y + x^7 + x^5 + x + 1) = 3 ],
             [ Gauss valuation induced by Valuation on rational function field induced by [ Gauss valuation induced by Trivial valuation, v(x^3 + x^2 + 1) = 1 ], v(y^3 + (x^8 + x^5 + x^4 + x^3 + x + 1)*y^2 + (x^7 + x^6 + x^5)*y + x^8 + x^5 + x^4 + x^3 + 1) = 3 ],
             [ Gauss valuation induced by Valuation on rational function field induced by [ Gauss valuation induced by Trivial valuation, v(x^3 + x^2 + 1) = 1 ], v(y^3 + (x^8 + x^4 + x^3 + x + 1)*y^2 + (x^4 + x^3 + 1)*y + x^8 + x^7 + x^4 + x + 1) = 3 ],
             [ Gauss valuation induced by Valuation on rational function field induced by [ Gauss valuation induced by Trivial valuation, v(x^3 + x^2 + 1) = 1 ], v(y^4 + (x^8 + x^7 + x^6 + x^5 + x^4 + x^3 + x^2 + x + 1)*y^3 + (x^8 + x^5 + x^4 + x^3 + x^2 + x + 1)*y^2 + (x^8 + x^7 + x^6 + x^5 + x^3 + x^2 + 1)*y + x^8 + x^7 + x^6 + x^5 + x^3 + 1) = 3 ],
             [ Gauss valuation induced by Valuation on rational function field induced by [ Gauss valuation induced by Trivial valuation, v(x^3 + x^2 + 1) = 1 ], v(y^7 + (x^8 + x^5 + x^4 + x)*y^6 + (x^7 + 1)*y^5 + (x^4 + x^2)*y^4 + (x^8 + x^3 + x + 1)*y^3 + (x^7 + x^6 + x^4 + x^2 + x + 1)*y^2 + (x^8 + x^7 + x^5 + x^3 + 1)*y + x^7 + x^6 + x^5 + x^4 + x^3 + x^2) = 3 ]]

        Cases with trivial residue field extensions::

            sage: K.<x> = FunctionField(QQ)
            sage: S.<y> = K[]
            sage: F = y^2 - x^2 - x^3 - 3
            sage: v0 = GaussValuation(K._ring,pAdicValuation(QQ,3))
            sage: v1 = v0.augmentation(K._ring.gen(),1/3)
            sage: mu0 = FunctionFieldValuation(K, v1)
            sage: mu0.mac_lane_approximants(F)
            [[ Gauss valuation induced by Valuation on rational function field induced by [ Gauss valuation induced by 3-adic valuation, v(x) = 1/3 ], v(y + x) = 2/3 ],
             [ Gauss valuation induced by Valuation on rational function field induced by [ Gauss valuation induced by 3-adic valuation, v(x) = 1/3 ], v(y + 2*x) = 2/3 ]]

        """
        R = G.parent()
        if R.base_ring() is not self.domain():
            raise ValueError("G must be defined over the domain of this valuation")
        if not assume_squarefree and not G.is_squarefree():
            raise ValueError("G must be squarefree")

        from sage.rings.all import infinity
        from gauss_valuation import GaussValuation

        leaves = [ GaussValuation(R, self)]
        while True:
            ef = [ v.E()*v.F() for v in leaves]
            if sum(ef) == G.degree():
                if precision_cap is None or all([v(v.phi())>precision_cap for v in leaves]):
                    return leaves

            expandables = []
            new_leaves = []
            for v in leaves:
                if v(G) is infinity:
                    new_leaves.append(v)
                else:
                    expandables.append(v)
            leaves = new_leaves

            if not expandables:
                return leaves

            for v in expandables:
                leaves.extend(v.mac_lane_step(G))

    def mac_lane_approximant(self, G, valuation, approximants = None):
        r"""
        sage: TODO
        """
        if valuation.constant_valuation() != self:
            raise ValueError

        # Check thet valuation is an approximant for a valuation
        # on domain that extends its restriction to the base field.
        from sage.rings.all import infinity
        if valuation(G) != infinity:
            G_integral = valuation._make_monic_integral(G)
            v = valuation
            while not v.is_gauss_valuation():
                if v(G_integral) <= v._base_valuation(G_integral):
                    raise ValueError("The valuation %r is not an approximant for a valuation on %r since the valuation of %r does not increase in every step"%(valuation, domain, G_integral))
                v = v._base_valuation

        if approximants is None:
            approximants = self.mac_lane_approximants(G)

        assert all(approximant.domain() is valuation.domain() for approximant in approximants)

        greater_approximants = [w for w in approximants if w <= valuation]
        if len(greater_approximants) > 1:
            raise ValueError("valuation %r does not uniquely describe an extension of %r to %r"%(valuation, self, domain))
        if len(greater_approximants) == 1:
            return greater_approximants[0]
        
        smaller_approximants = [w for w in approximants if w >= valuation]
        assert len(smaller_approximants) <= 1
        if len(smaller_approximants) == 0:
            raise ValueError("valuation %r does not describe an extension of %r with respect to %r"%(valuation, valuation.constant_valuation(), G))
        assert len(smaller_approximants) == 1
        return smaller_approximants[0]
