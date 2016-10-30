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

        """
        raise NotImplementedError("No total order for valuation.");

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
    def mac_lane_approximants(self, G, precision_cap=None, assume_squarefree=False):
        r"""
        sage: TODO
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

    def mac_lane_approximant(self, G, valuation):
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

        approximants = self.mac_lane_approximants(G)

        greater_approximants = [w for w in approximants if w <= valuation]
        if len(greater_approximants) > 1:
            raise ValueError("valuation %r does not uniquely describe an extension of %r to %r"%(valuation, self, domain))
        if len(greater_approximants) == 1:
            return greater_approximants[0]
        
        smaller_approximants = [w for w in approximants if w >= valuation]
        assert len(smaller_approximants) <= 1
        if len(smaller_approximants) == 0:
            raise ValueError("valuation %r does not describe an extension of %r to %r"%(valuation, valuation.constant_valuation(), domain))
        assert len(smaller_approximants) == 1
        return smaller_approximants[0]

