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
        from sage.rings.all import infinity
        if self(f) == infinity:
            return self(g) == infinity

        return self(f-g) > self(f)

    def __hash__(self):
        r"""
        The hash value of this valuation.

        We redirect to :meth:`_hash_`, so that subclasses can only override
        :meth:`_hash_` and :meth:`_eq_` if they want to provide a different
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
        ``id``. This works fine for objects which use unique representation.

        Note that the vast majority of valuations come out of a
        :class:`sage.structure.factory.UniqueFactory` and therefore override
        our implementation of :meth:`__hash__` and :meth:`__eq__`.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: v = pAdicValuation(QQ, 2)
            sage: hash(v) == hash(v) # indirect doctest
            True
            
        """
        return id(self)

    def _cmp_(self, other):
        r"""
        Compare this element to ``other``.

        Since there is no reasonable total order on valuations, this method
        just throws an exception.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
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

            sage: from mac_lane import * # optional: standalone
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
            return self._le_(other)
        if op == 2: # ==
            return self._eq_(other)
        if op == 3: # !=
            return not self == other
        if op == 5: # >=
            return self._ge_(other)
        raise NotImplementedError("Operator not implemented for this valuation.")

    def _eq_(self, other):
        r"""
        Return whether this valuation and ``other`` are indistinguishable.

        We override the strange default provided by
        ``sage.categories.marphism.Morphism`` here and implement equality by
        ``id``. This is the right behaviour in many cases.

        Note that the vast majority of valuations come out of a
        :class:`sage.structure.factory.UniqueFactory` and therefore override
        our implementation of :meth:`__hash__` and :meth:`__eq__`.

        When overriding this method, you can assume that ``other`` is a
        (pseudo-)valuation on the same domain.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: v = TrivialValuation(QQ)
            sage: v == v
            True

        """
        return self is other

    def _le_(self, other):
        r"""
        Return whether this valuation is less than or equal to ``other``
        pointwise.

        When overriding this method, you can assume that ``other`` is a
        (pseudo-)valuation on the same domain.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: v = TrivialValuation(QQ)
            sage: w = pAdicValuation(QQ, 2)
            sage: v <= w
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
        return other >= self

    def _ge_(self, other):
        r"""
        Return whether this valuation is greater than or equal to ``other``
        pointwise.

        When overriding this method, you can assume that ``other`` is a
        (pseudo-)valuation on the same domain.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: v = TrivialValuation(QQ)
            sage: w = pAdicValuation(QQ, 2)
            sage: v >= w
            False

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
        if self == other: return True
        raise NotImplementedError("Operator not implemented for this valuation.")

    # Remove the default implementation of Map.__reduce__ that does not play
    # nice with factories (a factory, does not override Map.__reduce__ because
    # it is not the generic reduce of object) and that does not match equality
    # by id.
    __reduce__ = object.__reduce__

class InfiniteDiscretePseudoValuation(DiscretePseudoValuation):
    r"""
    Abstract base class for infinite discrete pseudo-valuations, i.e., discrete
    pseudo-valuations which are not discrete valuations.

    EXAMPLES::

        sage: from mac_lane import * # optional: standalone
        sage: v = pAdicValuation(QQ, 2)
        sage: R.<x> = QQ[]
        sage: v = GaussValuation(R, v)
        sage: w = v.augmentation(x, infinity); w # indirect doctest
        [ Gauss valuation induced by 2-adic valuation, v(x) = +Infinity ]
    
    TESTS::

        sage: isinstance(w, InfiniteDiscretePseudoValuation)
        True
        sage: TestSuite(w).run()

    """
    def is_discrete_valuation(self):
        r"""
        Return whether this valuation is a discrete valuation.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: v = pAdicValuation(QQ, 2)
            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, v)
            sage: w = v.augmentation(x, infinity)
            sage: w.is_discrete_valuation()
            False

        """
        return False

class DiscreteValuation(DiscretePseudoValuation):
    r"""
    Abstract base class for discrete valuations.

    EXAMPLES::

        sage: from mac_lane import * # optional: standalone
        sage: v = pAdicValuation(QQ, 2)
        sage: R.<x> = QQ[]
        sage: v = GaussValuation(R, v)
        sage: w = v.augmentation(x, 1337); w # indirect doctest
        [ Gauss valuation induced by 2-adic valuation, v(x) = 1337 ]

    TESTS::

        sage: isinstance(w, DiscreteValuation)
        True
        sage: TestSuite(w).run()

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
        Return approximants on `K[x]` for the extensions of this valuation to
        `L=K[x]/(G)`.

        If `G` is an irreducible polynomial, then this corresponds to the
        extensions of this valuation to the completion of `L`.

        INPUT:

        - ``G`` -- a square-free polynomial defined over a univariate
          polynomial ring over the :meth:`domain` of this valuation.

        - ``precision_cap`` -- a number, infinity, or ``None`` (default:
          ``None``); the approximants are always determined such that they are in
          one-to-one correspondance to the extensions of this valuation to `L`
          and such that the approximants have the ramification index and
          residual degree of these extensions.
          If ``precision_cap`` is not ``None``, then the approximants are
          determined such that they last key polynomial also has valuation at
          least ``precision_cap``.

        - ``assume_squarefree`` -- a boolean (default: ``False``), whether or
          not to assume that ``G`` is squarefree.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: v = pAdicValuation(QQ, 2)
            sage: R.<x> = QQ[]
            sage: v.mac_lane_approximants(x^2 + 1)
            [[ Gauss valuation induced by 2-adic valuation, v(x + 1) = 1/2 ]]
            sage: v.mac_lane_approximants(x^2 + 1, precision_cap=infinity)
            [[ Gauss valuation induced by 2-adic valuation, v(x + 1) = 1/2, v(x^2 + 1) = +Infinity ]]
            sage: v.mac_lane_approximants(x^2 + x + 1)
            [[ Gauss valuation induced by 2-adic valuation, v(x^2 + x + 1) = +Infinity ]]

        Note that ``G`` does not need to be irreducible. Here, we detect a
        factor of the polynomial `x + 1` and an approximate factor `x + 1`
        (which is an approximation to `x - 1`)::

            sage: v.mac_lane_approximants(x^2 - 1)
            [[ Gauss valuation induced by 2-adic valuation, v(x + 1) = +Infinity ],
             [ Gauss valuation induced by 2-adic valuation, v(x + 3) = 2 ]]

        However, it needs to be squarefree::

            sage: v.mac_lane_approximants(x^2)
            Traceback (most recent call last):
            ...
            ValueError: G must be squarefree

        TESTS:

        Some difficult cases provided by Mark van Hoeij::

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

        Over a complete base field::

            sage: k=Qp(2,10)
            sage: v = pAdicValuation(k)

            sage: R.<x>=k[]
            sage: G = x
            sage: v.mac_lane_approximants(G)
            [Gauss valuation induced by 2-adic valuation]
            sage: v.mac_lane_approximants(G, precision_cap = infinity)
            [[ Gauss valuation induced by 2-adic valuation, v((1 + O(2^10))*x) = +Infinity ]]

            sage: G = x^2 + 1
            sage: v.mac_lane_approximants(G) # optional: integrated
            [[ Gauss valuation induced by 2-adic valuation, v((1 + O(2^10))*x + (1 + O(2^10))) = 1/2 ]]
            sage: v.mac_lane_approximants(G, precision_cap = infinity) # optional: integrated
            [[ Gauss valuation induced by 2-adic valuation, v((1 + O(2^10))*x + (1 + O(2^10))) = 1/2, v((1 + O(2^10))*x^2 + (1 + O(2^10))) = +Infinity ]]

            sage: G = x^4 + 2*x^3 + 2*x^2 - 2*x + 2
            sage: v.mac_lane_approximants(G)
            [[ Gauss valuation induced by 2-adic valuation, v((1 + O(2^10))*x) = 1/4 ]]
            sage: v.mac_lane_approximants(G,infinity)
            [[ Gauss valuation induced by 2-adic valuation, v((1 + O(2^10))*x) = 1/4, v((1 + O(2^10))*x^4 + (2 + O(2^11))*x^3 + (2 + O(2^11))*x^2 + (2 + 2^2 + 2^3 + 2^4 + 2^5 + 2^6 + 2^7 + 2^8 + 2^9 + 2^10 + O(2^11))*x + (2 + O(2^11))) = +Infinity ]]

        The factorization of primes in the Gaussian integers can be read off
        the Mac Lane approximants::

            sage: v0 = pAdicValuation(QQ, 2)
            sage: R.<x> = QQ[]
            sage: G = x^2 + 1
            sage: v0.mac_lane_approximants(G)
            [[ Gauss valuation induced by 2-adic valuation, v(x + 1) = 1/2 ]]

            sage: v0 = pAdicValuation(QQ, 3)
            sage: v0.mac_lane_approximants(G)
            [[ Gauss valuation induced by 3-adic valuation, v(x^2 + 1) = +Infinity ]]

            sage: v0 = pAdicValuation(QQ, 5)
            sage: v0.mac_lane_approximants(G)
            [[ Gauss valuation induced by 5-adic valuation, v(x + 2) = 1 ],
             [ Gauss valuation induced by 5-adic valuation, v(x + 3) = 1 ]]
            sage: v0.mac_lane_approximants(G, precision_cap = 10) # long time
            [[ Gauss valuation induced by 5-adic valuation, v(x + 6139557) = 10 ],
             [ Gauss valuation induced by 5-adic valuation, v(x + 3626068) = 10 ]]

        The same example over the 5-adic numbers. In the quadratic extension
        `\QQ[x]/(x^2+1)`, 5 factors `-(x - 2)(x + 2)`, this behaviour can be
        read off the Mac Lane approximants::

            sage: k=Qp(5,4)
            sage: v = pAdicValuation(k)
            sage: R.<x>=k[]
            sage: G = x^2 + 1
            sage: v1,v2 = v.mac_lane_approximants(G); v1,v2
            ([ Gauss valuation induced by 5-adic valuation, v((1 + O(5^4))*x + (2 + O(5^4))) = 1 ],
             [ Gauss valuation induced by 5-adic valuation, v((1 + O(5^4))*x + (3 + O(5^4))) = 1 ])
            sage: w1, w2 = v.mac_lane_approximants(G,precision_cap=2); w1,w2 # long time
            ([ Gauss valuation induced by 5-adic valuation, v((1 + O(5^4))*x + (2 + 5 + O(5^4))) = 2 ],
             [ Gauss valuation induced by 5-adic valuation, v((1 + O(5^4))*x + (3 + 3*5 + O(5^4))) = 2 ])

        Note how the latter give a better approximation to the factors of `x^2 + 1`::

            sage: v1.phi() * v2.phi() - G # optional: integrated
            (5 + O(5^4))*x + 5 + O(5^4)
            sage: w1.phi() * w2.phi() - G # optional: integrated
            (5^3 + O(5^4))*x + 5^3 + O(5^4)

        In this example, the process stops with a factorization of `x^2 + 1`::

            sage: v.mac_lane_approximants(G, precision_cap=infinity) # long time
            [[ Gauss valuation induced by 5-adic valuation, v((1 + O(5^4))*x + (2 + 5 + 2*5^2 + 5^3 + O(5^4))) = +Infinity ],
             [ Gauss valuation induced by 5-adic valuation, v((1 + O(5^4))*x + (3 + 3*5 + 2*5^2 + 3*5^3 + O(5^4))) = +Infinity ]]

        This obviously cannot happen over the rationals where we only get an
        approximate factorization::

            sage: v = pAdicValuation(QQ, 5)
            sage: R.<x>=QQ[]
            sage: G = x^2 + 1
            sage: v.mac_lane_approximants(G)
            [[ Gauss valuation induced by 5-adic valuation, v(x + 2) = 1 ], [ Gauss valuation induced by 5-adic valuation, v(x + 3) = 1 ]]
            sage: v.mac_lane_approximants(G, precision_cap=5)
            [[ Gauss valuation induced by 5-adic valuation, v(x + 2057) = 5 ],
             [ Gauss valuation induced by 5-adic valuation, v(x + 1068) = 6 ]]

        Initial versions ran into problems with the trivial residue field
        extensions in this case::

            sage: K = Qp(3,20)
            sage: R.<T> = K[]

            sage: alpha = T^3/4
            sage: G = 3^3*T^3*(alpha^4 - alpha)^2 - (4*alpha^3 - 1)^3
            sage: G = G/G.leading_coefficient()
            sage: pAdicValuation(K).mac_lane_approximants(G) # optional: integrated, long time
            [[ Gauss valuation induced by 3-adic valuation, v((1 + O(3^20))*T + (2 + O(3^20))) = 1/9, v((1 + O(3^20))*T^9 + (2*3 + 2*3^2 + O(3^21))*T^8 + (3 + 3^5 + O(3^21))*T^7 + (2*3 + 2*3^2 + 3^3 + 2*3^4 + 2*3^5 + 3^6 + O(3^21))*T^6 + (2*3 + 2*3^2 + 3^4 + 3^6 + 2*3^7 + O(3^21))*T^5 + (3 + 3^2 + 3^3 + 2*3^6 + 2*3^7 + 3^8 + O(3^21))*T^4 + (2*3 + 2*3^2 + 3^3 + 2*3^5 + 2*3^6 + 2*3^7 + 2*3^8 + O(3^21))*T^3 + (2*3 + 2*3^2 + 3^3 + 2*3^4 + 3^5 + 2*3^6 + 2*3^7 + 2*3^8 + O(3^21))*T^2 + (3 + 2*3^2 + 2*3^3 + 2*3^4 + 2*3^7 + 3^8 + O(3^21))*T + (2 + 2*3 + 2*3^2 + 2*3^4 + 2*3^5 + 3^7 + O(3^20))) = 55/27 ]]

        A similar example::

          sage: R.<x> = QQ[]
          sage: v = pAdicValuation(QQ, 3)
          sage: G = (x^3 + 3)^3 - 81
          sage: v.mac_lane_approximants(G) # optional: integrated
          [[ Gauss valuation induced by 3-adic valuation, v(x) = 1/3, v(x^3 + 3*x + 3) = 13/9 ]]

        Another problematic case::

            sage: R.<x> = QQ[] 
            sage: Delta = x^12 + 20*x^11 + 154*x^10 + 664*x^9 + 1873*x^8 + 3808*x^7 + 5980*x^6 + 7560*x^5 + 7799*x^4 + 6508*x^3 + 4290*x^2 + 2224*x + 887 
            sage: K.<theta> = NumberField(x^6 + 108) 
            sage: K.is_galois()
            True
            sage: vK = pAdicValuation(QQ, 2).extension(K)
            sage: vK(2) 
            1 
            sage: vK(theta) 
            1/3
            sage: G=Delta.change_ring(K) 
            sage: V=vK.mac_lane_approximants(G); V # long time
            [[ Gauss valuation induced by 2-adic valuation, v(x + 1) = 1/4, v(x^4 + 4*x^3 + 6*x^2 + 4*x + theta^4 + theta^3 + 1) = 5/3 ],
             [ Gauss valuation induced by 2-adic valuation, v(x + 1) = 1/4, v(x^4 + 4*x^3 + 6*x^2 + 4*x + 1/2*theta^4 + theta^3 - 27*theta + 1) = 5/3 ],
             [ Gauss valuation induced by 2-adic valuation, v(x + 1) = 1/4, v(x^4 + 4*x^3 + 6*x^2 + 4*x + 3/2*theta^4 + theta^3 - 27*theta + 1) = 5/3 ]]

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
                # the ramification indexes and residual degrees are final
                if precision_cap is None or all([v(v.phi()) >= precision_cap for v in leaves]):
                    # the precision to which we have determined the approximants is sufficient
                    if not [v for v in leaves if [w for w in leaves if w != v and w <= v]]:
                        # the approximants do not approximate each other
                        break

            expandables = []
            new_leaves = []
            for v in leaves:
                if v(G) is infinity:
                    new_leaves.append(v)
                else:
                    expandables.append(v)
            leaves = new_leaves

            assert expandables

            for v in expandables:
                leaves.extend(v.mac_lane_step(G))

        return leaves

    def mac_lane_approximant(self, G, valuation, approximants = None):
        r"""
        Return the approximant from :meth:`mac_lane_approximants` for ``G``
        which is approximated by or approximates ``valuation``.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: v = pAdicValuation(QQ, 2)
            sage: R.<x> = QQ[]
            sage: G = x^2 + 1

        We can select an approximant by approximating it::

            sage: w = GaussValuation(R, v).augmentation(x + 1, 1/2)
            sage: v.mac_lane_approximant(G, w)
            [ Gauss valuation induced by 2-adic valuation, v(x + 1) = 1/2 ]

        As long as this is the only matching approximant, the approximation can
        be very coarse::

            sage: w = GaussValuation(R, v)
            sage: v.mac_lane_approximant(G, w)
            [ Gauss valuation induced by 2-adic valuation, v(x + 1) = 1/2 ]

        Or it can be very specific::

            sage: w = GaussValuation(R, v).augmentation(x + 1, 1/2).augmentation(G, infinity)
            sage: v.mac_lane_approximant(G, w)
            [ Gauss valuation induced by 2-adic valuation, v(x + 1) = 1/2 ]

        But it must be an approximation of an approximant::

            sage: w = GaussValuation(R, v).augmentation(x, 1/2)
            sage: v.mac_lane_approximant(G, w)
            Traceback (most recent call last):
            ...
            ValueError: The valuation [ Gauss valuation induced by 2-adic valuation, v(x) = 1/2 ] is not an approximant for a valuation which extends 2-adic valuation with respect to x^2 + 1 since the valuation of x^2 + 1 does not increase in every step

        The ``valuation`` must single out one approximant::

            sage: G = x^2 - 1
            sage: w = GaussValuation(R, v)
            sage: v.mac_lane_approximant(G, w)
            Traceback (most recent call last):
            ...
            ValueError: The valuation Gauss valuation induced by 2-adic valuation does not approximate a unique extension of 2-adic valuation with respect to x^2 - 1

            sage: w = GaussValuation(R, v).augmentation(x + 1, 1)
            sage: v.mac_lane_approximant(G, w)
            Traceback (most recent call last):
            ...
            ValueError: The valuation [ Gauss valuation induced by 2-adic valuation, v(x + 1) = 1 ] does not approximate a unique extension of 2-adic valuation with respect to x^2 - 1

            sage: w = GaussValuation(R, v).augmentation(x + 1, 2)
            sage: v.mac_lane_approximant(G, w)
            [ Gauss valuation induced by 2-adic valuation, v(x + 1) = +Infinity ]

            sage: w = GaussValuation(R, v).augmentation(x + 3, 2)
            sage: v.mac_lane_approximant(G, w)
            [ Gauss valuation induced by 2-adic valuation, v(x + 3) = 2 ]

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
                    raise ValueError("The valuation %r is not an approximant for a valuation which extends %r with respect to %r since the valuation of %r does not increase in every step"%(valuation, self, G, G_integral))
                v = v._base_valuation

        if approximants is None:
            approximants = self.mac_lane_approximants(G)

        assert all(approximant.domain() is valuation.domain() for approximant in approximants)

        greater_approximants = [w for w in approximants if w >= valuation]
        if len(greater_approximants) > 1:
            raise ValueError("The valuation %r does not approximate a unique extension of %r with respect to %r"%(valuation, self, G))
        if len(greater_approximants) == 1:
            return greater_approximants[0]
        
        smaller_approximants = [w for w in approximants if w <= valuation]
        if len(smaller_approximants) > 1:
            raise ValueError("The valuation %r is not approximated by a unique extension of %r with respect to %r"%(valuation, valuation.constant_valuation(), G))
        if len(smaller_approximants) == 0:
            raise ValueError("The valuation %r is not related to an extension of %r with respect to %r"%(valuation, valuation.constant_valuation(), G))
        assert len(smaller_approximants) == 1
        return smaller_approximants[0]

    def _ge_(self, other):
        r"""
        Return whether this valuation is greater than or equal to ``other``
        pointwise.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: v = TrivialValuation(QQ)
            sage: w = pAdicValuation(QQ, 2)
            sage: v >= w
            False

        """
        if other.is_trivial():
            return other.is_discrete_valuation()
        return super(DiscreteValuation, self)._ge_(other)
