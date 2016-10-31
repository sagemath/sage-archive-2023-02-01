# -*- coding: utf-8 -*-
r"""
Valuations which are defined as limits of valuations.

The discrete valuation of a complete field extends uniquely to a finite field
extension. This is not the case anymore for fields which are not complete with
respect to their discrete valuation. In this case, the extensions essentially
correspond to the factors of the defining polynomial of the extension over the
completion. However, these factors only exist over the completion and this
makes it difficult to write down such valuations with a representation of
finite length.

More specifically, let `v` be a discrete valuation on `K` and let `L=K[x]/(G)`
a finite extension thereof. An extension of `v` to `L` can be represented as a
discrete pseudo-valuation `w'` on `K[x]` which sends `G` to infinity.
However, such `w'` might not be described by an :class:`AugmentedValuation`
over a :class:`GaussValuation` anymore. Instead, we may need to write is as a
limit of augmented valuations.

The classes in this module provide the means of writing down such limits and
resulting valuations on quotients.

EXAMPLES:

In this function field, the unique place of ``K`` which corresponds to the zero
point has two extensions to ``L``. The valuations corresponding to these
extensions can only be approximated::

    sage: from mac_lane import * # optional: standalone
    sage: K.<x> = FunctionField(QQ)
    sage: R.<y> = K[]
    sage: L.<y> = K.extension(y^2 - x)

    sage: v = FunctionFieldValuation(K, 1)
    sage: w = v.extensions(L); w
    [[ (x - 1)-adic valuation, v(y - 1) = 1 ]-adic valuation,
     [ (x - 1)-adic valuation, v(y + 1) = 1 ]-adic valuation]

The same phenomenon can be observed for valuations on number fields::

    sage: K = QQ
    sage: R.<t> = K[]
    sage: L.<t> = K.extension(t^2 + 1)
    sage: v = pAdicValuation(QQ, 5)
    sage: w = v.extensions(L); w
    [[ 5-adic valuation, v(t + 2) = 1 ]-adic valuation,
     [ 5-adic valuation, v(t + 3) = 1 ]-adic valuation]

.. NOTE::

    We often rely on approximations of valuations even if we could represent the
    valuation without using a limit. This is done to improve performance as many
    computations already can be done correctly with an approximation::

    sage: K.<x> = FunctionField(QQ)
    sage: R.<y> = K[]
    sage: L.<y> = K.extension(y^2 - x)

    sage: v = FunctionFieldValuation(K, 1/x)
    sage: w = v.extension(L); w
    Valuation at the infinite place
    sage: w._base_valuation._improve_approximation()
    sage: w._base_valuation._approximation
    [ Gauss valuation induced by Valuation at the infinite place, v(y) = 1/2, v(y^2 - 1/x) = +Infinity ]

AUTHORS:

- Julian Rueth (2016-10-19): initial version

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

from sage.misc.abstract_method import abstract_method
from valuation import DiscretePseudoValuation, InfiniteDiscretePseudoValuation, DiscreteValuation
from sage.structure.factory import UniqueFactory

class LimitValuationFactory(UniqueFactory):
    def create_key(self, base_valuation, G):
        return base_valuation, G

    def create_object(self, version, key):
        base_valuation, G = key
        from valuation_space import DiscretePseudoValuationSpace
        parent = DiscretePseudoValuationSpace(base_valuation.domain())
        return parent.__make_element_class__(MacLaneLimitValuation)(parent, base_valuation, G)

LimitValuation = LimitValuationFactory("LimitValuation")

class LimitValuation_generic(DiscretePseudoValuation):
    r"""
    Base class for limit valuations.

    A limit valuation is realized as an approximation of a valuation and means
    to improve that approximation when necessary.

    EXAMPLES::

        sage: from mac_lane import * # optional: standalone
        sage: K.<x> = FunctionField(QQ)
        sage: R.<y> = K[]
        sage: L.<y> = K.extension(y^2 - x)

        sage: v = FunctionFieldValuation(K, 0)
        sage: w = v.extension(L)
        sage: w._base_valuation
        [ Gauss valuation induced by (x)-adic valuation, v(y) = 1/2 , … ]

    The currently used approximation can be found in the ``_approximation``
    field::

        sage: w._base_valuation._approximation
        [ Gauss valuation induced by (x)-adic valuation, v(y) = 1/2 ]

    TESTS::

        sage: isinstance(w._base_valuation, LimitValuation_generic)
        True
        sage: TestSuite(w._base_valuation).run()

    """
    def __init__(self, parent, approximation):
        r"""
        TESTS::
    
            sage: from mac_lane import * # optional: standalone
            sage: R.<x> = QQ[]
            sage: K.<i> = QQ.extension(x^2 + 1)
            sage: v = pAdicValuation(K, 2)
            sage: isinstance(v._base_valuation, LimitValuation_generic)
            True
    
        """
        DiscretePseudoValuation.__init__(self, parent)
        self._initial_approximation = approximation # kept for consistent printing and equality tests
        self._approximation = approximation

    def reduce(self, f):
        r"""
        Return the reduction of ``f`` as an element of :meth:`residue_ring`.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - (x - 1))

            sage: v = FunctionFieldValuation(K, 0)
            sage: w = v.extension(L)
            sage: w.reduce(y) # indirect doctest
            u1

        """
        f = self.domain().coerce(f)
        self._improve_approximation_for_reduce(f)
        F = self._approximation.reduce(f)
        return self.residue_ring()(F)

    def _call_(self, f):
        r"""
        Return the valuation of ``f``.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)

            sage: v = FunctionFieldValuation(K, 0)
            sage: w = v.extension(L)
            sage: w(y) # indirect doctest
            1/2

        """
        self._improve_approximation_for_call(f)
        return self._approximation(f)

    @abstract_method
    def _improve_approximation_for_reduce(self, f):
        r"""
        Replace our approximation with a sufficiently precise approximation to
        correctly compute the reduction of ``f``.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - (x - 1337))

        For the unique extension over the place at 1337, the initial
        approximation is sufficient to compute the reduction of ``y``::

            sage: v = FunctionFieldValuation(K, 1337)
            sage: w = v.extension(L)
            sage: u = w._base_valuation
            sage: u._approximation
            [ Gauss valuation induced by (x - 1337)-adic valuation, v(y) = 1/2 ]
            sage: w.reduce(y)
            0
            sage: u._approximation
            [ Gauss valuation induced by (x - 1337)-adic valuation, v(y) = 1/2 ]

        However, at a place over 1341, the initial approximation is not sufficient
        for some values (note that 1341-1337 is a square)::

            sage: v = FunctionFieldValuation(K, 1341)
            sage: w = v.extensions(L)[0]
            sage: u = w._base_valuation
            sage: u._approximation
            [ Gauss valuation induced by (x - 1341)-adic valuation, v(y - 2) = 1 ]
            sage: w.reduce((y - 2) / (x - 1341)) # indirect doctest
            1/4
            sage: u._approximation
            [ Gauss valuation induced by (x - 1341)-adic valuation, v(y - 1/4*x + 1333/4) = 2 ]
            sage: w.reduce((y - 1/4*x + 1333/4) / (x - 1341)^2) # indirect doctest
            -1/64
            sage: u._approximation
            [ Gauss valuation induced by (x - 1341)-adic valuation, v(y + 1/64*x^2 - 1349/32*x + 1819609/64) = 3 ]

        """

    @abstract_method
    def _improve_approximation_for_call(self, f):
        r"""
        Replace our approximation with a sufficiently precise approximation to
        correctly compute the valuation of ``f``.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - (x - 23))

        For the unique extension over the place at 23, the initial
        approximation is sufficient to compute all valuations::

            sage: v = FunctionFieldValuation(K, 23)
            sage: w = v.extension(L)
            sage: u = w._base_valuation
            sage: u._approximation
            [ Gauss valuation induced by (x - 23)-adic valuation, v(y) = 1/2 ]
            sage: w(x - 23)
            1
            sage: u._approximation
            [ Gauss valuation induced by (x - 23)-adic valuation, v(y) = 1/2 ]

        However, due to performance reasons, sometimes we improve the
        approximation though it would not have been necessary (performing the
        improvement step is faster in this case than checking whether the
        approximation is sufficient)::

            sage: w(y) # indirect doctest
            1/2
            sage: u._approximation
            [ Gauss valuation induced by (x - 23)-adic valuation, v(y) = 1/2, v(y^2 - x + 23) = +Infinity ]

        """

    def _repr_(self):
        r"""
        Return a printable representation of this valuation.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: K = QQ
            sage: R.<t> = K[]
            sage: L.<t> = K.extension(t^2 + 1)
            sage: v = pAdicValuation(QQ, 2)
            sage: w = v.extension(L)
            sage: w._base_valuation # indirect doctest
            [ Gauss valuation induced by 2-adic valuation, v(t + 1) = 1/2 , … ]

        """
        from sage.rings.all import infinity
        from augmented_valuation import AugmentedValuation_generic
        if self._initial_approximation(self._G) < infinity:
            if isinstance(self._initial_approximation, AugmentedValuation_generic):
                return repr(self._initial_approximation)[:-1] + ", … ]"
        return repr(self._initial_approximation)


class MacLaneLimitValuation(LimitValuation_generic, InfiniteDiscretePseudoValuation):
    r"""
    A limit valuation that is a pseudo-valuation on polynomial ring `K[x]`
    which sends an irreducible polynomial `G` to infinity.

    This uses the MacLane algorithm to compute the next element in the limit.

    It starts from a first valuation ``approximation`` whose uniformizer must
    be a uniformizer of the limit and whose residue field must contain the
    residue field of the limit.

    EXAMPLES::

        sage: from mac_lane import * # optional: standalone
        sage: R.<x> = QQ[]
        sage: K.<i> = QQ.extension(x^2 + 1)

        sage: v = pAdicValuation(K, 2)
        sage: u = v._base_valuation; u
        [ Gauss valuation induced by 2-adic valuation, v(x + 1) = 1/2 , … ]

    """
    def __init__(self, parent, approximation, G):
        r"""
        TESTS::
    
            sage: from mac_lane import * # optional: standalone
            sage: R.<x> = QQ[]
            sage: K.<i> = QQ.extension(x^2 + 1)
            sage: v = pAdicValuation(K, 2)
            sage: u = v._base_valuation
            sage: isinstance(u, MacLaneLimitValuation)
            True
    
        """
        LimitValuation_generic.__init__(self, parent, approximation)
        self._G = G

    def lift(self, F):
        r"""
        Return a lift of ``F`` from the :meth:`residue_ring` to the
        :meth:`domain` of this valuatiion.

        EXAMPLES;;

            sage: from mac_lane import * # optional: standalone
            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^4 - x^2 - 2*x - 1)

            sage: v = FunctionFieldValuation(K, 1)
            sage: w = v.extensions(L)[0]; w
            [ (x - 1)-adic valuation, v(y^2 - 2) = 1 ]-adic valuation
            sage: s = w.reduce(y); s
            u1
            sage: w.lift(s) # indirect doctest
            y

        """
        F = self.residue_ring().coerce(F)
        return self._approximation.lift(F)

    def uniformizer(self):
        r"""
        Return a uniformizing element for this valuation.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)

            sage: v = FunctionFieldValuation(K, 0)
            sage: w = v.extension(L)
            sage: w.uniformizer() # indirect doctest
            y

        """
        return self._approximation.uniformizer()

    def _improve_approximation(self):
        r"""
        Perform one step of the Mac Lane algorithm to improve our approximation.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: K = QQ
            sage: R.<t> = K[]
            sage: L.<t> = K.extension(t^2 + 1)
            sage: v = pAdicValuation(QQ, 2)
            sage: w = v.extension(L)
            sage: u = w._base_valuation
            sage: u._approximation
            [ Gauss valuation induced by 2-adic valuation, v(t + 1) = 1/2 ]
            sage: u._improve_approximation()
            sage: u._approximation
            [ Gauss valuation induced by 2-adic valuation, v(t + 1) = 1/2, v(t^2 + 1) = +Infinity ]

        This method has no effect, if the approximation is already an infinite
        valuation::

            sage: u._improve_approximation()
            sage: u._approximation
            [ Gauss valuation induced by 2-adic valuation, v(t + 1) = 1/2, v(t^2 + 1) = +Infinity ]

        """
        from sage.rings.all import infinity
        if self._approximation(self._G) == infinity:
            # an infinite valuation can not be improved further
            return

        approximations = self._approximation.mac_lane_step(self._G)
        assert(len(approximations)==1)
        self._approximation = approximations[0]

    def _improve_approximation_for_call(self, f):
        r"""
        Replace our approximation with a sufficiently precise approximation to
        correctly compute the valuation of ``f``.

        EXAMPLES:

        In this examples, the approximation is increased unnecessarily. The
        first approximation would have been precise enough to compute the
        valuation of ``t + 2``. However, it is faster to improve the
        approximation (perform one step of the Mac Lane algorithm) than to
        check for this::

            sage: from mac_lane import * # optional: standalone
            sage: K = QQ
            sage: R.<t> = K[]
            sage: L.<t> = K.extension(t^2 + 1)
            sage: v = pAdicValuation(QQ, 5)
            sage: w = v.extensions(L)[0]
            sage: u = w._base_valuation
            sage: u._approximation
            [ Gauss valuation induced by 5-adic valuation, v(t + 2) = 1 ]
            sage: w(t + 2) # indirect doctest
            1
            sage: u._approximation
            [ Gauss valuation induced by 5-adic valuation, v(t + 7) = 2 ]

        ALGORITHM:

            Write `L=K[x]/(G)` and consider `g` a representative of the class
            of ``f`` in `K[x]` (of minimal degree.) Write `v` for
            ``self._approximation` and `\phi` for the last key polynomial of
            `v`. With repeated quotient and remainder `g` has a unique
            expansion as `g=\sum a_i\phi^i`. Suppose that `g` is an
            equivalence-unit with respect to ``self._approximation``, i.e.,
            `v(a_0) < v(a_i\phi^i)` for all `i\ne 0`. If we denote the limit
            valuation as `w`, then `v(a_i\phi^i)=w(a_i\phi^i)` since the
            valuation of key polynomials does not change during augmentations
            (Theorem 6.4 in [ML1936'].) By the strict triangle inequality,
            `w(g)=v(g)`.
            Note that any `g` which  is not in `(G)` is an equivalence-unit
            after finitely many steps of the Mac Lane algorithm. Indeed,
            otherwise the valuation of `g` would be infinite (follows from
            Theorem 5.1 in [ML1936']) since the valuation of the key
            polynomials increases.

        """
        from sage.rings.all import infinity
        if self._approximation._mu == infinity:
            # an infinite valuation can not be improved further
            return

        if f == 0:
            # zero always has infinite valuation (actually, this might
            # not be desirable for inexact zero elements with leading
            # zero coefficients.)
            return

        if self._approximation.is_equivalence_unit(f):
            # see ALGORITHM above
            return

        self._improve_approximation()
        return self._improve_approximation_for_call(f)

    def _improve_approximation_for_reduce(self, f):
        r"""
        Replace our approximation with a sufficiently precise approximation to
        correctly compute the reduction of ``f``.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: K = QQ
            sage: R.<t> = K[]
            sage: L.<t> = K.extension(t^2 + 1)
            sage: v = pAdicValuation(QQ, 13)
            sage: w = v.extensions(L)[0]
            sage: u = w._base_valuation
            sage: u._approximation
            [ Gauss valuation induced by 13-adic valuation, v(t + 5) = 1 ]
            sage: w.reduce((t + 5) / 13) # indirect doctest
            8
            sage: u._approximation
            [ Gauss valuation induced by 13-adic valuation, v(t + 70) = 2 ]

        ALGORITHM:

            The reduction produced by the approximation is correct for an
            equivalence-unit, see :meth:`_improve_approximation_for_call`.

        """
        if self._approximation(f) > 0:
            return
        self._improve_approximation_for_call(f)

    def residue_ring(self):
        r"""
        Return the residue field of this valuation.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: K = QQ
            sage: R.<t> = K[]
            sage: L.<t> = K.extension(t^2 + 1)
            sage: v = pAdicValuation(QQ, 2)
            sage: w = v.extension(L)
            sage: w.residue_ring()
            Finite Field of size 2

        """
        R = self._approximation.residue_ring()
        from sage.categories.fields import Fields
        if R in Fields():
            # the approximation ends in v(phi)=infty
            return R
        else:
            from sage.rings.polynomial.polynomial_ring import is_PolynomialRing
            assert(is_PolynomialRing(R))
            return R.base_ring()

class FiniteExtensionFromInfiniteValuation(DiscreteValuation):
    r"""
    A valuation on a quotient of the form `L=K[x]/(G)` with an irreducible `G`
    which is internally backed by a pseudo-valuations on `K[x]` which sends `G`
    to infinity.

    INPUT:

    - ``parent`` -- the containing valuation space (usually the space of
      discrete valuations on `L`)

    - ``base_valuation`` -- an infinite valuation on `K[x]` which takes `G` to
      infinity.

    EXAMPLES::

        sage: from mac_lane import * # optional: standalone
        sage: K.<x> = FunctionField(QQ)
        sage: R.<y> = K[]
        sage: L.<y> = K.extension(y^2 - x)

        sage: v = FunctionFieldValuation(K, 0)
        sage: w = v.extension(L); w
        (x)-adic valuation

    TESTS::

        sage: TestSuite(w).run()

    """
    def __init__(self, parent, base_valuation):
        r"""
        TESTS::

            sage: from mac_lane import * # optional: standalone
            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x^2 + 1)

            sage: v = FunctionFieldValuation(K, 0)
            sage: w = v.extension(L); w
            (x)-adic valuation
            sage: isinstance(w, FiniteExtensionFromInfiniteValuation)
            True

        """
        DiscretePseudoValuation.__init__(self, parent)
        self._base_valuation = base_valuation 

    def _eq_(self, other):
        r"""
        Return whether this valuation is indistinguishable from ``other``.

        EXAMPLES:

            sage: from mac_lane import * # optional: standalone
            sage: K = QQ
            sage: R.<t> = K[]
            sage: L.<t> = K.extension(t^2 + 1)
            sage: v = pAdicValuation(QQ, 2)
            sage: w = v.extension(L)
            sage: ww = v.extension(L)
            sage: w == ww # indirect doctest
            True

        """
        return isinstance(other, FiniteExtensionFromInfiniteValuation) and self._base_valuation == other._base_valuation

    def _repr_(self):
        r"""
        Return a printable representation of this valuation.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: K = QQ
            sage: R.<t> = K[]
            sage: L.<t> = K.extension(t^2 + 1)
            sage: v = pAdicValuation(QQ, 2)
            sage: v.extension(L) # indirect doctest
            2-adic valuation

        """
        return repr(self._base_valuation)

    def residue_ring(self):
        r"""
        Return the residue field of this valuation.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: K = QQ
            sage: R.<t> = K[]
            sage: L.<t> = K.extension(t^2 + 1)
            sage: v = pAdicValuation(QQ, 2)
            sage: v.extension(L).residue_ring()
            Finite Field of size 2

        """
        return self._base_valuation.residue_ring()

    def uniformizer(self):
        r"""
        Return a uniformizing element of this valuation.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: K = QQ
            sage: R.<t> = K[]
            sage: L.<t> = K.extension(t^2 + 1)
            sage: v = pAdicValuation(QQ, 2)
            sage: v.extension(L).uniformizer()
            t + 1

        """
        return self._from_base_domain(self._base_valuation.uniformizer())

    def _to_base_domain(self, f):
        r"""
        Return ``f`` as an element in the domain of ``_base_valuation``.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)

            sage: v = FunctionFieldValuation(K, 0)
            sage: w = v.extensions(L)[0]
            sage: w._to_base_domain(y).parent()
            Univariate Polynomial Ring in y over Rational function field in x over Rational Field

        """
        return self._base_valuation.domain().coerce(f)

    def _from_base_domain(self, f):
        r"""
        Return ``f`` as an element in the domain of this valuation.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)

            sage: v = FunctionFieldValuation(K, 0)
            sage: w = v.extension(L)
            sage: w._from_base_domain(w._base_valuation.domain().gen()).parent()
            Function field in y defined by y^2 - x

        """
        return self.domain().coerce(f)

    def _call_(self, f):
        r"""
        Evaluate this valuation at ``f``.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)

            sage: v = FunctionFieldValuation(K, 0)
            sage: w = v.extension(L)
            sage: w(y) # indirect doctest
            1/2

        """
        return self._base_valuation(self._to_base_domain(f))

    def reduce(self, f):
        r"""
        Return the reduction of ``f`` in the :meth:`residue_field` of this valuation.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - (x - 2))

            sage: v = FunctionFieldValuation(K, 0)
            sage: w = v.extension(L)
            sage: w.reduce(y)
            u1

        """
        return self._base_valuation.reduce(self._to_base_domain(f))

    def lift(self, F):
        r"""
        Lift ``F`` from the :meth;`residue_field` of this valuation into its
        domain.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)

            sage: v = FunctionFieldValuation(K, 2)
            sage: w = v.extension(L)
            sage: w.lift(w.residue_field().gen())
            y

        """
        F = self.residue_ring().coerce(F)
        F = self._base_valuation.residue_ring().coerce(F)
        f = self._base_valuation.lift(F)
        return self._from_base_domain(f)

class FiniteExtensionFromLimitValuation(FiniteExtensionFromInfiniteValuation):
    r"""
    An extension of a valuation on a finite field extensions `L=K[x]/(G)` which
    is induced by an infinite limit valuation on `K[x]`.

    EXAMPLES::

        sage: from mac_lane import * # optional: standalone
        sage: K.<x> = FunctionField(QQ)
        sage: R.<y> = K[]
        sage: L.<y> = K.extension(y^2 - x)
        sage: v = FunctionFieldValuation(K, 1)
        sage: w = v.extensions(L); w
        [[ (x - 1)-adic valuation, v(y - 1) = 1 ]-adic valuation,
         [ (x - 1)-adic valuation, v(y + 1) = 1 ]-adic valuation]

    TESTS::

        sage: TestSuite(w[0]).run()
        sage: TestSuite(w[1]).run()

    """
    def __init__(self, parent, approximant, G, approximants):
        r"""
        EXAMPLES:

        Note that this implementation is also used when the underlying limit is
        only taken over a finite sequence of valuations::

            sage: from mac_lane import * # optional: standalone
            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)
            sage: v = FunctionFieldValuation(K, 2)
            sage: w = v.extension(L); w
            (x - 2)-adic valuation
            sage: isinstance(w, FiniteExtensionFromLimitValuation)
            True

        """
        # keep track of all extensions to this field extension so we can print
        # this valuation nicely, dropping any unnecessary information
        self._approximants = approximants

        from valuation_space import DiscretePseudoValuationSpace
        from limit_valuation import LimitValuation
        limit = LimitValuation(approximant, G)
        FiniteExtensionFromInfiniteValuation.__init__(self, parent, limit)

    def _repr_(self):
        """
        Return a printable representation of this valuation.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: pAdicValuation(GaussianIntegers().fraction_field(), 2) # indirect doctest
            2-adic valuation

        """
        if isinstance(self._base_valuation, MacLaneLimitValuation):
            # print the minimal information that singles out this valuation from all approximants
            assert(self._base_valuation._initial_approximation in self._approximants)
            approximants = [v._augmentations() for v in self._approximants]
            augmentations = self._base_valuation._approximation._augmentations()
            unique_approximant = self._base_valuation._initial_approximation
            for l in range(len(augmentations)):
                if len([a for a in approximants if a[:l+1] == augmentations[:l+1]]) == 1:
                    unique_approximant = augmentations[:l+1]
                    break
            if unique_approximant[0].is_gauss_valuation():
                unique_approximant[0] = unique_approximant[0].constant_valuation()
            if len(unique_approximant) == 1:
                return repr(unique_approximant[0])
            from augmented_valuation import AugmentedValuation_generic
            return "[ %s ]-adic valuation"%(", ".join("v(%r) = %r"%(v._phi, v._mu) if (isinstance(v, AugmentedValuation_generic) and v.domain() == self._base_valuation.domain()) else repr(v) for v in unique_approximant))
        return "%s-adic valuation"%(self._base_valuation)

