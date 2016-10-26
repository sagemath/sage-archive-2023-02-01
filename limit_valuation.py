# -*- coding: utf-8 -*-
r"""
Valuations which are defined as limits of valuations.

The discrete valuation of a complete field extends uniquely to a finite field
extension. This is not the case anymore for fields which are not complete with
respect to their discrete valuation. In this case, the extensions essentially
correspond to the factors of the defining polynomial of the extension over the
completion. However, these factors only exist over the completion and so the
valuation can not written down finitely as a number of
:class:`AugmentedValuation`s but only as a limit thereof.

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
    [[ Gauss valuation induced by (x - 1)-adic valuation, v(y - 1) = 1 , … ],
     [ Gauss valuation induced by (x - 1)-adic valuation, v(y + 1) = 1 , … ]]

The same phenomenon can be observed for valuations on number fields::

    sage: K = QQ
    sage: R.<t> = K[]
    sage: L.<t> = K.extension(t^2 + 1)
    sage: v = pAdicValuation(QQ, 5)
    sage: w = v.extensions(L); w
    [[ Gauss valuation induced by 5-adic valuation, v(t + 2) = 1 , … ],
     [ Gauss valuation induced by 5-adic valuation, v(t + 3) = 1 , … ]]

.. NOTE::

    We often rely on approximations of valuations even if we could represent the
    valuation without using a limit. This is done to improve performance as many
    computations already can be done correctly with an approximation::

    sage: K.<x> = FunctionField(QQ)
    sage: R.<y> = K[]
    sage: L.<y> = K.extension(y^2 - x)

    sage: v = FunctionFieldValuation(K, 1/x)
    sage: w = v.extension(L); w
    [ Gauss valuation induced by Valuation at the infinite place, v(y) = 1/2 , … ]
    sage: w._improve_approximation()
    sage: w._approximation
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
from valuation import DiscretePseudoValuation

class LimitValuation_base(DiscretePseudoValuation):
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

    The currently used approximation can be found in the ``_approximation``
    field::

        sage: w._approximation
        [ Gauss valuation induced by (x)-adic valuation, v(y) = 1/2 ]

    Note that the approximation might be defined on a different domain::

        sage: w.domain()
        Function field in y defined by y^2 - x
        sage: w._approximation.domain()
        Univariate Polynomial Ring in y over Rational function field in x over Rational Field

    The methods :meth:`_to_approximation_domain` and
    :meth:`_from_approximation_domain` move items back and forth between the
    domains of definition.

    TESTS::

        sage: isinstance(w, LimitValuation_base)
        True

        sage: TestSuite(w).run()

    """
    def reduce(self, f):
        r"""
        Return the reduction of ``f`` as an element of :meth:`residue_ring`.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)

            sage: v = FunctionFieldValuation(K, 1)
            sage: w = v.extensions(L)[0]
            sage: w.reduce(y)
            1

        """
        f = self.domain().coerce(f)
        self._improve_approximation_for_reduce(f)
        F = self._approximation.reduce(self._to_approximation_domain(f))
        return self.residue_ring()(F)

    @abstract_method
    def _improve_approximation_for_reduce(self, f):
        r"""
        Replace our approximation with a sufficiently precise approximation to
        correctly compute the reduction of ``f``.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)

        For the unique extension over the place at zero, the initial
        approximation is sufficient to compute the reduction of ``y``::

            sage: v = FunctionFieldValuation(K, 0)
            sage: w = v.extension(L)
            sage: w._approximation
            [ Gauss valuation induced by (x)-adic valuation, v(y) = 1/2 ]
            sage: w.reduce(y)
            0
            sage: w._approximation
            [ Gauss valuation induced by (x)-adic valuation, v(y) = 1/2 ]

        However, at a places over 1, the initial approximation is not
        sufficient for some values::

            sage: v = FunctionFieldValuation(K, 1)
            sage: w = v.extensions(L)[0]
            sage: w._approximation
            [ Gauss valuation induced by (x - 1)-adic valuation, v(y - 1) = 1 ]
            sage: w.reduce((y - 1) / (x - 1)) # indirect doctest
            1/2
            sage: w._approximation
            [ Gauss valuation induced by (x - 1)-adic valuation, v(y - 1/2*x - 1/2) = 2 ]
            sage: w.reduce((y - 1/2*x - 1/2) / (x - 1)^2) # indirect doctest
            -1/8
            sage: w._approximation
            [ Gauss valuation induced by (x - 1)-adic valuation, v(y + 1/8*x^2 - 3/4*x - 3/8) = 3 ]

        """

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
            sage: w(y)
            1/2

        """
        f = self.domain().coerce(f)
        self._improve_approximation_for_call(f)
        return self._approximation(self._to_approximation_domain(f))

    @abstract_method
    def _improve_approximation_for_call(self, f):
        r"""
        Replace our approximation with a sufficiently precise approximation to
        correctly compute the valuation of ``f``.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)

        For the unique extension over the place at zero, the initial
        approximation is sufficient to compute all valuations::

            sage: v = FunctionFieldValuation(K, 0)
            sage: w = v.extension(L)
            sage: w._approximation
            [ Gauss valuation induced by (x)-adic valuation, v(y) = 1/2 ]
            sage: w(x)
            1
            sage: w._approximation
            [ Gauss valuation induced by (x)-adic valuation, v(y) = 1/2 ]

        However, due to performance reasons, sometimes we improve the
        approximation though it would not have been necessary (performing the
        improvement step is faster in this case than checking whether the
        approximation is sufficient)::

            sage: w(y) # indirect doctest
            1/2
            sage: w._approximation
            [ Gauss valuation induced by (x)-adic valuation, v(y) = 1/2, v(y^2 - x) = +Infinity ]

        """

    def _to_approximation_domain(self, f):
        r"""
        Return ``f`` as an element in the domain of ``_approximation``.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)

            sage: v = FunctionFieldValuation(K, 0)
            sage: w = v.extensions(L)[0]
            sage: w._to_approximation_domain(y).parent()
            Univariate Polynomial Ring in y over Rational function field in x over Rational Field

        """
        return self._approximation.domain().coerce(f)

    def _from_approximation_domain(self, f):
        r"""
        Return ``f`` as an element in the domain of this valuation.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)

            sage: v = FunctionFieldValuation(K, 0)
            sage: w = v.extensions(L)[0]
            sage: w._from_approximation_domain(w._approximation.domain().gen()).parent()
            Function field in y defined by y^2 - x

        """
        return self.domain().coerce(f)

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
            [ Gauss valuation induced by (x - 1)-adic valuation, v(y^2 - 2) = 1 , … ]
            sage: u = w.reduce(y); u
            u1
            sage: y == w.lift(u)
            True

        """
        F = self.residue_ring().coerce(F)
        f = self._approximation.lift(F)
        return self._from_approximation_domain(f)

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
            sage: w(w.uniformizer())
            1/2

        """
        return self._from_approximation_domain(self._approximation.uniformizer())

class LimitValuationFiniteExtension(LimitValuation_base):
    r"""
    A limit valuation that comes from a finite simple extension of the form
    `L=K[x]/(G)`.

    Starting from ``approximation``, the approximation used to perform
    computations is improved using the Mac Lane algorithm. The initial value of
    ``approximation`` must already have the final residue field and
    ramification index and it must approximate a unique extension on `L`.

    TESTS::

        sage: from mac_lane import * # optional: standalone
        sage: K.<x> = FunctionField(QQ)
        sage: R.<y> = K[]
        sage: L.<y> = K.extension(y^2 - x)

        sage: v = FunctionFieldValuation(K, 0)
        sage: w = v.extension(L)

        sage: TestSuite(w).run()

    """
    def __init__(self, parent, approximation, G):
        r"""
        TESTS::
    
            sage: from mac_lane import * # optional: standalone
            sage: K = QQ
            sage: R.<t> = K[]
            sage: L.<t> = K.extension(t^2 + 1)
            sage: v = pAdicValuation(QQ, 5)
            sage: w = v.extensions(L); w
            [[ Gauss valuation induced by 5-adic valuation, v(t + 2) = 1 , … ],
             [ Gauss valuation induced by 5-adic valuation, v(t + 3) = 1 , … ]]
    
            sage: isinstance(w[0], LimitValuation_base)
            True
            sage: isinstance(w[1], LimitValuation_base)
            True

        """
        DiscretePseudoValuation.__init__(self, parent)
        self._initial_approximation = approximation # kept for consistent printing
        self._approximation = approximation
        self._G = G

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
            sage: w._approximation
            [ Gauss valuation induced by 2-adic valuation, v(t + 1) = 1/2 ]
            sage: w._improve_approximation()
            sage: w._approximation
            [ Gauss valuation induced by 2-adic valuation, v(t + 1) = 1/2, v(t^2 + 1) = +Infinity ]
            sage: w._improve_approximation()
            sage: w._approximation
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
            sage: w._approximation
            [ Gauss valuation induced by 5-adic valuation, v(t + 2) = 1 ]
            sage: w(t + 2) # indirect doctest
            1
            sage: w._approximation
            [ Gauss valuation induced by 5-adic valuation, v(t + 7) = 2 ]

        """
        from sage.rings.all import infinity
        if self._approximation._mu == infinity:
            return
        if f == 0:
            return
        g = self._to_approximation_domain(f)
        if self._approximation.is_equivalence_unit(g):
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
            sage: v = pAdicValuation(QQ, 5)
            sage: w = v.extensions(L)[0]
            sage: w._approximation
            [ Gauss valuation induced by 5-adic valuation, v(t + 2) = 1 ]
            sage: w.reduce((t + 2) / 5) # indirect doctest
            4
            sage: w._approximation
            [ Gauss valuation induced by 5-adic valuation, v(t + 7) = 2 ]

        """
        g = self._to_approximation_domain(f)
        if self._approximation(g) > 0:
            return
        self._improve_approximation_for_call(f)

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
            [ Gauss valuation induced by 2-adic valuation, v(t + 1) = 1/2 , … ]

        """
        from sage.rings.all import infinity
        from augmented_valuation import AugmentedValuation
        if self._initial_approximation(self._G) < infinity:
            if isinstance(self._initial_approximation, AugmentedValuation):
                return repr(self._initial_approximation)[:-1] + ", … ]"
        return repr(self._initial_approximation)

    def _eq_(self, other):
        r"""
        Return whether this valuation is indistinguishable from ``other``.

        EXAMPLES:

        We deem two valuations indistinguishable if they can not be
        distinguished without considering hidden fields::

            sage: from mac_lane import * # optional: standalone
            sage: K = QQ
            sage: R.<t> = K[]
            sage: L.<t> = K.extension(t^2 + 1)
            sage: v = pAdicValuation(QQ, 2)
            sage: w = v.extension(L)
            sage: ww = v.extension(L)
            sage: w == ww # indirect doctest
            True
            sage: w._improve_approximation()
            sage: w == ww
            True
            sage: w._approximation == ww._approximation
            False

        """
        return isinstance(other, LimitValuationFiniteExtension) and self._G == other._G and self._initial_approximation == other._initial_approximation

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

        TESTS::

            sage: w._improve_approximation()
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
