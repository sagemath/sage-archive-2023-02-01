"""
Affine Lie Algebras

AUTHORS:

- Travis Scrimshaw (2013-05-03): Initial version

EXAMPLES::
"""

#*****************************************************************************
#       Copyright (C) 2013-2017 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.misc import repr_lincomb
from sage.structure.element import RingElement
from sage.categories.lie_algebras import LieAlgebras

from sage.algebras.lie_algebras.lie_algebra import LieAlgebra, FinitelyGeneratedLieAlgebra
from sage.algebras.lie_algebras.lie_algebra_element import UntwistedAffineLieAlgebraElement
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.root_system.cartan_matrix import CartanMatrix
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.categories.cartesian_product import cartesian_product
from sage.rings.integer_ring import ZZ
from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
from sage.sets.family import Family

class AffineLieAlgebra(FinitelyGeneratedLieAlgebra):
    r"""
    An (untwisted) affine Lie algebra.

    Given a finite dimensional simple Lie algebra `\mathfrak{g}` over `R`,
    we construct an affine Lie algebra `\widehat{\mathfrak{g}}^{\prime}` as

    .. MATH::

        \widehat{\mathfrak{g}}^{\prime} = \left( \mathfrak{g} \otimes
        R[t, t^{-1}] \right) \oplus R c,

    where `c` is the canonical central element and `R[t, t^{-1}]` is the
    Laurent polynomial ring over `R`. We define the Lie bracket as

    .. MATH::

        [a \otimes t^n + \alpha c, b \otimes t^m + \beta c] =
        [a, b] \otimes t^{n+m} + \delta_{n+m,0} ( a | b ) n c

    where `( a | b )` is the Killing form on `\mathfrak{g}`.

    There is a canonical derivative on `\widehat{\mathfrak{g}}^{\prime}`
    which is known as the *Lie derivative* and is denoted by `\delta`.
    The Lie derivative is defined as

    .. MATH::

        \delta(a \otimes t^m + \alpha c) = a \otimes m t^m,

    or equivalently by `\delta = t \frac{d}{dt}`.

    We can form the affine Kac-Moody algebra `\widehat{\mathfrak{g}}`
    by adding the additional generator `d` such that `[d, x] = \delta(x)`
    where `\delta` is the Lie derivative. We note that the derived subalgebra
    of the Kac-Moody algebra is the affine Lie algebra.

    .. NOTE::

        Our terminology is following that of :wikipedia:`Affine_Lie_algebra`.
    """
    @staticmethod
    def __classcall_private__(cls, arg0, kac_moody=True, cartan_type=None):
        """
        Parse input to ensure a unique representation.

        INPUT:

        - ``arg0`` -- a simple Lie algebra or a base ring
        - ``cartan_type`` -- a Cartan type

        EXAMPLES::
        """
        if isinstance(arg0, LieAlgebra):
            ct = arg0.cartan_type()
            if not ct.is_finite():
                raise ValueError("the base Lie algebra is not simple")
            cartan_type = ct.affine()
            g = arg0
        else:
            # arg0 is the base ring
            cartan_type = CartanType(cartan_type)
            if not cartan_type.is_affine():
                raise ValueError("the Cartan type must be affine")
            g = LieAlgebra(arg0, cartan_type=cartan_type.classical())

        if not cartan_type.is_untwisted_affine():
            raise NotImplementedError("only currently implemented for untwisted affine types")
        return super(AffineLieAlgebra, cls).__classcall__(cls, g, kac_moody)

    def __init__(self, g, kac_moody):
        """
        Initalize ``self``.
        """
        self._g = g
        self._cartan_type = g.cartan_type()
        R = g.base_ring()
        names = list(g.variable_names()) + ['e0', 'f0', 'c']

        if kac_moody:
            names += ['delta']
        self._kac_moody = kac_moody

        names = tuple(names)
        cat = LieAlgebras(R).WithBasis()
        FinitelyGeneratedLieAlgebra.__init__(self, R, names, names, category=cat)

    def _repr_(self):
        """
        Return a string representation of ``self``.
        """
        base = "Affine "
        rep = repr(self._g)
        if self._kac_moody:
            old_len = len(rep)
            rep = rep.replace("Lie", "Kac-Moody")
            if len(rep) == old_len: # We did not replace anything
                base += "Kac-Moody "
        return base + rep

    @cached_method
    def basis(self):
        """
        Return the basis of ``self``.
        """
        K = cartesian_product([self._g.basis().keys(), ZZ])
        if self._kac_moody:
            keys = DisjointUnionEnumeratedSets([('c',), ('delta',), K])
        else:
            keys = DisjointUnionEnumeratedSets([('c',), K])
        return Family(keys, self.monomial)

    def derived_subalgebra(self):
        """
        Return the derived subalgebra of ``self``.
        """
        if self._kac_moody:
            return AffineLieAlgebra(self._g, False)
        raise NotImplementedError # I think this is self...

    def cartan_type(self):
        """
        Return the Cartan type of ``self``.
        """
        return self._cartan_type

    def classical(self):
        """
        Return the classical Lie algebra of ``self``.
        """
        return self._g

    def _construct_UEA(self):
        """
        Construct the universal enveloping algebra of ``self``.
        """
        # These are placeholders and will require something more complex
        if self._kac_moody:
            PolynomialRing(self._g.universal_enveloping_algebra(), 't,c,delta')
        return PolynomialRing(self._g.universal_enveloping_algebra(), 't,c')

    @cached_method
    def zero(self):
        """
        Return the element `0`.

        EXAMPLES::

            sage: L.<x,y> = LieAlgebra(QQ, representation="polynomial")
            sage: L.zero()
            0
        """
        zero = self.base_ring().zero()
        return self.element_class(self, {}, zero, zero)

    @cached_method
    def lie_algebra_generators(self):
        """
        Return the Lie algebra generators of ``self``.
        """
        zero = self.base_ring().zero()
        one = self.base_ring().one()
        d = {}
        if self._kac_moody:
            d['delta'] = self.element_class(self, {}, zero, one)
        d['c'] = self.element_class(self, {}, one, zero)
        try:
            finite_gens = dict(self._g.lie_algebra_generators(True))
        except (TypeError):
            finite_gens = dict(self._g.lie_algebra_generators())
        for k,g in finite_gens.items():
            d[k] = self.element_class(self, {0: g}, zero, zero)
        # e_0 = f_{\theta} t
        d['e0'] = self.element_class(self, {1: self._g.highest_root_basis_elt(False)},
                                     zero, zero)
        # f_0 = e_{\theta} t^-1
        d['f0'] = self.element_class(self, {-1: self._g.highest_root_basis_elt(True)},
                                     zero, zero)
        return Family(self.variable_names(), d.__getitem__)

    def monomial(self, m):
        """
        Construct the monomial indexed by ``m``.
        """
        if m == 'c' or m == 'delta':
            return self.lie_algebra_generators()[m]
        G = self._g.lie_algebra_generators()
        zero = self.base_ring().zero()
        return self.element_class(self, {m[0]: G[m[1]]}, zero, zero)

    Element = UntwistedAffineLieAlgebraElement

