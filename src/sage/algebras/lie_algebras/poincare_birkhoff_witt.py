"""
The Poincare-Birkhoff-Witt Basis For A Universal Enveloping Algebra

AUTHORS:

- Travis Scrimshaw (2013-11-03): Initial version
"""

#*****************************************************************************
#  Copyright (C) 2013 Travis Scrimshaw <tscrim@ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty
#    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#  See the GNU General Public License for more details; the full text
#  is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.misc_c import prod
from sage.categories.algebras import Algebras
from sage.monoids.indexed_free_monoid import IndexedFreeAbelianMonoid
from sage.combinat.free_module import CombinatorialFreeModule
from sage.sets.family import Family
from sage.rings.all import ZZ

class PoincareBirkhoffWittBasis(CombinatorialFreeModule):
    r"""
    The Poincare-Birkhoff-Witt (PBW) basis of the universal enveloping
    algebra of a Lie algebra.

    Consider a Lie algebra `\mathfrak{g}` with basis `(b_i)_{i=1}^n`,
    so the universal enveloping algebra `U(\mathfrak{g})` is generated
    by `b_i` and subject to the relations

    .. MATH::

        [b_i, b_j] = \sum_{k \in I} c_{ij}^k b_k

    where `c_{ij}^k` are the structure coefficients of `\mathfrak{g}`. The
    Poincare-Birkhoff-Witt (PBW) basis is given by picking a particular
    ordering of `(b_i)_{i \in I}` and is given by the monomials
    `b_1^{e_1} b_2^{e_2} \cdots b_n^{e_n}`. Specifically, we can rewrite
    `b_j b_i = b_i b_j + [b_j, b_i]` where `j > i`, and we can repeat
    this to sort any monmoial into

    .. MATH::

        b_{i_1} \cdots b_{i_k} = b_1^{e_1} \cdots b_n^{e_n} + LOT

    where `LOT` are lower order terms. Thus the PBW is a filtered basis
    for `U(\mathfrak{g})`.

    EXAMPLES:

    We construct the PBW basis of `\mathfrak{sl}_2`::

        sage: L = lie_algebras.sl(QQ, 2)
        sage: PBW = L.pbw_basis()

    We then do some computations, in particular, we check that `[E, F] = H`::

        sage: E,F,H = PBW.algebra_generators()
        sage: E*F
        PBW['E']*PBW['F']
        sage: F*E
        PBW['E']*PBW['F'] - PBW['H']
        sage: E*F - F*E
        PBW['H']

    Next we construct another instance of the PBW basis, but sorted in the
    reverse order::

        sage: PBW2 = L.pbw_basis(prefix='PBW2', basis_cmp=lambda x,y: -cmp(x,y))

    We then check the multiplication is preserved::

        sage: PBW2(E) * PBW2(F)
        PBW2['F']*PBW2['E'] + PBW2['H']
        sage: PBW2(E*F)
        PBW2['F']*PBW2['E'] + PBW2['H']
        sage: F * E + H
        PBW['E']*PBW['F']

    We now construct the PBW basis for Lie algebra of regular
    vector fields on `\CC^{\times}`::

        sage: L = lie_algebras.regular_vector_fields(QQ)
        sage: PBW = L.pbw_basis()
        sage: G = PBW.algebra_generators()
        sage: G[2] * G[3]
        PBW[2]*PBW[3]
        sage: G[3] * G[2]
        PBW[2]*PBW[3] - PBW[5]
        sage: G[-2] * G[3] * G[2]
        PBW[-2]*PBW[2]*PBW[3] - PBW[-2]*PBW[5]
    """
    @staticmethod
    def __classcall_private__(cls, g, basis_cmp=None, prefix='PBW', **kwds):
        """
        Normalize input to ensure a unique representation.

        TESTS::

            sage: from sage.algebras.lie_algebras.poincare_birkhoff_witt import PoincareBirkhoffWittBasis
            sage: L = lie_algebras.sl(QQ, 2)
            sage: P1 = PoincareBirkhoffWittBasis(L)
            sage: P2 = PoincareBirkhoffWittBasis(L, prefix='PBW')
            sage: P1 is P2
            True
        """
        return super(PoincareBirkhoffWittBasis, cls).__classcall__(cls,
                            g, basis_cmp, prefix, **kwds)

    def __init__(self, g, basis_cmp, prefix, **kwds):
        """
        Initialize ``self``.

        TESTS::

            sage: L = lie_algebras.sl(QQ, 2)
            sage: PBW = L.pbw_basis()
            sage: E,F,H = PBW.algebra_generators()
            sage: TestSuite(PBW).run(elements=[E, F, H])
            sage: TestSuite(PBW).run(elements=[E, F, H, E*F + H]) # long time
        """
        if basis_cmp is not None:
            self._basis_cmp = basis_cmp

        R = g.base_ring()
        self._g = g
        monomials = IndexedFreeAbelianMonoid(g.basis().keys(), prefix,
                                             generator_cmp=self._monoid_cmp, **kwds)
        CombinatorialFreeModule.__init__(self, R, monomials,
                                         prefix='', bracket=False, latex_bracket=False,
                                         generator_cmp=self._monomial_cmp,
                                         category=Algebras(R).WithBasis())

    def _basis_cmp(self, x, y):
        """
        Compare the indices of ``x`` and ``y``.

        TESTS::

            sage: L = lie_algebras.sl(QQ, 2)
            sage: PBW = L.pbw_basis()
            sage: PBW._basis_cmp('E', 'H')
            -1

        ::

            sage: L = lie_algebras.sl(QQ, 2)
            sage: PBW = L.pbw_basis(basis_cmp=lambda x,y: -cmp(x,y))
            sage: prod(PBW.gens())
            PBW['H']*PBW['F']*PBW['E'] + PBW['H']^2
        """
        K = self._g.basis().keys()
        if K.cardinality() == float('inf'):
            return cmp(x, y)
        lst = list(K)
        return cmp(lst.index(x), lst.index(y))

    def _monoid_cmp(self, x, y):
        """
        Comparison function for the underlying monoid.

        EXAMPLES::

            sage: L = lie_algebras.sl(QQ, 2)
            sage: PBW = L.pbw_basis(basis_cmp=lambda x,y: -cmp(x,y))
            sage: M = PBW.basis().keys()
            sage: prod(M.gens()) # indirect doctest
            PBW['H']*PBW['F']*PBW['E']
        """
        return self._basis_cmp(x[0], y[0])

    def _monomial_cmp(self, x, y):
        """
        Compare the monomials ``x`` and ``y`` of ``self`` by reverse
        degree lexicographic order.

        EXAMPLES::

            sage: L = lie_algebras.sl(QQ, 2)
            sage: PBW = L.pbw_basis()
            sage: E,F,H = sorted(PBW.algebra_generators(), key=str)
            sage: H*F*F*E # indirect doctest
            PBW['E']*PBW['F']^2*PBW['H'] - 2*PBW['E']*PBW['F']^2
             - 2*PBW['F']*PBW['H']^2 + 6*PBW['F']*PBW['H'] - 4*PBW['F']

            sage: PBW = L.pbw_basis(basis_cmp=lambda x,y: -cmp(x,y))
            sage: E,F,H = sorted(PBW.algebra_generators(), key=str)
            sage: E*F*F*H # indirect doctest
            PBW['H']*PBW['F']^2*PBW['E'] + 2*PBW['H']^2*PBW['F']
             + 2*PBW['F']^2*PBW['E'] + 6*PBW['H']*PBW['F'] + 4*PBW['F']
        """
        c = cmp(len(y), len(x))
        if c:
            return c
        w = y.to_word_list()
        for i,a in enumerate(x.to_word_list()):
            c = self._basis_cmp(a, w[i])
            if c:
                return c
        return 0

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: L = lie_algebras.sl(QQ, 2)
            sage: L.pbw_basis()
            Universal enveloping algebra of sl2 over Rational Field
             in the Poincare-Birkhoff-Witt basis
        """
        return "Universal enveloping algebra of {} in the Poincare-Birkhoff-Witt basis".format(self._g)

    def _coerce_map_from_(self, R):
        """
        Return ``True`` if there is a coercion map from ``R`` to ``self``.

        EXAMPLES:

        We lift from the Lie algebra::

            sage: L = lie_algebras.sl(QQ, 2)
            sage: PBW = L.pbw_basis()
            sage: PBW.has_coerce_map_from(L)
            True
            sage: [PBW(g) for g in L.gens()]
            [PBW['E'], PBW['F'], PBW['H']]

        We can go between PBW bases under different sorting orders::

            sage: PBW2 = L.pbw_basis(basis_cmp=lambda x,y: -cmp(x,y))
            sage: E,F,H = sorted(PBW.algebra_generators(), key=str)
            sage: PBW2(E*F*H)
            PBW['H']*PBW['F']*PBW['E'] + PBW['H']^2
        """
        if R == self._g:
            # Make this into the lift map
            I = self._indices
            basis_function = lambda x: self.monomial(I.gen(x))
            # TODO: this diagonal, but with a smaller indexing set...
            return self._g.module_morphism(basis_function, codomain=self,
                                           triangular='upper', unitriangular=True)

        if isinstance(R, PoincareBirkhoffWittBasis) and self._g == R._g:
            I = self._indices
            def basis_function(x):
                return self.prod(self.monomial(I.gen(g)**e) for g,e in x._sorted_items())
            # TODO: this diagonal, but with a smaller indexing set...
            return R.module_morphism(basis_function, codomain=self)

        return super(PoincareBirkhoffWittBasis, self)._coerce_map_from_(R)

    def lie_algebra(self):
        """
        Return the underlying Lie algebra of ``self``.

        EXAMPLES::

            sage: L = lie_algebras.sl(QQ, 2)
            sage: PBW = L.pbw_basis()
            sage: PBW.lie_algebra() is L
            True
        """
        return self._g

    def algebra_generators(self):
        """
        Return the algebra generators of ``self``.

        EXAMPLES::

            sage: L = lie_algebras.sl(QQ, 2)
            sage: PBW = L.pbw_basis()
            sage: PBW.algebra_generators()
            Finite family {'H': PBW['H'], 'E': PBW['E'], 'F': PBW['F']}
        """
        G = self._indices.gens()
        return Family(self._indices._indices, lambda x: self.monomial(G[x]),
                      name="generator map")

    gens = algebra_generators

    @cached_method
    def one_basis(self):
        """
        Return the basis element indexing `1`.

        EXAMPLES::

            sage: L = lie_algebras.sl(QQ, 2)
            sage: PBW = L.pbw_basis()
            sage: ob = PBW.one_basis(); ob
            1
            sage: ob.parent()
            Free abelian monoid indexed by {'E', 'F', 'H'}
        """
        return self._indices.one()

    def product_on_basis(self, lhs, rhs):
        """
        Return the product of the two basis elements ``lhs`` and ``rhs``.

        EXAMPLES::

            sage: L = lie_algebras.sl(QQ, 2)
            sage: PBW = L.pbw_basis()
            sage: I = PBW.indices()
            sage: PBW.product_on_basis(I.gen('E'), I.gen('F'))
            PBW['E']*PBW['F']
            sage: PBW.product_on_basis(I.gen('E'), I.gen('H'))
            PBW['E']*PBW['H']
            sage: PBW.product_on_basis(I.gen('H'), I.gen('E'))
            PBW['E']*PBW['H'] + 2*PBW['E']
            sage: PBW.product_on_basis(I.gen('F'), I.gen('E'))
            PBW['E']*PBW['F'] - PBW['H']
            sage: PBW.product_on_basis(I.gen('F'), I.gen('H'))
            PBW['F']*PBW['H']
            sage: PBW.product_on_basis(I.gen('H'), I.gen('F'))
            PBW['F']*PBW['H'] - 2*PBW['F']
            sage: PBW.product_on_basis(I.gen('H')**2, I.gen('F')**2)
            PBW['F']^2*PBW['H']^2 - 8*PBW['F']^2*PBW['H'] + 16*PBW['F']^2

            sage: E,F,H = sorted(PBW.algebra_generators(), key=str)
            sage: E*F - F*E
            PBW['H']
            sage: H * F * E
            PBW['E']*PBW['F']*PBW['H'] - PBW['H']^2
            sage: E * F * H * E 
            PBW['E']^2*PBW['F']*PBW['H'] + 2*PBW['E']^2*PBW['F']
             - PBW['E']*PBW['H']^2 - 2*PBW['E']*PBW['H']
        """
        # Some trivial base cases
        if lhs == self.one_basis():
            return self.monomial(rhs)
        if rhs == self.one_basis():
            return self.monomial(lhs)

        I = self._indices
        trail = lhs.trailing_support()
        lead = rhs.leading_support()
        if self._basis_cmp(trail, lead) <= 0:
            return self.monomial(lhs * rhs)

        # Create the commutator
        # We have xy - yx = [x, y] -> xy = yx + [x, y] and we have x > y
        terms = self._g.monomial(trail).bracket(self._g.monomial(lead))
        lead = I.gen(lead)
        trail = I.gen(trail)
        terms = self.sum_of_terms((I.gen(t), c) for t,c in terms)
        terms += self.monomial(lead * trail)
        return self.monomial(lhs // trail) * terms * self.monomial(rhs // lead)

