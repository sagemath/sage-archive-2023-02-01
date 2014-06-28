r"""
Coalgebras with basis
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#  Copyright (C) 2008-2011 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.abstract_method import abstract_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.categories.all import ModulesWithBasis, tensor, Hom

class CoalgebrasWithBasis(CategoryWithAxiom_over_base_ring):
    """
    The category of coalgebras with a distinguished basis.

    EXAMPLES::

        sage: CoalgebrasWithBasis(ZZ)
        Category of coalgebras with basis over Integer Ring
        sage: sorted(CoalgebrasWithBasis(ZZ).super_categories(), key=str)
        [Category of coalgebras over Integer Ring,
         Category of modules with basis over Integer Ring]

    TESTS::

        sage: TestSuite(CoalgebrasWithBasis(ZZ)).run()
    """

    class ParentMethods:

        @abstract_method(optional = True)
        def coproduct_on_basis(self, i):
            """
            The coproduct of the algebra on the basis (optional).

            INPUT:

            - ``i`` -- the indices of an element of the basis of ``self``

            Returns the coproduct of the corresponding basis elements
            If implemented, the coproduct of the algebra is defined
            from it by linearity.

            EXAMPLES::

                sage: A = HopfAlgebrasWithBasis(QQ).example(); A
                An example of Hopf algebra with basis: the group algebra of the Dihedral group of order 6 as a permutation group over Rational Field
                sage: (a, b) = A._group.gens()
                sage: A.coproduct_on_basis(a)
                B[(1,2,3)] # B[(1,2,3)]
            """

        @lazy_attribute
        def coproduct(self):
            """
            If :meth:`coproduct_on_basis` is available, construct the
            coproduct morphism from ``self`` to ``self`` `\otimes`
            ``self`` by extending it by linearity. Otherwise, use
            :meth:`~Coalgebras.Realizations.ParentMethods.coproduct_by_coercion`,
            if available.

            EXAMPLES::

                sage: A = HopfAlgebrasWithBasis(QQ).example(); A
                An example of Hopf algebra with basis: the group algebra of the Dihedral group of order 6 as a permutation group over Rational Field
                sage: [a,b] = A.algebra_generators()
                sage: a, A.coproduct(a)
                (B[(1,2,3)], B[(1,2,3)] # B[(1,2,3)])
                sage: b, A.coproduct(b)
                (B[(1,3)], B[(1,3)] # B[(1,3)])

            """
            if self.coproduct_on_basis is not NotImplemented:
                # TODO: if self is a coalgebra, then one would want
                # to create a morphism of algebras with basis instead
                # should there be a method self.coproduct_hom_category?
                return Hom(self, tensor([self, self]), ModulesWithBasis(self.base_ring()))(on_basis = self.coproduct_on_basis)
            elif hasattr(self, "coproduct_by_coercion"):
                return self.coproduct_by_coercion

        @abstract_method(optional = True)
        def counit_on_basis(self, i):
            """
            The counit of the algebra on the basis (optional).

            INPUT:

            - ``i`` -- the indices of an element of the basis of ``self``

            Returns the counit of the corresponding basis elements
            If implemented, the counit of the algebra is defined
            from it by linearity.

            EXAMPLES::

                sage: A = HopfAlgebrasWithBasis(QQ).example(); A
                An example of Hopf algebra with basis: the group algebra of the Dihedral group of order 6 as a permutation group over Rational Field
                sage: (a, b) = A._group.gens()
                sage: A.counit_on_basis(a)
                1
            """

        @lazy_attribute
        def counit(self):
            """
            If :meth:`counit_on_basis` is available, construct the
            counit morphism from ``self`` to ``self`` `\otimes`
            ``self`` by extending it by linearity

            EXAMPLES::

                sage: A = HopfAlgebrasWithBasis(QQ).example(); A
                An example of Hopf algebra with basis: the group algebra of the Dihedral group of order 6 as a permutation group over Rational Field
                sage: [a,b] = A.algebra_generators()
                sage: a, A.counit(a)
                (B[(1,2,3)], 1)
                sage: b, A.counit(b)
                (B[(1,3)], 1)

            """
            if self.counit_on_basis is not NotImplemented:
                return self.module_morphism(self.counit_on_basis,codomain=self.base_ring())

    class ElementMethods:
        pass

