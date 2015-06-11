r"""
Graded Hopf algebras with basis
"""
#*****************************************************************************
#  Copyright (C) 2008      Teresa  Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2008-2011 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.categories.tensor import tensor
from sage.categories.graded_modules import GradedModulesCategory
from sage.categories.with_realizations import WithRealizationsCategory
from sage.misc.cachefunc import cached_method

class GradedHopfAlgebrasWithBasis(GradedModulesCategory):
    """
    The category of graded Hopf algebras with a distinguished basis.

    EXAMPLES::

        sage: C = GradedHopfAlgebrasWithBasis(ZZ); C
        Category of graded hopf algebras with basis over Integer Ring
        sage: C.super_categories()
        [Category of hopf algebras with basis over Integer Ring,
         Category of graded algebras with basis over Integer Ring]

        sage: C is HopfAlgebras(ZZ).WithBasis().Graded()
        True
        sage: C is HopfAlgebras(ZZ).Graded().WithBasis()
        False

    TESTS::

        sage: TestSuite(C).run()
    """

    class ParentMethods:
        pass

    class ElementMethods:
        pass

    class WithRealizations(WithRealizationsCategory):

        @cached_method
        def super_categories(self):
            """
            EXAMPLES::

                sage: GradedHopfAlgebrasWithBasis(QQ).WithRealizations().super_categories()
                [Join of Category of hopf algebras over Rational Field
                     and Category of graded algebras over Rational Field]

            TESTS::

                sage: TestSuite(GradedHopfAlgebrasWithBasis(QQ).WithRealizations()).run()
            """
            from sage.categories.graded_hopf_algebras import GradedHopfAlgebras
            R = self.base_category().base_ring()
            return [GradedHopfAlgebras(R)]


    class Connected(CategoryWithAxiom_over_base_ring):

            class ParentMethods:

                def antipode_on_basis(self, indice):
                    """
                    MATH::

                        S(x) := -\sum_{x^L \neq x} S(x^L) \times x^R

                    in general or `x` if `\mid x \mid = 0`.

                    TESTS::

                        sage: F = FQSym(QQ).F()
                        sage: F.antipode_on_basis(Permutation([1]))
                        -F[1]

                    """
                    if self.monomial(indice) == self.one():
                        return self.one()
                    else:
                        S = self.antipode_on_basis
                        x__S_Id = tensor([self, self]).module_morphism(
                            lambda (a, b): S(a) * self.monomial(b),
                            codomain=self)
                        return -x__S_Id(
                            self.monomial(indice).coproduct()
                            - tensor([self(indice), self.one()])
                        )

                def antipode(self, elem):
                    r"""
                    TESTS::

                        sage: F = FQSym(QQ).F()
                        sage: F.antipode(F[3,1,2])
                        -F[1, 3, 2] + F[2, 1, 3] + F[2, 3, 1] - 2*F[3, 1, 2]
                    """
                    import itertools
                    return self.linear_combination(itertools.imap(
                        lambda (mon, coeff): \
                            (self.antipode_on_basis(mon), coeff),
                        elem.monomial_coefficients().iteritems()
                    ))

            class ElementMethods:

                def antipode(self):
                    r"""
                    TESTS::

                        sage: F = FQSym(QQ).F()
                        sage: F[1].antipode()
                        -F[1]
                    """
                    return self.parent().antipode(self)
