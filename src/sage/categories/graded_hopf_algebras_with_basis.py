r"""
GradedHopfAlgebrasWithBasis
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from category_types import Category_over_base_ring
from sage.categories.all import AbstractCategory, GradedHopfAlgebras, HopfAlgebrasWithBasis, GradedBialgebrasWithBasis
from sage.misc.cachefunc import cached_method

class GradedHopfAlgebrasWithBasis(Category_over_base_ring):
    """
    The category of graded Hopf algebras with a distinguished basis

    EXAMPLES::

        sage: GradedHopfAlgebrasWithBasis(ZZ)
        Category of graded hopf algebras with basis over Integer Ring
        sage: GradedHopfAlgebrasWithBasis(ZZ).super_categories()
        [Category of graded bialgebras with basis over Integer Ring, Category of graded hopf algebras over Integer Ring, Category of hopf algebras with basis over Integer Ring]

    TESTS::

        sage: TestSuite(GradedHopfAlgebrasWithBasis(ZZ)).run()
    """

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: GradedHopfAlgebrasWithBasis(QQ).super_categories()
            [Category of graded bialgebras with basis over Rational Field, Category of graded hopf algebras over Rational Field, Category of hopf algebras with basis over Rational Field]
        """
        R = self.base_ring()
        return [GradedBialgebrasWithBasis(R), GradedHopfAlgebras(R), HopfAlgebrasWithBasis(R)]

    class ParentMethods:
        pass

    class ElementMethods:
        pass

    class AbstractCategory(AbstractCategory):
        @cached_method
        def super_categories(self):
            """
            EXAMPLES::

                sage: GradedHopfAlgebrasWithBasis(QQ).abstract_category().super_categories()
                [Category of graded hopf algebras over Rational Field]

            TESTS::

                sage: TestSuite(GradedHopfAlgebrasWithBasis(QQ).abstract_category()).run()
            """
            from sage.categories.graded_hopf_algebras import GradedHopfAlgebras
            R = self.base_category.base_ring()
            return [GradedHopfAlgebras(R)]
