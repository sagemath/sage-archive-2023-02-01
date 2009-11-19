r"""
GradedAlgebrasWithBasis
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from category_types import Category_over_base_ring
from sage.categories.all import GradedAlgebras, GradedModulesWithBasis, AlgebrasWithBasis
from sage.misc.cachefunc import cached_method

class GradedAlgebrasWithBasis(Category_over_base_ring):
    """
    The category of graded algebras with a distinguished basis

    EXAMPLES::

        sage: GradedAlgebrasWithBasis(ZZ)
        Category of graded algebras with basis over Integer Ring
        sage: GradedAlgebrasWithBasis(ZZ).super_categories()
        [Category of graded modules with basis over Integer Ring, Category of graded algebras over Integer Ring, Category of algebras with basis over Integer Ring]

    TESTS::

        sage: TestSuite(GradedAlgebrasWithBasis(ZZ)).run()
    """

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: GradedAlgebrasWithBasis(QQ).super_categories()
            [Category of graded modules with basis over Rational Field, Category of graded algebras over Rational Field, Category of algebras with basis over Rational Field]
        """
        R = self.base_ring()
        return [GradedModulesWithBasis(R),GradedAlgebras(R), AlgebrasWithBasis(R)]

    class ParentMethods:
        pass

    class ElementMethods:
        pass
