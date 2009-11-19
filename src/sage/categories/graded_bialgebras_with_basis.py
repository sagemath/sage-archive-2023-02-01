r"""
GradedBialgebrasWithBasis
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from category_types import Category_over_base_ring
from sage.categories.all import GradedBialgebras, GradedAlgebrasWithBasis, GradedCoalgebrasWithBasis
from sage.misc.cachefunc import cached_method

class GradedBialgebrasWithBasis(Category_over_base_ring):
    """
    The category of graded bialgebras with a distinguished basis

    EXAMPLES::

        sage: GradedBialgebrasWithBasis(ZZ)
        Category of graded bialgebras with basis over Integer Ring
        sage: GradedBialgebrasWithBasis(ZZ).super_categories()
        [Category of graded algebras with basis over Integer Ring, Category of graded coalgebras with basis over Integer Ring, Category of graded bialgebras over Integer Ring]

    TESTS::

        sage: TestSuite(GradedBialgebrasWithBasis(ZZ)).run()
    """

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: GradedBialgebrasWithBasis(QQ).super_categories()
            [Category of graded algebras with basis over Rational Field, Category of graded coalgebras with basis over Rational Field, Category of graded bialgebras over Rational Field]
        """
        R = self.base_ring()
        return [GradedAlgebrasWithBasis(R), GradedCoalgebrasWithBasis(R), GradedBialgebras(R)]

    class ParentMethods:
        pass

    class ElementMethods:
        pass
