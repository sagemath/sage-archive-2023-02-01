r"""
GradedCoalgebrasWithBasis
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from category_types import Category_over_base_ring
from sage.categories.all import GradedCoalgebras, GradedModulesWithBasis, CoalgebrasWithBasis
from sage.misc.cachefunc import cached_method

class GradedCoalgebrasWithBasis(Category_over_base_ring):
    """
    The category of graded coalgebras with a distinguished basis

    EXAMPLES::

        sage: GradedCoalgebrasWithBasis(ZZ)
        Category of graded coalgebras with basis over Integer Ring
        sage: GradedCoalgebrasWithBasis(ZZ).super_categories()
        [Category of graded modules with basis over Integer Ring, Category of graded coalgebras over Integer Ring, Category of coalgebras with basis over Integer Ring]

    TESTS::

        sage: TestSuite(GradedCoalgebrasWithBasis(ZZ)).run()
    """

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: GradedCoalgebrasWithBasis(QQ).super_categories()
            [Category of graded modules with basis over Rational Field, Category of graded coalgebras over Rational Field, Category of coalgebras with basis over Rational Field]
        """
        R = self.base_ring()
        return [GradedModulesWithBasis(R), GradedCoalgebras(R), CoalgebrasWithBasis(R)]

    class ParentMethods:
        pass

    class ElementMethods:
        pass
