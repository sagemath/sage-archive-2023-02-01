r"""
FiniteDimensionalCoalgebrasWithBasis
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from category_types import Category_over_base_ring
from sage.categories.all import FiniteDimensionalModulesWithBasis, CoalgebrasWithBasis
from sage.misc.cachefunc import cached_method

class FiniteDimensionalCoalgebrasWithBasis(Category_over_base_ring):
    """
    The category of finite dimensional coalgebras with a distinguished basis

    EXAMPLES::

      sage: FiniteDimensionalCoalgebrasWithBasis(QQ)
      Category of finite dimensional coalgebras with basis over Rational Field
      sage: FiniteDimensionalCoalgebrasWithBasis(QQ).super_categories()
      [Category of finite dimensional modules with basis over Rational Field, Category of coalgebras with basis over Rational Field]

    TESTS::

        sage: TestSuite(FiniteDimensionalCoalgebrasWithBasis(ZZ)).run()
    """

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: FiniteDimensionalCoalgebrasWithBasis(QQ).super_categories()
            [Category of finite dimensional modules with basis over Rational Field, Category of coalgebras with basis over Rational Field]
        """
        R = self.base_ring()
        return [FiniteDimensionalModulesWithBasis(R), CoalgebrasWithBasis(R)]

    class ParentMethods:
        pass

    class ElementMethods:
        pass
