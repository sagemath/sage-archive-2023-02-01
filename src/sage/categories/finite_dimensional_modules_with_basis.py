r"""
FiniteDimensionalModulesWithBasis
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from category_types import Category_over_base_ring
from sage.categories.all import ModulesWithBasis
from sage.misc.cachefunc import cached_method

class FiniteDimensionalModulesWithBasis(Category_over_base_ring):
    """
    The category of finite dimensional modules with a distinguished basis

    EXAMPLES::

      sage: FiniteDimensionalModulesWithBasis(QQ)
      Category of finite dimensional modules with basis over Rational Field
      sage: FiniteDimensionalModulesWithBasis(QQ).super_categories()
      [Category of modules with basis over Rational Field]

    TESTS::

        sage: TestSuite(FiniteDimensionalModulesWithBasis(ZZ)).run()
    """

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: FiniteDimensionalModulesWithBasis(QQ).super_categories()
            [Category of modules with basis over Rational Field]
        """
        R = self.base_ring()
        return [ModulesWithBasis(R)]

    class ParentMethods:
        pass

    class ElementMethods:
        pass
