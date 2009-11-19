r"""
GradedModulesWithBasis
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from category_types import Category_over_base_ring
from sage.categories.all import GradedModules, ModulesWithBasis
from sage.misc.cachefunc import cached_method

class GradedModulesWithBasis(Category_over_base_ring):
    """
    The category of graded modules with a distinguished basis

    EXAMPLES::

        sage: GradedModulesWithBasis(ZZ)
        Category of graded modules with basis over Integer Ring
        sage: GradedModulesWithBasis(ZZ).super_categories()
        [Category of graded modules over Integer Ring, Category of modules with basis over Integer Ring]

    TESTS::

        sage: TestSuite(GradedModulesWithBasis(ZZ)).run()
    """

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: GradedModulesWithBasis(QQ).super_categories()
            [Category of graded modules over Rational Field, Category of modules with basis over Rational Field]
        """
        R = self.base_ring()
        return [GradedModules(R), ModulesWithBasis(R)]

    class ParentMethods:
        pass

    class ElementMethods:
        pass
