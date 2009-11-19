r"""
GradedModules
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                     Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from category_types import Category_over_base_ring
from sage.categories.all import Modules
from sage.misc.cachefunc import cached_method

class GradedModules(Category_over_base_ring):
    """
    The category of graded modules

    EXAMPLES::

        sage: GradedModules(ZZ)
        Category of graded modules over Integer Ring
        sage: GradedModules(ZZ).super_categories()
        [Category of modules over Integer Ring]

    TESTS::

        sage: TestSuite(GradedModules(ZZ)).run()
    """

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: GradedModules(QQ).super_categories()
            [Category of modules over Rational Field]
        """
        R = self.base_ring()
        return [Modules(R)]

    class ParentMethods:
        pass

    class ElementMethods:
        pass
