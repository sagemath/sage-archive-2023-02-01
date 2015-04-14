r"""
Super algebras with basis
"""
#*****************************************************************************
#  Copyright (C) 2015 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.super_modules import SuperModulesCategory
from sage.categories.algebras import Algebras
from sage.categories.modules import Modules

class SuperAlgebrasWithBasis(SuperModulesCategory):
    """
    The category of super algebras with a distinguished basis

    EXAMPLES::

        sage: C = Algebras(ZZ).WithBasis().Super(); C
        Category of super algebras with basis over Integer Ring
        sage: sorted(C.super_categories(), key=str)
        [Category of graded algebras with basis over Integer Ring,
         Category of super modules with basis over Integer Ring]

    TESTS::

        sage: TestSuite(C).run()
    """
    def extra_super_categories(self):
        """
        EXAMPLES::

            sage: Algebras(ZZ).WithBasis().Super().super_categories()
            [Category of graded algebras with basis over Integer Ring,
             Category of super modules with basis over Integer Ring]
        """
        return [self.base_category().Graded()]

