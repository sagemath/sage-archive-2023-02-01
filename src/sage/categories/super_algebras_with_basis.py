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

class SuperAlgebrasWithBasis(SuperModulesCategory):
    """
    The category of super algebras with a distinguished basis

    EXAMPLES::

        sage: C = Algebras(ZZ).WithBasis().Super(); C
        Category of super algebras with basis over Integer Ring
        sage: sorted(C.super_categories(), key=str)
        [Category of graded algebras with basis over Integer Ring,
         Category of super algebras over Integer Ring,
         Category of super modules with basis over Integer Ring]

    TESTS::

        sage: TestSuite(C).run()
    """
    class ParentMethods:
        pass

    class ElementMethods:
        pass

