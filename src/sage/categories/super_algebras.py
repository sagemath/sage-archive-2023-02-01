r"""
Super Algebras
"""
#*****************************************************************************
#  Copyright (C) 2015 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.super_modules import SuperModulesCategory

class SuperAlgebras(SuperModulesCategory):
    """
    The category of super algebras.

    EXAMPLES::

        sage: Algebras(ZZ).Super()
        Category of super algebras over Integer Ring
        sage: Algebras(ZZ).Super().super_categories()
        [Category of graded algebras over Integer Ring,
         Category of super modules over Integer Ring]

    TESTS::

        sage: TestSuite(Algebras(ZZ).Super()).run()
    """
    class ParentMethods:
        pass

    class ElementMethods:
        pass

