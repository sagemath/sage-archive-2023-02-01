r"""
Graded Algebras
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                     Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from category_types import Category_over_base_ring
from sage.categories.all import Algebras, GradedModules
from sage.misc.cachefunc import cached_method

class GradedAlgebras(Category_over_base_ring):
    """
    The category of graded algebras

    EXAMPLES::

        sage: GradedAlgebras(ZZ)
        Category of graded algebras over Integer Ring
        sage: GradedAlgebras(ZZ).super_categories()
        [Category of graded modules over Integer Ring, Category of algebras over Integer Ring]

    TESTS::

        sage: TestSuite(GradedAlgebras(ZZ)).run()
    """

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: GradedAlgebras(QQ).super_categories()
            [Category of graded modules over Rational Field, Category of algebras over Rational Field]
        """
        R = self.base_ring()
        return [GradedModules(R), Algebras(R)]

    class ParentMethods:
        pass

    class ElementMethods:
        pass
