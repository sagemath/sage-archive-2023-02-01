r"""
Graded Coalgebras
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                     Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from category_types import Category_over_base_ring
from sage.categories.all import Coalgebras, GradedModules
from sage.misc.cachefunc import cached_method

class GradedCoalgebras(Category_over_base_ring):
    """
    The category of graded coalgebras

    EXAMPLES::

        sage: GradedCoalgebras(ZZ)
        Category of graded coalgebras over Integer Ring
        sage: GradedCoalgebras(ZZ).super_categories()
        [Category of graded modules over Integer Ring, Category of coalgebras over Integer Ring]

    TESTS::

        sage: TestSuite(GradedCoalgebras(ZZ)).run()
    """

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: GradedCoalgebras(QQ).super_categories()
            [Category of graded modules over Rational Field, Category of coalgebras over Rational Field]
        """
        R = self.base_ring()
        return [GradedModules(R), Coalgebras(R)]

    class ParentMethods:
        pass

    class ElementMethods:
        pass
