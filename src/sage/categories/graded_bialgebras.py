r"""
GradedBialgebras
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                     Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from category_types import Category_over_base_ring
from sage.categories.all import Bialgebras, GradedAlgebras, GradedCoalgebras
from sage.misc.cachefunc import cached_method

class GradedBialgebras(Category_over_base_ring):
    """
    The category of bialgebras with several bases

    EXAMPLES::

        sage: GradedBialgebras(ZZ)
        Category of graded bialgebras over Integer Ring
        sage: GradedBialgebras(ZZ).super_categories()
        [Category of graded algebras over Integer Ring, Category of graded coalgebras over Integer Ring, Category of bialgebras over Integer Ring]

    TESTS::

        sage: TestSuite(GradedBialgebras(ZZ)).run()
    """

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: GradedBialgebras(QQ).super_categories()
            [Category of graded algebras over Rational Field, Category of graded coalgebras over Rational Field, Category of bialgebras over Rational Field]
        """
        R = self.base_ring()
        return [GradedAlgebras(R), GradedCoalgebras(R), Bialgebras(R)]

    class ParentMethods:
        pass

    class ElementMethods:
        pass
