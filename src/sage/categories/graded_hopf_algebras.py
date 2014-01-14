r"""
Graded Hopf algebras
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                     Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from category_types import Category_over_base_ring
from sage.categories.all import HopfAlgebras, GradedBialgebras
from sage.misc.cachefunc import cached_method

class GradedHopfAlgebras(Category_over_base_ring):
    """
    The category of GradedHopf algebras with several bases

    EXAMPLES::

        sage: GradedHopfAlgebras(ZZ)
        Category of graded hopf algebras over Integer Ring
        sage: GradedHopfAlgebras(ZZ).super_categories()
        [Category of graded bialgebras over Integer Ring, Category of hopf algebras over Integer Ring]

    TESTS::

        sage: TestSuite(GradedHopfAlgebras(ZZ)).run()
    """

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: GradedHopfAlgebras(QQ).super_categories()
            [Category of graded bialgebras over Rational Field, Category of hopf algebras over Rational Field]
        """
        R = self.base_ring()
        return [GradedBialgebras(R), HopfAlgebras(R)]

    class ParentMethods:
        pass

    class ElementMethods:
        pass
