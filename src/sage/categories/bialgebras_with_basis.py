r"""
BialgebrasWithBasis
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#  Copyright (C) 2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from category_types import Category_over_base_ring
from sage.categories.all import AlgebrasWithBasis, CoalgebrasWithBasis, Bialgebras

class BialgebrasWithBasis(Category_over_base_ring):
    """
    The category of bialgebras with a distinguished basis

    EXAMPLES::

        sage: BialgebrasWithBasis(ZZ)
        Category of bialgebras with basis over Integer Ring
        sage: BialgebrasWithBasis(ZZ).super_categories()
        [Category of algebras with basis over Integer Ring, Category of coalgebras with basis over Integer Ring, Category of bialgebras over Integer Ring]

    TESTS::

        sage: TestSuite(BialgebrasWithBasis(ZZ)).run()
    """

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: BialgebrasWithBasis(QQ).super_categories()
            [Category of algebras with basis over Rational Field, Category of coalgebras with basis over Rational Field, Category of bialgebras over Rational Field]
        """
        R = self.base_ring()
        return [AlgebrasWithBasis(R), CoalgebrasWithBasis(R), Bialgebras(R)]

    class ParentMethods:
        pass

    class ElementMethods:
        pass
