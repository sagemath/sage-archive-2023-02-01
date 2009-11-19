r"""
FiniteDimensionalBialgebrasWithBasis
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from category_types import Category_over_base_ring
from sage.categories.all import FiniteDimensionalAlgebrasWithBasis, FiniteDimensionalCoalgebrasWithBasis, Bialgebras
from sage.misc.cachefunc import cached_method

class FiniteDimensionalBialgebrasWithBasis(Category_over_base_ring):
    """
    The category of finite dimensional bialgebras with a distinguished basis

    EXAMPLES::

        sage: FiniteDimensionalBialgebrasWithBasis(QQ)
        Category of finite dimensional bialgebras with basis over Rational Field
        sage: FiniteDimensionalBialgebrasWithBasis(QQ).super_categories()
        [Category of finite dimensional algebras with basis over Rational Field, Category of finite dimensional coalgebras with basis over Rational Field, Category of bialgebras over Rational Field]

    TESTS::

        sage: TestSuite(FiniteDimensionalBialgebrasWithBasis(ZZ)).run()
    """

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: FiniteDimensionalBialgebrasWithBasis(QQ).super_categories()
            [Category of finite dimensional algebras with basis over Rational Field, Category of finite dimensional coalgebras with basis over Rational Field, Category of bialgebras over Rational Field]
        """
        R = self.base_ring()
        return [FiniteDimensionalAlgebrasWithBasis(R), FiniteDimensionalCoalgebrasWithBasis(R), Bialgebras(R)]

    class ParentMethods:
        pass

    class ElementMethods:
        pass
