r"""
FiniteDimensionalHopfAlgebrasWithBasis
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from category_types import Category_over_base_ring
from sage.categories.all import FiniteDimensionalBialgebrasWithBasis, HopfAlgebrasWithBasis
from sage.misc.cachefunc import cached_method

class FiniteDimensionalHopfAlgebrasWithBasis(Category_over_base_ring):
    """
    The category of finite dimensional Hopf algebras with a distinguished basis

    EXAMPLES::

        sage: FiniteDimensionalHopfAlgebrasWithBasis(QQ) # fixme: Hopf should be capitalized
        Category of finite dimensional hopf algebras with basis over Rational Field
        sage: FiniteDimensionalHopfAlgebrasWithBasis(QQ).super_categories()
        [Category of finite dimensional bialgebras with basis over Rational Field, Category of hopf algebras with basis over Rational Field]

    TESTS::

        sage: TestSuite(FiniteDimensionalHopfAlgebrasWithBasis(ZZ)).run()
    """

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: FiniteDimensionalHopfAlgebrasWithBasis(QQ).super_categories()
            [Category of finite dimensional bialgebras with basis over Rational Field, Category of hopf algebras with basis over Rational Field]
        """
        R = self.base_ring()
        return [FiniteDimensionalBialgebrasWithBasis(R), HopfAlgebrasWithBasis(R)]

    class ParentMethods:
        pass

    class ElementMethods:
        pass
