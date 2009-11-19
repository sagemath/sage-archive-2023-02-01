r"""
LeftModules
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.all import CommutativeAdditiveGroups
from category_types import Category_over_base_ring
from sage.misc.cachefunc import cached_method

#?class LeftModules(Category_over_base_rng):
class LeftModules(Category_over_base_ring):
    """
    The category of left modules
    left modules over an rng (ring not necessarily with unit), i.e.
    an abelian group with left multiplation by elements of the rng

    EXAMPLES::

        sage: LeftModules(ZZ)
        Category of left modules over Integer Ring
        sage: LeftModules(ZZ).super_categories()
        [Category of commutative additive groups]

    TESTS::

        sage: TestSuite(LeftModules(ZZ)).run()
    """

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: LeftModules(QQ).super_categories()
            [Category of commutative additive groups]
        """
        return [CommutativeAdditiveGroups()]

    class ParentMethods:
        pass

    class ElementMethods:
        ## r * x
        pass
