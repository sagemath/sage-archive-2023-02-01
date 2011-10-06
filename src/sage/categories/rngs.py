r"""
Rngs
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category import Category
from sage.categories.category_singleton import Category_singleton
from sage.misc.cachefunc import cached_method

class Rngs(Category_singleton):
    """
    The category of rngs
    associative rings, not necessarily commutative, and not necessarily with  1
    This is a combination of an abelian group (+) and a semigroup (*),
    with * distributing over +

    EXAMPLES::

      sage: Rngs()
      Category of rngs
      sage: Rngs().super_categories()
      [Category of commutative additive groups, Category of semigroups]

    TESTS::

        sage: TestSuite(Rngs()).run()
    """

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: Rngs().super_categories()
            [Category of commutative additive groups, Category of semigroups]
        """
        from commutative_additive_groups import CommutativeAdditiveGroups
        from semigroups import Semigroups
        return [CommutativeAdditiveGroups(), Semigroups()]

    class ParentMethods:
        pass

    class ElementMethods:
        pass
