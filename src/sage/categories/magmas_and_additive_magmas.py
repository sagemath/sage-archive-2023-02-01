r"""
Magmas and Additive Magmas
"""
#*****************************************************************************
#  Copyright (C) 2010 Nicolas Borie <nicolas.borie@math.u-psud.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import LazyImport
from sage.categories.category_singleton import Category_singleton
from sage.categories.additive_magmas import AdditiveMagmas
from sage.categories.magmas import Magmas

class MagmasAndAdditiveMagmas(Category_singleton):
    """
    The category of sets `(S,+,*)` with an additive operation '+' and a multiplicative operation `*`

    EXAMPLES::

        sage: from sage.categories.magmas_and_additive_magmas import MagmasAndAdditiveMagmas
        sage: C = MagmasAndAdditiveMagmas(); C
        Category of magmas and additive magmas

    This is the base category for the categories of rings and their variants::

        sage: C.Distributive()
        Category of distributive magmas and additive magmas
        sage: C.Distributive().Associative().AdditiveAssociative().AdditiveCommutative().AdditiveUnital().AdditiveInverse()
        Category of rngs
        sage: C.Distributive().Associative().AdditiveAssociative().AdditiveCommutative().AdditiveUnital().Unital()
        Category of semirings
        sage: C.Distributive().Associative().AdditiveAssociative().AdditiveCommutative().AdditiveUnital().AdditiveInverse().Unital()
        Category of rings

    This category is really meant to represent the intersection of the
    categories of :class:`Magmas` and :class:`AdditiveMagmas`; however
    Sage's infrastructure does not allow yet to model this::

        sage: Magmas() & AdditiveMagmas()
        Join of Category of magmas and Category of additive magmas

        sage: Magmas() & AdditiveMagmas()        # todo: not implemented
        Category of magmas and additive magmas

    TESTS::

        sage: TestSuite(MagmasAndAdditiveMagmas()).run()
    """

    class SubcategoryMethods:

        @cached_method
        def Distributive(self):
            """
            Return the full subcategory of the objects of ``self`` where `*` is distributive on `+`.

            EXAMPLES::

                sage: from sage.categories.magmas_and_additive_magmas import MagmasAndAdditiveMagmas
                sage: C = MagmasAndAdditiveMagmas().Distributive()

            Given that `Sage` does not know that
            :class:`MagmasAndAdditiveMagmas` is the intersection of
            :class:`Magmas` and :class:`AdditiveMagmas`, this method
            is not available for

                sage: Magmas() & AdditiveMagmas()
                Join of Category of magmas and Category of additive magmas

            Still, the natural syntax works::

                sage: (Magmas() & AdditiveMagmas()).Distributive()
                Category of distributive magmas and additive magmas

            thanks to a workaround implemented in
            :meth:`Magmas.SubcategoryMethods.Distributive`::

                sage: (Magmas() & AdditiveMagmas()).Distributive.__module__
                'sage.categories.magmas'

            TESTS::

                sage: TestSuite(C).run()
                sage: Fields().Distributive.__module__
                'sage.categories.magmas_and_additive_magmas'
            """
            return self._with_axiom('Distributive')

    def super_categories(self):
        """
        EXAMPLES::

            sage: from sage.categories.magmas_and_additive_magmas import MagmasAndAdditiveMagmas
            sage: MagmasAndAdditiveMagmas().super_categories()
            [Category of magmas, Category of additive magmas]
        """
        return [Magmas(), AdditiveMagmas()]

    Distributive = LazyImport('sage.categories.distributive_magmas_and_additive_magmas', 'DistributiveMagmasAndAdditiveMagmas', at_startup=True)
