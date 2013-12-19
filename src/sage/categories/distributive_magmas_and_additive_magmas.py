r"""
Semirngs
"""
#*****************************************************************************
#  Copyright (C) 2010 Nicolas Borie <nicolas.borie@math.u-psud.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.lazy_import import LazyImport
from sage.categories.category_with_axiom import CategoryWithAxiom
from sage.categories.category_singleton import Category_singleton
from sage.categories.additive_magmas import AdditiveMagmas
from sage.categories.magmas import Magmas

class DistributiveMagmasAndAdditiveMagmas(Category_singleton):
    """
    The category of sets `(S,+,*)` with `*` distributing on `+`

    This is similar to a ring, but `+` and `*` are only required to be
    (additive) magmas.

    EXAMPLES::

        sage: DistributiveMagmasAndAdditiveMagmas()
        Category of distributive magmas and additive magmas
        sage: DistributiveMagmasAndAdditiveMagmas().super_categories()
        [Category of magmas, Category of additive magmas]

        sage: DistributiveMagmasAndAdditiveMagmas().Associative().AdditiveAssociative().AdditiveCommutative().AdditiveUnital().AdditiveInverse()
        Category of rngs
        sage: DistributiveMagmasAndAdditiveMagmas().Associative().AdditiveAssociative().AdditiveCommutative().AdditiveUnital().Unital()
        Category of semirings
        sage: DistributiveMagmasAndAdditiveMagmas().Associative().AdditiveAssociative().AdditiveCommutative().AdditiveUnital().AdditiveInverse().Unital()
        Category of rings

    TESTS::

        sage: TestSuite(DistributiveMagmasAndAdditiveMagmas()).run()
    """

    def super_categories(self):
        """
        EXAMPLES::

            sage: DistributiveMagmasAndAdditiveMagmas().super_categories()
            [Category of magmas, Category of additive magmas]
        """
        return [Magmas(), AdditiveMagmas()]

    class AdditiveAssociative(CategoryWithAxiom):
        class AdditiveCommutative(CategoryWithAxiom):
            class AdditiveUnital(CategoryWithAxiom):
                class AdditiveInverse(CategoryWithAxiom):
                    Associative = LazyImport('sage.categories.rngs', 'Rngs', at_startup=True)
                class Associative(CategoryWithAxiom):
                    Unital          = LazyImport('sage.categories.semirings', 'Semirings', at_startup=True)

    class ParentMethods:

        def _test_distributivity(self, **options):
            r"""
            Test the distributivity of `*` on `+` on (not necessarily
            all) elements of this set.

            INPUT::

             - ``options`` -- any keyword arguments accepted by :meth:`_tester`.

            EXAMPLES:

            By default, this method runs the tests only on the
            elements returned by ``self.some_elements()``::

                sage: NN.some_elements()
                [0, 1, 3, 42]
                sage: NN._test_distributivity()

            However, the elements tested can be customized with the
            ``elements`` keyword argument::

                sage: CC._test_distributivity(elements=[CC(0),CC(1),CC(3),CC(I)])

            See the documentation for :class:`TestSuite` for more information.
            """
            tester = self._tester(**options)
            S = tester.some_elements()
            from sage.combinat.cartesian_product import CartesianProduct
            for x,y,z in tester.some_elements(CartesianProduct(S,S,S)):
                # left distributivity
                tester.assert_(x * (y + z) == (x * y) + (x * z))
                # right distributivity
                tester.assert_((x + y) * z == (x * z) + (y * z))
