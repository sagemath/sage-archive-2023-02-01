r"""
Semirings
"""
#*****************************************************************************
#  Copyright (C) 2010 Nicolas Borie <nicolas.borie@math.u-psud.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category import Category
from sage.misc.cachefunc import cached_method

class Semirings(Category):
    """
    The category of semirings.

    A semiring `(S,+,*)` is similar to a ring, but without the
    requirement that each element must have an additive inverse. In
    other words, it is a combination of a commutative additive monoid
    `(S,+)` and a multiplicative monoid `(S,*)`, where `*` distributes
    over `+`.

    See: http://en.wikipedia.org/wiki/Semiring

    EXAMPLES::

        sage: Semirings()
        Category of semirings
        sage: Semirings().super_categories()
        [Category of commutative additive monoids, Category of monoids]

    TESTS::

        sage: TestSuite(Semirings()).run()
    """

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: Semirings().super_categories()
            [Category of commutative additive monoids, Category of monoids]
        """
        from sage.categories.commutative_additive_monoids import CommutativeAdditiveMonoids
        from sage.categories.monoids import Monoids
        return [CommutativeAdditiveMonoids(), Monoids()]

    class ParentMethods:

        def _test_distributivity(self, **options):
            r"""
            Test the distributivity of `*` on `+` on (not necessarily
            all) elements of this semiring.

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
            for x in tester.some_elements():
                for y in tester.some_elements():
                    for z in tester.some_elements():
                        # left distributivity
                        tester.assert_(x * (y + z) == (x * y) + (x * z))
                        # right distributivity
                        tester.assert_((x + y) * z == (x * z) + (y * z))
