r"""
CommutativeAdditiveMonoids
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category import Category
from sage.misc.cachefunc import cached_method
from sage.categories.commutative_additive_semigroups import CommutativeAdditiveSemigroups

# CHANGE: AbelianMonoid does not inherit any more from Monoids
class CommutativeAdditiveMonoids(Category):
    """
    The category of abelian monoids
    semigroups with an additive identity element

    EXAMPLES::

        sage: CommutativeAdditiveMonoids()
        Category of commutative additive monoids
        sage: CommutativeAdditiveMonoids().super_categories()
        [Category of commutative additive semigroups]
        sage: CommutativeAdditiveMonoids().all_super_categories()
        [Category of commutative additive monoids, Category of commutative additive semigroups, Category of additive magmas, Category of sets, Category of sets with partial maps, Category of objects]

    TESTS::

        sage: TestSuite(CommutativeAdditiveMonoids()).run()
    """

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: CommutativeAdditiveMonoids().super_categories()
            [Category of commutative additive semigroups]
        """
        return [CommutativeAdditiveSemigroups()]

    class ParentMethods:

        def _test_zero(self, **options):
            r"""
            Test that ``self.zero()`` is an element of self and is neutral for the addition

            INPUT::

             - ``options`` -- any keyword arguments accepted by :meth:`_tester`.

            EXAMPLES::

            By default, this method tests only the elements returned by
            ``self.some_elements()``::

                sage: S = CommutativeAdditiveMonoids().example()
                sage: S._test_zero()

            However, the elements tested can be customized with the
            ``elements`` keyword argument::

                sage: (a,b,c,d) = S.additive_semigroup_generators()
                sage: S._test_zero(elements = (a, a+c))

            See the documentation for :class:`TestSuite` for more information.
            """
            tester = self._tester(**options)
            zero = self.zero()
            # TODO: also call is_zero once it will work
            tester.assert_(zero.parent() == self)
            for x in tester.some_elements():
                tester.assert_(x + zero == x)

        @cached_method
        def zero(self):
            """
            Returns the zero of the abelian monoid, that is the unique neutral element for `+`.

            The default implementation is to coerce 0 into self.

            It is recommended to override this method because the
            coercion from the integers:

             - is not always meaningful (except for `0`),
             - often uses self.zero() otherwise.

            EXAMPLES::

                sage: S = CommutativeAdditiveMonoids().example()
                sage: S.zero()
                0
            """
            # TODO: add a test that actually exercise this default implementation
            return self(0)

        def zero_element(self):
            """
            Backward compatibility alias for self.zero()

            TESTS::

                sage: S = CommutativeAdditiveMonoids().example()
                sage: S.zero_element()
                0
            """
            return self.zero()

        def sum(self, args):
            r"""
            n-ary sum

            INPUT:

             - ``args`` -- a list (or iterable) of elements of ``self``

            Returns the sum of the elements in `args`, as an element of `self`.

            EXAMPLES::

                sage: S = CommutativeAdditiveMonoids().example()
                sage: (a,b,c,d) = S.additive_semigroup_generators()
                sage: S.sum((a,b,a,c,a,b))
                3*a + c + 2*b
                sage: S.sum(())
                0
                sage: S.sum(()).parent() == S
                True
            """
            return sum(args, self.zero())

    class ElementMethods:
        pass

#         TODO: merge with the implementation in Element which currently
#         overrides this one, and further requires self.parent()(0) to work.
#
#         def is_zero(self):
#             """
#             Returns whether self is the zero of the monoid

#             The default implementation, is to compare with ``self.zero()``.

#             TESTS::

#                 sage: S = CommutativeAdditiveMonoids().example()
#                 sage: S.zero().is_zero()
#                 True
#                 sage: S("aa").is_zero()
#                 False
#             """
#             return self == self.parent().zero()

