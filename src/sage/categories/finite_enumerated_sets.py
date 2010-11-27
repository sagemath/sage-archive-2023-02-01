r"""
Finite Enumerated Sets
"""
#*****************************************************************************
#  Copyright (C) 2009 Florent Hivert <Florent.Hivert@univ-rouen.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************


from category_types import Category
from sage.categories.enumerated_sets import EnumeratedSets
from sage.categories.isomorphic_objects   import IsomorphicObjectsCategory
from sage.rings.integer import Integer
from sage.misc.cachefunc import cached_method

class FiniteEnumeratedSets(Category):
    """
    The category of finite enumerated sets

    EXAMPLES::

        sage: FiniteEnumeratedSets()
        Category of finite enumerated sets
        sage: FiniteEnumeratedSets().super_categories()
        [Category of enumerated sets]
        sage: FiniteEnumeratedSets().all_super_categories()
        [Category of finite enumerated sets,
         Category of enumerated sets,
         Category of sets,
         Category of sets with partial maps,
         Category of objects]

    TESTS::

        sage: C = FiniteEnumeratedSets()
        sage: TestSuite(C).run()
    """

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: FiniteEnumeratedSets().super_categories()
            [Category of enumerated sets]
        """
        return [EnumeratedSets()]


    def _call_(self, X):
        """
        Construct an object in this category from the data in ``X``.

        EXAMPLES::

            sage: FiniteEnumeratedSets()(GF(3))
            Finite Field of size 3
            sage: FiniteEnumeratedSets()(Partitions(3)) # todo: not implemented: Partitions(3) is not yet in FiniteEnumeratedSets()
            Partitions of 3

        For now, lists, tuples, sets, Sets are coerced into finite
        enumerated sets::

            sage: FiniteEnumeratedSets()([1, 2, 3])
            {1, 2, 3}
            sage: FiniteEnumeratedSets()((1, 2, 3))
            {1, 2, 3}
            sage: FiniteEnumeratedSets()(set([1, 2, 3]))
            {1, 2, 3}
            sage: FiniteEnumeratedSets()(Set([1, 2, 3]))
            {1, 2, 3}
        """
        return EnumeratedSets()._call_(X)

    class ParentMethods:

        def is_finite(self):
            """
            Returns ``True`` since self is finite.

            EXAMPLES::

                sage: C = FiniteEnumeratedSets().example()
                sage: C.is_finite()
                True
            """
            return True

        def _cardinality_from_iterator(self):
            """
            The cardinality of ``self``.

            ``self.cardinality()`` returns the cardinality of the set ``self``
            as a sage ``Integer`` or as ``+Infinity``.

            This is the default (brute force) implementation from the
            category ``FiniteEnumeratedSet()``. It iterates through
            the elements of ``self`` to count them.

            EXAMPLES::

                sage: C = FiniteEnumeratedSets().example()
                sage: C.cardinality() # indirect doctest
                3
            """
            c = Integer(0)
            one = Integer(1)
            for _ in self:
                c += one
            return c
        #Set cardinality to the default implementation
        cardinality = _cardinality_from_iterator

        def _list_from_iterator(self):
            """
            The list of elements ``self``.

            ``self.list()`` returns the list of the element of the set
            ``self``. This is the default implementation from the
            category ``EnumeratedSet()`` which builds the list from
            the iterator.

            EXAMPLES::

                sage: C = FiniteEnumeratedSets().example()
                sage: C._list_from_iterator()
                [1, 2, 3]
                sage: C.list() # indirect doctest
                [1, 2, 3]
            """
            return [x for x in self]
        #Set list to the default implementation
        _list_default = _list_from_iterator # needed by the check mechanism.
        list  = _list_default


        def _random_element_from_unrank(self):
            """
            A random element in ``self``.

            ``self.random_element()`` returns a random element in
            ``self`` with uniform probability.

            This is the default implementation from the category
            ``EnumeratedSet()`` which uses the method ``unrank``.

            EXAMPLES::

                sage: C = FiniteEnumeratedSets().example()
                sage: C.random_element()
                1
                sage: C._random_element_from_unrank()
                2

            TODO: implement _test_random which checks uniformness
            """
            from sage.misc.prandom import randint
            c = self.cardinality()
            r = randint(0, c-1)
            return self.unrank(r)
        #Set the default implementation of random
        random_element = _random_element_from_unrank

        @cached_method
        def _last_from_iterator(self):
            """
            The last element of ``self``.

            ``self.last()`` returns the last element of ``self``.

            This is the default (brute force) implementation from the
            category ``FiniteEnumeratedSet()`` which can be used when
            the method ``__iter__`` is provided. Its complexity is
            `O(n)` where `n` is the size of ``self``.

            EXAMPLES::

                sage: C = FiniteEnumeratedSets().example()
                sage: C.last()
                3
                sage: C._last_from_iterator()
                3
            """
            for i in self:
                pass
            return i
        last = _last_from_iterator

        def _last_from_unrank(self):
            """
            The last element of ``self``.

            ``self.last()`` returns the last element of ``self``

            This is a generic implementation from the category
            ``FiniteEnumeratedSet()`` which can be used when the
            method ``unrank`` is provided.

            EXAMPLES::

                sage: C = FiniteEnumeratedSets().example()
                sage: C._last_from_unrank()
                3
            """
            return self.unrank(self.cardinality() -1)

        def _test_enumerated_set_iter_cardinality(self, **options):
            """
            Checks that the methods :meth:`.cardinality` and
            :meth:`.__iter__` are consistent.

            See also: :class:`TestSuite`.

            EXAMPLES::

                sage: C = FiniteEnumeratedSets().example()
                sage: C._test_enumerated_set_iter_cardinality()

            Let us now break the class::

                sage: from sage.categories.examples.finite_enumerated_sets import Example
                sage: class CCls(Example):
                ...       def cardinality(self):
                ...           return 4
                sage: CC = CCls()
                sage: CC._test_enumerated_set_iter_cardinality()
                Traceback (most recent call last):
                ...
                AssertionError: 4 != 3
            """
            tester = self._tester(**options)
            if self.cardinality != self._cardinality_from_iterator and self.cardinality() <= self.max_test_enumerated_set_loop:
                tester.assertEqual(self.cardinality(),
                                   self._cardinality_from_iterator())


    class IsomorphicObjects(IsomorphicObjectsCategory):

        def example(self):
            """
            Returns an example of isomorphic object of a finite
            enumerated set, as per :meth:`Category.example
            <sage.categories.category.Category.example>`.

            EXAMPLES::

                sage: FiniteEnumeratedSets().IsomorphicObjects().example()
                The image by some isomorphism of An example of a finite enumerated set: {1,2,3}
            """
            from sage.categories.examples.finite_enumerated_sets import IsomorphicObjectOfFiniteEnumeratedSet
            return IsomorphicObjectOfFiniteEnumeratedSet()

        class ParentMethods:

            def cardinality(self):
                r"""
                Returns the cardinality of ``self`` which is the same
                as that of the ambient set ``self`` is isomorphic to.

                EXAMPLES::

                    sage: A = FiniteEnumeratedSets().IsomorphicObjects().example(); A
                    The image by some isomorphism of An example of a finite enumerated set: {1,2,3}
                    sage: A.cardinality()
                    3
                """
                return self.ambient().cardinality()

            def __iter__(self):
                r"""
                Returns an iterator over ``self``, using the bijection
                with the ambient space.

                EXAMPLES::

                    sage: A = FiniteEnumeratedSets().IsomorphicObjects().example(); A
                    The image by some isomorphism of An example of a finite enumerated set: {1,2,3}
                    sage: list(A)                  # indirect doctest
                    [1, 4, 9]
                """
                for x in self.ambient():
                    yield self.retract(x)
