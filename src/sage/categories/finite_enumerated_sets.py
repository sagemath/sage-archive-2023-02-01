r"""
Finite Enumerated Sets
"""
#*****************************************************************************
#  Copyright (C) 2009 Florent Hivert <Florent.Hivert@univ-rouen.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************
from __future__ import print_function

import itertools

from sage.categories.category_with_axiom import CategoryWithAxiom
from sage.categories.enumerated_sets import EnumeratedSets
from sage.categories.sets_cat import Sets
from sage.categories.cartesian_product import CartesianProductsCategory
from sage.categories.isomorphic_objects   import IsomorphicObjectsCategory
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import lazy_import
lazy_import("sage.rings.integer", "Integer")

class FiniteEnumeratedSets(CategoryWithAxiom):
    """
    The category of finite enumerated sets

    EXAMPLES::

        sage: FiniteEnumeratedSets()
        Category of finite enumerated sets
        sage: FiniteEnumeratedSets().super_categories()
        [Category of enumerated sets, Category of finite sets]
        sage: FiniteEnumeratedSets().all_super_categories()
        [Category of finite enumerated sets,
         Category of enumerated sets,
         Category of finite sets,
         Category of sets,
         Category of sets with partial maps,
         Category of objects]

    TESTS::

        sage: C = FiniteEnumeratedSets()
        sage: TestSuite(C).run()
        sage: sorted(C.Algebras(QQ).super_categories(), key=str)
        [Category of finite dimensional modules with basis over Rational Field,
         Category of set algebras over Rational Field]

    .. TODO::

        :class:`sage.combinat.debruijn_sequence.DeBruijnSequences` should
        not inherit from this class. If that is solved, then
        :class:`FiniteEnumeratedSets` shall be turned into a subclass of
        :class:`~sage.categories.category_singleton.Category_singleton`.
    """

    def _call_(self, X):
        """
        Construct an object in this category from the data in ``X``.

        EXAMPLES::

            sage: FiniteEnumeratedSets()(GF(3))
            Finite Field of size 3
            sage: Partitions(3)
            Partitions of the integer 3

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

        @cached_method
        def cardinality(self, *ignored_args, **ignored_kwds):
            """
            Return the number of elements in self.

            EXAMPLES::
                sage: p=Partitions(3)
                sage: p.cardinality()
                3
            """
            return Integer(len(self.list()))

        def unrank_range(self, start=None, stop=None, step=None):
            r"""
            Return the range of elements of ``self`` starting at ``start``,
            ending at ``stop``, and stepping by ``step``.

            See also ``unrank()``.

            EXAMPLES::

                sage: F = FiniteEnumeratedSet([1,2,3])
                sage: F.unrank_range(1)
                [2, 3]
                sage: F.unrank_range(stop=2)
                [1, 2]
                sage: F.unrank_range(stop=2, step=2)
                [1]
                sage: F.unrank_range(start=1, step=2)
                [2]
                sage: F.unrank_range(stop=-1)
                [1, 2]

                sage: F = FiniteEnumeratedSet([1,2,3,4])
                sage: F.unrank_range(stop=10)
                [1, 2, 3, 4]
            """
            try:
                return self._list[start:stop:step]
            except AttributeError:
                pass
            card = self.cardinality() # This may set the list
            try:
                return self._list[start:stop:step]
            except AttributeError:
                pass
            if start is None and stop is not None and stop >= 0 and step is None:
                if stop < card:
                    it = self.__iter__()
                    return [next(it) for j in range(stop)]
                return self.list()
            return self.list()[start:stop:step]

        def iterator_range(self, start=None, stop=None, step=None):
            r"""
            Iterate over the range of elements of ``self`` starting
            at ``start``, ending at ``stop``, and stepping by ``step``.

            .. SEEALSO::

                ``unrank()``, ``unrank_range()``

            EXAMPLES::

                sage: F = FiniteEnumeratedSet([1,2,3])
                sage: list(F.iterator_range(1))
                [2, 3]
                sage: list(F.iterator_range(stop=2))
                [1, 2]
                sage: list(F.iterator_range(stop=2, step=2))
                [1]
                sage: list(F.iterator_range(start=1, step=2))
                [2]
                sage: list(F.iterator_range(start=1, stop=2))
                [2]
                sage: list(F.iterator_range(start=0, stop=1))
                [1]
                sage: list(F.iterator_range(start=0, stop=3, step=2))
                [1, 3]
                sage: list(F.iterator_range(stop=-1))
                [1, 2]

                sage: F = FiniteEnumeratedSet([1,2,3,4])
                sage: list(F.iterator_range(start=1, stop=3))
                [2, 3]
                sage: list(F.iterator_range(stop=10))
                [1, 2, 3, 4]
            """
            L = None
            try:
                L = self._list
            except AttributeError:
                pass
            card = self.cardinality() # This may set the list
            try:
                L = self._list
            except AttributeError:
                pass
            if L is None and start is None and stop is not None and stop >= 0 and step is None:
                if stop < card:
                    it = self.__iter__()
                    for j in range(stop):
                        yield next(it)
                    return
                for x in self:
                    yield x
                return
            if L is None:
                L = self.list()
            for x in L[start:stop:step]:
                yield x

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

    class CartesianProducts(CartesianProductsCategory):

        def extra_super_categories(self):
            """
            A Cartesian product of finite enumerated sets is a finite
            enumerated set.

            EXAMPLES::

                sage: C = FiniteEnumeratedSets().CartesianProducts()
                sage: C.extra_super_categories()
                [Category of finite enumerated sets]
            """
            return [FiniteEnumeratedSets()]

        class ParentMethods:
            r"""
            TESTS:

            Ideally, these tests should be just after the declaration of the
            associated attributes. But doing this way, Sage will not consider
            them as a doctest.

            We check that Cartesian products of finite enumerated sets
            inherit various methods from `Sets.CartesianProducts`
            and not from :class:`EnumeratedSets.Finite`::

                sage: C = cartesian_product([Partitions(10), Permutations(20)])
                sage: C in EnumeratedSets().Finite()
                True

                sage: C.random_element.__module__
                'sage.categories.sets_cat'

                sage: C.cardinality.__module__
                'sage.categories.sets_cat'

                sage: C.__iter__.__module__
                'sage.categories.sets_cat'
            """

            # Ambiguity resolution between methods inherited from
            # Sets.CartesianProducts and from EnumeratedSets.Finite.
            random_element = Sets.CartesianProducts.ParentMethods.random_element.__func__
            cardinality = Sets.CartesianProducts.ParentMethods.cardinality.__func__
            __iter__ = Sets.CartesianProducts.ParentMethods.__iter__.__func__

            def last(self):
                r"""
                Return the last element

                EXAMPLES::

                    sage: C = cartesian_product([Zmod(42), Partitions(10), IntegerRange(5)])
                    sage: C.last()
                    (41, [1, 1, 1, 1, 1, 1, 1, 1, 1, 1], 4)
                """
                return self._cartesian_product_of_elements(
                        tuple(c.last() for c in self.cartesian_factors()))

            def rank(self, x):
                r"""
                Return the rank of an element of this Cartesian product.

                The *rank* of ``x`` is its position in the
                enumeration. It is an integer between ``0`` and
                ``n-1`` where ``n`` is the cardinality of this set.

                .. SEEALSO::

                    - :meth:`EnumeratedSets.ParentMethods.rank`
                    - :meth:`unrank`

                EXAMPLES::

                    sage: C = cartesian_product([GF(2), GF(11), GF(7)])
                    sage: C.rank(C((1,2,5)))
                    96
                    sage: C.rank(C((0,0,0)))
                    0

                    sage: for c in C: print(C.rank(c))
                    0
                    1
                    2
                    3
                    4
                    5
                    ...
                    150
                    151
                    152
                    153

                    sage: F1 = FiniteEnumeratedSet('abcdefgh')
                    sage: F2 = IntegerRange(250)
                    sage: F3 = Partitions(20)
                    sage: C = cartesian_product([F1, F2, F3])
                    sage: c = C(('a', 86, [7,5,4,4]))
                    sage: C.rank(c)
                    54213
                    sage: C.unrank(54213)
                    ('a', 86, [7, 5, 4, 4])
                """
                from builtins import zip
                from sage.rings.integer_ring import ZZ
                x = self(x)
                b = ZZ.one()
                rank = ZZ.zero()
                for f, c in zip(reversed(x.cartesian_factors()),
                                reversed(self.cartesian_factors())):
                    rank += b * c.rank(f)
                    b *= c.cardinality()
                return rank

            def unrank(self, i):
                r"""
                Return the ``i``-th element of this Cartesian product.

                INPUT:

                - ``i`` -- integer between ``0`` and ``n-1`` where
                  ``n`` is the cardinality of this set.

                .. SEEALSO::

                    - :meth:`EnumeratedSets.ParentMethods.unrank`
                    - :meth:`rank`

                EXAMPLES::

                    sage: C = cartesian_product([GF(3), GF(11), GF(7), GF(5)])
                    sage: c = C.unrank(123); c
                    (0, 3, 3, 3)
                    sage: C.rank(c)
                    123

                    sage: c = C.unrank(857); c
                    (2, 2, 3, 2)
                    sage: C.rank(c)
                    857

                    sage: C.unrank(2500)
                    Traceback (most recent call last):
                    ...
                    IndexError: index i (=2) is greater than the cardinality
                """
                from sage.rings.integer_ring import ZZ
                i = ZZ(i)
                if i < 0:
                    raise IndexError("i (={}) must be a non-negative integer")
                elt = []
                for c in reversed(self.cartesian_factors()):
                    card = c.cardinality()
                    elt.insert(0, c.unrank(i % card))
                    i //= card
                if i:
                    raise IndexError("index i (={}) is greater than the cardinality".format(i))
                return self._cartesian_product_of_elements(elt)

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
