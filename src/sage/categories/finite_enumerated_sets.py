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
from sage.misc.decorators import sage_wraps

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

    TODO:

    :class:`sage.combinat.debruijn_sequence.DeBruijnSequences` should
    not inherit from this class. If that is solved, then
    :class:`FiniteEnumeratedSets` shall be turned into a subclass of
    :class:`~sage.categories.category_singleton.Category_singleton`.
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

        def _cardinality_from_iterator(self, *ignored_args, **ignored_kwds):
            """
            The cardinality of ``self``.

            OUTPUT: an ``Integer``

            This brute force implementation of :meth:`cardinality`
            iterates through the elements of ``self`` to count them.

            EXAMPLES::

                sage: C = FiniteEnumeratedSets().example(); C
                An example of a finite enumerated set: {1,2,3}
                sage: C._cardinality_from_iterator()
                3

            This is the default implementation of :meth:`cardinality`
            from the category ``FiniteEnumeratedSet()``. To test this,
            we need a fresh example::

                sage: from sage.categories.examples.finite_enumerated_sets import Example
                sage: class FreshExample(Example): pass
                sage: C = FreshExample(); C.rename("FreshExample")
                sage: C.cardinality
                <bound method FreshExample_with_category._cardinality_from_iterator of FreshExample>

            TESTS:

            This method shall return an ``Integer``; we test this
            here, because :meth:`_test_enumerated_set_iter_cardinality`
            does not do it for us::

                sage: type(C._cardinality_from_iterator())
                <type 'sage.rings.integer.Integer'>

            We ignore additional inputs since during doctests classes which
            override ``cardinality()`` call up to the category rather than
            their own ``cardinality()`` method (see :trac:`13688`)::

                sage: C = FiniteEnumeratedSets().example()
                sage: C._cardinality_from_iterator(algorithm='testing')
                3

            Here is a more complete example::

                sage: class TestParent(Parent):
                ...     def __init__(self):
                ...         Parent.__init__(self, category=FiniteEnumeratedSets())
                ...     def __iter__(self):
                ...         yield 1
                ...         return
                ...     def cardinality(self, dummy_arg):
                ...         return 1 # we don't want to change the semantics of cardinality()
                sage: P = TestParent()
                sage: P.cardinality(-1)
                1
                sage: v = P.list(); v
                [1]
                sage: P.cardinality()
                1
                sage: P.cardinality('use alt algorithm') # Used to break here: see :trac:`13688`
                1
                sage: P.cardinality(dummy_arg='use alg algorithm') # Used to break here: see :trac:`13688`
                1
            """
            c = 0
            for _ in self:
                c += 1
            return Integer(c)
        #Set cardinality to the default implementation
        cardinality = _cardinality_from_iterator

        def _list_from_iterator(self):
            """
            The list of the elements of ``self``.

            This implementation computes this list from the iterator
            of ``self``. This is used by the default implementation of
            :meth:`list`.

            EXAMPLES::

                sage: C = FiniteEnumeratedSets().example()
                sage: C._list_from_iterator()
                [1, 2, 3]
                sage: C.list() # indirect doctest
                [1, 2, 3]
            """
            return [x for x in self]

        def _cardinality_from_list(self, *ignored_args, **ignored_kwds):
            """
            The cardinality of ``self``.

            This implementation of :meth:`cardinality` computes the
            cardinality from :meth:`list` (which is
            cached). Reciprocally, calling ``self.list()`` makes this
            method the default implementation of :meth:`cardinality`.

            EXAMPLES::

                sage: C = FiniteEnumeratedSets().example()
                sage: C._cardinality_from_list()
                3

            We ignore additional inputs since during doctests classes which
            override ``cardinality()`` call up to the category rather than
            their own ``cardinality()`` method (see :trac:`13688`)::

                sage: C = FiniteEnumeratedSets().example()
                sage: C._cardinality_from_list(algorithm='testing')
                3
            """
            # We access directly the cache self._list to bypass the
            # copy that self.list() currently does each time.
            try:
                lst = self._list
            except AttributeError:
                lst = self.list()
            return Integer(len(lst))

        def _unrank_from_list(self, r):
            """
            The ``r``-th element of ``self``

            INPUT:

              - ``r`` -- an integer between ``0`` and ``n-1``,
                where ``n`` is the cardinality of ``self``.

            OUTPUT: the ``r``-th element of ``self``

            This implementation of :meth:`unrank` uses the method
            :meth:`list` (which is cached). Reciprocally, calling
            ``self.list()`` makes this method the default
            implementation of :meth:`unrank`.

            EXAMPLES::

                sage: C = FiniteEnumeratedSets().example()
                sage: C._unrank_from_list(1)
                2
            """
            # We access directly the cache self._list to bypass the
            # copy that self.list() currently does each time.
            try:
                lst = self._list
            except AttributeError:
                lst = self.list()
            try:
                return lst[r]
            except IndexError:
                raise ValueError, "the value must be between %s and %s inclusive"%(0,len(lst)-1)


        def list(self):
            """
            The list of the elements of ``self``.

            This default implementation from the category
            ``FiniteEnumeratedSet()`` computes the list of the
            elements of ``self`` from the iterator of ``self`` and
            caches the result. It moreover overrides the following
            methods to use this cache:

            - ``self.cardinality()``
            - ``self.__iter__()``    (but see below)
            - ``self.unrank()``

            .. seealso:: :meth:`_list_from_iterator`, :meth:`_cardinality_from_list`,
                :meth:`_iterator_from_list`, and :meth:`_unrank_from_list`

            EXAMPLES::

                sage: C = FiniteEnumeratedSets().example()
                sage: C.list()
                [1, 2, 3]

            .. warning::

                The overriding of ``self.__iter__`` to use the cache
                is ignored upon calls such as ``for x in C:`` or
                ``list(C)`` (which essentially ruins its purpose).
                Indeed, Python looks up the ``__iter__`` method
                directly in the class of ``C``, bypassing ``C``'s
                dictionary (see the Python reference manual,
                `Special method lookup for new-style classes <http://docs.python.org/reference/datamodel.html#special-method-lookup-for-new-style-classes>`_)

                Let's take an example::

                    sage: class Example(Parent):
                    ...       def __init__(self):
                    ...           Parent.__init__(self, category = FiniteEnumeratedSets())
                    ...       def __iter__(self):
                    ...           print "hello!"
                    ...           for x in [1,2,3]: yield x
                    sage: C = Example()
                    sage: list(C)
                    hello!
                    hello!
                    [1, 2, 3]
                    sage: list(C)
                    hello!
                    [1, 2, 3]

                Note that ``hello!`` actually gets printed twice in
                the first call to ``list(C)``. That's because of the
                current (dubious) implementation of
                :meth:`Parent.__len__`. Let's call :meth:`list`::

                    sage: C.list()
                    [1, 2, 3]

                Now we would want the original iterator of ``C`` not
                to be called anymore, but that's not the case::

                    sage: list(C)
                    hello!
                    [1, 2, 3]


            TESTS:

            To test if the caching and overriding works, we need a
            fresh finite enumerated set example, because the caching
            mechanism has already been triggered::

                sage: from sage.categories.examples.finite_enumerated_sets import Example
                sage: class FreshExample(Example): pass
                sage: C = FreshExample(); C.rename("FreshExample")
                sage: C.list
                <bound method FreshExample_with_category.list of FreshExample>
                sage: C.unrank
                <bound method FreshExample_with_category._unrank_from_iterator of FreshExample>
                sage: C.cardinality
                <bound method FreshExample_with_category._cardinality_from_iterator of FreshExample>

                sage: l1 = C.list(); l1
                [1, 2, 3]
                sage: C.list
                <bound method FreshExample_with_category.list of FreshExample>
                sage: C.unrank
                <bound method FreshExample_with_category._unrank_from_list of FreshExample>
                sage: C.cardinality
                <bound method FreshExample_with_category._cardinality_from_list of FreshExample>
                sage: C.__iter__
                <bound method FreshExample_with_category._iterator_from_list of FreshExample>

            We finally check that nothing breaks before and after
            calling explicitly the method ``.list()``::

                sage: class FreshExample(Example): pass
                sage: import __main__; __main__.FreshExample = FreshExample # Fake FreshExample being defined in a python module
                sage: C = FreshExample()
                sage: TestSuite(C).run()
                sage: C.list()
                [1, 2, 3]
                sage: TestSuite(C).run()
            """
            try:
                return self._list[:]
            except AttributeError:
                self._list = self._list_from_iterator()
                self.cardinality = self._cardinality_from_list
                self.__iter__ = self._iterator_from_list
                self.unrank = self._unrank_from_list
            return self._list[:]

        _list_default  = list # needed by the check system.

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
            :meth:`.__iter__` are consistent. Also checks that
            :meth:`.cardinality` returns an ``Integer``.

            For efficiency reasons, those tests are not run if
            :meth:`.cardinality` is
            :meth:`._cardinality_from_iterator`, or if ``self`` is too
            big.

            .. seealso:: :class:`TestSuite`.

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
                sage: class CCls(Example):
                ...       def cardinality(self):
                ...           return int(3)
                sage: CC = CCls()
                sage: CC._test_enumerated_set_iter_cardinality()
                Traceback (most recent call last):
                ...
                AssertionError: False is not true
            """
            tester = self._tester(**options)
            if self.cardinality != self._cardinality_from_iterator:
                card = self.cardinality()
                tester.assert_(isinstance(card, Integer))
                if card <= tester._max_runs:
                    tester.assertEqual(card,
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
