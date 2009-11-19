r"""
Enumerated Sets
"""
#*****************************************************************************
#  Copyright (C) 2009 Florent Hivert <Florent.Hivert@univ-rouen.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from category_types import Category
from sage.categories.sets_cat import Sets

class EnumeratedSets(Category):
    """
    The category of enumerated sets

    An *enumerated set* is a *finite* or *countable* set or multiset `S`
    together with a canonical enumeration of its elements;
    conceptually, this is very similar to an immutable list. The main
    difference lies in the names and the return type of the methods,
    and of course the fact that the list of element is not supposed to
    be expanded in memory. Whenever possible one should use one of the
    two sub-categories :class:`FiniteEnumeratedSets` or
    :class:`InfiniteEnumeratedSets`.

    The purpose of this category is threefold:

     - to fix a common interface for all these sets;
     - to provide a bunch of default implementations;
     - to provide consistency tests.

    The standard methods for an enumerated set ``S`` are:

       - ``S.cardinality()``: the number of elements of the set. This
         is the equivalent for ``len`` on a list except that the
         return value is specified to be a Sage :class:`Integer` or
         ``infinity``, instead of a Python ``int``;

       - ``iter(S)``: an iterator for the elements of the set;

       - ``S.list()``: the list of the elements of the set, when
         possible; raises a NotImplementedError if the list is
         predictably too large to be expanded in memory.

       - ``S.unrank(n)``: the  ``n-th`` element of the set when ``n`` is a sage
         ``Integer``. This is the equivanlent for ``l[n]`` on a list.

       - ``S.rank(e)``: the position of the element ``e`` in the set;
         This is equivalent to ``l.index(e)`` for a list except that
         the return value is specified to be a Sage :class:`Integer`,
         instead of a Python ``int``;

       - ``S.first()``: the first object of the set; it is equivalent to
         ``S.unrank(0)``;

       - ``S.next(e)``: the object of the set which follows ``e``; It is
         equivalent to ``S.unrank(S.rank(e)+1)``.

       - ``S.random_element()``: a random generator for an element of
         the set. Unless otherwise stated, and for finite enumerated
         sets, the probability is uniform.

    For examples and tests see:

       - ``FiniteEnumeratedSets().example()``
       - ``InfiniteEnumeratedSets().example()``


    EXAMPLES::

        sage: EnumeratedSets()
        Category of enumerated sets
        sage: EnumeratedSets().super_categories()
        [Category of sets]
        sage: EnumeratedSets().all_super_categories()
        [Category of enumerated sets, Category of sets, Category of objects]

    TESTS::

        sage: C = EnumeratedSets()
        sage: TestSuite(C).run()
    """

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: EnumeratedSets().super_categories()
            [Category of sets]
        """
        return [Sets()]


    class ParentMethods:

        def __iter__(self):
            """
            An iterator for the enumerated set.

            ``iter(self)`` allows the combinatorial class to be treated as an
            iterable. This if the default implementation from the category
            ``EnumeratedSet()`` it just goes through the iterator of the set
            to count the number of objects.

            By decreasing order of priority, the second column of the
            following array shows which methods is used to define
            ``__iter__``, when the methods of the first column are overloaded:

            +------------------------+---------------------------------+
            | Needed methods         | Default ``__iterator`` provided |
            +========================+=================================+
            | ``first`` and ``next`` | ``_iterator_from_next``         |
            +------------------------+---------------------------------+
            | ``unrank``             | ``_iterator_from_unrank``       |
            +------------------------+---------------------------------+
            | ``list`                | ``_iterator_from_next``         |
            +------------------------+---------------------------------+

            If non of these are provided raise a ``NotImplementedError``

            EXAMPLES::

            We start with an example where nothing is implemented::

                sage: class broken(UniqueRepresentation, Parent):
                ...    def __init__(self):
                ...        Parent.__init__(self, category = EnumeratedSets())
                ...
                sage: it = iter(broken()); [it.next(), it.next(), it.next()]
                Traceback (most recent call last):
                ...
                NotImplementedError: iterator called but not implemented

            Here is what happends when ``first`` and ``next`` are implemeted::

                sage: class set_first_next(UniqueRepresentation, Parent):
                ...    def __init__(self):
                ...        Parent.__init__(self, category = EnumeratedSets())
                ...    def first(self):
                ...        return 0
                ...    def next(self, elt):
                ...        return elt+1
                ...
                sage: it = iter(set_first_next()); [it.next(), it.next(), it.next()]
                [0, 1, 2]

            Let us try with ``unrank``::

                sage: class set_unrank(UniqueRepresentation, Parent):
                ...    def __init__(self):
                ...        Parent.__init__(self, category = EnumeratedSets())
                ...    def unrank(self, i):
                ...        return i + 5
                ...
                sage: it = iter(set_unrank()); [it.next(), it.next(), it.next()]
                [5, 6, 7]

            Let us finally try with ``list``::

                sage: class set_list(UniqueRepresentation, Parent):
                ...    def __init__(self):
                ...        Parent.__init__(self, category = EnumeratedSets())
                ...    def list(self):
                ...        return [5, 6, 7]
                ...
                sage: it = iter(set_list()); [it.next(), it.next(), it.next()]
                [5, 6, 7]

            """
            #Check to see if .first() and .next() are overridden in the subclass
            if ( self.first != self._first_from_iterator and
                 self.next  != self._next_from_iterator ):
                return self._iterator_from_next()
            #Check to see if .unrank() is overridden in the subclass
            elif self.unrank != self._unrank_from_iterator:
                return self._iterator_from_unrank()
            #Finally, check to see if .list() is overridden in the subclass
            elif self.list != self._list_default:
                return self._iterator_from_list()
            else:
                raise NotImplementedError, "iterator called but not implemented"

        def cardinality(self):
            """
            The cardinality of ``self``.

            ``self.cardinality()`` should return the cardinality of the set
            ``self`` as a sage :class:`Integer` or as ``infinity``.

            This if the default implementation from the category
            ``EnumeratedSet()`` it returns `NotImplementedError` since one does
            not know whether the set is finite or not.

            EXAMPLES::

                sage: class broken(UniqueRepresentation, Parent):
                ...    def __init__(self):
                ...        Parent.__init__(self, category = EnumeratedSets())
                ...
                sage: broken().cardinality()
                Traceback (most recent call last):
                ...
                NotImplementedError: unknown cardinality
            """
            raise NotImplementedError, "unknown cardinality"

        def list(self):
            """
            Returns an error since the cardinality of self is not known.

            EXAMPLES::

                sage: class broken(UniqueRepresentation, Parent):
                ...    def __init__(self):
                ...        Parent.__init__(self, category = EnumeratedSets())
                ...
                sage: broken().list()
                Traceback (most recent call last):
                ...
                NotImplementedError: unknown cardinality
            """
            raise NotImplementedError, "unknown cardinality"
        _list_default  = list # needed by the check system.


        def _first_from_iterator(self):
            """
            The "first" element of ``self``.

            ``self.first()`` returns the first element of the set
            ``self``. This is a generic implementation from the category
            ``EnumeratedSet()`` which can be used when the method ``__iter__`` is
            provided.

            EXAMPLES::

                sage: C = FiniteEnumeratedSets().example()
                sage: C.first() # indirect doctest
                1
            """
            it = self.__iter__()
            return it.next()
        first = _first_from_iterator

        def _next_from_iterator(self, obj):
            """
            The "next" element after ``obj`` in ``self``.

            ``self.next(e)`` returns the element of the set ``self`` which
            follows ``e``. This is a generic implementation from the category
            ``EnumeratedSet()`` which can be used when the method ``__iter__``
            is provided.

            Remark: this is the default (brute force) implementation
            of the category ``EnumeratedSets()``. Its complexity is
            `O(r)`, where `r` is the rank of ``obj``.

            EXAMPLES::

                sage: C = InfiniteEnumeratedSets().example()
                sage: C._next_from_iterator(10) # indirect doctest
                11

            TODO: specify the behavior when ``obj`` is not in ``self``.
            """
            it = iter(self)
            el = it.next()
            while el != obj:
                el = it.next()
            return it.next()
        next = _next_from_iterator


        def _unrank_from_iterator(self, r):
            """
            The ``r``-th element of ``self``

            ``self.unrank(r)`` returns the ``r``-th element of self where
            ``r`` is an integer between ``0`` and ``n-1`` where ``n`` is the
            cardinality of ``self``.

            This is the default (brute force) implementation from the
            category ``EnumeratedSets()`` which can be used when the
            method ``__iter__`` is provided. Its complexity is `O(r)`,
            where `r` is the rank of ``obj``.

            EXAMPLES::

                sage: C = FiniteEnumeratedSets().example()
                sage: C.unrank(2) # indirect doctest
                3
                sage: C.unrank(5) # indirect doctest
                Traceback (most recent call last):
                ...
                ValueError: the value must be between 0 and 2 inclusive
            """
            counter = 0
            for u in self:
                if counter == r:
                    return u
                counter += 1
            raise ValueError, "the value must be between %s and %s inclusive"%(0,counter-1)
        unrank = _unrank_from_iterator

        def _rank_from_iterator(self, x):
            """
            The rank of an element of ``self``

            ``self.unrank(x)`` returns the rank of `x`, that is its
            position in the enumeration of ``self``. This is an
            integer between ``0`` and ``n-1`` where ``n`` is the
            cardinality of ``self``, or None if `x` is not in `self`.

            This is the default (brute force) implementation from the
            category ``EnumeratedSet()`` which can be used when the
            method ``__iter__`` is provided. Its complexity is `O(r)`,
            where `r` is the rank of ``obj``. For infinite enumerated
            sets, this won't terminate when `x` is not in ``self``

            EXAMPLES::

                sage: C = FiniteEnumeratedSets().example()
                sage: list(C)
                [1, 2, 3]
                sage: C.rank(3) # indirect doctest
                2
                sage: C.rank(5) # indirect doctest
            """
            counter = 0
            for u in self:
                if u == x:
                    return counter
                counter += 1
            return None

        rank = _rank_from_iterator

        def _iterator_from_list(self):
            """
            An iterator for the elements of ``self``.

            ``iter(self)`` returns an iterator for the elements
            of ``self``. This is a generic implementation from the
            category ``EnumeratedSet()`` which can be used when the
            method ``list`` is provided.

            EXAMPLES::

                sage: C = FiniteEnumeratedSets().example()
                sage: it = C._iterator_from_list()
                sage: [it.next(), it.next(), it.next()]
                [1, 2, 3]
            """
            for x in self.list():
                yield x

        def _iterator_from_next(self):
            """
            An iterator for the elements of ``self``.

            ``iter(self)`` returns an iterator for the element of
            the set ``self``. This is a generic implementation from
            the category ``EnumeratedSet()`` which can be used when
            the methods ``first`` and ``next`` are provided.

            EXAMPLES::

                sage: C = InfiniteEnumeratedSets().example()
                sage: it = C._iterator_from_next()
                sage: [it.next(), it.next(), it.next(), it.next(), it.next()]
                [0, 1, 2, 3, 4]
            """
            f = self.first()
            yield f
            while True:
                try:
                    f = self.next(f)
                except (TypeError, ValueError ):
                    break

                if f is None or f is False :
                    break
                else:
                    yield f

        def _iterator_from_unrank(self):
            """
            An iterator for the elements of ``self``.

            ``iter(self)`` returns an iterator for the elements
            of the set ``self``. This is a generic implementation from
            the category ``EnumeratedSet()`` which can be used when
            the method ``unrank`` is provided.

            EXAMPLES::

                sage: C = InfiniteEnumeratedSets().example()
                sage: it = C._iterator_from_unrank()
                sage: [it.next(), it.next(), it.next(), it.next(), it.next()]
                [0, 1, 2, 3, 4]
            """
            r = 0
            u = self.unrank(r)
            yield u
            while True:
                r += 1
                try:
                    u = self.unrank(r)
                except (TypeError, ValueError, IndexError):
                    break

                if u == None:
                    break
                else:
                    yield u

        # This @cached_method also allows for the temporary assignments
        # an_element = EnumeratedSets.ParentMethods._an_element_
        # in parents in this category (see e.g. examples/FiniteEnumeratedSets
        @cached_method
        def _an_element_from_iterator(self):
            """
            An element in ``self``.

            ``self.an_element()`` returns a particular element of the set
            ``self``. This is a generic implementation from the category
            ``EnumeratedSet()`` which can be used when the method ``__iter__``
            is provided.

            EXAMPLES::

                sage: C = FiniteEnumeratedSets().example()
                sage: C.an_element() # indirect doctest
                1
            """
            it = self.__iter__()
            return it.next()
        # Should this implemented from first instead?
        _an_element_ = _an_element_from_iterator

        #FIXME: use combinatorial_class_from_iterator once class_from_iterator.patch is in
        def _some_elements_from_iterator(self):
            """
            Returns some elements in ``self``.

            See :class:`TestSuite` for a typical use case.

            This is a generic implementation from the category
            ``EnumeratedSet()`` which can be used when the method
            ``__iter__`` is provided. It returns an iterator for up to
            the first 100 elements of ``self``

            EXAMPLES::

                sage: C = FiniteEnumeratedSets().example()
                sage: list(C.some_elements()) # indirect doctest
                [1, 2, 3]
            """
            nb = 0
            for i in self:
                yield i
                nb += 1
                if nb >= 100:
                    break
        some_elements = _some_elements_from_iterator

        def random_element(self):
            """
            Returns a random element in ``self``.

            Unless otherwise stated, and for finite enumerated sets,
            the probability is uniform.

            This is a generic implementation from the category
            ``EnumeratedSet()``. It raise a ``NotImplementedError``
            since one does not know whether the set is finite.

            EXAMPLES::

                sage: class broken(UniqueRepresentation, Parent):
                ...    def __init__(self):
                ...        Parent.__init__(self, category = EnumeratedSets())
                ...
                sage: broken().random_element()
                Traceback (most recent call last):
                ...
                NotImplementedError: unknown cardinality
                """
            raise NotImplementedError, "unknown cardinality"


#
#  Consistency test suite for an enumerated set:
#
        # If the enumerated set is large, one can stop some coeherency tests
        # after looking at the first element by setting the following variable:
        max_test_enumerated_set_loop=100 # TODO: Devise a sensible bound !!
        ##########
        def _test_enumerated_set_contains(self, **options):
            """
            Checks that the methods :meth:`.__contains__` and :meth:`.__iter__` are consistent.

            See also :class:`TestSuite`.

            TESTS::

                sage: C = FiniteEnumeratedSets().example()
                sage: C._test_enumerated_set_contains()
                sage: TestSuite(C).run()

            Let us now break the class::

                sage: from sage.categories.examples.finite_enumerated_sets import Example
                sage: class CCls(Example):
                ...       def __contains__(self, obj):
                ...           if obj == 3:
                ...               return False
                ...           else:
                ...               return obj in C
                sage: CC = CCls()
                sage: TestSuite(CC).run()
                Traceback (most recent call last):
                ...
                AssertionError
            """
            tester = self._tester(**options)
            i = 0
            for w in self:
                tester.assertTrue(w in self)
                i += 1
                if i > self.max_test_enumerated_set_loop:
                    return

        def _test_enumerated_set_iter_list(self, **options):
            """
            Checks that the methods :meth:`.list` and :meth:`.__iter__` are consistent.

            See also: :class:`TestSuite`.

            EXAMPLES::

                sage: C = FiniteEnumeratedSets().example()
                sage: C._test_enumerated_set_iter_list()
                sage: TestSuite(C).run()

            Let us now break the class::

                sage: from sage.categories.examples.finite_enumerated_sets import Example
                sage: class CCls(Example):
                ...       def list(self):
                ...           return [1,2,3,4]
                sage: CC = CCls()
                sage: TestSuite(CC).run()
                Traceback (most recent call last):
                ...
                AssertionError: 3 != 4
            """
            tester = self._tester(**options)
            if self.list != self._list_default:
                ls = self.list()
                i = 0
                for obj in self:
                    tester.assertEqual(obj, ls[i])
                    i += 1
                    if i > self.max_test_enumerated_set_loop:
                        return
                tester.assertEqual(i, len(ls))

    class ElementMethods:

        def rank(self):
            """
            Returns the rank of ``self`` in its parent.

            See also :meth:`EnumeratedSets.ElementMethods.rank`

            EXAMPLES::

                sage: F = FiniteSemigroups().example(('a','b','c'))
                sage: L = list(F); L
                ['a', 'c', 'ac', 'b', 'ba', 'bc', 'cb', 'ca', 'bca', 'ab', 'bac', 'cab', 'acb', 'cba', 'abc']
                sage: L[7].rank()
                7
            """
            return self.parent().rank(self)
