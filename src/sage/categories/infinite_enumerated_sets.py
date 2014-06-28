r"""
Infinite Enumerated Sets

AUTHORS:

 - Florent Hivert (2009-11): initial revision.

"""
#*****************************************************************************
#  Copyright (C) 2009 Florent Hivert <Florent.Hivert@univ-rouen.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************


from sage.categories.category_with_axiom import CategoryWithAxiom

class InfiniteEnumeratedSets(CategoryWithAxiom):
    """
    The category of infinite enumerated sets

    An infinite enumerated sets is a countable set together with a
    canonical enumeration of its elements.

    EXAMPLES::

        sage: InfiniteEnumeratedSets()
        Category of infinite enumerated sets
        sage: InfiniteEnumeratedSets().super_categories()
        [Category of enumerated sets, Category of infinite sets]
        sage: InfiniteEnumeratedSets().all_super_categories()
        [Category of infinite enumerated sets,
         Category of enumerated sets,
         Category of infinite sets,
         Category of sets,
         Category of sets with partial maps,
         Category of objects]

    TESTS::

        sage: C = InfiniteEnumeratedSets()
        sage: TestSuite(C).run()
    """

    class ParentMethods:

        def random_element(self):
            """
            Returns an error since self is an infinite enumerated set.

            EXAMPLES::

                sage: NN = InfiniteEnumeratedSets().example()
                sage: NN.random_element()
                Traceback (most recent call last):
                ...
                NotImplementedError: infinite set

            TODO: should this be an optional abstract_method instead?
            """
            raise NotImplementedError("infinite set")

        def list(self):
            """
            Returns an error since self is an infinite enumerated set.

            EXAMPLES::

                sage: NN = InfiniteEnumeratedSets().example()
                sage: NN.list()
                Traceback (most recent call last):
                ...
                NotImplementedError: infinite list
            """
            raise NotImplementedError("infinite list")
        _list_default  = list # needed by the check system.

        def _test_enumerated_set_iter_cardinality(self, **options):
            """
            Check that the methods :meth:`.cardinality` and
            :meth:`.__iter__`. are consistent.

            See also :class:`TestSuite`.

            For infinite enumerated sets:

               * :meth:`.cardinality` is supposed to return `infinity`

               * :meth:`.list`` is supposed to raise a ``NotImplementedError``.

            EXAMPLES::

                sage: NN = InfiniteEnumeratedSets().example()
                sage: NN._test_enumerated_set_iter_cardinality()
            """
            tester = self._tester(**options)
            from sage.rings.infinity import infinity
            tester.assertEqual(self.cardinality(), infinity)
            tester.assertRaises(NotImplementedError, self.list)
