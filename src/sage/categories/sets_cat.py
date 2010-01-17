r"""
Sets
"""
#*****************************************************************************
#  Copyright (C) 2005      David Kohel <kohel@maths.usyd.edu>
#                          William Stein <wstein@math.ucsd.edu>
#                2008      Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.sage_unittest import TestSuite
from sage.misc.abstract_method import abstract_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.categories.category import Category, HomCategory
# Do not use sage.categories.all here to avoid initialization loop
from sage.categories.objects import Objects

class Sets(Category):
    """
    The category of sets

    The base category for collections of elements with = (equality)

    This is also the category whose objects are all parents.

    EXAMPLES::

        sage: Sets()
        Category of sets
        sage: Sets().super_categories()
        [Category of objects]
        sage: Sets().all_super_categories()
        [Category of sets, Category of objects]

    Let us consider an example of set::

        sage: P = Sets().example("inherits")
        sage: P
        Set of prime numbers

    See ``P??`` for the code.


    P is in the category of sets::

        sage: P.category()
        Category of sets

    and therefore gets its methods from the following classes::

        sage: for cl in P.__class__.mro(): print(cl)
        <class 'sage.categories.examples.sets_cat.PrimeNumbers_Inherits_with_category'>
        <class 'sage.categories.examples.sets_cat.PrimeNumbers_Inherits'>
        <class 'sage.categories.examples.sets_cat.PrimeNumbers_Abstract'>
        <class 'sage.structure.unique_representation.UniqueRepresentation'>
        <type 'sage.structure.parent.Parent'>
        <type 'sage.structure.category_object.CategoryObject'>
        <type 'sage.structure.sage_object.SageObject'>
        <class 'sage.categories.sets_cat.Sets.parent_class'>
        <class 'sage.categories.objects.Objects.parent_class'>
        <type 'object'>

    We run some generic checks on P::

        sage: TestSuite(P).run(verbose=True)
        running ._test_an_element() . . . pass
        running ._test_category() . . . pass
        running ._test_elements() . . .
          Running the test suite of self.an_element()
          running ._test_category() . . . pass
          running ._test_not_implemented_methods() . . . pass
          running ._test_pickling() . . . pass
          pass
        running ._test_not_implemented_methods() . . . pass
        running ._test_pickling() . . . pass
        running ._test_some_elements() . . . pass

    Now, we manipulate some elements of P::

        sage: P.an_element()
        47
        sage: x = P(3)
        sage: x.parent()
        Set of prime numbers
        sage: x in P, 4 in P
        (True, False)
        sage: x.is_prime()
        True

    They get their methods from the following classes::

        sage: for cl in x.__class__.mro(): print(cl)
        <class 'sage.categories.examples.sets_cat.PrimeNumbers_Inherits_with_category.element_class'>
        <class 'sage.categories.examples.sets_cat.PrimeNumbers_Inherits.Element'>
        <type 'sage.rings.integer.IntegerWrapper'>
        <type 'sage.rings.integer.Integer'>
        <type 'sage.structure.element.EuclideanDomainElement'>
        <type 'sage.structure.element.PrincipalIdealDomainElement'>
        <type 'sage.structure.element.DedekindDomainElement'>
        <type 'sage.structure.element.IntegralDomainElement'>
        <type 'sage.structure.element.CommutativeRingElement'>
        <type 'sage.structure.element.RingElement'>
        <type 'sage.structure.element.ModuleElement'>
        <class 'sage.categories.examples.sets_cat.PrimeNumbers_Abstract.Element'>
        <type 'sage.structure.element.Element'>
        <type 'sage.structure.sage_object.SageObject'>
        <class 'sage.categories.sets_cat.Sets.element_class'>
        <class 'sage.categories.objects.Objects.element_class'>
        <type 'object'>

    FIXME: Objects.element_class is not very meaningfull ...


    TESTS::

          sage: TestSuite(Sets()).run()

    """

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: Sets().super_categories()
            [Category of objects]
        """
        return [Objects()]

    def __call__(self, X):
        """
        Construct an object in this category from the data in ``X``.

        EXAMPLES::

            sage: Sets()(ZZ)
            Set of elements of Integer Ring

        FIXME: the above behavior dates back from the first category
        writeup. It is not consistent with :meth:`Category.__call__`.
        Should we change it to just return ``ZZ`` instead?
        """
        import sage.sets.all
        return sage.sets.all.Set(X)

    def example(self, choice = None):
        """
        Returns examples of objects of ``Sets()``, as per :meth:`Category.example`.

        EXAMPLES::

            sage: Sets().example()
            Set of prime numbers (basic implementation)

            sage: Sets().example("inherits")
            Set of prime numbers

            sage: Sets().example("facade")
            Set of prime numbers (facade implementation)

            sage: Sets().example("wrapper")
            Set of prime numbers (wrapper implementation)
        """
        if choice is None:
            from sage.categories.examples.sets_cat import PrimeNumbers
            return PrimeNumbers()
        elif choice == "inherits":
            from sage.categories.examples.sets_cat import PrimeNumbers_Inherits
            return PrimeNumbers_Inherits()
        elif choice == "facade":
            from sage.categories.examples.sets_cat import PrimeNumbers_Facade
            return PrimeNumbers_Facade()
        elif choice == "wrapper":
            from sage.categories.examples.sets_cat import PrimeNumbers_Wrapper
            return PrimeNumbers_Wrapper()
        else:
            raise ValueError, "Unkown choice"

    class ParentMethods:
#         # currently overriden by the default implementation in sage.structure.Parent
#         def __call__(self, *args, **options):
#             return self.element_class(*args, **options)

        # Todo: simplify the _element_constructor definition logic
        # Todo: find a nicer mantra for conditionaly defined methods
        @lazy_attribute
        def _element_constructor_(self):
            r"""
            TESTS::

                sage: S = Sets().example()
                sage: S._element_constructor_(17)
                17
                sage: S(17) # indirect doctest
                17

            Caveat: for some parents, element_class is a method, and
            not an attribute. We do not provide a default
            implementation of _element_constructor for those.

                sage: FreeModule(QQ,3).element_class
                <bound method FreeModule_ambient_field_with_category.element_class of Vector space of dimension 3 over Rational Field>
                sage: FreeModule(QQ,3)._element_constructor
            """
            if hasattr(self, "element_class") and issubclass(self.element_class, object):
                return self._element_constructor_from_element_class
            else:
                return NotImplemented

        def _element_constructor_from_element_class(self, *args, **keywords):
            """
            The default constructor for elements of this parent

            Among other things, it is called upon my_parent(data) when
            the coercion model did not find a way to coerce data into
            this parent.

            This default implementation for
            :meth:`_element_constructor_` calls the constructor of the
            element class.

            Caveat: ``self`` is passed to the constructor of the
            element class as a keyword argument ``parent``. Many
            element classes in Sage, in particular those implemented
            by mean of extension types, take ``parent`` as first
            mandatory argument instead.

            This incompatibility will be fixed soon (Fall 2009?) by
            having all element classes take ``parent`` as first
            mandatory argument, and updating this default
            implementation of :meth:`_element_constructor_`.

            EXAMPLES::

                sage: S = Sets().example("inherits")
                sage: S._element_constructor_from_element_class(17)
                17
            """
            return self.element_class(parent = self, *args, **keywords)

        @abstract_method
        def __contains__(self, x):
            """
            Tests whether the set contains the object ``x``.

            All parents in the category ``Sets()`` should implement this method.

            TESTS::

                sage: P = sage.categories.examples.sets_cat.PrimeNumbers()
                sage: 12 in P
                False
                sage: P(5) in P
                True
            """


        @cached_method
        def an_element(self):
            r"""
            Returns a (preferably typical) element of this parent.

            This is used both for illustration and testing purposes.

            This default implementation calls :meth:`_an_element_` and
            cache the result. Any parent should implement either
            :meth:`an_element` or :meth:`_an_element_`.

            EXAMPLES::

               sage: CDF.an_element()
               1.0*I
               sage: ZZ[['t']].an_element()
               t
            """
            return self._an_element_()

        def _test_an_element(self, **options):
            """
            Run generic tests on the method :meth:`.an_element`.

            See also: :class:`TestSuite`.

            EXAMPLES::

                sage: C = Sets().example()
                sage: C._test_an_element()

            Let us now write a broken :meth:`.an_element` method::

                sage: from sage.categories.examples.sets_cat import PrimeNumbers
                sage: class CCls(PrimeNumbers):
                ...       def an_element(self):
                ...           return 18
                sage: CC = CCls()
                sage: CC._test_an_element()
                Traceback (most recent call last):
                ...
                AssertionError: self.an_element() is not in self
            """
            tester = self._tester(**options)
            an_element = self.an_element()
            tester.assert_(an_element in self, "self.an_element() is not in self")

        def _test_elements(self, tester = None, **options):
            """
            Run generic tests on element(s) of ``self``.

            See also: :class:`TestSuite`.

            EXAMPLES::

                sage: C = Sets().example()
                sage: C._test_elements(verbose = True)
                <BLANKLINE>
                  Running the test suite of self.an_element()
                  running ._test_category() . . . pass
                  running ._test_not_implemented_methods() . . . pass
                  running ._test_pickling() . . . pass
                <BLANKLINE>

            Debugging tip: in case of failure of this test, run instead:

                sage: TestSuite(C.an_element()).run()

            Let us now implement a parent whose elements cannot be pickled::

                sage: from sage.categories.examples.sets_cat import PrimeNumbers
                sage: class Bla(SageObject): pass
                sage: class CCls(PrimeNumbers):
                ...       def an_element(self):
                ...           return Bla()
                sage: CC = CCls()
                sage: CC._test_elements()
                  Failure in _test_pickling:
                  ...
                  PicklingError: Can't pickle <class '__main__.Bla'>: attribute lookup __main__.Bla failed
                  ...
                  The following tests failed: _test_pickling
            """
            # TODO: add native support for nested test suites to TestSuite

            # The intention is to raise an exception only if this is
            # run as a sub-testsuite of a larger testsuite.
            is_sub_testsuite = (tester is not None)
            tester = self._tester(tester = tester, **options)
            # Or do we want to run the test on some_elements?
            an_element = self.an_element()
            tester.info("\n  Running the test suite of self.an_element()")
            TestSuite(an_element).run(verbose = tester._verbose, prefix = tester._prefix+"  ",
                                      raise_on_failure = is_sub_testsuite)
            tester.info(tester._prefix+"  ", newline = False)

        def some_elements(self):
            """
            Returns a list (or iterable) of elements of self.

            This is typically used for running generic tests (see :class:`TestSuite`).

            This default implementation calls :meth:`.an_element`.

            EXAMPLES::

                sage: S = Sets().example(); S
                Set of prime numbers (basic implementation)
                sage: S.an_element()
                47
                sage: S.some_elements()
                [47]

            This method should return an iterable, *not* an iterator.
            """
            return [ self.an_element() ]

        def _test_some_elements(self, **options):
            """
            Run generic tests on the method :meth:`.some_elements`.

            See also: :class:`TestSuite`.

            EXAMPLES::

                sage: C = Sets().example()
                sage: C._test_some_elements()

            Let us now write a broken :meth:`.some_elements` method::

                sage: from sage.categories.examples.sets_cat import PrimeNumbers
                sage: class CCls(PrimeNumbers):
                ...       def some_elements(self):
                ...           return [self(17), 32]
                sage: CC = CCls()
                sage: CC._test_some_elements()
                Traceback (most recent call last):
                ...
                AssertionError: the object 32 in self.some_elements() is not in self
            """
            tester = self._tester(**options)
            elements = self.some_elements()
            # Todo: enable this once
            #tester.assert_(elements != iter(elements),
            #               "self.some_elements() should return an iterable, not an iterator")
            for x in elements:
                tester.assert_(x in self, "the object %s in self.some_elements() is not in self"%(x,))

    class ElementMethods:
        ##def equal(x,y):
        ##def =(x,y):

        # Used by Element._test_category
        _dummy_attribute = None

    class HomCategory(HomCategory):
        pass
