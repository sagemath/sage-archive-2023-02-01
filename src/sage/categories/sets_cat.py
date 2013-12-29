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

from sage.misc.lazy_import import lazy_import, LazyImport
from sage.misc.cachefunc import cached_method
from sage.misc.sage_unittest import TestSuite
from sage.misc.abstract_method import abstract_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.lazy_format import LazyFormat
from sage.categories.category import Category, HomCategory
from sage.categories.category_singleton import Category_singleton
# Do not use sage.categories.all here to avoid initialization loop
from sage.categories.sets_with_partial_maps import SetsWithPartialMaps
from sage.categories.subquotients import SubquotientsCategory
from sage.categories.quotients    import QuotientsCategory
from sage.categories.subobjects   import SubobjectsCategory
from sage.categories.isomorphic_objects   import IsomorphicObjectsCategory
from sage.categories.algebra_functor import AlgebrasCategory
from sage.categories.cartesian_product import cartesian_product, CartesianProductsCategory
from sage.categories.realizations import RealizationsCategory, Category_realization_of_parent
from sage.categories.with_realizations import WithRealizationsCategory

lazy_import('sage.sets.cartesian_product', 'CartesianProduct')

def print_compare(x, y):
    """
    Helper method used in
    :meth:`Sets.ParentMethods._test_elements_eq_symmetric`,
    :meth:`Sets.ParentMethods._test_elements_eq_tranisitive`.

    INPUT:

    - ``x`` -- an element

    - ``y`` -- an element

    EXAMPLES::

        sage: from sage.categories.sets_cat import print_compare
        sage: print_compare(1,2)
        1 != 2
        sage: print_compare(1,1)
        1 == 1

    """
    if x == y:
        return LazyFormat("%s == %s")%(x, y)
    else:
        return LazyFormat("%s != %s")%(x, y)

class EmptySetError(ValueError):
    """
    Exception raised when some operation can't be performed on the empty set.

    EXAMPLES::

        sage: def first_element(st):
        ...    if not st: raise EmptySetError, "no elements"
        ...    else: return st[0]
        sage: first_element(Set((1,2,3)))
        1
        sage: first_element(Set([]))
        Traceback (most recent call last):
        ...
        EmptySetError: no elements
    """
    pass

class Sets(Category_singleton):
    r"""
    The category of sets

    The base category for collections of elements with = (equality)

    This is also the category whose objects are all parents.

    EXAMPLES::

        sage: Sets()
        Category of sets
        sage: Sets().super_categories()
        [Category of sets with partial maps]
        sage: Sets().all_super_categories()
        [Category of sets, Category of sets with partial maps, Category of objects]

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
        <class 'sage.structure.unique_representation.CachedRepresentation'>
        <type 'sage.misc.fast_methods.WithEqualityById'>
        <type 'sage.structure.parent.Parent'>
        <type 'sage.structure.category_object.CategoryObject'>
        <type 'sage.structure.sage_object.SageObject'>
        <class 'sage.categories.sets_cat.Sets.parent_class'>
        <class 'sage.categories.category.SetsWithPartialMaps.parent_class'>
        <class 'sage.categories.objects.Objects.parent_class'>
        <type 'object'>

    We run some generic checks on P::

        sage: TestSuite(P).run(verbose=True)
        running ._test_an_element() . . . pass
        running ._test_category() . . . pass
        running ._test_elements() . . .
          Running the test suite of self.an_element()
          running ._test_category() . . . pass
          running ._test_eq() . . . pass
          running ._test_not_implemented_methods() . . . pass
          running ._test_pickling() . . . pass
          pass
        running ._test_elements_eq_reflexive() . . . pass
        running ._test_elements_eq_symmetric() . . . pass
        running ._test_elements_eq_transitive() . . . pass
        running ._test_elements_neq() . . . pass
        running ._test_eq() . . . pass
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
        <class 'sage.categories.category.SetsWithPartialMaps.element_class'>
        <class 'sage.categories.objects.Objects.element_class'>
        <type 'object'>

    FIXME: Objects.element_class is not very meaningful ...


    TESTS::

          sage: TestSuite(Sets()).run()

    """

    def super_categories(self):
        r"""
        We include SetsWithPartialMaps between Sets and Objects so that we
        can define morphisms between sets that are only partially defined.
        This is also to have the Homset constructor not complain that
        SetsWithPartialMaps is not a supercategory of Fields, for example.

        EXAMPLES::

            sage: Sets().super_categories()
            [Category of sets with partial maps]
        """
        return [SetsWithPartialMaps()]

    def _call_(self, X, enumerated_set_if_possible = False):
        r"""
        Construct an object in this category from the data in ``X``.

        EXAMPLES::

            sage: Sets()(ZZ)
            Integer Ring
            sage: Sets()([1, 2, 3])
            {1, 2, 3}

        .. note::

           Using ``Sets()(A)`` used to implement some sort of
           forgetful functor into the ``Sets()`` category. This
           feature has been removed, because it was not consistent
           with the semantic of :meth:`Category.__call__`. Proper
           forgetful functors will eventually be implemented, with
           another syntax.

        - ``enumerated_set_if_possible`` -- an option to ask Sage to
          try to build an ``EnumeratedSets()`` rather that a
          ``Sets()`` if possible.  This is experimental an may change
          in the future::

              sage: S = Sets()([1, 2, 3]); S.category()
              Category of sets
              sage: S = Sets()([1, 2, 3], True); S.category()
              Category of facade finite enumerated sets
        """
        import sage.sets.all
        if enumerated_set_if_possible:
            from sage.categories.enumerated_sets import EnumeratedSets
            try:
                return EnumeratedSets()(X)
            except NotImplementedError:
                pass
        return sage.sets.all.Set(X)

    def example(self, choice = None):
        """
        Returns examples of objects of ``Sets()``, as per
        :meth:`Category.example()
        <sage.categories.category.Category.example>`.

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

            Caveat: For some parents, element_class is a method, and
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
            return self.element_class(self, *args, **keywords)

        def is_parent_of(self, element):
            """
            Returns whether ``self`` is the parent of ``element``

            INPUT:
             - ``element`` -- any object

            EXAMPLES::

                sage: S = ZZ
                sage: S.is_parent_of(1)
                True
                sage: S.is_parent_of(2/1)
                False

            This method differs from :meth:`__contains__` because it
            does not attempt any coercion::

                sage: 2/1 in S, S.is_parent_of(2/1)
                (True, False)
                sage: int(1) in S, S.is_parent_of(int(1))
                (True, False)
            """
            from sage.structure.element import parent
            return parent(element) == self

        @abstract_method
        def __contains__(self, x):
            """
            Tests whether the set contains the object ``x``.

            All parents in the category ``Sets()`` should implement this method.

            TESTS::

                sage: from sage.categories.examples.sets_cat import PrimeNumbers
                sage: P = PrimeNumbers()
                sage: 12 in P
                False
                sage: P(5) in P
                True
            """

        @cached_method
        def an_element(self):
            r"""
            Returns a (preferably typical) element of this parent.

            This is used both for illustration and testing purposes. If the
            set ``self`` is empty, :meth:`an_element` should raise the exception
            :class:`EmptySetError`.

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

            TESTS::

                sage: FiniteEnumeratedSet([])._test_an_element()
            """
            tester = self._tester(**options)
            try:
                an_element = self.an_element()
            except EmptySetError:
                return
            tester.assertTrue(an_element in self, "self.an_element() is not in self")
#            tester.assertTrue(self.is_parent_of(an_element), "self is not the parent of self.an_element()")
#            tester.assertEqual(self(an_element), an_element, "element construction is not idempotent")
            if self.is_parent_of(an_element):
                tester.assertEqual(self(an_element), an_element, "element construction is not idempotent")
            else: # Allows self(an_element) to fails for facade parent.
                try:
                    rebuilt_element = self(an_element)
                except NotImplementedError:
                    tester.info("\n  The set doesn't seems to implement __call__; skipping test of construction idempotency")
                    pass
                else:
                    tester.assertEqual(rebuilt_element, an_element, "element construction is not idempotent")


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
                  running ._test_eq() . . . pass
                  running ._test_nonzero_equal() . . . pass
                  running ._test_not_implemented_methods() . . . pass
                  running ._test_pickling() . . . pass
                <BLANKLINE>

            Debugging tip: in case of failure of this test, run instead::

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
            try:
                an_element = self.an_element()
            except EmptySetError:
                return
            tester.info("\n  Running the test suite of self.an_element()")
            TestSuite(an_element).run(verbose = tester._verbose, prefix = tester._prefix+"  ",
                                      raise_on_failure = is_sub_testsuite)
            tester.info(tester._prefix+" ", newline = False)

        def _test_elements_eq_reflexive(self, **options):
            """
            Run generic tests on the equality of elements.

            Test that ``==`` is reflexive.

            See also: :class:`TestSuite`.

            EXAMPLES::

                sage: C = Sets().example()
                sage: C._test_elements_eq_reflexive()

            We try a non-reflexive equality::

                sage: P = Sets().example("wrapper")
                sage: P._test_elements_eq_reflexive()
                sage: eq = P.element_class.__eq__

                sage: P.element_class.__eq__ = (lambda x, y:
                ...        False if eq(x, P(47)) and eq(y, P(47)) else eq(x, y))
                sage: P._test_elements_eq_reflexive()
                Traceback (most recent call last):
                ...
                AssertionError: 47 != 47

            We restore ``P.element_class`` in a proper state for further tests::

                sage: P.element_class.__eq__ = eq

            """
            tester = self._tester(**options)
            S = list(tester.some_elements()) + [None, 0]
            for x in S:
                tester.assertEqual(x, x)

        def _test_elements_eq_symmetric(self, **options):
            """
            Run generic tests on the equality of elements.

            This tests that ``==`` is symmetric.

            See also: :class:`TestSuite`.

            EXAMPLES::

                sage: C = Sets().example()
                sage: C._test_elements_eq_symmetric()

            We test a non symmetric equality::

                sage: P = Sets().example("wrapper")
                sage: P._test_elements_eq_symmetric()
                sage: eq = P.element_class.__eq__

                sage: def non_sym_eq(x, y):
                ...      if not y in P:                      return False
                ...      elif eq(x, P(47)) and eq(y, P(53)): return True
                ...      else:                               return eq(x, y)
                sage: P.element_class.__eq__ = non_sym_eq
                sage: P._test_elements_eq_symmetric()
                Traceback (most recent call last):
                ...
                AssertionError: non symmetric equality: 47 == 53 but 53 != 47

            We restore ``P.element_class`` in a proper state for further tests::

                sage: P.element_class.__eq__ = eq

            """
            tester = self._tester(**options)
            S = list(tester.some_elements()) + [None, 0]
            n = tester._max_runs
            from sage.combinat.cartesian_product import CartesianProduct
            S = CartesianProduct(S,S)
            if len(S) > n:
                from random import sample
                S = sample(S, n)

            for x, y in S:
                tester.assertEqual(x==y, y==x,
                    LazyFormat("non symmetric equality: %s but %s")%(
                        print_compare(x, y), print_compare(y, x)))

        def _test_elements_eq_transitive(self, **options):
            """
            Run generic tests on the equality of elements.

            Test that ``==`` is transitive.

            See also: :class:`TestSuite`.

            EXAMPLES::

                sage: C = Sets().example()
                sage: C._test_elements_eq_transitive()

            We test a non transitive equality::

                sage: R = Zp(3)
                sage: Sets().ParentMethods._test_elements_eq_transitive.im_func(R,elements=[R(3,2),R(3,1),R(0)])
                Traceback (most recent call last):
                ...
                AssertionError: non transitive equality:
                3 + O(3^2) == O(3) and O(3) == 0 but 3 + O(3^2) != 0

            """
            tester = self._tester(**options)
            S = list(tester.some_elements())
            n = tester._max_runs
            if (len(S)+2)**3 <= n:
                S = list(S) + [None, 0]
            else:
                from random import sample
                from sage.rings.integer import Integer
                S = sample(S, Integer(n).nth_root(3,truncate_mode=1)[0] - 2) + [None, 0]

            for x in S:
                for y in S:
                    if not x == y: continue
                    for z in S:
                        if not y == z: continue
                        tester.assertTrue(x == z,
                            LazyFormat("non transitive equality:\n"
                                       "%s and %s but %s")%(
                                print_compare(x, y),
                                print_compare(y, z),
                                print_compare(x, z)))

        def _test_elements_neq(self, **options):
            """
            Run generic tests on the equality of elements.

            Test that ``==`` and ``!=`` are consistent.

            See also: :class:`TestSuite`.

            EXAMPLES::

                sage: C = Sets().example()
                sage: C._test_elements_neq()

            We try a broken inequality::

                sage: P = Sets().example("wrapper")
                sage: P._test_elements_neq()
                sage: ne = P.element_class.__ne__
                sage: eq = P.element_class.__eq__

                sage: P.element_class.__ne__ = lambda x, y: False
                sage: P._test_elements_neq()
                Traceback (most recent call last):
                ...
                AssertionError: __eq__ and __ne__ inconsistency:
                  47 == 53 returns False  but  47 != 53 returns False

                sage: P.element_class.__ne__ = lambda x, y: not(x == y)

            We restore ``P.element_class`` in a proper state for further tests::

                sage: P.element_class.__ne__ = ne
                sage: P.element_class.__eq__ = eq
            """
            tester = self._tester(**options)
            S = list(tester.some_elements()) + [None, 0]
            n = tester._max_runs
            from sage.combinat.cartesian_product import CartesianProduct
            S = CartesianProduct(S,S)
            if len(S) > n:
                from random import sample
                S = sample(S, n)

            for x,y in S:
                tester.assertNotEqual(x == y, x != y,
                    LazyFormat("__eq__ and __ne__ inconsistency:\n"
                        "  %s == %s returns %s  but  %s != %s returns %s")%(
                            x, y, (x == y), x, y, (x != y)))

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
                sage: S = Set([])
                sage: S.some_elements()
                []

            This method should return an iterable, *not* an iterator.
            """
            try:
                return [ self.an_element() ]
            except EmptySetError:
                return []

        def _test_some_elements(self, **options):
            """
            Run generic tests on the method :meth:`.some_elements`.

            See also: :class:`TestSuite`.

            EXAMPLES::

                sage: C = Sets().example()
                sage: C._test_some_elements()

            Let us now write a broken :meth:`.some_elements` method::

                sage: from sage.categories.examples.sets_cat import *
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
                tester.assertTrue(x in self, LazyFormat(
                    "the object %s in self.some_elements() is not in self")%(x,))

        # Functorial constructions

        CartesianProduct = CartesianProduct
        def cartesian_product(*parents):
            """
            Returns the cartesian product of the parents

            EXAMPLES::

                sage: C = AlgebrasWithBasis(QQ)
                sage: A = C.example(); A.rename("A")
                sage: A.cartesian_product(A,A)
                A (+) A (+) A
            """
            return parents[0].__class__.CartesianProduct(parents, category = cartesian_product.category_from_parents(parents))

        def algebra(self, base_ring, category = None):
            """
            Returns the algebra of self

            EXAMPLES::

                sage: A = Sets().example().algebra(QQ); A
                Free module generated by Set of prime numbers (basic implementation) over Rational Field
                sage: A.category()
                Category of set algebras over Rational Field

                sage: A = Monoids().example().algebra(QQ); A
                Free module generated by An example of a monoid: the free monoid generated by ('a', 'b', 'c', 'd') over Rational Field
                sage: A.category()
                Category of monoid algebras over Rational Field

                sage: A = Monoids().example().algebra(QQ, category = Sets()); A
                Free module generated by An example of a monoid: the free monoid generated by ('a', 'b', 'c', 'd') over Rational Field
                sage: A.category()
                Category of set algebras over Rational Field
            """
            if category is None:
                category = self.category()
            from sage.combinat.free_module import CombinatorialFreeModule
            return CombinatorialFreeModule(base_ring, self, category = category.Algebras(base_ring))


    class ElementMethods:
        ##def equal(x,y):
        ##def =(x,y):

        # Used by Element._test_category
        _dummy_attribute = None

        def cartesian_product(*elements):
            """
            Returns the cartesian product of its arguments, as an element of
            the cartesian product of the parents of those elements.

            EXAMPLES::

                sage: C = AlgebrasWithBasis(QQ)
                sage: A = C.example()
                sage: (a,b,c) = A.algebra_generators()
                sage: a.cartesian_product(b, c)
                B[(0, word: a)] + B[(1, word: b)] + B[(2, word: c)]

            FIXME: is this a policy that we want to enforce on all parents?
            """
            from sage.structure.element import parent, Element
            assert(all(isinstance(element, Element) for element in elements))
            parents = [parent(element) for element in elements]
            return cartesian_product(parents)._cartesian_product_of_elements(elements) # good name???


    class HomCategory(HomCategory):
        pass

    Facades = LazyImport('sage.categories.facade_sets', 'FacadeSets')

    class Subquotients(SubquotientsCategory):
        """
        A category for subquotients of sets

        EXAMPLES::

            sage: Sets().Subquotients()
            Category of subquotients of sets
            sage: Sets().Subquotients().all_super_categories()
            [Category of subquotients of sets, Category of sets,
             Category of sets with partial maps,
             Category of objects]
        """

        class ParentMethods:

            def _repr_(self):
                """
                EXAMPLES::

                    sage: from sage.categories.examples.semigroups import IncompleteSubquotientSemigroup
                    sage: S = IncompleteSubquotientSemigroup()
                    sage: S._repr_()
                    'A subquotient of An example of a semigroup: the left zero semigroup'
                """
                return "A subquotient of %s"%(self.ambient())

            @abstract_method
            def ambient(self):
                """
                Returns the ambient space for `self`.

                EXAMPLES::

                    sage: Semigroups().Subquotients().example().ambient()
                    An example of a semigroup: the left zero semigroup

                See :class:`~sage.category.subquotients.Subquotients` for
                specification and also :meth:`.lift` and :meth:`.retract`
                """

            # Should lift and retract be declared as conversions to the coercion mechanism ?
            # Compatibility issue: in IntegerModRing, lift is a method returning the actual
            # lifting morphism::
            #
            #   sage: F = IntegerModRing(3)
            #   sage: F.lift()
            #   Set-theoretic ring morphism:
            #   From: Ring of integers modulo 3
            #   To:   Integer Ring
            #   Defn: Choice of lifting map
            @abstract_method
            def lift(self, x):
                """
                INPUT:
                 - ``x`` -- an element of ``self``

                Lifts `x` to the ambient space for ``self``.

                EXAMPLES::

                    sage: S = Semigroups().Subquotients().example()
                    sage: s = S.an_element()
                    sage: s, s.parent()
                    (42, An example of a (sub)quotient semigroup: a quotient of the left zero semigroup)
                    sage: S.lift(s), S.lift(s).parent()
                    (42, An example of a semigroup: the left zero semigroup)
                    sage: s.lift(), s.lift().parent()
                    (42, An example of a semigroup: the left zero semigroup)

                See :class:`~sage.category.subquotients.Subquotients` for
                specification, :meth:`.ambient`, :meth:`.retract`, and also
                :meth:`Sets.Subquotients.ElementMethods.lift`
                """

            @abstract_method
            def retract(self, x):
                """
                INPUT:
                 - ``x`` -- an element of the ambient space for ``self``

                Retracts ``x`` to ``self``.

                See :class:`~sage.category.subquotients.Subquotients` for
                specification, :meth:`.ambient`, :meth:`.retract`, and also
                :meth:`Sets.Subquotients.ElementMethods.retract`

                EXAMPLES::

                    sage: S = Semigroups().Subquotients().example()
                    sage: s = S.ambient().an_element()
                    sage: s, s.parent()
                    (42, An example of a semigroup: the left zero semigroup)
                    sage: S.retract(s), S.retract(s).parent()
                    (42, An example of a (sub)quotient semigroup: a quotient of the left zero semigroup)
                """

        class ElementMethods:

            def lift(self):
                """
                Lifts ``self`` to the ambient space for its parent

                EXAMPLES::

                    sage: S = Semigroups().Subquotients().example()
                    sage: s = S.an_element()
                    sage: s, s.parent()
                    (42, An example of a (sub)quotient semigroup: a quotient of the left zero semigroup)
                    sage: S.lift(s), S.lift(s).parent()
                    (42, An example of a semigroup: the left zero semigroup)
                    sage: s.lift(), s.lift().parent()
                    (42, An example of a semigroup: the left zero semigroup)
                """
                return self.parent().lift(self)

    class Quotients(QuotientsCategory):
        """
        A category for subquotients of sets

        EXAMPLES::

            sage: Sets().Quotients()
            Category of quotients of sets
            sage: Sets().Quotients().all_super_categories()
            [Category of quotients of sets, Category of subquotients of sets, Category of sets,
             Category of sets with partial maps, Category of objects]
        """

        class ParentMethods:

            def _repr_(self):
                """
                EXAMPLES::

                    sage: from sage.categories.examples.semigroups import IncompleteSubquotientSemigroup
                    sage: S = IncompleteSubquotientSemigroup(category=Semigroups().Quotients())
                    sage: S._repr_()
                    'A quotient of An example of a semigroup: the left zero semigroup'
                """
                return "A quotient of %s"%(self.ambient())

            def _an_element_(self):
                """
                Returns an element of ``self``, as per
                :meth:`Sets.ParentMethods.an_element`

                EXAMPLES::

                    sage: S = FiniteEnumeratedSets().IsomorphicObjects().example()
                    sage: S.an_element()  # indirect doctest
                    1
                """
                return self.retract(self.ambient().an_element())

    class Subobjects(SubobjectsCategory):
        """
        A category for subquotients of sets

        EXAMPLES::

            sage: Sets().Subobjects()
            Category of subobjects of sets
            sage: Sets().Subobjects().all_super_categories()
            [Category of subobjects of sets, Category of subquotients of sets, Category of sets,
             Category of sets with partial maps, Category of objects]
        """

        class ParentMethods:

            def _repr_(self):
                """
                EXAMPLES::

                    sage: from sage.categories.examples.semigroups import IncompleteSubquotientSemigroup
                    sage: S = IncompleteSubquotientSemigroup(category = Semigroups().Subobjects())
                    sage: S._repr_()
                    'A subobject of An example of a semigroup: the left zero semigroup'
                """
                return "A subobject of %s"%(self.ambient())

    class IsomorphicObjects(IsomorphicObjectsCategory):
        """
        A category for isomorphic objects of sets

        EXAMPLES::

            sage: Sets().IsomorphicObjects()
            Category of isomorphic objects of sets
            sage: Sets().IsomorphicObjects().all_super_categories()
            [Category of isomorphic objects of sets,
             Category of subobjects of sets, Category of quotients of sets,
             Category of subquotients of sets,
             Category of sets,
             Category of sets with partial maps,
             Category of objects]
        """

        class ParentMethods:

            def _repr_(self):
                """
                EXAMPLES::

                    sage: S = FiniteEnumeratedSets().IsomorphicObjects().example()
                    sage: S._repr_()
                    'The image by some isomorphism of An example of a finite enumerated set: {1,2,3}'
                """
                return "The image by some isomorphism of %s"%(self.ambient())

    class CartesianProducts(CartesianProductsCategory):
        """
        EXAMPLES::

            sage: C = Sets().CartesianProducts().example()
            sage: C
            The cartesian product of (Set of prime numbers (basic implementation), An example of an infinite enumerated set: the non negative integers, An example of a finite enumerated set: {1,2,3})
            sage: C.category()
            Category of Cartesian products of sets
            sage: C.categories()
            [Category of Cartesian products of sets, Category of sets,
             Category of sets with partial maps,
             Category of objects]
            sage: TestSuite(C).run()
        """

        def extra_super_categories(self):
            """
            A cartesian product of sets is a set

            EXAMPLES::

                sage: Sets().CartesianProducts().extra_super_categories()
                [Category of sets]
                sage: Sets().CartesianProducts().super_categories()
                [Category of sets]
            """
            return [Sets()]

        def example(self):
            """
            EXAMPLES::

                sage: Sets().CartesianProducts().example()
                The cartesian product of (Set of prime numbers (basic implementation), An example of an infinite enumerated set: the non negative integers, An example of a finite enumerated set: {1,2,3})
            """
            from finite_enumerated_sets import FiniteEnumeratedSets
            from infinite_enumerated_sets import InfiniteEnumeratedSets
            from cartesian_product import cartesian_product
            S1 = Sets().example()
            S2 = InfiniteEnumeratedSets().example()
            S3 = FiniteEnumeratedSets().example()
            return cartesian_product([S1, S2, S3])


        class ParentMethods:

            @cached_method
            def an_element(self):
                """
                EXAMPLES::

                    sage: C = Sets().CartesianProducts().example(); C
                    The cartesian product of (Set of prime numbers (basic implementation), An example of an infinite enumerated set: the non negative integers, An example of a finite enumerated set: {1,2,3})
                    sage: C.an_element()
                    (47, 42, 1)
                """
                return self._cartesian_product_of_elements(s.an_element() for s in self._sets)

            @abstract_method
            def _sets_keys(self):
                """
                Returns the indices of the summands of ``self``.

                EXAMPLES::

                    sage: cartesian_product([QQ, ZZ, ZZ])._sets_keys()
                    [0, 1, 2]
                """

            @abstract_method
            def summand_projection(self, i):
                """
                Returns the natural projection onto the `i`-th summand of self

                INPUTS:

                 - ``i`` -- the index of a summand of ``self``

                EXAMPLES::

                    sage: C = Sets().CartesianProducts().example(); C
                    The cartesian product of (Set of prime numbers (basic implementation), An example of an infinite enumerated set: the non negative integers, An example of a finite enumerated set: {1,2,3})
                    sage: x = C.an_element(); x
                    (47, 42, 1)
                    sage: pi = C.summand_projection(1)
                    sage: pi(x)
                    42
                """

            @abstract_method
            def _cartesian_product_of_elements(self, elements):
                """
                Returns the cartesian product of the given ``elements``

                INPUT:

                 - ``elements`` - a tuple with one element of each summand of ``self``

                EXAMPLES::

                    sage: S1 = Sets().example()
                    sage: S2 = InfiniteEnumeratedSets().example()
                    sage: C = cartesian_product([S2, S1, S2])
                    sage: C._cartesian_product_of_elements([S2.an_element(), S1.an_element(), S2.an_element()])
                    (42, 47, 42)
                """

        class ElementMethods:

            def summand_projection(self, i):
                """
                Returns the projection of ``self`` on the `i`-th
                summand of the cartesian product.

                INPUTS:

                 - ``i`` -- the index of a summand of the cartesian product

                EXAMPLES::

                    sage: F = CombinatorialFreeModule(ZZ, [4,5]); F.__custom_name = "F"
                    sage: G = CombinatorialFreeModule(ZZ, [4,6]); G.__custom_name = "G"
                    sage: S = cartesian_product([F, G])
                    sage: x = S.monomial((0,4)) + 2 * S.monomial((0,5)) + 3 * S.monomial((1,6))
                    sage: x.summand_projection(0)
                    B[4] + 2*B[5]
                    sage: x.summand_projection(1)
                    3*B[6]
                """
                return self.parent().summand_projection(i)(self)

            def summand_split(self):
                """
                Splits ``self`` into its summands

                EXAMPLES::

                    sage: F = CombinatorialFreeModule(ZZ, [4,5]); F.__custom_name = "F"
                    sage: G = CombinatorialFreeModule(ZZ, [4,6]); G.__custom_name = "G"
                    sage: H = CombinatorialFreeModule(ZZ, [4,7]); H.__custom_name = "H"
                    sage: S = cartesian_product([F, G, H])
                    sage: x = S.monomial((0,4)) + 2 * S.monomial((0,5)) + 3 * S.monomial((1,6)) + 4 * S.monomial((2,4)) + 5 * S.monomial((2,7))
                    sage: x.summand_split()
                    (B[4] + 2*B[5], 3*B[6], 4*B[4] + 5*B[7])
                    sage: [s.parent() for s in x.summand_split()]
                    [F, G, H]
                    sage: S.zero().summand_split()
                    (0, 0, 0)
                    sage: [s.parent() for s in S.zero().summand_split()]
                    [F, G, H]
                """
                # TODO: optimize
                return tuple(self.summand_projection(i) for i in self.parent()._sets_keys())
                #return Family(self._sets.keys(), self.projection)


    class Algebras(AlgebrasCategory):

        def extra_super_categories(self):
            """
            EXAMPLES::

                sage: Sets().Algebras(QQ).extra_super_categories()
                [Category of modules with basis over Rational Field]

                sage: Sets().Algebras(QQ).super_categories()
                [Category of modules with basis over Rational Field]

                sage: Sets().example().algebra(ZZ).categories()
                [Category of set algebras over Integer Ring,
                 Category of modules with basis over Integer Ring,
                 ...
                 Category of objects]

            """
            from sage.categories.modules_with_basis import ModulesWithBasis
            return [ModulesWithBasis(self.base_ring())]

    class WithRealizations(WithRealizationsCategory):

        def extra_super_categories(self):
            """
            A set with multiple realizations is a facade parent

            EXAMPLES::

                sage: Sets().WithRealizations().extra_super_categories()
                [Category of facade sets]
                sage: Sets().WithRealizations().super_categories()
                [Category of facade sets]
            """
            return [Sets().Facades()]

        def example(self, base_ring = None, set = None):
            r"""
            Returns an example of set with multiple realizations, as
            per :meth:`Category.example`.

            EXAMPLES::

                sage: Sets().WithRealizations().example()
                The subset algebra of {1, 2, 3} over Rational Field

                sage: Sets().WithRealizations().example(ZZ, Set([1,2]))
                The subset algebra of {1, 2} over Integer Ring
            """
            from sage.rings.rational_field import QQ
            from sage.sets.set import Set
            if base_ring is None:
                base_ring = QQ
            if set is None:
                set = Set([1,2,3])
            from sage.categories.examples.with_realizations import SubsetAlgebra
            return SubsetAlgebra(base_ring, set)


        class ParentMethods:

            def _test_with_realizations(self, **options):
                r"""
                Test that this parent with realizations is properly implemented

                INPUT::

                 - ``options`` -- any keyword arguments accepted by :meth:`_tester`.

                EXAMPLES::

                    sage: A = Sets().WithRealizations().example()
                    sage: A._test_with_realizations()

                See the documentation for :class:`TestSuite` for more information.
                """
                tester = self._tester(**options)
                for R in self.realizations():
                    tester.assert_(R in self.Realizations())
                # Could check that there are coerce maps between any two realizations

            @lazy_attribute
            def _realizations(self):
                """
                This lazily initializes the attribute
                ``_realizations`` the first time it is needed.

                TESTS::

                    sage: class MyParent(Parent):
                    ...      pass
                    sage: P = MyParent(category = Sets().WithRealizations())
                    sage: P._realizations
                    []
                """
                return []

            def _register_realization(self, realization):
                """
                EXAMPLES::

                    sage: A = Sets().WithRealizations().example(QQ[x]); A
                    The subset algebra of {1, 2, 3} over Univariate Polynomial Ring in x over Rational Field
                    sage: class ANewRealizationOfA(CombinatorialFreeModule):
                    ...       pass
                    sage: R = ANewRealizationOfA(A.base_ring(), A.F().basis().keys(), category = A.Realizations())
                    sage: R in A.realizations()  # indirect doctest
                    True

                Note: the test above uses ``QQ[x]`` to not interfer with other tests
                """
                assert realization.realization_of() is self
                self._realizations.append(realization)

            def inject_shorthands(self, verbose=True):
                """
                EXAMPLES::

                    sage: A = Sets().WithRealizations().example(QQ); A
                    The subset algebra of {1, 2, 3} over Rational Field
                    sage: A.inject_shorthands()
                    Injecting F as shorthand for The subset algebra of {1, 2, 3} over Rational Field in the Fundamental basis
                    Injecting In as shorthand for The subset algebra of {1, 2, 3} over Rational Field in the In basis
                    ...
                """
                from sage.misc.misc import inject_variable
                if not hasattr(self, "_shorthands"):
                    raise NotImplementedError("no shorthands defined for %s" % self)
                for shorthand in self._shorthands:
                    realization = getattr(self, shorthand)()
                    if verbose:
                        print 'Injecting %s as shorthand for %s'%(shorthand, realization)
                    inject_variable(shorthand, realization)

            def realizations(self):
                """
                Returns all the realizations of ``self`` that ``self`` is aware of.

                EXAMPLES::

                    sage: A = Sets().WithRealizations().example(); A
                    The subset algebra of {1, 2, 3} over Rational Field
                    sage: A.realizations()
                    [The subset algebra of {1, 2, 3} over Rational Field in the Fundamental basis, The subset algebra of {1, 2, 3} over Rational Field in the In basis, The subset algebra of {1, 2, 3} over Rational Field in the Out basis]

                .. note::

                    Constructing a parent ``P`` in the category
                    ``A.Realizations()`` automatically adds ``P`` to
                    this list by calling ``A._register_realization(A)``
                """
                return self._realizations

            def facade_for(self):
                """
                Returns the parents ``self`` is a facade for, that is
                the realizations of ``self``

                EXAMPLES::

                    sage: A = Sets().WithRealizations().example(); A
                    The subset algebra of {1, 2, 3} over Rational Field
                    sage: A.facade_for()
                    [The subset algebra of {1, 2, 3} over Rational Field in the Fundamental basis, The subset algebra of {1, 2, 3} over Rational Field in the In basis, The subset algebra of {1, 2, 3} over Rational Field in the Out basis]

                    sage: A = Sets().WithRealizations().example(); A
                    The subset algebra of {1, 2, 3} over Rational Field
                    sage: f = A.F().an_element(); f
                    F[{}] + 2*F[{1}] + 3*F[{2}] + F[{1, 2}]
                    sage: i = A.In().an_element(); i
                    In[{}] + 2*In[{1}] + 3*In[{2}] + In[{1, 2}]
                    sage: o = A.Out().an_element(); o
                    Out[{}] + 2*Out[{1}] + 3*Out[{2}] + Out[{1, 2}]
                    sage: f in A, i in A, o in A
                    (True, True, True)
                """
                return self.realizations()

            # Do we really want this feature?
            class Realizations(Category_realization_of_parent):

                def super_categories(self):
                    """
                    EXAMPLES::

                        sage: A = Sets().WithRealizations().example(); A
                        The subset algebra of {1, 2, 3} over Rational Field
                        sage: A.Realizations().super_categories()
                        [Category of realizations of sets]
                    """
                    return [Sets().Realizations()]

            def _an_element_(self):
                """
                Returns an element of some realization of ``self``

                EXAMPLES::

                    sage: A = Sets().WithRealizations().example(); A
                    The subset algebra of {1, 2, 3} over Rational Field
                    sage: A.an_element()        # indirect doctest
                    F[{}] + 2*F[{1}] + 3*F[{2}] + F[{1, 2}]
                """
                return self.realizations()[0].an_element()

            # TODO: maybe this could be taken care of by Sets.Facades()?
            def __contains__(self, x):
                r"""
                Test whether ``x`` is in ``self``, that is if it is an
                element of some realization of ``self``.

                EXAMPLES::

                    sage: A = Sets().WithRealizations().example(); A
                    The subset algebra of {1, 2, 3} over Rational Field
                    sage: A.an_element() in A
                    True
                    sage: A.In().an_element() in A
                    True
                    sage: A.F().an_element() in A
                    True
                    sage: A.Out().an_element() in A
                    True
                    sage: 1 in A
                    True
                    sage: QQ['x'].an_element() in A
                    False
                """
                return any(x in realization for realization in self.realizations())

    class Realizations(RealizationsCategory):

        class ParentMethods:

            def __init_extra__(self):
                """
                Registers ``self`` as a realization of ``self.realization_of``

                TESTS::

                    sage: A = Sets().WithRealizations().example()
                    sage: A.realizations()    # indirect doctest
                    [The subset algebra of {1, 2, 3} over Rational Field in the Fundamental basis,
                     The subset algebra of {1, 2, 3} over Rational Field in the In basis,
                     The subset algebra of {1, 2, 3} over Rational Field in the Out basis]
                """
                self.realization_of()._register_realization(self)

            @cached_method
            def realization_of(self):
                """
                Returns the parent this is a realization of

                EXAMPLES::

                    sage: A = Sets().WithRealizations().example(); A
                    The subset algebra of {1, 2, 3} over Rational Field
                    sage: In = A.In(); In
                    The subset algebra of {1, 2, 3} over Rational Field in the In basis
                    sage: In.realization_of()
                    The subset algebra of {1, 2, 3} over Rational Field
                """
                for category in self.categories():
                    if isinstance(category, Category_realization_of_parent):
                        return category.base()

            def _realization_name(self):
                """
                Returns the name of this realization

                In this default implementation, this is guessed from
                the name of its class.

                EXAMPLES::

                    sage: A = Sets().WithRealizations().example(); A
                    The subset algebra of {1, 2, 3} over Rational Field
                    sage: In = A.In(); In
                    The subset algebra of {1, 2, 3} over Rational Field in the In basis
                    sage: In._realization_name()
                    'In'
                """
                # The __base__ gets rid of the with_category
                # The split adds support for nested classes
                return self.__class__.__base__.__name__.split('.')[-1]

            def _repr_(self):
                """
                EXAMPLES::

                    sage: A = Sets().WithRealizations().example(); A
                    The subset algebra of {1, 2, 3} over Rational Field
                    sage: In = A.In(); In
                    The subset algebra of {1, 2, 3} over Rational Field in the In basis

                In the example above, :meth:`repr` was overriden by
                the category ``A.Realizations()``. We now add a new
                (fake) realization which is not in
                ``A.Realizations()`` to actually exercise this
                method::

                    sage: from sage.categories.realizations import Realizations
                    sage: class Blah(Parent):
                    ...       pass
                    sage: P = Blah(category = Sets.WithRealizations.ParentMethods.Realizations(A))
                    sage: P     # indirect doctest
                    The subset algebra of {1, 2, 3} over Rational Field in the realization Blah
                """
                return "%s in the realization %s"%(self.realization_of(), self._realization_name())
