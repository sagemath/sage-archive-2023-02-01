r"""
Unit testing for Sage objects
"""
#*****************************************************************************
#  Copyright (C) 2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

##############################################################################

import unittest
import sys
import traceback

class TestSuite(object):
    """
    Test suites for Sage objects.

    EXAMPLES::

        sage: TestSuite(ZZ).run()

    No output means that all tests passed. Which tests?
    In practice this calls all the methods ``._test_*`` of this
    object, in alphabetic order::

        sage: TestSuite(ZZ).run(verbose = True)
        running ._test_not_implemented_methods() . . . pass
        running ._test_pickling() . . . pass

    Those methods are typically implemented by abstract
    super classes, in particular via categories, in order to
    enforce standard behavior and API, or provide mathematical
    sanity checks. For example if ``self`` is in the category of
    finite semigroups, this checks that the multiplication is
    associative (at least on some elements)::

        sage: S = FiniteSemigroups().example(alphabet = ('a', 'b'))
        sage: TestSuite(S).run(verbose = True)
        running ._test_an_element() . . . pass
        running ._test_associativity() . . . pass
        running ._test_element_pickling() . . . pass
        running ._test_enumerated_set_contains() . . . pass
        running ._test_enumerated_set_iter_cardinality() . . . pass
        running ._test_enumerated_set_iter_list() . . . pass
        running ._test_not_implemented_methods() . . . pass
        running ._test_pickling() . . . pass
        running ._test_some_elements() . . . pass

    The different test methods can be called independently::

        sage: S._test_associativity()

    When meaningful, one can further customize on which elements
    the tests are run. Here, we use it to *prove* that the
    multiplication is indeed associative, by running the test on
    all the elements::

        sage: S._test_associativity(elements = S)

    Adding a new test boils down to adding a new method in the class
    of the object or any super class (e.g. in a category). This method
    should use the utility :meth:`._tester` to handle standard options
    and report test failures. See the code of
    :meth:`._test_an_element` for an example. Note: Python's testunit
    convention is to look for methods called ``.test*``; we use instead
    ``._test_*`` so as not to pollute the object's interface.

    Eventually, every implementation of a :class:`SageObject` should
    run a :class:`TestSuite` on one of its instances in its doctest
    (replacing the current ``loads(dumps(x))`` tests).

    TODO:

     - allow for customized behavior in case of failing assertion
       (warning, error, statistic accounting).
       This involves reimplementing the methods fail / failIf / ...
       of unittest.TestCase in InstanceTester

     - Don't catch the exceptions if ``TestSuite(..).run()`` is called
       under the debugger, or with ``%pdb`` on (how to detect this? see
       ``IPython.ipapi.get()``, ``IPython.Magic.shell.call_pdb``, ...)

     - Run the tests according to the inheritance order, from most
       generic to most specific, rather than alphabetically. Then, the
       first failure will be the most relevant, the others being
       usually consequences.

     - Improve integration with doctests (statistics on failing/passing tests)

     - Integration with unittest:
       Make TestSuite inherit from unittest.TestSuite?
       Make ``.run(...)`` accept a result object

     - Add some standard option ``proof = True``, asking for the
       test method to choose appropriately the elements so as to
       prove the desired property. The test method may assume that
       a parent implements properly all the super categories. For
       example, the ``_test_commutative`` method of the category
       ``CommutativeSemigroups()`` may just check that the
       provided generators commute, implicitly assuming that
       generators indeed generate the semigroup (as required by
       ``Semigroups()``).
    """

    def __init__(self, instance):
        """
        TESTS::

            sage: TestSuite(ZZ)
            Test suite for Integer Ring
        """
        self._instance = instance

    def __repr__(self):
        """
        TESTS::

            sage: TestSuite(ZZ)
            Test suite for Integer Ring
        """
        return "Test suite for %s"%self._instance


    def run(self, category = None, skip = [], **options):
        """
        Run all the tests from this test suite:

        INPUT:

         - ``category`` - a category; reserved for future use
         - ``skip` ` - a string or list (or iterable) of strings

        All other options are passed down to the individual tests.

        EXAMPLES::

            sage: TestSuite(ZZ).run()

        We now use the ``verbose`` option::

            sage: TestSuite(ZZ).run(verbose = True)
            running ._test_not_implemented_methods() . . . pass
            running ._test_pickling() . . . pass

        Some tests may be skipped using the ``skip`` option::

            sage: TestSuite(ZZ).run(verbose = True, skip ="_test_pickling")
            running ._test_not_implemented_methods() . . . pass
            sage: TestSuite(ZZ).run(verbose = True, skip =["_test_pickling"])
            running ._test_not_implemented_methods() . . . pass

            sage: class Blah(SageObject):
            ...       def _test_a(self, tester): pass
            ...       def _test_b(self, tester): tester.fail()
            ...       def _test_c(self, tester): pass
            ...       def _test_d(self, tester): tester.fail()

            sage: TestSuite(Blah()).run()
            Failure in _test_b:
            Traceback (most recent call last):
              ...
            AssertionError
            ------------------------------------------------------------
            Failure in _test_d:
            Traceback (most recent call last):
              ...
            AssertionError
            ------------------------------------------------------------
            Failure in _test_pickling:
            Traceback (most recent call last):
              ...
            PicklingError: Can't pickle <class '__main__.Blah'>: attribute lookup __main__.Blah failed
            ------------------------------------------------------------
            The following tests failed: _test_b, _test_d, _test_pickling

            sage: TestSuite(Blah()).run(verbose = True)
            running ._test_a() . . . pass
            running ._test_b() . . . fail
            Traceback (most recent call last):
              ...
            AssertionError
            ------------------------------------------------------------
            running ._test_c() . . . pass
            running ._test_d() . . . fail
            Traceback (most recent call last):
              ...
            AssertionError
            ------------------------------------------------------------
            running ._test_not_implemented_methods() . . . pass
            running ._test_pickling() . . . fail
            Traceback (most recent call last):
              ...
            PicklingError: Can't pickle <class '__main__.Blah'>: attribute lookup __main__.Blah failed
            ------------------------------------------------------------
            The following tests failed: _test_b, _test_d, _test_pickling

            File "/opt/sage/local/lib/python/site-packages/sage/misc/sage_unittest.py", line 183, in run
            test_method(tester = tester)

        """
        if type(skip) == str:
            skip = [skip]
        else:
            skip = tuple(skip)
        tester = self._instance._tester(**options)
        failed = []
        for method_name in dir(self._instance):
            if method_name[0:6] == "_test_" and method_name not in skip:
                # TODO: improve pretty printing
                # could use the doc string of the test method?
                tester.info("running .%s() . . . "%method_name, newline = False)
                test_method = getattr(self._instance, method_name)
                try:
                    test_method(tester = tester)
                    tester.info("pass")
                except:
                    if tester._verbose:
                        tester.info("fail")
                    else:
                        print "Failure in %s:"%method_name
                    traceback.print_exc(file = sys.stdout)
                    print "-" * 60
                    failed.append(method_name)
        if len(failed) > 0:
            print "The following tests failed: %s"%(", ".join(failed))

class InstanceTester(unittest.TestCase):
    """
    A gadget attached to an instance providing it with testing utilities.

    EXAMPLES::

        sage: from sage.misc.sage_unittest import InstanceTester
        sage: InstanceTester(instance = ZZ, verbose = True, elements = [1,2,3])
        Testing utilities for Integer Ring

    This is used by ``SageObject._tester``, which see::

        sage: ZZ._tester()
        Testing utilities for Integer Ring
    """

    def __init__(self, instance, elements = None, verbose = False, **options):
        """
        A gadget attached to an instance providing it with testing utilities.

        EXAMPLES::

            sage: from sage.misc.sage_unittest import InstanceTester
            sage: InstanceTester(instance = ZZ, verbose = True, elements = [1,2,3])
            Testing utilities for Integer Ring

        This is used by ``SageObject._tester``, which see::

            sage: ZZ._tester()
            Testing utilities for Integer Ring
        """
        self._instance = instance
        self._verbose = verbose
        self._elements = elements

    def runTest(self):
        """
        Trivial implementation of :meth:`unittest.TestCase.runTest` to
        please the super class :class:`TestCase`. That's the price to
        pay for abusively inheriting from it.

            sage: from sage.misc.sage_unittest import InstanceTester
            sage: tester = InstanceTester(ZZ, verbose = True)
            sage: tester.runTest()
        """
        pass

    def info(self, message, newline = True):
        """
        Displays user information

        EXAMPLES::

            sage: from sage.misc.sage_unittest import InstanceTester
            sage: tester = InstanceTester(ZZ, verbose = True)

            sage: tester.info("hello"); tester.info("world")
            hello
            world

            sage: tester = InstanceTester(ZZ, verbose = False)
            sage: tester.info("hello"); tester.info("world")

            sage: tester = InstanceTester(ZZ, verbose = True)
            sage: tester.info("hello", newline = False); tester.info(" world")
            hello world
        """
        if self._verbose:
            if newline:
                sys.stdout.write(message+"\n")
            else:
                sys.stdout.write(message)
                sys.stdout.flush()

    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.misc.sage_unittest import InstanceTester
            sage: InstanceTester(ZZ, verbose = True)
            Testing utilities for Integer Ring

        """
        return "Testing utilities for %s"%self._instance


    def some_elements(self):
        """
        Returns a list (or iterable) of elements of ``self`` on which
        the tests should be run. This is only meaningful for container
        objects like parents.

        By default, this calls :meth:`.some_elements`::

            sage: from sage.misc.sage_unittest import InstanceTester
            sage: class MyParent(Parent):
            ...       def some_elements(self):
            ...           return [1,2,3,4,5]
            ...
            sage: tester = InstanceTester(MyParent())
            sage: list(tester.some_elements())
            [1, 2, 3, 4, 5]

            sage: tester = InstanceTester(ZZ, elements = [1,3,5])
            sage: list(tester.some_elements())
            [1, 3, 5]

        """
        if self._elements is None:
            return self._instance.some_elements()
        else:
            return self._elements
