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

        sage: TestSuite(1).run(verbose = True)
        running ._test_category() . . . pass
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
        running ._test_category() . . . pass
        running ._test_elements() . . .
          Running the test suite of self.an_element()
          running ._test_category() . . . pass
          running ._test_not_implemented_methods() . . . pass
          running ._test_pickling() . . . pass
          pass
        running ._test_enumerated_set_contains() . . . pass
        running ._test_enumerated_set_iter_cardinality() . . . pass
        running ._test_enumerated_set_iter_list() . . . pass
        running ._test_not_implemented_methods() . . . pass
        running ._test_pickling() . . . pass
        running ._test_some_elements() . . . pass

    The different test methods can be called independently::

        sage: S._test_associativity()

    Debugging tip: in case of failure of some test, use `%pdb on` to
    turn on automatic debugging on error. Run the failing test
    independtly: the debugger will stop right where the first
    assertion fails. Then, introspection can be used to analyse what
    exactly the problem is. See also the ``catch = False`` option to
    :meth:`.run`.

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

    Finally, running ``TestSuite`` on a standard Python object does
    some basic sanity checks:

        sage: TestSuite(int(1)).run(verbose = True)
        running ._test_pickling() . . . pass

    TODO:

     - Allow for customized behavior in case of failing assertion
       (warning, error, statistic accounting).
       This involves reimplementing the methods fail / failIf / ...
       of unittest.TestCase in InstanceTester

     - Don't catch the exceptions if ``TestSuite(..).run()`` is called
       under the debugger, or with ``%pdb`` on (how to detect this? see
       ``IPython.ipapi.get()``, ``IPython.Magic.shell.call_pdb``, ...)
       In the mean time, see the ``catch=False`` option.

     - Run the tests according to the inheritance order, from most
       generic to most specific, rather than alphabetically. Then, the
       first failure will be the most relevant, the others being
       usually consequences.

     - Improve integration with doctests (statistics on failing/passing tests)

     - Add proper support for nested testsuites.

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
        from sage.structure.sage_object import SageObject
        if not isinstance(instance, (SageObject,PythonObjectWithTests)):
            instance = PythonObjectWithTests(instance)
        self._instance = instance

    def __repr__(self):
        """
        TESTS::

            sage: TestSuite(ZZ)
            Test suite for Integer Ring
        """
        return "Test suite for %s"%self._instance


    def run(self, category = None, skip = [], catch = True, raise_on_failure = False, **options):
        """
        Run all the tests from this test suite:

        INPUT:

         - ``category``         - a category; reserved for future use
         - ``skip``             - a string or list (or iterable) of strings
         - ``raise_on_failure`` - a boolean (default: False)
         - ``catch``            - a boolean (default: True)

        All other options are passed down to the individual tests.

        EXAMPLES::

            sage: TestSuite(ZZ).run()

        We now use the ``verbose`` option::

            sage: TestSuite(1).run(verbose = True)
            running ._test_category() . . . pass
            running ._test_not_implemented_methods() . . . pass
            running ._test_pickling() . . . pass

        Some tests may be skipped using the ``skip`` option::

            sage: TestSuite(1).run(verbose = True, skip ="_test_pickling")
            running ._test_category() . . . pass
            running ._test_not_implemented_methods() . . . pass
            sage: TestSuite(1).run(verbose = True, skip =["_test_pickling", "_test_category"])
            running ._test_not_implemented_methods() . . . pass

        We now show (and test) some standard error reports::

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
            running ._test_category() . . . pass
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

        The ``catch=False`` option prevents ``TestSuite`` from
        catching exceptions::

            sage: TestSuite(Blah()).run(catch=False)
            Traceback (most recent call last):
              ...
              File ..., in _test_b
                def _test_b(self, tester): tester.fail()
              ...
            AssertionError

        In conjonction with ``%pdb on``, this allows for the debbuger
        to jump directly to the first failure location.
        """
        if type(skip) == str:
            skip = [skip]
        else:
            skip = tuple(skip)

        # The class of exceptions that will be catched and reported;
        # other exceptions will get trough. None catches nothing.
        catch_exception = Exception if catch else None

        tester = instance_tester(self._instance, **options)
        failed = []
        for method_name in dir(self._instance):
            if method_name[0:6] == "_test_" and method_name not in skip:
                # TODO: improve pretty printing
                # could use the doc string of the test method?
                tester.info(tester._prefix+"running .%s() . . . "%method_name, newline = False)
                test_method = getattr(self._instance, method_name)
                try:
                    test_method(tester = tester)
                    tester.info("pass")
                except catch_exception as e:
                    failed.append(method_name)
                    if isinstance(e, TestSuiteFailure):
                        # The failure occured in a nested testsuite
                        # which has already reported the details of
                        # that failure
                        if not tester._verbose:
                            print tester._prefix+"Failure in %s"%method_name
                    else:
                        if tester._verbose:
                            tester.info("fail")
                        else:
                            print tester._prefix+"Failure in %s:"%method_name
                        s = traceback.format_exc()
                        print tester._prefix + s.strip().replace("\n", "\n"+tester._prefix)
                        print tester._prefix + "-" * 60
        if len(failed) > 0:
            print tester._prefix+"The following tests failed: %s"%(", ".join(failed))
            if raise_on_failure:
                raise TestSuiteFailure

class TestSuiteFailure(AssertionError):
    pass

def instance_tester(instance, tester = None, **options):
    """
    Returns a gadget attached to ``instance`` providing testing utilities.

    EXAMPLES::

        sage: from sage.misc.sage_unittest import instance_tester
        sage: tester = instance_tester(ZZ)

        sage: tester.assert_(1 == 1)
        sage: tester.assert_(1 == 0)
        Traceback (most recent call last):
        ...
        AssertionError
        sage: tester.assert_(1 == 0, "this is expected to fail")
        Traceback (most recent call last):
        ...
        AssertionError: this is expected to fail

        sage: tester.assertEquals(1, 1)
        sage: tester.assertEquals(1, 0)
        Traceback (most recent call last):
        ...
        AssertionError: 1 != 0

    The available assertion testing facilities are the same as in
    :class:`unittest.TestCase`, which see (actually, by a slight
    abuse, tester is currently an instance of this class).

    TESTS::

        sage: instance_tester(ZZ, tester = tester) is tester
        True
    """
    if tester is None:
        return InstanceTester(instance, **options)
    else:
        assert len(options) == 0
        assert tester._instance is instance
        return tester

class InstanceTester(unittest.TestCase):
    """
    A gadget attached to an instance providing it with testing utilities.

    EXAMPLES::

        sage: from sage.misc.sage_unittest import InstanceTester
        sage: InstanceTester(instance = ZZ, verbose = True, elements = [1,2,3])
        Testing utilities for Integer Ring

    This is used by ``SageObject._tester``, which see::

        sage: QQ._tester()
        Testing utilities for Rational Field
    """

    def __init__(self, instance, elements = None, verbose = False, prefix = "", **options):
        """
        A gadget attached to an instance providing it with testing utilities.

        EXAMPLES::

            sage: from sage.misc.sage_unittest import InstanceTester
            sage: InstanceTester(instance = ZZ, verbose = True, elements = [1,2,3])
            Testing utilities for Integer Ring

        This is used by ``SageObject._tester``, which see::

            sage: QQ._tester()
            Testing utilities for Rational Field
        """
        self._instance = instance
        self._verbose = verbose
        self._elements = elements
        self._prefix = prefix

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

class PythonObjectWithTests(object):
    """
    Utility class for running basis tests on a plain Python object
    (that is not in SageObject). More test methods can be added here.

    EXAMPLES::

            sage: TestSuite("bla").run()

    """
    def __init__(self, instance):
        """
        EXAMPLES::

            sage: from sage.misc.sage_unittest import PythonObjectWithTests
            sage: x = PythonObjectWithTests(int(1)); x
            <sage.misc.sage_unittest.PythonObjectWithTests object at ...>
            sage: TestSuite(x).run()
        """

        self._instance = instance

    def _test_pickling(self, **options):
        """
        Checks that the instance in self can be pickled and unpickled properly.

        EXAMPLES::

            sage: from sage.misc.sage_unittest import PythonObjectWithTests
            sage: PythonObjectWithTests(int(1))._test_pickling()

        SEE ALSO: :func:`dumps` :func:`loads`
        """
        tester = instance_tester(self, **options)
        from sage.misc.all import loads, dumps
        tester.assertEqual(loads(dumps(self._instance)), self._instance)
