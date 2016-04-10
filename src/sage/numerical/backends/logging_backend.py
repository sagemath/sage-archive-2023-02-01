r"""
Logging Backend

It records all calls to backend methods.
"""

#*****************************************************************************
#       Copyright (C) 2016 Matthias Koeppe <mkoeppe@math.ucdavis.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.numerical.backends.generic_backend import GenericBackend

def format_function_call(fn_name, *v, **k):
    args = [ repr(a) for a in v ] + [ "%s=%r" % (arg,val) for arg, val in k.items() ]
    return "{}({})".format(fn_name, ", ".join(args))

def _make_wrapper(attr):
    def m(self, *args, **kwdargs):
        funcall = format_function_call("p." + attr, *args, **kwdargs)
        a = getattr(self._backend, attr)
        if self._printing:
            print "# {}".format(funcall)
        if self._doctest:
            self._doctest.write("        sage: {}\n".format(funcall))
        try:
            result = a(*args, **kwdargs)
        except Exception as e:
            if self._printing:
                print "# exception: {}".format(e)
            if self._doctest:
                self._doctest.write("        {}\n".format(e)) # TODO: print "Traceback" etc.
            if self._test_method:
                self._test_method.write(("        with tester.assertRaises({}) as cm:\n"+
                                         "            {}\n").format(type(e).__name__, funcall))
            raise
        else:
            if self._printing:
                print "# result: {}".format(result)
            if self._doctest:
                self._doctest.write("        {}\n".format(result))
            if self._test_method:
                if result is None:
                    self._test_method.write("        tester.assertIsNone({})\n".format(funcall))
                elif type(result) is float:
                    # TODO: by default assertAlmostEqual does 7 decimal places (not significant digits)
                    # better perhaps to compute an appropriate 'places' or 'delta' parameter from result.
                    self._test_method.write("        tester.assertAlmostEqual({}, {})\n".format(funcall, result))
                else:
                    self._test_method.write("        tester.assertEqual({}, {})\n".format(funcall, result))
        return result
    return m

class LoggingBackend (GenericBackend):

    """
    EXAMPLES::

        sage: import sage.numerical.backends.logging_backend
        sage: from sage.numerical.backends.logging_backend import LoggingBackend
        sage: from sage.numerical.backends.generic_backend import get_solver
        sage: b = get_solver(solver = "GLPK")
        sage: lb = LoggingBackend(backend=b)
        sage: lb.add_variable(obj=42, name='Helloooooo')
        # p.add_variable(obj=42, name='Helloooooo')
        # result: 0
        0
        sage: lb.add_variable(obj=1789)
        # p.add_variable(obj=1789)
        # result: 1
        1
    """

    def __init__(self, backend, printing=True, doctest=None, test_method=None):
        self._backend = backend
        self._printing = printing
        self._doctest = doctest
        self._test_method = test_method

    # This getattr is there to create delegating method for all methods
    # that are not part of the GenericBackend interface
    def __getattr__(self, attr):
        a = getattr(self._backend, attr)
        if callable(a):
            # make a bound method
            import types
            mm = types.MethodType(_make_wrapper(attr), self)
            # cache it
            setattr(self, attr, mm)
            return mm
        else:
            return a

# Override all methods that we inherited from GenericBackend
# by delegating methods
for attr in dir(LoggingBackend):
    if not attr.startswith("_") and attr not in ("zero", "base_ring"):
        a = getattr(LoggingBackend, attr)
        if callable(a):
            # make an unbound method
            import types
            mm = types.MethodType(_make_wrapper(attr), None, LoggingBackend)
            setattr(LoggingBackend, attr, mm)

test_method_template = \
r'''
    @classmethod
    def _test_{name}(cls, tester=None, **options):
        """
        Run tests on ...

        TEST::

            SAGE: from sage.numerical.backends.generic_backend import GenericBackend
            SAGE: p = GenericBackend()
            SAGE: p._test_{name}()
            Traceback (most recent call last):
            ...
            NotImplementedError

        """
        p = cls()                         # fresh instance of the backend
        if tester is None:
            tester = p._tester(**options)
'''.replace("SAGE:", "sage:") # so that the above test does not get picked up by the doctester

def LoggingBackendFactory(solver=None, printing=True, doctest_file=None, test_method_file=None, method_name='CHANGE'):

    """
    The result of this can be passed as an argument to `get_solver`.

    EXAMPLES:

    In-terminal logging::

        sage: import sage.numerical.backends.logging_backend
        sage: from sage.numerical.backends.logging_backend import LoggingBackendFactory
        sage: from sage.numerical.backends.generic_backend import get_solver
        sage: lb = get_solver(solver = LoggingBackendFactory(solver='GLPK'))
        sage: lb.add_variable(obj=42, name='Helloooooo')
        # p.add_variable(obj=42, name='Helloooooo')
        # result: 0
        0
        sage: lb.add_variable(obj=1789)
        # p.add_variable(obj=1789)
        # result: 1
        1

    Write doctests to file::

        sage: fname = tmp_filename()
        sage: logging_solver = LoggingBackendFactory(solver='GLPK', printing=False, doctest_file=fname)
        sage: lb = get_solver(solver = logging_solver)
        sage: lb.add_variable(obj=42, name='Helloooooo')
        0
        sage: lb.add_variable(obj=1789)
        1
        sage: with open(fname) as f:
        ....:     for line in f.readlines(): print '|{}'.format(line),
        |        sage: p.add_variable(obj=42, name='Helloooooo')
        |        0
        |        sage: p.add_variable(obj=1789)
        |        1

    Write test method to file::

        sage: fname = tmp_filename()
        sage: logging_solver = LoggingBackendFactory(solver='GLPK', printing=False, test_method_file=fname, method_name='something')
        sage: lb = get_solver(solver = logging_solver)
        sage: lb.add_variable(obj=42, name='Helloooooo')
        0
        sage: lb.add_variable(obj=1789)
        1
        sage: lb.solve()
        Traceback (most recent call last):
        ...
        MIPSolverException...
        sage: with open(fname) as f:
        ....:     for line in f.readlines(): print '|{}'.format(line),

   """

    if doctest_file is not None:
        doctest = open(doctest_file, "w", 1) #line-buffered
    else:
        doctest = None
    if test_method_file is not None:
        test_method = open(test_method_file, "w", 1) #line-buffered
    else:
        test_method = None

    def logging_solver(**kwds):
        if test_method is not None:
            test_method.write(test_method_template.format(name=method_name))
        from sage.numerical.backends.generic_backend import get_solver
        return LoggingBackend(backend=get_solver(solver),
                              printing=printing, doctest=doctest, test_method=test_method)

    return logging_solver
