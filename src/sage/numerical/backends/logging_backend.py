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

#from sage.numerical.backends.generic_backend import GenericBackend

def format_function_call(fn_name, *v, **k):
    args = [ repr(a) for a in v ] + [ "%s=%r" % (arg,val) for arg, val in k.items() ]
    return "{}({})".format(fn_name, ", ".join(args))

class LoggingBackend:

    """
    EXAMPLES::

        sage: import sage.numerical.backends.logging_backend
        sage: from sage.numerical.backends.logging_backend import LoggingBackend
        sage: from sage.numerical.backends.generic_backend import get_solver
        sage: b = get_solver(solver = "GLPK")
        sage: lb = LoggingBackend(backend=b)
        sage: lb.add_variable(obj=42, name='Helloooooo')
        # b.add_variable(obj=42, name='Helloooooo')
        # result: 0
        0
        sage: lb.add_variable(obj=1789)
        # b.add_variable(obj=1789)
        # result: 1
        1
    """
    
    def __init__(self, backend, printing=True, doctest=None, test_method=None):
        self._backend = backend
        self._printing = printing
        self._doctest = doctest
        self._test_method = test_method

    def __getattr__(self, attr):
        a = getattr(self._backend, attr)
        if callable(a):
            # a method
            def m(backend, *args, **kwdargs):
                if self._printing:
                    print "# {}".format(format_function_call("b." + attr, *args, **kwdargs))
                if self._doctest:
                    self._doctest.write("        sage: {}\n".format(format_function_call("b." + attr, *args, **kwdargs)))
                result = a(*args, **kwdargs)
                if self._printing:
                    print "# result: {}".format(result)
                if self._doctest:
                    self._doctest.write("        {}".format(result))
                return result
            import types
            mm = types.MethodType(m,self)
            # cache it
            setattr(self, attr, mm)
            return mm
        else:
            return a

### does not work because get_solver needs to return a subclass of GenericBackend.
### but if we make LoggingBackend a subclass of GenericBackend, then __getattr__ does not get run for all those methods.
### Need to override methods by iterating through the dictionary, I guess...

def LoggingBackendFactory(solver=None, printing=True, doctest_file=None, test_method_file=None):

    """
    The result of this can be passed as an argument to `get_solver`.

    EXAMPLES:

    In-terminal logging::

        sage: import sage.numerical.backends.logging_backend
        sage: from sage.numerical.backends.logging_backend import LoggingBackendFactory
        sage: from sage.numerical.backends.generic_backend import get_solver
        sage: lb = get_solver(solver = LoggingBackendFactory(solver='GLPK'))
        sage: lb.add_variable(obj=42, name='Helloooooo')
        # b.add_variable(obj=42, name='Helloooooo')
        # result: 0
        0
        sage: lb.add_variable(obj=1789)
        # b.add_variable(obj=1789)
        # result: 1
        1

    Write doctests to file::

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
        from sage.numerical.backends.generic_backend import get_solver
        return LoggingBackend(backend=get_solver(solver),
                              printing=printing, doctest=doctest, test_method=test_method)

    return logging_solver
