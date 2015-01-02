r"""
Fixtures to help testing functionality

Utilities which modify or replace code to help with doctesting functionality.
Wrappers, proxies and mockups are typical examples of fixtures.

AUTHORS:

- Martin von Gagern (2014-12-15): PropertyAccessTracerProxy and trace_method

EXAMPLES:

You can use :func:`trace_method` to see how a method
communicates with its surroundings::

    sage: class Foo(object):
    ....:     def f(self):
    ....:         self.y = self.g(self.x)
    ....:     def g(self, arg):
    ....:         return arg + 1
    ....:
    sage: foo = Foo()
    sage: foo.x = 3
    sage: from sage.doctest.fixtures import trace_method
    sage: trace_method(foo, "f")
    sage: foo.f()
    enter f()
      read x = 3
      call g(3) -> 4
      write y = 4
    exit f -> None
"""

#*****************************************************************************
#       Copyright (C) 2014 Martin von Gagern <Martin.vGagern@gmx.net>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from functools import wraps

class PropertyAccessTracerHelper(object):

    def __init__(self, delegate, prefix="  ", reads=True):
        r"""
        Helper to print proxied access to properties.

        This class does the actual printing of access traces.
        The fact that it's not a proxy at the same time
        helps avoiding complicated attribute access syntax.

        INPUT:

        - ``delegate``: The actual object to be proxied.

        - ``prefix``: String to prepend to each printed output.

        - ``reads``: Whether to trace read access as well.

        EXAMPLE::

            sage: class Foo(object):
            ....:     def f(self, *args):
            ....:         return self.x*self.x
            ....:
            sage: foo = Foo()
            sage: from sage.doctest.fixtures import PropertyAccessTracerHelper
            sage: pat = PropertyAccessTracerHelper(foo)
            sage: pat.set("x", 2)
              write x = 2
            sage: pat.get("x")
              read x = 2
            2
            sage: pat.get("f")(3)
              call f(3) -> 4
            4
        """
        self.delegate = delegate
        self.prefix = prefix
        self.reads = reads

    def get(self, name):
        r"""
        EXAMPLE::

            sage: class Foo(object):
            ....:     def f(self, *args):
            ....:         return self.x*self.x
            ....:
            sage: foo = Foo()
            sage: foo.x = 2
            sage: from sage.doctest.fixtures import PropertyAccessTracerHelper
            sage: pat = PropertyAccessTracerHelper(foo)
            sage: pat.get("x")
              read x = 2
            2
            sage: pat.get("f")(3)
              call f(3) -> 4
            4
        """
        val = getattr(self.delegate, name)
        if callable(val) and name not in self.delegate.__dict__:
            @wraps(val)
            def wrapper(*args, **kwds):
                arglst = [self.fmt(arg) for arg in args]
                arglst.extend("{}={}".format(k, self.fmt(v))
                              for k, v in sorted(kwds.items()))
                res = val(*args, **kwds)
                print("{}call {}({}) -> {}"
                      .format(self.prefix, name, ", ".join(arglst),
                              self.fmt(res)))
                return res
            return wrapper
        else:
            if self.reads:
                print("{}read {} = {}".format(self.prefix, name, self.fmt(val)))
            return val

    def set(self, name, val):
        r"""
        EXAMPLE::

            sage: class Foo(object):
            ....:     pass
            ....:
            sage: foo = Foo()
            sage: from sage.doctest.fixtures import PropertyAccessTracerHelper
            sage: pat = PropertyAccessTracerHelper(foo)
            sage: pat.set("x", 2)
              write x = 2
            sage: foo.x
            2
        """
        print("{}write {} = {}".format(self.prefix, name, self.fmt(val)))
        setattr(self.delegate, name, val)

    @classmethod
    def fmt(cls, val):
        r"""
        Format a value to be printed.

        This can be used to introduce normalization,
        such that the printed value does not depend on factors
        outside the control of the doctest.
        One example is the order of elements in a hash-based structure.
        For most objects, this is simply the ``repr`` of the object.

        EXAMPLE::

            sage: from sage.doctest.fixtures import PropertyAccessTracerHelper
            sage: fmt = PropertyAccessTracerHelper.fmt
            sage: print(fmt(set(["a", "c", "b", "d"])))
            set(['a', 'b', 'c', 'd'])
            sage: print(fmt(frozenset(["a", "c", "b", "d"])))
            frozenset(['a', 'b', 'c', 'd'])
            sage: print(fmt("foo\nbar"))
            'foo\nbar'
        """
        if isinstance(val, frozenset):
            return ("frozenset([{}])".format
                    (", ".join(map(cls.fmt, sorted(val)))))
        if isinstance(val, set):
            return ("set([{}])".format
                    (", ".join(map(cls.fmt, sorted(val)))))
        r = repr(val)
        return r


class PropertyAccessTracerProxy(object):

    def __init__(self, delegate, **kwds):
        r"""
        Proxy object which prints all property and method access to an object.

        The implementation is kept lean since all access to properties of
        the proxy itself requires complicated syntax.
        For this reason, the actual handling of property access
        is delegated to a :class:`PropertyAccessTracerHelper`.

        INPUT:

        - ``delegate``: The actual object to be proxied.

        - ``prefix``: String to prepend to each printed output.
          (Default: ``"  "``)

        - ``reads``: Whether to trace read access as well.
          (Default: ``True)

        EXAMPLE::

            sage: class Foo(object):
            ....:     def f(self, *args):
            ....:         return self.x*self.x
            ....:
            sage: foo = Foo()
            sage: from sage.doctest.fixtures import PropertyAccessTracerProxy
            sage: pat = PropertyAccessTracerProxy(foo)
            sage: pat.x = 2
              write x = 2
            sage: pat.x
              read x = 2
            2
            sage: pat.f(3)
              call f(3) -> 4
            4
        """
        helper = PropertyAccessTracerHelper(delegate, **kwds)
        object.__setattr__(self, "helper", helper)

    def __getattribute__(self, name):
        r"""
        EXAMPLE::

            sage: class Foo(object):
            ....:     def f(self, *args):
            ....:         return self.x*self.x
            ....:
            sage: foo = Foo()
            sage: foo.x = 2
            sage: from sage.doctest.fixtures import PropertyAccessTracerProxy
            sage: pat = PropertyAccessTracerProxy(foo)
            sage: pat.x
              read x = 2
            2
            sage: pat.f(3)
              call f(3) -> 4
            4
        """
        helper = object.__getattribute__(self, "helper")
        return helper.get(name)

    def __setattr__(self, name, val):
        r"""
        EXAMPLE::

            sage: class Foo(object):
            ....:     pass
            ....:
            sage: foo = Foo()
            sage: from sage.doctest.fixtures import PropertyAccessTracerProxy
            sage: pat = PropertyAccessTracerProxy(foo)
            sage: pat.x = 2
              write x = 2
            sage: foo.x
            2
        """
        helper = object.__getattribute__(self, "helper")
        return helper.set(name, val)


def trace_method(obj, meth, **kwds):
    r"""
    Trace the doings of a given method.
    It prints method entry with arguments,
    access to members and other methods during method execution
    as well as method exit with return value.

    INPUT:

    - ``obj``: The object containing the method.

    - ``meth``: The name of the method to be traced.

    - ``prefix``: String to prepend to each printed output.
      (Default: ``"  "``)

    - ``reads``: Whether to trace read access as well.
      (Default: ``True)

    EXAMPLE::

        sage: class Foo(object):
        ....:     def f(self):
        ....:         self.y = self.g(self.x)
        ....:     def g(self, arg):
        ....:         return arg + 1
        ....:
        sage: foo = Foo()
        sage: foo.x = 3
        sage: from sage.doctest.fixtures import trace_method
        sage: trace_method(foo, "f")
        sage: foo.f()
        enter f()
          read x = 3
          call g(3) -> 4
          write y = 4
        exit f -> None
    """
    f = getattr(obj, meth).__func__
    t = PropertyAccessTracerProxy(obj, **kwds)
    fmt = PropertyAccessTracerHelper.fmt
    @wraps(f)
    def g(*args, **kwds):
        arglst = [fmt(arg) for arg in args]
        arglst.extend("{}={}".format(k, fmt(v))
                      for k, v in sorted(kwds.items()))
        print("enter {}({})".format(meth, ", ".join(arglst)))
        res = f(t, *args, **kwds)
        print("exit {} -> {}".format(meth, fmt(res)))
        return res
    setattr(obj, meth, g)
