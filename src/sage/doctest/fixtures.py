from functools import wraps

class PropertyAccessTracerProxy(object):
    r"""
    Proxy object which prints all property and method access to an object.

    INPUT:

    -  ``delegate``: The actual object to be proxied.

    -  ``prefix``: String to prepend to each printed output.

    -  ``reads``: Whether to trace read access as well

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

    def __init__(self, delegate, prefix="  ", reads=True):
        object.__setattr__(self, "delegate", delegate)
        object.__setattr__(self, "prefix", prefix)
        object.__setattr__(self, "reads", reads)

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
        delegate = object.__getattribute__(self, "delegate")
        prefix = object.__getattribute__(self, "prefix")
        fmt = object.__getattribute__(self, "fmt")
        val = getattr(delegate, name)
        if callable(val) and name not in delegate.__dict__:
            @wraps(val)
            def wrapper(*args, **kwds):
                arglst = [fmt(arg) for arg in args]
                arglst.extend("{}={}".format(k, fmt(v))
                              for k, v in sorted(kwds.items()))
                res = val(*args, **kwds)
                print("{}call {}({}) -> {}"
                      .format(prefix, name, ", ".join(arglst), fmt(res)))
                return res
            return wrapper
        else:
            if object.__getattribute__(self, "reads"):
                print("{}read {} = {}".format(prefix, name, fmt(val)))
            return val

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
        delegate = object.__getattribute__(self, "delegate")
        prefix = object.__getattribute__(self, "prefix")
        fmt = object.__getattribute__(self, "fmt")
        print("{}write {} = {}".format(prefix, name, fmt(val)))
        setattr(delegate, name, val)

    @classmethod
    def fmt(cls, val):
        if isinstance(val, frozenset):
            return ("frozenset([{}])".format
                    (", ".join(map(cls.fmt, sorted(val)))))
        if isinstance(val, set):
            return ("set([{}])".format
                    (", ".join(map(cls.fmt, sorted(val)))))
        r = repr(val)
        return r


def traceMethod(obj, meth, **kwds):
    r"""
    Trace the doings of a given method.
    It prints method entry with arguments,
    access to members and other methods during method execution
    as well as method exit with return value.

    INPUT:

    -  ``obj``: The object containing the method.

    -  ``meth``: The name of the method to be traced.

    EXAMPLE::

        sage: class Foo(object):
        ....:     def f(self):
        ....:         self.y = self.g(self.x)
        ....:     def g(self, arg):
        ....:         return arg + 1
        ....:
        sage: foo = Foo()
        sage: foo.x = 3
        sage: from sage.doctest.fixtures import traceMethod
        sage: traceMethod(foo, "f")
        sage: foo.f()
        enter f()
          read x = 3
          call g(3) -> 4
          write y = 4
        exit f -> None
    """
    f = getattr(obj, meth).__func__
    t = PropertyAccessTracerProxy(obj, **kwds)
    fmt = PropertyAccessTracerProxy.fmt
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
