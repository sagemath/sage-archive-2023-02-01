"""
Cached Functions

AUTHOR:
    -- William Stein (inspired by conversation with Justin Walker).
    -- Mike Hansen (added doctests and made it work with class methods).
"""
########################################################################
#       Copyright (C) 2008 William Stein <wstein@gmail.com>
#                          Mike Hansen <mhansen@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
########################################################################
from function_mangling import ArgumentFixer

class CachedFunction(object):
    def __init__(self, f, classmethod=False):
        """
        Create a cached version of a function, which only recomputes
        values it hasn't already computed.

        If f is a function, do either g = CachedFunction(f) to make
        a cached version of f, or put @CachedFunction right before
        the definition of f (i.e., use Python decorators):

        @CachedFunction
        def f(...):
            ....

        The inputs to the function must be hashable.

        EXAMPLES:
            sage: g = CachedFunction(number_of_partitions)
            sage: g.__name__
            'number_of_partitions'
            sage: 'partitions' in g.__doc__
            True
            sage: g(5)
            7
            sage: g.cache
            {((5, None, 'default'), ()): 7}
            sage: def sleep1(t=1): sleep(t)
            sage: h = CachedFunction(sleep1)
            sage: w = walltime()
            sage: h(); h(1); h(t=1)
            sage: walltime(w) < 2
            True

        """
        self.f = f
        self.cache = {}
        if hasattr(f, "func_doc"):
            self.__doc__ = f.func_doc
        if hasattr(f, "func_name"):
            self.__name__ = f.func_name
        self.__module__ = f.__module__
        self._argumentfixer = ArgumentFixer(f,classmethod=classmethod)

    def _sage_src_(self):
        """
        Returns the source code for the wrapped function.

        EXAMPLES:
            sage: from sage.misc.sageinspect import sage_getsource
            sage: g = CachedFunction(number_of_partitions)
            sage: 'bober' in sage_getsource(g)
            True

        """
        from sage.misc.sageinspect import sage_getsource
        return sage_getsource(self.f)


    def __call__(self, *args, **kwds):
        """
        EXAMPLES:
            sage: g = CachedFunction(number_of_partitions)
            sage: a = g(5)
            sage: g.get_cache()
            {((5, None, 'default'), ()): 7}
            sage: a = g(10^5)
            sage: a == number_of_partitions(10^5)
            True
        """
        cache = self.get_cache()
        k = self.get_key(*args, **kwds)
        if cache.has_key(k):
            return cache[k]
        w = self.f(*args, **kwds)
        cache[k] = w
        return w

    def get_cache(self):
        """
        Returns the cache dictionary.

        EXAMPLES:
            sage: g = CachedFunction(number_of_partitions)
            sage: a = g(5)
            sage: g.get_cache()
            {((5, None, 'default'), ()): 7}

        """
        return self.cache

    def is_in_cache(self, *args, **kwds):
        """
        EXAMPLES:
            sage: class Foo:
            ...       def __init__(self, x):
            ...           self._x = x
            ...       @cached_method
            ...       def f(self, z):
            ...           return self._x*z
            ...
            sage: a = Foo(2)
            sage: a.f.is_in_cache(3)
            False
            sage: a.f(3)
            6
            sage: a.f.is_in_cache(3)
            True
        """
        cache = self.get_cache()
        return self.get_key(*args, **kwds) in cache

    def set_cache(self, value, *args, **kwds):
        """
        Set the value for those args and keyword args
        Mind the unintuitive syntax (value first)
        Any idea on how to improve that welcome

        EXAMPLES:
            sage: g = CachedFunction(number_of_partitions)
            sage: a = g(5)
            sage: g.get_cache()
            {((5, None, 'default'), ()): 7}
            sage: g.set_cache(17, 5)
            sage: g.get_cache()
            {((5, None, 'default'), ()): 17}
            sage: g(5)
            17

          Is there a way to use the following intuitive syntax?
            sage: g(5) = 19    # todo: not implemented
            sage: g(5)         # todo: not implemented
            19
        """
        cache = self.get_cache()
        cache[self.get_key(*args, **kwds)] = value


    def get_key(self, *args, **kwds):
        """
        Returns the key in the cache to be used when args
        and kwds are passed in as parameters.

        EXAMPLES:
            sage: class Foo:
            ...       def __init__(self, x):
            ...           self._x = x
            ...       @cached_method
            ...       def f(self):
            ...           return self._x^2
            ...
            sage: a = Foo(2)
            sage: a.f.get_key()
            ((), ())
        """
        return self._argumentfixer.fix_to_pos(*args, **kwds)

    def __repr__(self):
        """
        EXAMPLES:
            sage: g = CachedFunction(number_of_partitions)
            sage: g
            Cached version of <function number_of_partitions at 0x...>
        """
        return "Cached version of %s"%self.f

    def clear_cache(self):
        """
        Clear the cache dictionary.

        EXAMPLES:
            sage: g = CachedFunction(number_of_partitions)
            sage: a = g(5)
            sage: g.get_cache()
            {((5, None, 'default'), ()): 7}
            sage: g.clear_cache()
            sage: g.get_cache()
            {}
        """
        cache = self.get_cache()
        for key in cache.keys():
            del cache[key]


cached_function = CachedFunction

class CachedMethod(CachedFunction):
    def __init__(self, f):
        """
        EXAMPLES:
            sage: class Foo:
            ...       def __init__(self, x):
            ...           self._x = x
            ...       @cached_method
            ...       def f(self):
            ...           return self._x^2
            ...
            sage: Foo.f._cache_name
            '_cache__f'
        """
        self._cache_name = '_cache__' + f.__name__
        CachedFunction.__init__(self, f, classmethod=True)

    def __call__(self, *args, **kwds):
        """
        EXAMPLES:
            sage: class Foo(object):
            ...       def __init__(self, x):
            ...           self._x = x
            ...       @cached_method
            ...       def f(self):
            ...           return self._x^2
            ...       def g(self):
            ...           return self._x^2
            ...
            sage: a = Foo(2)
            sage: a.f()
            4
            sage: a.f() is a.f()
            True
            sage: a.g()
            4
            sage: a.g() is a.g()
            False
            sage: b = Foo(3)
            sage: b.f()
            9

        """
        cache = self.get_cache()
        key = self.get_key(*args, **kwds)
        if cache.has_key(key):
            return cache[key]
        else:
            cache[key] = self.f(self._instance, *args, **kwds)
            return cache[key]

    def get_cache(self):
        """
        Returns the cache dictionary.

        EXAMPLES:
            sage: class Foo:
            ...       def __init__(self, x):
            ...           self._x = x
            ...       @cached_method
            ...       def f(self):
            ...           return self._x^2
            ...
            sage: a = Foo(2)
            sage: a.f()
            4
            sage: a.f.get_cache()
            {((), ()): 4}

        """
        return self._instance.__dict__.setdefault(self._cache_name, {})

    def __get__(self, inst, cls=None):
        """
        This is needed to allow CachedFunction to decorate
        methods.
        """
        self._instance = inst
        return self

cached_method = CachedMethod

class CachedInParentMethod(CachedMethod):
    def __init__(self, f):
        """
        Constructs a new method with cache stored in the parent of the instance.

        See also ``cached_method`` and ``cached_function``.

        EXAMPLES::

            sage: class Foo:
            ...       def __init__(self, x):
            ...           self._x = x
            ...       @cached_in_parent_method
            ...       def f(self):
            ...           return self._x^2
            ...
            sage: Foo.f._cache_name
            '_cache__element_f'
        """
        self._cache_name = '_cache__' + 'element_' + f.__name__
        CachedFunction.__init__(self, f, classmethod=True)

    def get_key(self, *args, **kwds):
        """
        Returns the key used to lookup in the cache dictionary.

        EXAMPLES::

            sage: class MyParent(Parent):
            ...       pass
            ...
            sage: class Foo:
            ...       def __init__(self, x):
            ...           self._x = x
            ...       _parent = MyParent()
            ...       def parent(self):
            ...           return self._parent
            ...       def __repr__(self):
            ...           return str(self._x)
            ...       @cached_in_parent_method
            ...       def f(self, *args, **keywords):
            ...           return self._x^2
            ...
            sage: a = Foo(2)
            sage: a.f.get_key()
            ((2,), ())
            sage: a = Foo(2)
            sage: a.f.get_key(1,3)
            ((2, 1, 3), ())
            sage: a.f.get_key(1,3, bla=4)
            ((2, 1, 3), (('bla', 4),))
        """
        return super(CachedInParentMethod, self).get_key(self._instance, *args, **kwds)

    def get_cache(self):
        """
        Returns the cache dictionary, which is stored in the parent.

        EXAMPLES::

            sage: class MyParent(Parent):
            ...       pass
            ...
            sage: class Foo:
            ...       def __init__(self, x):
            ...           self._x = x
            ...       _parent = MyParent()
            ...       def parent(self):
            ...           return self._parent
            ...       def __repr__(self):
            ...           return str(self._x)
            ...       @cached_in_parent_method
            ...       def f(self):
            ...           return self._x^2
            ...
            sage: a = Foo(2)
            sage: a.f()
            4
            sage: a.f.get_cache()
            {((2,), ()): 4}
            sage: b = Foo(2)
            sage: a is not b
            True
            sage: b.f.get_cache()
            {((2,), ()): 4}
            sage: c = Foo(3)
            sage: c.f()
            9
            sage: c.f.get_cache()
            {((2,), ()): 4, ((3,), ()): 9}

        """
        return self._instance.parent().__dict__.setdefault(self._cache_name, {})


cached_in_parent_method = CachedInParentMethod

class ClearCacheOnPickle(object):
    """
    This class implements an appropriate __getstate__ method that
    clears the cache of the methods (see @cached_method) before
    passing them on to the caller, typically the pickle and copy modules.

    The implemented __getstate__ method calls the __getstate__ methods
    of classes later in the method resolution order. Therefore,
    classes which wants this behaviour should inherit first from this
    one.

    """

    def __getstate__(self):
        """
        """
        return dict( (key, value) for (key, value) in super(ClearCacheOnPickle, self).__getstate__().iteritems() if not (type(key) == str and key[0:8] == '_cache__') )
