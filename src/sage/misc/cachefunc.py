"""
Cached Functions and Methods

AUTHORS:

- William Stein (inspired by conversation with Justin Walker).
- Mike Hansen (added doctests and made it work with class methods).
- Willem Jan Palenstijn (add CachedMethodCaller for binding cached
  methods to instances).
- Tom Boothby (added DiskCachedFunction).
- Simon King (improved performance, more doctests).

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
import os

class CachedFunction(object):
    """
    Create a cached version of a function, which only recomputes
    values it hasn't already computed. Synonyme: ``cached_function``

    If ``f`` is a function, do either ``g = CachedFunction(f)``
    or ``g = cached_function(f)`` to make a cached version of ``f``,
    or put ``@cached_function`` right before the definition of ``f``
    (i.e., use Python decorators)::

        @cached_function
        def f(...):
            ....

    The inputs to the function must be hashable.

    EXAMPLES::

        sage: @cached_function
        ... def mul(x, y=2):
        ...     return x*y
        ...
        sage: mul(3)
        6

    We demonstrate that the result is cached, and that, moreover,
    the cache takes into account the various ways of providing
    default arguments::

        sage: mul(3) is mul(3,2)
        True
        sage: mul(3,y=2) is mul(3,2)
        True

    The user can clear the cache::

        sage: a = mul(4)
        sage: mul.clear_cache()
        sage: a is mul(4)
        False

    It is also possible to explicitly override the cache with
    a different value::

        sage: mul.set_cache('foo',5)
        sage: mul(5,2)
        'foo'

    """
    def __init__(self, f, classmethod=False):
        """
        Create a cached version of a function, which only recomputes
        values it hasn't already computed.

        If f is a function, do either g = CachedFunction(f) to make
        a cached version of f, or put @CachedFunction right before
        the definition of f (i.e., use Python decorators)::

            @CachedFunction
            def f(...):
                ....

        The inputs to the function must be hashable.

        TESTS::

            sage: g = CachedFunction(number_of_partitions)
            sage: g.__name__
            'number_of_partitions'
            sage: 'partitions' in g.__doc__
            True
            sage: g(5)
            7
            sage: g.cache
            {((5, None, 'default'), ()): 7}
            sage: def f(t=1): print(t)
            sage: h = CachedFunction(f)
            sage: w = walltime()
            sage: h(); h(1); h(t=1)
            1
            sage: walltime(w) < 2
            True

        """
        self._common_init(f, ArgumentFixer(f,classmethod=classmethod))
        self.cache = {}

    def _common_init(self, f, argumentfixer):
        """
        Perform initialization common to CachedFunction and CachedMethodCaller.

        TESTS::

            sage: @cached_function
            ... def test_cache(x):
            ...     return -x
            sage: hasattr(test_cache, '_argumentfixer')  # indirect doctest
            True
        """
        self.f = f
        if hasattr(f, "func_doc"):
            self.__doc__ = f.func_doc
        if hasattr(f, "func_name"):
            self.__name__ = f.func_name
        self.__module__ = f.__module__
        self._argumentfixer = argumentfixer

    def _sage_src_(self):
        """
        Returns the source code for the wrapped function.

        TESTS::

            sage: from sage.misc.sageinspect import sage_getsource
            sage: g = CachedFunction(number_of_partitions)
            sage: 'bober' in sage_getsource(g)  # indirect doctest
            True

        """
        from sage.misc.sageinspect import sage_getsource
        return sage_getsource(self.f)

    def __call__(self, *args, **kwds):
        """
        Return value from cache or call the wrapped function,
        caching the output.

        TESTS::

            sage: g = CachedFunction(number_of_partitions)
            sage: a = g(5)
            sage: g.get_cache()
            {((5, None, 'default'), ()): 7}
            sage: a = g(10^5)
            sage: a == number_of_partitions(10^5)
            True
            sage: a is g(10^5)
            True
            sage: a is number_of_partitions(10^5)
            False

        """
        # We shortcut a common case of no arguments
        if args or kwds:
            k = self._argumentfixer.fix_to_pos(*args, **kwds)
        else:
            try:
                k = self._default_key
            except AttributeError:
                k = self._default_key = self._argumentfixer.fix_to_pos()
        try:
            return self.cache[k]
        except KeyError:
            w = self.f(*args, **kwds)
            self.cache[k] = w
            return w

    def get_cache(self):
        """
        Returns the cache dictionary.

        EXAMPLES::

            sage: g = CachedFunction(number_of_partitions)
            sage: a = g(5)
            sage: g.get_cache()
            {((5, None, 'default'), ()): 7}

        """
        return self.cache

    def is_in_cache(self, *args, **kwds):
        """
        Checks if the argument list is in the cache.

        EXAMPLES::

            sage: class Foo:
            ...       def __init__(self, x):
            ...           self._x = x
            ...       @cached_method
            ...       def f(self, z, y=0):
            ...           return self._x*z+y
            ...
            sage: a = Foo(2)
            sage: a.f.is_in_cache(3)
            False
            sage: a.f(3)
            6
            sage: a.f.is_in_cache(3,y=0)
            True
        """
        return self._argumentfixer.fix_to_pos(*args, **kwds) in self.cache

    def set_cache(self, value, *args, **kwds):
        """
        Set the value for those args and keyword args
        Mind the unintuitive syntax (value first).
        Any idea on how to improve that welcome!

        EXAMPLES::

            sage: g = CachedFunction(number_of_partitions)
            sage: a = g(5)
            sage: g.get_cache()
            {((5, None, 'default'), ()): 7}
            sage: g.set_cache(17, 5)
            sage: g.get_cache()
            {((5, None, 'default'), ()): 17}
            sage: g(5)
            17

        DEVELOPER NOTE:

        Is there a way to use the following intuitive syntax?

        ::

            sage: g(5) = 19    # todo: not implemented
            sage: g(5)         # todo: not implemented
            19
        """
        self.cache[self._argumentfixer.fix_to_pos(*args, **kwds)] = value

    def get_key(self, *args, **kwds):
        """
        Returns the key in the cache to be used when args
        and kwds are passed in as parameters.

        EXAMPLES::

            sage: @cached_function
            ... def foo(x):
            ...    return x^2
            ...
            sage: foo(2)
            4
            sage: foo.get_key(2)
            ((2,), ())
            sage: foo.get_key(x=3)
            ((3,), ())
        """
        return self._argumentfixer.fix_to_pos(*args, **kwds)

    def __repr__(self):
        """
        EXAMPLES::

            sage: g = CachedFunction(number_of_partitions)
            sage: g
            Cached version of <function number_of_partitions at 0x...>
        """
        try:
            return "Cached version of %s"%self.f
        except AttributeError:
            return "Cached version of a method (pending reassignment)"

    def clear_cache(self):
        """
        Clear the cache dictionary.

        EXAMPLES::

            sage: g = CachedFunction(number_of_partitions)
            sage: a = g(5)
            sage: g.get_cache()
            {((5, None, 'default'), ()): 7}
            sage: g.clear_cache()
            sage: g.get_cache()
            {}
        """
        cache = self.cache
        for key in cache.keys():
            del cache[key]

    def precompute(self, arglist, num_processes=1):
        """

        Cache values for a number of inputs.  Do the computation
        in parallel, and only bother to compute values that we
        haven't already cached.

        EXAMPLES::

            sage: @cached_function
            ... def oddprime_factors(n):
            ...     l = [p for p,e in factor(n) if p != 2]
            ...     return len(l)
            sage: oddprime_factors.precompute(range(1,100), 4)
            sage: oddprime_factors(25) is oddprime_factors(25)
            True
        """
        from sage.parallel.decorate import parallel, normalize_input
        P = parallel(num_processes)(self.f)
        cache = self.cache
        get_key = self._argumentfixer.fix_to_pos
        new = lambda x: not cache.has_key(get_key(*x[0],**x[1]))
        arglist = filter(new, map(normalize_input, arglist))
        for ((args,kwargs), val) in P(arglist):
            self.set_cache(val, *args, **kwargs)


cached_function = CachedFunction

class CachedMethodPickle:
    """
    This class helps to unpickle cached methods.

    NOTE:

    Since trac ticket #8611, a cached method is an attribute
    of the instance (provided that it has a ``__dict__``).
    Hence, when pickling the instance, it would be attempted
    to pickle that attribute as well, but this is a problem,
    since functions can not be pickled, currently. Therefore,
    we replace the actual cached method by a place holder,
    that kills itself as soon as any attribute is requested.
    Then, the original cached attribute is reinstated.

    EXAMPLE::

        sage: R.<x, y, z> = PolynomialRing(QQ, 3)
        sage: I = R*(x^3 + y^3 + z^3,x^4-y^4)
        sage: I.groebner_basis()
        [y^5*z^3 - 1/4*x^2*z^6 + 1/2*x*y*z^6 + 1/4*y^2*z^6, x^2*y*z^3 - x*y^2*z^3 + 2*y^3*z^3 + z^6, x*y^3 + y^4 + x*z^3, x^3 + y^3 + z^3]
        sage: I.groebner_basis
        Cached version of <function groebner_basis at 0x...>

    We now pickle and unpickle the ideal. The cached method
    ``groebner_basis`` is replaced by a placeholder::

        sage: J = loads(dumps(I))
        sage: J.groebner_basis
        Pickle of the cached method "groebner_basis"

    But as soon as any other attribute is requested from the
    placeholder, it replaces itself by the cached method, and
    the entries of the cache are actually preserved::

        sage: J.groebner_basis.is_in_cache()
        True
        sage: J.groebner_basis
        Cached version of <function groebner_basis at 0x...>
        sage: J.groebner_basis() == I.groebner_basis()
        True

    AUTHOR:

    - Simon King (2011-01)
    """
    def __init__(self, inst, name):
        """
        INPUT:

        - ``inst`` - some instance.
        - ``name`` (string) - usually the name of an attribute
          of ``inst`` to which ``self`` is assigned.

        TEST::

            sage: from sage.misc.cachefunc import CachedMethodPickle
            sage: P = CachedMethodPickle(1, 'foo')
            sage: P
            Pickle of the cached method "foo"

        """
        self._instance = inst
        self._name = name
    def __repr__(self):
        """
        TEST::

            sage: R.<x, y, z> = PolynomialRing(QQ, 3)
            sage: I = R*(x^3 + y^3 + z^3,x^4-y^4)
            sage: G = I.groebner_basis()
            sage: J = loads(dumps(I))
            sage: J.groebner_basis  #indirect doctest
            Pickle of the cached method "groebner_basis"
        """
        return 'Pickle of the cached method "%s"'%self._name
    def __getattr__(self,s):
        """
        TEST::

            sage: R.<x, y, z> = PolynomialRing(QQ, 3)
            sage: I = R*(x^3 + y^3 + z^3,x^4-y^4)
            sage: G = I.groebner_basis()
            sage: J = loads(dumps(I))
            sage: J.groebner_basis
            Pickle of the cached method "groebner_basis"

        If an attribute of name ``s`` is requested (say,
        ``is_in_cache``), the attribute ``self._name`` of
        ``self._instance`` is deleted. Then, the attribute
        of name ``s`` of the attribute ``self._name`` of
        ``self._instance`` is requested. Since ``self._name``
        is a cached method defined for the class of
        ``self._instance``, retrieving the just-deleted
        attribute ``self._name`` succeeds.

        In that way, the unpickling of the cached method is
        finally accomplished::

            sage: J.groebner_basis.is_in_cache()  #indirect doctest
            True
            sage: J.groebner_basis
            Cached version of <function groebner_basis at 0x...>

        """
        self._instance.__dict__.__delitem__(self._name)
        return getattr(getattr(self._instance,self._name),s)

class CachedMethodCaller(CachedFunction):
    """
    Utility class that is used by :class:`CachedMethod` to bind a
    cached method to an instance.

    EXAMPLE::

        sage: class A:
        ...    @cached_method
        ...    def bar(self,x):
        ...        return x^2
        sage: a = A()
        sage: a.bar
        Cached version of <function bar at 0x...>
        sage: type(a.bar)
        <class 'sage.misc.cachefunc.CachedMethodCaller'>
        sage: a.bar(2) is a.bar(x=2)
        True

    """
    def __init__(self, cachedmethod, inst, cache=None, inst_in_key=False):
        """
        EXAMPLES::

            sage: class Foo:
            ...       def __init__(self, x):
            ...           self._x = x
            ...       @cached_method
            ...       def f(self):
            ...           return self._x^2
            ...
            sage: a = Foo(2)
            sage: a.f.get_cache()
            {}
            sage: a.f()
            4
            sage: a.f.get_cache()
            {((), ()): 4}
        """
        # initialize CachedFunction, but re-use the ArgumentFixer
        self._common_init(cachedmethod._cachedfunc.f, cachedmethod._cachedfunc._argumentfixer)
        self.cache = {} if cache is None else cache
        self._instance = inst
        self._inst_in_key = inst_in_key
        self._cachedmethod = cachedmethod

    def __reduce__(self):
        """
        The pickle of a :class:`CachedMethodCaller` unpickles
        to a :class:`CachedMethodPickle`, that is able to replace
        itself by a copy of the original :class:`CachedMethodCaller`.

        TEST::

            sage: R.<x, y, z> = PolynomialRing(QQ, 3)
            sage: I = R*(x^3 + y^3 + z^3,x^4-y^4)
            sage: G = I.groebner_basis()
            sage: J = loads(dumps(I))  #indirect doctest
            sage: J.groebner_basis
            Pickle of the cached method "groebner_basis"
            sage: J.groebner_basis.is_in_cache()
            True
            sage: J.groebner_basis
            Cached version of <function groebner_basis at 0x...>

        """
        return CachedMethodPickle,(self._instance,self.__name__)

    def __call__(self, *args, **kwds):
        """
        Call the cached method.

        TESTS::

            sage: class Foo:
            ...       @cached_method
            ...       def f(self, x,y=1):
            ...           return x+y
            ...
            sage: a = Foo()
            sage: a.f(1)  #indirect doctest
            2

        The result is cached, taking into account
        the three ways of providing (named) arguments::

            sage: a.f(5) is a.f(5,1)
            True
            sage: a.f(5) is a.f(5,y=1)
            True
            sage: a.f(5) is a.f(y=1,x=5)
            True

        We test that #5843 is fixed::

            sage: class Foo:
            ...       def __init__(self, x):
            ...           self._x = x
            ...       @cached_method
            ...       def f(self, y):
            ...           return self._x
            ...
            sage: a = Foo(2)
            sage: b = Foo(3)
            sage: a.f(b.f)
            2
        """
        # We shortcut a common case of no arguments
        if args or kwds:
            k = self.get_key(*args, **kwds)
        else:
            try:
                k = self._default_key
            except AttributeError:
                k = self._default_key = self.get_key()
        try:
            return self.cache[k]
        except KeyError:
            w = self._cachedmethod._instance_call(self._instance, *args, **kwds)
            self.cache[k] = w
            return w

    def get_key(self, *args, **kwds):
        """
        Convert arguments to the key for this instance's cache.

        EXAMPLES::

            sage: class Foo:
            ...       def __init__(self, x):
            ...           self._x = x
            ...       @cached_method
            ...       def f(self, y, z=0):
            ...           return self._x * y + z
            ...
            sage: a = Foo(2)
            sage: z = a.f(37)
            sage: k = a.f.get_key(37); k
            ((37, 0), ())
            sage: a.f.get_cache()[k] is z
            True

        Note that the method does not test whether there are
        too many arguments, or wrong argument names::

            sage: a.f.get_key(1,2,3,x=4,y=5,z=6)
            ((1, 2, 3), (('x', 4), ('y', 5), ('z', 6)))

        It does, however, take into account the different
        ways of providing named arguments, possibly with a
        default value::

            sage: a.f.get_key(5)
            ((5, 0), ())
            sage: a.f.get_key(y=5)
            ((5, 0), ())
            sage: a.f.get_key(5,0)
            ((5, 0), ())
            sage: a.f.get_key(5,z=0)
            ((5, 0), ())
            sage: a.f.get_key(y=5,z=0)
            ((5, 0), ())

        """
        if self._inst_in_key:
            return (self._instance,self._argumentfixer.fix_to_pos(*args,**kwds))
        return self._argumentfixer.fix_to_pos(*args,**kwds)

    def __get__(self, inst, cls=None):
        """
        Get a :class:`CachedMethodCaller` bound to a specific
        instance of the class of the cached method.

        NOTE:

        :class:`CachedMethodCaller` has a separate ``__get__``
        since the categories framework creates and caches the
        return value of ``CachedMethod.__get__`` with
        ``inst==None``.

        TESTS:

        Due to the separate ``__get__`` method, it is possible
        to define a cached method in one class and use it as
        an attribute of another class.

            sage: class Foo:
            ...       @cached_method
            ...       def f(self, y):
            ...           return y - 1
            sage: class Bar:
            ...       f = Foo.f
            sage: b1 = Bar()
            sage: b2 = Bar()

        The :class:`CachedMethod` is replaced by an instance
        of :class:`CachedMethodCaller` that (by trac ticket
        #8611) is set as an attribute. Hence, we have::

            sage: b1.f is b1.f
            True

        Any instance of ``Bar`` gets its own instance of
        :class:`CachedMethodCaller``::

            sage: b1.f is b2.f
            False

        The method caller knows the instance that it belongs
        to::

            sage: Foo.f._instance is None
            True
            sage: b1.f._instance is b1
            True
            sage: b2.f._instance is b2
            True

        """
        Caller = CachedMethodCaller(self._cachedmethod, inst, cache=self._cachedmethod._get_instance_cache(inst), inst_in_key=self._inst_in_key)
        try:
            setattr(inst,self._cachedmethod._cachedfunc.f.__name__, Caller)
        except AttributeError:
            pass
        return Caller

class CachedMethod(object):
    """
    A decorator that creates a cached version of an instance
    method of a class.

    NOTE:

    For proper behavior, the method must be a pure function
    (no side effects). Arguments to the method must be hashable.

    EXAMPLES::

        sage: class Foo(object):
        ...       @cached_method
        ...       def f(self, t, x=2):
        ...           print 'computing'
        ...           return t**x
        sage: a = Foo()

    The example shows that the actual computation
    takes place only once, and that the result is
    identic for equivalent input::

        sage: res = a.f(3, 2); res
        computing
        9
        sage: a.f(t = 3, x = 2) is res
        True
        sage: a.f(3) is res
        True

    """
    def __init__(self, f):
        """
        EXAMPLES::

            sage: class Foo:
            ...       def __init__(self, x):
            ...           self._x = x
            ...       @cached_method
            ...       def f(self,n):
            ...           return self._x^n
            ...
            sage: a = Foo(2)
            sage: a.f(2)
            4

        The computations in method ``f`` are
        stored in a dictionary assigned to the
        instance ``a``::

            sage: hasattr(a, '_cache__f')
            True
            sage: a._cache__f
            {((2,), ()): 4}

        As a shortcut, useful to speed up internal computations,
        the same dictionary is also available as an attribute
        of the ``CachedMethodCaller``::

            sage: type(a.f)
            <class 'sage.misc.cachefunc.CachedMethodCaller'>
            sage: a.f.cache is a._cache__f
            True

        """
        self._cache_name = '_cache__' + f.__name__
        self._cachedfunc = CachedFunction(f, classmethod=True)
        self._argument_fix_to_pos = self._cachedfunc._argumentfixer.fix_to_pos

    def _instance_call(self, inst, *args, **kwds):
        """
        Call the cached method *without* using the cache.

        INPUT:

        - ``inst`` - an instance at which the method is to be called
        - Further positional or named arguments.

        EXAMPLES::

            sage: class Foo(object):
            ...       def __init__(self, x):
            ...           self._x = x
            ...       @cached_method
            ...       def f(self):
            ...           return self._x^2
            ...
            sage: a = Foo(2)
            sage: a.f()
            4

        Usually, a cached meth is indeed cached::

            sage: a.f() is a.f()
            True

        However, when it becomes necessary, one can call
        it without using the cache::

            sage: a.f._cachedmethod._instance_call(a) is a.f()
            False
            sage: a.f._cachedmethod._instance_call(a) == a.f()
            True

        It is also possible to call the method with a
        different instance::

            sage: b = Foo(3)
            sage: b.f()
            9
            sage: a.f._cachedmethod._instance_call(b)
            9

        """
        return self._cachedfunc.f(inst, *args, **kwds)

    def _get_instance_cache(self, inst):
        """
        Returns the cache dictionary.

        TESTS::

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
            sage: a.f._cachedmethod._get_instance_cache(a)
            {((), ()): 4}
        """
        try:
            return inst.__dict__.setdefault(self._cache_name, {})
        except AttributeError:
            return {}

    def __get__(self, inst, cls=None):
        """
        Get a CachedMethodCaller bound to this specific instance of
        the class of the cached method.

        TESTS::

            sage: class Foo:
            ...       @cached_method
            ...       def f(self):
            ...           return 1
            sage: a = Foo()
            sage: type(a.f)
            <class 'sage.misc.cachefunc.CachedMethodCaller'>

        By trac ticket #8611, the CachedMethodCaller is set as
        an attribute of the instance ``a``, replacing the original
        cached attribute. Therefore, the ``__get__`` method will
        be used only once, which saves much time. Hence, we have::

            sage: a.f is a.f
            True

        """
        Caller = CachedMethodCaller(self, inst, cache=self._get_instance_cache(inst))
        try:
            setattr(inst,self._cachedfunc.f.__name__, Caller)
        except AttributeError,msg:
            pass
        return Caller

        # Note: a simpler approach to this would be
        # def caller(*args, **kwds):
        #     return self._instance_call(inst, *args, **kwds)
        # return caller
        # The disadvantage to this is that it does not provide
        # is_in_cache(), set_cache(), clear_cache(), ... methods.


cached_method = CachedMethod

class CachedInParentMethod(CachedMethod):
    """
    A decorator that creates a cached version of an instance
    method of a class.

    In contrast to :class:`CachedMethod`,
    the cache dictionary is an attribute of the parent of
    the instance to which the method belongs.

    ASSUMPTION:

    This way of caching works only if

    - the instances *have* a parent,
    - the parent allows assignment of attributes, and
    - the instances are hashable (they are part of the cache key).

    NOTE:

    For proper behavior, the method must be a pure function
    (no side effects). Arguments to the method must be hashable.

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
        ...       def __eq__(self, other):
        ...           return self._x^2 == other._x^2
        ...       def __hash__(self):
        ...           return hash(self._x^2)
        ...       def __repr__(self):
        ...           return 'My %s'%self._x
        ...       @cached_in_parent_method
        ...       def f(self):
        ...           return self._x^3
        ...
        sage: a = Foo(2)
        sage: a.f()
        8
        sage: a.f.get_cache()   #indirect doctest
        {(My 2, ((), ())): 8}

    Since the key for the cache depends on equality of
    the instances, we obtain *identical* result for
    *equal* instances - even though in this particular
    example the result of the method should be different
    for ``a`` and ``b``::

        sage: b = Foo(-2)
        sage: a is not b
        True
        sage: a == b
        True
        sage: b.f() is a.f()
        True

    Non-equal instances do not share the result of
    the cached method, but they do share the cache::

        sage: c = Foo(3)
        sage: c.f()
        27
        sage: c.f.get_cache() is a.f.get_cache() #indirect doctest
        True

    Note that the cache is also available as an
    attribute of the cached method, which speeds
    up internal computations::

        sage: a.f.cache is b.f.get_cache() is c.f._cachedmethod._get_instance_cache(c)
        True
    """
    def __init__(self, f):
        """
        Constructs a new method with cache stored in the parent of the instance.

        See also ``cached_method`` and ``cached_function``.

        EXAMPLES::

            sage: class MyParent(Parent):
            ...       pass
            sage: class Foo:
            ...       def __init__(self, x):
            ...           self._x = x
            ...       _parent = MyParent()
            ...       def parent(self):
            ...           return self._parent
            ...       @cached_in_parent_method  #indirect doctest
            ...       def f(self):
            ...           return self._x^2
            ...
            sage: a = Foo(2)
            sage: a.f()
            4
            sage: hasattr(a.parent(), '_cache__element_f')
            True

        For speeding up internal computations, this dictionary
        is also accessible as an attribute of the CachedMethodCaller
        (by trac ticket #8611)::

            sage: a.parent()._cache__element_f is a.f.cache
            True
        """
        self._cache_name = '_cache__' + 'element_' + f.__name__
        self._cachedfunc = CachedFunction(f, classmethod=True)
        self._argument_fix_to_pos = self._cachedfunc._argumentfixer.fix_to_pos

    def _get_instance_cache(self, inst):
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
            ...       def __eq__(self, other):
            ...           return self._x^2 == other._x^2
            ...       def __hash__(self):
            ...           return hash(self._x^2)
            ...       def __repr__(self):
            ...           return 'My %s'%self._x
            ...       @cached_in_parent_method
            ...       def f(self):
            ...           return self._x^3
            ...
            sage: a = Foo(2)
            sage: a.f()
            8
            sage: a.f.get_cache()   #indirect doctest
            {(My 2, ((), ())): 8}

        Since the key for the cache depends on equality of
        the instances, we obtain *identical* result for
        *equal* instance - even though in this particular
        example the result is wrong::

            sage: b = Foo(-2)
            sage: a is not b
            True
            sage: a == b
            True
            sage: b.f() is a.f()
            True

        Non-equal instances do not share the result of
        the cached method, but they do share the cache::

            sage: c = Foo(3)
            sage: c.f()
            27
            sage: c.f.get_cache() is a.f.get_cache() #indirect doctest
            True

        Note that the cache is also available as an
        attribute of the cached method, which speeds
        up internal computations::

            sage: a.f.cache is b.f.get_cache() is c.f._cachedmethod._get_instance_cache(c)
            True

        """
        try:
            return inst.parent().__dict__.setdefault(self._cache_name, {})
        except AttributeError:
            return {}

    def __get__(self, inst, cls=None):
        """
        Get a CachedMethodCaller bound to this specific instance of
        the class of the cached-in-parent method.
        """
        Caller = CachedMethodCaller(self, inst, cache=self._get_instance_cache(inst), inst_in_key=True)
        try:
            setattr(inst,self._cachedfunc.f.__name__, Caller)
        except AttributeError:
            pass
        return Caller

cached_in_parent_method = CachedInParentMethod

class FileCache:
    """
    FileCache is a dictionary-like class which stores keys and
    values on disk.  The keys take the form of a tuple (A,K)

    - A is a tuple of objects t where each t is an exact
      object which is uniquely identified by a short string.

    - K is a tuple of tuples (s,v) where s is a valid
      variable name and v is an exact object which is uniquely
      identified by a short string with letters [a-zA-Z0-9-._]

    The primary use case is the DiskCachedFunction.  If
    ``memory_cache == True``, we maintain a cache of objects seen
    during this session in memory -- but we don't load them from
    disk until necessary.  The keys and values are stored in a
    pair of files:

    - ``prefix-argstring.key.sobj`` contains the ``key`` only,
    - ``prefix-argstring.sobj`` contains the tuple ``(key,val)``

    where ``self[key] == val``.

    NOTE:

    We assume that each FileCache lives in its own directory.
    Use **extreme** caution if you wish to break that assumption.
    """
    def __init__(self, dir, prefix = '', memory_cache = False):
        """
        EXAMPLES::

            sage: from sage.misc.cachefunc import FileCache
            sage: dir = tmp_dir()
            sage: FC = FileCache(dir, memory_cache = True)
            sage: FC[((),())] = 1
            sage: FC[((1,2),())] = 2
            sage: FC[((),())]
            1
        """
        from sage.misc.misc import sage_makedirs
        if len(dir) == 0 or dir[-1] != '/':
            dir += '/'
        self._dir = dir
        sage_makedirs(dir)

        self._prefix = prefix + '-'

        if memory_cache:
            self._cache = {}
        else:
            self._cache = None

    def file_list(self):
        """
        Returns the list of files corresponding to self.

        EXAMPLES::

            sage: from sage.misc.cachefunc import FileCache
            sage: dir = tmp_dir()
            sage: FC = FileCache(dir, memory_cache = True, prefix='t')
            sage: FC[((),())] = 1
            sage: FC[((1,2),())] = 2
            sage: FC[((1,),(('a',1),))] = 3
            sage: for f in sorted(FC.file_list()): print f[len(dir):]
            /t-.key.sobj
            /t-.sobj
            /t-1_2.key.sobj
            /t-1_2.sobj
            /t-a-1.1.key.sobj
            /t-a-1.1.sobj
        """
        files = []
        prefix = self._prefix
        dir = self._dir
        l = len(prefix)
        for f in os.listdir(dir):
            if f[:l] == prefix:
                files.append( dir + f )
        return files

    def items(self):
        """
        Returns a list of tuples ``(k,v)`` where ``self[k] = v``.

        EXAMPLES::

            sage: from sage.misc.cachefunc import FileCache
            sage: dir = tmp_dir()
            sage: FC = FileCache(dir, memory_cache = False)
            sage: FC[((),())] = 1
            sage: FC[((1,2),())] = 2
            sage: FC[((1,),(('a',1),))] = 3
            sage: I = FC.items()
            sage: I.sort(); print I
            [(((), ()), 1), (((1,), (('a', 1),)), 3), (((1, 2), ()), 2)]
        """
        return [(k,self[k]) for k in self]

    def values(self):
        """
        Returns a list of values that are stored in ``self``.

        EXAMPLES::

            sage: from sage.misc.cachefunc import FileCache
            sage: dir = tmp_dir()
            sage: FC = FileCache(dir, memory_cache = False)
            sage: FC[((),())] = 1
            sage: FC[((1,2),())] = 2
            sage: FC[((1,),(('a',1),))] = 3
            sage: FC[((),(('a',1),))] = 4
            sage: v = FC.values()
            sage: v.sort(); print v
            [1, 2, 3, 4]
        """
        return [self[k] for k in self]

    def __iter__(self):
        """
        Returns a list of keys of ``self``.

        EXAMPLES::

            sage: from sage.misc.cachefunc import FileCache
            sage: dir = tmp_dir()
            sage: FC = FileCache(dir, memory_cache = False)
            sage: FC[((),())] = 1
            sage: FC[((1,2),())] = 2
            sage: FC[((1,),(('a',1),))] = 3
            sage: for k in sorted(FC): print k
            ((), ())
            ((1,), (('a', 1),))
            ((1, 2), ())
        """
        from sage.structure.sage_object import load

        for f in self.file_list():
            if f[-9:] == '.key.sobj':
                 yield load(f)

    def keys(self):
        """
        Returns a list of keys ``k`` where ``self[k]`` is defined.

        EXAMPLES::

            sage: from sage.misc.cachefunc import FileCache
            sage: dir = tmp_dir()
            sage: FC = FileCache(dir, memory_cache = False)
            sage: FC[((),())] = 1
            sage: FC[((1,2),())] = 2
            sage: FC[((1,),(('a',1),))] = 3
            sage: K = FC.keys()
            sage: K.sort(); print K
            [((), ()), ((1,), (('a', 1),)), ((1, 2), ())]
        """
        return [k for k in self]

    def _filename(self, key):
        """
        Computes the filename associated with a certain key.

        EXAMPLES::

            sage: from sage.misc.cachefunc import FileCache
            sage: dir = tmp_dir() + '/'
            sage: FC = FileCache(dir, memory_cache = False, prefix='foo')
            sage: N = FC._filename(((1,2), (('a',1),('b',2))))
            sage: print N[len(dir):]
            foo-a-1_b-2.1_2
            sage: N = FC._filename(((), (('a',1),('b',2))))
            sage: print N[len(dir):]
            foo-a-1_b-2
            sage: N = FC._filename(((1,2), ()))
            sage: print N[len(dir):]
            foo-1_2
        """
        a,k = key
        kwdstr = '_'.join('%s-%s'%x for x in k)
        argstr = '_'.join('%s'%x for x in a)
        if kwdstr and argstr:
            keystr = kwdstr + '.' + argstr
        else:
            keystr = kwdstr + argstr
        return self._dir + self._prefix + keystr

    def has_key(self, key):
        """
        Returns ``True`` if ``self[key]`` is defined and False otherwise.

        EXAMPLES::

            sage: from sage.misc.cachefunc import FileCache
            sage: dir = tmp_dir() + '/'
            sage: FC = FileCache(dir, memory_cache = False, prefix='foo')
            sage: k = ((),(('a',1),))
            sage: FC[k] = True
            sage: FC.has_key(k)
            True
            sage: FC.has_key(((),()))
            False
        """
        return os.path.exists(self._filename(key) + '.key.sobj')

    def __getitem__(self, key):
        """
        Returns the value set by ``self[key] = val``, in this session
        or a previous one.

        EXAMPLES::

            sage: from sage.misc.cachefunc import FileCache
            sage: dir = tmp_dir() + '/'
            sage: FC1 = FileCache(dir, memory_cache = False, prefix='foo')
            sage: FC2 = FileCache(dir, memory_cache = False, prefix='foo')
            sage: k = ((),(('a',1),))
            sage: t = randint(0, 1000)
            sage: FC1[k] = t
            sage: FC2[k] == FC1[k] == t
            True
            sage: FC1[(1,2),(('a',4),('b',2))]
            Traceback (most recent call last):
            ...
            KeyError: ((1, 2), (('a', 4), ('b', 2)))

        """
        from sage.structure.sage_object import load

        cache = self._cache
        if cache is not None:
            if cache.has_key(key):
                return cache[key]

        f = self._filename(key) + '.sobj'
        try:
            k,v = load(f)
        except IOError:
            raise KeyError, key
        if k != key:
            raise RuntimeError, "cache corrupted"

        if cache is not None:
            cache[key] = v
        return v

    def __setitem__(self, key, value):
        """
        Sets ``self[key] = value`` and stores both key and value on
        disk.

        EXAMPLES::

            sage: from sage.misc.cachefunc import FileCache
            sage: dir = tmp_dir() + '/'
            sage: FC1 = FileCache(dir, memory_cache = False, prefix='foo')
            sage: FC2 = FileCache(dir, memory_cache = False, prefix='foo')
            sage: k = ((),(('a',1),))
            sage: t = randint(0, 1000)
            sage: FC1[k] = t
            sage: FC2[k] == t
            True
            sage: FC1[k] = 2000
            sage: FC2[k]!= t
            True
        """
        from sage.structure.sage_object import save

        f = self._filename(key)

        save(key, f+'.key.sobj')
        save((key,value), f + '.sobj')
        if self._cache is not None:
            self._cache[key] = value

    def __delitem__(self, key):
        """
        Delete the key,value pair from self and unlink the associated
        files from the file cache.

        EXAMPLES::

            sage: from sage.misc.cachefunc import FileCache
            sage: dir = tmp_dir() + '/'
            sage: FC1 = FileCache(dir, memory_cache = False, prefix='foo')
            sage: FC2 = FileCache(dir, memory_cache = False, prefix='foo')
            sage: k = ((),(('a',1),))
            sage: t = randint(0, 1000)
            sage: FC1[k] = t
            sage: del FC2[k]
            sage: FC1.has_key(k)
            False
       """
        f = self._filename(key)
        cache = self._cache
        if cache is not None and cache.has_key(key):
            del self._cache[key]
        if os.path.exists(f + '.sobj'):
            os.remove(f + '.sobj')
        if  os.path.exists(f + '.key.sobj'):
           os.remove(f + '.key.sobj')


class DiskCachedFunction(CachedFunction):
    """
    Works similar to CachedFunction, but instead, we keep the
    cache on disk (optionally, we keep it in memory too).

    EXAMPLES::

        sage: from sage.misc.cachefunc import DiskCachedFunction
        sage: dir = tmp_dir()
        sage: factor = DiskCachedFunction(factor, dir, memory_cache=True)
        sage: f = factor(2775); f
        3 * 5^2 * 37
        sage: f is factor(2775)
        True
    """
    def __init__(self, f, dir, memory_cache=False):
        """
        EXAMPLES::

            sage: from sage.misc.cachefunc import DiskCachedFunction
            sage: def foo(x): sleep(x)
            sage: dir = tmp_dir()
            sage: bar = DiskCachedFunction(foo, dir, memory_cache = False)
            sage: w = walltime()
            sage: for i in range(10): bar(1)
            sage: walltime(w) < 2
            True
        """
        CachedFunction.__init__(self, f)
        prefix = f.__name__
        self.cache = FileCache(dir, prefix=prefix, memory_cache = memory_cache)


class disk_cached_function:
    """
    Decorator for :class:`DiskCachedFunction`.

    EXAMPLES::

        sage: dir = tmp_dir()
        sage: @disk_cached_function(dir)
        ... def foo(x): return next_prime(2^x)%x
        sage: x = foo(200);x
        11
        sage: @disk_cached_function(dir)
        ... def foo(x): return 1/x
        sage: foo(200)
        11
        sage: foo.clear_cache()
        sage: foo(200)
        1/200
    """
    def __init__(self, dir, memory_cache = False):
        """
        EXAMPLES::

            sage: dir = tmp_dir()
            sage: @disk_cached_function(dir, memory_cache=True)
            ... def foo(x): return next_prime(2^x)
            sage: x = foo(200)
            sage: x is foo(200)
            True
            sage: @disk_cached_function(dir, memory_cache=False)
            ... def foo(x): return next_prime(2^x)
            sage: x is foo(200)
            False
        """
        self._dir = dir
        self._memory_cache = memory_cache

    def __call__(self, f):
        """
        EXAMPLES::

            sage: dir = tmp_dir()
            sage: @disk_cached_function(dir)
            ... def foo(x): return ModularSymbols(x)
            sage: foo(389)
            Modular Symbols space of dimension 65 for Gamma_0(389) of weight 2 with sign 0 over Rational Field
        """
        return DiskCachedFunction(f, self._dir, memory_cache = self._memory_cache)

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
