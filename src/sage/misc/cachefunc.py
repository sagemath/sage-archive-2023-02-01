"""
Cached Functions

AUTHOR:
    -- William Stein (inspired by conversation with Justin Walker).
    -- Mike Hansen (added doctests and made it work with class methods).
    -- Willem Jan Palenstijn (add CachedMethodCaller for binding
                              cached methods to instances)
    -- Tom Boothby (added DiskCachedFunction)
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

        EXAMPLES::

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
        cache = self.get_cache()
        cache[self.get_key(*args, **kwds)] = value

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
        return "Cached version of %s"%self.f

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
        cache = self.get_cache()
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
        cache = self.get_cache()
        new = lambda x: not cache.has_key(self.get_key(*x[0],**x[1]))
        arglist = filter(new, map(normalize_input, arglist))
        for ((args,kwargs), val) in P(arglist):
            self.set_cache(val, *args, **kwargs)


cached_function = CachedFunction

class CachedMethodCaller(CachedFunction):
    """
    Utility class that is used by CachedMethod to bind a
    cached method to an instance.
    """
    def __init__(self, cachedmethod, inst):
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
        self._instance = inst
        self._cachedmethod = cachedmethod

    def __call__(self, *args, **kwds):
        """
        Call the cached method.

        TESTS::

            sage: class Foo:
            ...       @cached_method
            ...       def f(self, x):
            ...           return x+1
            ...
            sage: a = Foo()
            sage: a.f(1)
            2

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
        return self._cachedmethod._instance_call(self._instance, *args, **kwds)

    def get_cache(self, *args, **kwds):
        """
        Retrieve the cache for the instance.

        EXAMPLES::

            sage: class Foo:
            ...       def __init__(self, x):
            ...           self._x = x
            ...       @cached_method
            ...       def f(self, y):
            ...           return self._x * y
            ...
            sage: a = Foo(2)
            sage: a.f.get_cache()
            {}
            sage: a.f(37)
            74
            sage: a.f.get_cache()
            {((37,), ()): 74}
        """
        return self._cachedmethod._get_instance_cache(self._instance)

    def get_key(self, *args, **kwds):
        """
        Convert arguments to the key for this instance's cache.

        EXAMPLES::

            sage: class Foo:
            ...       def __init__(self, x):
            ...           self._x = x
            ...       @cached_method
            ...       def f(self, y):
            ...           return self._x * y
            ...
            sage: a = Foo(2)
            sage: z = a.f(37)
            sage: k = a.f.get_key(37); k
            ((37,), ())
            sage: a.f.get_cache()[k] == z
            True

        """
        return self._cachedmethod._get_instance_key(self._instance, *args, **kwds)

    def __get__(self, inst, cls=None):
        """
        Get a CachedMethodCaller bound to this specific instance of
        the class of the cached method.

        CachedMethodCaller has a separate __get__ since
        the categories framework creates and caches the return
        value of CachedMethod.__get__ with inst==None.

        TESTS::

            sage: class Foo:
            ...       @cached_method
            ...       def f(self, y):
            ...           return y - 1
            sage: class Bar:
            ...       f = Foo.f
            sage: b = Bar()
            sage: b.f is b.f
            False
            sage: b.f._instance is None
            False
        """
        return CachedMethodCaller(self._cachedmethod, inst)


class CachedMethod(object):
    """
    A decorator that creates a cached version of an instance method of
    a class. For proper behavior, the method must be a pure function (no side effects).
    Arguments to the method must be hashable.

    EXAMPLES::

        sage: class Foo(object):
        ...       @cached_method
        ...       def f(self, t, x=2):
        ...           print x + 1
        sage: a = Foo()
        sage: w = walltime()
        sage: a.f(1, 2); a.f(t = 1, x = 2); a.f(1)
        3
        sage: walltime(w) < 2
        True
    """
    def __init__(self, f):
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
            sage: a.f()
            4
            sage: hasattr(a, '_cache__f')
            True
        """
        self._cache_name = '_cache__' + f.__name__
        self._cachedfunc = CachedFunction(f, classmethod=True)

    def _instance_call(self, inst, *args, **kwds):
        """
        EXAMPLES::

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
            sage: a.f()  # indirect doctest
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
        cache = self._get_instance_cache(inst)
        key = self._get_instance_key(inst, *args, **kwds)
        if cache.has_key(key):
            return cache[key]
        else:
            cache[key] = self._cachedfunc.f(inst, *args, **kwds)
            return cache[key]

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
            sage: a.f.get_cache()  # indirect doctest
            {((), ()): 4}
        """
        return inst.__dict__.setdefault(self._cache_name, {})

    def _get_instance_key(self, inst, *args, **kwds):
        """
        Returns the key in the cache to be used when args
        and kwds are passed in as parameters with the given instance.

        TESTS::

            sage: class Foo:
            ...       @cached_method
            ...       def f(self, y):
            ...           return y
            sage: a = Foo()
            sage: a.f.get_key(37)  # indirect doctest
            ((37,), ())
            sage: a.f.get_key(y=5)
            ((5,), ())
        """
        return self._cachedfunc.get_key(*args, **kwds)

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
            sage: a.f is a.f
            False
        """
        return CachedMethodCaller(self, inst)

        # Note: a simpler approach to this would be
        # def caller(*args, **kwds):
        #     return self._instance_call(inst, *args, **kwds)
        # return caller
        # The disadvantage to this is that it does not provide
        # is_in_cache(), set_cache(), clear_cache(), ... methods.


cached_method = CachedMethod

class CachedInParentMethod(CachedMethod):
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
            ...       @cached_in_parent_method
            ...       def f(self):
            ...           return self._x^2
            ...
            sage: a = Foo(2)
            sage: a.f()
            4
            sage: hasattr(a.parent(), '_cache__element_f')
            True
        """
        self._cache_name = '_cache__' + 'element_' + f.__name__
        self._cachedfunc = CachedFunction(f, classmethod=True)

    def _get_instance_key(self, inst, *args, **kwds):
        """
        Returns the key used to lookup in the cache dictionary.
        This key includes the specific element instance.

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
            sage: a.f.get_key()   # indirect doctest
            ((2,), ())
            sage: a = Foo(2)
            sage: a.f.get_key(1,3)
            ((2, 1, 3), ())
            sage: a.f.get_key(1,3, bla=4)
            ((2, 1, 3), (('bla', 4),))
        """
        return self._cachedfunc.get_key(inst, *args, **kwds)

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
            ...       def __repr__(self):
            ...           return str(self._x)
            ...       @cached_in_parent_method
            ...       def f(self):
            ...           return self._x^2
            ...
            sage: a = Foo(2)
            sage: a.f()
            4
            sage: a.f.get_cache()   # indirect doctest
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
        return inst.parent().__dict__.setdefault(self._cache_name, {})


cached_in_parent_method = CachedInParentMethod

class FileCache:
    """
    FileCache is a dictionary-like class which stores keys and
    values on disk.  The keys take the form of a tuple (A,K)

        * A is a tuple of objects t where each t is an exact
          object which is uniquely identified by a short string.

        * K is a tuple of tuples (s,v) where s is a valid
          variable name and v is an exact object which is uniquely
          identified by a short string with letters [a-zA-Z0-9-._]

    The primary use case is the DiskCachedFunction.  If
    memory_cache == True, we maintain a cache of objects seen
    during this session in memory -- but we don't load them from
    disk until necessary.  The keys and values are stored in a
    pair of files:

        prefix-argstring.key.sobj contains the key only
        prefix-argstring.sobj contains the tuple (key,val)

    where self[key] == val.

    NOTE: We assume that each FileCache lives in its own directory.
    Use EXTREME caution if you wish to break that assumption.
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
        if len(dir) == 0 or dir[-1] != '/':
            dir += '/'
        self._dir = dir
        if not os.path.exists(dir):
            os.mkdir(dir)

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
        Returns a list of tuples (k,v) where self[k] = v.

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
        Returns a list of values where self[k] = v.

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
        Returns a list of values where self[k] = v for some k.

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
        Returns a list of keys k where self[k] is defined.

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
        Returns True if self[key] is defined and False otherwise.

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
        Returns the value set by self[key] = val, in this session
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
        """
        from sage.structure.sage_object import load

        cache = self._cache
        if cache is not None:
            if cache.has_key(key):
                return cache[key]

        f = self._filename(key) + '.sobj'
        k,v = load(f)
        if k != key:
            raise RuntimeError, "cache corrupted"

        if cache is not None:
            cache[key] = v
        return v

    def __setitem__(self, key, value):
        """
        Sets self[key] = value and stores both key and value on
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
    Decorator for DiskCachedFunction.

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
