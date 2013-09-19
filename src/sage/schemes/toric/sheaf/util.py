"""
Utility classes
"""

from copy import copy


def ensure_immutable(x):
    """
    EXAMPLES::

        sage: from sage.schemes.toric.sheaf.util import ensure_immutable
        sage: ensure_immutable(123)
        123
        sage: ensure_immutable([1, 2, 3])
        (1, 2, 3)
        sage: v = vector([1, 2, 3])
        sage: ensure_immutable(v)
        (1, 2, 3)
        sage: ensure_immutable(v) == v
        True
        sage: ensure_immutable(v) is v
        False
    """
    if isinstance(x, list):
        return tuple(x)
    try:
        is_immutable = x.is_immutable()
        if not is_immutable:
            x = copy(x)
            x.set_immutable()
        return x
    except AttributeError:
        return x
        

class CachedCaller(object):
    
    def __init__(self, original_caller):
        self.original_caller = original_caller
        self.cache = {}

    def __call__(self, *args, **kwds):
        args = tuple(map(ensure_immutable, args))
        kwds = dict((key, ensure_immutable(value)) for key, value in kwds.iteritems())
        cache = self.cache
        key = (args, tuple(kwds))
        try:
            return cache[key]
        except TypeError:
            print 'not hashable'
            return self.original_caller(*args, **kwds)
        except KeyError:
            result = self.original_caller(*args, **kwds)
            cache[key] = result
            return result


class CachedMethodImmutable(object):

    def __init__(self, method, name=None):
        """
        
        EXAMPLES::
        
            sage: from sage.schemes.toric.sheaf.util import cached_method_immutable
            sage: class Foo(object):
            ....:     @cached_method_immutable
            ....:     def bar(self, x):
            ....:         return x
            sage: foo = Foo()
            sage: foo.bar(123)
            123
            sage: foo.bar(vector([1, 2, 3]))
            (1, 2, 3)
        """
        self.method = method
        self.name = method.__name__ if name is None else name
        self.__doc__ = method.__doc__

    def __get__(self, instance, cls, *args):
        if instance is None:
            return self
        instance_caller = self.method.__get__(instance, cls)
        cached_caller = CachedCaller(instance_caller)
        setattr(instance, self.name, cached_caller)
        return cached_caller


cached_method_immutable = CachedMethodImmutable
