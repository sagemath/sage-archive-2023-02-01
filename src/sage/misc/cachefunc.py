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
class CachedFunction(object):
    def __init__(self, f):
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

        EXAMPELES:
            sage: g = CachedFunction(number_of_partitions)
            sage: g.__name__
            'number_of_partitions'
            sage: 'partitions' in g.__doc__
            True
            sage: g(5)
            7
            sage: g.cache
            {((5,), ()): 7}

        """
        self.f = f
        self.cache = {}
        self.__doc__ = f.func_doc
        self.__name__ = f.func_name
        self.instance = None

    def __call__(self, *args, **kwds):
        """
        EXAMPLES:
            sage: g = CachedFunction(number_of_partitions)
            sage: a = g(5)
            sage: g.cache
            {((5,), ()): 7}
            sage: a = g(10^5)
            sage: a == number_of_partitions(10^5)
            True
        """
        k = (args, tuple(kwds))
        if self.cache.has_key(k):
            return self.cache[k]

        #Handle the case where self.f is a method of a class
        if self.instance is not None:
            w = self.f(self.instance, *args, **kwds)
        else:
            w = self.f(*args, **kwds)
        self.cache[k] = w
        return w

    def __repr__(self):
        """
        EXAMPLES:
            sage: g = CachedFunction(number_of_partitions)
            sage: g
            Cached version of <function number_of_partitions at 0x...>
        """
        return "Cached version of %s"%self.f

    def __get__(self, inst, cls=None):
        """
        This is needed to allow CachedFunction to decorate
        methods.

        EXAMPLES:
            sage: class Foo:
            ...       @CachedFunction
            ...       def f(self, x):
            ...           return x^2
            ...
            sage: a = Foo()
            sage: a.f(2)
            4
        """
        self.instance = inst
        return self
