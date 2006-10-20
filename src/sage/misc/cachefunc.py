"""
Cached Functions

AUTHOR:
    -- William Stein (inspired by conversation with Justin Walker).
"""

########################################################################
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
########################################################################

class CachedFunction:
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
    """
    def __init__(self, f):
        self.f = f
        self.cache = {}

    def __call__(self, *args, **kwds):
        k = (args, tuple(kwds))
        if self.cache.has_key(k):
            return self.cache[k]
        w = self.f(*args, **kwds)
        self.cache[k] = w
        return w

    def __repr__(self):
        return "Cached version of %s"%self.f
