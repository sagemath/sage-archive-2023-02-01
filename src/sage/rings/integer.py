"""
PYREX: sage.ext.integer
"""

from sage.ext.integer import *

def make_integer(s):
    r = Integer()
    r._reduce_set(s)
    return r
