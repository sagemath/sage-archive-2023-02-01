"""
PYREX: sage.ext.rational
"""

from sage.ext.rational import *

def make_rational(s):
    r = Rational()
    r._reduce_set(s)
    return r
