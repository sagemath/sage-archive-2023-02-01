"""
PYREX: sage.ext.mpfr
"""

from sage.ext.mpfr import RealField, RealNumber, mpfr_prec_min, mpfr_prec_max

def is_RealField(x):
    return isinstance(x, RealField)

def is_RealNumber(x):
    return isinstance(x, RealNumber)

def __reduce__RealField(prec, sci_not, rnd):
    return RealField(prec, sci_not, rnd)

def __reduce__RealNumber(parent, x, base=10):
    return RealNumber(parent, x, base=base)

