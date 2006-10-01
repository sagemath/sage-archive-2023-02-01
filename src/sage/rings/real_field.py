"""
PYREX: sage.rings.mpfr
"""

from sage.rings.mpfr import (RealField, RR, RealNumber as RealNumberClass,
                             mpfr_prec_min, mpfr_prec_max,
                             create_RealNumber as RealNumber)

def is_RealField(x):
    return isinstance(x, RealField)

def is_RealNumber(x):
    return isinstance(x, RealNumberClass)

def __reduce__RealField(prec, sci_not, rnd):
    return RealField(prec, sci_not, rnd)

def __reduce__RealNumber(parent, x, base=10):
    return RealNumberClass(parent, x, base=base)





