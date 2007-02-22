"""
Field of Arbitrary Precision Real Number Intervals
"""

from sage.rings.real_mpfi import (RealIntervalField, RealIntervalFieldElement as RealIntervalFieldElementClass,
                                  create_RealIntervalFieldElement as RealIntervalFieldElement)

def is_RealIntervalField(x):
    return isinstance(x, RealIntervalField)

def is_RealIntervalFieldElement(x):
    return isinstance(x, RealIntervalFieldElementClass)

def __reduce__RealIntervalField(prec, sci_not):
    return RealIntervalField(prec, sci_not)

def __reduce__RealIntervalFieldElement(parent, x, base=10):
    return RealIntervalFieldElementClass(parent, x, base=base)
