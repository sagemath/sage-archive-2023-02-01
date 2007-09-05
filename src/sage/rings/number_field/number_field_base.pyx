"""
Base class for all number fields.


TESTS:
    sage: k = NumberField(x^2 + 1, 'i'); k == loads(dumps(k))
    True
"""

def is_NumberField(x):
    """
    Return True if x is of number field type.

    EXAMPLES:
        sage: is_NumberField(NumberField(x^2+1,'a'))
        True
        sage: is_NumberField(QuadraticField(-97,'theta'))
        True
        sage: is_NumberField(CyclotomicField(97))
        True

    Note that the rational numbers QQ are a number field.
        sage: is_NumberField(QQ)
        True
        sage: is_NumberField(ZZ)
        False
    """
    return isinstance(x, NumberField)

from sage.rings.ring cimport Field

cdef class NumberField(Field):
    pass
