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
        sage: is_NumberField(QQ)
        True
    """
    return isinstance(x, NumberField)

from sage.rings.ring cimport Field

cdef class NumberField(Field):
    pass
