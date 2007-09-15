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

    def ring_of_integers(self):
        r"""
        Synomym for \code{self.maximal_order()}.
        """
        return self.maximal_order()

    def maximal_order(self):
        raise NotImplementedError

    def is_absolute(self):
        """
        Return True if self is viewed as a single extension over Q.

        EXAMPLES:
            sage: K.<a> = NumberField(x^3+2)
            sage: K.is_absolute()
            True
            sage: y = polygen(K)
            sage: L.<b> = NumberField(y^2+1)
            sage: L.is_absolute()
            False
            sage: QQ.is_absolute()
            True
        """
        return False

