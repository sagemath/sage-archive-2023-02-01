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

    def ring_of_integers(self, *args, **kwds):
        r"""
        Synomym for \code{self.maximal_order(...)}.

        EXAMPLES:
            sage: K.<a> = NumberField(x^2 + 1)
            sage: K.ring_of_integers()
            Maximal Order in Number Field in a with defining polynomial x^2 + 1
        """
        return self.maximal_order()

    def OK(self, *args, **kwds):
        r"""
        Synomym for \code{self.maximal_order(...)}.

        EXAMPLES:
            sage: NumberField(x^3 - 2,'a').OK()
            Maximal Order in Number Field in a with defining polynomial x^3 - 2
        """
        return self.maximal_order(*args, **kwds)

    def maximal_order(self):
        """
        Return the maximal order, i.e., the ring of integers of this
        number field.

        EXAMPLES:
            sage: NumberField(x^3 - 2,'b').maximal_order()
            Maximal Order in Number Field in b with defining polynomial x^3 - 2
        """
        raise NotImplementedError

    def is_finite(self):
        """
        Return False since number fields are not finite.

            sage: z = polygen(QQ)
            sage: K.<theta, beta> = NumberField([z^3 - 3, z^2 + 1])
            sage: K.is_finite()
            False
            sage: K.order()
            +Infinity
        """
        return False


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

