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
        sage: from sage.rings.number_field.number_field_base import is_NumberField
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

    def signature(self):
        """
        Return (r1, r2), where r1 and r2 are the number of real embeddings
        and pairs of complex embeddings of this field, respectively.
        """
        raise NotImplementedError

    def degree(self):
        """
        Return the degree of this number field.
        """
        raise NotImplementedError

    def discriminant(self):
        """
        Return the discriminant of this number field.
        """
        raise NotImplementedError

    def minkowski_bound(self):
        r"""
        Return the Minkowski bound associated to this number field.

        EXAMPLES:
        The Minkowski bound for $\QQ[i]$ tells us that the class
        number is 1:
            sage: K = QQ[I]
            sage: B = K.minkowski_bound(); B
            4/pi
            sage: B.n()
            1.27323954473516

        We compute the Minkowski bound for $\QQ[\sqrt[3]{2}]$:
            sage: K = QQ[2^(1/3)]
            sage: B = K.minkowski_bound(); B
            16*sqrt(3)/(3*pi)
            sage: B.n()
            2.94042077558289
            sage: int(B)
            2

        We compute the Minkowski bound for $\QQ[\sqrt{10}]$, which
        has class number $2$:
            sage: K = QQ[sqrt(10)]
            sage: B = K.minkowski_bound(); B
            sqrt(10)
            sage: int(B)
            3
            sage: K.class_number()
            2

        The bound of course also works for the rational numbers:
            sage: QQ.minkowski_bound()
            1
        """
        _, s = self.signature()
        n = self.degree()
        d = self.discriminant().abs().sqrt()
        from sage.functions.constants import pi
        if s > 0:
            return d * (4/pi)**s * n.factorial() / (n**n)
        else:
            return d * n.factorial() / (n**n)

