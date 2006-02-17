r"""
Field $\C$ of Complex Numbers

AUTHOR:
    -- William Stein (2006-01-26): complete rewrite
"""

#################################################################################
#       Copyright (C) 2006 William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import complex_number
import field
import real_field
import integer_ring
import integer
import weakref

def is_ComplexField(x):
    return isinstance(x, ComplexField_class)

cache = {}
def ComplexField(prec=53):
    global cache
    if cache.has_key(prec):
        X = cache[prec]
        C = X()
        if not C is None:
            return C
    C = ComplexField_class(prec)
    cache[prec] = weakref.ref(C)
    return C


class ComplexField_class(field.Field):
    """
    The field of complex numbers.

    EXAMPLES:
        sage: C = ComplexField(); C
        Complex Field with 53 bits of precision
        sage: Q = RationalField()
        sage: C(1/3)
        0.33333333333333331
        sage: C(1/3, 2)
        0.33333333333333331 + 2.0000000000000000*I

    Note that the second argument is the number of *bits* of precision,
    not the number of digits of precision:
        sage: C(1/3, 2)
        0.33333333333333331 + 2.0000000000000000*I

    We can also coerce rational numbers and integers into C, but
    coercing a polynomial in raising an exception.

        sage: Q = RationalField()
        sage: C(1/3)
        0.33333333333333331
        sage: S = PolynomialRing(Q)
        sage: C(S.gen())
        Traceback (most recent call last):
        ...
        TypeError: unable to coerce (x,0) to a ComplexNumber

    This illustrates precision.
        sage: CC = ComplexField(10); CC(1/3, 2/3)
        0.33350 + 0.66699*I
        sage: CC
        Complex Field with 10 bits of precision
        sage: CC = ComplexField(100); CC
        Complex Field with 100 bits of precision
        sage: z = CC(1/3, 2/3); z
        0.33333333333333333333333333333346 + 0.66666666666666666666666666666693*I

    We can load and save complex numbers and the complex field.
        sage: loads(z.dumps()) == z
        True
        sage: loads(CC.dumps()) == CC
        True

    This illustrates basic properties of a complex field.
        sage: CC = ComplexField(200)
        sage: CC.is_field()
        True
        sage: CC.characteristic()
        0
        sage: CC.precision()
        200
        sage: CC.variable_name()
        'x'
        sage: CC == ComplexField(200)
        True
        sage: CC == ComplexField(53)
        False
        sage: CC == 1.1
        False

    """
    def __init__(self, prec=53):
        self.__prec = int(prec)

    def prec(self):
        return self.__prec

    precision = prec


    # very useful to cache this.
    def _real_field(self):
        try:
            return self.__real_field
        except AttributeError:
            self.__real_field = real_field.RealField(self.__prec)
            return self.__real_field

    def __cmp__(self, other):
        if not isinstance(other, ComplexField_class):
            return -1
        return cmp(self.__prec, other.__prec)

    def __call__(self, x, im=None):
        if isinstance(x, complex_number.ComplexNumber) and x.parent() == self:
            return x
        return complex_number.ComplexNumber(self, x, im)

    def _repr_(self):
        return "Complex Field with %s bits of precision"%self.__prec

    def _latex_(self):
        return "\\C"

    def characteristic(self):
        return integer.Integer(0)

    def gen(self, n=0):
        if n != 0:
            raise IndexError, "n (=%s) must be 0"%n
        return complex_number.ComplexNumber(self, 0, 1)

    def is_field(self):
        return True

    def pi(self):
        return self(self._real_field().pi())

    def ngens(self):
        return 1

    def zeta(self, n=2):
        """
        Return a primitive $n$-th root of unity.

        INPUT:
            n -- an integer (default: 2)

        OUTPUT:
            a complex n-th root of unity.
        """
        from integer import Integer
        n = Integer(n)
        if n == 1:
            x = self(1)
        elif n == 2:
            x = self(-1)
        elif n >= 3:
            # Use De Moivre
            # e^(2*pi*i/n) = cos(2pi/n) + i *sin(2pi/n)
            RR = self._real_field()
            pi = RR.pi()
            z = 2*pi/n
            x = complex_number.ComplexNumber(self, z.cos(), z.sin())
        x._multiplicative_order = n
        return x

    def scientific_notation(self, status=None):
        return self._real_field().scientific_notation(status)

CC = ComplexField()

