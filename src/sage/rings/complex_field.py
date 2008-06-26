r"""
Field of Arbitrary Precision Complex Numbers

AUTHOR:
    -- William Stein (2006-01-26): complete rewrite
"""

#################################################################################
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import complex_number
import complex_double
import field
import integer_ring
import integer
import real_mpfr
import weakref
import ring
from sage.misc.sage_eval import sage_eval

from sage.structure.parent_gens import ParentWithGens

NumberFieldElement_quadratic = None
AlgebraicNumber_base = None
AlgebraicNumber = None
AlgebraicReal = None
def late_import():
    global NumberFieldElement_quadratic
    global AlgebraicNumber_base
    global AlgebraicNumber
    global AlgebraicReal
    if NumberFieldElement_quadratic is None:
        import sage.rings.number_field.number_field_element_quadratic as nfeq
        NumberFieldElement_quadratic = nfeq.NumberFieldElement_quadratic
        import sage.rings.qqbar
        AlgebraicNumber_base = sage.rings.qqbar.AlgebraicNumber_base
        AlgebraicNumber = sage.rings.qqbar.AlgebraicNumber
        AlgebraicReal = sage.rings.qqbar.AlgebraicReal

def is_ComplexField(x):
    return isinstance(x, ComplexField_class)

cache = {}
def ComplexField(prec=53, names=None):
    """
    Return the complex field with real and imaginary parts having prec
    *bits* of precision.

    EXAMPLES:
        sage: ComplexField()
        Complex Field with 53 bits of precision
        sage: ComplexField(100)
        Complex Field with 100 bits of precision
        sage: ComplexField(100).base_ring()
        Real Field with 100 bits of precision
        sage: i = ComplexField(200).gen()
        sage: i^2
        -1.0000000000000000000000000000000000000000000000000000000000
    """
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
    An approximation to the field of complex numbers using floating
    point numbers with any specified precision. Answers derived from
    calculations in this approximation may differ from what they would
    be if those calculations were performed in the true field of
    complex numbers. This is due to the rounding errors inherent to
    finite precision calculations.

    EXAMPLES:
        sage: C = ComplexField(); C
        Complex Field with 53 bits of precision
        sage: Q = RationalField()
        sage: C(1/3)
        0.333333333333333
        sage: C(1/3, 2)
        0.333333333333333 + 2.00000000000000*I
        sage: C(RQDF.pi())
        3.14159265358979
        sage: C(RQDF.log2(), RQDF.e())
        0.693147180559945 + 2.71828182845905*I

    We can also coerce rational numbers and integers into C, but
    coercing a polynomial will raise an exception.

        sage: Q = RationalField()
        sage: C(1/3)
        0.333333333333333
        sage: S = PolynomialRing(Q, 'x')
        sage: C(S.gen())
        Traceback (most recent call last):
        ...
        TypeError: unable to coerce to a ComplexNumber

    This illustrates precision.
        sage: CC = ComplexField(10); CC(1/3, 2/3)
        0.33 + 0.67*I
        sage: CC
        Complex Field with 10 bits of precision
        sage: CC = ComplexField(100); CC
        Complex Field with 100 bits of precision
        sage: z = CC(1/3, 2/3); z
        0.33333333333333333333333333333 + 0.66666666666666666666666666667*I

    We can load and save complex numbers and the complex field.
        sage: loads(z.dumps()) == z
        True
        sage: loads(CC.dumps()) == CC
        True
        sage: k = ComplexField(100)
        sage: loads(dumps(k)) == k
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
        'I'
        sage: CC == ComplexField(200)
        True
        sage: CC == ComplexField(53)
        False
        sage: CC == 1.1
        False

    """
    def __init__(self, prec=53):
        self._prec = int(prec)
        ParentWithGens.__init__(self, self._real_field(), ('I',), False)

    def __reduce__(self):
        return ComplexField, (self._prec, )

    def is_exact(self):
        return False

    def prec(self):
        return self._prec

    precision = prec


    # very useful to cache this.
    def _real_field(self):
        try:
            return self.__real_field
        except AttributeError:
            self.__real_field = real_mpfr.RealField(self._prec)
            return self.__real_field

    def __cmp__(self, other):
        if not isinstance(other, ComplexField_class):
            return cmp(type(self), type(other))
        return cmp(self._prec, other._prec)

    def __call__(self, x, im=None):
        """
        EXAMPLES:
            sage: CC(2)
            2.00000000000000
            sage: CC(CC.0)
            1.00000000000000*I
            sage: CC('1+I')
            1.00000000000000 + 1.00000000000000*I
            sage: CC(2,3)
            2.00000000000000 + 3.00000000000000*I
            sage: CC(QQ[I].gen())
            1.00000000000000*I
            sage: CC.gen() + QQ[I].gen()
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '+': 'Complex Field with 53 bits of precision' and 'Number Field in I with defining polynomial x^2 + 1'
        """
        if im is None:
            if isinstance(x, complex_number.ComplexNumber) and x.parent() is self:
                return x
            elif isinstance(x, complex_double.ComplexDoubleElement):
                return complex_number.ComplexNumber(self, x.real(), x.imag())
            elif isinstance(x, str):
                # TODO: this is probably not the best and most
                # efficient way to do this.  -- Martin Albrecht
                return complex_number.ComplexNumber(self,
                            sage_eval(x.replace(' ',''), locals={"I":self.gen(),"i":self.gen()}))

            late_import()
            if isinstance(x, NumberFieldElement_quadratic) and list(x.parent().polynomial()) == [1, 0, 1]:
                (re, im) = list(x)
                return complex_number.ComplexNumber(self, re, im)

            try:
                return x._complex_mpfr_field_( self )
            except AttributeError:
                pass
        return complex_number.ComplexNumber(self, x, im)

    def _coerce_impl(self, x):
        """
        Return the canonical coerce of x into this complex field, if it is defined,
        otherwise raise a TypeError.

        The rings that canonicaly coerce to the MPFS complex field are:
           * this MPFR complex field, or any other of higher precision
           * anything that canonically coerces to the mpfr real field with this prec

        EXAMPLES:
        sage: ComplexField(200)(1) + RealField(90)(1)
        2.0000000000000000000000000
        """
        try:
            K = x.parent()
            if is_ComplexField(K) and K._prec >= self._prec:
                return self(x)
        except AttributeError:
            pass
        late_import()
        if isinstance(x, AlgebraicNumber_base):
            return self(x)
        return self._coerce_try(x, self._real_field())

    def _repr_(self):
        return "Complex Field with %s bits of precision"%self._prec

    def _latex_(self):
        return "\\C"

    def characteristic(self):
        return integer.Integer(0)

    def gen(self, n=0):
        if n != 0:
            raise IndexError, "n must be 0"
        return complex_number.ComplexNumber(self, 0, 1)

    def is_field(self):
        """
        Return True, since the complex numbers are a field.

        EXAMPLES:
            sage: CC.is_field()
            True
        """
        return True

    def is_finite(self):
        """
        Return False, since the complex numbers are infinite.

        EXAMPLES:
            sage: CC.is_finite()
            False
        """
        return False

    def construction(self):
        """
        Returns the functorial construction of self, namely,
        algebraic closure of the real field with the same precision.

        EXAMPLES:
            sage: c, S = CC.construction(); S
            Real Field with 53 bits of precision
            sage: CC == c(S)
            True
        """
        from sage.categories.pushout import AlgebraicClosureFunctor
        return (AlgebraicClosureFunctor(), self._real_field())


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
        x._set_multiplicative_order( n )
        return x

    def scientific_notation(self, status=None):
        return self._real_field().scientific_notation(status)
