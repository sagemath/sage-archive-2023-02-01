"""
Field of Double-Precision Real Numbers

PYREX: sage.rings.real_double
"""

include '../ext/cdefs.pxi'
include '../ext/interrupt.pxi'
include '../gsl/gsl.pxi'

import operator

from sage.misc.sage_eval import sage_eval

cimport sage.structure.element
import  sage.structure.element

cimport sage.rings.ring
import  sage.rings.ring

import sage.misc.functional
#import real_number


cdef class RealDoubleField_class(sage.rings.ring.Field):
    """
    The field of real double precision numbers.

    ALGORITHM: Arithmetic is done through pyrex.
    """

    def __cmp__(self, other):
        """
        Returns True if and only if other is the unique real double field.

        EXAMPLES:
            sage: RR == RDF
            False
            sage: RDF == RealDoubleField     # RDF is the shorthand
            True
        """
        if other is RealDoubleField:
            return 0
        return -1

    def __repr__(self):
        """
        Print out this real double field.

        EXAMPLES:
            sage: RealDoubleField
            Real Double Field
            sage: RDF
            Real Double Field
        """
        return "Real Double Field"

    def __call__(self, x):
        """
        Create a real double using x.

        EXAMPLES:
            sage: RDF(1)
            1.0
            sage: RDF(2/3)
            0.666666666667

        A TypeError is raised if the coercion doesn't make sense:
            sage: RDF(QQ['x'].0)
            Traceback (most recent call last):
            ...
            TypeError: cannot coerce nonconstant polynomial to float

        One can convert back and forth between double precision real
        numbers and higher-precision ones, though of course there may
        be loss of precision:
            sage: a = RealField(200)(2).sqrt(); a
            1.4142135623730950488016887242096980785696718753769480731766796
            sage: b = RDF(a); b
            1.41421356237
            sage: a.parent()(b)
            1.4142135623700000000000000000000000000000000000000000000000002
        """
        return RealDoubleElement(x)

    def gen(self, n=0):
        """
        Return the generator of the real double field.
        EXAMPLES:
            sage: RDF.0
            1.0
            sage: RDF.gens()
            (1.0,)
        """
        if n != 0:
            raise ValueError, "only 1 generator"
        return RealDoubleElement(1)

    def ngens(self):
        return 1

    def is_atomic_repr(self):
        """
        Returns True, to signify that elements of this field print
        without sums, so parenthesis aren't required, e.g., in
        coefficients of polynomials.

        EXAMPLES:
            sage: RealDoubleField.is_atomic_repr()
            True
        """
        return True

    def is_finite(self):
        """
        Returns False, since the field of real numbers is not finite.
        Technical note:  There exists an upper bound on the double representation.

        EXAMPLES:
            sage: RealDoubleField.is_finite()
            False
        """
        return False

    def characteristic(self):
        """
        Returns 0, since the field of real numbers has characteristic 0.

        EXAMPLES:
            sage: RealDoubleField.characteristic()
            0
        """
        return 0

    def name(self):
        return "RealDoubleField"

    def __hash__(self):
        return 1455926870 #return hash(self.name())

    def pi(self):
        """
        Returns pi to double-precision.

        EXAMPLES:
            sage: R = RealField(100)
            sage: RDF.pi()
            3.14159265359
            sage: RDF.pi().sqrt()/2
            0.886226925453
        """
        return self(M_PI)


    def euler_constant(self):
        """
        Returns Euler's gamma constant to double precision

        EXAMPLES:
            sage: RDF.euler_constant()
            0.577215664902
        """
        return self(M_EULER)

    def log2(self):
        """
        Returns log(2) to the precision of this field.

        EXAMPLES:
            sage: RDF.log2()
            0.69314718056
            sage: RDF(2).log()
            0.69314718056
        """
        return self(M_LN2)

    def factorial(self, int n):
        """
        Return the factorial of the integer n as a real number.
        EXAMPLES:
            sage: RDF.factorial(100)
            9.33262154439e+157
        """
        if n < 0:
            raise ArithmeticError, "n must be nonnegative"
        return self(gsl_sf_fact(n))

    def zeta(self, n=2):
        """
        Return an $n$-th root of unity in the real field,
        if one exists, or raise a ValueError otherwise.

        EXAMPLES:
            sage: RDF.zeta()
            -1.0
            sage: RDF.zeta(1)
            1.0
            sage: RDF.zeta(5)
            Traceback (most recent call last):
            ...
            ValueError: No 5th root of unity in self
        """
        if n == 1:
            return self(1)
        elif n == 2:
            return self(-1)
        raise ValueError, "No %sth root of unity in self"%n





cdef class RealDoubleElement(sage.structure.element.FieldElement):
    cdef double _value
    def __init__(self, x):
        self._value = float(x)

    def real(self):
        """
        Returns itself -- we're already real.

        EXAMPLES:
            sage: a = RDF(3)
            sage: a.real()
            3.0
        """
        return self

    def imag(self):
        """
        Returns the imaginary part of this number.  (hint: it's zero.)

        EXAMPLES:
            sage: a = RDF(3)
            sage: a.imag()
            0.0
        """
        return RealDoubleElement(0)

    def __complex__(self):
        """
        EXAMPLES:
            sage: a = 2303
            sage: RDF(a)
            2303.0
            sage: complex(RDF(a))
            (2303+0j)
        """
        return complex(self._value,0)

    def parent(self):
        """
        Return the real double field, which is the parent of self.

        EXAMPLES:
            sage: a = RDF(2.3)
            sage: a.parent()
            Real Double Field
            sage: parent(a)
            Real Double Field
        """
        return RDF

    def __repr__(self):
        """
        Return print version of self.

        EXAMPLES:
            sage: a = RDF(2); a
            2.0
            sage: a^2
            4.0
        """
        return self.str()

    def _latex_(self):
        return self.str()

    def __hash__(self):
        return hash(self.str())

    def _im_gens_(self, codomain, im_gens):
        return codomain(self) # since 1 |--> 1

    def str(self, no_sci=None):
        if gsl_isnan(self._value):
            return "nan"
        else:
            v = gsl_isinf(self._value)
            if v == 1:
                return "inf"
            elif v == -1:
                return "-inf"
        if no_sci is not None and not no_sci:
            return "%e"%self._value
        else:
            return str(self._value)

    def copy(self):
        cdef RealDoubleElement z
        z = RealDoubleElement(self._value)
        return z

    def integer_part(self):
        """
        If in decimal this number is written n.defg, returns n.
        """
        return sage.rings.all.integer.Integer(int(self._value))


    ########################
    #   Basic Arithmetic
    ########################
    def _add_(RealDoubleElement self, RealDoubleElement other):
        """
        Add two real numbers with the same parent.

        EXAMPLES:
            sage: R = RealDoubleField
            sage: R(-1.5) + R(2.5)
            1.0
        """
        return RealDoubleElement(self._value + other._value)

    def __invert__(self):
        """
        Compute the multiplicative inverse of self.

        EXAMPLES:
            sage: a = RDF(-1.5)*RDF(2.5)
            sage: a.__invert__()
            -0.266666666667
        """
        return RealDoubleElement(1/self._value)

    def _sub_(RealDoubleElement self, RealDoubleElement other):
        """
        Subtract two real numbers with the same parent.

        EXAMPLES:
            sage: R = RealDoubleField
            sage: R(-1.5) - R(2.5)
            -4.0
        """
        return RealDoubleElement(self._value - other._value)

    def _mul_(RealDoubleElement self, RealDoubleElement other):
        """
        Multiply two real numbers with the same parent.

        EXAMPLES:
            sage: R = RealDoubleField
            sage: R(-1.5) * R(2.5)
            -3.75
        """
        return RealDoubleElement(self._value * other._value)

    def __div_(RealDoubleElement self, RealDoubleElement other):
        if not other:
            raise ZeroDivisionError, "RealDoubleElement division by zero"
        return RealDoubleElement(self._value / other._value)

    def __div__(x, y):
        """
        EXAMPLES:
            sage: R = RealDoubleField
            sage: R(-1.5) / R(2.5)
            -0.6
        """
        if isinstance(x, RealDoubleElement) and isinstance(y, RealDoubleElement):
            return x.__div_(y)
        return sage.rings.coerce.bin_op(x, y, operator.div)

    def __neg__(self):
        return RealDoubleElement(-self._value)

    def __pos__(self):
        return self

    def __abs__(self):
        return self.abs()

    cdef RealDoubleElement abs(RealDoubleElement self):
        return RealDoubleElement(abs(self._value))

    def __lshift__(x, y):
        """
        LShifting a double is not supported; nor is lshifting a RealDoubleElement.
        """
        raise TypeError, "unsupported operand type(s) for <<: '%s' and '%s'"%(typeof(self), typeof(n))

    def __rshift__(x, y):
        """
        RShifting a double is not supported; nor is rshifting a RealDoubleElement.
        """
        raise TypeError, "unsupported operand type(s) for >>: '%s' and '%s'"%(typeof(self), typeof(n))

    def multiplicative_order(self):
        if self == 1:
            return 1
        elif self == -1:
            return -1
        return sage.rings.infinity.infinity

    def sign(self):
        """
        Returns -1,0, or 1 if self is negative, zero, or positive; respectively.

        Examples:
            sage: RDF(-1.5).sign()
            -1
            sage: RDF(0).sign()
            0
            sage: RDF(2.5).sign()
            1
        """
        if not self._value:
            return 0
        if self._value > 0:
            return 1
        return -1


    ###################
    # Rounding etc
    ###################

    def round(self):
        """
        Given real number x, rounds up if fractional part is greater than .5,
        rounds down if fractional part is lesser than .5.
        EXAMPLES:
            sage: RDF(0.49).round()
            0.0
            sage: RDF(0.51).round()
            1.0
        """
        return RealDoubleElement(round(self._value))

    def floor(self):
        """
        Returns the floor of this number

        EXAMPLES:
            sage: RDF(2.99).floor()
            2
            sage: RDF(2.00).floor()
            2
            sage: RDF(-5/2).floor()
            -3
        """
        return sage.misc.functional.floor(self._value)

    def ceil(self):
        """
        Returns the ceiling of this number

        OUTPUT:
            integer

        EXAMPLES:
            sage: RDF(2.99).ceil()
            3
            sage: RDF(2.00).ceil()
            2
            sage: RDF(-5/2).ceil()
            -2
        """
        return sage.misc.functional.ceil(self._value)

    def ceiling(self):
        return self.ceil()

    def trunc(self):
        """
        Truncates this number (returns integer part).

        EXAMPLES:
            sage: RDF(2.99).trunc()
            2.0
            sage: RDF(-2.00).trunc()
            -2.0
            sage: RDF(0.00).trunc()
            0.0
        """
        return RealDoubleElement(int(self._value))

    def frac(self):
        """
        frac returns a real number > -1 and < 1. it satisfies the
        relation:
            x = x.trunc() + x.frac()

        EXAMPLES:
            sage: RDF(2.99).frac()
            0.99
            sage: RDF(2.50).frac()
            0.5
            sage: RDF(-2.79).frac()
            -0.79
        """
        return RealDoubleElement(self._value - int(self._value))

    ###########################################
    # Conversions
    ###########################################

    def __float__(self):
        return self._value

    def __int__(self):
        """
        Returns integer truncation of this real number.
        """
        return int(self._value)

    def __long__(self):
        """
        Returns long integer truncation of this real number.
        """
        return long(self._value)

    def _complex_number_(self):
        return sage.rings.complex_field.ComplexField()(self)

    def _complex_double_(self):
         return sage.rings.complex_double.ComplexDoubleField(self)

    def _pari_(self):
        return sage.libs.pari.all.pari.new_with_bits_prec("%.15e"%self._value, 64)


    ###########################################
    # Comparisons: ==, !=, <, <=, >, >=
    ###########################################

    def is_NaN(self):
        return bool(gsl_isnan(self._value))

    cdef int cmp(RealDoubleElement self, RealDoubleElement x):
        if self._value < x._value:
            return -1
        elif self._value > x._value:
            return 1
        return 0

    def __cmp__(RealDoubleElement self, RealDoubleElement x):
        return self.cmp(x)

    def __richcmp__(RealDoubleElement self, x, int op):
        cdef int n
        if not isinstance(x, RealDoubleElement):
            try:
                x = RealDoubleElement(x)
            except TypeError:
                n = sage.rings.coerce.cmp(self, x)
            else:
                n = self.cmp(x)
        else:
            n = self.cmp(x)
        if op == 0:
            return bool(n < 0)
        elif op == 1:
            return bool(n <= 0)
        elif op == 2:
            return bool(n == 0)
        elif op == 3:
            return bool(n != 0)
        elif op == 4:
            return bool(n > 0)
        elif op == 5:
            return bool(n >= 0)


    ############################
    # Special Functions
    ############################

    def sqrt(self):
        """
        Return a square root of self.

        If self is negative a complex number is returned.

        If you use self.square_root() then a real number will always
        be returned (though it will be NaN if self is negative).

        EXAMPLES:
            sage: r = 4.0
            sage: r.sqrt()
            2.0000000000000000
            sage: r.sqrt()^2 == r
            True

            sage: r = 4344
            sage: r.sqrt()
            65.909028213136324
            sage: r.sqrt()^2 == r
            True

            sage: r = -2.0
            sage: r.sqrt()
            1.4142135623730951*I
            """
        if self >= 0:
            return self.square_root()
        return self._complex_double_().sqrt()


    def square_root(self):
        """
        Return a square root of self.  A real number will always be
        returned (though it will be NaN if self is negative).

        Use self.sqrt() to get a complex number if self is negative.

        EXAMPLES:
            sage: r = -2.0
            sage: r.square_root()
            NaN
            sage: r.sqrt()
            1.4142135623730951*I
        """
        return RealDoubleElement(sqrt(self._value))

    def cube_root(self):
        """
        Return the cubic root (defined over the real numbers) of self.

        EXAMPLES:
            sage: r = RDF(125.0); r.cube_root()
            5.0
            sage: r = RDF(-119.0)
            sage: r.cube_root()^3 - r         # output is random, depending on arch.
            0.0
        """
        return self.nth_root(3)


    def nth_root(self, int n):
        """
        Returns the $n^{th}$ root of self.
        EXAMPLES:
            sage: r = RDF(-125.0); r.nth_root(3)
            -5.0
            sage: r.nth_root(5)
            -2.6265278044
        """
        if n == 0:
            return RealDoubleElement(float('nan'))
        if self._value < 0 and GSL_IS_EVEN(n):
            pass #return self._complex_double_().pow(1.0/n)
        else:
            return RealDoubleElement(self.__nth_root(n))

    cdef double __nth_root(RealDoubleElement self, int n):
        cdef int m
        cdef double x
        cdef double x0
        cdef double dx
        cdef double dx0
        m  = n-1
        x  = ( m + self._value ) / n
        x0 = 0
        dx = abs(x - x0)
        dx0= dx + 1
        while dx < dx0:
            x0= x
            dx0 = dx
            x = ( m*x + self._value / gsl_pow_int(x,m) ) / n
            dx=abs(x - x0)
        return x

    def __pow(self, RealDoubleElement exponent):
        return RealDoubleElement(self._value**exponent._value)

    def __pow_int(self, int exponent):
        return RealDoubleElement(gsl_pow_int(self._value, exponent))

    def __pow__(self, exponent, modulus):
        """
        Compute self raised to the power of exponent, rounded in
        the direction specified by the parent of self.

        If the result is not a real number, self and the exponent are
        both coerced to complex numbers (with sufficient precision),
        then the exponentiation is computed in the complex numbers.
        Thus this function can return either a real or complex number.

        EXAMPLES:
            sage: a = RDF('1.23456')
            sage: a^20
            67.6462977039
            sage: a^a
            1.29711148178
        """
        cdef RealDoubleElement x
        if isinstance(self, RealDoubleElement):
            return self.__pow(RealDoubleElement(exponent))
        if isinstance(exponent, (int,Integer)):
            return self.__pow_int(int(exponent))
        elif not isinstance(exponent, RealDoubleElement):
            x = RealDoubleElement(exponent)
        else:
            x = exponent
        return self.__pow(x)


    def __log_(self, double log_of_base):
        if self._value < 2:
            if self._value == 0:
                return -1./0
            if self._value < 0:
                return 0./0
            return RealDoubleElement(gsl_sf_log_1plusx(self._value - 1) / log_of_base)
        return RealDoubleElement(gsl_sf_log(self._value) / log_of_base)

    def log(self, base='e'):
        """
        EXAMPLES:
            sage: RDF(2).log()
            0.69314718056
            sage: RDF(2).log(2)
            1.0
            sage: RDF(2).log(pi)
            0.605511561398
            sage: RDF(2).log(10)
            0.301029995664
            sage: RDF(2).log(1.5)
            1.70951129135
            sage: RDF(0).log()
            -inf
            sage: RDF(-1).log()
            nan
        """
        if base == 'e':
            return self.__log_(1)
        elif base == 'pi':
            return self.logpi()
        elif base == 2:
            return self.log2()
        elif base == 10:
            return self.log10()
        else:
            if isinstance(base, RealDoubleElement):
                return self.__log_(base.__log_(1))
            else:
                return self.__log_(gsl_sf_log(float(base)))

    def log2(self):
        """
        Returns log to the base 2 of self

        EXAMPLES:
            sage: r = RDF(16.0)
            sage: r.log2()
            4.0

            sage: r = RDF(31.9); r.log2()
            4.99548451888

        """
        return RealDoubleElement(gsl_sf_log(self._value) / M_LN2)


    def log10(self):
        """
        Returns log to the base 10 of self

        EXAMPLES:
            sage: r = RDF('16.0'); r.log10()
            1.20411998266
            sage: r.log() / log(10)
            1.20411998266
            sage: r = RDF('39.9'); r.log10()
            1.60097289569
        """
        return RealDoubleElement(gsl_sf_log(self._value) / M_LN10)

    def logpi(self):
        """
        Returns log to the base pi of self

        EXAMPLES:
            sage: r = RDF(16); r.logpi()
            2.42204624559
            sage: r.log() / log(pi)
            2.42204624559
            sage: r = RDF('39.9'); r.logpi()
            3.22030233461
        """
        return RealDoubleElement(gsl_sf_log(self._value) / M_LNPI)

    def exp(self):
        r"""
        Returns $e^\code{self}$

        EXAMPLES:
            sage: r = RDF(0.0)
            sage: r.exp()
            1.0

            sage: r = RDF('32.3')
            sage: a = r.exp(); a
            1.06588847275e+14
            sage: a.log()
            32.3

            sage: r = RDF('-32.3')
            sage: r.exp()
            9.3818445885e-15
        """
        return RealDoubleElement(gsl_sf_exp(self._value))

    def exp2(self):
        """
        Returns $2^\code{self}$

        EXAMPLES:
            sage: r = RDF(0.0)
            sage: r.exp2()
            1.0

            sage: r = RDF(32.0)
            sage: r.exp2()
            4294967296.0

            sage: r = RDF(-32.3)
            sage: r.exp2()
            1.89117248253e-10

        """
        return RealDoubleElement(gsl_sf_exp(self._value * M_LN2))

    def exp10(self):
        r"""
        Returns $10^\code{self}$

        EXAMPLES:
            sage: r = RDF(0.0)
            sage: r.exp10()
            1.0

            sage: r = RDF(32.0)
            sage: r.exp10()
            1e+32

            sage: r = RDF(-32.3)
            sage: r.exp10()
            5.01187233627e-33
        """
        return RealDoubleElement(gsl_sf_exp(self._value * M_LN10))

    def cos(self):
        """
        Returns the cosine of this number

        EXAMPLES:
            sage: t=RDF.pi()/2
            sage: t.cos()
            6.12323399574e-17
        """
        return RealDoubleElement(gsl_sf_cos(self._value))

    def sin(self):
        """
        Returns the sine of this number

        EXAMPLES:
            sage: RDF(2).sin()
            0.909297426826
        """
        return RealDoubleElement(gsl_sf_sin(self._value))

    def tan(self):
        """
        Returns the tangent of this number

        EXAMPLES:
            sage: q = RDF.pi()/3
            sage: q.tan()
            1.73205080757
            sage: q = RDF.pi()/6
            sage: q.tan()
            0.57735026919
        """
        return RealDoubleElement(tan(self._value))

    def sincos(self):
        """
        Returns a pair consisting of the sine and cosine.

        EXAMPLES:
            sage: t = RDF.pi()/6
            sage: t.sincos()
            (0.5, 0.866025403784)
        """
        return self.sin(), self.cos()

    def hypot(self, other):
        return RealDoubleElement(gsl_sf_hypot(self._value, float(other)))

    def acos(self):
        """
        Returns the inverse cosine of this number

        EXAMPLES:
            sage: q = RDF.pi()/3
            sage: i = q.cos()
            sage: i.acos() == q
            True
        """
        return RealDoubleElement(acos(self._value))

    def asin(self):
        """
        Returns the inverse sine of this number

        EXAMPLES:
            sage: q = RDF.pi()/5
            sage: i = q.sin()
            sage: i.asin() == q
            True
        """
        return RealDoubleElement(asin(self._value))

    def atan(self):
        """
        Returns the inverse tangent of this number

        EXAMPLES:
            sage: q = RDF.pi()/5
            sage: i = q.tan()
            sage: i.atan() == q
            True
        """
        return RealDoubleElement(atan(self._value))


    def cosh(self):
        """
        Returns the hyperbolic cosine of this number

        EXAMPLES:
            sage: q = RDF.pi()/12
            sage: q.cosh()
            1.0344656401
        """
        return RealDoubleElement(cosh(self._value))

    def sinh(self):
        """
        Returns the hyperbolic sine of this number

        EXAMPLES:
            sage: q = RDF.pi()/12
            sage: q.sinh()
            0.264800227602

        """
        return RealDoubleElement(sinh(self._value))

    def tanh(self):
        """
        Returns the hyperbolic tangent of this number

        EXAMPLES:
            sage: q = RDF.pi()/12
            sage: q.tanh()
            0.255977789246
        """
        return RealDoubleElement(tanh(self._value))

    def acosh(self):
        """
        Returns the hyperbolic inverse cosine of this number

        EXAMPLES:
            sage: q = RDF.pi()/2
            sage: i = q.cosh() ; i
            2.50917847866
            sage: i.acosh() == q
            True
        """
        return RealDoubleElement(gsl_acosh(self._value))

    def asinh(self):
        """
        Returns the hyperbolic inverse sine of this number

        EXAMPLES:
            sage: q = RDF.pi()/2
            sage: i = q.sinh() ; i
            2.30129890231
            sage: i.asinh() == q
            True
        """
        return RealDoubleElement(gsl_asinh(self._value))

    def atanh(self):
        """
        Returns the hyperbolic inverse tangent of this number

        EXAMPLES:
            sage: q = RDF.pi()/2
            sage: i = q.tanh() ; i
            0.917152335667
            sage: i.atanh() - q      # output is random, depending on arch.
            -4.4408920985e-16
        """
        return RealDoubleElement(gsl_atanh(self._value))

    def agm(self, other):
        """
        Return the arithmetic-geometric mean of self and other. The
        arithmetic-geometric mean is the common limit of the sequences
        $u_n$ and $v_n$, where $u_0$ is self, $v_0$ is other,
        $u_{n+1}$ is the arithmetic mean of $u_n$ and $v_n$, and
        $v_{n+1}$ is the geometric mean of u_n and v_n. If any operand
        is negative, the return value is \code{NaN}.
        """
        return RealDoubleElement(sage.rings.all.RR(self).agm(sage.rings.all.RR(other)))

    def erf(self):
        """
        Returns the value of the error function on self.

        EXAMPLES:
           sage: RDF(6).erf()
           1.0
        """
        return RealDoubleElement(gsl_sf_erf(self._value))

    def gamma(self):
        """
        The Euler gamma function. Return gamma of self.

        EXAMPLES:
           sage: RDF(6).gamma()
           120.0
           sage: RDF(1.5).gamma()
           0.886226925453
        """
        return RealDoubleElement(gsl_sf_gamma(self._value))

    def zeta(self):
        r"""
        Return the Riemann zeta function evaluated at this real number.

        \note{PARI is vastly more efficient at computing the Riemann zeta
        function.   See the example below for how to use it.}

        EXAMPLES:
            sage: RDF(2).zeta()
            1.64493406685
            sage: RDF.pi()^2/6
            1.64493406685
            sage: RDF(-2).zeta()       # slightly random-ish arch dependent output
            -2.37378795339e-18
            sage: RDF(1).zeta()
            inf
        """
        if self._value == 1:
            return RealDoubleElement(1)/RealDoubleElement(0)
        return RealDoubleElement(gsl_sf_zeta(self._value))

    def algdep(self, n):
        """
        Returns a polynomial of degree at most $n$ which is approximately
        satisfied by this number.  Note that the returned polynomial
        need not be irreducible, and indeed usually won't be if this number
        is a good approximation to an algebraic number of degree less than $n$.

        ALGORITHM: Uses the PARI C-library algdep command.

        EXAMPLE:
            sage: r = RDF(2).sqrt(); r
            1.41421356237
            sage: r.algdep(5)
            x^4 - 2*x^2
        """
        return sage.rings.arith.algdep(self,n)

    def algebraic_dependency(self, n):
        """
        Returns a polynomial of degree at most $n$ which is approximately
        satisfied by this number.  Note that the returned polynomial
        need not be irreducible, and indeed usually won't be if this number
        is a good approximation to an algebraic number of degree less than $n$.

        ALGORITHM: Uses the PARI C-library algdep command.

        EXAMPLE:
            sage: r = sqrt(RDF(2)); r
            1.41421356237
            sage: r.algdep(5)
            x^4 - 2*x^2
        """
        return sage.rings.arith.algdep(self,n)


#####################################################
# unique objects
#####################################################
RealDoubleField = RealDoubleField_class()
RDF = RealDoubleField





