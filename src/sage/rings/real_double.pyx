r"""
Double Precision Real Numbers

EXAMPLES:

We create the real double vector space of dimension $3$:
    sage: V = RDF^3; V
    Vector space of dimension 3 over Real Double Field

Notice that this space is unique.
    sage: V is RDF^3
    True
    sage: V is FreeModule(RDF, 3)
    True
    sage: V is VectorSpace(RDF, 3)
    True

Also, you can instantly create a space of large dimension.
    sage: V = RDF^10000
"""

include '../ext/cdefs.pxi'
include '../ext/stdsage.pxi'
include '../ext/random.pxi'
include '../ext/interrupt.pxi'
include '../gsl/gsl.pxi'

gsl_set_error_handler_off()

import math, operator

cimport sage.libs.pari.gen
import sage.libs.pari.gen

from sage.misc.sage_eval import sage_eval

import sage.rings.complex_double
import sage.rings.complex_field

import sage.rings.integer
import sage.rings.rational

from sage.rings.integer cimport Integer

def is_RealDoubleField(x):
    return bool(PY_TYPE_CHECK(x, RealDoubleField_class))

cdef class RealDoubleField_class(Field):
    """
    The field of real double precision numbers.

    EXAMPLES:
        sage: RR == RDF
        False
        sage: RDF == RealDoubleField()    # RDF is the shorthand
        True
    """

    def is_exact(self):
        return False

    def _latex_(self):
        return "\\R"

    def __repr__(self):
        """
        Print out this real double field.

        EXAMPLES:
            sage: RealDoubleField()
            Real Double Field
            sage: RDF
            Real Double Field
        """
        return "Real Double Field"

    def __cmp__(self, x):
        """
        EXAMPLES:
            sage: RDF == 5
            False
            sage: loads(dumps(RDF)) == RDF
            True
        """
        if PY_TYPE_CHECK(x, RealDoubleField_class):
            return 0
        return cmp(type(self), type(x))

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
            1.4142135623730950488016887242096980785696718753769480731767
            sage: b = RDF(a); b
            1.41421356237
            sage: a.parent()(b)
            1.4142135623730951454746218587388284504413604736328125000000
            sage: a.parent()(b) == b
            True
            sage: b == RR(a)
            True
        """
        if hasattr(x, '_real_double_'):
            return x._real_double_(self)
        return RealDoubleElement(x)

    cdef _coerce_c_impl(self, x):
        """
        Canonical coercion of x to the real double field.

        The rings that canonically coerce to the real double field are:
             * the real double field itself
             * int, long, integer, and rational rings
             * real mathematical constants
             * the mpfr real field

        EXAMPLES:
            sage: RDF._coerce_(5)
            5.0
            sage: RDF._coerce_(9499294r)
            9499294.0
            sage: RDF._coerce_(61/3)
            20.3333333333
            sage: parent(RDF(3) + CDF(5))
            Complex Double Field
            sage: parent(CDF(5) + RDF(3))
            Complex Double Field
        """
        if isinstance(x, (int, long, sage.rings.integer.Integer,
                          sage.rings.rational.Rational)):
            return self(x)
        import real_mpfr
        return self._coerce_try(x, [sage.functions.constants.ConstantRing,
                                    real_mpfr.RR])


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
            sage: RDF.is_atomic_repr()
            True
        """
        return True

    def is_finite(self):
        """
        Returns False, since the field of real numbers is not finite.
        Technical note:  There exists an upper bound on the double representation.

        EXAMPLES:
            sage: RDF.is_finite()
            False
        """
        return False

    def characteristic(self):
        """
        Returns 0, since the field of real numbers has characteristic 0.

        EXAMPLES:
            sage: RDF.characteristic()
            0
        """
        return 0

    cdef _new_c(self, double value):
        cdef RealDoubleElement x
        x = PY_NEW(RealDoubleElement)
        x._value = value
        return x

    def random_element(self, double min=-1, double max=1):
        """
        Return a random element of this real double field in the interval [min, max].

        EXAMPLES:
	    sage: RDF.random_element()
	    -0.657233364114
	    sage: RDF.random_element(min=100, max=110)
	    106.592535785
        """
        return self._new_c((max-min)*(<double>random())/RAND_MAX + min)

    def name(self):
        return "RealDoubleField"

    def __hash__(self):
        return 1455926870 #return hash(self.name())

    def pi(self):
        """
        Returns pi to double-precision.

        EXAMPLES:
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



def new_RealDoubleElement():
    cdef RealDoubleElement x
    x = PY_NEW(RealDoubleElement)
    return x

cdef class RealDoubleElement(FieldElement):
    def __new__(self, x=None):
        (<Element>self)._parent = _RDF

    def __init__(self, x):
        self._value = float(x)

    def __reduce__(self):
        """
        EXAMPLES:
            sage: a = RDF(-2.7)
            sage: loads(dumps(a)) == a
            True
        """
        return RealDoubleElement, (self._value, )

    cdef _new_c(self, double value):
        cdef RealDoubleElement x
        x = PY_NEW(RealDoubleElement)
        x._value = value
        return x

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
        return self._parent

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
        s = self.str()
        parts = s.split('e')
        if len(parts) > 1:
            # scientific notation
            if parts[1][0] == '+':
                parts[1] = parts[1][1:]
            s = "%s \\times 10^{%s}" % (parts[0], parts[1])
        return s

    def __hash__(self):
        return hash(float(self))

    def _im_gens_(self, codomain, im_gens):
        return codomain(self) # since 1 |--> 1

    def str(self):
        """
        Return string representation of self.

        EXAMPLES:
            sage: a = RDF('4.5'); a.str()
            '4.5'
            sage: a = RDF('49203480923840.2923904823048'); a.str()
            '4.92034809238e+13'
            sage: a = RDF(1)/RDF(0); a.str()
            'inf'
            sage: a = -RDF(1)/RDF(0); a.str()
            '-inf'
            sage: a = RDF(0)/RDF(0); a.str()
            'nan'
        """
        if gsl_isnan(self._value):
            return "nan"
        else:
            v = gsl_isinf(self._value)
            if v == 1:
                return "inf"
            elif v == -1:
                return "-inf"
        return str(self._value)

    def __copy__(self):
        """
        Return copy of self, which since self is immutable, is just self.

        EXAMPLES:
            sage: r = RDF('-1.6')
            sage: r.__copy__() is r
            True
        """
        return self

    def integer_part(self):
        """
        If in decimal this number is written n.defg, returns n.

        EXAMPLES:
            sage: r = RDF('-1.6')
            sage: a = r.integer_part(); a
            -1
            sage: type(a)
            <type 'sage.rings.integer.Integer'>
        """
        return sage.rings.integer.Integer(int(self._value))


    ########################
    #   Basic Arithmetic
    ########################
    def __invert__(self):
        """
        Compute the multiplicative inverse of self.

        EXAMPLES:
            sage: a = RDF(-1.5)*RDF(2.5)
            sage: a.__invert__()
            -0.266666666667
	    sage: ~a
            -0.266666666667
        """
        cdef RealDoubleElement x = <RealDoubleElement>PY_NEW(RealDoubleElement)
        x._value = 1.0 / self._value
        return x

    cdef ModuleElement _add_c_impl(self, ModuleElement right):
        """
        Add two real numbers with the same parent.

        EXAMPLES:
            sage: RDF('-1.5') + RDF('2.5')
            1.0
        """
        cdef RealDoubleElement x = <RealDoubleElement>PY_NEW(RealDoubleElement)
        x._value = self._value + (<RealDoubleElement>right)._value
        return x

    cdef ModuleElement _sub_c_impl(self, ModuleElement right):
        """
        Subtract two real numbers with the same parent.

        EXAMPLES:
            sage: RDF('-1.5') - RDF('2.5')
            -4.0
        """
        cdef RealDoubleElement x = <RealDoubleElement>PY_NEW(RealDoubleElement)
        x._value = self._value - (<RealDoubleElement>right)._value
        return x

    cdef RingElement _mul_c_impl(self, RingElement right):
        """
        Multiply two real numbers with the same parent.

        EXAMPLES:
            sage: RDF('-1.5') * RDF('2.5')
            -3.75
        """
        cdef RealDoubleElement x = <RealDoubleElement>PY_NEW(RealDoubleElement)
        x._value = self._value * (<RealDoubleElement>right)._value
        return x

    cdef RingElement _div_c_impl(self, RingElement right):
        """
        EXAMPLES:
            sage: RDF('-1.5') / RDF('2.5')
            -0.6
            sage: RDF(1)/RDF(0)
            inf
        """
        cdef RealDoubleElement x = <RealDoubleElement>PY_NEW(RealDoubleElement)
        x._value = self._value / (<RealDoubleElement>right)._value
        return x

    cdef ModuleElement _neg_c_impl(self):
        """
        Negates a real number.

        EXAMPLES:
            sage: -RDF('-1.5')
            1.5
        """
        cdef RealDoubleElement x = <RealDoubleElement>PY_NEW(RealDoubleElement)
        x._value = -self._value
        return x

    def __abs__(self):
        return self.abs()

    cdef RealDoubleElement abs(RealDoubleElement self):
        return RealDoubleElement(abs(self._value))

    def __lshift__(x, y):
        """
        LShifting a double is not supported; nor is lshifting a RealDoubleElement.
        """
        raise TypeError, "unsupported operand type(s) for <<"

    def __rshift__(x, y):
        """
        RShifting a double is not supported; nor is rshifting a RealDoubleElement.
        """
        raise TypeError, "unsupported operand type(s) for >>"

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
        return sage.rings.integer.Integer(int(math.floor(self._value)))

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
        return sage.rings.integer.Integer(int(math.ceil(self._value)))

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
        return self._new_c(self._value - int(self._value))

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
        return sage.rings.complex_double.ComplexDoubleField()(self)

    def _pari_(self):
        cdef sage.libs.pari.gen.PariInstance P = sage.libs.pari.gen.pari
        return P.double_to_gen_c(self._value)


    ###########################################
    # Comparisons: ==, !=, <, <=, >, >=
    ###########################################

    def is_NaN(self):
        return bool(gsl_isnan(self._value))

    def is_positive_infinity(self):
        return bool(gsl_isinf(self._value) > 0)

    def is_negative_infinity(self):
        return bool(gsl_isinf(self._value) < 0)

    def is_infinity(self):
        return bool(gsl_isinf(self._value))

    def __richcmp__(left, right, int op):
        return (<Element>left)._richcmp(right, op)

    cdef int _cmp_c_impl(left, Element right) except -2:
        if left._value < (<RealDoubleElement>right)._value:
            return -1
        elif left._value > (<RealDoubleElement>right)._value:
            return 1
        return 0


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
            sage: r = RDF(4.0)
            sage: r.sqrt()
            2.0
            sage: r.sqrt()^2 == r
            True

            sage: r = RDF(4344)
            sage: r.sqrt()
            65.9090282131
            sage: r.sqrt()^2 - r             # random low order bits
            0.0

            sage: r = RDF(-2.0)
            sage: r.sqrt()
            1.41421356237*I
            """
        if self._value >= 0:
            return self.square_root()
        return self._complex_double_().sqrt()


    def square_root(self):
        """
        Return a square root of self.  A real number will always be
        returned (though it will be NaN if self is negative).

        Use self.sqrt() to get a complex number if self is negative.

        EXAMPLES:
            sage: r = RDF(-2.0)
            sage: r.square_root()
            nan
            sage: r.sqrt()
            1.41421356237*I
        """
        return self._new_c(sqrt(self._value))

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
        return self._new_c(gsl_sf_exp(gsl_sf_log(self._value) * exponent._value))

    def __pow_int(self, int exponent):
        return self._new_c(gsl_pow_int(self._value, exponent))

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

        Symbolic examples:
            sage: RDF('-2.3')^(x+y^3+sin(x))
            -2.30000000000000^(y^3 + sin(x) + x)
            sage: RDF('-2.3')^x
            -2.30000000000000^x
        """
        cdef RealDoubleElement x
        if PY_TYPE_CHECK(exponent, RealDoubleElement):
            return self.__pow(RealDoubleElement(exponent))
        elif PY_TYPE_CHECK(exponent, int):
            return self.__pow_int(exponent)
        elif PY_TYPE_CHECK(exponent, Integer) and exponent < INT_MAX:
            return self.__pow_int(int(exponent))
        try:
            x = self.parent()(exponent)
        except TypeError:
            try:
                return exponent.parent()(self)**exponent
            except AttributeError:
                raise TypeError
        return self.__pow(x)


    def _log_base(self, double log_of_base):
        if self._value < 2:
            if self._value == 0:
                return -1./0
            if self._value < 0:
                return 0./0
            _sig_on
            a = self._new_c(gsl_sf_log_1plusx(self._value - 1) / log_of_base)
            _sig_off
            return a
        _sig_on
        a = self._new_c(gsl_sf_log(self._value) / log_of_base)
        _sig_off
        return a

    def log(self, base=None):
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
        if base is None:
            return self._log_base(1)
        else:
            if isinstance(base, RealDoubleElement):
                return self._log_base(base._log_base(1))
            else:
                return self._log_base(gsl_sf_log(float(base)))

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
        _sig_on
        a = self._new_c(gsl_sf_log(self._value) / M_LN2)
        _sig_off
        return a


    def log10(self):
        """
        Returns log to the base 10 of self

        EXAMPLES:
            sage: r = RDF('16.0'); r.log10()
            1.20411998266
            sage: r.log() / RDF(log(10))
            1.20411998266
            sage: r = RDF('39.9'); r.log10()
            1.60097289569
        """
        _sig_on
        a = self._new_c(gsl_sf_log(self._value) / M_LN10)
        _sig_off
        return a

    def logpi(self):
        """
        Returns log to the base pi of self

        EXAMPLES:
            sage: r = RDF(16); r.logpi()
            2.42204624559
            sage: r.log() / RDF(log(pi))
            2.42204624559
            sage: r = RDF('39.9'); r.logpi()
            3.22030233461
        """
        _sig_on
        a = self._new_c(gsl_sf_log(self._value) / M_LNPI)
        _sig_off
        return a

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

            sage: RDF(1000).exp()
            inf
        """
        _sig_on
        a = self._new_c(gsl_sf_exp(self._value))
        _sig_off
        return a

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
        _sig_on
        a = self._new_c(gsl_sf_exp(self._value * M_LN2))
        _sig_off
        return a

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
        _sig_on
        a = self._new_c(gsl_sf_exp(self._value * M_LN10))
        _sig_off
        return a

    def cos(self):
        """
        Returns the cosine of this number

        EXAMPLES:
            sage: t=RDF.pi()/2
            sage: t.cos()
            6.12323399574e-17
        """
        return self._new_c(gsl_sf_cos(self._value))

    def sin(self):
        """
        Returns the sine of this number

        EXAMPLES:
            sage: RDF(2).sin()
            0.909297426826
        """
        return self._new_c(gsl_sf_sin(self._value))

    def restrict_angle(self):
        return self._new_c(gsl_sf_angle_restrict_symm(self._value))

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
        cdef double denom
        cos = gsl_sf_cos(self._value)
        a = self._new_c(gsl_sf_sin(self._value) / cos)
        return a

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
        _sig_on
        a = self._new_c(gsl_sf_hypot(self._value, float(other)))
        _sig_off
        return a

    def acos(self):
        """
        Returns the inverse cosine of this number

        EXAMPLES:
            sage: q = RDF.pi()/3
            sage: i = q.cos()
            sage: i.acos() == q
            True
        """
        return self._new_c(acos(self._value))

    def asin(self):
        """
        Returns the inverse sine of this number

        EXAMPLES:
            sage: q = RDF.pi()/5
            sage: i = q.sin()
            sage: i.asin() == q
            True
        """
        return self._new_c(asin(self._value))

    def atan(self):
        """
        Returns the inverse tangent of this number

        EXAMPLES:
            sage: q = RDF.pi()/5
            sage: i = q.tan()
            sage: i.atan() == q
            True
        """
        return self._new_c(atan(self._value))


    def cosh(self):
        """
        Returns the hyperbolic cosine of this number

        EXAMPLES:
            sage: q = RDF.pi()/12
            sage: q.cosh()
            1.0344656401
        """
        return self._new_c(gsl_ldexp( gsl_sf_exp(self._value) + gsl_sf_exp(-self._value), -1)) # (e^x + x^-x)/2

    def sinh(self):
        """
        Returns the hyperbolic sine of this number

        EXAMPLES:
            sage: q = RDF.pi()/12
            sage: q.sinh()
            0.264800227602

        """
        return self._new_c(gsl_ldexp( gsl_sf_expm1(self._value) - gsl_sf_expm1(-self._value), -1)) # (e^x - x^-x)/2

    def tanh(self):
        """
        Returns the hyperbolic tangent of this number

        EXAMPLES:
            sage: q = RDF.pi()/12
            sage: q.tanh()
            0.255977789246
        """
        return self.sinh() / self.cosh()

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
        return self._new_c(gsl_acosh(self._value))

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
        return self._new_c(gsl_asinh(self._value))

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
        return self._new_c(gsl_atanh(self._value))

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
        return self._new_c(gsl_sf_erf(self._value))

    def gamma(self):
        """
        The Euler gamma function. Return gamma of self.

        EXAMPLES:
           sage: RDF(6).gamma()
           120.0
           sage: RDF(1.5).gamma()
           0.886226925453
        """
        _sig_on
        a = self._new_c(gsl_sf_gamma(self._value))
        _sig_off
        return a

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
            return self._new_c(1)/self._new_c(0)
        return self._new_c(gsl_sf_zeta(self._value))

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
            x^2 - 2
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
            x^2 - 2
        """
        return sage.rings.arith.algdep(self,n)


#####################################################
# unique objects
#####################################################
cdef RealDoubleField_class _RDF
_RDF = RealDoubleField_class()

RDF = _RDF   # external interface

def RealDoubleField():
    global _RDF
    return _RDF


def is_RealDoubleElement(x):
    return PY_TYPE_CHECK(x, RealDoubleElement)






################# FAST CREATION CODE ######################
########### Based on fast integer creation code   #########
######## There is nothing to see here, move along   #######

cdef extern from *:

    ctypedef struct RichPyObject "PyObject"

    # We need a PyTypeObject with elements so we can
    # get and set tp_new, tp_dealloc, tp_flags, and tp_basicsize
    ctypedef struct RichPyTypeObject "PyTypeObject":

        # We replace this one
        PyObject*      (*    tp_new) ( RichPyTypeObject*, PyObject*, PyObject*)

        # Not used, may be useful to determine correct memory management function
        RichPyObject *(*   tp_alloc) ( RichPyTypeObject*, size_t )

        # We replace this one
        void           (*tp_dealloc) ( PyObject*)

        # Not used, may be useful to determine correct memory management function
        void          (*    tp_free) ( PyObject* )

        # We set a flag here to circumvent the memory manager
        long tp_flags

    cdef long Py_TPFLAGS_HAVE_GC

    # We need a PyObject where we can get/set the refcnt directly
    # and access the type.
    ctypedef struct RichPyObject "PyObject":
        int ob_refcnt
        RichPyTypeObject* ob_type

    # Allocation
    RichPyObject* PyObject_MALLOC(int)

    # Useful for debugging, see below
    void PyObject_INIT(RichPyObject *, RichPyTypeObject *)

    # Free
    void PyObject_FREE(PyObject*)

# We use a global element to steal all the references
# from.  DO NOT INITIALIZE IT AGAIN and DO NOT REFERENCE IT!
cdef RealDoubleElement global_dummy_element
global_dummy_element = RealDoubleElement(0)

# A global pool for performance when elements are rapidly created and destroyed.
# It operates on the following principles:
#
# - The pool starts out empty.
# - When an new element is needed, one from the pool is returned
#   if available, otherwise a new Integer object is created
# - When an element is collected, it will add it to the pool
#   if there is room, otherwise it will be deallocated.

cdef enum:
    element_pool_size = 50 # Pyrex has no way of defining constants

cdef PyObject* element_pool[element_pool_size]
cdef int element_pool_count = 0

# used for profiling the pool
cdef int total_alloc = 0
cdef int use_pool = 0

# The signature of tp_new is
# PyObject* tp_new(RichPyTypeObject *t, PyObject *a, PyObject *k).
# However we don't actually use any of it.
#
# t in this case is the RealDoubleElement TypeObject.

cdef PyObject* fast_tp_new(RichPyTypeObject *t, PyObject *a, PyObject *k):

    global element_pool, element_pool_count, total_alloc, use_pool

    cdef RichPyObject* new

    # for profiling pool usage
    # total_alloc += 1

    # If there is a ready integer in the pool, we will
    # decrement the counter and return that.

    if element_pool_count > 0:

        # for profiling pool usage
        # use_pool += 1

        element_pool_count -= 1
        new = <RichPyObject *> element_pool[element_pool_count]

    # Otherwise, we have to create one.

    else:

        # allocate enough room for the Integer, sizeof_Integer is
        # sizeof(Integer). The use of PyObject_MALLOC directly
        # assumes that Integers are not garbage collected, i.e.
        # they do not pocess references to other Python
        # objects (Aas indicated by the Py_TPFLAGS_HAVE_GC flag).
        # See below for a more detailed description.

        new = PyObject_MALLOC( sizeof(RealDoubleElement) )

        # Now set every member as set in z, the global dummy Integer
        # created before this tp_new started to operate.

        memcpy(new, (<void*>global_dummy_element), sizeof(RealDoubleElement) )

        # This line is only needed if Python is compiled in debugging
        # mode './configure --with-pydebug'. If that is the case a Python
        # object has a bunch of debugging fields which are initialized
        # with this macro. For speed reasons, we don't call it if Python
        # is not compiled in debug mode. So uncomment the following line
        # if you are debugging Python.

        #PyObject_INIT(new, (<RichPyObject*>global_dummy_element).ob_type)

    # The global_dummy_element may have a reference count larger than
    # one, but it is expected that newly created objects have a
    # reference count of one. This is potentially unneeded if
    # everybody plays nice, because the gobal_dummy_Integer has only
    # one reference in that case.

    # Objects from the pool have reference count zero, so this
    # needs to be set in this case.

    new.ob_refcnt = 1

    return new

cdef void fast_tp_dealloc(PyObject* o):

    # If there is room in the pool for a used integer object,
    # then put it in rather than deallocating it.

    global element_pool, element_pool_count

    if element_pool_count < element_pool_size:

        # And add it to the pool.
        element_pool[element_pool_count] = o
        element_pool_count += 1
        return

    # Free the object. This assumes that Py_TPFLAGS_HAVE_GC is not
    # set. If it was set another free function would need to be
    # called.

    PyObject_FREE(o)

hook_fast_tp_functions()

def hook_fast_tp_functions():
    """
    """
    global global_dummy_element

    cdef long flag

    cdef RichPyObject* o
    o = <RichPyObject*>global_dummy_element

    # By default every object created in Pyrex is garbage
    # collected. This means it may have references to other objects
    # the Garbage collector has to look out for. We remove this flag
    # as the only reference an Integer has is to the global Integer
    # ring. As this object is unique we don't need to garbage collect
    # it as we always have a module level reference to it. If another
    # attribute is added to the Integer class this flag removal so as
    # the alloc and free functions may not be used anymore.
    # This object will still be reference counted.
    flag = Py_TPFLAGS_HAVE_GC
    o.ob_type.tp_flags = <long>(o.ob_type.tp_flags & (~flag))

    # Finally replace the functions called when an Integer needs
    # to be constructed/destructed.
    o.ob_type.tp_new = &fast_tp_new
    o.ob_type.tp_dealloc = &fast_tp_dealloc

def time_alloc_list(n):
    cdef int i
    l = []
    for i from 0 <= i < n:
        l.append(PY_NEW(RealDoubleElement))

    return l

def time_alloc(n):
    cdef int i
    for i from 0 <= i < n:
        z = PY_NEW(RealDoubleElement)

def pool_stats():
    print "Used pool %s / %s times" % (use_pool, total_alloc)
    print "Pool contains %s / %s items" % (integer_pool_count, integer_pool_size)
