"""nodoctest
Field of quad double numbers
"""
include '../ext/cdefs.pxi'
include '../ext/interrupt.pxi'
#include '../ext/stdsage.pxi'

import operator
from random import random

from sage.misc.sage_eval import sage_eval

import sage.rings.complex_double
import sage.rings.complex_field

from sage.rings.integer import Integer
from sage.rings.rational import Rational
from real_mpfr import RealNumber
from real_double import RealDoubleElement

import sage.structure.coerce
from sage.structure.parent_base cimport ParentWithBase
from sage.structure.parent_gens cimport ParentWithGens

cdef class RealQuadDoubleField_class(Field):
    """
    Quad Double Field
    """

    def is_exact(self):
        return False

    def _latex_(self):
        return "\\R"

    def __repr__(self):
        """
        Print out this real double field.

        EXAMPLES:
            sage: RQDF = RealQuadDoubleField ()
            sage: RQDF
            Quad Double Real Field

        """
        return "Quad Double Real Field"

    def __cmp__(self, x):
        """
        EXAMPLES:
            sage: RQDF = RealQuadDoubleField ()
            sage: RQDF == 1
            False
            sage: RQDF == RealQuadDoubleField ()
            True

        """
        if PY_TYPE_CHECK(x, RealQuadDoubleField_class):
            return 0
        return cmp(type(self), type(x))

    def __call__(self, x):
        """
        Create a real quad double using x.

        EXAMPLES:
            sage: RQDF = RealQuadDoubleField ()
            sage: RQDF ('-1')
            -1.00000e+00

            sage: RQDF (-1/3)
            -3.33333e-01

            sage: RQDF (3.01)
            3.01000e+00
        """
        return QuadDoubleElement(x)

    def gen(self, i=0):
        if i == 0:
            # TODO:  this currently
            # does not accept an int
            return self('1')
        else:
            raise IndexError

    def ngens(self):
        return 1

    def gens(self):
        return [self.gen()]

    cdef _coerce_c_impl(self, x):
        if isinstance(x, (int, long, RealNumber, Integer, Rational,RealDoubleElement)):
            return self(x)

        import sage.functions.constants
        return self._coerce_try(x, [sage.functions.constants.ConstantRing])

    def name(self):
        return "QuadDoubleField"

    def __hash__(self):
        return hash(self.name())

    def pi(self):
        """
        Returns pi

        EXAMPLES:
            sage: R = RealField(100)
            sage: RQDF.pi()
            3.14159265359
            sage: RQDF.pi().sqrt()/2
            0.886226925453
        """
        cdef qd z
        cdef char *s
        s = <char*>PyMem_Malloc(sizeof(char)+10+8) # See docs for write()
        _sig_on
        z._pi.write(s,10,0,1)
        _sig_off
        return QuadDoubleElement(str(s))

    def log2(self):
        """
        Returns log(2)

        EXAMPLES:

        """
        cdef qd z
        cdef char *s
        s = <char*>PyMem_Malloc(sizeof(char)+10+8) # See docs for write()
        _sig_on
        z._log2.write(s,10,0,1)
        _sig_off
        return QuadDoubleElement(str(s))

cdef class QuadDoubleElement(FieldElement):
    """
    A quad double real number
    """

    cdef _new(self):
        cdef QuadDoubleElement q
        q = PY_NEW(QuadDoubleElement)
        q.initptr = new_qd_real()
        return q

    cdef _new_c(self, qd a):
        cdef QuadDoubleElement q
        q = PY_NEW(QuadDoubleElement)
        q.initptr =  qd_from_qd(a.x[0],a.x[1],a.x[2],a.x[3])
        return q

    def __new__(self, x=None):
        # explicit cast required for C++
        self._parent = <ParentWithBase> _RQDF

    def __init__(self, number):
        self._set(number)

    cdef _set(self, x):
        if PY_TYPE_CHECK(x, Rational) or \
               PY_TYPE_CHECK(x, RealNumber) or \
               PY_TYPE_CHECK(x, RealDoubleElement):

            value = float(x)
            self.initptr = qd_from_double(value)
        s = str(x)
        try:
            _sig_on
            self.initptr = qd_from_str(s)
            _sig_off
        except RuntimeError:
            raise TypeError

    def real(self):
        """
        Returns itself

        EXAMPLES:

        """
        return self

    def imag(self):
        """
        Returns the imaginary part of this number.

        EXAMPLES:

        """
        return QuadDoubleElement(0)

    def __complex__(self):
        """
        Returns self as a complex number
        EXAMPLES:

        """
        return complex(float(self))

    def __reduce__(self):
        """
        EXAMPLES:

        """
        s = str(self)
        return __create__QuadDoubleElement_version0, (s )

    def __str__(self):
        return self.str()

    def __repr__(self):
        return self.str()

    def str(self, precision=60):
        """
        Returns the string representation of self
        """
        # TODO: customize precision, sci_not, etc
        cdef char *s
        s = <char*>PyMem_Malloc(sizeof(char)+precision+8) # See docs for write()
        _sig_on
        self.initptr.write(s,precision,0,1)
        _sig_off
        return str(s)

    def parent(self):
        """
        Returns the parent of this number

        EXAMPLES
        """
        return RQDF

    def __copy__(self):
        """
        Return copy of self, which since self is immutable, is just self.

        EXAMPLES:

        """
        return self

    def integer_part(self):
        """
        If in decimal this number is written n.defg, returns n.

        EXAMPLES:

        """
        return Integer(int(self))

    ########################
    #   Basic Arithmetic
    ########################
    def __invert__(self):
        """
        Compute the multiplicative inverse of self.

        EXAMPLES:

        """
        cdef QuadDoubleElement res
        res = self._new()
        _sig_on
        c_qd_npwr(self.initptr.x,-1,res.initptr.x)
        _sig_off
        return res

    cdef ModuleElement _add_c_impl(self, ModuleElement right):
        """
        Add two quad double numbers
        """

        cdef QuadDoubleElement res
        res = self._new()
        c_qd_add(self.initptr.x,(<QuadDoubleElement>right).initptr.x,res.initptr.x)
        return res

    cdef ModuleElement _sub_c_impl(self, ModuleElement right):
        """
        Substract two quad double numbers
        """
        cdef QuadDoubleElement res
        res = self._new()
        c_qd_sub(self.initptr.x,(<QuadDoubleElement>right).initptr.x,res.initptr.x)
        return res

    cdef RingElement _mul_c_impl(self, RingElement right):
        """
        Multiply two quad double numbers
        """
        cdef QuadDoubleElement res
        res = self._new()
        #res = PY_NEW(QuadDoubleElement)
        #res.initptr = new_qd_real()
        c_qd_mul(self.initptr.x,(<QuadDoubleElement>right).initptr.x,res.initptr.x)
        return res

    cdef RingElement _div_c_impl(self, RingElement right):
        """
        Divide two quad double numbers
        """
        cdef QuadDoubleElement res
        res = self._new()
        c_qd_div(self.initptr.x,(<QuadDoubleElement>right).initptr.x,res.initptr.x)
        return res

    cdef ModuleElement _neg_c_impl(self):
        """
        Negates a real number.

        EXAMPLES:

        """
        cdef QuadDoubleElement res
        res = self._new()
        c_qd_neg(self.initptr.x,res.initptr.x)
        return res

    def __abs__(self):
        cdef QuadDoubleElement res
        res = self._new()
        c_qd_abs(self.initptr.x,res.initptr.x)
        return res

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
        if self == 1: return 1
        if self == -1: return -1
        return sage.rings.infinity.infinity


    ###################
    # Rounding etc
    ###################

    def round(self):
        """
        Given real number x, rounds up if fractional part is greater than .5,
        rounds down if fractional part is lesser than .5.
        EXAMPLES:
            sage: RQDF(0.49).round()
            0.0
            sage: RQDF(0.51).round()
            1.0
        """
        return -1 #TODO

    def floor(self):
        """
        Returns the floor of this number

        EXAMPLES:
            sage: RQDF(2.99).floor()
            2
            sage: RQDF(2.00).floor()
            2
            sage: RQDF(-5/2).floor()
            -3
        """
        cdef QuadDoubleElement res
        res = self._new()
        c_qd_floor(self.initptr.x,res.initptr.x)
        return res.integer_part()

    def ceil(self):
        """
        Returns the ceiling of this number

        EXAMPLES:
            sage: RQDF(2.99).ceil()
            3
            sage: RQDF(2.00).ceil()
            2
            sage: RQDF(-5/2).ceil()
            -2
        """
        cdef QuadDoubleElement res
        res = self._new()
        c_qd_ceil(self.initptr.x,res.initptr.x)
        return res.integer_part()

    def ceiling(self):
        return self.ceil()

    def trunc(self):
        """
        Truncates this number (returns integer part).

        EXAMPLES:
            sage: RQDF(2.99).trunc()
            2.000000000E0
            sage: RQDF(-2.00).trunc()
            -2.000000000E0
            sage: RQDF(0.00).trunc()
            0.000000000E0
        """
        return QuadDoubleElement(self.floor())


    def frac(self):
        """
        frac returns a real number > -1 and < 1. it satisfies the
        relation:
            x = x.trunc() + x.frac()

        EXAMPLES:
            sage: RQDF(2.99).frac()
            0.99
            sage: RQDF(2.50).frac()
            0.5
            sage: RQDF(-2.79).frac()
            -0.79
        """
        return 0 #TODO

    ###########################################
    # Conversions
    ###########################################

    def __float__(self):
        cdef double d
        d = qd_to_double(qd_deref(self.initptr))
        return d

    def __int__(self):
        """
        Returns integer truncation of this real number.
        """
        cdef int i
        i = qd_to_int(qd_deref(self.initptr))
        return i

    def __long__(self):
        """
        Returns long integer truncation of this real number.
        """
        return long(self.__int__())


    def _complex_number_(self):
        return  sage.rings.complex_field.ComplexField()(float(self))

    def _complex_double_(self):
        return  sage.rings.complex_double.ComplexDoubleField(float(self))

    def _pari_(self):
        return  sage.libs.pari.all.pari.new_with_bits_prec("%.15e"%float(self), 64)


    ###########################################
    # Comparisons: ==, !=, <, <=, >, >=
    ###########################################

    def is_NaN(self):
        return bool(qd_is_nan(qd_deref(self.initptr)))

    def is_infinity(self):
        return bool(qd_is_inf(qd_deref(self.initptr)))

    def __richcmp__(left, right, int op):
        return (<Element>left)._richcmp(right, op)

    cdef int _cmp_c_impl(left, Element right) except -2:
        cdef int i
        c_qd_comp(left.initptr.x,(<QuadDoubleElement>right).initptr.x,&i)
        return i

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
            2.00000000000000
            sage: r.sqrt()^2 == r
            True

            sage: r = 4344
            sage: r.sqrt()
            65.9090282131363
            sage: r.sqrt()^2 - r
            0.000000000000000

            sage: r = -2.0
            sage: r.sqrt()
            1.41421356237309*I
            """
        return self.square_root()
    #return self._complex_number_().sqrt()

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
            1.41421356237309*I
        """
        cdef QuadDoubleElement res
        res = self._new()
        c_qd_sqrt(self.initptr.x, res.initptr.x)
        return res

    def cube_root(self):
        """
        Return the cubic root (defined over the real numbers) of self.

        EXAMPLES:
            sage: r = RQDF(125.0); r.cube_root()
            5.000000000E0
            sage: r = RQDF(-119.0)
            sage: r.cube_root()^3 - r         # output is random, depending on arch.
            0.0
        """
        return self.nth_root(3)


    def nth_root(self, int n):
        """
        Returns the $n^{th}$ root of self.
        EXAMPLES:
            sage: r = RQDF(-125.0); r.nth_root(3)
            -5.0
            sage: r.nth_root(5)
            -2.6265278044
        """

        cdef QuadDoubleElement res
        res = self._new()
        _sig_on
        c_qd_nroot(self.initptr.x, n, res.initptr.x)
        _sig_off
        return res

    def __pow(self, n, modulus):
        cdef QuadDoubleElement res
        res = self._new()
        _sig_on
        c_qd_npwr(self.initptr.x, n, res.initptr.x)
        _sig_off
        return res

    def __pow__(self,n,d):
        """
        Compute self raised to the power of exponent, rounded in
        the direction specified by the parent of self.

        If the result is not a real number, self and the exponent are
        both coerced to complex numbers (with sufficient precision),
        then the exponentiation is computed in the complex numbers.
        Thus this function can return either a real or complex number.

        EXAMPLES:
            sage: a = RQDF('1.23456')
            sage: a^20
            6.764629770E1

        """
        return self.__pow(n,d)

    def log(self):
        """
        EXAMPLES:
            sage: RQDF(2).log()
            6.931471806E-1
            sage: RQDF(0).log()
            -inf
            sage: RQDF(-1).log()
            nan
        """
        cdef QuadDoubleElement res
        res = self._new()
        _sig_on
        c_qd_log(self.initptr.x,res.initptr.x)
        _sig_off
        return res

    def log10(self):
        """
        Returns log to the base 10 of self

        EXAMPLES:
            sage: r = RQDF('16.0'); r.log10()
            1.204119983E0
            sage: r.log() / log(10)
            1.20411998266
            sage: r = RQDF('39.9'); r.log10()
            1.600972896E0
        """
        cdef QuadDoubleElement res
        res = self._new()
        _sig_on
        c_qd_log10(self.initptr.x,res.initptr.x)
        _sig_off
        return res

    def exp(self):
        r"""
        Returns $e^\code{self}$

        EXAMPLES:
            sage: r = RQDF(0.0)
            sage: r.exp()
            1.000000000E0

            sage: r = RQDF('32.3')
            sage: a = r.exp(); a
            1.065888473E14
            sage: a.log()
            3.230000000E1

            sage: RQDF('-32.3').exp()
            9.381844588E-15

        """
        cdef QuadDoubleElement res
        res = self._new()
        _sig_on
        c_qd_exp(self.initptr.x,res.initptr.x)
        _sig_off
        return res

    def cos(self):
        """
        Returns the cosine of this number

        EXAMPLES:
            sage: t=RQDF.pi()/2
            sage: t.cos()
            6.12323399574e-17
        """
        cdef QuadDoubleElement res
        res = self._new()
        _sig_on
        c_qd_cos(self.initptr.x,res.initptr.x)
        _sig_off
        return res

    def sin(self):
        """
        Returns the sine of this number

        EXAMPLES:
            sage: RQDF(2).sin()
            0.909297426826
        """
        cdef QuadDoubleElement res
        res = self._new()
        _sig_on
        c_qd_sin(self.initptr.x,res.initptr.x)
        _sig_off
        return res

    def tan(self):
        """
        Returns the tangent of this number

        EXAMPLES:
            sage: q = RQDF.pi()/3
            sage: q.tan()
            1.73205080757
            sage: q = RQDF.pi()/6
            sage: q.tan()
            0.57735026919
        """
        cdef QuadDoubleElement res
        res = self._new()
        _sig_on
        c_qd_tan(self.initptr.x,res.initptr.x)
        _sig_off
        return res

    def sincos(self):
        """
        Returns a pair consisting of the sine and cosine.

        EXAMPLES:
            sage: t = RQDF.pi()/6
            sage: t.sincos()
            (0.5, 0.866025403784)
        """
        return self.sin(), self.cos()

    def acos(self):
        """
        Returns the inverse cosine of this number

        EXAMPLES:
            sage: q = RQDF.pi()/3
            sage: i = q.cos()
            sage: i.acos() == q
            True
        """
        cdef QuadDoubleElement res
        res = self._new()
        _sig_on
        c_qd_acos(self.initptr.x,res.initptr.x)
        _sig_off
        return res

    def asin(self):
        """
        Returns the inverse sine of this number

        EXAMPLES:
            sage: q = RQDF.pi()/5
            sage: i = q.sin()
            sage: i.asin() == q
            True
        """
        cdef QuadDoubleElement res
        res = self._new()
        _sig_on
        c_qd_asin(self.initptr.x,res.initptr.x)
        _sig_off
        return res

    def atan(self):
        """
        Returns the inverse tangent of this number

        EXAMPLES:
            sage: q = RQDF.pi()/5
            sage: i = q.tan()
            sage: i.atan() == q
            True
        """
        cdef QuadDoubleElement res
        res = self._new()
        _sig_on
        c_qd_atan(self.initptr.x,res.initptr.x)
        _sig_off
        return res

    def cosh(self):
        """
        Returns the hyperbolic cosine of this number

        EXAMPLES:
            sage: q = RQDF.pi()/12
            sage: q.cosh()
            1.0344656401
        """
        cdef QuadDoubleElement res
        res = self._new()
        _sig_on
        c_qd_cosh(self.initptr.x,res.initptr.x)
        _sig_off
        return res

    def sinh(self):
        """
        Returns the hyperbolic sine of this number

        EXAMPLES:
            sage: q = RQDF.pi()/12
            sage: q.sinh()
            0.264800227602

        """
        cdef QuadDoubleElement res
        res = self._new()
        _sig_on
        c_qd_sinh(self.initptr.x,res.initptr.x)
        _sig_off
        return res

    def tanh(self):
        """
        Returns the hyperbolic tangent of this number

        EXAMPLES:
            sage: q = RQDF.pi()/12
            sage: q.tanh()
            0.255977789246
        """
        cdef QuadDoubleElement res
        res = self._new()
        _sig_on
        c_qd_tanh(self.initptr.x,res.initptr.x)
        _sig_off
        return res

    def acosh(self):
        """
        Returns the hyperbolic inverse cosine of this number

        EXAMPLES:
            sage: q = RQDF.pi()/2
            sage: i = q.cosh() ; i
            2.50917847866
            sage: i.acosh() == q
            True
        """
        cdef QuadDoubleElement res
        res = self._new()
        _sig_on
        c_qd_acosh(self.initptr.x,res.initptr.x)
        _sig_off
        return res

    def asinh(self):
        """
        Returns the hyperbolic inverse sine of this number

        EXAMPLES:
            sage: q = RQDF.pi()/2
            sage: i = q.sinh() ; i
            2.30129890231
            sage: i.asinh() == q
            True
        """
        cdef QuadDoubleElement res
        res = self._new()
        _sig_on
        c_qd_asinh(self.initptr.x,res.initptr.x)
        _sig_off
        return res

    def atanh(self):
        """
        Returns the hyperbolic inverse tangent of this number

        EXAMPLES:
            sage: q = RQDF.pi()/2
            sage: i = q.tanh() ; i
            0.917152335667
            sage: i.atanh() - q      # output is random, depending on arch.
            -4.4408920985e-16
        """
        cdef QuadDoubleElement res
        res = self._new()
        _sig_on
        c_qd_atanh(self.initptr.x,res.initptr.x)
        _sig_off
        return res

    def agm(self, other):
        """
        Return the arithmetic-geometric mean of self and other. The
        arithmetic-geometric mean is the common limit of the sequences
        $u_n$ and $v_n$, where $u_0$ is self, $v_0$ is other,
        $u_{n+1}$ is the arithmetic mean of $u_n$ and $v_n$, and
        $v_{n+1}$ is the geometric mean of u_n and v_n. If any operand
        is negative, the return value is \code{NaN}.
        """
        return QuadDoubleElement(sage.rings.all.RR(self).agm(sage.rings.all.RR(other)))


cdef RealQuadDoubleField_class _RQDF
_RQDF = RealQuadDoubleField_class()

RQDF = _RQDF   # external interface

def RealQuadDoubleField():
    global _RQDF
    return _RQDF

def is_QuadDoubleElement(x):
    return PY_TYPE_CHECK(x, QuadDoubleElement)

def __create__QuadDoubleField_version0():
    return QuadDoubleField()

def __create__QuadDoubleElement_version0(x):
    return QuadDoubleElement(x)
