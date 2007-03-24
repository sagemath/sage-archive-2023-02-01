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
            -1.0000000000000000000E0

            sage: RQDF (-1/3)
            -3.3333333333333331483E-1

            sage: RQDF('-0.33333333333333333333')
            -3.3333333333333333333E-1
        """
        return QuadDoubleElement(x)

    def gen(self, i=0):
        """
        Return the generator of the field.
        EXAMPLES:
            sage: RQDF.gen()
            1.0000000000000000000E0

        """
        if i == 0:
            return self(int(1))
        else:
            raise IndexError

    def ngens(self):
        """
        Return the number of generators of the field.
        EXAMPLES:
            sage: RQDF.ngens()
            1
        """
        return 1

    def gens(self):
        """
        Returns a list of the generators of this field

        EXAMPLES:
            sage: RQDF.gens()
            [1.0000000000000000000E0]
        """
        return [self.gen()]

    cdef _coerce_c_impl(self, x):
        if isinstance(x, (int, long, RealNumber, Integer, Rational,RealDoubleElement)):
            return self(x)

        import sage.functions.constants
        return self._coerce_try(x, [sage.functions.constants.ConstantRing])

    def name(self):
        return "QuadDoubleField"

    def __hash__(self):
        return 219746133

    def pi(self):
        """
        Returns pi

        EXAMPLES:
            sage: RQDF = RealQuadDoubleField ()
            sage: RQDF.pi()
            3.1415926535897932385E0
            sage: RQDF.pi().sqrt()/2
            8.8622692545275801365E-1

        """
        cdef qd z
        cdef char *s
        s = <char*>PyMem_Malloc(sizeof(char)+20+8) # See docs for write()
        _sig_on
        z._pi.write(s,20,0,1)
        _sig_off
        return QuadDoubleElement(str(s))

    def log2(self):
        """
        Returns log(2)

        EXAMPLES:
            sage: RQDF = RealQuadDoubleField ()
            sage: RQDF.log2()
            6.9314718055994530942E-1
        """
        cdef qd z
        cdef char *s
        s = <char*>PyMem_Malloc(sizeof(char)+20+8) # See docs for write()
        _sig_on
        z._log2.write(s,20,0,1)
        _sig_off
        return QuadDoubleElement(str(s))

    def random_element(self,x=0,y=1):
        """
        Generate a random quad double between x and y
        RQDF.random_element() -- random real number between 0 and 1
        RQDF.random_element(n) -- return an real number between 0 and n-1, inclusive.
        RQDF.random_element(min, max) -- return a real number between min and max-1, inclusive.

        EXAMPLES:
            sage: RQDF = RealQuadDoubleField ()
            sage: RQDF.random_element(-10,10)
            4.0887067399587205430E0

            sage: RQDF.random_element(10)
            3.6374177750879467821E0

            sage: [RQDF.random_element(10) for _ in range(5)]
            [9.7759052462604299129E0,
            9.2496342244458553030E0,
            2.9211666756880459942E0,
            1.4497166816703933368E0,
            1.5131396931839003209E0]
        """
        cdef QuadDoubleElement res, upper, lower
        res = QuadDoubleElement(0)
        upper = self(x)
        lower = self(y)
        c_qd_rand(res.initptr.x)
        return (upper-lower)*res + lower

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

        else:
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
           sage: RQDF = RealQuadDoubleField ()
           sage: w=RQDF(2) ; w
           2.0000000000000000000E0
        """
        return self

    def imag(self):
        """
        Returns the imaginary part of this number.

        EXAMPLES:
           sage: RQDF = RealQuadDoubleField ()
           sage: w=RQDF(2)
           sage: w.imag()
           0.0000000000000000000E0
        """
        return QuadDoubleElement(0)

    def __complex__(self):
        """
        Returns self as a complex number
        EXAMPLES:
           sage: RQDF = RealQuadDoubleField ()
           sage: w=RQDF(2)
           sage: complex(w)
           (2+0j)

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

    def str(self, precision=20):
        """
        Returns the string representation of self
        EXAMPLES:
           sage: RQDF = RealQuadDoubleField ()
           sage: w=RQDF(-21.2) ; str(w)
           '-2.1199999999999999289E1'
        """
        cdef char *s
        s = <char*>PyMem_Malloc(sizeof(char)+precision+8) # See docs for write()
        _sig_on
        self.initptr.write(s,precision,0,1)
        _sig_off
        return str(s)

    def parent(self):
        """
        Returns the parent of this number

        EXAMPLES:
           sage: RQDF = RealQuadDoubleField ()
           sage: w=RQDF(-21.2) ; w.parent()
           Quad Double Real Field
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
            sage: RQDF = RealQuadDoubleField ()
            sage: test = '253536646425647436353675786864535364746'
            sage: RQDF(test).integer_part()
            253536646425647436353675786864535364745

        """
        # TODO : because of rounding errors, this
        # hard to get for large integers
        s = self.str(300)
        digits = int(s.split('E')[1])
        # self is < 0
        if digits < 0: return Integer(0)
        num = s.split('E')[0].replace('.','')
        if s[0] == '-': # negative
            return Integer(num[:digits+2])
        else:
            return Integer(num[:digits+1])


    ########################
    #   Basic Arithmetic
    ########################
    def __invert__(self):
        """
        Compute the multiplicative inverse of self.

        EXAMPLES:
            sage: RQDF = RealQuadDoubleField ()
            sage: w=RQDF(1/3)
            sage: 1/w
            3.0000000000000001665E0

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

        EXAMPLES:
            sage: RQDF = RealQuadDoubleField ()
            sage: RQDF(1/3) + RQDF('1')
            1.3333333333333333148E0
        """

        cdef QuadDoubleElement res
        res = self._new()
        c_qd_add(self.initptr.x,(<QuadDoubleElement>right).initptr.x,res.initptr.x)
        return res

    cdef ModuleElement _sub_c_impl(self, ModuleElement right):
        """
        Substract two quad double numbers

        EXAMPLES:
            sage: RQDF = RealQuadDoubleField ()
            sage: RQDF(1/3) - RQDF('1')
            -6.6666666666666668517E-1
        """
        cdef QuadDoubleElement res
        res = self._new()
        c_qd_sub(self.initptr.x,(<QuadDoubleElement>right).initptr.x,res.initptr.x)
        return res

    cdef RingElement _mul_c_impl(self, RingElement right):
        """
        Multiply two quad double numbers

        EXAMPLES:
            sage: RQDF = RealQuadDoubleField ()
            sage: RQDF(1.3) * RQDF(10)
            1.3000000000000000444E1
        """
        cdef QuadDoubleElement res
        res = self._new()
        c_qd_mul(self.initptr.x,(<QuadDoubleElement>right).initptr.x,res.initptr.x)
        return res

    cdef RingElement _div_c_impl(self, RingElement right):
        """
        Divide two quad double numbers

        EXAMPLES:
            sage: RQDF = RealQuadDoubleField ()
            sage: RQDF(1/3) / RQDF(100)
            3.3333333333333331483E-3
        """
        cdef QuadDoubleElement res
        res = self._new()
        c_qd_div(self.initptr.x,(<QuadDoubleElement>right).initptr.x,res.initptr.x)
        return res

    cdef ModuleElement _neg_c_impl(self):
        """
        Negates a real number.

        EXAMPLES:
            sage: RQDF = RealQuadDoubleField ()
            sage: -RQDF('0.056')
            -5.6000000000000000000E-2
        """
        cdef QuadDoubleElement res
        res = self._new()
        c_qd_neg(self.initptr.x,res.initptr.x)
        return res

    def __abs__(self):
        """
        Negates a real number.

        EXAMPLES:
            sage: RQDF = RealQuadDoubleField ()
            sage: abs(RQDF('-0.45'))
            4.5000000000000000000E-1
        """
        cdef QuadDoubleElement res
        res = self._new()
        c_qd_abs(self.initptr.x,res.initptr.x)
        return res

    def __lshift__(x, y):
        """
        LShifting a quad double is not supported; nor is lshifting a RealDoubleElement.
        """
        raise TypeError, "unsupported operand type(s) for <<: '%s' and '%s'"%(typeof(self), typeof(n))

    def __rshift__(x, y):
        """
        RShifting a quad double is not supported; nor is rshifting a RealDoubleElement.
        """
        raise TypeError, "unsupported operand type(s) for >>: '%s' and '%s'"%(typeof(self), typeof(n))

    def multiplicative_order(self):
        """
        Returns the multiplicative order of self

        EXAMPLES:
            sage: RQDF = RealQuadDoubleField ()
            sage: w=RQDF(-1)
            sage: w.multiplicative_order()
            -1
            sage: w=RQDF(1)
            sage: w.multiplicative_order()
            1
            sage: w=RQDF(0)
            sage: w.multiplicative_order()
            +Infinity

        """
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
            sage: RQDF = RealQuadDoubleField ()
            sage: RQDF(0.49).round()
            0.0000000000000000000E0
            sage: RQDF(0.51).round()
            1.0000000000000000000E0
        """
        cdef QuadDoubleElement res
        res = self._new()
        if self.frac() < 0.5:
            _sig_on
            c_qd_floor(self.initptr.x,res.initptr.x)
            _sig_off
            return res

        _sig_on
        c_qd_ceil(self.initptr.x,res.initptr.x)
        _sig_off
        return res


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
            2.0000000000000000000E0
            sage: RQDF(-2.00).trunc()
            -2.0000000000000000000E0
            sage: RQDF(0.00).trunc()
            0.0000000000000000000E0
        """
        return QuadDoubleElement(self.floor())


    def frac(self):
        """
        frac returns a real number > -1 and < 1. that satisfies the
        relation:
            x = x.trunc() + x.frac()

        EXAMPLES:
            sage: RQDF(2.99).frac()
            9.9000000000000021316E-1
            sage: RQDF(2.50).frac()
            5.0000000000000000000E-1
            sage: RQDF(-2.79).frac()
            -7.9000000000000003553E-1
        """
        return self - self.integer_part()


    ###########################################
    # Conversions
    ###########################################

    def __float__(self):
        """
        Returns the floating-point value of this number

        EXAMPLES:
            sage: RQDF = RealQuadDoubleField ()
            sage: w=RQDF(-23.79)
            sage: float(w)
            -23.789999999999999

        """
        cdef double d
        d = qd_to_double(qd_deref(self.initptr))
        return d

    def __int__(self):
        """
        Returns integer truncation of this real number.
        EXAMPLES:
            sage: RQDF = RealQuadDoubleField ()
            sage: w=RQDF(-23.79)
            sage: int(w)
            -23

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
        """
        Returns True if self is NaN
        EXAMPLES:
            sage: RQDF = RealQuadDoubleField () ; w=RQDF(9)
            sage: w.is_NaN()
            False
            sage: w=RQDF(1<<3000)
            sage: w.is_NaN()
            True
        """
        return bool(qd_is_nan(qd_deref(self.initptr)))

    def is_infinity(self):
        """
        Returns True if self is infinifty
        EXAMPLES:

        """
        return bool(qd_is_inf(qd_deref(self.initptr)))

    def __richcmp__(left, right, int op):
        return (<Element>left)._richcmp(right, op)

    cdef int _cmp_c_impl(left, Element right) except -2:
        """
        Compares 2 quad double numbers

        Returns True if self is NaN
        EXAMPLES:
            sage: RQDF = RealQuadDoubleField ()
            sage: RQDF('3233') > RQDF('323')
            True
            sage: RQDF('-3233') > RQDF('323')
            False
            sage: RQDF('-3233') == RQDF('-3233')
            True
        """
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
            sage: RQDF = RealQuadDoubleField ()
            sage: RQDF('-3233').sqrt()
            56.8594759033180*I
            sage: RQDF('-4').sqrt()
            2.00000000000000*I
            sage: RQDF('4').sqrt()
            2.0000000000000000000E0
        """
        if self >=0:
            return self.square_root()
        return self._complex_number_().sqrt()


    def square_root(self):
        """
        Return a square root of self.  A real number will always be
        returned (though it will be NaN if self is negative).

        Use self.sqrt() to get a complex number if self is negative.

        EXAMPLES:
            sage: RQDF('4').square_root()
            2.0000000000000000000E0
            sage: RQDF('-4').square_root()
            NAN
        """
        cdef QuadDoubleElement res
        res = self._new()
        if self > 0:
            c_qd_sqrt(self.initptr.x, res.initptr.x)
            return res
        else:
            res = self._new_c(self.initptr._nan)
            return res

    def cube_root(self):
        """
        Return the cubic root (defined over the real numbers) of self.

        EXAMPLES:
            sage: RQDF = RealQuadDoubleField ()
            sage: r = RQDF(125.0); r.cube_root()
            5.0000000000000000000E0
            sage: RQDF('-4').cube_root()
            NAN
        """
        return self.nth_root(3)


    def nth_root(self, int n):
        """
        Returns the $n^{th}$ root of self.
        Returns NaN if self is negative
        EXAMPLES:
            sage: RQDF = RealQuadDoubleField ()
            sage: r = RQDF(125.0); r.nth_root(3)
            5.0000000000000000000E0
            sage: r.nth_root(5)
            2.6265278044037672365E0
            sage: RQDF('-4987').nth_root(2)
            NAN
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
            sage: RQDF = RealQuadDoubleField ()
            sage: a = RQDF('1.23456')
            sage: a^20
            6.7646297703853969093E1

        """
        return self.__pow(n,d)

    def log(self):
        """
        Returns the log of this number.
        Returns NaN is self < 0
        EXAMPLES:
            sage: RQDF = RealQuadDoubleField ()
            sage: RQDF(2).log()
            6.9314718055994530942E-1
            sage: RQDF(0).log()
            0.0000000000000000000E0
            sage: RQDF(-1).log()
            0.0000000000000000000E0
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
        Outputs an error if self < 0

        EXAMPLES:
            sage: RQDF = RealQuadDoubleField ()
            sage: r = RQDF('16.0'); r.log10()
            1.2041199826559247809E0
            sage: r.log() / log(10)
            1.2041199826559246673E0
            sage: r = RQDF('-16.0'); r.log10()
            0.0000000000000000000E0

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
            sage: RQDF = RealQuadDoubleField ()
            sage: r = RQDF(0.0) ; r.exp()
            1.0000000000000000000E0
            sage: RQDF('-32.3').exp()
            9.3818445884986577852E-15
            sage: r = RQDF('16.0');
            sage: r.log().exp() == r
            True
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
            sage: RQDF = RealQuadDoubleField ()
            sage: t=RQDF.pi()/2
            sage: t.cos()
            -1.8678308360248557901E-20
            sage: t.cos()^2 + t.sin()^2
            1.0000000000000000000E0
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
            sage: RQDF = RealQuadDoubleField ()
            sage: RQDF.pi().sin()
            -3.7356616720497115803E-20
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
            sage: RQDF = RealQuadDoubleField ()
            sage: q = RQDF.pi()/3
            sage: q.tan()
            1.7320508075688772936E0
            sage: q = RQDF.pi()/6
            sage: q.tan()
            5.7735026918962576452E-1
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
            sage: RQDF = RealQuadDoubleField ()
            sage: t = RQDF.pi()/6
            sage: t.sincos()
            (5.0000000000000000001E-1, 8.6602540378443864676E-1)
        """
        return self.sin(), self.cos()

    def acos(self):
        """
        Returns the inverse cosine of this number

        EXAMPLES:
            sage: RQDF = RealQuadDoubleField ()
            sage: q = RQDF.pi()/3
            sage: i = q.cos()
            sage: q
            1.0471975511965977462E0
            sage: i.acos()
            1.0471975511965977462E0

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
            sage: RQDF = RealQuadDoubleField ()
            sage: q = RQDF.pi()/3
            sage: i = q.sin()
            sage: q
            1.0471975511965977462E0
            sage: i.asin()
            1.0471975511965977462E0
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
            sage: RQDF = RealQuadDoubleField ()
            sage: q = RQDF.pi()/3
            sage: i = q.tan()
            sage: q
            1.0471975511965977462E0
            sage: i.atan()
            1.0471975511965977462E0

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
            sage: RQDF = RealQuadDoubleField ()
            sage: q = RQDF.pi()/12
            sage: q.cosh()
            1.0344656400955105653E0
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
            sage: RQDF = RealQuadDoubleField ()
            sage: q = -RQDF.pi()/12
            sage: q.sinh()
            -2.6480022760227075770E-1

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
            sage: RQDF = RealQuadDoubleField ()
            sage: q = RQDF.pi()/12
            sage: q.tanh()
            2.5597778924568453946E-1
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
            sage: RQDF = RealQuadDoubleField ()
            sage: q = RQDF.pi()/2
            sage: i = q.cosh() ; i
            2.5091784786580567821E0
            sage: i.acosh()
            1.5707963267948966193E0
            sage: q
            1.5707963267948966193E0
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
            sage: RQDF = RealQuadDoubleField ()
            sage: q = RQDF.pi()/2
            sage: i = q.sinh() ; i
            2.3012989023072948735E0
            sage: i.asinh() ; q
            1.5707963267948966193E0
            1.5707963267948966193E0
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
            sage: RQDF = RealQuadDoubleField ()
            sage: q = RQDF.pi()/2
            sage: i = q.tanh() ; i
            9.1715233566727434638E-1
            sage: i.atanh() ; q
            1.5707963267948966193E0
            1.5707963267948966193E0

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
