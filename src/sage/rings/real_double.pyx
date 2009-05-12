r"""
Double Precision Real Numbers

EXAMPLES:

We create the real double vector space of dimension `3`::

    sage: V = RDF^3; V
    Vector space of dimension 3 over Real Double Field

Notice that this space is unique.

::

    sage: V is RDF^3
    True
    sage: V is FreeModule(RDF, 3)
    True
    sage: V is VectorSpace(RDF, 3)
    True

Also, you can instantly create a space of large dimension.

::

    sage: V = RDF^10000
"""

include '../ext/python_float.pxi'
include '../ext/cdefs.pxi'
include '../ext/stdsage.pxi'
include '../ext/random.pxi'
include '../ext/interrupt.pxi'
include '../gsl/gsl.pxi'

gsl_set_error_handler_off()

import math, operator

cimport sage.libs.pari.gen
import sage.libs.pari.gen


import sage.rings.integer
import sage.rings.rational

from sage.rings.integer cimport Integer

from sage.categories.morphism cimport Morphism

def is_RealDoubleField(x):
    """
    Returns True if x is the field of real double precision numbers.

    EXAMPLE::

        sage: from sage.rings.real_double import is_RealDoubleField
        sage: is_RealDoubleField(RDF)
        True
        sage: is_RealDoubleField(RealField(53))
        False
    """
    return PY_TYPE_CHECK(x, RealDoubleField_class)

cdef class RealDoubleField_class(Field):
    """
    An approximation to the field of real numbers using double
    precision floating point numbers. Answers derived from calculations
    in this approximation may differ from what they would be if those
    calculations were performed in the true field of real numbers. This
    is due to the rounding errors inherent to finite precision
    calculations.

    EXAMPLES::

        sage: RR == RDF
        False
        sage: RDF == RealDoubleField()    # RDF is the shorthand
        True

    ::

        sage: RDF(1)
        1.0
        sage: RDF(2/3)
        0.666666666667

    A TypeError is raised if the coercion doesn't make sense::

        sage: RDF(QQ['x'].0)
        Traceback (most recent call last):
        ...
        TypeError: cannot coerce nonconstant polynomial to float
        sage: RDF(QQ['x'](3))
        3.0

    One can convert back and forth between double precision real
    numbers and higher-precision ones, though of course there may be
    loss of precision::

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
    def __init__(self):
        Field.__init__(self, self)
        self._populate_coercion_lists_(element_constructor=RealDoubleElement,
                                       init_no_parent=True,
                                       convert_method_name='_real_double_')

    def __reduce__(self):
        """
        EXAMPLES::

            sage: loads(dumps(RDF)) is RDF
            True
        """
        return RealDoubleField, ()

    cpdef bint is_exact(self) except -2:
        """
        Returns False, because doubles are not exact.

        EXAMPLE::

            sage: RDF.is_exact()
            False
        """
        return False

    def _latex_(self):
        return "\\Bold{R}"

    def _sage_input_(self, sib, coerced):
        r"""
        Produce an expression which will reproduce this value when evaluated.

        EXAMPLES:
            sage: sage_input(RDF, verify=True)
            # Verified
            RDF
            sage: from sage.misc.sage_input import SageInputBuilder
            sage: RDF._sage_input_(SageInputBuilder(), False)
            {atomic:RDF}
        """
        return sib.name('RDF')

    def __repr__(self):
        """
        Print out this real double field.

        EXAMPLES::

            sage: RealDoubleField()
            Real Double Field
            sage: RDF
            Real Double Field
        """
        return "Real Double Field"

    def __cmp__(self, x):
        """
        EXAMPLES::

            sage: RDF == 5
            False
            sage: loads(dumps(RDF)) == RDF
            True
        """
        if PY_TYPE_CHECK(x, RealDoubleField_class):
            return 0
        return cmp(type(self), type(x))

    def construction(self):
        """
        Returns the functorial construction of self, namely, completion of
        the rational numbers with respect to the prime at
        `\infty`.

        Also preserves other information that makes this field unique (i.e.
        the Real Double Field).

        EXAMPLES::

            sage: c, S = RDF.construction(); S
            Rational Field
            sage: RDF == c(S)
            True
        """
        from sage.categories.pushout import CompletionFunctor
        return (CompletionFunctor(sage.rings.infinity.Infinity,
                                  53,
                                  {'type': 'RDF'}),
               sage.rings.rational_field.QQ)

    def complex_field(self):
        """
        Returns the complex field with the same precision as self, ie, the
        complex double field.

        EXAMPLES::

            sage: RDF.complex_field()
            Complex Double Field
        """
        from sage.rings.complex_double import CDF
        return CDF

    def algebraic_closure(self):
        """
        Returns the algebraic closure of self, ie, the complex double
        field.

        EXAMPLES::

            sage: RDF.algebraic_closure()
            Complex Double Field
        """
        from sage.rings.complex_double import CDF
        return CDF

    cpdef _coerce_map_from_(self, S):
        """
        Canonical coercion of x to the real double field.

        The rings that canonically coerce to the real double field are:

        - the real double field itself

        - int, long, integer, and rational rings

        - real mathematical constants

        - the MPFR real field with <= 53 bits of precision

        EXAMPLES::

            sage: RDF.coerce(5)
            5.0
            sage: RDF.coerce(9499294r)
            9499294.0
            sage: RDF.coerce(61/3)
            20.3333333333
            sage: parent(RDF(3) + CDF(5))
            Complex Double Field
            sage: parent(CDF(5) + RDF(3))
            Complex Double Field
            sage: CDF.gen(0) + 5.0
            5.0 + 1.0*I
            sage: RLF(2/3) + RDF(1)
            1.66666666667
        """
        from integer_ring import ZZ
        from rational_field import QQ
        from real_lazy import RLF
        from real_mpfr import RR, RealField
        if S in [int, float, ZZ, QQ, RLF] or isinstance(S, RealField) and S.prec() >= 53:
            return ToRDF(S)
        connecting = RR.coerce_map_from(S)
        if connecting is not None:
            return ToRDF(RR) * connecting

    def _magma_init_(self, magma):
        r"""
        Return a string representation of self in the Magma language.

        EXAMPLES:

        Magma handles precision in decimal digits, so we lose a
        bit::

            sage: magma(RDF) # optional - magma
            Real field of precision 15
            sage: 10^15 < 2^53 < 10^16
            True

        When we convert back from Magma, we convert to a generic real field
        that has 53 bits of precision::

            sage: magma(RDF).sage() # optional - magma
            Real Field with 53 bits of precision
        """
        return "RealField(%s : Bits := true)" % self.prec()

    def prec(self):
        """
        Return the precision of this real double field in bits.

        Always returns 53.

        EXAMPLES::

            sage: RDF.prec()
            53
        """
        return 53

    def to_prec(self, prec):
        """
        Returns the real field to the specified precision. As doubles have
        fixed precision, this will only return a real double field if prec
        is exactly 53.

        EXAMPLES::

            sage: RDF.to_prec(52)
            Real Field with 52 bits of precision
            sage: RDF.to_prec(53)
            Real Double Field
        """
        if prec == 53:
            return self
        else:
            from real_mpfr import RealField_constructor
            return RealField_constructor(prec)


    def gen(self, n=0):
        """
        Return the generator of the real double field.

        EXAMPLES::

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
        Returns True, to signify that elements of this field print without
        sums, so parenthesis aren't required, e.g., in coefficients of
        polynomials.

        EXAMPLES::

            sage: RDF.is_atomic_repr()
            True
        """
        return True

    def is_finite(self):
        """
        Returns False, since the field of real numbers is not finite.
        Technical note: There exists an upper bound on the double
        representation.

        EXAMPLES::

            sage: RDF.is_finite()
            False
        """
        return False

    def characteristic(self):
        """
        Returns 0, since the field of real numbers has characteristic 0.

        EXAMPLES::

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
        Return a random element of this real double field in the interval
        [min, max].

        EXAMPLES::

            sage: RDF.random_element()
            0.736945423566
            sage: RDF.random_element(min=100, max=110)
            102.815947352
        """
        cdef randstate rstate = current_randstate()

        return self._new_c((max-min)*rstate.c_rand_double() + min)

    def name(self):
        return "RealDoubleField"

    def __hash__(self):
        """
        TEST::

            sage: hash(RDF) % 2^32 == hash(str(RDF)) % 2^32
            True
        """
        return 1157042230 #return hash(str(self))

    def pi(self):
        """
        Returns pi to double-precision.

        EXAMPLES::

            sage: RDF.pi()
            3.14159265359
            sage: RDF.pi().sqrt()/2
            0.886226925453
        """
        return self(M_PI)


    def euler_constant(self):
        """
        Returns Euler's gamma constant to double precision

        EXAMPLES::

            sage: RDF.euler_constant()
            0.577215664902
        """
        return self(M_EULER)

    def log2(self):
        """
        Returns log(2) to the precision of this field.

        EXAMPLES::

            sage: RDF.log2()
            0.69314718056
            sage: RDF(2).log()
            0.69314718056
        """
        return self(M_LN2)

    def factorial(self, int n):
        """
        Return the factorial of the integer n as a real number.

        EXAMPLES::

            sage: RDF.factorial(100)
            9.33262154439e+157
        """
        if n < 0:
            raise ArithmeticError, "n must be nonnegative"
        return self(gsl_sf_fact(n))

    def zeta(self, n=2):
        """
        Return an `n`-th root of unity in the real field, if one
        exists, or raise a ValueError otherwise.

        EXAMPLES::

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

    def NaN(self):
        """
        EXAMPLES::

            sage: RDF.NaN()
            NaN
        """
        return self(0)/self(0)

    def nan(self):
        """
        EXAMPLES::

            sage: RDF.nan()
            NaN
        """
        return self(0)/self(0)

cdef class RealDoubleElement(FieldElement):
    """
    An approximation to a real number using double precision floating
    point numbers. Answers derived from calculations with such
    approximations may differ from what they would be if those
    calculations were performed with true real numbers. This is due to
    the rounding errors inherent to finite precision calculations.
    """
    def __new__(self, x=None):
        (<Element>self)._parent = _RDF

    def __init__(self, x):
        """
        Create a new RealDoubleElement with value x.

        EXAMPLES::

            sage: RDF(10^100)
            1e+100
        """
        self._value = float(x)

    def _magma_init_(self, magma):
        r"""
        Return a string representation of self in the Magma language.

        EXAMPLES::

            sage: RDF(10.5)
            10.5
            sage: magma(RDF(10.5)) # optional - magma
            10.5000000000000
        """
        return "%s!%s" % (self.parent()._magma_init_(magma), self)

    def __reduce__(self):
        """
        EXAMPLES::

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

    def prec(self):
        """
        Returns the precision of this number in bits.

        Always returns 53.

        EXAMPLES::

            sage: RDF(0).prec()
            53
        """

        return 53

    def ulp(self):
        """
        Returns the unit of least precision of self, which is the weight of
        the least significant bit of self. Unless self is exactly a power of
        two, it is gap between this number and the next closest distinct
        number that can be represented.

        EXAMPLES::

            sage: a = RDF(1)
            sage: a - a.ulp() == a
            False
            sage: a - a.ulp()/2 == a
            True

            sage: a = RDF.pi()
            sage: b = a + a.ulp()
            sage: (a+b)/2 in [a,b]
            True

            sage: a = RDF(1)/RDF(0); a
            +infinity
            sage: a.ulp()
            +infinity
            sage: (-a).ulp()
            +infinity
            sage: a = RR('nan')
            sage: a.ulp() is a
            True
        """
        cdef int e, v = gsl_isinf(self._value)
        if gsl_isnan(self._value):
            return self
        elif self._value == 0:
            return RealDoubleElement(ldexp(1.0, e-1082))
        elif v == 1:
            return self
        elif v == -1:
            return -self
        else:
            frexp(self._value, &e)
            return RealDoubleElement(ldexp(1.0, e-54))

    def real(self):
        """
        Returns itself - we're already real.

        EXAMPLES::

            sage: a = RDF(3)
            sage: a.real()
            3.0
        """
        return self

    def imag(self):
        """
        Returns the imaginary part of this number. (hint: it's zero.)

        EXAMPLES::

            sage: a = RDF(3)
            sage: a.imag()
            0.0
        """
        return RealDoubleElement(0)

    def __complex__(self):
        """
        EXAMPLES::

            sage: a = 2303
            sage: RDF(a)
            2303.0
            sage: complex(RDF(a))
            (2303+0j)
        """
        return complex(self._value,0)

    def _integer_(self, ZZ=None):
        """
        If this floating-point number is actually an integer, return
        that integer.  Otherwise, raise an exception.

        EXAMPLES::

            sage: ZZ(RDF(237.0))
            237
            sage: ZZ(RDF(0.0/0.0))
            Traceback (most recent call last):
            ...
            TypeError: Cannot convert non-integral float to integer
            sage: ZZ(RDF(1.0/0.0))
            Traceback (most recent call last):
            ...
            OverflowError: cannot convert float infinity to long
            sage: ZZ(RDF(-123456789.0))
            -123456789
            sage: ZZ(RDF((2.0))^290)
            1989292945639146568621528992587283360401824603189390869761855907572637988050133502132224
            sage: ZZ(RDF(-2345.67))
            Traceback (most recent call last):
            ...
            TypeError: Cannot convert non-integral float to integer
        """
        return Integer(self._value)

    def parent(self):
        """
        Return the real double field, which is the parent of self.

        EXAMPLES::

            sage: a = RDF(2.3)
            sage: a.parent()
            Real Double Field
            sage: parent(a)
            Real Double Field
        """
        return self._parent

    def _interface_init_(self):
        """
        Returns self formatted as a string, suitable as input to another
        computer algebra system. (This is the default function used for
        exporting to other computer algebra systems.)

        EXAMPLES::

            sage: s1 = RDF(sin(1)); s1
            0.841470984808
            sage: s1._interface_init_()
            '0.8414709848078965'
            sage: s1 == RDF(gp(s1))
            True
        """
        return repr(self._value)

    def _sage_input_(self, sib, coerced):
        r"""
        Produce an expression which will reproduce this value when evaluated.

        EXAMPLES::

            sage: sage_input(RDF(NaN), verify=True)
            # Verified
            RDF(NaN)
            sage: sage_input(RDF(-infinity), verify=True)
            # Verified
            -RDF(infinity)
            sage: sage_input(RDF(-infinity)*polygen(RDF), verify=True)
            # Verified
            R.<x> = RDF[]
            -RDF(infinity)*x + RDF(NaN)
            sage: sage_input(RDF(pi), verify=True)
            # Verified
            RDF(3.1415926535897931)
            sage: sage_input(RDF(-e), verify=True, preparse=False)
            # Verified
            -RDF(2.718281828459045...)
            sage: sage_input(RDF(pi)*polygen(RDF), verify=True, preparse=None)
            # Verified
            R = RDF['x']
            x = R.gen()
            3.1415926535897931*x
            sage: from sage.misc.sage_input import SageInputBuilder
            sage: sib = SageInputBuilder()
            sage: RDF(22/7)._sage_input_(sib, True)
            {atomic:3.1428571428571428}
            sage: RDF(22/7)._sage_input_(sib, False)
            {call: {atomic:RDF}({atomic:3.1428571428571428})}
        """
        cdef int isinf = gsl_isinf(self._value)
        cdef bint isnan = gsl_isnan(self._value)
        if isinf or isnan:
            if isnan:
                v = sib.name('NaN')
            else:
                v = sib.name('infinity')
            v = sib(self.parent())(v)
            if isinf < 0:
                v = -v
            return v

        from sage.rings.all import ZZ, RR

        cdef bint negative = self._value < 0
        if negative:
            self = -self

        # There are five possibilities for printing this floating-point
        # number, ordered from prettiest to ugliest (IMHO).
        # 1) An integer: 42
        # 2) A simple literal: 3.14159
        # 3) A coerced integer: RDF(42)
        # 4) A coerced literal: RDF(3.14159)
        # 5) A coerced RR value: RDF(RR('3.14159'))

        # str(self) works via libc, which we don't necessarily trust
        # to produce the best possible answer.  So this function prints
        # via RR/MPFR.  Without the preparser, input works via libc as
        # well, but we don't have a choice about that.

        # To use choice 1 or choice 3, this number must be an integer.
        cdef bint can_use_int_literal = \
            self.abs() < (Integer(1) << self.prec()) and self in ZZ

        self_str = RR(self._value).str(truncate=False, skip_zeroes=True)

        # To use choice 2 or choice 4, we must be able to read
        # numbers of this precision as a literal.
        cdef bint can_use_float_literal = \
            (sib.preparse() or float(self_str) == self)

        if can_use_int_literal or can_use_float_literal:
            if can_use_int_literal:
                v = sib.int(self._integer_())
            else:
                v = sib.float_str(self_str)
        else:
            v = sib(RR(self))
        if not coerced:
            v = sib(self.parent())(v)

        if negative:
            v = -v

        return v

    def __repr__(self):
        """
        Return print version of self.

        EXAMPLES::

            sage: a = RDF(2); a
            2.0
            sage: a^2
            4.0
        """
        return self.str()

    def _latex_(self):
        """
        EXAMPLE::

            sage: RDF(2e-100)._latex_()
            '2 \\times 10^{-100}'
        """
        s = self.str()
        parts = s.split('e')
        if len(parts) > 1:
            # scientific notation
            if parts[1][0] == '+':
                parts[1] = parts[1][1:]
            s = "%s \\times 10^{%s}" % (parts[0], parts[1])
        return s

    def __hash__(self):
        """
        Returns the hash of self, which coincides with the python float
        (and often int) type.

        EXAMPLE::

            sage: hash(RDF(1.2)) == hash(1.2r)
            True
        """
        return hash(self._value)

    def _im_gens_(self, codomain, im_gens):
        return codomain(self) # since 1 |--> 1

    def str(self):
        """
        Return string representation of self.

        EXAMPLES::

            sage: a = RDF('4.5'); a.str()
            '4.5'
            sage: a = RDF('49203480923840.2923904823048'); a.str()
            '4.92034809238e+13'
            sage: a = RDF(1)/RDF(0); a.str()
            '+infinity'
            sage: a = -RDF(1)/RDF(0); a.str()
            '-infinity'
            sage: a = RDF(0)/RDF(0); a.str()
            'NaN'

        We verify consistency with RR (mpfr reals)::

            sage: str(RR(RDF(1)/RDF(0))) == str(RDF(1)/RDF(0))
            True
            sage: str(RR(-RDF(1)/RDF(0))) == str(-RDF(1)/RDF(0))
            True
            sage: str(RR(RDF(0)/RDF(0))) == str(RDF(0)/RDF(0))
            True
        """
        cdef int v
        if gsl_isnan(self._value):
            return "NaN"
        else:
            v = gsl_isinf(self._value)
            if v == 1:
                return "+infinity"
            elif v == -1:
                return "-infinity"
        return str(self._value)

    def __copy__(self):
        """
        Return copy of self, which since self is immutable, is just self.

        EXAMPLES::

            sage: r = RDF('-1.6')
            sage: r.__copy__() is r
            True
        """
        return self

    def integer_part(self):
        """
        If in decimal this number is written n.defg, returns n.

        EXAMPLES::

            sage: r = RDF('-1.6')
            sage: a = r.integer_part(); a
            -1
            sage: type(a)
            <type 'sage.rings.integer.Integer'>
            sage: r = RDF(0.0/0.0)
            sage: a = r.integer_part()
            Traceback (most recent call last):
            ...
            TypeError: Attempt to get integer part of NaN
        """
        if gsl_isnan(self._value):
            raise TypeError, "Attempt to get integer part of NaN"
        else:
            return Integer(int(self._value))


    ########################
    #   Basic Arithmetic
    ########################
    def __invert__(self):
        """
        Compute the multiplicative inverse of self.

        EXAMPLES::

            sage: a = RDF(-1.5)*RDF(2.5)
            sage: a.__invert__()
            -0.266666666667
            sage: ~a
            -0.266666666667
        """
        cdef RealDoubleElement x = <RealDoubleElement>PY_NEW(RealDoubleElement)
        x._value = 1.0 / self._value
        return x

    cpdef ModuleElement _add_(self, ModuleElement right):
        """
        Add two real numbers with the same parent.

        EXAMPLES::

            sage: RDF('-1.5') + RDF('2.5')
            1.0
        """
        cdef RealDoubleElement x = <RealDoubleElement>PY_NEW(RealDoubleElement)
        x._value = self._value + (<RealDoubleElement>right)._value
        return x

    cpdef ModuleElement _iadd_(self, ModuleElement right):
        """
        EXAMPLES::

            sage: a = RDF(0.5)
            sage: a += RDF(3); a
            3.5
        """
        # self and right are guaranteed to be Integers
        self._value += (<RealDoubleElement>right)._value
        return self

    cpdef ModuleElement _sub_(self, ModuleElement right):
        """
        Subtract two real numbers with the same parent.

        EXAMPLES::

            sage: RDF('-1.5') - RDF('2.5')
            -4.0
        """
        cdef RealDoubleElement x = <RealDoubleElement>PY_NEW(RealDoubleElement)
        x._value = self._value - (<RealDoubleElement>right)._value
        return x

    cpdef ModuleElement _isub_(self, ModuleElement right):
        """
        EXAMPLES::

            sage: a = RDF(0.5)
            sage: a -= RDF(3); a
            -2.5
        """
        self._value -= (<RealDoubleElement>right)._value
        return self

    cpdef RingElement _mul_(self, RingElement right):
        """
        Multiply two real numbers with the same parent.

        EXAMPLES::

            sage: RDF('-1.5') * RDF('2.5')
            -3.75
        """
        cdef RealDoubleElement x = <RealDoubleElement>PY_NEW(RealDoubleElement)
        x._value = self._value * (<RealDoubleElement>right)._value
        return x

    cpdef RingElement _imul_(self, RingElement right):
        """
        EXAMPLES::

            sage: a = RDF(2.5)
            sage: a *= RDF(3); a
            7.5
        """
        self._value *= (<RealDoubleElement>right)._value
        return self

    cpdef RingElement _div_(self, RingElement right):
        """
        EXAMPLES::

            sage: RDF('-1.5') / RDF('2.5')
            -0.6
            sage: RDF(1)/RDF(0)
            +infinity
        """
        cdef RealDoubleElement x = <RealDoubleElement>PY_NEW(RealDoubleElement)
        x._value = self._value / (<RealDoubleElement>right)._value
        return x

    cpdef RingElement _idiv_(self, RingElement right):
        """
        EXAMPLES::

            sage: a = RDF(1.5)
            sage: a /= RDF(2); a
            0.75
            sage: a /= RDF(0); a
            +infinity
        """
        self._value /= (<RealDoubleElement>right)._value
        return self

    def __neg__(self):
        """
        Negates a real number.

        EXAMPLES::

            sage: -RDF('-1.5')
            1.5
        """
        cdef RealDoubleElement x = <RealDoubleElement>PY_NEW(RealDoubleElement)
        x._value = -self._value
        return x

    def __abs__(self):
        """
        Returns the absolute value of self.

        EXAMPLES::

            sage: abs(RDF(1.5))
            1.5
            sage: abs(RDF(-1.5))
            1.5
        """
        if self._value >= 0:
            return self
        else:
            return self._new_c(-self._value)

    cpdef RealDoubleElement abs(RealDoubleElement self):
        """
        Returns the absolute value of self.

        EXAMPLES::

            sage: RDF(1e10).abs()
            10000000000.0
            sage: RDF(-1e10).abs()
            10000000000.0
        """
        if self._value >= 0:
            return self
        else:
            return self._new_c(-self._value)

    def __lshift__(x, y):
        """
        LShifting a double is not supported; nor is lshifting a
        RealDoubleElement.
        """
        raise TypeError, "unsupported operand type(s) for <<"

    def __rshift__(x, y):
        """
        RShifting a double is not supported; nor is rshifting a
        RealDoubleElement.
        """
        raise TypeError, "unsupported operand type(s) for >>"

    def multiplicative_order(self):
        r"""
        Returns `n` such that self^n == 1.

        Only `\pm 1` have finite multiplicative order.

        EXAMPLES::

            sage: RDF(1).multiplicative_order()
            1
            sage: RDF(-1).multiplicative_order()
            2
            sage: RDF(3).multiplicative_order()
            +Infinity
        """
        if self._value == 1:
            return 1
        elif self._value == -1:
            return 2
        return sage.rings.infinity.infinity

    def sign(self):
        """
        Returns -1,0, or 1 if self is negative, zero, or positive;
        respectively.

        EXAMPLES::

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
        Given real number x, rounds up if fractional part is greater than
        .5, rounds down if fractional part is less than .5.

        EXAMPLES::

            sage: RDF(0.49).round()
            0
            sage: a=RDF(0.51).round(); a
            1
        """
        return Integer(round(self._value))

    def floor(self):
        """
        Returns the floor of this number

        EXAMPLES::

            sage: RDF(2.99).floor()
            2
            sage: RDF(2.00).floor()
            2
            sage: RDF(-5/2).floor()
            -3
        """
        return Integer(math.floor(self._value))

    def ceil(self):
        """
        Returns the ceiling of this number

        OUTPUT: integer

        EXAMPLES::

            sage: RDF(2.99).ceil()
            3
            sage: RDF(2.00).ceil()
            2
            sage: RDF(-5/2).ceil()
            -2
        """
        return Integer(math.ceil(self._value))

    ceiling = ceil

    def trunc(self):
        """
        Truncates this number (returns integer part).

        EXAMPLES::

            sage: RDF(2.99).trunc()
            2
            sage: RDF(-2.00).trunc()
            -2
            sage: RDF(0.00).trunc()
            0
        """
        return Integer(int(self._value))

    def frac(self):
        """
        frac returns a real number > -1 and < 1. it satisfies the relation:
        x = x.trunc() + x.frac()

        EXAMPLES::

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
        """
        Return self as a python float.

        EXAMPLES::

            sage: float(RDF(1.5))
            1.5
            sage: type(float(RDF(1.5)))
            <type 'float'>
        """
        return self._value

    def _rpy_(self):
        """
        Returns self.__float__() for rpy to convert into the
        appropriate R object.

        EXAMPLES::

            sage: n = RDF(2.0)
            sage: n._rpy_()
            2.0
            sage: type(n._rpy_())
            <type 'float'>
        """
        return self.__float__()

    def __int__(self):
        """
        Returns integer truncation of this real number.

        EXAMPLES::

            sage: int(RDF(2.99))
            2
            sage: int(RDF(-2.99))
            -2
        """
        return int(self._value)

    def __long__(self):
        """
        Returns long integer truncation of this real number.

        EXAMPLES::

            sage: int(RDF(10e15))
            10000000000000000L                   # 32-bit
            10000000000000000                    # 64-bit
            sage: long(RDF(2^100)) == 2^100
            True
        """
        return long(self._value)

    def _complex_mpfr_field_(self, CC):
        """
        EXAMPLES::

            sage: a = RDF(1/3)
            sage: CC(a)
            0.333333333333333
            sage: a._complex_mpfr_field_(CC)
            0.333333333333333

        If we coerce to a higher-precision field the extra bits appear
        random; they are actualy 0's in base 2.

        ::

            sage: a._complex_mpfr_field_(ComplexField(100))
            0.33333333333333331482961625625
            sage: a._complex_mpfr_field_(ComplexField(100)).str(2)
            '0.01010101010101010101010101010101010101010101010101010100000000000000000000000000000000000000000000000'
        """
        return CC(self._value)

    def _complex_double_(self, CDF):
        """
        EXAMPLES::

            sage: CDF(RDF(1/3))
            0.333333333333
        """
        return CDF(self._value)

    def _pari_(self):
        """
        EXAMPLES::

            sage: RDF(1.5)._pari_()
            1.50000000000000
        """
        cdef sage.libs.pari.gen.PariInstance P = sage.libs.pari.gen.pari
        return P.double_to_gen_c(self._value)


    ###########################################
    # Comparisons: ==, !=, <, <=, >, >=
    ###########################################

    def is_NaN(self):
        """
        EXAMPLES::

            sage: RDF(1).is_NaN()
            False
            sage: a = RDF(0)/RDF(0)
            sage: a.is_NaN()
            True
        """
        return gsl_isnan(self._value)

    def is_positive_infinity(self):
        """
        EXAMPLES::

            sage: a = RDF(1)/RDF(0)
            sage: a.is_positive_infinity()
            True
            sage: a = RDF(-1)/RDF(0)
            sage: a.is_positive_infinity()
            False
        """
        return gsl_isinf(self._value) > 0

    def is_negative_infinity(self):
        """
        EXAMPLES::

            sage: a = RDF(2)/RDF(0)
            sage: a.is_negative_infinity()
            False
            sage: a = RDF(-3)/RDF(0)
            sage: a.is_negative_infinity()
            True
        """
        return gsl_isinf(self._value) < 0

    def is_infinity(self):
        """
        EXAMPLES::

            sage: a = RDF(2); b = RDF(0)
            sage: (a/b).is_infinity()
            True
            sage: (b/a).is_infinity()
            False
        """
        return gsl_isinf(self._value)

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

    def NaN(self):
        """
        EXAMPLES::

            sage: RDF.NaN()
            NaN
        """
        return self(0)/self(0)

    def nan(self):
        """
        EXAMPLES::


            sage: RDF.nan()
            NaN
        """
        return self(0)/self(0)

    def sqrt(self, extend=True, all=False):
        """
        The square root function.

        INPUT:


        -  ``extend`` - bool (default: True); if True, return a
           square root in a complex field if necessary if self is negative;
           otherwise raise a ValueError

        -  ``all`` - bool (default: False); if True, return a
           list of all square roots.


        EXAMPLES::

            sage: r = RDF(4.0)
            sage: r.sqrt()
            2.0
            sage: r.sqrt()^2 == r
            True

        ::

            sage: r = RDF(4344)
            sage: r.sqrt()
            65.9090282131
            sage: r.sqrt()^2 - r
            0.0

        ::

            sage: r = RDF(-2.0)
            sage: r.sqrt()
            1.41421356237*I

        ::

            sage: RDF(2).sqrt(all=True)
            [1.41421356237, -1.41421356237]
            sage: RDF(0).sqrt(all=True)
            [0.0]
            sage: RDF(-2).sqrt(all=True)
            [1.41421356237*I, -1.41421356237*I]
        """
        if self._value >= 0:
            x = self._new_c(sqrt(self._value))
            if all:
                if x.is_zero():
                    return [x]
                else:
                    return [x, -x]
            else:
                return x
        if not extend:
            raise ValueError, "negative number %s does not have a square root in the real field"%self
        import sage.rings.complex_double
        return self._complex_double_(sage.rings.complex_double.CDF).sqrt(all=all)

    def is_square(self):
        """
        Returns whether or not this number is a square in this field. For
        the real numbers, this is True if and only if self is
        non-negative.

        EXAMPLES::

            sage: RDF(3.5).is_square()
            True
            sage: RDF(0).is_square()
            True
            sage: RDF(-4).is_square()
            False
        """
        return self._value >= 0

    def cube_root(self):
        """
        Return the cubic root (defined over the real numbers) of self.

        EXAMPLES::

            sage: r = RDF(125.0); r.cube_root()
            5.0
            sage: r = RDF(-119.0)
            sage: r.cube_root()^3 - r         # output is random, depending on arch.
            0.0
        """
        return self.nth_root(3)


    def nth_root(self, int n):
        """
        Returns the `n^{th}` root of self.

        INPUT:


        -  ``n`` - an integer


        OUTPUT: an real or complex double

        The output is complex if self is negative and n is even.

        EXAMPLES::

            sage: r = RDF(-125.0); r.nth_root(3)
            -5.0
            sage: r.nth_root(5)
            -2.6265278044
            sage: RDF(-2).nth_root(5)^5
            -2.0
            sage: RDF(-1).nth_root(5)^5
            -1.0
            sage: RDF(3).nth_root(10)^10
            3.0
            sage: RDF(-1).nth_root(2)
            6.12323399574e-17 + 1.0*I
            sage: RDF(-1).nth_root(4)
            0.707106781187 + 0.707106781187*I
        """
        if n == 0:
            return RealDoubleElement(float('nan'))
        if self._value < 0:
            if GSL_IS_EVEN(n):
                return self._complex_double_(sage.rings.complex_double.CDF).nth_root(n)
            else:
                return - ( (-self).__pow__(float(1)/n) )
        else:
            return self.__pow__(float(1)/n)

    cdef RealDoubleElement __pow_float(self, double exponent):
        """
        TESTS:

              sage: RDF(0)^.5
              0.0
              sage: RDF(0)^(1/2)
              0.0
              sage: RDF(0)^RDF(0)
              1.0
        """
        if exponent == 0:
            return self._new_c(1)
        elif self._value == 0 or self._value == 1:
            return self
        else:
            return self._new_c(gsl_sf_exp(gsl_sf_log(self._value) * exponent))

    cdef RealDoubleElement __pow_int(self, int exponent):
        return self._new_c(gsl_pow_int(self._value, exponent))

    def __pow__(self, exponent, modulus):
        """
        Compute self raised to the power of exponent, rounded in the
        direction specified by the parent of self.

        If the result is not a real number, self and the exponent are both
        coerced to complex numbers (with sufficient precision), then the
        exponentiation is computed in the complex numbers. Thus this
        function can return either a real or complex number.

        EXAMPLES::

            sage: a = RDF('1.23456')
            sage: a^20
            67.6462977039
            sage: a^a
            1.29711148178

        Symbolic examples::

            sage: x, y = var('x,y')
            sage: RDF('-2.3')^(x+y^3+sin(x))
            (-2.3)^(y^3 + sin(x) + x)
            sage: RDF('-2.3')^x
            (-2.3)^x
        """
        cdef RealDoubleElement base, exp
        if PY_TYPE_CHECK(self, RealDoubleElement):
            base = self
            if PY_TYPE_CHECK(exponent, RealDoubleElement):
                return base.__pow_float((<RealDoubleElement>exponent)._value)
            elif PY_TYPE_CHECK(exponent, float):
                return base.__pow_float(exponent)
            elif PY_TYPE_CHECK(exponent, int):
                return base.__pow_int(exponent)
            elif PY_TYPE_CHECK(exponent, Integer) and exponent < INT_MAX:
                return base.__pow_int(exponent)
            try:
                exp = base._parent(exponent)
                return base.__pow_float(exp._value)
            except TypeError:
                return exponent.parent()(self) ** exponent # neither operand is RealDoubleElement
        else:
            try:
                base = exp.parent()(self)
                return base.__pow_float((<RealDoubleElement>exponent)._value)
            except TypeError:
                return self ** self.parent()(exponent) # neither operand is RealDoubleElement


    cdef _log_base(self, double log_of_base):
        if self._value < 2:
            if self._value == 0:
                return RDF(-1)/RDF(0)
            if self._value < 0:
                return RDF(0)/RDF(0)
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
        EXAMPLES::

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
            -infinity
            sage: RDF(-1).log()
            NaN
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

        EXAMPLES::

            sage: r = RDF(16.0)
            sage: r.log2()
            4.0

        ::

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

        EXAMPLES::

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

        EXAMPLES::

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
        Returns `e^\mathtt{self}`

        EXAMPLES::

            sage: r = RDF(0.0)
            sage: r.exp()
            1.0

        ::

            sage: r = RDF('32.3')
            sage: a = r.exp(); a
            1.06588847275e+14
            sage: a.log()
            32.3

        ::

            sage: r = RDF('-32.3')
            sage: r.exp()
            9.3818445885e-15

        ::

            sage: RDF(1000).exp()
            +infinity
        """
        _sig_on
        a = self._new_c(gsl_sf_exp(self._value))
        _sig_off
        return a

    def exp2(self):
        """
        Returns `2^\mathtt{self}`

        EXAMPLES::

            sage: r = RDF(0.0)
            sage: r.exp2()
            1.0

        ::

            sage: r = RDF(32.0)
            sage: r.exp2()
            4294967296.0

        ::

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
        Returns `10^\mathtt{self}`

        EXAMPLES::

            sage: r = RDF(0.0)
            sage: r.exp10()
            1.0

        ::

            sage: r = RDF(32.0)
            sage: r.exp10()
            1e+32

        ::

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

        EXAMPLES::

            sage: t=RDF.pi()/2
            sage: t.cos()
            6.12323399574e-17
        """
        return self._new_c(gsl_sf_cos(self._value))

    def sin(self):
        """
        Returns the sine of this number

        EXAMPLES::

            sage: RDF(2).sin()
            0.909297426826
        """
        return self._new_c(gsl_sf_sin(self._value))

    def restrict_angle(self):
        r"""
        Returns a number congruent to self mod `2\pi` that lies in
        the interval `(-\pi, \pi]`.

        Specifically, it is the unique `x \in (-\pi, \pi]` such
        that ```self`` = x + 2\pi n` for some
        `n \in \ZZ`.

        EXAMPLES::

            sage: RDF(pi).restrict_angle()
            3.14159265359
            sage: RDF(pi + 1e-10).restrict_angle()
            -3.14159265349
            sage: RDF(1+10^10*pi).restrict_angle()
            0.9999977606...
        """
        return self._new_c(gsl_sf_angle_restrict_symm(self._value))

    def tan(self):
        """
        Returns the tangent of this number

        EXAMPLES::

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

        EXAMPLES::

            sage: t = RDF.pi()/6
            sage: t.sincos()
            (0.5, 0.866025403784)
        """
        return self.sin(), self.cos()

    def hypot(self, other):
        r"""
        Computes the value `\sqrt(self^2 + other^2)` in such a way
        as to avoid overflow.

        EXAMPLES::

            sage: x = RDF(4e300); y = RDF(3e300);
            sage: x.hypot(y)
            5e+300
            sage: sqrt(x^2+y^2) # overflow
            +infinity
        """
        _sig_on
        a = self._new_c(gsl_sf_hypot(self._value, float(other)))
        _sig_off
        return a

    def arccos(self):
        """
        Returns the inverse cosine of this number

        EXAMPLES::

            sage: q = RDF.pi()/3
            sage: i = q.cos()
            sage: i.arccos() == q
            True
        """
        return self._new_c(acos(self._value))

    def arcsin(self):
        """
        Returns the inverse sine of this number

        EXAMPLES::

            sage: q = RDF.pi()/5
            sage: i = q.sin()
            sage: i.arcsin() == q
            True
        """
        return self._new_c(asin(self._value))

    def arctan(self):
        """
        Returns the inverse tangent of this number

        EXAMPLES::

            sage: q = RDF.pi()/5
            sage: i = q.tan()
            sage: i.arctan() == q
            True
        """
        return self._new_c(atan(self._value))


    def cosh(self):
        """
        Returns the hyperbolic cosine of this number

        EXAMPLES::

            sage: q = RDF.pi()/12
            sage: q.cosh()
            1.0344656401
        """
        return self._new_c(gsl_ldexp( gsl_sf_exp(self._value) + gsl_sf_exp(-self._value), -1)) # (e^x + x^-x)/2

    def sinh(self):
        """
        Returns the hyperbolic sine of this number

        EXAMPLES::

            sage: q = RDF.pi()/12
            sage: q.sinh()
            0.264800227602
        """
        return self._new_c(gsl_ldexp( gsl_sf_expm1(self._value) - gsl_sf_expm1(-self._value), -1)) # (e^x - x^-x)/2

    def tanh(self):
        """
        Returns the hyperbolic tangent of this number

        EXAMPLES::

            sage: q = RDF.pi()/12
            sage: q.tanh()
            0.255977789246
        """
        return self.sinh() / self.cosh()

    def acosh(self):
        """
        Returns the hyperbolic inverse cosine of this number

        EXAMPLES::

            sage: q = RDF.pi()/2
            sage: i = q.cosh() ; i
            2.50917847866
            sage: abs(i.acosh()-q) < 1e-15
            True
        """
        return self._new_c(gsl_acosh(self._value))

    def arcsinh(self):
        """
        Returns the hyperbolic inverse sine of this number

        EXAMPLES::

            sage: q = RDF.pi()/2
            sage: i = q.sinh() ; i
            2.30129890231
            sage: abs(i.arcsinh()-q) < 1e-15
            True
        """
        return self._new_c(gsl_asinh(self._value))

    def arctanh(self):
        """
        Returns the hyperbolic inverse tangent of this number

        EXAMPLES::

            sage: q = RDF.pi()/2
            sage: i = q.tanh() ; i
            0.917152335667
            sage: i.arctanh() - q      # output is random, depending on arch.
            -4.4408920985e-16
        """
        return self._new_c(gsl_atanh(self._value))

    def sech(self):
        r"""
        This function returns the hyperbolic secant.

        EXAMPLES::

            sage: RDF(pi).sech()
            0.0862667383341
            sage: CDF(pi).sech()
            0.0862667383341
        """
        return 1/self.cosh()

    def csch(self):
        r"""
        This function returns the hyperbolic cosecant.

        EXAMPLES::

            sage: RDF(pi).csch()
            0.08658953753
            sage: CDF(pi).csch()
            0.08658953753
        """
        return 1/self.sinh()

    def coth(self):
        r"""
        This function returns the hyperbolic cotangent.

        EXAMPLES::

            sage: RDF(pi).coth()
            1.0037418732
            sage: CDF(pi).coth()
            1.0037418732
        """
        return self.cosh() / self.sinh()

    def agm(self, other):
        """
        Return the arithmetic-geometric mean of self and other. The
        arithmetic-geometric mean is the common limit of the sequences
        `u_n` and `v_n`, where `u_0` is self,
        `v_0` is other, `u_{n+1}` is the arithmetic mean
        of `u_n` and `v_n`, and `v_{n+1}` is the
        geometric mean of u_n and v_n. If any operand is negative, the
        return value is ``NaN``.

        EXAMPLES::

            sage: a = RDF(1.5)
            sage: b = RDF(2.3)
            sage: a.agm(b)
            1.87864845581

        The arithmetic-geometric mean always lies between the geometric and
        arithmetic mean.

        ::

            sage: sqrt(a*b) < a.agm(b) < (a+b)/2
            True
        """
        return RealDoubleElement(sage.rings.all.RR(self).agm(sage.rings.all.RR(other)))

    def erf(self):
        """
        Returns the value of the error function on self.

        EXAMPLES::

            sage: RDF(6).erf()
            1.0
        """
        return self._new_c(gsl_sf_erf(self._value))

    def gamma(self):
        """
        The Euler gamma function. Return gamma of self.

        EXAMPLES::

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

        .. note::

           PARI is vastly more efficient at computing the Riemann zeta
           function. See the example below for how to use it.

        EXAMPLES::

            sage: RDF(2).zeta()
            1.64493406685
            sage: RDF.pi()^2/6
            1.64493406685
            sage: RDF(-2).zeta()       # slightly random-ish arch dependent output
            -2.37378795339e-18
            sage: RDF(1).zeta()
            +infinity
        """
        if self._value == 1:
            return self._new_c(1)/self._new_c(0)
        return self._new_c(gsl_sf_zeta(self._value))

    def algdep(self, n):
        """
        Returns a polynomial of degree at most `n` which is
        approximately satisfied by this number. Note that the returned
        polynomial need not be irreducible, and indeed usually won't be if
        this number is a good approximation to an algebraic number of
        degree less than `n`.

        ALGORITHM: Uses the PARI C-library algdep command.

        EXAMPLE::

            sage: r = RDF(2).sqrt(); r
            1.41421356237
            sage: r.algdep(5)
            x^2 - 2
        """
        return sage.rings.arith.algdep(self,n)

    def algebraic_dependency(self, n):
        """
        Returns a polynomial of degree at most `n` which is
        approximately satisfied by this number. Note that the returned
        polynomial need not be irreducible, and indeed usually won't be if
        this number is a good approximation to an algebraic number of
        degree less than `n`.

        ALGORITHM: Uses the PARI C-library algdep command.

        EXAMPLE::

            sage: r = sqrt(RDF(2)); r
            1.41421356237
            sage: r.algdep(5)
            x^2 - 2
        """
        return sage.rings.arith.algdep(self,n)

cdef class ToRDF(Morphism):
    def __init__(self, R):
        """
        Fast morphism from anything with a __float__ method to an RDF
        element.

        EXAMPLES::

            sage: f = RDF.coerce_map_from(ZZ); f
            Native morphism:
              From: Integer Ring
              To:   Real Double Field
            sage: f(4)
            4.0
            sage: f = RDF.coerce_map_from(QQ); f
            Native morphism:
              From: Rational Field
              To:   Real Double Field
            sage: f(1/2)
            0.5
            sage: f = RDF.coerce_map_from(int); f
            Native morphism:
              From: Set of Python objects of type 'int'
              To:   Real Double Field
            sage: f(3r)
            3.0
            sage: f = RDF.coerce_map_from(float); f
            Native morphism:
              From: Set of Python objects of type 'float'
              To:   Real Double Field
            sage: f(3.5)
            3.5
        """
        from sage.categories.homset import Hom
        if isinstance(R, type):
            from sage.structure.parent import Set_PythonType
            R = Set_PythonType(R)
        Morphism.__init__(self, Hom(R, RDF))
    cpdef Element _call_(self, x):
        cdef RealDoubleElement r = <RealDoubleElement>PY_NEW(RealDoubleElement)
        r._value = PyFloat_AsDouble(x)
        return r
    def _repr_type(self):
        return "Native"


#####################################################
# unique objects
#####################################################
cdef RealDoubleField_class _RDF
_RDF = RealDoubleField_class()

RDF = _RDF   # external interface

def RealDoubleField():
    """
    Return the unique instance of the Real Double Field.

    EXAMPLES::

        sage: RealDoubleField() is RealDoubleField()
        True
    """
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

cdef hook_fast_tp_functions():
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
    print "Pool contains %s / %s items" % (element_pool_count, element_pool_size)
