"""
Real Numbers

AUTHORS: Kyle Schalm <kschalm@math.utexas.edu> (2005-09)
         William Stein <wstein@ucsd.edu>
"""

#*****************************************************************************
#
#   SAGE: System for Algebra and Geometry Experimentation
#
#       Copyright (C) 2005 William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

#*****************************************************************************
# general TODOs:
#
# more type conversions and coercion. examples:
# sage: R(1/2)
# TypeError: Unable to convert x (='1/2') to mpfr.
#
# sage: 1 + R(42)
# _49 = 1
#
#*****************************************************************************

include 'interrupt.pxi'

cimport ring
import ring
cimport element
import element
import sage.rings.coerce
import operator
import rational

import sage.rings.complex_number

import integer

import sage.rings.infinity


#*****************************************************************************
# Headers.  When you past things in here from mpfr, be sure
# to remove const's, since those aren't allowed in pyrex.  Also, it can be
# challenging figuring out how to modify things from mpfr.h to be valid pyrex
# code.    Note that what is here is only used for generating the C code.
# The C compiler doesn't see any of this -- it only sees mpfr.h and stdlib.h
#*****************************************************************************

cdef extern from "stdlib.h":
    ctypedef int size_t
    void free(void *ptr)

cdef extern from "mpfr.h":
    ctypedef struct __mpfr_struct:
        pass
    #ctypedef __mpfr_struct mpfr_t[1]
    ctypedef __mpfr_struct* mpfr_t
    ctypedef mpfr_t mpfr_ptr
    ctypedef mpfr_t mpfr_srcptr
    ctypedef enum mpfr_rnd_t:
        GMP_RNDN = 0
        GMP_RNDZ = 1
        GMP_RNDU = 2
        GMP_RNDD = 3
        GMP_RND_MAX = 4
        GMP_RNDNA = -1

    ctypedef mpfr_rnd_t mp_rnd_t
    ctypedef long int mp_exp_t
    ctypedef long mp_prec_t

    int MPFR_PREC_MIN, MPFR_PREC_MAX

    #mp_rnd_t GMP_RNDZ, GMP_RNDN, GMP_RNDU, GMP_RNDD

    void mpfr_init (mpfr_t x)
    void mpfr_init2 (mpfr_t x, mp_prec_t prec)
    void mpfr_clear (mpfr_t x)

    int mpfr_set (mpfr_t rop, mpfr_t op, mp_rnd_t rnd)
    int mpfr_set_str (mpfr_t rop, char *s, int base, mp_rnd_t rnd)
    void mpfr_set_inf (mpfr_t x, int sign)
    void mpfr_set_nan (mpfr_t x)

    char * mpfr_get_str (char *str, mp_exp_t *expptr, int base, size_t n, mpfr_t op, mp_rnd_t rnd)
    size_t mpfr_out_str (int *stream, int base, size_t n, mpfr_t op, mp_rnd_t rnd)

    # Arithmetic
    int mpfr_add (mpfr_t rop, mpfr_t op1, mpfr_t op2, mp_rnd_t rnd)
    int mpfr_sub (mpfr_t rop, mpfr_t op1, mpfr_t op2, mp_rnd_t rnd)
    int mpfr_mul (mpfr_t rop, mpfr_t op1, mpfr_t op2, mp_rnd_t rnd)
    int mpfr_div (mpfr_t rop, mpfr_t op1, mpfr_t op2, mp_rnd_t rnd)

    # constants
    int mpfr_const_log2 (mpfr_t rop, mp_rnd_t rnd)
    int mpfr_const_pi (mpfr_t rop, mp_rnd_t rnd)
    int mpfr_const_euler (mpfr_t rop, mp_rnd_t rnd)
    int mpfr_const_catalan (mpfr_t rop, mp_rnd_t rnd)

    # Special functions
    int mpfr_sqrt (mpfr_t rop, mpfr_t op, mp_rnd_t rnd)
    #int mpfr_sqrt_ui _PROTO ((mpfr_ptr, unsigned long, mp_rnd_t));
    int mpfr_cbrt (mpfr_t rop, mpfr_t op, mp_rnd_t rnd)

    int mpfr_log (mpfr_t rop, mpfr_t op, mp_rnd_t rnd)
    int mpfr_log2 (mpfr_t rop, mpfr_t op, mp_rnd_t rnd)
    int mpfr_log10 (mpfr_t rop, mpfr_t op, mp_rnd_t rnd)

    int mpfr_exp (mpfr_t rop, mpfr_t op, mp_rnd_t rnd)
    int mpfr_exp2 (mpfr_t rop, mpfr_t op, mp_rnd_t rnd)
    int mpfr_exp10 (mpfr_t rop, mpfr_t op, mp_rnd_t rnd)

    int mpfr_cos (mpfr_t rop, mpfr_t op, mp_rnd_t rnd)
    int mpfr_sin (mpfr_t rop, mpfr_t op, mp_rnd_t rnd)
    int mpfr_tan (mpfr_t rop, mpfr_t op, mp_rnd_t rnd)
    int mpfr_sin_cos (mpfr_t rop, mpfr_t op, mpfr_t, mp_rnd_t rnd)

    int mpfr_acos (mpfr_ptr, mpfr_srcptr, mp_rnd_t)
    int mpfr_asin (mpfr_ptr, mpfr_srcptr, mp_rnd_t)
    int mpfr_atan (mpfr_ptr, mpfr_srcptr, mp_rnd_t)

    int mpfr_cosh (mpfr_ptr, mpfr_srcptr, mp_rnd_t)
    int mpfr_sinh (mpfr_ptr, mpfr_srcptr, mp_rnd_t)
    int mpfr_tanh (mpfr_ptr, mpfr_srcptr, mp_rnd_t)

    int mpfr_atanh (mpfr_ptr, mpfr_srcptr, mp_rnd_t)
    int mpfr_acosh (mpfr_ptr, mpfr_srcptr, mp_rnd_t)
    int mpfr_asinh (mpfr_ptr, mpfr_srcptr, mp_rnd_t)

    int mpfr_agm (mpfr_t rop, mpfr_t op1, mpfr_t op2, mp_rnd_t rnd)
    int mpfr_gamma (mpfr_t rop, mpfr_t op, mp_rnd_t rnd)
    int mpfr_zeta (mpfr_t rop, mpfr_t op, mp_rnd_t rnd)
    int mpfr_erf (mpfr_t rop, mpfr_t op, mp_rnd_t rnd)

    int mpfr_pow (mpfr_t rop, mpfr_t op1, mpfr_t op2, mp_rnd_t rnd)
    #int mpfr_ui_pow _PROTO ((mpfr_ptr, unsigned long int, mpfr_srcptr, mp_rnd_t));
    #int mpfr_pow_si _PROTO ((mpfr_ptr, mpfr_srcptr, long int, mp_rnd_t));

    #int mpfr_min _PROTO ((mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mp_rnd_t));
    #int mpfr_max _PROTO ((mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mp_rnd_t));

    int mpfr_fac_ui (mpfr_t rop, unsigned long int op, mp_rnd_t rnd)

    int mpfr_abs(mpfr_ptr rop, mpfr_srcptr op, mp_rnd_t rnd)
    int mpfr_sgn(mpfr_t op)
    #define mpfr_abs(a,b,r) mpfr_set4(a,b,r,1)
    #define mpfr_sgn(x) mpfr_cmp_ui(x,0)

    int mpfr_round (mpfr_ptr rop, mpfr_srcptr op)
    int mpfr_trunc (mpfr_ptr rop, mpfr_srcptr op)
    int mpfr_ceil (mpfr_ptr rop, mpfr_srcptr op)
    int mpfr_floor (mpfr_ptr rop, mpfr_srcptr op)
    int mpfr_frac (mpfr_t rop, mpfr_t op, mp_rnd_t rnd)

    # Status functions
    int mpfr_nan_p (mpfr_t op)
    int mpfr_inf_p (mpfr_t op)
    int mpfr_number_p (mpfr_t op)
    int mpfr_zero_p (mpfr_t op)

    double mpfr_get_d (mpfr_t op, mp_rnd_t rnd)


    # Operators

    int mpfr_neg (mpfr_ptr rop, mpfr_srcptr op, mp_rnd_t rnd)
    # int mpfr_eq (mpfr_srcptr rop, mpfr_srcptr op, unsigned long i)
    # int mpfr_less_p (mpfr_t op1, mpfr_t op2)
    int mpfr_cmp (mpfr_t op1, mpfr_t op2)


cdef class RealNumber(element.RingElement)

#*****************************************************************************
#
#       Implementation
#
#*****************************************************************************

#*****************************************************************************
#
#       External Python access to constants
#
#*****************************************************************************

def mpfr_prec_min():
    """
    Return the mpfr variable MPFR_PREC_MIN.
    """
    return MPFR_PREC_MIN

def mpfr_prec_max():
    return MPFR_PREC_MAX

#*****************************************************************************
#
#       Real Field
#
#*****************************************************************************
# The real field is in Pyrex, so mpfr elements will have access to
# their parent via direct C calls, which will be faster.

_rounding_modes = ['RNDN', 'RNDZ', 'RNDU', 'RNDD']

cdef class RealField(ring.Field):
    """
    RealField(prec, sci_not, rnd):

    INPUT:
        prec -- (integer) precision; default = 53
                prec is the number of bits used to represent the
                mantissa of a floating-point number.  The
                precision can be any integer between mpfr_prec_min()
                and mpfr_prec_max(). In the current implementation,
                mpfr_prec_min() is equal to 2.

        sci_not -- (default: True) whether or not to display
                using scientific notation

        rnd -- (string) the rounding mode
                RNDN -- round to nearest
                RNDZ -- round towards zero
                RNDU -- round towards plus infinity
                RNDD -- round towards minus infinity

    EXAMPLES:
        sage: RealField(10)
        Real Field with 10 bits of precision
        sage: RealField()
        Real Field with 53 bits of precision
        sage: RealField(100000)
        Real Field with 100000 bits of precision

    NOTE: The default precision is 53, since according to the GMP
       manual: "mpfr should be able to exactly reproduce all
       computations with double-precision machine floating-point
       numbers (double type in C), except the default exponent
       range is much wider and subnormal numbers are not
       implemented."
    """
    cdef int prec, sci_not
    cdef mp_rnd_t rnd
    cdef object rnd_str

    def __init__(self, int prec=53, int sci_not=0, rnd="RNDN"):
        if prec < MPFR_PREC_MIN or prec > MPFR_PREC_MAX:
            raise ValueError, "prec (=%s) must be >= %s and <= %s."%(
                prec, MPFR_PREC_MIN, MPFR_PREC_MAX)
        self.prec = prec
        if not isinstance(rnd, str):
            raise TypeError, "rnd must be a string"
        self.sci_not = sci_not
        n = _rounding_modes.index(rnd)
        if n == -1:
            raise ValueError, "rnd (=%s) must be one of RNDN, RNDZ, RNDU, or RNDD"%rnd
        self.rnd = n
        self.rnd_str = rnd

    def _repr_(self):
        s = "Real Field with %s bits of precision"%self.prec
        if self.rnd != GMP_RNDN:
            s = s + " and rounding %s"%(self.rnd_str)
        return s

    def _latex_(self):
        return "\\R"

    def __call__(self, x, base=10):
        """
        Coerce x into this real field.

        EXAMPLES:
            sage: R = RealField(10)
            sage: R('1.234')
            1.2344
            sage: R('2', base=2)
            Traceback (most recent call last):
            ...
            TypeError: Unable to convert x (='2') to real number.
            sage: a = R('1.1001', base=2); a
            1.5625
            sage: a.str(2)
            '1.100100000'
        """
        if hasattr(x, '_mpfr_'):
            return x._mpfr_(self)
        return RealNumber(self, x, base)

    def _coerce_(self, x):
        cdef RealField K
        if isinstance(x, RealNumber):
            K = x.parent()
            if K is self:
                return x
            elif K.prec >= self.prec:
                return self(x)
            else:
                raise TypeError
        if isinstance(x, (int, long, integer.Integer, rational.Rational)):
            return self(x)
        raise TypeError

    def __cmp__(self, other):
        """
        EXAMPLES:
            sage: RealField(10) == RealField(11)
            False
            sage: RealField(10) == RealField(10)
            True
            sage: RealField(10,rnd='RNDN') == RealField(10,rnd='RNDZ')
            False
            sage: RealField(10,sci_not=True) == RealField(10,sci_not=False)
            False
            sage: RealField(10) == IntegerRing()
            False
        """
        if not isinstance(other, RealField):
            return -1
        cdef RealField _other
        _other = other  # to access C structure
        if self.prec == _other.prec and self.rnd == _other.rnd \
               and self.sci_not == _other.sci_not:
            return 0
        return 1

    def __reduce__(self):
        """
        EXAMPLES:
            sage: R = RealField(sci_not=1, prec=200, rnd='RNDU')
            sage: loads(dumps(R)) == R
            True
        """
        return sage.rings.real_field.__reduce__RealField, \
                (self.prec, self.sci_not, self.rnd_str)

    def gen(self, i=0):
        if i == 0:
            return self(1)
        else:
            raise IndexError

    def ngens(self):
        return 1

    def gens(self):
        return [self.gen()]

    def _is_valid_homomorphism_(self, codomain, im_gens):
        try:
            s = codomain._coerce_(self(1))
        except TypeError:
            return False
        return s == im_gens[0]

    def is_atomic_repr(self):
        """
        Returns True, to signify that elements of this field print
        without sums, so parenthesis aren't required, e.g., in
        coefficients of polynomials.

        EXAMPLES:
            sage: RealField(10).is_atomic_repr()
            True
        """
        return True

    def is_finite(self):
        """
        Returns False, since the field of real numbers is not finite.

        EXAMPLES:
            sage: RealField(10).is_finite()
            False
        """
        return False

    def characteristic(self):
        """
        Returns 0, since the field of real numbers has characteristic 0.

        EXAMPLES:
            sage: RealField(10).characteristic()
            0
        """
        return 0

    def name(self):
        return "RealField%s_%s"%(self.prec,self.rnd)

    def __hash__(self):
        return hash(self.name())

    def precision(self):
        return self.prec

    # int mpfr_const_pi (mpfr_t rop, mp_rnd_t rnd)
    def pi(self):
        """
        Returns pi to the precision of this field.

        EXAMPLES:
            sage: R = RealField(100)
            sage: R.pi()
            3.1415926535897932384626433832793
            sage: R.pi().sqrt()/2
            0.88622692545275801364908374167063
        """
        cdef RealNumber x
        x = RealNumber(self)
        mpfr_const_pi(x.value, self.rnd)
        return x


    # int mpfr_const_euler (mpfr_t rop, mp_rnd_t rnd)
    def euler_constant(self):
        """
        Returns Euler's gamma constant to the precision of this field.

        EXAMPLES:
            sage: RealField(100).euler_constant()
            0.57721566490153286060651209008234
        """
        cdef RealNumber x
        x = RealNumber(self)
        mpfr_const_euler(x.value, self.rnd)
        return x

    # int mpfr_const_catalan (mpfr_t rop, mp_rnd_t rnd)
    def catalan_constant(self):
        """
        Returns Catalan's constant to the precision of this field.

        EXAMPLES:
            sage: RealField(100).catalan_constant()
            0.91596559417721901505460351493252
        """
        cdef RealNumber x
        x = RealNumber(self)
        mpfr_const_catalan(x.value, self.rnd)
        return x

    # int mpfr_const_log2 (mpfr_t rop, mp_rnd_t rnd)
    def log2(self):
        """
        Returns log(2) to the precision of this field.

        EXAMPLES:
            sage: R=RealField(100)
            sage: R.log2()
            0.69314718055994530941723212145798
            sage: R(2).log()
            0.69314718055994530941723212145798
        """
        cdef RealNumber x
        x = RealNumber(self)
        mpfr_const_log2(x.value, self.rnd)
        return x

    def factorial(self, int n):
        """
        Return the factorial of the integer n as a real number.
        """
        cdef RealNumber x
        if n < 0:
            raise ArithmeticError, "n must be nonnegative"
        x = RealNumber(self)
        mpfr_fac_ui(x.value, n, self.rnd)
        return x

    def rounding_mode(self):
        return _rounding_modes[self.rnd]

    def scientific_notation(self, status=None):
        """
        Set or return the scientific notation printing flag.  If this flag
        is True then real numbers with this space as parent print using
        scientific notation.

        INPUT:
            status -- (bool --) optional flag
        """
        if status is None:
            return bool(self.sci_not)
        else:
            self.sci_not = status

    def zeta(self, n=2):
        """
        Return an $n$-th root of unity in the real field,
        if one exists, or raise a ValueError otherwise.

        EXAMPLES:
            sage: R = RealField()
            sage: R.zeta()
            -1.0000000000000000
            sage: R.zeta(1)
            1.0000000000000000
            sage: R.zeta(5)
            Traceback (most recent call last):
            ...
            ValueError: No 5th root of unity in self
        """
        if n == 1:
            return self(1)
        elif n == 2:
            return self(-1)
        raise ValueError, "No %sth root of unity in self"%n

R = RealField()

#*****************************************************************************
#
#     RealNumber -- element of Real Field
#
#
#
#*****************************************************************************
cdef class RealNumber(element.RingElement):
    """
    A real number.
    """
    cdef char init
    cdef mpfr_t value
    cdef RealField _parent

    def __new__(self, RealField parent, x=0, int base=10):
        self.init = 0

    def __init__(self, RealField parent, x=0, int base=10, special=None):
        """
        Create a real number.  Should be called by first creating
        a RealField, as illustrated in the examples.

        EXAMPLES:
            sage: R = RealField()
            sage: R('1.2456')
            1.2456000000000000
            sage: R = RealField(3)
            sage: R('1.2456')
            1.2

        EXAMPLE: Rounding Modes
            sage: w = RealField(3)(5/2)
            sage: RealField(2, rnd="RNDZ")(w).str(2)
            '10'
            sage: RealField(2, rnd="RNDD")(w).str(2)
            '10'
            sage: RealField(2, rnd="RNDU")(w).str(2)
            '11'
            sage: RealField(2, rnd="RNDN")(w).str(2)
            '10'

        NOTES: A real number is an arbitrary precision mantissa with a
        limited precision exponent.  A real number can have three
        special values: Not-a-Number (NaN) or plus or minus
        Infinity. NaN represents an uninitialized object, the result
        of an invalid operation (like 0 divided by 0), or a value that
        cannot be determined (like +Infinity minus
        +Infinity). Moreover, like in the IEEE 754-1985 standard, zero
        is signed, i.e. there are both +0 and ?0; the behavior is the
        same as in the IEEE 754-1985 standard and it is generalized to
        the other functions supported by MPFR.

        """
        if parent is None:
            raise TypeError
        self._parent = parent
        mpfr_init2(self.value, parent.prec)
        self.init = 1
        if x is None: return
        cdef RealNumber _x, n, d
        if isinstance(x, RealNumber):
            _x = x  # so we can get at x.value
            mpfr_set(self.value, _x.value, parent.rnd)
        elif isinstance(x, rational.Rational):
            n = parent(x.numerator())
            d = parent(x.denominator())
            mpfr_div(self.value, n.value, d.value, parent.rnd)
        else:
            s = str(x)
            if mpfr_set_str(self.value, s, base, parent.rnd):
                if s == 'NaN' or s == '@NaN@':
                    mpfr_set_nan(self.value)
                elif s == '+infinity':
                    mpfr_set_inf(self.value, 1)
                elif s == '-infinity':
                    mpfr_set_inf(self.value, -1)
                else:
                    raise TypeError, "Unable to convert x (='%s') to real number."%s

    def __reduce__(self):
        """
        EXAMPLES:
            sage: R = RealField(sci_not=1, prec=200, rnd='RNDU')
            sage: b = R('393.39203845902384098234098230948209384028340')
            sage: loads(dumps(b)) == b
            True
            sage: b = R(1)/R(0); b
            +infinity
            sage: loads(dumps(b)) == b
            True
            sage: b = R(-1)/R(0); b
            -infinity
            sage: loads(dumps(b)) == b
            True
            sage: b = R(-1).sqrt(); b
            NaN
            sage: loads(dumps(b)) == b
            True
        """
        s = self.str(32, no_sci=False, e='@')
        return (sage.rings.real_field.__reduce__RealNumber, (self._parent, s, 32))

    def  __dealloc__(self):
        if self.init:
            mpfr_clear(self.value)

    def __repr__(self):
        return self.str(10)

    def _latex_(self):
        return str(self)

    def _interface_init_(self):
        """
        Return string representation of self in base 10 with
        no scientific notation.

        This is most likely to make sense in other computer algebra
        systems (this function is the default for exporting to other
        computer algebra systems).

        EXAMPLES:
            sage: n = 1.3939494594
            sage: n._interface_init_()
            '1.3939494593999999'
        """
        return self.str(10, no_sci=True)

    def __hash__(self):
        return hash(self.str(16))

    def _im_gens_(self, codomain, im_gens):
        return codomain(self) # since 1 |--> 1

    def real(self):
        """
        Return the real part of self.

        (Since self is a real number, this simply returns self.)
        """
        return self

    def parent(self):
        """
        EXAMPLES:
            sage: R = RealField()
            sage: a = R('1.2456')
            sage: a.parent()
            Real Field with 53 bits of precision
        """
        return self._parent

    def str(self, int base=10, no_sci=None, e='e'):
        if base < 2 or base > 36:
            raise ValueError, "the base (=%s) must be between 2 and 36"%base
        if mpfr_nan_p(self.value):
            if base >= 24:
                return "@NaN@"
            else:
                return "NaN"
        elif mpfr_inf_p(self.value):
            if mpfr_sgn(self.value) > 0:
                return "+infinity"
            else:
                return "-infinity"

        cdef char *s
        cdef mp_exp_t exponent

        _sig_on
        s = mpfr_get_str(<char*>0, &exponent, base, 0,
                         self.value, self._parent.rnd)
        _sig_off
        if s == <char*> 0:
            raise RuntimeError, "Unable to convert an mpfr number to a string."
        t = str(s)
        free(s)

        if no_sci==False or (self._parent.sci_not and not (no_sci==True)):
            if t[0] == "-":
                return "-%s.%s%s%s"%(t[1:2], t[2:], e, exponent-1)
            return "%s.%s%s%s"%(t[0], t[1:], e, exponent-1)

        n = abs(exponent)
        lpad = ''
        if exponent <= 0:
            n = len(t)
            lpad = '0.' + '0'*abs(exponent)
        if t[0] == '-':
            lpad = '-' + lpad
            t = t[1:]
        z = lpad + str(t[:n])
        w = t[n:]
        if len(w) > 0:
            z = z + ".%s"%w
        elif lpad == '':
            z = z + '0'*(n-len(t))
        return z

    def copy(self):
        cdef RealNumber z
        z = RealNumber(self._parent)
        mpfr_set(z.value, self.value, self._parent.rnd)
        return z

    def integer_part(self):
        """
        If in decimal this number is written n.defg, returns n.
        """
        s = self.str(base=10, no_sci=True)
        i = s.find(".")
        return integer.Integer(s[:i])

    ########################
    #   Basic Arithmetic
    ########################

    def _add_(RealNumber self, RealNumber other):
        """
        Add two real numbers with the same parent.

        EXAMPLES:
            sage: R = RealField()
            sage: R(-1.5) + R(2.5)
            1.0000000000000000
        """
        cdef RealNumber x
        x = RealNumber(self._parent, None)
        mpfr_add(x.value, self.value, other.value, self._parent.rnd)
        return x

    def __invert__(self):
        return self.parent()(1) / self

    def _sub_(RealNumber self, RealNumber other):
        """
        Subtract two real numbers with the same parent.

        EXAMPLES:
            sage: R = RealField()
            sage: R(-1.5) - R(2.5)
            -4.0000000000000000
        """
        cdef RealNumber x
        x = RealNumber(self._parent, None)
        mpfr_sub(x.value, self.value, other.value, self._parent.rnd)
        return x

    def _mul_(RealNumber self, RealNumber other):
        """
        Multiply two real numbers with the same parent.

        EXAMPLES:
            sage: R = RealField()
            sage: R(-1.5) * R(2.5)
            -3.7500000000000000

        If two elements have different precision, arithmetic
        operations are performed after coercing to the lower
        precision.

            sage: R10 = RealField(10)
            sage: R100 = RealField(100)
            sage: a = R10('393.3902834028345')
            sage: b = R100('393.3902834028345')
            sage: a
            393.50
            sage: b
            393.39028340283450000000000000011
            sage: a*b
            154880
            sage: b*a
            154880
            sage: parent(b*a)
            Real Field with 10 bits of precision
        """
        cdef RealNumber x
        x = RealNumber(self._parent, None)
        mpfr_mul(x.value, self.value, other.value, self._parent.rnd)
        return x

    def __div_(RealNumber self, RealNumber other):
        if not other:
            raise ZeroDivisionError, "RealNumber division by zero"
        cdef RealNumber x
        x = RealNumber(self._parent, None)
        mpfr_div(x.value, self.value, other.value, self._parent.rnd)
        return x

    def __div__(x, y):
        """
        EXAMPLES:
            sage: R = RealField()
            sage: R(-1.5) / R(2.5)
            -0.59999999999999998
        """
        if isinstance(x, RealNumber) and isinstance(y, RealNumber):
            return x.__div_(y)
        return sage.rings.coerce.bin_op(x, y, operator.div)

    def __neg__(self):
        cdef RealNumber x
        x = RealNumber(self._parent, None)
        mpfr_neg(x.value, self.value, self._parent.rnd)
        return x

    def __pos__(self):
        return self

    def __abs__(self):
        return self.abs()

    cdef RealNumber abs(RealNumber self):
        cdef RealNumber x
        x = RealNumber(self._parent, None)
        mpfr_abs(x.value, self.value, self._parent.rnd)
        return x

    def multiplicative_order(self):
        if self == 1:
            return 1
        elif self == -1:
            return -1
        return sage.rings.infinity.infinity

    def sign(self):
        return mpfr_sgn(self.value)

    def prec(self):
        return self._parent.prec

    ###################
    # Rounding etc
    ###################

    def round(self):
        cdef RealNumber x
        x = RealNumber(self._parent, None)
        mpfr_round(x.value, self.value)
        return x

    def floor(self):
        cdef RealNumber x
        x = RealNumber(self._parent, None)
        mpfr_floor(x.value, self.value)
        return x

    def ceil(self):
        cdef RealNumber x
        x = RealNumber(self._parent, None)
        mpfr_ceil(x.value, self.value)
        return x

    def trunc(self):
        cdef RealNumber x
        x = RealNumber(self._parent, None)
        mpfr_trunc(x.value, self.value)
        return x

    def frac(self):
        """
        frac returns a real number > -1 and < 1. it satisfies the
        relation:
            x = x.trunc() + x.frac()
        """
        cdef RealNumber x
        x = RealNumber(self._parent, None)
        mpfr_frac(x.value, self.value, self._parent.rnd)
        return x

    ###########################################
    # Conversions
    ###########################################

    def __float__(self):
        return mpfr_get_d(self.value, self._parent.rnd)

    def __int__(self):
        """
        Returns integer truncation of this real number.
        """
        s = self.str(32)
        i = s.find('.')
        return int(s[:i], 32)

    def __long__(self):
        """
        Returns long integer truncation of this real number.
        """
        s = self.str(32)
        i = s.find('.')
        return long(s[:i], 32)

    def __complex__(self):
        return complex(float(self))

    ###########################################
    # Comparisons: ==, !=, <, <=, >, >=
    ###########################################

    cdef int cmp(RealNumber self, RealNumber x):
        cdef int a,b
        a = mpfr_nan_p(self.value)
        b = mpfr_nan_p(x.value)
        if a != b:
            return -1    # nothing equal to Nan
        cdef int i
        i = mpfr_cmp(self.value, x.value)
        if i < 0:
            return -1
        elif i == 0:
            return 0
        else:
            return 1

    def __cmp__(RealNumber self, RealNumber x):
        return self.cmp(x)

    def __richcmp__(RealNumber self, x, int op):
        cdef int n
        if not isinstance(x, RealNumber):
            try:
                x = RealNumber(self._parent, x)
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
        Return the square root of self.
        """
        cdef RealNumber x
        x = RealNumber(self._parent, None)
        _sig_on
        mpfr_sqrt(x.value, self.value, self._parent.rnd)
        _sig_off
        return x

    def cube_root(self):
        """
        Return the cubic root (defined over the real numbers) of self.
        """
        cdef RealNumber x
        x = RealNumber(self._parent, None)
        _sig_on
        mpfr_cbrt(x.value, self.value, self._parent.rnd)
        _sig_off
        return x

    def __pow(self, RealNumber exponent):
        cdef RealNumber x
        x = RealNumber(self._parent, None)
        _sig_on
        mpfr_pow(x.value, self.value, exponent.value, self._parent.rnd)
        _sig_off
        return x

    def __pow__(self, exponent, modulus):
        """
        Compute self raised to the power of exponent, rounded in
        the direction specified by the parent of self.

        EXAMPLES:
            sage: R = RealField(30)
            sage: a = R('1.23456')
            sage: a^20
            67.646297455
            sage: a^a
            1.2971114814
            sage: b = R(-1)
            sage: b^0.5
            NaN
        """
        cdef RealNumber x
        if not isinstance(self, RealNumber):
            return self.__pow__(float(exponent))
        if not isinstance(exponent, RealNumber):
            x = self
            exponent = x._parent(exponent)
        return self.__pow(exponent)

    def log(self):
        """
        EXAMPLES:
            sage: R = RealField()
            sage: R(2).log()
            0.69314718055994529
        """
        cdef RealNumber x
        x = RealNumber(self._parent, None)
        _sig_on
        mpfr_log(x.value, self.value, self._parent.rnd)
        _sig_off
        return x

    def log2(self):
        cdef RealNumber x
        x = RealNumber(self._parent, None)
        _sig_on
        mpfr_log2(x.value, self.value, self._parent.rnd)
        _sig_off
        return x

    def log10(self):
        cdef RealNumber x
        x = RealNumber(self._parent, None)
        _sig_on
        mpfr_log10(x.value, self.value, self._parent.rnd)
        _sig_off
        return x

    def exp(self):
        cdef RealNumber x
        x = RealNumber(self._parent, None)
        _sig_on
        mpfr_exp(x.value, self.value, self._parent.rnd)
        _sig_off
        return x

    def exp2(self):
        cdef RealNumber x
        x = RealNumber(self._parent, None)
        _sig_on
        mpfr_exp2(x.value, self.value, self._parent.rnd)
        _sig_off
        return x

    def exp10(self):
        cdef RealNumber x
        x = RealNumber(self._parent, None)
        _sig_on
        mpfr_exp10(x.value, self.value, self._parent.rnd)
        _sig_off
        return x

    def cos(self):
        cdef RealNumber x
        x = RealNumber(self._parent, None)
        _sig_on
        mpfr_cos(x.value, self.value, self._parent.rnd)
        _sig_off
        return x

    ##########################################################
    # it would be nice to get zero back here:
    # sage: R(-1).acos().sin()
    # _57 = -0.50165576126683320234e-19
    # i think this could be "fixed" by using MPFI. (put on to-do list.)
    #
    # this seems to work ok:
    # sage: R(-1).acos().cos()
    # _58 = -0.10000000000000000000e1
    def sin(self):
        """
        EXAMPLES:
            sage: R = RealField(100)
            sage: R(2).sin()
            0.90929742682568169539601986591150
        """
        cdef RealNumber x
        x = RealNumber(self._parent, None)
        _sig_on
        mpfr_sin(x.value, self.value, self._parent.rnd)
        _sig_off
        return x

    def tan(self):
        cdef RealNumber x
        x = RealNumber(self._parent, None)
        _sig_on
        mpfr_tan(x.value, self.value, self._parent.rnd)
        _sig_off
        return x

    def sincos(self):
        """
        Returns a pair consisting of the sine and cosine.

        EXAMPLES:
            sage: R = RealField()
            sage: t = R.pi()/6
            sage: t.sincos()
            (0.49999999999999994, 0.86602540378443871)
        """
        cdef RealNumber x,y
        x = RealNumber(self._parent, None)
        y = RealNumber(self._parent, None)
        _sig_on
        mpfr_sin_cos(x.value, y.value, self.value, self._parent.rnd)
        _sig_off
        return x,y


    # int mpfr_sin_cos (mpfr_t rop, mpfr_t op, mpfr_t, mp_rnd_t rnd)

    def acos(self):
        cdef RealNumber x
        x = RealNumber(self._parent, None)
        _sig_on
        mpfr_acos(x.value, self.value, self._parent.rnd)
        _sig_off
        return x

    def asin(self):
        cdef RealNumber x
        x = RealNumber(self._parent, None)
        _sig_on
        mpfr_asin(x.value, self.value, self._parent.rnd)
        _sig_off
        return x

    def atan(self):
        cdef RealNumber x
        x = RealNumber(self._parent, None)
        _sig_on
        mpfr_atan(x.value, self.value, self._parent.rnd)
        _sig_off
        return x

    #int mpfr_acos _PROTO ((mpfr_ptr, mpfr_srcptr, mp_rnd_t));
    #int mpfr_asin _PROTO ((mpfr_ptr, mpfr_srcptr, mp_rnd_t));
    #int mpfr_atan _PROTO ((mpfr_ptr, mpfr_srcptr, mp_rnd_t));

    def cosh(self):
        cdef RealNumber x
        x = RealNumber(self._parent, None)
        _sig_on
        mpfr_cosh(x.value, self.value, self._parent.rnd)
        _sig_off
        return x

    def sinh(self):
        cdef RealNumber x
        x = RealNumber(self._parent, None)
        _sig_on
        mpfr_sinh(x.value, self.value, self._parent.rnd)
        _sig_off
        return x

    def tanh(self):
        cdef RealNumber x
        x = RealNumber(self._parent, None)
        _sig_on
        mpfr_tanh(x.value, self.value, self._parent.rnd)
        _sig_off
        return x

    def acosh(self):
        cdef RealNumber x
        x = RealNumber(self._parent, None)
        _sig_on
        mpfr_acosh(x.value, self.value, self._parent.rnd)
        _sig_off
        return x

    def asinh(self):
        cdef RealNumber x
        x = RealNumber(self._parent, None)
        _sig_on
        mpfr_asinh(x.value, self.value, self._parent.rnd)
        _sig_off
        return x

    def atanh(self):
        cdef RealNumber x
        x = RealNumber(self._parent, None)
        _sig_on
        mpfr_atanh(x.value, self.value, self._parent.rnd)
        _sig_off
        return x

    def agm(self, other):
        """
        Return the arithmetic-geometric mean of self and other. The
        arithmetic-geometric mean is the common limit of the sequences
        $u_n$ and $v_n$, where $u_0$ is self, $v_0$ is other,
        $u_{n+1}$ is the arithmetic mean of $u_n$ and $v_n$, and
        $v_{n+1}$ is the geometric mean of u_n and v_n. If any operand
        is negative, the return value is \code{NaN}.
        """
        cdef RealNumber x, _other
        if not isinstance(other, RealNumber) or other.parent() != self._parent:
            _other = self._parent(other)
        else:
            _other = other
        x = RealNumber(self._parent, None)
        _sig_on
        mpfr_agm(x.value, self.value, _other.value, self._parent.rnd)
        _sig_off
        return x


    def erf(self):
        """
        Returns the value of the error function on self.

        EXAMPLES:
           sage: R = RealField()
           sage: R(6).erf()
           1.0000000000000000
        """
        cdef RealNumber x
        x = RealNumber(self._parent, None)
        _sig_on
        mpfr_erf(x.value, self.value, self._parent.rnd)
        _sig_off
        return x


    def gamma(self):
        """
        The Euler gamma function. Return gamma of self.

        EXAMPLES:
           sage: R = RealField()
           sage: R(6).gamma()
           120.00000000000000
           sage: R(1.5).gamma()
           0.88622692545275805
        """
        cdef RealNumber x
        x = RealNumber(self._parent, None)
        _sig_on
        mpfr_gamma(x.value, self.value, self._parent.rnd)
        _sig_off
        return x

    def zeta(self):
        r"""
        Return the Riemann zeta function evaluated at this real number.

        \note{PARI is vastly more efficient at computing the Riemann zeta
        function.   See the example below for how to use it.}

        EXAMPLES:
            sage: R = RealField()
            sage: R(2).zeta()
            1.6449340668482264
            sage: R.pi()^2/6
            1.6449340668482264
            sage: R(-2).zeta()
            0.00000000000000000
            sage: R(1).zeta()
            +infinity

        Computing zeta using PARI is much more efficient in difficult cases.
        Here's how to compute zeta with at least a given precision:

             sage: z = pari.new_with_bits_prec(2, 53).zeta(); z
             1.644934066848226436472415167              # 32-bit
             1.6449340668482264364724151666460251892    # 64-bit

        Note that the number of bits of precision in the constructor only
        effects the internel precision of the pari number, not the number
        of digits that gets displayed.  To increase that you must
        use \code{pari.set_real_precision}.

             sage: type(z)
             <type 'gen.gen'>
             sage: R(z)
             1.6449340668482264
        """
        cdef RealNumber x
        x = RealNumber(self._parent, None)
        _sig_on
        mpfr_zeta(x.value, self.value, self._parent.rnd)
        _sig_off
        return x

