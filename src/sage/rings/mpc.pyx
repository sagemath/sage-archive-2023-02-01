"""nodoctest
SAGE Interface to MPC.

AUTHORS: William Stein

DATE: 2005-09-26
"""

#*****************************************************************************
#
#   SAGE: System for Algebra and Geometry Experimentation
#
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
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

include 'interrupt.pxi'

cimport ring
import ring
cimport sage.structure.element
import  sage.structure.element

import sage.rings.coerce
import operator
import sage.rings.rational

#*****************************************************************************
# Headers.  When you past things in here from the mpc docs, be sure to
# remove const's, since those aren't allowed in pyrex.  Also, it can
# be challenging figuring out how to modify things from mpc.h to be
# valid pyrex code.  Note that what is here is only used for
# generating the C code.  The C compiler doesn't see any of this -- it
# only sees mpc.h and stdlib.h, and what Pyrex generates.
#*****************************************************************************

cdef extern from "stdlib.h":
    ctypedef int size_t
    void free(void *ptr)

cdef extern from "gmp.h":
    pass

cdef extern from "mpfr.h":
    int MPFR_PREC_MIN, MPFR_PREC_MAX

cdef extern from "mpc.h":
    ctypedef struct __mpc_struct:
        pass
    ctypedef __mpc_struct* mpc_t
    ctypedef enum mpc_rnd_t:
        GMP_RNDN = 0
        GMP_RNDZ = 1
        GMP_RNDU = 2
        GMP_RNDD = 3
        GMP_RND_MAX = 4
        GMP_RNDNA = -1

    ctypedef mpc_rnd_t mp_rnd_t
    ctypedef long int mp_exp_t
    ctypedef long mp_prec_t

    #mp_rnd_t GMP_RNDZ, GMP_RNDN, GMP_RNDU, GMP_RNDD

    void mpc_init (mpc_t x)
    void mpc_init2 (mpc_t x, mp_prec_t prec)
    void mpc_clear (mpc_t x)

    int mpc_set (mpc_t rop, mpc_t op, mp_rnd_t rnd)
    int mpc_set_str (mpc_t rop, char *s, int base, mp_rnd_t rnd)


    # Arithmetic
    int mpc_add (mpc_t rop, mpc_t op1, mpc_t op2, mp_rnd_t rnd)
    int mpc_sub (mpc_t rop, mpc_t op1, mpc_t op2, mp_rnd_t rnd)
    int mpc_mul (mpc_t rop, mpc_t op1, mpc_t op2, mp_rnd_t rnd)
    int mpc_div (mpc_t rop, mpc_t op1, mpc_t op2, mp_rnd_t rnd)

    # constants
    int mpc_const_log2 (mpc_t rop, mp_rnd_t rnd)
    int mpc_const_pi (mpc_t rop, mp_rnd_t rnd)
    int mpc_const_euler (mpc_t rop, mp_rnd_t rnd)


    # Status functions
    bint mpc_nan_p (mpc_t op)
    bint mpc_inf_p (mpc_t op)
    bint mpc_number_p (mpc_t op)
    bint mpc_zero_p (mpc_t op)

    # Operators
    int mpc_cmp (mpc_t op1, mpc_t op2)


cdef class ComplexNumber(sage.structure.element.RingElement)

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

def mpc_prec_min():
    return MPC_PREC_MIN

def mpc_prec_max():
    return MPC_PREC_MAX

#*****************************************************************************
#
#       Real Field
#
#*****************************************************************************
# The real field is in Pyrex, so mpc elements will have access to
# their parent via direct C calls, which will be faster.

_rounding_modes = ['RNDN', 'RNDZ', 'RNDU', 'RNDD']

cdef class ComplexField(ring.Field):
    cdef int prec
    cdef bint sci_not
    cdef mp_rnd_t rnd

    def __init__(self, int prec=53, rnd="RNDN"):
        """
        ComplexField(prec, rnd):

        INPUT:
            prec -- (integer) precision; default = 53
                    prec is the number of bits used to represent the
                    mantissa of a floating-point number.  The
                    precision can be any integer between mpc_prec_min()
                    and mpc_prec_max(). In the current implementation,
                    mpc_prec_min() is equal to 2.

            rnd -- (string) the rounding mode
                    RNDN -- round to nearest
                    RNDZ -- round towards zero
                    RNDU -- round towards plus infinity
                    RNDD -- round towards minus infinity

        EXAMPLES:
            sage: ComplexField(10)
            Real Field with 10 bits of precision
            sage: ComplexField()
            Real Field with 53 bits of precision
            sage: ComplexField(100000)
            Real Field with 100000 bits of precision

        NOTE: The default precision is 53, since according to the GMP
           manual: "mpc should be able to exactly reproduce all
           computations with double-precision machine floating-point
           numbers (double type in C), except the default exponent
           range is much wider and subnormal numbers are not
           implemented."
        """
        if prec < MPFR_PREC_MIN or prec > MPFR_PREC_MAX:
            raise ValueError, "prec (=%s) must be >= %s and <= %s."%(
                prec, MPFR_PREC_MIN, MPFR_PREC_MAX)
        self.prec = prec
        if not isinstance(rnd, str):
            raise TypeError, "rnd must be a string"
        self.sci_not = 0
        n = _rounding_modes.index(rnd)
        if n == -1:
            raise ValueError, "rnd (=%s) must be one of RNDN, RNDZ, RNDU, or RNDD"%rnd
        self.rnd = n

    def __repr__(self):
        s = "Complex Field with %s bits of precision"%self.prec
        if self.rnd != GMP_RNDN:
            s = s + " and rounding %s"%(self.rnd)
        return s

    def __call__(self, x):
        """
        Coerce x into this real field.

        EXAMPLES:
            sage: C = ComplexField(10)
            sage: R('1.234')
            1.2324
        """
        return ComplexNumber(self, x)

    def __contains__(self, x):
        """
        True if x is in this field.

        EXAMPLES:
            sage: C = ComplexField(10)
            sage: 1 in R
            False
            sage: R(1.0) in R
            True
        """
        return isinstance(x, ComplexNumber) and x.parent() == self

    def __cmp__(self, other):
        """
        EXAMPLES:
            sage: ComplexField(10) == ComplexField(11)
            False
            sage: ComplexField(10) == ComplexField(10)
            True
            sage: ComplexField(10,'RNDN') == ComplexField(10,'RNDZ')
            False
            sage: ComplexField(10) == IntegerRing()
            False
        """
        if not isinstance(other, ComplexField):
            return -1
        cdef ComplexField _other
        _other = other  # to access C structure
        if self.prec == _other.prec and self.rnd == _other.rnd:
            return 0
        return 1

    def __reduce__(self):
        """
        Needed for object persistence.
        """
        return (sage.rings.real_field.__init, (self.prec, ))

    def is_atomic_repr(self):
        """
        Returns True, to signify that elements of this field print
        without sums, so parenthesis aren't required, e.g., in
        coefficients of polynomials.

        EXAMPLES:
            sage: ComplexField(10).is_atomic_repr()
            True
        """
        return True

    def is_finite(self):
        """
        Returns False, since the field of real numbers is not finite.

        EXAMPLES:
            sage: ComplexField(10).is_finite()
            False
        """
        return False

    def characteristic(self):
        """
        Returns 0, since the field of real numbers has characteristic 0.

        EXAMPLES:
            sage: ComplexField(10).characteristic()
            0
        """
        return 0

    def __hash__(self):
        return hash("ComplexField%s_%s"%(self.prec,self.rnd))

    def rounding_mode(self):
        return _rounding_modes[self.rnd]

    def scientific_notation(self, status=None):
        """
        Set or return the scientific notation printing flag.  If this flag
        is true then real numbers with this space as parent print using
        scientific notation.

        INPUT:
            status -- optional flag
        """
        if status is None:
            return self.sci_not
        else:
            self.sci_not = status


C = ComplexField()

#*****************************************************************************
#
#     ComplexNumber -- element of Real Field
#
#
#
#*****************************************************************************
cdef class ComplexNumber(sage.structure.element.RingElement):
    """
    A real number.
    """
    cdef char init
    cdef mpc_t value
    cdef ComplexField _parent

    def __new__(self, ComplexField parent, x=0):
        self.init = 0

    def __init__(self, ComplexField parent, x=0):
        """
        Create a real number.  Should be called by first creating
        a ComplexField, as illustrated in the examples.

        EXAMPLES:
            sage: C = ComplexField()
            sage: R('1.2456')
            1.2455999999999998
            sage: C = ComplexField(3)
            sage: R('1.2456')
            1.0

        NOTES: A real number is an arbitrary precision mantissa with a
        limited precision exponent.  A real number can have three
        special values: Not-a-Number (NaN) or plus or minus
        Infinity. NaN represents an uninitialized object, the result
        of an invalid operation (like 0 divided by 0), or a value that
        cannot be determined (like +Infinity minus
        +Infinity). Moreover, like in the IEEE 754-1985 standard, zero
        is signed, i.e. there are both +0 and ?0; the behavior is the
        same as in the IEEE 754-1985 standard and it is generalized to
        the other functions supported by MPC.

        """
        if parent is None:
            raise TypeError
        self._parent = parent
        mpc_init2(self.value, parent.prec)
        self.init = 1
        if x is None: return
        cdef ComplexNumber _x, n, d
        if isinstance(x, ComplexNumber):
            _x = x  # so we can get at x.value
            mpc_set(self.value, _x.value, GMP_RNDZ)
        elif isinstance(x, sage.rings.rational.Rational):
            n = parent(x.numerator())
            d = parent(x.denominator())
            mpc_div_(self.value, n.value, d.value, parent.rnd)
        else:
            s = str(x)
            raise NotImplementedError

    cdef ComplexField c_parent(ComplexNumber self):
        """
        Returns the parent of self as a ComplexField.

        TODO: we *really* need to work out a way to make this insanely
        fast. This is called extremely frequently. All this overhead is
        probably already unacceptable.
        """
        return <ComplexField> self._parent

    def  __dealloc__(self):
        if self.init:
            mpc_clear(self.value)

    def __repr__(self):
        return self.str(10)

    def real(self):
        """
        Return the real part of self.

        (Since self is a real number, this simply returns self.)
        """
        return self

    def parent(self):
        """
        EXAMPLES:
            sage: C = ComplexField()
            sage: a = R('1.2456')
            sage: a.parent()
            Real Field with 53 bits of precision
        """
        return self._parent

    def str(self, int base=10):
        raise NotImplementedError

    def copy(self):
        cdef ComplexNumber z
        z = ComplexNumber(self._parent)
        mpc_set(z.value, self.value, self._parent.rnd)
        return z

    ########################
    #   Basic Arithmetic
    ########################

    def __add_(ComplexNumber self, ComplexNumber other):
        cdef ComplexNumber x
        x = ComplexNumber(self._parent, None)
        mpc_add_(x.value, self.value, other.value, self._parent.rnd)
        return x

    def __add__(x, y):
        """
        EXAMPLES:
            sage: C = ComplexField()
            sage: R(-1.5) + R(2.5)
            1.0000000000000000
        """
        if isinstance(x, ComplexNumber) and isinstance(y, ComplexNumber):
            return x.__add_(y)
        return sage.rings.coerce.bin_op(x, y, operator.add)

    def __sub_(ComplexNumber self, ComplexNumber other):
        cdef ComplexNumber x
        x = ComplexNumber(self._parent, None)
        mpc_sub_(x.value, self.value, other.value, self._parent.rnd)
        return x

    def __sub__(x, y):
        """
        EXAMPLES:
            sage: C = ComplexField()
            sage: R(-1.5) - R(2.5)
            -4.0000000000000000
        """
        if isinstance(x, ComplexNumber) and isinstance(y, ComplexNumber):
            return x.__sub_(y)
        return sage.rings.coerce.bin_op(x, y, operator.sub)

    def __mul_(ComplexNumber self, ComplexNumber other):
        cdef ComplexNumber x
        x = ComplexNumber(self._parent, None)
        mpc_mul_(x.value, self.value, other.value, self._parent.rnd)
        return x

    def __mul__(x, y):
        """
        EXAMPLES:
            sage: C = ComplexField()
            sage: R(-1.5) * R(2.5)
            -3.7500000000000000
        """
        if isinstance(x, ComplexNumber) and isinstance(y, ComplexNumber):
            return x.__mul_(y)
        return sage.rings.coerce.bin_op(x, y, operator.mul)

    def __div_(ComplexNumber self, ComplexNumber other):
        if not other:
            raise ZeroDivisionError, "ComplexNumber division by zero"
        cdef ComplexNumber x
        x = ComplexNumber(self._parent, None)
        mpc_div_(x.value, self.value, other.value, self._parent.rnd)
        return x

    def __div__(x, y):
        """
        EXAMPLES:
            sage: C = ComplexField()
            sage: R(-1.5) / R(2.5)
            -0.59999999999999998
        """
        if isinstance(x, ComplexNumber) and isinstance(y, ComplexNumber):
            return x.__div_(y)
        return sage.rings.coerce.bin_op(x, y, operator.div)

    def __neg__(self):
        raise NotImplementedError

    def __pos__(self):
        return self

    def __abs__(self):
        return self.abs()

    cdef ComplexNumber abs(ComplexNumber self):
        raise NotImplementedError

    def prec(self):
        return self._parent.prec

    ###########################################
    # Conversions
    ###########################################

    def __float__(self):
        raise NotImplementedError

    def __int__(self):
        return int(float(self))

    def __long__(self):
        return long(float(self))

    def __complex__(self):
        return complex(float(self))

    ###########################################
    # Comparisons: ==, !=, <, <=, >, >=
    ###########################################

    cdef int cmp(ComplexNumber self, ComplexNumber x):
        cdef int i
        i = mpc_cmp(self.value, x.value)
        if i < 0:
            return -1
        elif i == 0:
            return 0
        else:
            return 1

    def __cmp__(ComplexNumber self, ComplexNumber x):
        return self.cmp(x)

    def __richcmp__(ComplexNumber self, x, int op):
        cdef int n
        if not isinstance(x, ComplexNumber):
            try:
                x = ComplexNumber(self.__parent, x)
            except TypeError:
                n = sage.rings.coerce.cmp(self, x)
            else:
                n = self.cmp(x)
        else:
            n = self.cmp(x)
        if op == 0:
            return (n < 0)
        elif op == 1:
            return (n <= 0)
        elif op == 2:
            return (n == 0)
        elif op == 3:
            return (n != 0)
        elif op == 4:
            return (n > 0)
        elif op == 5:
            return (n >= 0)


    ############################
    # Special Functions
    ############################





