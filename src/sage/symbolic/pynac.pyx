cdef extern from "math.h":
    double sin(double)
    double cos(double)

include "../ext/cdefs.pxi"
include "../ext/stdsage.pxi"

from sage.rings.integer cimport Integer

#################################################################
# Binomial Coefficients
#################################################################

cdef extern from "gmp.h":
    void mpz_bin_uiui(mpz_t, unsigned int, unsigned int)

cdef public object py_binomial_int(int n, unsigned int k):
    cdef bint sign
    if n < 0:
        n = -n + (k-1)
        sign = k%2
    else:
        sign = 0
    cdef Integer ans = PY_NEW(Integer)
    # Compute the binomial coefficient using GMP.
    mpz_bin_uiui(ans.value, n, k)
    # Return the answer or the negative of it (only if k is odd and n is negative).
    if sign:
        return -ans
    else:
        return ans

cdef public object py_binomial(object n, object k):
    # Keep track of the sign we should use.
    cdef bint sign
    if n < 0:
        n = k-n-1
        sign = k%2
    else:
        sign = 0
    # Convert n and k to unsigned ints.
    cdef unsigned int n_ = n, k_ = k
    cdef Integer ans = PY_NEW(Integer)
    # Compute the binomial coefficient using GMP.
    mpz_bin_uiui(ans.value, n_, k_)
    # Return the answer or the negative of it (only if k is odd and n is negative).
    if sign:
        return -ans
    else:
        return ans

def test_binomial(n, k):
    """
    The Binomial coefficients.  It computes the binomial coefficients.  For
    integer n and k and positive n this is the number of ways of choosing k
    objects from n distinct objects.  If n is negative, the formula
    binomial(n,k) == (-1)^k*binomial(k-n-1,k) is used to compute the result.

    INPUT:
        n, k -- integers, with k >= 0.

    OUTPUT:
        integer

    EXAMPLES:
        sage: import sage.symbolic.pynac
        sage: sage.symbolic.pynac.test_binomial(5,2)
        10
        sage: sage.symbolic.pynac.test_binomial(-5,3)
        -35
        sage: -sage.symbolic.pynac.test_binomial(3-(-5)-1, 3)
        -35
    """
    return py_binomial(n, k)

#################################################################
# GCD
#################################################################
import sage.rings.arith
cdef public object py_gcd(object n, object k):
    try:
        return sage.rings.arith.gcd(n,k)
    except (TypeError, ValueError, AttributeError):
        # some strange meaning in case of weird things with no usual lcm.
        return 1


#################################################################
# LCM
#################################################################
cdef public object py_lcm(object n, object k):
    try:
        return sage.rings.arith.lcm(n,k)
    except (TypeError, ValueError, AttributeError):
        # some strange meaning in case of weird things with no usual lcm, e.g.,
        # elements of finite fields.
        return 1


#################################################################
# Real Part
#################################################################
cdef public object py_real(object x):
    try:
        return x.real()
    except AttributeError:
        if isinstance(x, complex):
            return x.real
        return x  # assume is real

#################################################################
# Imaginary Part
#################################################################
cdef public object py_imag(object x):
    try:
        return x.imag()
    except AttributeError:
        if isinstance(x, complex):
            return x.imag
        return 0 # assume is real since it isn't of type complex and doesn't have an imag attribute.


from sage.rings.rational cimport Rational
cdef public bint py_is_rational(object x):
    return PY_TYPE_CHECK_EXACT(x, Rational) or  PY_TYPE_CHECK_EXACT(x, Integer) or\
           IS_INSTANCE(x, int) or IS_INSTANCE(x, long)


cdef public bint py_is_integer(object x):
    return PY_TYPE_CHECK_EXACT(x, Integer) or\
           IS_INSTANCE(x, int) or IS_INSTANCE(x, long)


#cdef public object py_Rational(object x):
#    return Rational(x)


#cdef public object py_Integer(object x):
#    return Integer(x)

cdef public bint py_is_real(object a):
    return py_imag(a) == 0



import sage.rings.arith
cdef public bint py_is_prime(object n):
    try:
        return n.is_prime()
    except:  # yes, I'm doing this on purpose.
        pass
    try:
        return sage.rings.arith.is_prime(n)
    except:
        pass
    return False

cdef public object py_numer(object n):
    if isinstance(n, (int, long, Integer)):
        return n
    try:
        return n.numerator()
    except AttributeError:
        return 1

cdef public object py_denom(object n):
    if isinstance(n, (int, long, Integer)):
        return 1
    try:
        return n.denominator()
    except AttributeError:
        return 1

#################################################################
# SPECIAL FUNCTIONS
#################################################################
from sage.rings.arith import factorial
cdef public object py_factorial(object x):
    return factorial(x)

cdef public object py_doublefactorial(object x):
    raise NotImplementedError #TODO!

from sage.libs.pari.all import pari
cdef public object py_fibonacci(object n):
    return Integer(pari(n).fibonacci())

from sage.rings.arith import bernoulli
cdef public object py_bernoulli(object x):
    return bernoulli(x)

cdef public object py_sin(object x):
    try:
        return x.sin()
    except AttributeError:
        return sin(float(x))

cdef public object py_cos(object x):
    try:
        return x.cos()
    except AttributeError:
        return cos(float(x))

from sage.rings.real_mpfr import RR, RealField

cdef public object py_zeta(object x):
    try:
        return x.zeta()
    except AttributeError:
        return RR(x).zeta()

cdef public object py_exp(object x):
    try:
        return x.exp()
    except AttributeError:
        return RR(x).exp()

cdef public object py_log(object x):
    try:
        return x.log()
    except AttributeError:
        return RR(x).log()


cdef public object py_tan(object x):
    try:
        return x.tan()
    except AttributeError:
        return RR(x).tan()

cdef public object py_asin(object x):
    try:
        return x.arcsin()
    except AttributeError:
        return RR(x).arcsin()

cdef public object py_acos(object x):
    try:
        return x.arccos()
    except AttributeError:
        return RR(x).arccos()

cdef public object py_atan(object x):
    try:
        return x.arctan()
    except AttributeError:
        return RR(x).arctan()

cdef public object py_atan2(object x, object y):
    raise NotImplementedError


cdef public object py_sinh(object x):
    try:
        return x.sinh()
    except AttributeError:
        return RR(x).sinh()


cdef public object py_cosh(object x):
    try:
        return x.cosh()
    except AttributeError:
        return RR(x).cosh()

cdef public object py_tanh(object x):
    try:
        return x.tanh()
    except AttributeError:
        return RR(x).tanh()


cdef public object py_asinh(object x):
    try:
        return x.arcsinh()
    except AttributeError:
        return RR(x).arcsinh()

cdef public object py_acosh(object x):
    try:
        return x.arccosh()
    except AttributeError:
        return RR(x).arccosh()


cdef public object py_atanh(object x):
    try:
        return x.arctanh()
    except AttributeError:
        return RR(x).arctanh()


cdef public object py_li2(object x):
    raise NotImplementedError

cdef public object py_lgamma(object x):
    try:
        return x.lngamma()
    except AttributeError:
        return RR(x).lngamma()

cdef public object py_tgamma(object x):
    raise NotImplementedError

cdef public object py_psi(object x):
    raise NotImplementedError

cdef public object py_psi2(object x, object y):
    raise NotImplementedError

cdef public object py_isqrt(object x):
    raise NotImplementedError

cdef public object py_sqrt(object x):
    try:
        # TODO: worry about recursive call into symbolic system!
        # What if Integer's sqrt calls symbolic one and we go in circle?
        return x.sqrt()
    except AttributeError:
        return RR(x).sqrt()

cdef public object py_abs(object x):
    return abs(x)

cdef public object py_mod(object x, object n):
    return x % n

cdef public object py_smod(object x, object n):
    raise NotImplementedError

cdef public object py_irem(object x, object n):
    raise NotImplementedError

cdef public object py_irem2(object x, object n):
    raise NotImplementedError

cdef public object py_iquo(object x, object n):
    raise NotImplementedError

cdef public object py_iquo2(object x, object n):
    raise NotImplementedError



##################################################################
# Constants
##################################################################

cdef int prec(long ndigits):
    return <int>((ndigits+1) * 3.32192 + 1)

cdef public object py_eval_pi(long ndigits):
    return RealField(prec(ndigits)).pi()

cdef public object py_eval_euler_gamma(long ndigits):
    return RealField(prec(ndigits)).euler_constant()

cdef public object py_eval_catalan(long ndigits):
    return RealField(prec(ndigits)).catalan_constant()



##################################################################
# Constructors
##################################################################
cdef public object py_rational_long(long numer, long denom):
    return Rational((numer, denom))
