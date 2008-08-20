include "../ext/cdefs.pxi"
include "../ext/stdsage.pxi"

from sage.rings.integer cimport Integer

#################################################################
# Binomial Coefficients
#################################################################

cdef extern from "gmp.h":
    void mpz_bin_uiui(mpz_t, unsigned int, unsigned int)

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
    mpz_init(ans.value)
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
        return x  # assume is imag


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

#################################################################
# Factorial
#################################################################
from sage.rings.arith import factorial
cdef public object py_factorial(object x):
    return factorial(x)

from sage.rings.arith import bernoulli
cdef public object py_bernoulli(object x):
    return bernoulli(x)

cdef public object py_sin(object x):
    try:
        return x.sin()
    except AttributeError:
        return float(x).sin()

cdef public object py_cos(object x):
    try:
        return x.cos()
    except AttributeError:
        return float(x).cos()

from sage.rings.real_mpfr import RR

cdef public object py_zeta(object x):
    try:
        return x.zeta()
    except AttributeError:
        return RR(x).zeta()
