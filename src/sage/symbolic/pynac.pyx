###############################################################################
#   SAGE: Open Source Mathematical Software
#       Copyright (C) 2008 William Stein <wstein@gmail.com>
#       Copyright (C) 2008 Burcin Erocal <burcin@erocal.org>
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or any later version.  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
###############################################################################

cdef extern from "math.h":
    double sin(double)
    double cos(double)
    double sqrt(double)

include "../ext/cdefs.pxi"
include "../ext/stdsage.pxi"
include "../libs/ginac/decl.pxi"

from sage.structure.element import Element
from sage.rings.integer_ring import ZZ
from sage.rings.integer cimport Integer
from sage.rings.real_mpfr import RR, RealField
from sage.rings.all import CC
from sage.calculus.all import SR

from sage.symbolic.expression cimport Expression, new_Expression_from_GEx

import constants
import ring

#################################################################
# Symbolic function helpers
#################################################################

cdef extern from *:
    object exvector_to_PyTuple(GExVector seq)
    GEx pyExpression_to_ex(object res) except *
    object ex_to_pyExpression(GEx juice)

cdef public object ex_to_pyExpression(GEx juice):
    """
    Convert given GiNaC::ex object to a python Expression class.

    Used to pass parameters to custom power and series functions.
    """
    cdef Expression nex
    nex = <Expression>PY_NEW(Expression)
    GEx_construct_ex(&nex._gobj, juice)
    nex._parent = ring.NSR
    return nex

cdef public object exvector_to_PyTuple(GExVector seq):
    """
    Converts arguments list given to a function to a PyTuple.

    Used to pass arguments to python functions assigned to be custom
    evaluation, derivative, etc. functions of symbolic functions.
    """

    return tuple([new_Expression_from_GEx(seq.at(i)) \
            for i from 0 <= i < seq.size()])

cdef public GEx pyExpression_to_ex(object res) except *:
    """
    Converts an Expression object to a GiNaC::ex.

    Used to pass returen values of custom python evaluation, derivation
    functions back to C++ level.
    """
    if res is None:
        raise TypeError, "eval function returned None, expected return value of type sage.symbolic.expression.Expression"
    try:
        t = ring.NSR._coerce_c(res)
    except TypeError, err:
        raise TypeError, "eval function did not return a symbolic expression or an element that can be coerced into a symbolic expression"
    return (<Expression>t)._gobj

cdef int GINAC_FN_SERIAL = 0

cdef set_ginac_fn_serial():
    """
    Initialize the GINAC_FN_SERIAL variable to the number of functions
    defined by GiNaC. This allows us to prevent collisions with C++ level
    special functions when a user asks to construct a symbolic function
    with the same name.
    """
    global GINAC_FN_SERIAL
    GINAC_FN_SERIAL = g_registered_functions().size()

def get_ginac_serial():
    """
    Returns the number of C++ level functions defined by GiNaC.

    EXAMPLES:
        sage: from sage.symbolic.pynac import get_ginac_serial
        sage: get_ginac_serial() >= 40
        True
    """
    return GINAC_FN_SERIAL

#################################################################
# Printing helpers
#################################################################

cdef extern from *:
    char* py_latex(object o) except +

cdef public char* py_latex(object o) except +:
    from sage.misc.latex import latex
    s = latex(o)
    return s

#################################################################
# Modular helpers
#################################################################

from sage.structure.element cimport Element

cdef extern from *:
    int py_get_parent_char(object o) except *

cdef public int py_get_parent_char(object o) except *:
    if isinstance(o, Element):
        return (<Element>o)._parent.characteristic()
    else:
        return 0

#################################################################
# Binomial Coefficients
#################################################################

# We declare the functions defined below as extern here, to prevent Cython
# from generating separate declarations for them which confuse g++
cdef extern from *:
    object py_binomial_int(int n, unsigned int k)
    object py_binomial(object n, object k)
    object py_gcd(object n, object k)
    object py_lcm(object n, object k)
    object py_real(object x)
    object py_imag(object x)
    object py_conjugate(object x)
    bint py_is_rational(object x)
    bint py_is_integer(object x)
    bint py_is_real(object x)
    bint py_is_prime(object x)
    object py_numer(object x)
    object py_denom(object x)
    object py_float(object x)
    object py_RDF_from_double(double x)
    object py_factorial(object x)
    object py_doublefactorial(object x)
    object py_fibonacci(object x)
    object py_step(object x)
    object py_bernoulli(object x)
    object py_sin(object x)
    object py_cos(object x)
    object py_zeta(object x)
    object py_exp(object x)
    object py_log(object x)
    object py_tan(object x)
    object py_asin(object x)
    object py_acos(object x)
    object py_atan(object x)
    object py_atan2(object x, object y)
    object py_sinh(object x)
    object py_cosh(object x)
    object py_tanh(object x)
    object py_asinh(object x)
    object py_acosh(object x)
    object py_atanh(object x)
    object py_tgamma(object x)
    object py_lgamma(object x)
    object py_isqrt(object x)
    object py_sqrt(object x)
    object py_abs(object x)
    object py_mod(object x, object y)
    object py_smod(object x, object y)
    object py_irem(object x, object y)
    object py_iquo(object x, object y)
    object py_iquo2(object x, object y)
    int py_int_length(object x) except -1
    object py_li2(object x)
    object py_psi(object x)
    object py_psi2(object x, object y)
    object py_eval_pi(long int)
    object py_eval_euler_gamma(long int)
    object py_eval_catalan(long int)
    object py_integer_from_long(long int)
    object py_integer_from_python_obj(object x)

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
    if PY_TYPE_CHECK_EXACT(n, Rational) and PY_TYPE_CHECK_EXACT(k, Rational):
        return n.content(k)
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
    """
    Returns the real part of x.

    TESTS::

        sage: from sage.symbolic.pynac import py_real_for_doctests as py_real
        sage: py_real(I)
        0.0
        sage: py_real(CC(1,5))
        1.00000000000000
        sage: py_real(CC(1))
        1.00000000000000
        sage: py_real(RR(1))
        1.00000000000000

        sage: py_real(Mod(2,7))
        2

        sage: py_real(QQ['x'].gen())
        x
    """
    try:
        return x.real()
    except AttributeError:
        pass
    if isinstance(x, complex):
        return x.real
    from sage.misc.functional import real
    try:
        return real(x)
    except TypeError:
        return x # assume x is real

def py_real_for_doctests(x):
    """
    Used for doctesting py_real.

    TESTS::

        sage: from sage.symbolic.pynac import py_real_for_doctests
        sage: py_real_for_doctests(I)
        0.0
    """
    return py_real(x)

#################################################################
# Imaginary Part
#################################################################
cdef public object py_imag(object x):
    """
    Return the imaginary part of x.

    TESTS::

        sage: from sage.symbolic.pynac import py_imag_for_doctests as py_imag
        sage: py_imag(I)
        1.0
        sage: py_imag(CC(1,5))
        5.00000000000000
        sage: py_imag(CC(1))
        0.000000000000000
        sage: py_imag(RR(1))
        0.0

        sage: py_imag(Mod(2,7))
        0.0

        sage: py_imag(QQ['x'].gen())
        0
    """
    try:
        return x.imag()
    except AttributeError:
        pass

    if isinstance(x, complex):
        return x.imag
    from sage.misc.functional import imag
    try:
        return imag(x)
    except TypeError:
        return 0 # assume x is real

def py_imag_for_doctests(x):
    """
    Used for doctesting py_imag.

    TESTS::

        sage: from sage.symbolic.pynac import py_imag_for_doctests
        sage: py_imag_for_doctests(I)
        1.0
    """
    return py_imag(x)


#################################################################
# Conjugate
#################################################################
cdef public object py_conjugate(object x):
    try:
        return x.conjugate()
    except AttributeError:
        return x # assume is real since it doesn't have an imag attribute.

from sage.rings.rational cimport Rational
cdef public bint py_is_rational(object x):
    return PY_TYPE_CHECK_EXACT(x, Rational) or \
           PY_TYPE_CHECK_EXACT(x, Integer) or\
           IS_INSTANCE(x, int) or IS_INSTANCE(x, long)


cdef public bint py_is_integer(object x):
    r"""
    Returns True if pynac should treat this object as an integer.

    EXAMPLES:
        sage: from sage.symbolic.pynac import py_is_integer_for_doctests
        sage: py_is_integer = py_is_integer_for_doctests

        sage: py_is_integer(1r)
        True
        sage: py_is_integer(long(1))
        True
        sage: py_is_integer(3^57)
        True
        sage: py_is_integer(SR(5))
        True
        sage: py_is_integer(4/2)
        True
        sage: py_is_integer(QQbar(sqrt(2))^2)
        True
        sage: py_is_integer(3.0)
        False
        sage: py_is_integer(3.0r)
        False
    """
    return IS_INSTANCE(x, int) or IS_INSTANCE(x, long) or \
           (IS_INSTANCE(x, Element) and
            (x.parent().is_exact() or x.parent() == SR) and
            (x in ZZ))

def py_is_integer_for_doctests(x):
    return py_is_integer(x)

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

cdef public object py_float(object n):
    return float(n)

# TODO: Optimize this
from sage.rings.real_double import RDF
cdef public object py_RDF_from_double(double x):
    return RDF(x)

#################################################################
# SPECIAL FUNCTIONS
#################################################################
from sage.rings.arith import factorial
cdef public object py_factorial(object x):
    return factorial(x)

cdef public object py_doublefactorial(object x):
    n = Integer(x)
    if n < -1:
        raise ValueError, "argument must be >= -1"
    from sage.misc.misc_c import prod  # fast balanced product
    return prod([n - 2*i for i in range(n//2)])

def doublefactorial(n):
    """
    The double factorial combinatorial function:
        n!! == n * (n-2) * (n-4) * ... * ({1|2}) with 0!! == (-1)!! == 1.

    INPUT:
        n -- an integer > = 1

    EXAMPLES:
        sage: from sage.symbolic.pynac import doublefactorial
        sage: doublefactorial(-1)
        1
        sage: doublefactorial(0)
        1
        sage: doublefactorial(1)
        1
        sage: doublefactorial(5)
        15
        sage: doublefactorial(20)
        3715891200
        sage: prod( [20,18,..,2] )
        3715891200
    """
    return py_doublefactorial(n)


from sage.libs.pari.all import pari
cdef public object py_fibonacci(object n):
    return Integer(pari(n).fibonacci())

cdef public object py_step(object n):
    """
    Return step function of n.
    """
    cdef int c = cmp(n, 0)
    if c < 0:
        return ZERO
    elif c > 0:
        return ONE
    return ONE_HALF

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
    pi = constants.pi
    cdef int sgn_y = cmp(y, 0)
    cdef int sgn_x = cmp(x, 0)
    if sgn_y:
        if sgn_x > 0:
            return py_atan(abs(y/x)) * sgn_y
        elif sgn_x == 0:
            return pi/2 * sgn_y
        else:
            return (pi - py_atan(abs(y/x))) * sgn_y
    else:
        if sgn_x > 0:
            return 0
        elif x == 0:
            raise ValueError, "arctan2(0,0) undefined"
        else:
            return pi

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
        return CC(x).arcsinh()

cdef public object py_acosh(object x):
    try:
        return x.arccosh()
    except AttributeError:
        return CC(x).arccosh()


cdef public object py_atanh(object x):
    try:
        return x.arctanh()
    except AttributeError:
        return CC(x).arctanh()


cdef public object py_tgamma(object x):
    try:
        return x.gamma()
    except AttributeError:
        return RR(x).gamma()

cdef public object py_lgamma(object x):
    try:
        return x.lngamma()
    except AttributeError:
        return RR(x).lngamma()

cdef public object py_isqrt(object x):
    return Integer(x).isqrt()

cdef public object py_sqrt(object x):
    try:
        # WORRY: What if Integer's sqrt calls symbolic one and we go in circle?
        return x.sqrt()
    except AttributeError, msg:
        return sqrt(float(x))

cdef public object py_abs(object x):
    return abs(x)

cdef public object py_mod(object x, object n):
    return x % n

cdef public object py_smod(object a, object b):
    # Modulus (in symmetric representation).
    # Equivalent to Maple's mods.
    # returns a mod b in the range [-iquo(abs(b)-1,2), iquo(abs(b),2)]
    a = Integer(a); b = Integer(b)
    b = abs(b)
    c = a % b
    if c > b//2:
        c -= b
    return c

cdef public object py_irem(object x, object n):
    return Integer(x) % Integer(n)

cdef public object py_iquo(object x, object n):
    return Integer(x)//Integer(n)

cdef public object py_iquo2(object x, object n):
    x = Integer(x); n = Integer(n)
    try:
        q = x//n
        r = x - q*n
        return q, r
    except (TypeError, ValueError):
        return 0, 0

cdef public int py_int_length(object x) except -1:
    # Size in binary notation.  For integers, this is the smallest n >= 0 such
    # that -2^n <= x < 2^n. If x > 0, this is the unique n > 0 such that
    # 2^(n-1) <= x < 2^n.  This returns 0 if x is not an integer.
    return Integer(x).nbits()


##################################################################
# Not yet implemented
##################################################################
cdef public object py_li2(object x):
    raise NotImplementedError

cdef public object py_psi(object x):
    raise NotImplementedError

cdef public object py_psi2(object x, object y):
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
cdef Integer z = Integer(0)
cdef public object py_integer_from_long(long x):
    cdef Integer z = PY_NEW(Integer)
    mpz_init_set_si(z.value, x)
    return z

cdef public object py_integer_from_python_obj(object x):
    return Integer(x)


ZERO = ring.NSR(0)
ONE = ring.NSR(1)
ONE_HALF = ring.NSR(Rational((1,2)))


import sage.rings.integer
ginac_pyinit_Integer(sage.rings.integer.Integer)

import sage.rings.real_double
ginac_pyinit_Float(sage.rings.real_double.RDF)

from sage.rings.complex_double import CDF
ginac_pyinit_I(CDF.gen())


set_ginac_fn_serial()

