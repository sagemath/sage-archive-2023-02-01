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
from sage.rings.rational cimport Rational
from sage.rings.real_mpfr import RR, RealField_constructor
from sage.rings.complex_field import ComplexField
from sage.rings.all import CC

from sage.symbolic.expression cimport Expression, new_Expression_from_GEx
from sage.symbolic.function import get_sfunction_from_serial
from sage.symbolic.function cimport SFunction
from sage.symbolic.constants_c cimport PynacConstant

import ring

from sage.rings.integer cimport Integer

#################################################################
# Symbolic function helpers
#################################################################

# We declare the functions defined below as extern here, to prevent Cython
# from generating separate declarations for them which confuse g++
cdef extern from *:
    object exvector_to_PyTuple(GExVector seq)
    GEx pyExpression_to_ex(object res) except *
    object ex_to_pyExpression(GEx juice)
    object paramset_to_PyTuple(GParamSet juice)
    int py_get_ginac_serial()

cdef public object ex_to_pyExpression(GEx juice):
    """
    Convert given GiNaC::ex object to a python Expression instance.

    Used to pass parameters to custom power and series functions.
    """
    cdef Expression nex
    nex = <Expression>PY_NEW(Expression)
    GEx_construct_ex(&nex._gobj, juice)
    nex._parent = ring.SR
    return nex

cdef public object exvector_to_PyTuple(GExVector seq):
    """
    Converts arguments list given to a function to a PyTuple.

    Used to pass arguments to python functions assigned to be custom
    evaluation, derivative, etc. functions of symbolic functions.
    """
    from sage.symbolic.ring import SR
    return tuple([new_Expression_from_GEx(SR, seq.at(i))
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
        t = ring.SR._coerce_c(res)
    except TypeError, err:
        raise TypeError, "eval function did not return a symbolic expression or an element that can be coerced into a symbolic expression"
    return (<Expression>t)._gobj

cdef public object paramset_to_PyTuple(GParamSet s):
    """
    Converts a std::multiset<unsigned> to a PyTuple.

    Used to pass a list of parameter numbers with respect to which a function
    is differentiated to the printing functions py_print_fderivative and
    py_latex_fderivative.
    """
    cdef GParamSetIter itr = s.begin()
    res = []
    while itr.is_not_equal(s.end()):
        res.append(itr.obj())
        itr.inc()
    return res

def paramset_from_Expression(Expression e):
    """
    EXAMPLES::

        sage: from sage.symbolic.pynac import paramset_from_Expression
        sage: f = function('f')
        sage: paramset_from_Expression(f(x).diff(x))
        [0L]
    """
    return paramset_to_PyTuple(ex_to_fderivative(e._gobj).get_parameter_set())

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

cdef public int py_get_ginac_serial():
    """
    Returns the number of C++ level functions defined by GiNaC.

    EXAMPLES::

        sage: from sage.symbolic.pynac import get_ginac_serial
        sage: get_ginac_serial() >= 40
        True
    """
    global GINAC_FN_SERIAL
    return GINAC_FN_SERIAL

def get_ginac_serial():
    """
    Number of C++ level functions defined by GiNaC. (Defined mainly for testing.)

    EXAMPLES::

        sage: sage.symbolic.pynac.get_ginac_serial() >= 40
        True
    """
    return py_get_ginac_serial()

#################################################################
# Printing helpers
#################################################################

# We declare the functions defined below as extern here, to prevent Cython
# from generating separate declarations for them which confuse g++
cdef extern from *:
    stdstring* py_repr(object o, int level) except +
    stdstring* py_latex(object o, int level) except +
    stdstring* py_latex_variable(char* var_name) except +
    stdstring* py_print_function(unsigned id, object args) except +
    stdstring* py_latex_function(unsigned id, object args) except +
    stdstring* py_print_fderivative(unsigned id, object params, object args) except +
    stdstring* py_latex_fderivative(unsigned id, object params, object args) except +

##########################################################################
# Pynac's precedence levels, as extracted from the raw source code on
# 2009-05-15.  If this changes in Pynac it could cause a bug in
# printing.  But it's hardcoded in Pynac now, so there's not much to
# be done (at present).
#    Container: 10
#    Expairseq: 10
#    Relational: 20
#    Numeric: 30
#    Pseries: 38
#    Addition: 40
#    Integral: 45
#    Multiplication: 50
#    Noncummative mult: 50
#    Index: 55
#    Power: 60
#    Clifford: 65
#    Function: 70
#    Structure: 70
##########################################################################

cdef public stdstring* py_repr(object o, int level) except +:
    """
    Return string representation of o.  If level > 0, possibly put
    paranthesis around the string.
    """
    s = repr(o)
    if level >= 20:
        # s may need parens (e.g., is in an exponent), so decide if we
        # have to put paranthesis around s:
        # A regexp might seem better, but I don't think it's really faster.
        # It would be more readable. Python does the below (with in) very quickly.
        if level <= 50:
            t = s[1:]   # ignore leading minus
        else:
            t = s
        if ' ' in t or '/' in t or '+' in t or '-' in t or '*' in t or '^' in t:
            s = '(%s)'%s
    return string_from_pystr(s)

cdef public stdstring* py_latex(object o, int level) except +:
    """
    Return latex string representation of o.  If level > 0, possibly
    put paranthesis around the string.
    """
    from sage.misc.latex import latex
    s = latex(o)
    if level > 0:
        if ' ' in s or '/' in s or '+' in s or '-' in s or '*' in s or '^' in s or '\\frac' in s:
            s = '\\left(%s\\right)'%s
    return string_from_pystr(s)

cdef stdstring* string_from_pystr(object py_str):
    """
    Creates a C++ string with the same contents as the given python string.

    Used when passing string output to pynac for printing, since we don't want
    to mess with reference counts of the python objects and we cannot guarantee
    they won't be garbage collected before the output is printed.

    WARNING: You *must* call this with py_str a str type, or it will segfault.
    """
    cdef char *t_str = PyString_AsString(py_str)
    cdef Py_ssize_t slen = len(py_str)
    cdef stdstring* sout = stdstring_construct_cstr(t_str, slen)
    return sout

cdef public stdstring* py_latex_variable(char* var_name) except +:
    """
    Returns a c++ string containing the latex representation of the given
    variable name.

    Real work is done by the function sage.misc.latex.latex_variable_name.

    EXAMPLES::

        sage: from sage.symbolic.pynac import py_latex_variable_for_doctests
        sage: py_latex_variable = py_latex_variable_for_doctests

        sage: py_latex_variable('a')
        a
        sage: py_latex_variable('abc')
        \mbox{abc}
        sage: py_latex_variable('a_00')
        a_{00}
        sage: py_latex_variable('sigma_k')
        \sigma_{k}
        sage: py_latex_variable('sigma389')
        \sigma_{389}
        sage: py_latex_variable('beta_00')
        \beta_{00}
    """
    cdef Py_ssize_t slen
    from sage.misc.latex import latex_variable_name
    py_vlatex = latex_variable_name(var_name)
    return string_from_pystr(py_vlatex)

def py_latex_variable_for_doctests(x):
    """
    Internal function used so we can doctest a certain cdef'd method.

    EXAMPLES::

        sage: sage.symbolic.pynac.py_latex_variable_for_doctests('x')
        x
        sage: sage.symbolic.pynac.py_latex_variable_for_doctests('sigma')
        \sigma
    """
    assert isinstance(x, str)
    cdef stdstring* ostr = py_latex_variable(PyString_AsString(x))
    print(ostr.c_str())
    stdstring_delete(ostr)

def py_print_function_pystring(id, args, fname_paren=False):
    """
    Return a string with the representation of the symbolic function specified
    by the given id applied to args.

    INPUT:

        id --   serial number of the corresponding symbolic function
        params -- Set of parameter numbers with respect to which to take
                    the derivative.
        args -- arguments of the function.

    EXAMPLES::

        sage: from sage.symbolic.pynac import py_print_function_pystring
        sage: from sage.symbolic.function import function, get_sfunction_from_serial, get_ginac_serial
        sage: var('x,y,z')
        (x, y, z)
        sage: foo = function('foo', 2)
        sage: for i in range(get_ginac_serial(), get_ginac_serial()+50):
        ...     if get_sfunction_from_serial(i) == foo: break

        sage: get_sfunction_from_serial(i) == foo
        True
        sage: py_print_function_pystring(i, (x,y))
        'foo(x, y)'
        sage: py_print_function_pystring(i, (x,y), True)
        '(foo)(x, y)'
        sage: def my_print(*args): return "my args are: " + ', '.join(map(repr, args))
        sage: foo = function('foo', 2, print_func=my_print)
        sage: for i in range(get_ginac_serial(), get_ginac_serial()+50):
        ...     if get_sfunction_from_serial(i) == foo: break

        sage: get_sfunction_from_serial(i) == foo
        True
        sage: py_print_function_pystring(i, (x,y))
        'my args are: x, y'
    """
    cdef SFunction func = get_sfunction_from_serial(id)
    # This function is called from two places, from function::print in pynac
    # and from py_print_fderivative. function::print checks if the serial
    # belongs to a function defined at the C++ level. There are no C++ level
    # functions that return an instance of fderivative when derivated. Hence,
    # func will never be None.
    assert(func is not None)

    # if function has a custom print function call it
    if func._print_ is not None:
        res = func._print_(*args)
        # make sure the output is a string
        if res is None:
            return ""
        if not isinstance(res, str):
            return str(res)
        return res

    # otherwise use default output
    if fname_paren:
        olist = ['(', func._name, ')']
    else:
        olist = [func._name]
    olist.extend(['(', ', '.join(map(repr, args)), ')'])
    return ''.join(olist)

cdef public stdstring* py_print_function(unsigned id, object args) except +:
    return string_from_pystr(py_print_function_pystring(id, args))

def py_latex_function_pystring(id, args, fname_paren=False):
    """
    Return a string with the latex representation of the symbolic function
    specified by the given id applied to args.

    See documentation of py_print_function_pystring for more information.

    EXAMPLES::

        sage: from sage.symbolic.pynac import py_latex_function_pystring
        sage: from sage.symbolic.function import function, get_sfunction_from_serial, get_ginac_serial
        sage: var('x,y,z')
        (x, y, z)
        sage: foo = function('foo', 2)
        sage: for i in range(get_ginac_serial(), get_ginac_serial()+50):
        ...     if get_sfunction_from_serial(i) == foo: break

        sage: get_sfunction_from_serial(i) == foo
        True
        sage: py_latex_function_pystring(i, (x,y^z))
        '{\\rm foo}\\left(x, y^{z}\\right)'
        sage: py_latex_function_pystring(i, (x,y^z), True)
         '\\left({\\rm foo}\\right)\\left(x, y^{z}\\right)'

    Test latex_name::

        sage: foo = function('foo', 2, latex_name=r'\mathrm{bar}')
        sage: for i in range(get_ginac_serial(), get_ginac_serial()+50):
        ...     if get_sfunction_from_serial(i) == foo: break

        sage: get_sfunction_from_serial(i) == foo
        True
        sage: py_latex_function_pystring(i, (x,y^z))
        '\\mathrm{bar}\\left(x, y^{z}\\right)'

    Test custom func::

        sage: def my_print(*args): return "my args are: " + ', '.join(map(repr, args))
        sage: foo = function('foo', 2, print_latex_func=my_print)
        sage: for i in range(get_ginac_serial(), get_ginac_serial()+50):
        ...     if get_sfunction_from_serial(i) == foo: break

        sage: get_sfunction_from_serial(i) == foo
        True
        sage: py_latex_function_pystring(i, (x,y^z))
        'my args are: x, y^z'


    """
    cdef SFunction func = get_sfunction_from_serial(id)
    # This function is called from two places, from function::print in pynac
    # and from py_latex_fderivative. function::print checks if the serial
    # belongs to a function defined at the C++ level. There are no C++ level
    # functions that return an instance of fderivative when derivated. Hence,
    # func will never be None.
    assert(func is not None)

    # if function has a custom print function call it
    if func._print_latex_:
        res = func._print_latex_(*args)
        # make sure the output is a string
        if res is None:
            return ""
        if not isinstance(res, str):
            return str(res)
        return res

    # otherwise, use the latex name if defined
    if func._latex_name:
        name = func._latex_name
    else:
        # if latex_name is not defined, then call
        # latex_variable_name with "is_fname=True" flag
        from sage.misc.latex import latex_variable_name
        name = latex_variable_name(func._name, is_fname=True)
    if fname_paren:
        olist = [r'\left(', name, r'\right)']
    else:
        olist = [name]
    # print the arguments
    olist.extend([r'\left(', ', '.join([x._latex_() for x in args]),
        r'\right)'] )
    return ''.join(olist)

cdef public stdstring* py_latex_function(unsigned id, object args) except +:
    return string_from_pystr(py_latex_function_pystring(id, args))

cdef public stdstring* py_print_fderivative(unsigned id, object params,
        object args) except +:
    """
    Return a string with the representation of the derivative of the symbolic
    function specified by the given id, lists of params and args.

    INPUT:

        id --   serial number of the corresponding symbolic function
        params -- Set of parameter numbers with respect to which to take
                    the derivative.
        args -- arguments of the function.


    """
    ostr = ''.join(['D[', ', '.join([repr(int(x)) for x in params]), ']'])
    fstr = py_print_function_pystring(id, args, True)
    py_res = ostr + fstr
    return string_from_pystr(py_res)

def py_print_fderivative_for_doctests(id, params, args):
    """
    Used for testing a cdef'd function.

    EXAMPLES::

        sage: from sage.symbolic.pynac import py_print_fderivative_for_doctests as py_print_fderivative
        sage: var('x,y,z',ns=1)
        (x, y, z)
        sage: from sage.symbolic.function import function, get_sfunction_from_serial, get_ginac_serial
        sage: foo = function('foo', 2)
        sage: for i in range(get_ginac_serial(), get_ginac_serial()+50):
        ...     if get_sfunction_from_serial(i) == foo: break

        sage: get_sfunction_from_serial(i) == foo
        True
        sage: py_print_fderivative(i, (0, 1, 0, 1), (x, y^z))
        D[0, 1, 0, 1](foo)(x, y^z)

    Test custom print function::

        sage: def my_print(*args): return "func_with_args(" + ', '.join(map(repr, args)) +')'
        sage: foo = function('foo', 2, print_func=my_print)
        sage: for i in range(get_ginac_serial(), get_ginac_serial()+50):
        ...     if get_sfunction_from_serial(i) == foo: break

        sage: get_sfunction_from_serial(i) == foo
        True
        sage: py_print_fderivative(i, (0, 1, 0, 1), (x, y^z))
        D[0, 1, 0, 1]func_with_args(x, y^z)

    """
    cdef stdstring* ostr = py_print_fderivative(id, params, args)
    print(ostr.c_str())
    stdstring_delete(ostr)

cdef public stdstring* py_latex_fderivative(unsigned id, object params,
        object args) except +:
    """
    Return a string with the latex representation of the derivative of the
    symbolic function specified by the given id, lists of params and args.

    See documentation of py_print_fderivative for more information.

    """
    ostr = ''.join(['D[', ', '.join([repr(int(x)) for x in params]), ']'])
    fstr = py_latex_function_pystring(id, args, True)
    py_res = ostr + fstr
    return string_from_pystr(py_res)

def py_latex_fderivative_for_doctests(id, params, args):
    """
    Used internally for writing doctests for certain cdef'd functions.

    EXAMPLES::

        sage: from sage.symbolic.pynac import py_latex_fderivative_for_doctests as py_latex_fderivative
        sage: var('x,y,z',ns=1)
        (x, y, z)
        sage: from sage.symbolic.function import function, get_sfunction_from_serial, get_ginac_serial
        sage: foo = function('foo', 2)
        sage: for i in range(get_ginac_serial(), get_ginac_serial()+50):
        ...     if get_sfunction_from_serial(i) == foo: break

        sage: get_sfunction_from_serial(i) == foo
        True
        sage: py_latex_fderivative(i, (0, 1, 0, 1), (x, y^z))
        D[0, 1, 0, 1]\left({\rm foo}\right)\left(x, y^{z}\right)

    Test latex_name::

        sage: foo = function('foo', 2, latex_name=r'\mathrm{bar}')
        sage: for i in range(get_ginac_serial(), get_ginac_serial()+50):
        ...     if get_sfunction_from_serial(i) == foo: break

        sage: get_sfunction_from_serial(i) == foo
        True
        sage: py_latex_fderivative(i, (0, 1, 0, 1), (x, y^z))
        D[0, 1, 0, 1]\left(\mathrm{bar}\right)\left(x, y^{z}\right)

    Test custom func::

        sage: def my_print(*args): return "func_with_args(" + ', '.join(map(repr, args)) +')'
        sage: foo = function('foo', 2, print_latex_func=my_print)
        sage: for i in range(get_ginac_serial(), get_ginac_serial()+50):
        ...     if get_sfunction_from_serial(i) == foo: break

        sage: get_sfunction_from_serial(i) == foo
        True
        sage: py_latex_fderivative(i, (0, 1, 0, 1), (x, y^z))
        D[0, 1, 0, 1]func_with_args(x, y^z)
    """
    cdef stdstring* ostr = py_latex_fderivative(id, params, args)
    print(ostr.c_str())
    stdstring_delete(ostr)

#################################################################
# Archive helpers
#################################################################

cdef extern from *:
    stdstring* py_dumps(object o) except +
    object py_loads(object s) except +
    object py_get_sfunction_from_serial(unsigned s) except +
    unsigned py_get_serial_from_sfunction(SFunction f) except +

from sage.structure.sage_object import loads, dumps
cdef public stdstring* py_dumps(object o) except +:
    s = dumps(o, compress=False)
    # pynac archive format terminates atoms with zeroes.
    # since pickle output can break the archive format
    # we use the base64 data encoding
    import base64
    s = base64.b64encode(s)
    return string_from_pystr(s)

cdef public object py_loads(object s) except +:
    import base64
    s = base64.b64decode(s)
    return loads(s)

cdef public object py_get_sfunction_from_serial(unsigned s) except +:
    return get_sfunction_from_serial(s)

cdef public unsigned py_get_serial_from_sfunction(SFunction f) except +:
    return f._serial

#################################################################
# Modular helpers
#################################################################

from sage.structure.element cimport Element

# We declare the functions defined below as extern here, to prevent Cython
# from generating separate declarations for them which confuse g++
cdef extern from *:
    int py_get_parent_char(object o) except -1

cdef public int py_get_parent_char(object o) except -1:
    if isinstance(o, Element):
        return (<Element>o)._parent.characteristic()
    else:
        return 0

#################################################################
# power helpers
#################################################################

cdef extern from *:
    object py_rational_power_parts(object base, object exp) except +

from sage.rings.rational cimport rational_power_parts
cdef public object py_rational_power_parts(object base, object exp) except +:
    if not PY_TYPE_CHECK_EXACT(base, Rational):
        base = Rational(base)
    if not PY_TYPE_CHECK_EXACT(exp, Rational):
        exp = Rational(exp)
    res= rational_power_parts(base, exp)
    return res + (bool(res[0] == 1),)

#################################################################
# Binomial Coefficients
#################################################################


# We declare the functions defined below as extern here, to prevent Cython
# from generating separate declarations for them which confuse g++.
cdef extern from *:
    object py_binomial_int(int n, unsigned int k)
    object py_binomial(object n, object k)
    object py_gcd(object n, object k)
    object py_lcm(object n, object k)
    object py_real(object x)
    object py_imag(object x)
    object py_conjugate(object x)
    bint py_is_rational(object x)
    bint py_is_crational(object x)
    bint py_is_integer(object x)
    bint py_is_equal(object x, object y)
    bint py_is_even(object x)
    bint py_is_cinteger(object x)
    bint py_is_real(object x)
    bint py_is_prime(object x)
    object py_numer(object x)
    object py_denom(object x)
    object py_float(object x, int prec) except +
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
    object py_li(object x, object n, object prec)
    object py_li2(object x)
    object py_psi(object x)
    object py_psi2(object x, object y)
    object py_eval_constant(unsigned serial, int prec)
    object py_eval_unsigned_infinity()
    object py_eval_infinity()
    object py_eval_neg_infinity()
    object py_integer_from_long(long int)
    object py_integer_from_python_obj(object x)

    GConstant py_get_constant(char* name)


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

    EXAMPLES::

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
    if PY_TYPE_CHECK(n, Integer) and PY_TYPE_CHECK(k, Integer):
        if mpz_cmp_si((<Integer>n).value,1) == 0:
            return n
        elif mpz_cmp_si((<Integer>k).value,1) == 0:
            return k
        return n.gcd(k)

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
    if PY_TYPE_CHECK(n, Integer) and PY_TYPE_CHECK(k, Integer):
        if mpz_cmp_si((<Integer>n).value,1) == 0:
            return k
        elif mpz_cmp_si((<Integer>k).value,1) == 0:
            return n
        return n.lcm(k)
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
        0
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
    try:
        return x.real_part()
    except AttributeError:
        pass

    if isinstance(x, complex):
        return x.real

    return x # assume x is real

def py_real_for_doctests(x):
    """
    Used for doctesting py_real.

    TESTS::

        sage: from sage.symbolic.pynac import py_real_for_doctests
        sage: py_real_for_doctests(I)
        0
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
        1
        sage: py_imag(CC(1,5))
        5.00000000000000
        sage: py_imag(CC(1))
        0.000000000000000
        sage: py_imag(RR(1))
        0
        sage: py_imag(Mod(2,7))
        0

        sage: py_imag(QQ['x'].gen())
        0
    """
    try:
        return x.imag()
    except AttributeError:
        pass
    try:
        return x.imag_part()
    except AttributeError:
        pass

    if isinstance(x, complex):
        return x.imag

    return 0 # assume x is real

def py_imag_for_doctests(x):
    """
    Used for doctesting py_imag.

    TESTS::

        sage: from sage.symbolic.pynac import py_imag_for_doctests
        sage: py_imag_for_doctests(I)
        1
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

cdef public bint py_is_rational(object x):
    return PY_TYPE_CHECK_EXACT(x, Rational) or \
           PY_TYPE_CHECK_EXACT(x, Integer) or\
           IS_INSTANCE(x, int) or IS_INSTANCE(x, long)

cdef public bint py_is_equal(object x, object y):
    """
    Return True precisely if x and y are equal.
    """
    return bool(x==y)

cdef public bint py_is_integer(object x):
    r"""
    Returns True if pynac should treat this object as an integer.

    EXAMPLES::

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
    return IS_INSTANCE(x, int) or IS_INSTANCE(x, long) or PY_TYPE_CHECK(x, Integer) or \
           (IS_INSTANCE(x, Element) and
            ((<Element>x)._parent.is_exact() or (<Element>x)._parent == ring.SR) and
            (x in ZZ))

def py_is_integer_for_doctests(x):
    """
    Used internally for doctesting purposes.

    TESTS::

        sage: sage.symbolic.pynac.py_is_integer_for_doctests(1r)
        True
        sage: sage.symbolic.pynac.py_is_integer_for_doctests(1/3)
        False
        sage: sage.symbolic.pynac.py_is_integer_for_doctests(2)
        True
    """
    return py_is_integer(x)

cdef public bint py_is_even(object x):
    try:
        return not(x%2)
    except:
        try:
            return not(ZZ(x)%2)
        except:
            pass
    return 0


cdef public bint py_is_crational(object x):
    if py_is_rational(x):
        return True
    elif isinstance(x, Element) and (<Element>x)._parent is pynac_I._parent:
        return True
    else:
        return False

def py_is_crational_for_doctest(x):
    """
    Returns True if pynac should treat this object as an element of `\QQ(i)`.

    TESTS::

        sage: from sage.symbolic.pynac import py_is_crational_for_doctest
        sage: py_is_crational_for_doctest(1)
        True
        sage: py_is_crational_for_doctest(-2r)
        True
        sage: py_is_crational_for_doctest(1.5)
        False
        sage: py_is_crational_for_doctest(I.pyobject())
        True
        sage: py_is_crational_for_doctest(I.pyobject()+1/2)
        True
    """
    return py_is_crational(x)

cdef public bint py_is_real(object a):
    if PyInt_CheckExact(a) or PY_TYPE_CHECK(a, Integer) or PyLong_CheckExact(a):
        return True
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

cdef public bint py_is_cinteger(object x):
    return py_is_integer(x) or (py_is_crational(x) and py_denom(x) == 1)

def py_is_cinteger_for_doctest(x):
    """
    Returns True if pynac should treat this object as an element of `\ZZ(i)`.

    TESTS::

        sage: from sage.symbolic.pynac import py_is_cinteger_for_doctest
        sage: py_is_cinteger_for_doctest(1)
        True
        sage: py_is_cinteger_for_doctest(long(-3))
        True
        sage: py_is_cinteger_for_doctest(I.pyobject())
        True
        sage: py_is_cinteger_for_doctest(I.pyobject() - 3)
        True
        sage: py_is_cinteger_for_doctest(I.pyobject() + 1/2)
        False
    """
    return py_is_cinteger(x)

cdef public object py_float(object n, int prec) except +:
    """
    Evaluate pynac numeric objects numerically.

    TESTS::

        sage: from sage.symbolic.pynac import py_float_for_doctests as py_float
        sage: py_float(I, 10)
        1.0*I
        sage: py_float(pi, 100)
        3.1415926535897932384626433833
        sage: py_float(CDF(10), 53)
        10.0
        sage: type(py_float(CDF(10), 53))
        <type 'sage.rings.complex_double.ComplexDoubleElement'>
        sage: py_float(CC(1/2), 53)
        0.500000000000000
        sage: type(py_float(CC(1/2), 53))
        <type 'sage.rings.complex_number.ComplexNumber'>
    """
    if isinstance(n, (int, long, float)):
        return RealField_constructor(prec)(n)
    from sage.structure.coerce import parent
    if hasattr(parent(n), 'to_prec'):
        return parent(n).to_prec(prec)(n)
    from sage.misc.functional import numerical_approx
    return numerical_approx(n, prec)

def py_float_for_doctests(n, prec):
    """
    This function is for testing py_float.

    EXAMPLES::

        sage: from sage.symbolic.pynac import py_float_for_doctests
        sage: py_float_for_doctests(pi, 80)
        3.1415926535897932384626
    """
    return py_float(n, prec)

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

    EXAMPLES::

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
    """
    Return the value of the zeta function at the given value.

    The value is expected to be a numerical object, in RR, CC, RDF or CDF,
    different from 1.

    TESTS::

        sage: from sage.symbolic.pynac import py_zeta_for_doctests as py_zeta
        sage: py_zeta(CC.0)
        0.00330022368532410 - 0.418155449141322*I
        sage: py_zeta(CDF(5))
        1.03692775514
        sage: py_zeta(RealField(100)(5))
        1.0369277551433699263313654865
    """
    try:
        return x.zeta()
    except AttributeError:
        pass
    # zeta() evaluates its argument before calling this, so we should never
    # end up here. Precision is not a problem here, since zeta() passes that
    # on when it calls evalf() on the argument.
    raise RuntimeError, "py_zeta() received non evaluated argument"

def py_zeta_for_doctests(x):
    """
    This function is for testing py_zeta().

    EXAMPLES::

        sage: from sage.symbolic.pynac import py_zeta_for_doctests
        sage: py_zeta_for_doctests(CC.0)
        0.00330022368532410 - 0.418155449141322*I
    """
    return py_zeta(x)

cdef public object py_exp(object x):
    """
    Return the value of the exp function at the given value.

    The value is expected to be a numerical object, in RR, CC, RDF or CDF.

    TESTS::

        sage: from sage.symbolic.pynac import py_exp_for_doctests as py_exp
        sage: py_exp(CC(1))
        2.71828182845905
        sage: py_exp(CC(.5*I))
        0.877582561890373 + 0.479425538604203*I
    """
    try:
        return x.exp()
    except AttributeError:
        pass
    # exp() evaluates its argument before calling this, so we should never
    # end up here. Precision is not a problem here, since exp() passes that
    # on when it calls evalf() on the argument.
    raise RuntimeError, "py_exp() received non evaluated argument"

def py_exp_for_doctests(x):
    """
    This function tests py_exp.

    EXAMPLES::

        sage: from sage.symbolic.pynac import py_exp_for_doctests
        sage: py_exp_for_doctests(CC(2))
        7.38905609893065
    """
    return py_exp(x)

cdef public object py_log(object x):
    """
    Return the value of the log function at the given value.

    The value is expected to be a numerical object, in RR, CC, RDF or CDF.

    TESTS::

        sage: from sage.symbolic.pynac import py_log_for_doctests as py_log
        sage: py_log(CC(e))
        1.00000000000000
        sage: py_log(CC.0)
        1.57079632679490*I
    """
    try:
        return x.log()
    except AttributeError:
        pass
    # log() evaluates its argument before calling this, so we should never
    # end up here. Precision is not a problem here, since log() passes that
    # on when it calls evalf() on the argument.
    raise RuntimeError, "py_log() received non evaluated argument"

def py_log_for_doctests(x):
    """
    This function tests py_log.

    EXAMPLES::

        sage: from sage.symbolic.pynac import py_log_for_doctests
        sage: py_log_for_doctests(CC(e))
        1.00000000000000
    """
    return py_log(x)

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
    from sage.symbolic.constants import pi
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
    """
    Return x mod n. Both x and n are assumed to be integers.

    EXAMPLES::

        sage: from sage.symbolic.pynac import py_mod_for_doctests as py_mod
        sage: py_mod(I.parent(5), 4)
        1
        sage: py_mod(3, -2)
        -1
        sage: py_mod(3/2, 2)
        Traceback (most recent call last):
        ...
        TypeError: no conversion of this rational to integer


    Note: The original code for this function in GiNaC checks if the arguments
    are integers, and returns 0 otherwise. We omit this check, since all the
    calls to py_mod are preceeded by an integer check. We also raise an error
    if converting the arguments to integers fails, since silently returning 0
    would hide possible misuses of this function.

    Regarding the sign of the return value, the CLN reference manual says:

        If x and y are both >= 0, mod(x,y) = rem(x,y) >= 0. In general,
        mod(x,y) has the sign of y or is zero, and rem(x,y) has the sign of
        x or is zero.

    This matches the behavior of the % operator for integers in Sage.
    """
    return Integer(x) % Integer(n)

def py_mod_for_doctests(x, n):
    """
    This function is a python wrapper so py_mod can be tested. The real tests
    are in the docstring for py_mod.

    EXAMPLES::

        sage: from sage.symbolic.pynac import py_mod_for_doctests
        sage: py_mod_for_doctests(5, 2)
        1
    """
    return py_mod(x, n)

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


cdef public object py_li(object x, object n, object prec):
    """
    Returns a numerical approximation of polylog(n, x) with *prec*
    bits of precision.
    """
    #This is really awful and slow.  We should really improve our
    #interface with mpmath! :-)
    from sympy import mpmath
    old_prec, mpmath.mp.prec = mpmath.mp.prec, prec
    n = int(n)
    res = mpmath.polylog(n, str(x))
    try:
        res = RealField_constructor(prec)(str(res))
    except TypeError:
        res = ComplexField(prec)(str(res))
    mpmath.mp.prec = old_prec
    return res

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

cdef public GConstant py_get_constant(char* name):
    """
    Returns a constant given its name. This is called by
    constant::unarchive in constant.cpp in Pynac and is used for
    pickling.
    """
    from sage.symbolic.constants import constants_name_table
    cdef PynacConstant pc
    c = constants_name_table.get(name, None)
    if c is None:
        raise RuntimeError
    else:
        pc = c._pynac
        return pc.object

cdef public object py_eval_constant(unsigned serial, int prec):
    from sage.symbolic.constants import constants_table
    constant = constants_table[serial]
    try:
        return RealField_constructor(prec)(constant)
    except TypeError:
        return ComplexField(prec)(constant)

cdef public object py_eval_unsigned_infinity():
    """
    Returns unsigned_infinity.
    """
    from sage.rings.infinity import unsigned_infinity
    return unsigned_infinity

def py_eval_unsigned_infinity_for_doctests():
    """
    This function tests py_eval_unsigned_infinity.

    TESTS::
        sage: from sage.symbolic.pynac import py_eval_unsigned_infinity_for_doctests as py_eval_unsigned_infinity
        sage: py_eval_unsigned_infinity()
        Infinity
    """
    return py_eval_unsigned_infinity()

cdef public object py_eval_infinity():
    """
    Returns positive infinity, i.e., oo.
    """
    from sage.rings.infinity import infinity
    return infinity

def py_eval_infinity_for_doctests():
    """
    This function tests py_eval_infinity.

    TESTS::
        sage: from sage.symbolic.pynac import py_eval_infinity_for_doctests as py_eval_infinity
        sage: py_eval_infinity()
        +Infinity
    """
    return py_eval_infinity()

cdef public object py_eval_neg_infinity():
    """
    Returns minus_infinity.
    """
    from sage.rings.infinity import minus_infinity
    return minus_infinity

def py_eval_neg_infinity_for_doctests():
    """
    This function tests py_eval_neg_infinity.

    TESTS::
        sage: from sage.symbolic.pynac import py_eval_neg_infinity_for_doctests as py_eval_neg_infinity
        sage: py_eval_neg_infinity()
        -Infinity
    """
    return py_eval_neg_infinity()

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


ZERO = ring.SR(0)
ONE = ring.SR(1)
ONE_HALF = ring.SR(Rational((1,2)))

symbol_table = {'functions':{}}
def register_symbol(obj, conversions):
    """
    Add an object to the symbol table, along with how to convert it to
    other systems such as Maxima, Mathematica, etc.  This table is used
    to convert *from* other systems back to Sage.

    INPUTS:

        - `obj` -- a symbolic object or function.

        - `conversions` -- a dictionary of conversions, where the keys
                           are the names of interfaces (e.g.,
                           'maxima'), and the values are the string
                           representation of obj in that system.



    EXAMPLES::

        sage: sage.symbolic.pynac.register_symbol(SR(5),{'maxima':'five'})
        sage: SR(maxima_calculus('five'))
        5
    """
    conversions = dict(conversions)
    try:
        conversions['sage'] = obj.name()
    except AttributeError:
        pass
    for system, value in conversions.iteritems():
        system_table = symbol_table.get(system, None)
        if system_table is None:
            symbol_table[system] = system_table = {}
        system_table[value] = obj



import sage.rings.integer
ginac_pyinit_Integer(sage.rings.integer.Integer)

import sage.rings.real_double
ginac_pyinit_Float(sage.rings.real_double.RDF)

cdef Element pynac_I
I = None

def init_pynac_I():
    """
    Initialize the numeric I object in pynac. We use the generator of QQ(i).

    EXAMPLES::

        sage: sage.symbolic.pynac.init_pynac_I()
        sage: type(sage.symbolic.pynac.I)
        <type 'sage.symbolic.expression.Expression'>
        sage: type(sage.symbolic.pynac.I.pyobject())
        <type 'sage.rings.number_field.number_field_element_quadratic.NumberFieldElement_quadratic'>
    """
    global pynac_I, I
    from sage.rings.number_field.number_field import QuadraticField
    K = QuadraticField(-1, 'I')
    pynac_I = K.gen()
    ginac_pyinit_I(pynac_I)
    I = new_Expression_from_GEx(ring.SR, g_I)

init_pynac_I()


set_ginac_fn_serial()
