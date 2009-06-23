r"""
Symbolic Computation.

AUTHORS:

- Bobby Moretti and William Stein (2006-2007)

The Sage calculus module is loosely based on the Sage Enhancement
Proposal found at: http://www.sagemath.org:9001/CalculusSEP.

EXAMPLES:

The basic units of the calculus package are symbolic expressions which
are elements of the symbolic expression ring (SR). To create a
symbolic variable object in Sage, use the :func:`var` function, whose
argument is the text of that variable. Note that Sage is intelligent
about LaTeXing variable names.

::

    sage: x1 = var('x1'); x1
    x1
    sage: latex(x1)
    x_{1}
    sage: theta = var('theta'); theta
    theta
    sage: latex(theta)
    \theta

Sage predefines ``x`` to be a global indeterminate.
Thus the following works::

    sage: x^2
    x^2
    sage: type(x)
    <type 'sage.symbolic.expression.Expression'>

More complicated expressions in Sage can be built up using ordinary
arithmetic. The following are valid, and follow the rules of Python
arithmetic: (The '=' operator represents assignment, and not
equality)

::

    sage: var('x,y,z')
    (x, y, z)
    sage: f = x + y + z/(2*sin(y*z/55))
    sage: g = f^f; g
    (x + y + 1/2*z/sin(1/55*y*z))^(x + y + 1/2*z/sin(1/55*y*z))

Differentiation and integration are available, but behind the
scenes through maxima::

    sage: f = sin(x)/cos(2*y)
    sage: f.derivative(y)
    2*sin(x)*sin(2*y)/cos(2*y)^2
    sage: g = f.integral(x); g
    -cos(x)/cos(2*y)

Note that these methods require an explicit variable name. If none
is given, Sage will try to find one for you.

::

    sage: f = sin(x); f.derivative()
    cos(x)

However when this is ambiguous, Sage will raise an exception::

    sage: f = sin(x+y); f.derivative()
    Traceback (most recent call last):
    ...
    ValueError: No differentiation variable specified.

Substitution works similarly. We can substitute with a python
dict::

    sage: f = sin(x*y - z)
    sage: f({x: var('t'), y: z})
    sin(t*z - z)

Also we can substitute with keywords::

    sage: f = sin(x*y - z)
    sage: f(x = t, y = z)
    sin(t*z - z)

It was formerly the case that if there was no ambiguity of variable
names, we didn't have to specify them; that still works for the moment,
but the behavior is deprecated::

    sage: f = sin(x)
    sage: f(y)
    doctest:...: DeprecationWarning: Substitution using function-call syntax and unnamed arguments is deprecated and will be removed from a future release of Sage; you can use named arguments instead, like EXPR(x=..., y=...)
    sin(y)
    sage: f(pi)
    0

However if there is ambiguity, we should explicitly state what
variables we're substituting for::

    sage: f = sin(2*pi*x/y)
    sage: f(x=4)
    sin(8*pi/y)

We can also make a ``CallableSymbolicExpression``,
which is a ``SymbolicExpression`` that is a function of
specified variables in a fixed order. Each
``SymbolicExpression`` has a
``function(...)`` method that is used to create a
``CallableSymbolicExpression``, as illustrated below::

    sage: u = log((2-x)/(y+5))
    sage: f = u.function(x, y); f
    (x, y) |--> log(-(x - 2)/(y + 5))

There is an easier way of creating a
``CallableSymbolicExpression``, which relies on the
Sage preparser.

::

    sage: f(x,y) = log(x)*cos(y); f
    (x, y) |--> log(x)*cos(y)

Then we have fixed an order of variables and there is no ambiguity
substituting or evaluating::

    sage: f(x,y) = log((2-x)/(y+5))
    sage: f(7,t)
    log(-5/(t + 5))

Some further examples::

    sage: f = 5*sin(x)
    sage: f
    5*sin(x)
    sage: f(x=2)
    5*sin(2)
    sage: f(x=pi)
    0
    sage: float(f(x=pi))
    0.0

Another example::

    sage: f = integrate(1/sqrt(9+x^2), x); f
    arcsinh(1/3*x)
    sage: f(x=3)
    arcsinh(1)
    sage: f.derivative(x)
    1/3/sqrt(1/9*x^2 + 1)

We compute the length of the parabola from 0 to 2::

    sage: x = var('x')
    sage: y = x^2
    sage: dy = derivative(y,x)
    sage: z = integral(sqrt(1 + dy^2), x, 0, 2)
    sage: z
    sqrt(17) + 1/4*arcsinh(4)
    sage: n(z,200)
    4.6467837624329358733826155674904591885104869874232887508703
    sage: float(z)
    4.6467837624329356

We test pickling::

    sage: x, y = var('x,y')
    sage: f = -sqrt(pi)*(x^3 + sin(x/cos(y)))
    sage: bool(loads(dumps(f)) == f)
    True

Coercion examples:

We coerce various symbolic expressions into the complex numbers::

    sage: CC(I)
    1.00000000000000*I
    sage: CC(2*I)
    2.00000000000000*I
    sage: ComplexField(200)(2*I)
    2.0000000000000000000000000000000000000000000000000000000000*I
    sage: ComplexField(200)(sin(I))
    1.1752011936438014568823818505956008151557179813340958702296*I
    sage: f = sin(I) + cos(I/2); f
    sin(I) + cos(1/2*I)
    sage: CC(f)
    1.12762596520638 + 1.17520119364380*I
    sage: ComplexField(200)(f)
    1.1276259652063807852262251614026720125478471180986674836290 + 1.1752011936438014568823818505956008151557179813340958702296*I
    sage: ComplexField(100)(f)
    1.1276259652063807852262251614 + 1.1752011936438014568823818506*I

We illustrate construction of an inverse sum where each denominator
has a new variable name::

    sage: f = sum(1/var('n%s'%i)^i for i in range(10))
    sage: f
    1/n1 + 1/n2^2 + 1/n3^3 + 1/n4^4 + 1/n5^5 + 1/n6^6 + 1/n7^7 + 1/n8^8 + 1/n9^9 + 1

Note that after calling var, the variables are immediately
available for use::

    sage: (n1 + n2)^5
    (n1 + n2)^5

We can, of course, substitute::

    sage: f(n9=9,n7=n6)
    1/n1 + 1/n2^2 + 1/n3^3 + 1/n4^4 + 1/n5^5 + 1/n6^6 + 1/n6^7 + 1/n8^8 + 387420490/387420489

TESTS:

Substitution::

    sage: f = x
    sage: f(x=5)
    5

Simplifying expressions involving scientific notation::

    sage: k = var('k')
    sage: a0 = 2e-06; a1 = 12
    sage: c = a1 + a0*k; c
    (2.00000000000000e-6)*k + 12
    sage: sqrt(c)
    sqrt((2.00000000000000e-6)*k + 12)
    sage: sqrt(c^3)
    sqrt(((2.00000000000000e-6)*k + 12)^3)

The symbolic Calculus package uses its own copy of maxima for
simplification, etc., which is separate from the default
system-wide version::

    sage: maxima.eval('[x,y]: [1,2]')
    '[1,2]'
    sage: maxima.eval('expand((x+y)^3)')
    '27'

If the copy of maxima used by the symbolic calculus package were
the same as the default one, then the following would return 27,
which would be very confusing indeed!

::

    sage: x, y = var('x,y')
    sage: expand((x+y)^3)
    x^3 + 3*x^2*y + 3*x*y^2 + y^3

Set x to be 5 in maxima::

    sage: maxima('x: 5')
    5
    sage: maxima('x + x + %pi')
    %pi+10

This simplification is done using maxima (behind the scenes)::

    sage: x + x + pi
     pi + 2*x

Note that ``x`` is still ``x``, since the
maxima used by the calculus package is different than the one in
the interactive interpreter.

Check to see that the problem with the variables method mentioned
in Trac ticket #3779 is actually fixed::

    sage: f = function('F',x)
    sage: diff(f*SR(1),x)
    D[0](F)(x)
"""

import weakref
import re
from sage.rings.all import (CommutativeRing, RealField, is_Polynomial,
                            is_MPolynomial, is_MPolynomialRing, is_FractionFieldElement,
                            is_RealNumber, is_ComplexNumber, RR, is_NumberFieldElement,
                            Integer, Rational, CC, QQ, CDF, ZZ,
                            QuadDoubleElement,
                            PolynomialRing, ComplexField, RealDoubleElement,
                            algdep, Integer, RealNumber, RealIntervalField)

from sage.rings.real_mpfr import create_RealNumber

from sage.structure.element import RingElement, is_Element
from sage.structure.parent_base import ParentWithBase

import operator
from sage.misc.latex import latex, latex_variable_name
from sage.misc.misc import uniq as unique
from sage.structure.sage_object import SageObject

from sage.interfaces.maxima import Maxima

import sage.numerical.optimize

from sage.misc.parser import Parser

import math

import sage.ext.fast_eval as fast_float

from sage.symbolic.ring import var, SR, is_SymbolicVariable
from sage.symbolic.expression import Expression
from sage.symbolic.function import SFunction, PrimitiveFunction
from sage.symbolic.pynac import symbol_table

arc_functions =  ['asin', 'acos', 'atan', 'asinh', 'acosh', 'atanh', 'acoth', 'asech', 'acsch', 'acot', 'acsc', 'asec']

maxima = Maxima(init_code = ['display2d:false; domain: complex; keepfloat: true; load(topoly_solver)'],
                script_subdirectory=None)

########################################################


###################################################################
# integration
###################################################################
def integral(expression, v=None, a=None, b=None, algorithm='maxima'):
    r"""
    Returns the indefinite integral with respect to the variable
    `v`, ignoring the constant of integration. Or, if endpoints
    `a` and `b` are specified, returns the definite
    integral over the interval `[a, b]`.

    If ``self`` has only one variable, then it returns the
    integral with respect to that variable.

    INPUT:


    -  ``v`` - (optional) a variable or variable name

    -  ``a`` - (optional) lower endpoint of definite
       integral

    -  ``b`` - (optional) upper endpoint of definite
       integral

    - ``algorithm`` - (default: 'maxima')  one of

              - 'maxima' - use maxima (the default)

              - 'sympy' - use sympy (also in Sage)

              - 'mathematica_free' - use http://integrals.wolfram.com/

    EXAMPLES::

        sage: x = var('x')
        sage: h = sin(x)/(cos(x))^2
        sage: h.integral(x)
        1/cos(x)

    ::

        sage: f = x^2/(x+1)^3
        sage: f.integral()
        1/2*(4*x + 3)/(x^2 + 2*x + 1) + log(x + 1)

    ::

        sage: f = x*cos(x^2)
        sage: f.integral(x, 0, sqrt(pi))
        0
        sage: f.integral(a=-pi, b=pi)
        0

    ::

        sage: f(x) = sin(x)
        sage: f.integral(x, 0, pi/2)
        1

    The variable and endpoints are both optional::

        sage: integral(sin(x))
        -cos(x)
        sage: integral(sin(x), var('y'))
        y*sin(x)
        sage: integral(sin(x), pi, 2*pi)
        -2
        sage: integral(sin(x), var('y'), pi, 2*pi)
        pi*sin(x)

    Constraints are sometimes needed::

        sage: var('x, n')
        (x, n)
        sage: integral(x^n,x)
        Traceback (most recent call last):
        ...
        TypeError: Computation failed since Maxima requested additional constraints (try the command 'assume(n+1>0)' before integral or limit evaluation, for example):
        Is  n+1  zero or nonzero?
        sage: assume(n > 0)
        sage: integral(x^n,x)
        x^(n + 1)/(n + 1)
        sage: forget()

    Usually the constraints are of sign, but others are possible::

        sage: assume(n==-1)
        sage: integral(x^n,x)
        log(x)

    Note that an exception is raised when a definite integral is
    divergent.

    ::

        sage: forget()

        sage: integrate(1/x^3,x,0,1)
        Traceback (most recent call last):
        ...
        ValueError: Integral is divergent.
        sage: integrate(1/x^3,x,-1,3)
        Traceback (most recent call last):
        ...
        ValueError: Integral is divergent.

    .. note::

       Above, putting assume(n == -1) does not yield the right
       behavior.

    The examples in the Maxima documentation::

        sage: var('x, y, z, b')
        (x, y, z, b)
        sage: integral(sin(x)^3)
        1/3*cos(x)^3 - cos(x)
        sage: integral(x/sqrt(b^2-x^2))
        x*log(2*b + 2*sqrt(b^2 - x^2))
        sage: integral(x/sqrt(b^2-x^2), x)
        -sqrt(b^2 - x^2)
        sage: integral(cos(x)^2 * exp(x), x, 0, pi)
        3/5*e^pi - 3/5
        sage: integral(x^2 * exp(-x^2), x, -oo, oo)
        1/2*sqrt(pi)

    We integrate the same function in both Mathematica and Sage (via
    Maxima)::

        sage: _ = var('x, y, z')
        sage: f = sin(x^2) + y^z
        sage: g = mathematica(f)                           # optional  -- requires mathematica
        sage: print g                                      # optional -- requires mathematica
                  z        2
                 y  + Sin[x ]
        sage: print g.Integrate(x)                         # optional -- requires mathematica
                    z        Pi                2
                 x y  + Sqrt[--] FresnelS[Sqrt[--] x]
                             2                 Pi
        sage: print f.integral(x)
        y^z*x + 1/8*((I - 1)*sqrt(2)*erf((1/2*I - 1/2)*sqrt(2)*x) + (I + 1)*sqrt(2)*erf((1/2*I + 1/2)*sqrt(2)*x))*sqrt(pi)

    Alternatively, just use algorithm='mathematica_free' to integrate via Mathematica
    over the internet (deos NOT require a mathematica license!)::

        sage: _ = var('x, y, z')
        sage: f = sin(x^2) + y^z
        sage: f.integrate(algorithm="mathematica_free")       # optional -- requires internet
        sqrt(pi)*sqrt(1/2)*fresnels(sqrt(2)*x/sqrt(pi)) + y^z*x

    We can also use Sympy::

        sage: _ = var('x, y, z')
        sage: (x^y-z).integrate(y)
        -y*z + x^y/log(x)
        sage: (x^y-z).integrate(y,algorithm="sympy")
        -y*z + x^y/log(x)


    We integrate the above function in maple now::

        sage: g = maple(f); g                             # optional -- requires maple
        sin(x^2)+y^z
        sage: g.integrate(x)                              # optional -- requires maple
        1/2*2^(1/2)*Pi^(1/2)*FresnelS(2^(1/2)/Pi^(1/2)*x)+y^z*x

    We next integrate a function with no closed form integral. Notice
    that the answer comes back as an expression that contains an
    integral itself.

    ::

        sage: A = integral(1/ ((x-4) * (x^3+2*x+1)), x); A
        1/73*log(x - 4) - 1/73*integrate((x^2 + 4*x + 18)/(x^3 + 2*x + 1), x)

    We now show that floats are not converted to rationals
    automatically since we by default have keepfloat: true in maxima.

    ::

        sage: integral(e^(-x^2),x, 0, 0.1)
        0.0562314580091*sqrt(pi)

    ALIASES: integral() and integrate() are the same.

    EXAMPLES: Here is example where we have to use assume::

        sage: a,b = var('a,b')
        sage: integrate(1/(x^3 *(a+b*x)^(1/3)), x)
        Traceback (most recent call last):
        ...
        TypeError: Computation failed since Maxima requested additional constraints (try the command 'assume(a>0)' before integral or limit evaluation, for example):
        Is  a  positive or negative?

    So we just assume that `a>0` and the integral works::

        sage: assume(a>0)
        sage: integrate(1/(x^3 *(a+b*x)^(1/3)), x)
        2/9*sqrt(3)*b^2*arctan(1/3*(2*(b*x + a)^(1/3) + a^(1/3))*sqrt(3)/a^(1/3))/a^(7/3) + 2/9*b^2*log((b*x + a)^(1/3) - a^(1/3))/a^(7/3) - 1/9*b^2*log((b*x + a)^(2/3) + (b*x + a)^(1/3)*a^(1/3) + a^(2/3))/a^(7/3) + 1/6*(4*(b*x + a)^(5/3)*b^2 - 7*(b*x + a)^(2/3)*a*b^2)/((b*x + a)^2*a^2 - 2*(b*x + a)*a^3 + a^4)

    TESTS:

    The following integral was broken prior to Maxima 5.15.0 -
    see #3013

    ::

        sage: integrate(sin(x)*cos(10*x)*log(x))
        1/18*log(x)*cos(9*x) - 1/22*log(x)*cos(11*x) - 1/18*integrate(cos(9*x)/x, x) + 1/22*integrate(cos(11*x)/x, x)


    It is no longer possible to use certain functions without an
    explicit variable.  Instead, evaluate the function at a variable,
    and then take the integral::

        sage: integrate(sin)
        Traceback (most recent call last):
        ...
        TypeError

        sage: integrate(sin(x))
        -cos(x)
        sage: integrate(sin(x), 0, 1)
        -cos(1) + 1
    """
    if b is None and a is not None:
        # two arguments, must be endpoints
        a, b = v, a
        v = None

    if v is None:
        v = expression.default_variable()
        if isinstance(expression, SFunction):
            # a bare function like sin
            expression = expression(v)

    elif not is_SymbolicVariable(v):
        v = var(repr(v))
        #raise TypeError, 'must integrate with respect to a variable'

    if (a is None) ^ (b is None):
        raise TypeError, 'only one endpoint given'

    if algorithm == 'maxima':
        if a is None:
            result = expression._maxima_().integrate(v)
        else:
            try:
                result = expression._maxima_().integrate(v, a, b)
            except TypeError, error:
                s = str(error)
                if "divergent" in s or 'Principal Value' in s:
                    raise ValueError, "Integral is divergent."
                else:
                    raise

    elif algorithm == 'mathematica_free':
        import urllib, re
        # We need to integrate against x
        vars = [str(x) for x in expression.variables()]
        if any(len(x)>1 for x in vars):
            raise NotImplementedError, "Mathematica online integrator can only handle single letter variables."
        x = var('x')
        if repr(v) != 'x':
            for i in range(ord('a'), ord('z')+1):
                if chr(i) not in vars:
                    shadow_x = var(chr(i))
                    break
            expression = expression.subs({x:shadow_x}).subs({dvar: x})
        params = urllib.urlencode({'expr': expression._mathematica_init_(), 'random': 'false'})
        page = urllib.urlopen("http://integrals.wolfram.com/index.jsp", params).read()
        page = page[page.index('"inputForm"'):page.index('"outputForm"')]
        page = re.sub("\s", "", page)
        mexpr = re.match(r".*Integrate.*==</em><br/>(.*)</p>", page).groups()[0]
        try:
            ans = SR(mexpr.lower().replace('[', '(').replace(']', ')'))
            if repr(v) != 'x':
                ans = ans.subs({x:v}).subs({shadow_x:x})
            return ans
        except TypeError:
            raise ValueError, "Unable to parse: %s" % mexpr

    elif algorithm == 'sympy':
        import sympy
        ex = expression._sympy_()
        v = v._sympy_()
        if a is None:
            result = sympy.integrate(ex, v)
        else:
            result = sympy.integrate(ex, (v, a._sympy_(), b._sympy_()))

    else:
        raise ValueError, "unknown algorithm: %s" % algorithm

    if a is None:
        return expression.parent()(result)
    else:
        return SR(result)

integrate = integral

def nintegral(ex, x, a, b,
              desired_relative_error='1e-8',
              maximum_num_subintervals=200):
    r"""
    Return a floating point machine precision numerical approximation
    to the integral of ``self`` from `a` to
    `b`, computed using floating point arithmetic via maxima.

    INPUT:


    -  ``x`` - variable to integrate with respect to

    -  ``a`` - lower endpoint of integration

    -  ``b`` - upper endpoint of integration

    -  ``desired_relative_error`` - (default: '1e-8') the
       desired relative error

    -  ``maximum_num_subintervals`` - (default: 200)
       maxima number of subintervals


    OUTPUT:


    -  float: approximation to the integral

    -  float: estimated absolute error of the
       approximation

    -  the number of integrand evaluations

    -  an error code:

       -  ``0`` - no problems were encountered

       -  ``1`` - too many subintervals were done

       -  ``2`` - excessive roundoff error

       -  ``3`` - extremely bad integrand behavior

       -  ``4`` - failed to converge

       -  ``5`` - integral is probably divergent or slowly
          convergent

       -  ``6`` - the input is invalid


    ALIAS: nintegrate is the same as nintegral

    REMARK: There is also a function
    ``numerical_integral`` that implements numerical
    integration using the GSL C library. It is potentially much faster
    and applies to arbitrary user defined functions.

    Also, there are limits to the precision to which Maxima can compute
    the integral to due to limitations in quadpack.

    ::

        sage: f = x
        sage: f = f.nintegral(x,0,1,1e-14)
        Traceback (most recent call last):
        ...
        ValueError: Maxima (via quadpack) cannot compute the integral to that precision

    EXAMPLES::

        sage: f(x) = exp(-sqrt(x))
        sage: f.nintegral(x, 0, 1)
        (0.52848223531423055, 4.163...e-11, 231, 0)

    We can also use the ``numerical_integral`` function,
    which calls the GSL C library.

    ::

        sage: numerical_integral(f, 0, 1)
        (0.52848223225314706, 6.83928460...e-07)

    Note that in exotic cases where floating point evaluation of the
    expression leads to the wrong value, then the output can be
    completely wrong::

        sage: f = exp(pi*sqrt(163)) - 262537412640768744

    Despite appearance, `f` is really very close to 0, but one
    gets a nonzero value since the definition of
    ``float(f)`` is that it makes all constants inside the
    expression floats, then evaluates each function and each arithmetic
    operation using float arithmetic::

        sage: float(f)
        -480.0

    Computing to higher precision we see the truth::

        sage: f.n(200)
        -7.4992740280181431112064614366622348652078895136533593355718e-13
        sage: f.n(300)
        -7.49927402801814311120646143662663009137292462589621789352095066181709095575681963967103004e-13

    Now numerically integrating, we see why the answer is wrong::

        sage: f.nintegrate(x,0,1)
        (-480.00000000000011, 5.3290705182007538e-12, 21, 0)

    It is just because every floating point evaluation of return -480.0
    in floating point.

    Important note: using GP/PARI one can compute numerical integrals
    to high precision::

        sage: gp.eval('intnum(x=17,42,exp(-x^2)*log(x))')
        '2.565728500561051482917356396 E-127'        # 32-bit
        '2.5657285005610514829173563961304785900 E-127'    # 64-bit
        sage: old_prec = gp.set_real_precision(50)
        sage: gp.eval('intnum(x=17,42,exp(-x^2)*log(x))')
        '2.5657285005610514829173563961304785900147709554020 E-127'
        sage: gp.set_real_precision(old_prec)
        57

    Note that the input function above is a string in PARI syntax.
    """
    try:
        v = ex._maxima_().quad_qags(x, a, b,
                                    epsrel=desired_relative_error,
                                    limit=maximum_num_subintervals)
    except TypeError, err:
        if "ERROR" in str(err):
            raise ValueError, "Maxima (via quadpack) cannot compute the integral to that precision"
        else:
            raise TypeError, err

    #This is just a work around until there is a response to
    #http://www.math.utexas.edu/pipermail/maxima/2008/012975.html
    if 'quad_qags' in str(v):
        raise ValueError, "Maxima (via quadpack) cannot compute the integral to that precision"

    return float(v[0]), float(v[1]), Integer(v[2]), Integer(v[3])

nintegrate = nintegral

def minpoly(ex, var='x', algorithm=None, bits=None, degree=None, epsilon=0):
    r"""
    Return the minimal polynomial of self, if possible.

    INPUT:

    -  ``var`` - polynomial variable name (default 'x')

    -  ``algorithm`` - 'algebraic' or 'numerical' (default
       both, algebraic first)

    -  ``bits`` - the number of bits to use in numerical
       approx

    -  ``degree`` - the expected algebraic degree

    -  ``epsilon`` - return without error as long as
       f(self) epsilon, in the case that the result cannot be proven.

       All of the above parameters are optional, with epsilon=0, bits and
       degree tested up to 1000 and 24 by default respectively. The
       numerical algorithm will be faster if bits and/or degree are given
       explicitly. The algebraic algorithm ignores the last three
       parameters.


    OUTPUT: The minimal polynomial of self. If the numerical algorithm
    is used then it is proved symbolically when epsilon=0 (default).

    If the minimal polynomial could not be found, two distinct kinds of
    errors are raised. If no reasonable candidate was found with the
    given bit/degree parameters, a ``ValueError`` will be
    raised. If a reasonable candidate was found but (perhaps due to
    limits in the underlying symbolic package) was unable to be proved
    correct, a ``NotImplementedError`` will be raised.

    ALGORITHM: Two distinct algorithms are used, depending on the
    algorithm parameter. By default, the algebraic algorithm is
    attempted first, then the numerical one.

    Algebraic: Attempt to evaluate this expression in QQbar, using
    cyclotomic fields to resolve exponential and trig functions at
    rational multiples of pi, field extensions to handle roots and
    rational exponents, and computing compositums to represent the full
    expression as an element of a number field where the minimal
    polynomial can be computed exactly. The bits, degree, and epsilon
    parameters are ignored.

    Numerical: Computes a numerical approximation of
    ``self`` and use PARI's algdep to get a candidate
    minpoly `f`. If `f(\mathtt{self})`,
    evaluated to a higher precision, is close enough to 0 then evaluate
    `f(\mathtt{self})` symbolically, attempting to prove
    vanishing. If this fails, and ``epsilon`` is non-zero,
    return `f` if and only if
    `f(\mathtt{self}) < \mathtt{epsilon}`.
    Otherwise raise a ``ValueError`` (if no suitable
    candidate was found) or a ``NotImplementedError`` (if a
    likely candidate was found but could not be proved correct).

    EXAMPLES: First some simple examples::

        sage: sqrt(2).minpoly()
        x^2 - 2
        sage: minpoly(2^(1/3))
        x^3 - 2
        sage: minpoly(sqrt(2) + sqrt(-1))
        x^4 - 2*x^2 + 9
        sage: minpoly(sqrt(2)-3^(1/3))
        x^6 - 6*x^4 + 6*x^3 + 12*x^2 + 36*x + 1

    Works with trig and exponential functions too.

    ::

        sage: sin(pi/3).minpoly()
        x^2 - 3/4
        sage: sin(pi/7).minpoly()
        x^6 - 7/4*x^4 + 7/8*x^2 - 7/64
        sage: minpoly(exp(I*pi/17))
        x^16 - x^15 + x^14 - x^13 + x^12 - x^11 + x^10 - x^9 + x^8 - x^7 + x^6 - x^5 + x^4 - x^3 + x^2 - x + 1

    Here we verify it gives the same result as the abstract number
    field.

    ::

        sage: (sqrt(2) + sqrt(3) + sqrt(6)).minpoly()
        x^4 - 22*x^2 - 48*x - 23
        sage: K.<a,b> = NumberField([x^2-2, x^2-3])
        sage: (a+b+a*b).absolute_minpoly()
        x^4 - 22*x^2 - 48*x - 23

    Here we solve a cubic and then recover it from its complicated
    radical expansion.

    ::

        sage: f = x^3 - x + 1
        sage: a = f.solve(x)[0].rhs(); a
        -1/2*(I*sqrt(3) + 1)*(1/18*sqrt(3)*sqrt(23) - 1/2)^(1/3) - 1/6*(-I*sqrt(3) + 1)/(1/18*sqrt(3)*sqrt(23) - 1/2)^(1/3)
        sage: a.minpoly()
        x^3 - x + 1

    Note that simplification may be necessary to see that the minimal
    polynomial is correct.

    ::

        sage: a = sqrt(2)+sqrt(3)+sqrt(5)
        sage: f = a.minpoly(); f
        x^8 - 40*x^6 + 352*x^4 - 960*x^2 + 576
        sage: f(a)
        ((((sqrt(2) + sqrt(3) + sqrt(5))^2 - 40)*(sqrt(2) + sqrt(3) + sqrt(5))^2 + 352)*(sqrt(2) + sqrt(3) + sqrt(5))^2 - 960)*(sqrt(2) + sqrt(3) + sqrt(5))^2 + 576
        sage: f(a).expand()
        0

    Here we show use of the ``epsilon`` parameter. That
    this result is actually exact can be shown using the addition
    formula for sin, but maxima is unable to see that.

    ::

        sage: a = sin(pi/5)
        sage: a.minpoly(algorithm='numerical')
        Traceback (most recent call last):
        ...
        NotImplementedError: Could not prove minimal polynomial x^4 - 5/4*x^2 + 5/16 (epsilon 0.00000000000000e-1)
        sage: f = a.minpoly(algorithm='numerical', epsilon=1e-100); f
        x^4 - 5/4*x^2 + 5/16
        sage: f(a).numerical_approx(100)
        0.00000000000000000000000000000

    The degree must be high enough (default tops out at 24).

    ::

        sage: a = sqrt(3) + sqrt(2)
        sage: a.minpoly(algorithm='numerical', bits=100, degree=3)
        Traceback (most recent call last):
        ...
        ValueError: Could not find minimal polynomial (100 bits, degree 3).
        sage: a.minpoly(algorithm='numerical', bits=100, degree=10)
        x^4 - 10*x^2 + 1

    There is a difference between algorithm='algebraic' and
    algorithm='numerical'::

        sage: cos(pi/22).minpoly(algorithm='algebraic')
        x^10 - 11/4*x^8 + 11/4*x^6 - 77/64*x^4 + 55/256*x^2 - 11/1024
        sage: cos(pi/22).minpoly(algorithm='numerical')
        Traceback (most recent call last):
        NotImplementedError: Could not prove minimal polynomial x^10 - 11/4*x^8 + 11/4*x^6 - 77/64*x^4 + 55/256*x^2 - 11/1024 (epsilon ...)

    Sometimes it fails.

    ::

        sage: sin(1).minpoly()
        Traceback (most recent call last):
        ...
        ValueError: Could not find minimal polynomial (1000 bits, degree 24).

    .. note::

       Failure to produce a minimal polynomial does not
       necessarily indicate that this number is transcendental.

    AUTHORS:

    - Robert Bradshaw (2007-10): numerical algorithm

    - Robert Bradshaw (2008-10): algebraic algorithm
    """
    if algorithm is None or algorithm == 'algebraic':
        from sage.rings.all import QQbar
        try:
            return QQ[var](QQbar(ex).minpoly())
        except (TypeError, ValueError):
            if algorithm == 'algebraic':
                raise

    if algorithm is None or algorithm.startswith('numeric'):
        bits_list = [bits] if bits else [100,200,500,1000]
        degree_list = [degree] if degree else [2,4,8,12,24]

        for bits in bits_list:
            a = ex.numerical_approx(bits)
            check_bits = int(1.25 * bits + 80)
            aa = ex.numerical_approx(check_bits)

            for degree in degree_list:

                f = QQ[var](algdep(a, degree)) # TODO: use the known_bits parameter?
                # If indeed we have found a minimal polynomial,
                # it should be accurate to a much higher precision.
                error = abs(f(aa))
                dx = ~RR(Integer(1) << (check_bits - degree - 2))
                expected_error = abs(f.derivative()(CC(aa))) * dx

                if error < expected_error:
                    # Degree might have been an over-estimate, factor because we want (irreducible) minpoly.
                    ff = f.factor()
                    for g, e in ff:
                        lead = g.leading_coefficient()
                        if lead != 1:
                            g = g / lead
                        expected_error = abs(g.derivative()(CC(aa))) * dx
                        error = abs(g(aa))
                        if error < expected_error:
                            # See if we can prove equality exactly
                            if g(ex).simplify_trig().simplify_radical() == 0:
                                return g
                            # Otherwise fall back to numerical guess
                            elif epsilon and error < epsilon:
                                return g
                            else:
                                raise NotImplementedError, "Could not prove minimal polynomial %s (epsilon %s)" % (g, RR(error).str(no_sci=False))

        raise ValueError, "Could not find minimal polynomial (%s bits, degree %s)." % (bits, degree)

    raise ValueError, "Unknown algorithm: %s" % algorithm


###################################################################
# limits
###################################################################
def limit(ex, dir=None, taylor=False, algorithm='maxima', **argv):
    r"""
    Return the limit as the variable `v` approaches `a`
    from the given direction.

    ::

       expr.limit(x = a)
       expr.limit(x = a, dir='above')

    INPUT:

    -  ``dir`` - (default: None); dir may have the value
       'plus' (or 'above') for a limit from above, 'minus' (or 'below')
       for a limit from below, or may be omitted (implying a two-sided
       limit is to be computed).

    -  ``taylor`` - (default: False); if True, use Taylor
       series, which allows more integrals to be computed (but may also
       crash in some obscure cases due to bugs in Maxima).

    -  ``**argv`` - 1 named parameter

    .. note::

       The output may also use 'und' (undefined), 'ind'
       (indefinite but bounded), and 'infinity' (complex
       infinity).

    EXAMPLES::

        sage: x = var('x')
        sage: f = (1+1/x)^x
        sage: f.limit(x = oo)
        e
        sage: f.limit(x = 5)
        7776/3125
        sage: f.limit(x = 1.2)
        2.06961575467...
        sage: f.limit(x = I, taylor=True)
        (-I + 1)^I
        sage: f(x=1.2)
        2.0696157546720...
        sage: f(x=I)
        (-I + 1)^I
        sage: CDF(f(x=I))
        2.06287223508 + 0.74500706218*I
        sage: CDF(f.limit(x = I))
        2.06287223508 + 0.74500706218*I

    More examples::

        sage: limit(x*log(x), x = 0, dir='above')
        0
        sage: lim((x+1)^(1/x),x = 0)
        e
        sage: lim(e^x/x, x = oo)
        +Infinity
        sage: lim(e^x/x, x = -oo)
        0
        sage: lim(-e^x/x, x = oo)
        -Infinity
        sage: lim((cos(x))/(x^2), x = 0)
        +Infinity
        sage: lim(sqrt(x^2+1) - x, x = oo)
        0
        sage: lim(x^2/(sec(x)-1), x=0)
        2
        sage: lim(cos(x)/(cos(x)-1), x=0)
        -Infinity
        sage: lim(x*sin(1/x), x=0)
        0

    ::

        sage: f = log(log(x))/log(x)
        sage: forget(); assume(x<-2); lim(f, x=0, taylor=True)
        limit(log(log(x))/log(x), x, 0)

    Here ind means "indefinite but bounded"::

        sage: lim(sin(1/x), x = 0)
        ind
    """
    if not isinstance(ex, Expression):
        ex = SR(ex)

    if len(argv) != 1:
        raise ValueError, "call the limit function like this, e.g. limit(expr, x=2)."
    else:
        k = argv.keys()[0]
        v = var(k)
        a = argv[k]

    if taylor and algorithm == 'maxima':
        algorithm = 'maxima_taylor'

    if dir not in [None, 'plus', 'above', 'minus', 'below']:
        raise ValueError, "dir must be one of 'plus' or 'minus'"

    if algorithm == 'maxima':
        if dir is None:
            l = ex._maxima_().limit(v, a)
        elif dir == 'plus' or dir == 'above':
            l = ex._maxima_().limit(v, a, 'plus')
        elif dir == 'minus' or dir == 'below':
            l = ex._maxima_().limit(v, a, 'minus')
    elif algorithm == 'maxima_taylor':
        if dir is None:
            l = ex._maxima_().tlimit(v, a)
        elif dir == 'plus' or dir == 'above':
            l = ex._maxima_().tlimit(v, a, 'plus')
        elif dir == 'minus' or dir == 'below':
            l = ex._maxima_().tlimit(v, a, 'minus')
    elif algorithm == 'sympy':
        if dir is None:
            import sympy
            l = sympy.limit(ex._sympy_(), v._sympy_(), a._sympy_())
        else:
            raise NotImplementedError, "sympy does not support one-sided limits"

    return ex.parent()(l)

lim = limit

###################################################################
# Laplace transform
###################################################################
def laplace(ex, t, s):
    r"""
    Attempts to compute and return the Laplace transform of
    ``self`` with respect to the variable `t` and
    transform parameter `s`. If this function cannot find a
    solution, a formal function is returned.

    The function that is returned may be be viewed as a function of
    `s`.

    DEFINITION: The Laplace transform of a function `f(t)`,
    defined for all real numbers `t \geq 0`, is the function
    `F(s)` defined by

    .. math::

                      F(s) = \int_{0}^{\infty} e^{-st} f(t) dt.



    EXAMPLES: We compute a few Laplace transforms::

        sage: var('x, s, z, t, t0')
        (x, s, z, t, t0)
        sage: sin(x).laplace(x, s)
        1/(s^2 + 1)
        sage: (z + exp(x)).laplace(x, s)
        z/s + 1/(s - 1)
        sage: log(t/t0).laplace(t, s)
         -(euler_gamma + log(s) + log(t0))/s

    We do a formal calculation::

        sage: f = function('f', x)
        sage: g = f.diff(x); g
        D[0](f)(x)
        sage: g.laplace(x, s)
        s*laplace(f(x), x, s) - f(0)

    EXAMPLE: A BATTLE BETWEEN the X-women and the Y-men (by David
    Joyner): Solve

    .. math::

                   x' = -16y, x(0)=270,  y' = -x + 1, y(0) = 90.


    This models a fight between two sides, the "X-women" and the
    "Y-men", where the X-women have 270 initially and the Y-men have
    90, but the Y-men are better at fighting, because of the higher
    factor of "-16" vs "-1", and also get an occasional reinforcement,
    because of the "+1" term.

    ::

        sage: var('t')
        t
        sage: t = var('t')
        sage: x = function('x', t)
        sage: y = function('y', t)
        sage: de1 = x.diff(t) + 16*y
        sage: de2 = y.diff(t) + x - 1
        sage: de1.laplace(t, s)
        s*laplace(x(t), t, s) + 16*laplace(y(t), t, s) - x(0)
        sage: de2.laplace(t, s)
        s*laplace(y(t), t, s) - 1/s + laplace(x(t), t, s) - y(0)

    Next we form the augmented matrix of the above system::

        sage: A = matrix([[s, 16, 270],[1, s, 90+1/s]])
        sage: E = A.echelon_form()
        sage: xt = E[0,2].inverse_laplace(s,t)
        sage: yt = E[1,2].inverse_laplace(s,t)
        sage: xt
        629/2*e^(-4*t) - 91/2*e^(4*t) + 1
        sage: yt
        629/8*e^(-4*t) + 91/8*e^(4*t)
        sage: p1 = plot(xt,0,1/2,rgbcolor=(1,0,0))
        sage: p2 = plot(yt,0,1/2,rgbcolor=(0,1,0))
        sage: (p1+p2).save()

    Another example::

        sage: var('a,s,t')
        (a, s, t)
        sage: f = exp (2*t + a) * sin(t) * t; f
        t*e^(a + 2*t)*sin(t)
        sage: L = laplace(f, t, s); L
        2*(s - 2)*e^a/(s^2 - 4*s + 5)^2
        sage: inverse_laplace(L, s, t)
        t*e^(a + 2*t)*sin(t)

    Unable to compute solution::

        sage: laplace(1/s, s, t)
        laplace(1/s, s, t)

    """
    if not isinstance(ex, (Expression, SFunction)):
        ex = SR(ex)
    return ex.parent()(ex._maxima_().laplace(var(t), var(s)))

def inverse_laplace(ex, t, s):
    r"""
    Attempts to compute the inverse Laplace transform of
    ``self`` with respect to the variable `t` and
    transform parameter `s`. If this function cannot find a
    solution, a formal function is returned.

    The function that is returned may be be viewed as a function of
    `s`.

    DEFINITION: The inverse Laplace transform of a function
    `F(s)`, is the function `f(t)` defined by

    .. math::

                      F(s) = \frac{1}{2\pi i} \int_{\gamma-i\infty}^{\gamma + i\infty} e^{st} F(s) dt,


    where `\gamma` is chosen so that the contour path of
    integration is in the region of convergence of `F(s)`.

    EXAMPLES::

        sage: var('w, m')
        (w, m)
        sage: f = (1/(w^2+10)).inverse_laplace(w, m); f
        1/10*sqrt(10)*sin(sqrt(10)*m)
        sage: laplace(f, m, w)
        1/(w^2 + 10)

        sage: f(t) = t*cos(t)
        sage: s = var('s')
        sage: L = laplace(f, t, s); L
        t |--> 2*s^2/(s^2 + 1)^2 - 1/(s^2 + 1)
        sage: inverse_laplace(L, s, t)
        t |--> t*cos(t)
        sage: inverse_laplace(1/(s^3+1), s, t)
        1/3*(sqrt(3)*sin(1/2*sqrt(3)*t) - cos(1/2*sqrt(3)*t))*e^(1/2*t) + 1/3*e^(-t)

    No explicit inverse Laplace transform, so one is returned formally
    as a function ``ilt``::

        sage: inverse_laplace(cos(s), s, t)
        ilt(cos(s), s, t)

    """
    if not isinstance(ex, Expression):
        ex = SR(ex)
    return ex.parent()(ex._maxima_().ilt(var(t), var(s)))










#############################################3333
def var_cmp(x,y):
    """
    Return comparison of the two variables x and y, which is just the
    comparison of the underlying string representations of the
    variables. This is used internally by the Calculus package.

    INPUT:

    -  ``x, y`` - symbolic variables

    OUTPUT: Python integer; either -1, 0, or 1.

    EXAMPLES::

        sage: sage.calculus.calculus.var_cmp(x,x)
        0
        sage: sage.calculus.calculus.var_cmp(x,var('z'))
        -1
        sage: sage.calculus.calculus.var_cmp(x,var('a'))
        1
    """
    return cmp(repr(x), repr(y))

def function(s, *args, **kwds):
    """
    Create a formal symbolic function with the name *s*.

    EXAMPLES::

        sage: var('a, b')
        (a, b)
        sage: f = function('cr', a)
        sage: g = f.diff(a).integral(b)
        sage: g
        b*D[0](cr)(a)

    In Sage 4.0, you need to use :meth:`substitute_function` to
    replace all occurrences of a function with another::

        sage: g.substitute_function(cr, cos)
        -b*sin(a)

        sage: g.substitute_function(cr, (sin(x) + cos(x)).function(x))
        -(sin(a) - cos(a))*b

    In Sage 4.0, basic arithmetic with unevaluated functions is no
    longer supported::

        sage: x = var('x')
        sage: f = function('f')
        sage: 2*f
        Traceback (most recent call last):
        ...
        TypeError: unsupported operand parent(s) for '*': 'Integer Ring' and '<type 'sage.symbolic.function.SFunction'>'

    You now need to evaluate the function in order to do the arithmetic::

        sage: 2*f(x)
        2*f(x)
    """
    if isinstance(s, SFunction):
        return s

    if len(args) > 0:
        return function(s, **kwds)(*args)

    s = str(s)
    if ',' in s:
        return tuple([function(x.strip(), **kwds) for x in s.split(',')])
    elif ' ' in s:
        return tuple([function(x.strip(), **kwds) for x in s.split()])
    try:
        return symbol_table['functions'][s]
    except KeyError:
        pass
    f = SFunction(s, **kwds)
    symbol_table['functions'][s] = f
    return f

def clear_functions():
    """
    Clear all user-defined functions from the symbol tables.

    EXAMPLES::

        sage: f = function('foo')
        sage: id_f1 = id(f)
        sage: id_f1 == id(function('foo'))
        True
        sage: clear_functions()
        sage: id_f1 == id(function('foo'))
        False
    """
    functions = symbol_table['functions']
    for name, f in functions.items():
        if name in ['limit']:
            continue
        if not isinstance(f, PrimitiveFunction):
            del functions[name]

def dummy_limit(*args):
    """
    This function is called to create formal wrappers of limits that
    Maxima can't compute:

    EXAMPLES::

        sage: a = lim(exp(x^2)*(1-erf(x)), x=infinity); a
        limit(-e^(x^2)*erf(x) + e^(x^2), x, +Infinity)
        sage: a = sage.calculus.calculus.dummy_limit(sin(x)/x, x, 0);a
        limit(sin(x)/x, x, 0)
    """
    return _limit(args[0], var(repr(args[1])), SR(args[2]))

def dummy_diff(*args):
    """
    This function is called when 'diff' appears in a Maxima string.

    EXAMPLES::

        sage: from sage.calculus.calculus import dummy_diff
        sage: x,y = var('x,y')
        sage: dummy_diff(sin(x*y), x, SR(2), y, SR(1))
        -x*y^2*cos(x*y) - 2*y*sin(x*y)

    Here the function is used implicitly::

        sage: a = var('a')
        sage: f = function('cr', a)
        sage: g = f.diff(a); g
        D[0](cr)(a)
    """
    f = args[0]
    args = list(args[1:])
    for i in range(1, len(args), 2):
        args[i] = Integer(args[i])
    return f.diff(*args)

def dummy_integrate(*args):
    """
    This function is called to create formal wrappers of integrals that
    Maxima can't compute:

    EXAMPLES::

        sage: from sage.calculus.calculus import dummy_integrate
        sage: f(x) = function('f',x)
        sage: dummy_integrate(f(x), x)
        integrate(f(x), x)
        sage: a,b = var('a,b')
        sage: dummy_integrate(f(x), x, a, b)
        integrate(f(x), x, a, b)

    """
    if len(args) == 4:
        return _integrate(args[0], var(repr(args[1])), SR(args[2]), SR(args[3]))
    else:
        return _integrate(args[0], var(repr(args[1])))

def dummy_laplace(*args):
    """
    This function is called to create formal wrappers of laplace transforms
    that Maxima can't compute:

    EXAMPLES::

        sage: from sage.calculus.calculus import dummy_laplace
        sage: s,t = var('s,t')
        sage: f(t) = function('f',t)
        sage: dummy_laplace(f(t),t,s)
        laplace(f(t), t, s)
    """
    return _laplace(args[0], var(repr(args[1])), var(repr(args[2])))

def dummy_inverse_laplace(*args):
    """
    This function is called to create formal wrappers of inverse laplace
    transforms that Maxima can't compute:

    EXAMPLES::

        sage: from sage.calculus.calculus import dummy_inverse_laplace
        sage: s,t = var('s,t')
        sage: F(s) = function('F',s)
        sage: dummy_inverse_laplace(F(s),s,t)
        ilt(F(s), s, t)
    """
    return _inverse_laplace(args[0], var(repr(args[1])), var(repr(args[2])))

#######################################################
#
# Helper functions for printing latex expression
#
#######################################################

def _limit_latex_(*args):
    r"""
    Return latex expression for limit of a symbolic function.

    EXAMPLES::

        sage: from sage.calculus.calculus import _limit_latex_
        sage: var('x,a')
        (x, a)
        sage: f(x) = function('f',x)
        sage: _limit_latex_(f(x), x, a)
        '\\lim_{x \\to a}\\, f\\left(x\\right)'

    AUTHORS:

    - Golam Mortuza Hossain (2009-06-15)
    """
    # Read f,x,a from arguments
    f = args[0]
    x = args[1]
    a = args[2]
    return "\\lim_{%s \\to %s}\\, %s"%(latex(x), latex(a), latex(f))

def _integrate_latex_(*args):
    r"""
    Return LaTeX expression for integration of a symbolic function.

    EXAMPLES::

        sage: from sage.calculus.calculus import _integrate_latex_
        sage: var('x,a,b')
        (x, a, b)
        sage: f(x) = function('f',x)
        sage: _integrate_latex_(f(x),x)
        '\\int f\\left(x\\right)\\,{d x}'
        sage: _integrate_latex_(f(x),x,a,b)
        '\\int_{a}^{b} f\\left(x\\right)\\,{d x}'

    AUTHORS:

    - Golam Mortuza Hossain (2009-06-22)
    """
    f = args[0]
    x = args[1]
    # Check whether its a definite integral
    if len(args) == 4:
        a = args[2]
        b = args[3]
        return "\\int_{%s}^{%s} %s\\,{d %s}"%(latex(a), latex(b), latex(f), latex(x))
    # Typeset as indefinite integral
    return "\\int %s\\,{d %s}"%(latex(f), latex(x))

def _laplace_latex_(*args):
    r"""
    Return LaTeX expression for Laplace transform of a symbolic function.

    EXAMPLES::

        sage: from sage.calculus.calculus import _laplace_latex_
        sage: var('s,t')
        (s, t)
        sage: f(t) = function('f',t)
        sage: _laplace_latex_(f(t),t,s)
        '\\mathcal{L}\\left(f\\left(t\\right), t, s\\right)'

    AUTHORS:

    - Golam Mortuza Hossain (2009-06-22)
    """
    return "\\mathcal{L}\\left(%s\\right)"%(', '.join([latex(x) for x in args]))

def _inverse_laplace_latex_(*args):
    r"""
    Return LaTeX expression for inverse Laplace transform of a symbolic function.

    EXAMPLES::

        sage: from sage.calculus.calculus import _inverse_laplace_latex_
        sage: var('s,t')
        (s, t)
        sage: F(s) = function('F',s)
        sage: _inverse_laplace_latex_(F(s),s,t)
        '\\mathcal{L}^{-1}\\left(F\\left(s\\right), s, t\\right)'

    AUTHORS:

    - Golam Mortuza Hossain (2009-06-22)
    """
    return "\\mathcal{L}^{-1}\\left(%s\\right)"%(', '.join([latex(x) for x in args]))

# Return un-evaluated expression as instances of SFunction class
_limit = SFunction('limit', print_latex_func=_limit_latex_)
_integrate = SFunction('integrate', print_latex_func=_integrate_latex_)
_laplace = SFunction('laplace', print_latex_func=_laplace_latex_)
_inverse_laplace = SFunction('ilt', print_latex_func=_inverse_laplace_latex_)

######################################i################




#######################################################

symtable = {'%pi':'pi', '%e': 'e', '%i':'I', '%gamma':'euler_gamma',
            'li[2]':'polylog2', 'li[3]':'polylog3'}

from sage.symbolic.pynac import register_symbol
from sage.rings.infinity import infinity, minus_infinity
register_symbol(infinity, dict(maxima='inf'))
register_symbol(minus_infinity, dict(maxima='minf'))

from sage.misc.multireplace import multiple_replace

import re


maxima_tick = re.compile("'[a-z|A-Z|0-9|_]*")

maxima_qp = re.compile("\?\%[a-z|A-Z|0-9|_]*")  # e.g., ?%jacobi_cd

maxima_var = re.compile("\%[a-z|A-Z|0-9|_]*")  # e.g., ?%jacobi_cd

sci_not = re.compile("(-?(?:0|[1-9]\d*))(\.\d+)?([eE][-+]\d+)")

def symbolic_expression_from_maxima_string(x, equals_sub=False, maxima=maxima):
    """
    Given a string representation of a Maxima expression, parse it and
    return the corresponding Sage symbolic expression.

    INPUT:


    -  ``x`` - a string

    -  ``equals_sub`` - (default: False) if True, replace
       '=' by '==' in self

    -  ``maxima`` - (default: the calculus package's
       Maxima) the Maxima interpreter to use.


    EXAMPLES::

        sage: sage.calculus.calculus.symbolic_expression_from_maxima_string('x^%e + %e^%pi + %i + sin(0)')
        x^e + e^pi + I
    """
    syms = dict(_syms)

    if len(x) == 0:
        raise RuntimeError, "invalid symbolic expression -- ''"
    maxima.set('_tmp_',x)

    # This is inefficient since it so rarely is needed:
    #r = maxima._eval_line('listofvars(_tmp_);')[1:-1]

    s = maxima._eval_line('_tmp_;')

    formal_functions = maxima_tick.findall(s)
    if len(formal_functions) > 0:
        for X in formal_functions:
            syms[X[1:]] = function(X[1:])
        # You might think there is a potential very subtle bug if 'foo
        # is in a string literal -- but string literals should *never*
        # ever be part of a symbolic expression.
        s = s.replace("'","")

    delayed_functions = maxima_qp.findall(s)
    if len(delayed_functions) > 0:
        for X in delayed_functions:
            syms[X[2:]] = SFunction(X[2:])
        s = s.replace("?%","")

    s = multiple_replace(symtable, s)
    s = s.replace("%","")

    if equals_sub:
        s = s.replace('=','==')

    #replace all instances of Maxima's scientific notation
    #with regular notation
    search = sci_not.search(s)
    while not search is None:
        (start, end) = search.span()
        r = create_RealNumber(s[start:end]).str(no_sci=2, truncate=True)
        s = s.replace(s[start:end], r)
        search = sci_not.search(s)

    # have to do this here, otherwise maxima_tick catches it
    syms['limit'] = dummy_limit
    syms['diff'] = dummy_diff
    syms['integrate'] = dummy_integrate
    syms['laplace'] = dummy_laplace
    syms['ilt'] = dummy_inverse_laplace

    global is_simplified
    try:
        # use a global flag so all expressions obtained via
        # evaluation of maxima code are assumed pre-simplified
        is_simplified = True
        return symbolic_expression_from_string(s, syms, accept_sequence=True)
    except SyntaxError:
        raise TypeError, "unable to make sense of Maxima expression '%s' in Sage"%s
    finally:
        is_simplified = False


# External access used by restore
from sage.symbolic.pynac import symbol_table
_syms = syms_cur = symbol_table.get('maxima', {})
syms_default = dict(syms_cur)

# Comma format options for Maxima
def mapped_opts(v):
    """
    Used internally when creating a string of options to pass to
    Maxima.

    INPUT:


    -  ``v`` - an object


    OUTPUT: a string.

    The main use of this is to turn Python bools into lower case
    strings.

    EXAMPLES::

        sage: sage.calculus.calculus.mapped_opts(True)
        'true'
        sage: sage.calculus.calculus.mapped_opts(False)
        'false'
        sage: sage.calculus.calculus.mapped_opts('bar')
        'bar'
    """
    if isinstance(v, bool):
        return str(v).lower()
    return str(v)

def maxima_options(**kwds):
    """
    Used internally to create a string of options to pass to Maxima.

    EXAMPLES::

        sage: sage.calculus.calculus.maxima_options(an_option=True, another=False, foo='bar')
        'an_option=true,foo=bar,another=false'
    """
    return ','.join(['%s=%s'%(key,mapped_opts(val)) for key, val in kwds.iteritems()])


# Parser for symbolic ring elements

_augmented_syms = {}

def _find_var(name):
    try:
        return (_augmented_syms or _syms)[name]
    except KeyError:
        pass
    try:
        return SR(sage.all.__dict__[name])
    except (KeyError, TypeError):
        return var(name)

def _find_func(name):
    try:
        func = (_augmented_syms or _syms)[name]
        if not isinstance(func, Expression):
            return func
    except KeyError:
        pass
    try:
        func = SR(sage.all.__dict__[name])
        if not isinstance(func, Expression):
            return func
    except (KeyError, TypeError):
        return function(name)

SR_parser = Parser(make_int      = lambda x: SR(Integer(x)),
                   make_float    = lambda x: SR(RealDoubleElement(x)),
                   make_var      = _find_var,
                   make_function = _find_func)

def symbolic_expression_from_string(s, syms=None, accept_sequence=False):
    # from sage.functions.constants import I # can't import this at the top, but need it now
    # _syms['i'] = _syms['I'] = I
    parse_func = SR_parser.parse_sequence if accept_sequence else SR_parser.parse_expression
    if syms is None:
        return parse_func(s)
    else:
        try:
            global _augmented_syms
            _augmented_syms = syms
            return parse_func(s)
        finally:
            _augmented_syms = {}

def symbolic_expression_from_maxima_element(x, maxima=maxima):
    """
    Given an element of the calculus copy of the Maxima interface,
    create the corresponding Sage symbolic expression.

    EXAMPLES::

        sage: a = sage.calculus.calculus.maxima('x^(sqrt(y)+%pi) + sin(%e + %pi)')
        sage: sage.calculus.calculus.symbolic_expression_from_maxima_element(a)
        x^(pi + sqrt(y)) - sin(e)
        sage: var('x, y')
        (x, y)
        sage: v = sage.calculus.calculus.maxima.vandermonde_matrix([x, y, 1/2])
        sage: sage.calculus.calculus.symbolic_expression_from_maxima_element(v)
        [  1   x x^2]
        [  1   y y^2]
        [  1 1/2 1/4]
    """
    return symbolic_expression_from_maxima_string(x.name(), maxima=maxima)
