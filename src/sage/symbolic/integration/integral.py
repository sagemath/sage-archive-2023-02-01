"""
Symbolic Integration
"""

##############################################################################
#
#       Copyright (C) 2009 Golam Mortuza Hossain <gmhossain@gmail.com>
#       Copyright (C) 2010 Burcin Erocal <burcin@erocal.org>
#
#  Distributed under the terms of the GNU General Public License (GPL v2+)
#                  http://www.gnu.org/licenses/
#
##############################################################################

from sage.symbolic.ring import SR, is_SymbolicVariable
from sage.symbolic.function import BuiltinFunction, Function

##################################################################
#  Table of available integration routines
##################################################################

# Add new integration routines to the dictionary below. This will make them
# accessible with the 'algorithm' keyword parameter of top level integrate().
available_integrators = {}

import sage.symbolic.integration.external as external
available_integrators['maxima'] = external.maxima_integrator
available_integrators['sympy'] = external.sympy_integrator
available_integrators['mathematica_free'] = external.mma_free_integrator

######################################################
#
# Class implementing symbolic integration
#
######################################################

class IndefiniteIntegral(BuiltinFunction):
    def __init__(self):
        """
        Class to represent an indefinite integral.

        EXAMPLES::

            sage: from sage.symbolic.integration.integral import indefinite_integral
            sage: indefinite_integral(log(x), x) #indirect doctest
            x*log(x) - x
            sage: indefinite_integral(x^2, x)
            1/3*x^3
            sage: indefinite_integral(4*x*log(x), x)
            2*x^2*log(x) - x^2
            sage: indefinite_integral(exp(x), 2*x)
            2*e^x

        """
        # The automatic evaluation routine will try these integrators
        # in the given order. This is an attribute of the class instead of
        # a global variable in this module to enable customization by
        # creating a subclasses which define a different set of integrators
        self.integrators = [external.maxima_integrator]

        BuiltinFunction.__init__(self, "integrate", nargs=2, conversions={'sympy': 'Integral'})

    def _eval_(self, f, x):
        """
        EXAMPLES::

            sage: from sage.symbolic.integration.integral import indefinite_integral
            sage: indefinite_integral(exp(x), x) # indirect doctest
            e^x
            sage: indefinite_integral(exp(x), x^2)
            2*(x - 1)*e^x
        """
        # Check for x
        if not is_SymbolicVariable(x):
            if len(x.variables()) == 1:
                nx = x.variables()[0]
                f = f*x.diff(nx)
                x = nx
            else:
                return None

        # we try all listed integration algorithms
        for integrator in self.integrators:
            res = integrator(f, x)
            try:
                return integrator(f, x)
            except NotImplementedError:
                pass
        return None

    def _tderivative_(self, f, x, diff_param=None):
        """
        EXAMPLES::

            sage: from sage.symbolic.integration.integral import indefinite_integral
            sage: f = function('f'); a,b=var('a,b')
            sage: h = indefinite_integral(f(x), x)
            sage: h.diff(x) # indirect doctest
            f(x)
            sage: h.diff(a)
            0
        """
        if x.has(diff_param):
            return f*x.derivative(diff_param)
        else:
            return f.derivative(diff_param).integral(x)

    def _print_latex_(self, f, x):
        """
        EXAMPLES::

            sage: from sage.symbolic.integration.integral import indefinite_integral
            sage: print_latex = indefinite_integral._print_latex_
            sage: var('x,a,b')
            (x, a, b)
            sage: f = function('f')
            sage: print_latex(f(x),x)
            '\\int f\\left(x\\right)\\,{d x}'
        """
        from sage.misc.latex import latex
        if not is_SymbolicVariable(x):
            dx_str = "{d \\left(%s\\right)}"%(latex(x))
        else:
            dx_str = "{d %s}"%(latex(x))

        return "\\int %s\\,%s"%(latex(f), dx_str)

indefinite_integral = IndefiniteIntegral()

class DefiniteIntegral(BuiltinFunction):
    def __init__(self):
        """
        Symbolic function representing a definite integral.

        EXAMPLES::

            sage: from sage.symbolic.integration.integral import definite_integral
            sage: definite_integral(sin(x),x,0,pi)
            2
        """
        # The automatic evaluation routine will try these integrators
        # in the given order. This is an attribute of the class instead of
        # a global variable in this module to enable customization by
        # creating a subclasses which define a different set of integrators
        self.integrators = [external.maxima_integrator]

        BuiltinFunction.__init__(self, "integrate", nargs=4, conversions={'sympy': 'Integral'})

    def _eval_(self, f, x, a, b):
        """
        Returns the results of symbolic evaluation of the integral

        EXAMPLES::

            sage: from sage.symbolic.integration.integral import definite_integral
            sage: definite_integral(exp(x),x,0,1) # indirect doctest
            e - 1
        """
        # Check for x
        if not is_SymbolicVariable(x):
            if len(x.variables()) == 1:
                nx = x.variables()[0]
                f = f*x.diff(nx)
                x = nx
            else:
                return None

        args = (f,x,a,b)

        # we try all listed integration algorithms
        for integrator in self.integrators:
            try:
                return integrator(*args)
            except NotImplementedError:
                pass
        return None

    def _evalf_(self, f, x, a, b, parent=None, algorithm=None):
        """
        Returns numerical approximation of the integral

        EXAMPLES::

            sage: from sage.symbolic.integration.integral import definite_integral
            sage: h = definite_integral(sin(x)*log(x)/x^2, x, 1, 2); h
            integrate(log(x)*sin(x)/x^2, x, 1, 2)
            sage: h.n() # indirect doctest
            0.14839875208053...

        TESTS:

        Check if #3863 is fixed::

            sage: integrate(x^2.7 * e^(-2.4*x), x, 0, 3).n()
            0.154572952320790
        """
        from sage.gsl.integration import numerical_integral
        # The gsl routine returns a tuple, which also contains the error.
        # We only return the result.
        return numerical_integral(f, a, b)[0]

    def _tderivative_(self, f, x, a, b, diff_param=None):
        """
        Returns derivative of symbolic integration

        EXAMPLES::

            sage: from sage.symbolic.integration.integral import definite_integral
            sage: f = function('f'); a,b=var('a,b')
            sage: h = definite_integral(f(x), x,a,b)
            sage: h.diff(x) # indirect doctest
            0
            sage: h.diff(a)
            -f(a)
            sage: h.diff(b)
            f(b)
        """
        if not x.has(diff_param):
            # integration variable != differentiation variable
            ans = definite_integral(f.diff(diff_param), x, a, b)
        else:
            ans = SR(0)
        return ans + f.subs(x==b)*b.diff(diff_param) \
                    - f.subs(x==a)*a.diff(diff_param)

    def _print_latex_(self, f, x, a, b):
        r"""
        Returns LaTeX expression for integration of a symbolic function.

        EXAMPLES::

            sage: from sage.symbolic.integration.integral import definite_integral
            sage: print_latex = definite_integral._print_latex_
            sage: var('x,a,b')
            (x, a, b)
            sage: f = function('f')
            sage: print_latex(f(x),x,0,1)
            '\\int_{0}^{1} f\\left(x\\right)\\,{d x}'
            sage: latex(integrate(1/(1+sqrt(x)),x,0,1))
            \int_{0}^{1} \frac{1}{\sqrt{x} + 1}\,{d x}
        """
        from sage.misc.latex import latex
        if not is_SymbolicVariable(x):
            dx_str = "{d \\left(%s\\right)}"%(latex(x))
        else:
            dx_str = "{d %s}"%(latex(x))
        return "\\int_{%s}^{%s} %s\\,%s"%(latex(a), latex(b), latex(f), dx_str)

definite_integral = DefiniteIntegral()


def _normalize_integral_input(f, v=None, a=None, b=None):
    r"""
    Validate and return variable and endpoints for an integral.

    INPUT:

    - ``f`` -- an expression to integrate;

    - ``v`` -- a variable of integration or a triple;

    - ``a`` -- (optional) the left endpoint of integration;

    - ``b`` -- (optional) the right endpoint of integration.

    It is also possible to pass last three parameters in ``v`` as a triple.

    OUPUT:

    - a tuple of ``f``, ``v``, ``a``, and ``b``.

    EXAMPLES::

        sage: from sage.symbolic.integration.integral import \
        ...       _normalize_integral_input
        sage: _normalize_integral_input(x^2, x, 0, 3)
        (x^2, x, 0, 3)
        sage: _normalize_integral_input(x^2, [x, 0, 3], None, None)
        (x^2, x, 0, 3)
        sage: _normalize_integral_input(x^2, [0, 3], None, None)
        doctest:...: DeprecationWarning:
        Variable of integration should be specified explicitly.
        See http://trac.sagemath.org/12438 for details.
        (x^2, x, 0, 3)
        sage: _normalize_integral_input(x^2, [x], None, None)
        (x^2, x, None, None)
    """
    if isinstance(v, (list, tuple)) and a is None and b is None:
        if len(v) == 1: # bare variable in a tuple
            v = v[0]
        elif len(v) == 2: # endpoints only
            a, b = v
            v = None
        elif len(v) == 3: # variable and endpoints
            v, a, b = v
        else:
            raise ValueError("invalid input %s - please use variable, "
                             "with or without two endpoints" % repr(v))
    elif b is None and a is not None:
        # two arguments, must be endpoints
        v, a, b = None, v, a
    if v is None:
        from sage.misc.superseded import deprecation
        deprecation(12438, "Variable of integration should be specified explicitly.")
        v = f.default_variable()
        if isinstance(f, Function):  # a bare function like sin
            f = f(v)
    if (a is None) ^ (b is None):
        raise TypeError('only one endpoint was given!')
    return f, v, a, b

def integrate(expression, v=None, a=None, b=None, algorithm=None):
    r"""
    Returns the indefinite integral with respect to the variable
    `v`, ignoring the constant of integration. Or, if endpoints
    `a` and `b` are specified, returns the definite
    integral over the interval `[a, b]`.

    If ``self`` has only one variable, then it returns the
    integral with respect to that variable.

    If definite integration fails, it could be still possible to
    evaluate the definite integral using indefinite integration with
    the Newton - Leibniz theorem (however, the user has to ensure that the
    indefinite integral is continuous on the compact interval `[a,b]` and
    this theorem can be applied).

    INPUT:

    - ``v`` - a variable or variable name.  This can also be a tuple of
      the variable (optional) and endpoints (i.e., ``(x,0,1)`` or ``(0,1)``).

    - ``a`` - (optional) lower endpoint of definite integral

    - ``b`` - (optional) upper endpoint of definite integral

    - ``algorithm`` - (default: 'maxima') one of

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
        sage: f.integral(x)
        1/2*(4*x + 3)/(x^2 + 2*x + 1) + log(x + 1)

    ::

        sage: f = x*cos(x^2)
        sage: f.integral(x, 0, sqrt(pi))
        0
        sage: f.integral(x, a=-pi, b=pi)
        0

    ::

        sage: f(x) = sin(x)
        sage: f.integral(x, 0, pi/2)
        1

    The variable is required, but the endpoints are optional::

        sage: y=var('y')
        sage: integral(sin(x), x)
        -cos(x)
        sage: integral(sin(x), y)
        y*sin(x)
        sage: integral(sin(x), x, pi, 2*pi)
        -2
        sage: integral(sin(x), y, pi, 2*pi)
        pi*sin(x)
        sage: integral(sin(x), (x, pi, 2*pi))
        -2
        sage: integral(sin(x), (y, pi, 2*pi))
        pi*sin(x)

    Constraints are sometimes needed::

        sage: var('x, n')
        (x, n)
        sage: integral(x^n,x)
        Traceback (most recent call last):
        ...
        ValueError: Computation failed since Maxima requested additional
        constraints; using the 'assume' command before integral evaluation
        *may* help (example of legal syntax is 'assume(n>0)', see `assume?`
        for more details)
        Is n equal to -1?
        sage: assume(n > 0)
        sage: integral(x^n,x)
        x^(n + 1)/(n + 1)
        sage: forget()

    Usually the constraints are of sign, but others are possible::

        sage: assume(n==-1)
        sage: integral(x^n,x)
        log(x)

    Note that an exception is raised when a definite integral is
    divergent::

        sage: forget() # always remember to forget assumptions you no longer need
        sage: integrate(1/x^3,(x,0,1))
        Traceback (most recent call last):
        ...
        ValueError: Integral is divergent.
        sage: integrate(1/x^3,x,-1,3)
        Traceback (most recent call last):
        ...
        ValueError: Integral is divergent.

    But Sage can calculate the convergent improper integral of
    this function::

        sage: integrate(1/x^3,x,1,infinity)
        1/2

    The examples in the Maxima documentation::

        sage: var('x, y, z, b')
        (x, y, z, b)
        sage: integral(sin(x)^3, x)
        1/3*cos(x)^3 - cos(x)
        sage: integral(x/sqrt(b^2-x^2), b)
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
        sage: g = mathematica(f)                           # optional - mathematica
        sage: print g                                      # optional - mathematica
                  z        2
                 y  + Sin[x ]
        sage: print g.Integrate(x)                         # optional - mathematica
                    z        Pi                2
                 x y  + Sqrt[--] FresnelS[Sqrt[--] x]
                             2                 Pi
        sage: print f.integral(x)
        x*y^z + 1/8*sqrt(pi)*((I + 1)*sqrt(2)*erf((1/2*I + 1/2)*sqrt(2)*x) + (I - 1)*sqrt(2)*erf((1/2*I - 1/2)*sqrt(2)*x))

    Alternatively, just use algorithm='mathematica_free' to integrate via Mathematica
    over the internet (does NOT require a Mathematica license!)::

        sage: _ = var('x, y, z')   # optional - internet
        sage: f = sin(x^2) + y^z   # optional - internet
        sage: f.integrate(x, algorithm="mathematica_free")   # optional - internet
        x*y^z + sqrt(1/2)*sqrt(pi)*fresnels(sqrt(2)*x/sqrt(pi))

    We can also use Sympy::

        sage: integrate(x*sin(log(x)), x)
        -1/5*x^2*(cos(log(x)) - 2*sin(log(x)))
        sage: integrate(x*sin(log(x)), x, algorithm='sympy')
        -1/5*x^2*cos(log(x)) + 2/5*x^2*sin(log(x))
        sage: _ = var('y, z')
        sage: (x^y - z).integrate(y)
        -y*z + x^y/log(x)
        sage: (x^y - z).integrate(y, algorithm="sympy")  # see Trac #14694
        Traceback (most recent call last):
        ...
        AttributeError: 'Piecewise' object has no attribute '_sage_'

    We integrate the above function in Maple now::

        sage: g = maple(f); g                             # optional - maple
        sin(x^2)+y^z
        sage: g.integrate(x)                              # optional - maple
        1/2*2^(1/2)*Pi^(1/2)*FresnelS(2^(1/2)/Pi^(1/2)*x)+y^z*x

    We next integrate a function with no closed form integral. Notice
    that the answer comes back as an expression that contains an
    integral itself.

    ::

        sage: A = integral(1/ ((x-4) * (x^3+2*x+1)), x); A
        -1/73*integrate((x^2 + 4*x + 18)/(x^3 + 2*x + 1), x) + 1/73*log(x - 4)

    We now show that floats are not converted to rationals
    automatically since we by default have keepfloat: true in maxima.

    ::

        sage: integral(e^(-x^2),(x, 0, 0.1))
        0.0562314580091*sqrt(pi)

    ALIASES: integral() and integrate() are the same.

    EXAMPLES:

    Here is an example where we have to use assume::

        sage: a,b = var('a,b')
        sage: integrate(1/(x^3 *(a+b*x)^(1/3)), x)
        Traceback (most recent call last):
        ...
        ValueError: Computation failed since Maxima requested additional
        constraints; using the 'assume' command before integral evaluation
        *may* help (example of legal syntax is 'assume(a>0)', see `assume?`
        for more details)
        Is a positive or negative?

    So we just assume that `a>0` and the integral works::

        sage: assume(a>0)
        sage: integrate(1/(x^3 *(a+b*x)^(1/3)), x)
        2/9*sqrt(3)*b^2*arctan(1/3*sqrt(3)*(2*(b*x + a)^(1/3) + a^(1/3))/a^(1/3))/a^(7/3) - 1/9*b^2*log((b*x + a)^(2/3) + (b*x + a)^(1/3)*a^(1/3) + a^(2/3))/a^(7/3) + 2/9*b^2*log((b*x + a)^(1/3) - a^(1/3))/a^(7/3) + 1/6*(4*(b*x + a)^(5/3)*b^2 - 7*(b*x + a)^(2/3)*a*b^2)/((b*x + a)^2*a^2 - 2*(b*x + a)*a^3 + a^4)

    TESTS:

    The following integral was broken prior to Maxima 5.15.0 -
    see :trac:`3013`::

        sage: integrate(sin(x)*cos(10*x)*log(x), x)
        -1/198*(9*cos(11*x) - 11*cos(9*x))*log(x) + 1/44*Ei(11*I*x) - 1/36*Ei(9*I*x) - 1/36*Ei(-9*I*x) + 1/44*Ei(-11*I*x)

    It is no longer possible to use certain functions without an
    explicit variable.  Instead, evaluate the function at a variable,
    and then take the integral::

        sage: integrate(sin)
        Traceback (most recent call last):
        ...
        TypeError

        sage: integrate(sin(x), x)
        -cos(x)
        sage: integrate(sin(x), x, 0, 1)
        -cos(1) + 1

    Check if :trac:`780` is fixed::

        sage: _ = var('x,y')
        sage: f = log(x^2+y^2)
        sage: res = integral(f,x,0.0001414, 1.); res
        Traceback (most recent call last):
        ...
        ValueError: Computation failed since Maxima requested additional constraints; using the 'assume' command before integral evaluation *may* help (example of legal syntax is 'assume(50015104*y^2-50015103>0)', see `assume?` for more details)
        Is 50015104*y^2-50015103 positive, negative or zero?
        sage: assume(y>1)
        sage: res = integral(f,x,0.0001414, 1.); res
        -2*y*arctan(0.0001414/y) + 2*y*arctan(1/y) + log(y^2 + 1.0) - 0.0001414*log(y^2 + 1.999396e-08) - 1.9997172
        sage: nres = numerical_integral(f.subs(y=2), 0.0001414, 1.); nres
        (1.4638323264144..., 1.6251803529759...e-14)
        sage: res.subs(y=2).n()
        1.46383232641443
        sage: nres = numerical_integral(f.subs(y=.5), 0.0001414, 1.); nres
        (-0.669511708872807, 7.768678110854711e-15)
        sage: res.subs(y=.5).n()
        -0.669511708872807

    Check if :trac:`6816` is fixed::

        sage: var('t,theta')
        (t, theta)
        sage: integrate(t*cos(-theta*t),t,0,pi)
        (pi*theta*sin(pi*theta) + cos(pi*theta))/theta^2 - 1/theta^2
        sage: integrate(t*cos(-theta*t),(t,0,pi))
        (pi*theta*sin(pi*theta) + cos(pi*theta))/theta^2 - 1/theta^2
        sage: integrate(t*cos(-theta*t),t)
        (t*theta*sin(t*theta) + cos(t*theta))/theta^2
        sage: integrate(x^2,(x)) # this worked before
        1/3*x^3
        sage: integrate(x^2,(x,)) # this didn't
        1/3*x^3
        sage: integrate(x^2,(x,1,2))
        7/3
        sage: integrate(x^2,(x,1,2,3))
        Traceback (most recent call last):
        ...
        ValueError: invalid input (x, 1, 2, 3) - please use variable, with or without two endpoints

    Note that this used to be the test, but it is
    actually divergent (though Maxima as yet does
    not say so)::

        sage: integrate(t*cos(-theta*t),(t,-oo,oo))
        integrate(t*cos(t*theta), t, -Infinity, +Infinity)

    Check if :trac:`6189` is fixed::

        sage: n = N; n
        <function numerical_approx at ...>
        sage: F(x) = 1/sqrt(2*pi*1^2)*exp(-1/(2*1^2)*(x-0)^2)
        sage: G(x) = 1/sqrt(2*pi*n(1)^2)*exp(-1/(2*n(1)^2)*(x-n(0))^2)
        sage: integrate( (F(x)-F(x))^2, x, -infinity, infinity).n()
        0.000000000000000
        sage: integrate( ((F(x)-G(x))^2).expand(), x, -infinity, infinity).n()
        -6.26376265908397e-17
        sage: integrate( (F(x)-G(x))^2, x, -infinity, infinity).n()# abstol 1e-6
        0

    This was broken before Maxima 5.20::

        sage: exp(-x*i).integral(x,0,1)
        I*e^(-I) - I

    Test deprecation warning when variable is not specified::

        sage: x.integral()
        doctest:...: DeprecationWarning:
        Variable of integration should be specified explicitly.
        See http://trac.sagemath.org/12438 for details.
        1/2*x^2

    Test that :trac:`8729` is fixed::

        sage: t = var('t')
        sage: a = sqrt((sin(t))^2 + (cos(t))^2)
        sage: integrate(a, t, 0, 2*pi)
        2*pi
        sage: a.simplify_full().simplify_trig()
        1

    Maxima uses Cauchy Principal Value calculations to
    integrate certain convergent integrals.  Here we test
    that this does not raise an error message (see :trac:`11987`)::

        sage: integrate(sin(x)*sin(x/3)/x^2, x, 0, oo)
        1/6*pi

    Maxima returned a negative value for this integral prior to
    maxima-5.24 (:trac:`10923`). Ideally we would get an answer in terms
    of the gamma function; however, we get something equivalent::

        sage: actual_result = integral(e^(-1/x^2), x, 0, 1)
        sage: actual_result.simplify_radical()
        (sqrt(pi)*(erf(1)*e - e) + 1)*e^(-1)
        sage: ideal_result = 1/2*gamma(-1/2, 1)
        sage: error = actual_result - ideal_result
        sage: error.numerical_approx() # abs tol 1e-10
        0

    We will not get an evaluated answer here, which is better than
    the previous (wrong) answer of zero. See :trac:`10914`::

        sage: f = abs(sin(x))
        sage: integrate(f, x, 0, 2*pi)  # long time (4s on sage.math, 2012)
        integrate(abs(sin(x)), x, 0, 2*pi)

    Another incorrect integral fixed upstream in Maxima, from
    :trac:`11233`::

        sage: a,t = var('a,t')
        sage: assume(a>0)
        sage: assume(x>0)
        sage: f = log(1 + a/(x * t)^2)
        sage: F = integrate(f, t, 1, Infinity)
        sage: F(x=1, a=7).numerical_approx() # abs tol 1e-10
        4.32025625668262
        sage: forget()

    Verify that MinusInfinity works with sympy (:trac:`12345`)::

        sage: integral(1/x^2, x, -infinity, -1, algorithm='sympy')
        1

    Check that :trac:`11737` is fixed::

        sage: N(integrate(sin(x^2)/(x^2), x, 1, infinity))
        0.285736646322858

    Check that :trac:`14209` is fixed::

        sage: integral(e^(-abs(x))/cosh(x),x,-infinity,infinity)
        2*log(2)
        sage: integral(e^(-abs(x))/cosh(x),x,-infinity,infinity)
        2*log(2)
    """
    expression, v, a, b = _normalize_integral_input(expression, v, a, b)
    if algorithm is not None:
        integrator = available_integrators.get(algorithm)
        if not integrator:
            raise ValueError("Unknown algorithm: %s" % algorithm)
        return integrator(expression, v, a, b)
    if a is None:
        return indefinite_integral(expression, v)
    else:
        return definite_integral(expression, v, a, b)

integral = integrate
