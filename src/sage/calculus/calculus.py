r"""
Symbolic Computation

AUTHORS:

- Bobby Moretti and William Stein (2006-2007)

- Robert Bradshaw (2007-10): minpoly(), numerical algorithm

- Robert Bradshaw (2008-10): minpoly(), algebraic algorithm

- Golam Mortuza Hossain (2009-06-15): _limit_latex()

- Golam Mortuza Hossain (2009-06-22): _laplace_latex(), _inverse_laplace_latex()

- Tom Coates (2010-06-11): fixed Trac #9217

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
scenes through Maxima::

    sage: f = sin(x)/cos(2*y)
    sage: f.derivative(y)
    2*sin(x)*sin(2*y)/cos(2*y)^2
    sage: g = f.integral(x); g
    -cos(x)/cos(2*y)

Note that these methods usually require an explicit variable name. If none
is given, Sage will try to find one for you.

::

    sage: f = sin(x); f.derivative()
    cos(x)

If the expression is a callable symbolic expression (i.e., the
variable order is specified), then Sage can calculate the matrix
derivative (i.e., the gradient, Jacobian matrix, etc.) if no variables
are specified.  In the example below, we use the second derivative
test to determine that there is a saddle point at (0,-1/2).

::

    sage: f(x,y)=x^2*y+y^2+y
    sage: f.diff() # gradient
    (x, y) |--> (2*x*y, x^2 + 2*y + 1)
    sage: solve(list(f.diff()),[x,y])
    [[x == -I, y == 0], [x == I, y == 0], [x == 0, y == (-1/2)]]
    sage: H=f.diff(2); H  # Hessian matrix
    [(x, y) |--> 2*y (x, y) |--> 2*x]
    [(x, y) |--> 2*x   (x, y) |--> 2]
    sage: H(x=0,y=-1/2)
    [-1  0]
    [ 0  2]
    sage: H(x=0,y=-1/2).eigenvalues()
    [-1, 2]

Here we calculate the Jacobian for the polar coordinate transformation::

    sage: T(r,theta)=[r*cos(theta),r*sin(theta)]
    sage: T
    (r, theta) |--> (r*cos(theta), r*sin(theta))
    sage: T.diff() # Jacobian matrix
    [   (r, theta) |--> cos(theta) (r, theta) |--> -r*sin(theta)]
    [   (r, theta) |--> sin(theta)  (r, theta) |--> r*cos(theta)]
    sage: diff(T) # Jacobian matrix
    [   (r, theta) |--> cos(theta) (r, theta) |--> -r*sin(theta)]
    [   (r, theta) |--> sin(theta)  (r, theta) |--> r*cos(theta)]
    sage: T.diff().det() # Jacobian
    (r, theta) |--> r*cos(theta)^2 + r*sin(theta)^2

When the order of variables is ambiguous, Sage will raise an
exception when differentiating::

    sage: f = sin(x+y); f.derivative()
    Traceback (most recent call last):
    ...
    ValueError: No differentiation variable specified.

Simplifying symbolic sums is also possible, using the
sum command, which also uses Maxima in the background::

    sage: k, m = var('k, m')
    sage: sum(1/k^4, k, 1, oo)
    1/90*pi^4
    sage: sum(binomial(m,k), k, 0, m)
    2^m

Symbolic matrices can be used as well in various ways,
including exponentiation::

    sage: M = matrix([[x,x^2],[1/x,x]])
    sage: M^2
    [x^2 + x   2*x^3]
    [      2 x^2 + x]
    sage: e^M
    [          1/2*(e^(2*sqrt(x)) + 1)*e^(x - sqrt(x)) 1/2*(x*e^(2*sqrt(x)) - x)*sqrt(x)*e^(x - sqrt(x))]
    [  1/2*(e^(2*sqrt(x)) - 1)*e^(x - sqrt(x))/x^(3/2)           1/2*(e^(2*sqrt(x)) + 1)*e^(x - sqrt(x))]

And complex exponentiation works now::

    sage: M = i*matrix([[pi]])
    sage: e^M
    [-1]
    sage: M = i*matrix([[pi,0],[0,2*pi]])
    sage: e^M
    [-1  0]
    [ 0  1]
    sage: M = matrix([[0,pi],[-pi,0]])
    sage: e^M
    [-1  0]
    [ 0 -1]

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
    doctest:...: DeprecationWarning: Substitution using function-call
    syntax and unnamed arguments is deprecated and will be removed
    from a future release of Sage; you can use named arguments instead,
    like EXPR(x=..., y=...)
    See http://trac.sagemath.org/5930 for details.
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
    (x, y) |--> cos(y)*log(x)

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
    4.646783762432936

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
    cos(1/2*I) + sin(I)
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

The symbolic calculus package uses its own copy of Maxima for
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

Simplifications like these are now done using Pynac::

    sage: x + x + pi
    pi + 2*x

But this still uses Maxima::

    sage: (x + x + pi).simplify()
    pi + 2*x

Note that ``x`` is still ``x``, since the
maxima used by the calculus package is different than the one in
the interactive interpreter.

Check to see that the problem with the variables method mentioned
in Trac ticket #3779 is actually fixed::

    sage: f = function('F',x)
    sage: diff(f*SR(1),x)
    D[0](F)(x)

Doubly ensure that Trac #7479 is working::

    sage: f(x)=x
    sage: integrate(f,x,0,1)
    1/2

Check that the problem with Taylor expansions of the gamma function
(Trac #9217) is fixed::

    sage: taylor(gamma(1/3+x),x,0,3)
    -1/432*((72*euler_gamma^3 + 36*euler_gamma^2*(sqrt(3)*pi + 9*log(3)) +
    27*pi^2*log(3) + 243*log(3)^3 + 18*euler_gamma*(6*sqrt(3)*pi*log(3) + pi^2
    + 27*log(3)^2 + 12*psi(1, 1/3)) + 324*log(3)*psi(1, 1/3) + sqrt(3)*(pi^3 +
    9*pi*(9*log(3)^2 + 4*psi(1, 1/3))))*gamma(1/3) - 72*psi(2,
    1/3)*gamma(1/3))*x^3 + 1/24*(6*sqrt(3)*pi*log(3) + 12*euler_gamma^2 + pi^2
    + 4*euler_gamma*(sqrt(3)*pi + 9*log(3)) + 27*log(3)^2 + 12*psi(1,
    1/3))*x^2*gamma(1/3) - 1/6*(6*euler_gamma + sqrt(3)*pi +
    9*log(3))*x*gamma(1/3) + gamma(1/3)
    sage: map(lambda f:f[0].n(), _.coeffs())  # numerical coefficients to make comparison easier; Maple 12 gives same answer
    [2.6789385347..., -8.3905259853..., 26.662447494..., -80.683148377...]

Ensure that ticket #8582 is fixed::

    sage: k = var("k")
    sage: sum(1/(1+k^2), k, -oo, oo)
    -1/2*I*psi(I + 1) + 1/2*I*psi(-I + 1) - 1/2*I*psi(I) + 1/2*I*psi(-I)

Ensure that ticket #8624 is fixed::

    sage: integrate(abs(cos(x)) * sin(x), x, pi/2, pi)
    1/2
    sage: integrate(sqrt(cos(x)^2 + sin(x)^2), x, 0, 2*pi)
    2*pi
"""

import re
from sage.rings.all import RR, Integer, CC, QQ, RealDoubleElement, algdep
from sage.rings.real_mpfr import create_RealNumber

from sage.misc.latex import latex
from sage.misc.parser import Parser

from sage.symbolic.ring import var, SR, is_SymbolicVariable
from sage.symbolic.expression import Expression
from sage.symbolic.function import Function
from sage.symbolic.function_factory import function_factory
from sage.symbolic.integration.integral import indefinite_integral, \
        definite_integral
import sage.symbolic.pynac

"""
Check if maxima has redundant variables defined after initialization #9538::

    sage: maxima = sage.interfaces.maxima.maxima
    sage: maxima('f1')
    f1
    sage: sage.calculus.calculus.maxima('f1')
    f1
"""
from sage.misc.lazy_import import lazy_import
lazy_import('sage.interfaces.maxima_lib','maxima')
# This is not the same instance of Maxima as the general purpose one
#from sage.interfaces.maxima import Maxima
#maxima = Maxima(init_code = ['display2d : false', 'domain : complex',
#                             'keepfloat : true', 'load(to_poly_solver)',
#                             'load(simplify_sum)'],
#                script_subdirectory=None)

########################################################
def symbolic_sum(expression, v, a, b, algorithm='maxima'):
    r"""
    Returns the symbolic sum `\sum_{v = a}^b expression` with respect
    to the variable `v` with endpoints `a` and `b`.

    INPUT:

    - ``expression`` - a symbolic expression

    - ``v`` - a variable or variable name

    - ``a`` - lower endpoint of the sum

    - ``b`` - upper endpoint of the sum

    - ``algorithm`` - (default: 'maxima')  one of
      - 'maxima' - use Maxima (the default)
      - 'maple' - (optional) use Maple
      - 'mathematica' - (optional) use Mathematica
      - 'giac' - (optional) use Giac

    EXAMPLES::

        sage: k, n = var('k,n')
        sage: from sage.calculus.calculus import symbolic_sum
        sage: symbolic_sum(k, k, 1, n).factor()
        1/2*(n + 1)*n

    ::

        sage: symbolic_sum(1/k^4, k, 1, oo)
        1/90*pi^4

    ::

        sage: symbolic_sum(1/k^5, k, 1, oo)
        zeta(5)

    A well known binomial identity::

        sage: symbolic_sum(binomial(n,k), k, 0, n)
        2^n

    And some truncations thereof::

        sage: assume(n>1)
        sage: symbolic_sum(binomial(n,k),k,1,n)
        2^n - 1
        sage: symbolic_sum(binomial(n,k),k,2,n)
        2^n - n - 1
        sage: symbolic_sum(binomial(n,k),k,0,n-1)
        2^n - 1
        sage: symbolic_sum(binomial(n,k),k,1,n-1)
        2^n - 2

    The binomial theorem::

        sage: x, y = var('x, y')
        sage: symbolic_sum(binomial(n,k) * x^k * y^(n-k), k, 0, n)
        (x + y)^n

    ::

        sage: symbolic_sum(k * binomial(n, k), k, 1, n)
        2^(n - 1)*n

    ::

        sage: symbolic_sum((-1)^k*binomial(n,k), k, 0, n)
        0

    ::

        sage: symbolic_sum(2^(-k)/(k*(k+1)), k, 1, oo)
        -log(2) + 1

    Summing a hypergeometric term::

        sage: symbolic_sum(binomial(n, k) * factorial(k) / factorial(n+1+k), k, 0, n)
        1/2*sqrt(pi)/factorial(n + 1/2)

    We check a well known identity::

        sage: bool(symbolic_sum(k^3, k, 1, n) == symbolic_sum(k, k, 1, n)^2)
        True

    A geometric sum::

        sage: a, q = var('a, q')
        sage: symbolic_sum(a*q^k, k, 0, n)
        (a*q^(n + 1) - a)/(q - 1)

    For the geometric series, we will have to assume
    the right values for the sum to converge::

        sage: assume(abs(q) < 1)
        sage: symbolic_sum(a*q^k, k, 0, oo)
        -a/(q - 1)

    A divergent geometric series.  Don't forget
    to forget your assumptions::

        sage: forget()
        sage: assume(q > 1)
        sage: symbolic_sum(a*q^k, k, 0, oo)
        Traceback (most recent call last):
        ...
        ValueError: Sum is divergent.
        sage: forget()
        sage: assumptions() # check the assumptions were really forgotten
        []

    This summation only Mathematica can perform::

        sage: symbolic_sum(1/(1+k^2), k, -oo, oo, algorithm = 'mathematica')     # optional - mathematica
        pi*coth(pi)

    An example of this summation with Giac::

        sage: symbolic_sum(1/(1+k^2), k, -oo, oo, algorithm = 'giac')           # optional - giac
        -(pi*e^(-2*pi) - pi*e^(2*pi))/(e^(-2*pi) + e^(2*pi) - 2)

    Use Maple as a backend for summation::

        sage: symbolic_sum(binomial(n,k)*x^k, k, 0, n, algorithm = 'maple')      # optional - maple
        (x + 1)^n

    TESTS:

    Trac #10564 is fixed::

        sage: sum (n^3 * x^n, n, 0, infinity)
        (x^3 + 4*x^2 + x)/(x^4 - 4*x^3 + 6*x^2 - 4*x + 1)

    .. note::

       Sage can currently only understand a subset of the output of Maxima,
       Maple and Mathematica, so even if the chosen backend can perform
       the summation the result might not be convertable into a Sage
       expression.
    """
    if not is_SymbolicVariable(v):
        if isinstance(v, str):
            v = var(v)
        else:
            raise TypeError("need a summation variable")

    if v in SR(a).variables() or v in SR(b).variables():
        raise ValueError("summation limits must not depend on the summation variable")

    if algorithm == 'maxima':
        return maxima.sr_sum(expression,v,a,b)

    elif algorithm == 'mathematica':
        try:
            sum = "Sum[%s, {%s, %s, %s}]" % tuple([repr(expr._mathematica_()) for expr in (expression, v, a, b)])
        except TypeError:
            raise ValueError("Mathematica cannot make sense of input")
        from sage.interfaces.mathematica import mathematica
        try:
            result = mathematica(sum)
        except TypeError:
            raise ValueError("Mathematica cannot make sense of: %s" % sum)
        return result.sage()

    elif algorithm == 'maple':
        sum = "sum(%s, %s=%s..%s)" % tuple([repr(expr._maple_()) for expr in (expression, v, a, b)])
        from sage.interfaces.maple import maple
        try:
            result = maple(sum).simplify()
        except TypeError:
            raise ValueError("Maple cannot make sense of: %s" % sum)
        return result.sage()

    elif algorithm == 'giac':
        sum = "sum(%s, %s, %s, %s)" % tuple([repr(expr._giac_()) for expr in (expression, v, a, b)])
        from sage.interfaces.giac import giac
        try:
            result = giac(sum)
        except TypeError:
            raise ValueError("Giac cannot make sense of: %s" % sum)
        return result.sage()

    else:
        raise ValueError("unknown algorithm: %s" % algorithm)

def nintegral(ex, x, a, b,
              desired_relative_error='1e-8',
              maximum_num_subintervals=200):
    r"""
    Return a floating point machine precision numerical approximation
    to the integral of ``self`` from `a` to
    `b`, computed using floating point arithmetic via maxima.

    INPUT:

    - ``x`` - variable to integrate with respect to

    - ``a`` - lower endpoint of integration

    - ``b`` - upper endpoint of integration

    - ``desired_relative_error`` - (default: '1e-8') the
      desired relative error

    - ``maximum_num_subintervals`` - (default: 200)
      maxima number of subintervals

    OUTPUT:

    - float: approximation to the integral

    - float: estimated absolute error of the
      approximation

    - the number of integrand evaluations

    - an error code:

      - ``0`` - no problems were encountered

      - ``1`` - too many subintervals were done

      - ``2`` - excessive roundoff error

      - ``3`` - extremely bad integrand behavior

      - ``4`` - failed to converge

      - ``5`` - integral is probably divergent or slowly
        convergent

      - ``6`` - the input is invalid; this includes the case of
                desired_relative_error being too small to be achieved

    ALIAS: nintegrate is the same as nintegral

    REMARK: There is also a function
    ``numerical_integral`` that implements numerical
    integration using the GSL C library. It is potentially much faster
    and applies to arbitrary user defined functions.

    Also, there are limits to the precision to which Maxima can compute
    the integral due to limitations in quadpack.
    In the following example, remark that the last value of the returned
    tuple is ``6``, indicating that the input was invalid, in this case
    because of a too high desired precision.

    ::

        sage: f = x
        sage: f.nintegral(x,0,1,1e-14)
        (0.0, 0.0, 0, 6)

    EXAMPLES::

        sage: f(x) = exp(-sqrt(x))
        sage: f.nintegral(x, 0, 1)
        (0.5284822353142306, 4.163...e-11, 231, 0)

    We can also use the ``numerical_integral`` function,
    which calls the GSL C library.

    ::

        sage: numerical_integral(f, 0, 1)
        (0.528482232253147, 6.83928460...e-07)

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
        (-480.00000000000006, 5.329070518200754e-12, 21, 0)

    It is just because every floating point evaluation of return -480.0
    in floating point.

    Important note: using PARI/GP one can compute numerical integrals
    to high precision::

        sage: gp.eval('intnum(x=17,42,exp(-x^2)*log(x))')
        '2.565728500561051482917356396 E-127'        # 32-bit
        '2.5657285005610514829173563961304785900 E-127'    # 64-bit
        sage: old_prec = gp.set_real_precision(50)
        sage: gp.eval('intnum(x=17,42,exp(-x^2)*log(x))')
        '2.5657285005610514829173563961304785900147709554020 E-127'
        sage: gp.set_real_precision(old_prec)
        50

    Note that the input function above is a string in PARI syntax.
    """
    try:
        v = ex._maxima_().quad_qags(x, a, b,
                                    epsrel=desired_relative_error,
                                    limit=maximum_num_subintervals)
    except TypeError as err:
        if "ERROR" in str(err):
            raise ValueError("Maxima (via quadpack) cannot compute the integral")
        else:
            raise TypeError(err)

    # Maxima returns unevaluated expressions when the underlying library fails
    # to perfom numerical integration. See:
    # http://www.math.utexas.edu/pipermail/maxima/2008/012975.html
    if 'quad_qags' in str(v):
        raise ValueError("Maxima (via quadpack) cannot compute the integral")

    return float(v[0]), float(v[1]), Integer(v[2]), Integer(v[3])

nintegrate = nintegral

def minpoly(ex, var='x', algorithm=None, bits=None, degree=None, epsilon=0):
    r"""
    Return the minimal polynomial of self, if possible.

    INPUT:

    - ``var`` - polynomial variable name (default 'x')

    - ``algorithm`` - 'algebraic' or 'numerical' (default
      both, but with numerical first)

    - ``bits`` - the number of bits to use in numerical
      approx

    - ``degree`` - the expected algebraic degree

    - ``epsilon`` - return without error as long as
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
    algorithm parameter. By default, the numerical algorithm is
    attempted first, then the algebraic one.

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

    The minpoly function is used implicitly when creating
    number fields::

        sage: x = var('x')
        sage: eqn =  x^3 + sqrt(2)*x + 5 == 0
        sage: a = solve(eqn, x)[0].rhs()
        sage: QQ[a]
        Number Field in a with defining polynomial x^6 + 10*x^3 - 2*x^2 + 25

    Here we solve a cubic and then recover it from its complicated
    radical expansion.

    ::

        sage: f = x^3 - x + 1
        sage: a = f.solve(x)[0].rhs(); a
        -1/2*(1/18*sqrt(23)*sqrt(3) - 1/2)^(1/3)*(I*sqrt(3) + 1) - 1/6*(-I*sqrt(3) + 1)/(1/18*sqrt(23)*sqrt(3) - 1/2)^(1/3)
        sage: a.minpoly()
        x^3 - x + 1

    Note that simplification may be necessary to see that the minimal
    polynomial is correct.

    ::

        sage: a = sqrt(2)+sqrt(3)+sqrt(5)
        sage: f = a.minpoly(); f
        x^8 - 40*x^6 + 352*x^4 - 960*x^2 + 576
        sage: f(a)
        ((((sqrt(5) + sqrt(3) + sqrt(2))^2 - 40)*(sqrt(5) + sqrt(3) + sqrt(2))^2 + 352)*(sqrt(5) + sqrt(3) + sqrt(2))^2 - 960)*(sqrt(5) + sqrt(3) + sqrt(2))^2 + 576
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

        sage: cos(pi/33).minpoly(algorithm='algebraic')
        x^10 + 1/2*x^9 - 5/2*x^8 - 5/4*x^7 + 17/8*x^6 + 17/16*x^5 - 43/64*x^4 - 43/128*x^3 + 3/64*x^2 + 3/128*x + 1/1024
        sage: cos(pi/33).minpoly(algorithm='numerical')
        Traceback (most recent call last):
        ...
        NotImplementedError: Could not prove minimal polynomial x^10 + 1/2*x^9 - 5/2*x^8 - 5/4*x^7 + 17/8*x^6 + 17/16*x^5 - 43/64*x^4 - 43/128*x^3 + 3/64*x^2 + 3/128*x + 1/1024 (epsilon ...)

    Sometimes it fails, as it must given that some numbers aren't algebraic::

        sage: sin(1).minpoly(algorithm='numerical')
        Traceback (most recent call last):
        ...
        ValueError: Could not find minimal polynomial (1000 bits, degree 24).

    .. note::

       Of course, failure to produce a minimal polynomial does not
       necessarily indicate that this number is transcendental.
    """
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
                    # Degree might have been an over-estimate,
                    # factor because we want (irreducible) minpoly.
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
                            elif algorithm is not None:
                                raise NotImplementedError("Could not prove minimal polynomial %s (epsilon %s)" % (g, RR(error).str(no_sci=False)))

        if algorithm is not None:
            raise ValueError("Could not find minimal polynomial (%s bits, degree %s)." % (bits, degree))

    if algorithm is None or algorithm == 'algebraic':
        from sage.rings.all import QQbar
        return QQ[var](QQbar(ex).minpoly())

    raise ValueError("Unknown algorithm: %s" % algorithm)


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

    - ``dir`` - (default: None); dir may have the value
      'plus' (or '+' or 'right') for a limit from above,
      'minus' (or '-' or 'left') for a limit from below, or may be omitted
      (implying a two-sided limit is to be computed).

    - ``taylor`` - (default: False); if True, use Taylor
      series, which allows more limits to be computed (but may also
      crash in some obscure cases due to bugs in Maxima).

    - ``**argv`` - 1 named parameter

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

    Notice that Maxima may ask for more information::

        sage: var('a')
        a
        sage: limit(x^a,x=0)
        Traceback (most recent call last):
        ...
        ValueError: Computation failed since Maxima requested additional
        constraints; using the 'assume' command before limit evaluation
        *may* help (see `assume?` for more details)
        Is a positive, negative or zero?

    With this example, Maxima is looking for a LOT of information::

        sage: assume(a>0)
        sage: limit(x^a,x=0)
        Traceback (most recent call last):
        ...
        ValueError: Computation failed since Maxima requested additional
        constraints; using the 'assume' command before limit evaluation
        *may* help (see `assume?` for more details)
        Is a an integer?
        sage: assume(a,'integer')
        sage: limit(x^a,x=0)
        Traceback (most recent call last):
        ...
        ValueError: Computation failed since Maxima requested additional
        constraints; using the 'assume' command before limit evaluation
        *may* help (see `assume?` for more details)
        Is a an even number?
        sage: assume(a,'even')
        sage: limit(x^a,x=0)
        0
        sage: forget()

    More examples::

        sage: limit(x*log(x), x = 0, dir='+')
        0
        sage: lim((x+1)^(1/x), x = 0)
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
        sage: limit(e^(-1/x), x=0, dir='right')
        0
        sage: limit(e^(-1/x), x=0, dir='left')
        +Infinity

    ::

        sage: f = log(log(x))/log(x)
        sage: forget(); assume(x<-2); lim(f, x=0, taylor=True)
        0
        sage: forget()

    Here ind means "indefinite but bounded"::

        sage: lim(sin(1/x), x = 0)
        ind

    TESTS::

        sage: lim(x^2, x=2, dir='nugget')
        Traceback (most recent call last):
        ...
        ValueError: dir must be one of None, 'plus', '+', 'right',
        'minus', '-', 'left'

    We check that :trac:`3718` is fixed, so that
    Maxima gives correct limits for the floor function::

        sage: limit(floor(x), x=0, dir='-')
        -1
        sage: limit(floor(x), x=0, dir='+')
        0
        sage: limit(floor(x), x=0)
        und

    Maxima gives the right answer here, too, showing
    that :trac:`4142` is fixed::

        sage: f = sqrt(1-x^2)
        sage: g = diff(f, x); g
        -x/sqrt(-x^2 + 1)
        sage: limit(g, x=1, dir='-')
        -Infinity

    ::

        sage: limit(1/x, x=0)
        Infinity
        sage: limit(1/x, x=0, dir='+')
        +Infinity
        sage: limit(1/x, x=0, dir='-')
        -Infinity

    Check that :trac:`8942` is fixed::

        sage: f(x) = (cos(pi/4-x) - tan(x)) / (1 - sin(pi/4+x))
        sage: limit(f(x), x = pi/4, dir='minus')
        +Infinity
        sage: limit(f(x), x = pi/4, dir='plus')
        -Infinity
        sage: limit(f(x), x = pi/4)
        Infinity

    Check that we give deprecation warnings for 'above' and 'below',
    :trac:`9200`::

        sage: limit(1/x, x=0, dir='above')
        doctest:...: DeprecationWarning: the keyword
        'above' is deprecated. Please use 'right' or '+' instead.
        See http://trac.sagemath.org/9200 for details.
        +Infinity
        sage: limit(1/x, x=0, dir='below')
        doctest:...: DeprecationWarning: the keyword
        'below' is deprecated. Please use 'left' or '-' instead.
        See http://trac.sagemath.org/9200 for details.
        -Infinity

    Check that :trac:`12708` is fixed::

        sage: limit(tanh(x),x=0)
        0

    Check that :trac:`15386` is fixed::

        sage: n = var('n')
        sage: assume(n>0)
        sage: sequence = -(3*n^2 + 1)*(-1)^n/sqrt(n^5 + 8*n^3 + 8)
        sage: limit(sequence, n=infinity)
        0
    """
    if not isinstance(ex, Expression):
        ex = SR(ex)

    if len(argv) != 1:
        raise ValueError("call the limit function like this, e.g. limit(expr, x=2).")
    else:
        k = argv.keys()[0]
        v = var(k)
        a = argv[k]

    if taylor and algorithm == 'maxima':
        algorithm = 'maxima_taylor'

    if dir not in [None, 'plus', '+', 'right', 'minus', '-', 'left',
            'above', 'below']:
        raise ValueError("dir must be one of None, 'plus', '+', 'right', 'minus', '-', 'left'")

    if algorithm == 'maxima':
        if dir is None:
            l = maxima.sr_limit(ex, v, a)
        elif dir in ['plus', '+', 'right', 'above']:
            if dir == 'above':
                from sage.misc.superseded import deprecation
                deprecation(9200, "the keyword 'above' is deprecated. Please use 'right' or '+' instead.")
            l = maxima.sr_limit(ex, v, a, 'plus')
        elif dir in ['minus', '-', 'left', 'below']:
            if dir == 'below':
                from sage.misc.superseded import deprecation
                deprecation(9200, "the keyword 'below' is deprecated. Please use 'left' or '-' instead.")
            l = maxima.sr_limit(ex, v, a, 'minus')
    elif algorithm == 'maxima_taylor':
        if dir is None:
            l = maxima.sr_tlimit(ex, v, a)
        elif dir == 'plus' or dir == 'above' or dir == 'from_right':
            l = maxima.sr_tlimit(ex, v, a, 'plus')
        elif dir == 'minus' or dir == 'below' or dir == 'from_left':
            l = maxima.sr_tlimit(ex, v, a, 'minus')
    elif algorithm == 'sympy':
        if dir is None:
            import sympy
            l = sympy.limit(ex._sympy_(), v._sympy_(), a._sympy_())
        else:
            raise NotImplementedError("sympy does not support one-sided limits")

    #return l.sage()
    return ex.parent()(l)

# lim is alias for limit
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

    DEFINITION:

    The Laplace transform of a function `f(t)`,
    defined for all real numbers `t \geq 0`, is the function
    `F(s)` defined by

    .. math::

                      F(s) = \int_{0}^{\infty} e^{-st} f(t) dt.

    EXAMPLES:

    We compute a few Laplace transforms::

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

    EXAMPLES:

    A BATTLE BETWEEN the X-women and the Y-men (by David
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
        -91/2*e^(4*t) + 629/2*e^(-4*t) + 1
        sage: yt
        91/8*e^(4*t) + 629/8*e^(-4*t)
        sage: p1 = plot(xt,0,1/2,rgbcolor=(1,0,0))
        sage: p2 = plot(yt,0,1/2,rgbcolor=(0,1,0))
        sage: (p1+p2).save(os.path.join(SAGE_TMP, "de_plot.png"))

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
    if not isinstance(ex, (Expression, Function)):
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

###################################################################
# symbolic evaluation "at" a point
###################################################################
def at(ex, *args, **kwds):
    """
    Parses ``at`` formulations from other systems, such as Maxima.
    Replaces evaluation 'at' a point with substitution method of
    a symbolic expression.

    EXAMPLES:

    We do not import ``at`` at the top level, but we can use it
    as a synonym for substitution if we import it::

        sage: g = x^3-3
        sage: from sage.calculus.calculus import at
        sage: at(g, x=1)
        -2
        sage: g.subs(x=1)
        -2

    We find a formal Taylor expansion::

        sage: h,x = var('h,x')
        sage: u = function('u')
        sage: u(x + h)
        u(h + x)
        sage: diff(u(x+h), x)
        D[0](u)(h + x)
        sage: taylor(u(x+h),h,0,4)
        1/24*h^4*D[0, 0, 0, 0](u)(x) + 1/6*h^3*D[0, 0, 0](u)(x) + 1/2*h^2*D[0, 0](u)(x) + h*D[0](u)(x) + u(x)

    We compute a Laplace transform::

        sage: var('s,t')
        (s, t)
        sage: f=function('f', t)
        sage: f.diff(t,2)
        D[0, 0](f)(t)
        sage: f.diff(t,2).laplace(t,s)
        s^2*laplace(f(t), t, s) - s*f(0) - D[0](f)(0)

    We can also accept a non-keyword list of expression substitutions,
    like Maxima does, :trac:`12796`::

        sage: from sage.calculus.calculus import at
        sage: f = function('f')
        sage: at(f(x), [x == 1])
        f(1)

    TESTS:

    Our one non-keyword argument must be a list::

        sage: from sage.calculus.calculus import at
        sage: f = function('f')
        sage: at(f(x), x == 1)
        Traceback (most recent call last):
        ...
        TypeError: at can take at most one argument, which must be a list

    We should convert our first argument to a symbolic expression::

        sage: from sage.calculus.calculus import at
        sage: at(int(1), x=1)
        1

    """
    if not isinstance(ex, (Expression, Function)):
        ex = SR(ex)
    if len(args) == 1 and isinstance(args[0],list):
        for c in args[0]:
            kwds[str(c.lhs())]=c.rhs()
    else:
        if len(args) !=0:
            raise TypeError("at can take at most one argument, which must be a list")

    return ex.subs(**kwds)


#############################################3333
def var_cmp(x,y):
    """
    Return comparison of the two variables x and y, which is just the
    comparison of the underlying string representations of the
    variables. This is used internally by the Calculus package.

    INPUT:

    - ``x, y`` - symbolic variables

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

def dummy_limit(*args):
    """
    This function is called to create formal wrappers of limits that
    Maxima can't compute:

    EXAMPLES::

        sage: a = lim(exp(x^2)*(1-erf(x)), x=infinity); a
        -limit((erf(x) - 1)*e^(x^2), x, +Infinity)
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
        return definite_integral(*args, hold=True)
    else:
        return indefinite_integral(*args, hold=True)

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

def _limit_latex_(self, f, x, a, direction=None):
    r"""
    Return latex expression for limit of a symbolic function.

    EXAMPLES::

        sage: from sage.calculus.calculus import _limit_latex_
        sage: var('x,a')
        (x, a)
        sage: f = function('f')
        sage: _limit_latex_(0, f(x), x, a)
        '\\lim_{x \\to a}\\, f\\left(x\\right)'
        sage: latex(limit(f(x), x=oo))
        \lim_{x \to +\infty}\, f\left(x\right)

    TESTS:

    When one-sided limits are converted back from maxima, the direction
    argument becomes a symbolic variable. We check if typesetting these works::

        sage: var('minus,plus')
        (minus, plus)
        sage: _limit_latex_(0, f(x), x, a, minus)
        '\\lim_{x \\to a^-}\\, f\\left(x\\right)'
        sage: _limit_latex_(0, f(x), x, a, plus)
        '\\lim_{x \\to a^+}\\, f\\left(x\\right)'
        sage: latex(limit(f(x),x=a,dir='+'))
        \lim_{x \to a^+}\, f\left(x\right)
        sage: latex(limit(f(x),x=a,dir='right'))
        \lim_{x \to a^+}\, f\left(x\right)
        sage: latex(limit(f(x),x=a,dir='-'))
        \lim_{x \to a^-}\, f\left(x\right)
        sage: latex(limit(f(x),x=a,dir='left'))
        \lim_{x \to a^-}\, f\left(x\right)

    Check if :trac:`13181` is fixed::

        sage: t = var('t')
        sage: latex(limit(exp_integral_e(1/2, I*t - I*x)*sqrt(-t + x),t=x,dir='-'))
        \lim_{t \to x^-}\, \sqrt{-t + x} exp_integral_e\left(\frac{1}{2}, i \, t - i \, x\right)
        sage: latex(limit(exp_integral_e(1/2, I*t - I*x)*sqrt(-t + x),t=x,dir='+'))
        \lim_{t \to x^+}\, \sqrt{-t + x} exp_integral_e\left(\frac{1}{2}, i \, t - i \, x\right)
        sage: latex(limit(exp_integral_e(1/2, I*t - I*x)*sqrt(-t + x),t=x))
        \lim_{t \to x}\, \sqrt{-t + x} exp_integral_e\left(\frac{1}{2}, i \, t - i \, x\right)
    """
    if repr(direction) == 'minus':
        dir_str = '^-'
    elif repr(direction) == 'plus':
        dir_str = '^+'
    else:
        dir_str = ''
    return "\\lim_{%s \\to %s%s}\\, %s"%(latex(x), latex(a), dir_str, latex(f))

def _laplace_latex_(self, *args):
    r"""
    Return LaTeX expression for Laplace transform of a symbolic function.

    EXAMPLES::

        sage: from sage.calculus.calculus import _laplace_latex_
        sage: var('s,t')
        (s, t)
        sage: f = function('f',t)
        sage: _laplace_latex_(0,f,t,s)
        '\\mathcal{L}\\left(f\\left(t\\right), t, s\\right)'
        sage: latex(laplace(f, t, s))
        \mathcal{L}\left(f\left(t\right), t, s\right)

    """
    return "\\mathcal{L}\\left(%s\\right)"%(', '.join([latex(x) for x in args]))

def _inverse_laplace_latex_(self, *args):
    r"""
    Return LaTeX expression for inverse Laplace transform
    of a symbolic function.

    EXAMPLES::

        sage: from sage.calculus.calculus import _inverse_laplace_latex_
        sage: var('s,t')
        (s, t)
        sage: F = function('F',s)
        sage: _inverse_laplace_latex_(0,F,s,t)
        '\\mathcal{L}^{-1}\\left(F\\left(s\\right), s, t\\right)'
        sage: latex(inverse_laplace(F,s,t))
        \mathcal{L}^{-1}\left(F\left(s\right), s, t\right)
    """
    return "\\mathcal{L}^{-1}\\left(%s\\right)"%(', '.join([latex(x) for x in args]))

# Return un-evaluated expression as instances of SFunction class
_limit = function_factory('limit', print_latex_func=_limit_latex_)
_laplace = function_factory('laplace', print_latex_func=_laplace_latex_)
_inverse_laplace = function_factory('ilt',
        print_latex_func=_inverse_laplace_latex_)

######################################i################




#######################################################

symtable = {'%pi':'pi', '%e': 'e', '%i':'I', '%gamma':'euler_gamma'}

from sage.misc.multireplace import multiple_replace
import re

maxima_tick = re.compile("'[a-z|A-Z|0-9|_]*")

maxima_qp = re.compile("\?\%[a-z|A-Z|0-9|_]*")  # e.g., ?%jacobi_cd

maxima_var = re.compile("\%[a-z|A-Z|0-9|_]*")  # e.g., ?%jacobi_cd

sci_not = re.compile("(-?(?:0|[1-9]\d*))(\.\d+)?([eE][-+]\d+)")

polylog_ex = re.compile('li\[([0-9]+?)\]\(')

maxima_polygamma = re.compile("psi\[(\d*)\]\(")  # matches psi[n]( where n is a number

def symbolic_expression_from_maxima_string(x, equals_sub=False, maxima=maxima):
    """
    Given a string representation of a Maxima expression, parse it and
    return the corresponding Sage symbolic expression.

    INPUT:

    - ``x`` - a string

    - ``equals_sub`` - (default: False) if True, replace
      '=' by '==' in self

    - ``maxima`` - (default: the calculus package's
      Maxima) the Maxima interpreter to use.

    EXAMPLES::

        sage: from sage.calculus.calculus import symbolic_expression_from_maxima_string as sefms
        sage: sefms('x^%e + %e^%pi + %i + sin(0)')
        x^e + e^pi + I
        sage: f = function('f',x)
        sage: sefms('?%at(f(x),x=2)#1')
        f(2) != 1
        sage: a = sage.calculus.calculus.maxima("x#0"); a
        x#0
        sage: a.sage()
        x != 0

    TESTS:

    Trac #8459 fixed::

        sage: maxima('3*li[2](u)+8*li[33](exp(u))').sage()
        8*polylog(33, e^u) + 3*polylog(2, u)

    Check if #8345 is fixed::

        sage: assume(x,'complex')
        sage: t = x.conjugate()
        sage: latex(t)
        \overline{x}
        sage: latex(t._maxima_()._sage_())
        \overline{x}

    Check that we can understand maxima's not-equals (:trac:`8969`)::

        sage: from sage.calculus.calculus import symbolic_expression_from_maxima_string as sefms
        sage: sefms("x!=3") == (factorial(x) == 3)
        True
        sage: sefms("x # 3") == SR(x != 3)
        True
        sage: solve([x != 5], x)
        #0: solve_rat_ineq(ineq=x # 5)
        [[x - 5 != 0]]
        sage: solve([2*x==3, x != 5], x)
        [[x == (3/2), (-7/2) != 0]]

    """
    syms = sage.symbolic.pynac.symbol_table.get('maxima', {}).copy()

    if len(x) == 0:
        raise RuntimeError("invalid symbolic expression -- ''")
    maxima.set('_tmp_',x)

    # This is inefficient since it so rarely is needed:
    #r = maxima._eval_line('listofvars(_tmp_);')[1:-1]

    s = maxima._eval_line('_tmp_;')

    formal_functions = maxima_tick.findall(s)
    if len(formal_functions) > 0:
        for X in formal_functions:
            syms[X[1:]] = function_factory(X[1:])
        # You might think there is a potential very subtle bug if 'foo
        # is in a string literal -- but string literals should *never*
        # ever be part of a symbolic expression.
        s = s.replace("'","")

    delayed_functions = maxima_qp.findall(s)
    if len(delayed_functions) > 0:
        for X in delayed_functions:
            if X == '?%at': # we will replace Maxima's "at" with symbolic evaluation, not an SFunction
                pass
            else:
                syms[X[2:]] = function_factory(X[2:])
        s = s.replace("?%","")

    s = polylog_ex.sub('polylog(\\1,',s)
    s = multiple_replace(symtable, s)
    s = s.replace("%","")

    s = s.replace("#","!=") # a lot of this code should be refactored somewhere...

    s = maxima_polygamma.sub('psi(\g<1>,',s) # this replaces psi[n](foo) with psi(n,foo), ensuring that derivatives of the digamma function are parsed properly below

    if equals_sub:
        s = s.replace('=','==')
        # unfortunately, this will turn != into !==, which we correct
        s = s.replace("!==", "!=")

    #replace %union from to_poly_solve with a list
    if s[0:5]=='union':
        s = s[5:]
        s = s[s.find("(")+1:s.rfind(")")]
        s = "[" + s + "]" # turn it into a string that looks like a list

    #replace %solve from to_poly_solve with the expressions
    if s[0:5]=='solve':
        s = s[5:]
        s = s[s.find("(")+1:s.find("]")+1]

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

    syms['at'] = at

    global is_simplified
    try:
        # use a global flag so all expressions obtained via
        # evaluation of maxima code are assumed pre-simplified
        is_simplified = True
        return symbolic_expression_from_string(s, syms, accept_sequence=True)
    except SyntaxError:
        raise TypeError("unable to make sense of Maxima expression '%s' in Sage"%s)
    finally:
        is_simplified = False


# Comma format options for Maxima
def mapped_opts(v):
    """
    Used internally when creating a string of options to pass to
    Maxima.

    INPUT:

    - ``v`` - an object

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

# We keep two dictionaries syms_cur and syms_default to keep the current symbol
# table and the state of the table at startup respectively. These are used by
# the restore() function (see sage.misc.reset).
#
# The dictionary _syms is used as a lookup table for the system function
# registry by _find_func() below. It gets updated by
# symbolic_expression_from_string() before calling the parser.
from sage.symbolic.pynac import symbol_table
_syms = syms_cur = symbol_table.get('functions', {})
syms_default = dict(syms_cur)

# This dictionary is used to pass a lookup table other than the system registry
# to the parser. A global variable is necessary since the parser calls the
# _find_var() and _find_func() functions below without extra arguments.
_augmented_syms = {}

from sage.symbolic.ring import pynac_symbol_registry

def _find_var(name):
    """
    Function to pass to Parser for constructing
    variables from strings.  For internal use.

    EXAMPLES::

        sage: y = var('y')
        sage: sage.calculus.calculus._find_var('y')
        y
        sage: sage.calculus.calculus._find_var('I')
        I
    """
    try:
        res = _augmented_syms.get(name)
        if res is None:
            return pynac_symbol_registry[name]
        # _augmented_syms might contain entries pointing to functions if
        # previous computations polluted the maxima workspace
        if not isinstance(res, Function):
            return res
    except KeyError:
        pass

    # try to find the name in the global namespace
    # needed for identifiers like 'e', etc.
    try:
        return SR(sage.all.__dict__[name])
    except (KeyError, TypeError):
        return var(name)

def _find_func(name, create_when_missing = True):
    """
    Function to pass to Parser for constructing
    functions from strings.  For internal use.

    EXAMPLES::

        sage: sage.calculus.calculus._find_func('limit')
        limit
        sage: sage.calculus.calculus._find_func('zeta_zeros')
        zeta_zeros
        sage: f(x)=sin(x)
        sage: sage.calculus.calculus._find_func('f')
        f
        sage: sage.calculus.calculus._find_func('g', create_when_missing=False)
        sage: s = sage.calculus.calculus._find_func('sin')
        sage: s(0)
        0
    """
    try:
        func = _augmented_syms.get(name)
        if func is None:
            func = _syms[name]
        if not isinstance(func, Expression):
            return func
    except KeyError:
        pass
    try:
        func = SR(sage.all.__dict__[name])
        if not isinstance(func, Expression):
            return func
    except (KeyError, TypeError):
        if create_when_missing:
            return function_factory(name)
        else:
            return None

SR_parser = Parser(make_int      = lambda x: SR(Integer(x)),
                   make_float    = lambda x: SR(RealDoubleElement(x)),
                   make_var      = _find_var,
                   make_function = _find_func)

def symbolic_expression_from_string(s, syms=None, accept_sequence=False):
    """
    Given a string, (attempt to) parse it and return the
    corresponding Sage symbolic expression.  Normally used
    to return Maxima output to the user.

    INPUT:

    - ``s`` - a string

    - ``syms`` - (default: None) dictionary of
      strings to be regarded as symbols or functions

    - ``accept_sequence`` - (default: False) controls whether
      to allow a (possibly nested) set of lists and tuples
      as input

    EXAMPLES::

        sage: y = var('y')
        sage: sage.calculus.calculus.symbolic_expression_from_string('[sin(0)*x^2,3*spam+e^pi]',syms={'spam':y},accept_sequence=True)
        [0, 3*y + e^pi]
    """
    global _syms
    _syms = sage.symbolic.pynac.symbol_table['functions'].copy()
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
