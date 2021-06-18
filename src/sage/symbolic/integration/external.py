"""Symbolic Integration via External Software

TESTS::

    sage: from sage.symbolic.integration.external import sympy_integrator
    sage: sympy_integrator(sin(x), x)
    -cos(x)
"""
from sage.symbolic.expression import Expression
from sage.symbolic.ring import SR


def maxima_integrator(expression, v, a=None, b=None):
    """
    Integration using Maxima

    EXAMPLES::

        sage: from sage.symbolic.integration.external import maxima_integrator
        sage: maxima_integrator(sin(x), x)
        -cos(x)
        sage: maxima_integrator(cos(x), x)
        sin(x)
        sage: f(x) = function('f')(x)
        sage: maxima_integrator(f(x), x)
        integrate(f(x), x)

    TESTS:

    Check that :trac:`25817` is fixed::

        sage: maxima_integrator(log(e^x*log(x)*sin(x))/x^2, x)
        1/2*(x*(Ei(-log(x)) + conjugate(Ei(-log(x))))
        - 2*x*integrate(sin(x)/(x*cos(x)^2 + x*sin(x)^2
        + 2*x*cos(x) + x), x) + 2*x*integrate(sin(x)/(x*cos(x)^2
        + x*sin(x)^2 - 2*x*cos(x) + x), x) + 2*x*log(x) + 2*log(2)
        - log(cos(x)^2 + sin(x)^2 + 2*cos(x) + 1) - log(cos(x)^2
        + sin(x)^2 - 2*cos(x) + 1) - 2*log(log(x)))/x
    """
    from sage.calculus.calculus import maxima
    if not isinstance(expression, Expression):
        expression = SR(expression)
    if a is None:
        result = maxima.sr_integral(expression, v)
    else:
        result = maxima.sr_integral(expression, v, a, b)
    return result._sage_()


def sympy_integrator(expression, v, a=None, b=None):
    """
    Integration using SymPy

    EXAMPLES::

        sage: from sage.symbolic.integration.external import sympy_integrator
        sage: sympy_integrator(sin(x), x)
        -cos(x)
        sage: sympy_integrator(cos(x), x)
        sin(x)
    """
    import sympy
    ex = expression._sympy_()
    v = v._sympy_()
    if a is None:
        result = sympy.integrate(ex, v)
    else:
        result = sympy.integrate(ex, (v, a._sympy_(), b._sympy_()))
    return result._sage_()


def mma_free_integrator(expression, v, a=None, b=None):
    """
    Integration using Mathematica's online integrator

    EXAMPLES::

        sage: from sage.symbolic.integration.external import mma_free_integrator
        sage: mma_free_integrator(sin(x), x) # optional - internet
        -cos(x)

    A definite integral::

        sage: mma_free_integrator(e^(-x), x, a=0, b=oo) # optional - internet
        1

    TESTS:

    Check that :trac:`18212` is resolved::

        sage: var('y')   # optional - internet
        y
        sage: result = integral(sin(y)^2, y, algorithm='mathematica_free') # optional - internet
        sage: result.simplify_trig()               # optional - internet
        -1/2*cos(y)*sin(y) + 1/2*y

    ::

    Check that :trac:`14764` is resolved::

        sage: integrate(x^2, x, 0, 1, algorithm="mathematica_free") # optional - internet
        1/3
        sage: integrate(sin(x), [x, 0, pi], algorithm="mathematica_free") # optional - internet
        2
        sage: integrate(sqrt(x), (x, 0, 1), algorithm="mathematica_free") # optional - internet
        2/3

    ::

        sage: mma_free_integrator(exp(-x^2)*log(x), x) # optional - internet
        1/2*sqrt(pi)*erf(x)*log(x) - x*hypergeometric((1/2, 1/2), (3/2, 3/2), -x^2)


    """
    from sage.interfaces.mathematica import request_wolfram_alpha, parse_moutput_from_json, symbolic_expression_from_mathematica_string
    math_expr = expression._mathematica_init_()
    variable = v._mathematica_init_()
    if a is None and b is None:
        input = "Integrate[{},{}]".format(math_expr, variable)
    elif a is not None and b is not None:
        input = "Integrate[{},{{{},{},{}}}]".format(math_expr, variable,
                    a._mathematica_init_(), b._mathematica_init_())
    else:
        raise ValueError('a(={}) and b(={}) should be both None'
                         ' or both defined'.format(a, b))
    json_page_data = request_wolfram_alpha(input)
    all_outputs = parse_moutput_from_json(json_page_data)
    if not all_outputs:
        raise ValueError("no outputs found in the answer from Wolfram Alpha")
    first_output = all_outputs[0]
    return symbolic_expression_from_mathematica_string(first_output)


def fricas_integrator(expression, v, a=None, b=None, noPole=True):
    """
    Integration using FriCAS

    EXAMPLES::

        sage: from sage.symbolic.integration.external import fricas_integrator  # optional - fricas
        sage: fricas_integrator(sin(x), x)                                      # optional - fricas
        -cos(x)
        sage: fricas_integrator(cos(x), x)                                      # optional - fricas
        sin(x)
        sage: fricas_integrator(1/(x^2-2), x, 0, 1)                             # optional - fricas
        1/4*sqrt(2)*(log(3*sqrt(2) - 4) - log(sqrt(2)))
        sage: fricas_integrator(1/(x^2+6), x, -oo, oo)                          # optional - fricas
        1/6*sqrt(6)*pi

    TESTS:

    Check that :trac:`25220` is fixed::

        sage: integral(sqrt(1-cos(x)), x, 0, 2*pi, algorithm="fricas")          # optional - fricas
        4*sqrt(2)

    Check that in case of failure one gets unevaluated integral::

        sage: integral(cos(ln(cos(x))), x, 0, pi/8, algorithm='fricas')   # optional - fricas
        integrate(cos(log(cos(x))), x, 0, 1/8*pi)
        sage: integral(cos(ln(cos(x))), x, algorithm='fricas')   # optional - fricas
        integral(cos(log(cos(x))), x)
    """
    if not isinstance(expression, Expression):
        expression = SR(expression)

    from sage.interfaces.fricas import fricas
    ex = fricas(expression)

    if a is None:
        result = ex.integrate(v)
    else:
        seg = fricas.equation(v, fricas.segment(a, b))

        if noPole:
            result = ex.integrate(seg, '"noPole"')
        else:
            result = ex.integrate(seg)

    result = result.sage()

    if result == "failed":
        return expression.integrate(v, a, b, hold=True)

    if result == "potentialPole":
        raise ValueError("The integrand has a potential pole"
                         " in the integration interval")

    return result


def giac_integrator(expression, v, a=None, b=None):
    """
    Integration using Giac

    EXAMPLES::

        sage: from sage.symbolic.integration.external import giac_integrator
        sage: giac_integrator(sin(x), x)
        -cos(x)
        sage: giac_integrator(1/(x^2+6), x, -oo, oo)
        1/6*sqrt(6)*pi

    TESTS::

        sage: giac_integrator(e^(-x^2)*log(x), x)
        integrate(e^(-x^2)*log(x), x)
    """
    ex = expression._giac_()
    if a is None:
        result = ex.integrate(v._giac_())
    else:
        result = ex.integrate(v._giac_(), a._giac_(), b._giac_())
    if 'integrate' in format(result) or 'integration' in format(result):
        return expression.integrate(v, a, b, hold=True)
    else:
        return result._sage_()
