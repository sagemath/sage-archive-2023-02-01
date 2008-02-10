from equations import (is_SymbolicEquation,
                       forget, assume, assumptions,
                       solve, solve_mod)

from calculus import (SymbolicExpressionRing,
                      is_SymbolicExpressionRing,
                      is_SymbolicExpression,
                      is_SymbolicVariable,
                      CallableSymbolicExpressionRing,
                      is_CallableSymbolicExpressionRing,
                      is_CallableSymbolicExpression,
                      SR,
                      sin, cos, sec, csc, cot, tan, log, erf, sqrt,
                      tanh, sinh, cosh, coth, sech, csch, ln,
                      asin, acos, atan,
                      asinh, acosh, atanh, acoth, asech, acsch,
                      acot, acsc, asec,
                      arcsin, arccos, arctan,
                      arcsinh, arccosh, arctanh, arccoth, arcsech, arccsch,
                      arccot, arccsc, arcsec,

                      ceil, floor,
                      polylog,
                      abs_symbolic, exp,
                      is_SymbolicExpression,
                      is_SymbolicExpressionRing)

from calculus import maxima as maxima_calculus


from functional import (diff, derivative,
                        laplace, inverse_laplace,
                        expand,
                        integrate, limit, lim,
                        taylor, simplify)

from desolvers import (desolve, desolve_laplace, desolve_system,
                       eulers_method, eulers_method_2x2,
                       eulers_method_2x2_plot)

from var import (var, function, clear_vars)

def symbolic_expression(x):
    """
    Create a symbolic expression from x.

    INPUT:
        x -- an object
    OUTPUT:
        a symbolic expression.

    EXAMPLES:
        sage: a = symbolic_expression(3/2); a
        3/2
        sage: type(a)
        <class 'sage.calculus.calculus.SymbolicConstant'>
        sage: R.<x> = QQ[]; type(x)
        <class 'sage.rings.polynomial.polynomial_element_generic.Polynomial_rational_dense'>
        sage: a = symbolic_expression(2*x^2 + 3); a
        2*x^2 + 3
        sage: type(a)
        <class 'sage.calculus.calculus.SymbolicPolynomial'>
        sage: is_SymbolicExpression(a)
        True
        sage: a in SR
        True
        sage: a.parent()
        Symbolic Ring
    """
    return SR(x)

import desolvers
