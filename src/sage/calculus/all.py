## from calculus import (SymbolicExpressionRing,
##                       is_SymbolicExpressionRing,
##                       is_SymbolicExpression,
##                       is_SymbolicVariable,
##                       CallableSymbolicExpressionRing,
##                       is_CallableSymbolicExpressionRing,
##                       is_CallableSymbolicExpression,
##                       is_SymbolicExpression,
##                       is_SymbolicExpressionRing)

from calculus import maxima as maxima_calculus
from calculus import (laplace, inverse_laplace,
                      limit, lim, clear_functions)

from functional import (diff, derivative,
                        expand,
                        taylor, simplify)

from functions import (wronskian,jacobian)

from desolvers import (desolve, desolve_laplace, desolve_system,
                       eulers_method, eulers_method_2x2,
                       eulers_method_2x2_plot)

from var import (var, function, clear_vars)

def symbolic_expression(x):
    """
    Create a symbolic expression from x.

    INPUT:

    - ``x`` - an object

    OUTPUT:

    - a symbolic expression.

    EXAMPLES::

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
        sage: from sage.calculus.calculus import is_SymbolicExpression
        sage: is_SymbolicExpression(a)
        True
        sage: a in SR
        True
        sage: a.parent()
        Symbolic Ring

    This may not always return an element of SR.

    ::

        sage: E = EllipticCurve('15a'); E
        Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 - 10*x - 10 over Rational Field
        sage: symbolic_expression(E)
        y^2 + x*y + y == x^3 + x^2 - 10*x - 10
        sage: symbolic_expression(E) in SR
        False

    """
    from sage.symbolic.expression import Expression
    from sage.symbolic.ring import SR
    if isinstance(x, Expression):
        return x
    elif hasattr(x, '_symbolic_'):
        return x._symbolic_(SR)
    else:
        return SR(x)

import desolvers
