from equations import SymbolicEquation, forget, assume, assumptions, solve

from calculus import (SymbolicExpressionRing,
                      is_SymbolicExpressionRing,
                      is_SymbolicExpression,
                      CallableSymbolicExpressionRing,
                      is_CallableSymbolicExpressionRing,
                      is_CallableSymbolicExpression,
                      SR,
                      sin, cos, sec, tan, log, erf, sqrt, asin, acos, atan,
                      tanh, sinh, cosh, coth, sech, csch, ln,
                      ceil, floor,
                      polylog,
                      abs_symbolic, exp,
                      is_SymbolicExpression,
                      is_SymbolicExpressionRing)


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
    return SR(x)

import desolvers
