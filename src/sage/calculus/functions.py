r"""
Calculus functions.
"""

from sage.matrix.all import matrix
from calculus import SymbolicVariable
from functional import diff

def wronskian(*args):
    """Returns the Wronskian of the provided functions, differentiating with
    respect to the given variable. If no variable is provided,
    diff(f) is called for each function f.

    wronskian(f1,...,fn, x) returns the Wronskian of f1,...,fn, with
    derivatives taken with respect to x.

    wronskian(f1,...,fn) returns the Wronskian of f1,...,fn where
    k'th derivatives are computed by doing `.derivative(k)' on each
    function.

    The Wronskian of a list of functions is a determinant of derivatives.
    The nth row (starting from 0) is a list of the nth derivatives of the
    given functions.

    For two functions:

                              | f   g  |
                 W(f, g) = det|        | = f*g' - g*f'.
                              | f'  g' |

    EXAMPLES:
        sage: wronskian(e^x, x^2)
        2*x*e^x - x^2*e^x

        sage: var('x, y'); wronskian(x*y, log(x), x)
        (x, y)
        y - log(x)*y

      If your functions are in a list, you can use `*' to turn them into
      arguments to wronskian():
        sage: wronskian(*[x^k for k in range(1, 5)])
        12*x^4

      If you want to use 'x' as one of the functions in the Wronskian,
      you can't put it last or it will be interpreted as the variable
      with respect to which we differentiate. There are several ways to
      get around this.

      Two-by-two Wronskian of sin(x) and e^x:
        sage: wronskian(sin(x), e^x, x)
        e^x*sin(x) - e^x*cos(x)

      Three-by-three Wronskian of sin(x), e^x, and x:
        sage: wronskian(sin(x), cos(x), x+0)
        x*(-sin(x)^2 - cos(x)^2)

      Or don't put x last:
        sage: wronskian(x, sin(x), e^x)
        x*(e^x*sin(x) + e^x*cos(x)) - 2*e^x*sin(x)

      Example where one of the functions is constant:
        sage: wronskian(1, e^(-x), e^(2*x))
        -6*e^x

    NOTES:
        http://en.wikipedia.org/wiki/Wronskian
        http://planetmath.org/encyclopedia/WronskianDeterminant.html

    AUTHORS:
        - Dan Drake (2008-03-12)
    """
    if len(args) == 0:
        raise TypeError, 'wronskian() takes at least one argument (0 given)'
    elif len(args) == 1:
        # a 1x1 Wronskian is just its argument
        return args[0]
    else:
        if isinstance(args[-1], SymbolicVariable):
            # if last argument is a variable, peel it off and
            # differentiate the other args
            v = args[-1]
            fs = args[0:-1]
            row = lambda n: map(lambda f: diff(f, v, n), fs)
        else:
            # if the last argument isn't a variable, just run
            # .derivative on everything
            fs = args
            row = lambda n: map(lambda f: diff(f, n), fs)
        # NOTE: I rewrote the below as two lines to avoid a possible subtle
        # memory management problem on some platforms (only vmware as far
        # as we know?).  See trac #2990.
        # There may still be a real problem that this is just hiding for now.
        A = matrix(map(row, range(len(fs))))
        return A.determinant()
        #return matrix(map(row, range(len(fs)))).determinant()
