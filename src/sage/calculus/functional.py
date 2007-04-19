"""
Functional notation support for common calculus methods.
"""

from calculus import SER, SymbolicExpression, CallableFunction

def simplify(f):
    """
    Simplify the expression f.
    """
    try:
        return f.simplify()
    except AttributeError:
        return f

def diff(f, *args, **kwds):
    """
    Formally differentiate the function f.
    """
    if isinstance(f, CallableFunction):
        return f.derivative(*args, **kwds)
    if not isinstance(f, SymbolicExpression):
        f = SER(f)
    return f.derivative(*args, **kwds)

derivative = diff

def integrate(f, *args, **kwds):
    """
    Returns the integral of f.

    EXAMPLES:
        sage: integrate(sin(x), x)
        -cos(x)
        sage: integrate(sin(x)^2, x, pi, 123*pi/2)
        121*pi/4
    """
    if isinstance(f, CallableFunction):
        return f.derivative(*args, **kwds)
    if not isinstance(f, SymbolicExpression):
        f = SER(f)
    return f.integral(*args, **kwds)

def limit(f, v, a, dir=None):
    """
    Return the limit as the variable v approaches a from the
    given direction.

    INPUT:
        f -- symbolic expression
        v -- variable
        a -- number
        dir -- (default: None); dir may have the value `plus' (or 'above')
               for a limit from above, `minus' (or 'below') for a limit from
               below, or may be omitted (implying a two-sided
               limit is to be computed).
    """
    if not isinstance(f, SymbolicExpression):
        f = SER(f)
    return f.limit(v=v,a=a,dir=dir)

lim = limit

def taylor(f, v, a, n):
    """
    Expands self in a truncated Taylor or Laurent series in the
    variable v around the point a, containing terms through $(x - a)^n$.

    INPUT:
        v -- variable
        a -- number
        n -- integer

    EXAMPLES:
        sage: taylor (sqrt (1 - k^2*sin(x)^2), x, 0, 6)
        1 - (k^2*x^2/2) - ((3*k^4 - 4*k^2)*x^4/24) - ((45*k^6 - 60*k^4 + 16*k^2)*x^6/720)
        sage: taylor ((x + 1)^n, x, 0, 4)
        1 + n*x + (n^2 - n)*x^2/2 + (n^3 - 3*n^2 + 2*n)*x^3/6 + (n^4 - 6*n^3 + 11*n^2 - 6*n)*x^4/24

    """
    if not isinstance(f, SymbolicExpression):
        f = SER(f)
    return f.taylor(v=v,a=a,n=n)


def expand(x, *args, **kwds):
    """
    EXAMPLES:
        sage: a = (1+I)*(2-sqrt(3)*I); a
        (I + 1)*(2 - sqrt(3)*I)
        sage: expand(a)
        -sqrt(3)*I + 2*I + sqrt(3) + 2
        sage: a = (x-1)*(x^2 - 1); a
        (x - 1)*(x^2 - 1)
        sage: expand(a)
        x^3 - x^2 - x + 1

    You can also use expand on polynomial, integer, and other
    factorizations:
        sage: x = polygen(ZZ)
        sage: F = factor(x^12 - 1); F
        (x - 1) * (x + 1) * (x^2 - x + 1) * (x^2 + 1) * (x^2 + x + 1) * (x^4 - x^2 + 1)
        sage: expand(F)
        x^12 - 1
        sage: F.expand()
        x^12 - 1
        sage: F = factor(2007); F
        3^2 * 223
        sage: expand(F)
        2007

    Note: If you want to compute the expanded form of a polynomial
    arithmetic operation quickly and the coefficients of the polynomial
    all lie in some ring, e.g., the integers, it is vastly faster to
    create a polynomial ring and do the arithmetic there.

        sage: x = polygen(ZZ)      # polynomial over a given base ring.
        sage: f = sum(x^n for n in range(5))
        sage: f*f                  # much faster, even if the degree is huge
        x^8 + 2*x^7 + 3*x^6 + 4*x^5 + 5*x^4 + 4*x^3 + 3*x^2 + 2*x + 1

    """
    try:
        return x.expand(*args, **kwds)
    except AttributeError:
        return x
