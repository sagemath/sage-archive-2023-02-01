"""
Functional notation support for common calculus methods.

EXAMPLES:
We illustrate each of the calculus functional functions.
    sage: simplify(x - x)
    0
    sage: a = var('a')
    sage: derivative(x^a + sin(x), x)
    cos(x) + a*x^(a - 1)
    sage: diff(x^a + sin(x), x)
    cos(x) + a*x^(a - 1)
    sage: derivative(x^a + sin(x), x)
    cos(x) + a*x^(a - 1)
    sage: integral(a*x*sin(x), x)
    a*(sin(x) - x*cos(x))
    sage: integrate(a*x*sin(x), x)
    a*(sin(x) - x*cos(x))
    sage: limit(a*sin(x)/x, x=0)
    a
    sage: taylor(a*sin(x)/x, x, 0, 4)
    a - a*x^2/6 + a*x^4/120
    sage: expand( (x-a)^3 )
    x^3 - 3*a*x^2 + 3*a^2*x - a^3
    sage: laplace( e^(x+a), x, a)
    e^a/(a - 1)
    sage: inverse_laplace( e^a/(a-1), x, a)
    ilt(e^a/(a - 1), x, a)
"""

from calculus import SR, SymbolicExpression, CallableSymbolicExpression

def simplify(f):
    r"""
    Simplify the expression $f$.

    EXAMPLES:
    We simplify the expression $i + x - x$.
        sage: f = I + x - x; simplify(f)
        I

    In fact, printing $f$ yields the same thing -- i.e., the simplified form.
        sage: f
        I

    Nonetheless $f$ and \code{simplify(f)} have a different type; one remembers
    that it is constructed as a sum, and the other is really just the
    simplified expression:
        sage: type(f)
        <class 'sage.calculus.calculus.SymbolicArithmetic'>
        sage: type(simplify(f))
        <class 'sage.calculus.calculus.SymbolicConstant'>
    """
    try:
        return f.simplify()
    except AttributeError:
        return f

def derivative(f, *args, **kwds):
    """
    The derivative of $f$.

    Repeated differentation is supported by the syntax given in the
    examples below.

    ALIAS:
        diff

    EXAMPLES:
    We differentiate a callable symbolic function:
        sage: f(x,y) = x*y + sin(x^2) + e^(-x)
        sage: f
        (x, y) |--> x*y + sin(x^2) + e^(-x)
        sage: derivative(f, x)
        (x, y) |--> y + 2*x*cos(x^2) - e^(-x)
        sage: derivative(f, y)
        (x, y) |--> x

    We differentiate a polynomial:
        sage: t = polygen(QQ, 't')
        sage: f = (1-t)^5; f
        -t^5 + 5*t^4 - 10*t^3 + 10*t^2 - 5*t + 1
        sage: derivative(f)
        -5*t^4 + 20*t^3 - 30*t^2 + 20*t - 5
        sage: derivative(f, t)
        -5*t^4 + 20*t^3 - 30*t^2 + 20*t - 5
        sage: derivative(f, t, t)
        -20*t^3 + 60*t^2 - 60*t + 20
        sage: derivative(f, t, 2)
        -20*t^3 + 60*t^2 - 60*t + 20
        sage: derivative(f, 2)
        -20*t^3 + 60*t^2 - 60*t + 20

    We differentiate a symbolic expression:
        sage: var('a x')
        (a, x)
        sage: f = exp(sin(a - x^2))/x
        sage: derivative(f, x)
        -2*cos(x^2 - a)*e^(-sin(x^2 - a)) - e^(-sin(x^2 - a))/x^2
        sage: derivative(f, a)
        cos(x^2 - a)*e^(-sin(x^2 - a))/x

    Syntax for repeated differentiation:
        sage: R.<u, v> = PolynomialRing(QQ)
        sage: f = u^4*v^5
        sage: derivative(f, u)
        4*u^3*v^5
        sage: f.derivative(u)   # can always use method notation too
        4*u^3*v^5

        sage: derivative(f, u, u)
        12*u^2*v^5
        sage: derivative(f, u, u, u)
        24*u*v^5
        sage: derivative(f, u, 3)
        24*u*v^5

        sage: derivative(f, u, v)
        20*u^3*v^4
        sage: derivative(f, u, 2, v)
        60*u^2*v^4
        sage: derivative(f, u, v, 2)
        80*u^3*v^3
        sage: derivative(f, [u, v, v])
        80*u^3*v^3

    """
    try:
        return f.derivative(*args, **kwds)
    except AttributeError:
        pass
    if not isinstance(f, SymbolicExpression):
        f = SR(f)
    return f.derivative(*args, **kwds)

diff = derivative

def integral(f, *args, **kwds):
    r"""
    The integral of $f$.

    EXAMPLES:
        sage: integral(sin(x), x)
        -cos(x)
        sage: integral(sin(x)^2, x, pi, 123*pi/2)
        121*pi/4
        sage: integral( sin(x), x, 0, pi)
        2

    We integrate a symbolic function:
        sage: f(x,y,z) = x*y/z + sin(z)
        sage: integral(f, z)
        (x, y, z) |--> x*y*log(z) - cos(z)

        sage: var('a,b')
        (a, b)
        sage: assume(b-a>0)
        sage: integral( sin(x), x, a, b)
        cos(a) - cos(b)
        sage: forget()

        sage: print integral(x/(x^3-1), x)
                                         2 x + 1
                       2          arctan(-------)
                  log(x  + x + 1)        sqrt(3)    log(x - 1)
                - --------------- + ------------- + ----------
                         6             sqrt(3)          3

        sage: print integral( exp(-x^2), x )
                               sqrt( pi) erf(x)
                               ----------------
                                      2

    We define the Gaussian, plot and integrate it numerically and symbolically:
        sage: f(x) = 1/(sqrt(2*pi)) * e^(-x^2/2)
        sage: P = plot(f, -4, 4, hue=0.8, thickness=2)
        sage: P.show(ymin=0, ymax=0.4)
        sage: numerical_integral(f, -4, 4)                    # random output
        (0.99993665751633376, 1.1101527003413533e-14)
        sage: integrate(f, x)
        x |--> erf(x/sqrt(2))/2

    You can have \sage calculate multiple integrals.  For example,
    consider the function $exp(y^2)$ on the region between the lines
    $x=y$, $x=1$, and $y=0$. We find the value of the integral on
    this region using the command:
        sage: area = integral(integral(exp(y^2),x,0,y),y,0,1); area
        e/2 - 1/2
        sage: float(area)
        0.85914091422952255

    We compute the line integral of $\sin(x)$ along the arc of the curve
    $x=y^4$ from $(1,-1)$ to $(1,1)$:
        sage: t = var('t')
        sage: (x,y) = (t^4,t)
        sage: (dx,dy) = (diff(x,t), diff(y,t))
        sage: integral(sin(x)*dx, t,-1, 1)
        0
        sage: restore('x,y')   # restore the symbolic variables x and y

    \sage is unable to do anything with the following integral:

        sage: print integral( exp(-x^2)*log(x), x )
                              /      2
                              [   - x
                              I  e     log(x) dx
                              ]
                              /

    \sage does not know how to compute this integral either.
        sage: print integral( exp(-x^2)*ln(x), x, 0, oo)
                              inf
                             /         2
                             [      - x
                             I     e     log(x) dx
                             ]
                             /
                              0

    This definite integral is easy:
        sage: integral( ln(x)/x, x, 1, 2)
        log(2)^2/2

    \sage can't do this elliptic integral (yet):
        sage: integral(1/sqrt(2*t^4 - 3*t^2 - 2), t, 2, 3)
        integrate(1/sqrt(2*t^4 - 3*t^2 - 2), t, 2, 3)

    A double integral:
        sage: y = var('y')
        sage: integral(integral(x*y^2, x, 0, y), y, -2, 2)
        32/5

    This illustrates using assumptions:
        sage: integral(abs(x), x, 0, 5)
        25/2
        sage: integral(abs(x), x, 0, a)
        integrate(abs(x), x, 0, a)
        sage: assume(a>0)
        sage: integral(abs(x), x, 0, a)
        a^2/2
        sage: forget()      # forget the assumptions.

    We integrate and differentiate a huge mess:
        sage: f = (x^2-1+3*(1+x^2)^(1/3))/(1+x^2)^(2/3)*x/(x^2+2)^2
        sage: g = integral(f, x)
        sage: h = f - diff(g, x)

        sage: [float(h(i)) for i in range(5)]     # random low-order bits
        [0.0, -1.1102230246251565e-16, -8.3266726846886741e-17, -4.163336342344337e-17, -6.9388939039072284e-17]
        sage: bool(h == 0)
        True
    """
    try:
        return f.integral(*args, **kwds)
    except ValueError, err:
        raise err
    except AttributeError:
        pass


    if not isinstance(f, SymbolicExpression):
        f = SR(f)
    return f.integral(*args, **kwds)

integrate = integral

def limit(f, dir=None, taylor=False, **argv):
    r"""
    Return the limit as the variable $v$ approaches $a$ from the
    given direction.

        \begin{verbatim}
        limit(expr, x = a)
        limit(expr, x = a, dir='above')
        \end{verbatim}

    INPUT:
        dir -- (default: None); dir may have the value `plus' (or 'above')
               for a limit from above, `minus' (or 'below') for a limit from
               below, or may be omitted (implying a two-sided
               limit is to be computed).
        taylor -- (default: False); if True, use Taylor series, which
               allows more integrals to be computed (but may also crash
               in some obscure cases due to bugs in Maxima).
        **argv -- 1 named parameter

    ALIAS: You can also use lim instead of limit.

    EXAMPLES:
        sage: limit(sin(x)/x, x=0)
        1
        sage: limit(exp(x), x=oo)
        +Infinity
        sage: lim(exp(x), x=-oo)
        0
        sage: lim(1/x, x=0)
        und
        sage: limit(sqrt(x^2+x+1)+x, taylor=True, x=-oo)
        -1/2
        sage: limit((tan(sin(x)) - sin(tan(x)))/x^7, taylor=True, x=0)
        1/30

    \sage does not know how to do this limit (which is 0),
    so it returns it unevaluated:
        sage: lim(exp(x^2)*(1-erf(x)), x=infinity)
        limit(e^x^2 - e^x^2*erf(x), x, +Infinity)
    """
    if not isinstance(f, SymbolicExpression):
        f = SR(f)
    return f.limit(dir=dir, taylor=taylor, **argv)

lim = limit

def taylor(f, v, a, n):
    """
    Expands self in a truncated Taylor or Laurent series in the
    variable $v$ around the point $a$, containing terms through $(x - a)^n$.

    INPUT:
        v -- variable
        a -- number
        n -- integer

    EXAMPLES:
        sage: var('x,k,n')
        (x, k, n)
        sage: taylor (sqrt (1 - k^2*sin(x)^2), x, 0, 6)
        1 - k^2*x^2/2 - (3*k^4 - 4*k^2)*x^4/24 - (45*k^6 - 60*k^4 + 16*k^2)*x^6/720
        sage: taylor ((x + 1)^n, x, 0, 4)
        1 + n*x + (n^2 - n)*x^2/2 + (n^3 - 3*n^2 + 2*n)*x^3/6 + (n^4 - 6*n^3 + 11*n^2 - 6*n)*x^4/24

    """
    if not isinstance(f, SymbolicExpression):
        f = SR(f)
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


def laplace(f, t, s):
    r"""
    Attempts to compute and return the Laplace transform of \code{self}
    with respect to the variable $t$ and transform parameter $s$.  If
    this function cannot find a solution, a delayed function is returned.

    The function that is returned may be be viewed as a function of $s$.

    EXAMPLES:
        sage: var('a,s,t')
        (a, s, t)
        sage: f = exp (2*t + a) * sin(t) * t; f
        t*e^(2*t + a)*sin(t)
        sage: L = laplace(f, t, s); L
        e^a*(2*s - 4)/(s^2 - 4*s + 5)^2
        sage: inverse_laplace(L, s, t)
        t*e^(2*t + a)*sin(t)

    Unable to compute solution:
        sage: laplace(1/s, s, t)
        laplace(1/s, s, t)
    """
    if not isinstance(f, SymbolicExpression):
        f = SR(f)
    return f.laplace(t,s)

def inverse_laplace(f, t, s):
    r"""
    Attempts to compute and return the inverse Laplace transform of
    \code{self} with respect to the variable $t$ and transform parameter $s$.  If
    this function cannot find a solution, a delayed function is returned, which
    is called \code{ilt}.

    EXAMPLES:
        sage: f(t) = t*cos(t)
        sage: s = var('s')
        sage: L = laplace(f, t, s); L
        t |--> 2*s^2/(s^2 + 1)^2 - 1/(s^2 + 1)
        sage: inverse_laplace(L, s, t)
        t |--> t*cos(t)
        sage: print inverse_laplace(1/(s^3+1), s, t)
                           sqrt(3) t        sqrt(3) t
                       sin(---------)   cos(---------)      - t
                  t/2          2                2          e
                 e    (-------------- - --------------) + -----
                          sqrt(3)             3             3

    No explicit inverse Laplace transform, so one is returned formally as a function \code{ilt}.
        sage: inverse_laplace(cos(s), s, t)
        ilt(cos(s), s, t)
    """
    if not isinstance(f, SymbolicExpression):
        f = SR(f)
    return f.inverse_laplace(t,s)


