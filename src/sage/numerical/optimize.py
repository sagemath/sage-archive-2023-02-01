"""
Numerical Root Finding and Optimization

AUTHOR:
    -- William Stein (2007)
"""

def find_root(f, a, b, xtol=10e-13, rtol=4.5e-16, maxiter=100, full_output=False):
    """
    Numerically find a root of f on the closed interval [a,b]
    (or [b,a]) if possible, where f is a function in the one variable.

    INPUT:
        f -- a function of one variable or symbolic equality
        a, b -- endpoints of the interval
        xtol, rtol -- the routine converges when a root is known
                to lie within xtol of the value return. Should be
                >= 0.  The routine modifies this to take into
                account the relative precision of doubles.
        maxiter -- integer; if convergence is not achieved in
                maxiter iterations, an error is raised. Must be >= 0.
        full_output -- bool (default: False), if True, also return
                object that contains information about convergence.

    EXAMPLES:
    An example involving an algebraic polynomial function.
        sage: R.<x> = QQ[]
        sage: f = (x+17)*(x-3)*(x-1/8)^3
        sage: find_root(f, 0,4)
        2.9999999999999951
        sage: find_root(f, 0,1)  # note -- precision of answer isn't very good on some machines.
        0.124999...
        sage: find_root(f, -20,-10)
        -17.0

    In Pomerance book on primes he asserts that the famous Riemann
    Hypothesis is equivalent to the statement that the function f(x)
    defined below is positive for all $x\geq 2.01$
        sage: def f(x):
        ...       return sqrt(x) * log(x) - abs(Li(x) - prime_pi(x))

    We find where $f$ equals, i.e., what value that is slightly smaller
    than $2.01$ that could have been used in the formulation of the Riemann
    Hypothesis:
        sage: find_root(f, 2, 4, rtol=0.0001)
        2.0082590205656166

    This agrees with the plot:
        sage: show(plot(f,2,2.01),xmin=2,xmax=2.01, ymin=0.01,ymax=0.01)
    """
    try:
        return f.find_root(a=a,b=b,xtol=xtol,rtol=rtol,maxiter=maxiter,full_output=full_output)
    except AttributeError:
        pass
    a = float(a); b = float(b)
    if a > b:
        a, b = b, a
    left = f(a)
    right = f(b)
    if left > 0 and right > 0:
        # Refine further -- try to find a point where this
        # function is negative in the interval
        val, s = find_minimum_on_interval(f, a, b)
        if val > 0:
            if val < rtol:
                if full_output:
                    return s, "No extra data"
                else:
                    return s
            raise RuntimeError, "f appears to have no zero on the interval"
        # If we found such an s, then we just instead find
        # a root between left and s or s and right.
        a = s   # arbitrary choice -- maybe should try both and take one that works?

    elif left < 0 and right < 0:
        # Refine further
        val, s = find_maximum_on_interval(f, a, b)
        if val < 0:
            if abs(val) < rtol:
                if full_output:
                    return s, "No extra data"
                else:
                    return s
            raise RuntimeError, "f appears to have no zero on the interval"
        a = s

    import scipy.optimize
    return scipy.optimize.brentq(f, a, b,
                                 full_output=full_output, xtol=xtol, rtol=rtol, maxiter=maxiter)

def find_maximum_on_interval(f, a, b, tol=1.48e-08, maxfun=500):
    """
    Numerically find the maximum of the expression f on the interval
    [a,b] (or [b,a]) along with the point at which the maximum is attained.

    See the documentation for \code{find_minimum_on_interval}
    for more details.

    EXAMPLES:
        sage: f = lambda x: x*cos(x)
        sage: find_maximum_on_interval(f, 0,5)
        (0.561096338191, 0.8603335890...)
        sage: find_maximum_on_interval(f, 0,5, tol=0.1, maxfun=10)
        (0.561090323458, 0.857926501456)
    """
    def g(z):
        r"""
        Returns the negative of the input function f. Finding the maximum
        of f(z) on [a,b] is equivalent to finding th minimum of -f(z) on
        [a,b].

        EXAMPLES:

        """
        return -f(z)
    minval, x = find_minimum_on_interval(g, a=a, b=b, tol=tol, maxfun=maxfun)
    return -minval, x

def find_minimum_on_interval(f, a, b, tol=1.48e-08, maxfun=500):
    """
    Numerically find the minimum of the expression self on the
    interval [a,b] (or [b,a]) and the point at which it attains that
    minimum.  Note that self must be a function of (at most) one
    variable.

    INPUT:
        a,b -- endpoints of interval on which to minimize self.
        tol -- the convergence tolerance
        maxfun -- maximum function evaluations

    OUTPUT:
        minval -- (float) the minimum value that self takes on in the interval [a,b]
        x -- (float) the point at which self takes on the minimum value

    EXAMPLES:
        sage: f = lambda x: x*cos(x)
        sage: find_minimum_on_interval(f, 1, 5)
        (-3.28837139559, 3.42561846957)
        sage: find_minimum_on_interval(f, 1, 5, tol=1e-3)
        (-3.28837136189098, 3.42575079030572)
        sage: find_minimum_on_interval(f, 1, 5, tol=1e-2, maxfun=10)
        (-3.28837084598, 3.42508402203)
        sage: show(plot(f, 0, 20))
        sage: find_minimum_on_interval(f, 1, 15)
        (-9.47729425948, 9.52933441095)

    ALGORITHM: Uses scipy.optimize.fminbound which uses Brent's method.

    AUTHOR:
         -- William Stein (2007-12-07)
    """
    try:
        return f.find_minimum_on_interval(a=a, b=b, tol=tol,maxfun=maxfun)
    except AttributeError:
        pass
    a = float(a); b = float(b)
    import scipy.optimize
    xmin, fval, iter, funcalls = scipy.optimize.fminbound(f, a, b, full_output=1, xtol=tol, maxfun=maxfun)
    return fval, xmin
