"""
Further Examples.

Compute the Christoffel symbol.

    sage: var('r theta phi')
    (r, theta, phi)
    sage: m = matrix(SR, [[(1-1/r),0,0,0],[0,-(1-1/r)^(-1),0,0],[0,0,-r^2,0],[0,0,0,-r^2*(sin(theta))^2]])
    sage: print m
    [        1 - (1/r)                 0                 0                 0]
    [                0    -1/(1 - (1/r))                 0                 0]
    [                0                 0              -r^2                 0]
    [                0                 0                 0 -r^2*sin(theta)^2]

    sage: def christoffel(i,j,k,vars,g):
    ...   s = 0
    ...   ginv = g^(-1)
    ...   for l in range(g.nrows()):
    ...      s = s + (1/2)*ginv[k,l]*(g[j,l].diff(vars[i])+g[i,l].diff(vars[j])-g[i,j].diff(vars[l]))
    ...   return s

    sage: christoffel(3,3,2, [t,r,theta,phi], m)
    -cos(theta)*sin(theta)
    sage: X = christoffel(1,1,1,[t,r,theta,phi],m)
    sage: print X
                                     1
                                     - - 1
                                     r
                                 -------------
                                        1 2  2
                                 2 (1 - -)  r
                                        r
    sage: print X.rational_simplify()
                                       1
                                 - ----------
                                      2
                                   2 r  - 2 r

Some basic things:

    sage: f(x,y) = x^3 + sinh(1/y)
    sage: f
    (x, y) |--> sinh(1/y) + x^3
    sage: f^3
    (x, y) |--> (sinh(1/y) + x^3)^3
    sage: (f^3).expand()
    (x, y) |--> sinh(1/y)^3 + 3*x^3*sinh(1/y)^2 + 3*x^6*sinh(1/y) + x^9

A polynomial over a symbolic base ring:
    sage: R = SR[x]
    sage: f = R([1/sqrt(2), 1/(4*sqrt(2))])
    sage: f
    1/(4*sqrt(2))*x + 1/sqrt(2)
    sage: -f
    (-1/(4*sqrt(2)))*x + -1/sqrt(2)
    sage: (-f).degree()
    1

Something that was a printing bug.  This tests that we print
the simplified version using ASCII art:
    sage: A = exp(I*pi/5)
    sage: print A*A*A*A*A*A*A*A*A*A
                                           1

We check a statement made at the beginning of Friedlander and Joshi's book
on Distributions:
    sage: f = sin(x^2)
    sage: g = cos(x) + x^3
    sage: u = f(x+t) + g(x-t)
    sage: u
    sin((x + t)^2) + cos(x - t) + (x - t)^3
    sage: u.diff(t,2) - u.diff(x,2)
    0

Restoring variables after they have been turned into functions:
    sage: x = function('x')
    sage: sin(x).variables()
    ()
    sage: restore('x')
    sage: sin(x).variables()
    (x,)

MATHEMATICA:
Some examples of integration and differentiation taken from some
Mathematica docs:
    sage: diff(x^n, x)
    n*x^(n - 1)
    sage: diff(x^2 * log(x+a), x)
    2*x*log(x + a) + x^2/(x + a)
    sage: derivative(atan(x), x)
    1/(x^2 + 1)
    sage: derivative(x^n, x, 3)
    (n - 2)*(n - 1)*n*x^(n - 3)
    sage: derivative( function('f')(x), x)
    diff(f(x), x, 1)
    sage: diff( 2*x*f(x^2), x)
    2*x*diff(f(x^2), x, 1) + 2*f(x^2)
    sage: print integrate( 1/(x^4 - a^4), x)
                                                           x
                                                      atan(-)
                            log(x + a)   log(x - a)	       a
                          - ---------- + ---------- - -------
                                  3		   3	      3
                               4 a	        4 a	   2 a
    sage: expand(integrate(log(1-x^2), x))
    x*log(1 - x^2) + log(x + 1) - log(x - 1) - 2*x
    sage: integrate(log(1-x^2)/x, x)
    log(x)*log(1 - x^2) + polylog(2, 1 - x^2)/2
    sage: integrate(exp(1-x^2),x)
    sqrt(pi)*e*erf(x)/2
    sage: integrate(sin(x^2),x)
    sqrt(pi)*((sqrt(2)*I + sqrt(2))*erf((sqrt(2)*I + sqrt(2))*x/2) + (sqrt(2)*I - sqrt(2))*erf((sqrt(2)*I - sqrt(2))*x/2))/8
    sage: integrate((1-x^2)^n,x)
    integrate((1 - x^2)^n, x)
    sage: integrate(x^x,x)
    integrate(x^x, x)
    sage: print integrate(1/(x^3+1),x)
                                             2 x - 1
                           2	    atan(-------)
                      log(x  - x + 1)	 sqrt(3)    log(x + 1)
                    - --------------- + ------------- + ----------
                             6	       sqrt(3)	        3
    sage: integrate(1/(x^3+1), x, 0, 1)
    (6*log(2) + sqrt(3)*pi)/18 + sqrt(3)*pi/18

    sage: forget(); assume(c > 0)
    sage: integrate(exp(-c*x^2), x, -oo, oo)
    sqrt(pi)/sqrt(c)
    sage: forget()

The following are a bunch of examples of integrals that Mathematica
can do, but SAGE currently can't do:
    sage: integrate(sqrt(x + sqrt(x)), x)    # todo -- mathematica can do this
    integrate(sqrt(x + sqrt(x)), x)
    sage: integrate(log(x)*exp(-x^2))        # todo -- mathematica can do this
    integrate(e^(-x^2)*log(x), x)

    sage.: # Todo -- Mathematica can do this and gets pi^2/15
    sage.: # integrate(log(1+sqrt(1+4*x)/2)/x, x, 0, 1)
    [boom!]
    Integral is divergent

    sage: integrate(ceil(x^2 + floor(x)), x, 0, 5)    # todo: mathematica can do this
    integrate(ceil(x^2) + floor(x), x, 0, 5)


MAPLE:
The basic differentiation and integration examples in the
Maple documentation:

    sage: diff(sin(x), x)
    cos(x)
    sage: diff(sin(x), y)
    0
    sage: diff(sin(x), x, 3)
    -cos(x)
    sage: diff(x*sin(cos(x)), x)
    sin(cos(x)) - x*sin(x)*cos(cos(x))
    sage: diff(tan(x), x)
    sec(x)^2
    sage: f = function('f'); f
    f
    sage: diff(f(x), x)
    diff(f(x), x, 1)
    sage: diff(f(x,y), x, y)
    diff(f(x, y), x, 1, y, 1)
    sage: diff(f(x,y), x, y) - diff(f(x,y), y, x)
    0
    sage: g = function('g')
    sage: diff(g(x,y,z), x,z,z)
    diff(g(x, y, z), x, 1, z, 2)
    sage: integrate(sin(x), x)
    -cos(x)
    sage: integrate(sin(x), x, 0, pi)
    2
    sage: assume(b-a>0)      # annoying -- maple doesn't require this...
    sage: print integrate(sin(x), x, a, b)
                                    cos(a) - cos(b)
    sage: forget()


    sage: integrate( x/(x^3-1), x)
    (-log(x^2 + x + 1))/6 + atan((2*x + 1)/sqrt(3))/sqrt(3) + log(x - 1)/3
    sage: integrate(exp(-x^2), x)
    sqrt(pi)*erf(x)/2
    sage: integrate(exp(-x^2)*log(x), x)       # todo: maple can compute this exactly.
    integrate(e^(-x^2)*log(x), x)
    sage: f = exp(-x^2)*log(x)
    sage: f.nintegral(x, 0, 999)
    (-0.87005772672831549, 7.5584116743243612e-10, 567, 0)
    sage: integral(1/sqrt(2*t^4 - 3*t^2 - 2), t, 2, 3)     # todo: maple can do this
    integrate(1/(sqrt(2*t^4 - 3*t^2 - 2)), t, 2, 3)
    sage: integral(integral(x*y^2, x, 0, y), y, -2, 2)
    32/5

We verify several standard differentiation rules:
    sage: function('f, g')
    (f, g)
    sage: diff(f(t)*g(t),t)
    f(t)*diff(g(t), t, 1) + g(t)*diff(f(t), t, 1)
    sage: diff(f(t)/g(t), t)
    diff(f(t), t, 1)/g(t) - (f(t)*diff(g(t), t, 1)/g(t)^2)
    sage: diff(f(t) + g(t), t)
    diff(g(t), t, 1) + diff(f(t), t, 1)
    sage: diff(c*f(t), t)
    c*diff(f(t), t, 1)
"""
