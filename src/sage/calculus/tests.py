"""
TESTS
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

Restoring variables after they have been turned into functions:
    sage: x = function('x')
    sage: sin(x).variables()
    ()
    sage: restore('x')
    sage: sin(x).variables()
    (x,)
"""
