"""
TESTS
Compute the Christoffel symbol.

    sage: var('r theta phi')
    (r, theta, phi)
    sage: m = matrix(SER, [[(1-1/r),0,0,0],[0,-(1-1/r)^(-1),0,0],[0,0,-r^2,0],[0,0,0,-r^2*(sin(theta))^2]])
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

"""
