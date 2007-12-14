"""
Doctests for cvxopt, scipy

sage: from cvxopt.base import *
sage: from cvxopt import umfpack
sage: from scipy import optimize
sage: from scipy import special
sage: from scipy import integrate
sage: from scipy import linsolve
sage: from scipy import interpolate
sage: from scipy import sparse
sage: import arpack

#Test arpack
sage: import scipy
sage: n=scipy.zeros((3,100))
sage: n[0,:]=-1
sage: n[1,:]=1
sage: n[2,:]=-2
sage: A=sparse.spdiags(n,[-1,0,1],int(100),int(100))
sage: e,v=arpack.eigen(A,3)
sage: e
array([ 3.8270...+0.j,  3.8229...+0.j,  0.        +0.j,  0.
+0.j])


"""

