"""
Doctests for cvxopt, scipy

sage: from cvxopt.base import *
sage: from cvxopt import umfpack
sage: from scipy import optimize
sage: from scipy import special
sage: from scipy import integrate
sage: from scipy.sparse.linalg import dsolve
sage: from scipy import interpolate
sage: from scipy import sparse
sage: from scipy.sparse.linalg.eigen.arpack import arpack

#Test arpack
#This matrix is the finite difference approximation to
# the eigenvalue problem
#d^2f/dx^2=\lambda f, on [0,\pi], which boundary values 0
# The lowest eigenvalue calulated should be close to 1
sage: import scipy
sage: from scipy.sparse.linalg.interface import aslinearoperator
sage: n=scipy.zeros((3,500))
sage: n[0,:]=-1
sage: n[1,:]=2
sage: n[2,:]=-1
sage: A=sparse.spdiags(n,[-1,0,1],int(500),int(500))
sage: A=aslinearoperator(A)
sage: e,v=arpack.eigs(A,k=6,which='SM')  # long time (4s on sage.math, 2012)
sage: e[0]*float(501/pi)**2  # long time
(0.999...+0j)
"""

