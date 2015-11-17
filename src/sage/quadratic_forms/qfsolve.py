"""
Solving quadratic equations.

Interface to the PARI/GP quadratic forms code of Denis Simon.

AUTHORS:

- Denis Simon (GP code)

- Nick Alexander (Sage interface)

- Jeroen Demeyer (2014-09-23): use PARI instead of GP scripts,
  return vectors instead of tuples (:trac:`16997`).

- Tyler Gaona (2015-11-14): added the `solve` method
"""

#*****************************************************************************
#       Copyright (C) 2008 Nick Alexander
#       Copyright (C) 2014 Jeroen Demeyer
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.all import ZZ, QQ
from sage.modules.free_module_element import vector

def qfsolve(G):
    r"""
    Find a solution `x = (x_0,...,x_n)` to `x G x^t = 0` for an
    `n \times n`-matrix ``G`` over `\QQ`.

    OUTPUT:

    If a solution exists, return a vector of rational numbers `x`.
    Otherwise, returns `-1` if no solution exists over the reals or a
    prime `p` if no solution exists over the `p`-adic field `\QQ_p`.

    ALGORITHM:

    Uses PARI/GP function ``qfsolve``.

    EXAMPLES::

        sage: from sage.quadratic_forms.qfsolve import qfsolve
        sage: M = Matrix(QQ, [[0, 0, -12], [0, -12, 0], [-12, 0, -1]]); M
        [  0   0 -12]
        [  0 -12   0]
        [-12   0  -1]
        sage: sol = qfsolve(M); sol
        (1, 0, 0)
        sage: sol.parent()
        Vector space of dimension 3 over Rational Field

        sage: M = Matrix(QQ, [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        sage: ret = qfsolve(M); ret
        -1
        sage: ret.parent()
        Integer Ring

        sage: M = Matrix(QQ, [[1, 0, 0], [0, 1, 0], [0, 0, -7]])
        sage: qfsolve(M)
        7

        sage: M = Matrix(QQ, [[3, 0, 0, 0], [0, 5, 0, 0], [0, 0, -7, 0], [0, 0, 0, -11]])
        sage: qfsolve(M)
        (3, -4, -3, -2)
    """
    ret = G._pari_().qfsolve()
    if ret.type() == 't_COL':
        return vector(QQ, ret)
    return ZZ(ret)

def qfparam(G, sol):
    r"""
    Parametrizes the conic defined by the matrix ``G``.

    INPUT:

     - ``G`` -- a `3 \times 3`-matrix over `\QQ`.

     - ``sol`` -- a triple of rational numbers providing a solution
       to sol*G*sol^t = 0.

    OUTPUT:

    A triple of polynomials that parametrizes all solutions of
    x*G*x^t = 0 up to scaling.

    ALGORITHM:

    Uses PARI/GP function ``qfparam``.

    EXAMPLES::

        sage: from sage.quadratic_forms.qfsolve import qfsolve, qfparam
        sage: M = Matrix(QQ, [[0, 0, -12], [0, -12, 0], [-12, 0, -1]]); M
        [  0   0 -12]
        [  0 -12   0]
        [-12   0  -1]
        sage: sol = qfsolve(M)
        sage: ret = qfparam(M, sol); ret
        (-12*t^2 - 1, 24*t, 24)
        sage: ret.parent()
        Ambient free module of rank 3 over the principal ideal domain Univariate Polynomial Ring in t over Rational Field
    """
    R = QQ['t']
    t = R.gen()
    mat = G._pari_().qfparam(sol)
    # Interpret the rows of mat as coefficients of polynomials
    return vector(R, mat.Col())

def solve(self, c=0):
    r""" 
    Return a vector x such that self(x) == c.

    INPUT:
    
    - ``self`` -- a quadratic form.
    - ``c`` -- (default: 0) a rational number. 

    OUTPUT:
    
    - A vector x satisfying self(x) == c.

    ALGORITHM:
    
    Uses PARI's qfsolve(). Algorithm described by Jeroen Demeyer; see comments on :trac:`19112`

    EXAMPLES::
    
        sage: Q = QuadraticForm(QQ, 2, [17, 94, 130])
        sage: x = Q.solve(5); x
        (17, -6)
        sage: Q(x)
        5
    
        sage: Q.solve(6)
        Traceback (most recent call last):
        ...
        ArithmeticError: no solution found (local obstruction at 3)
    
        sage: G = DiagonalQuadraticForm(QQ, [5, -3, -2])
        sage: x = G.solve(10); x
        (3/2, -1/2, 1/2)
        sage: G(x)
        10

        sage: F = DiagonalQuadraticForm(QQ, [1, -4])
        sage: x = F.solve(); x
        (2, 1)
        sage: F(x)
        0
    """
    from sage.quadratic_forms.quadratic_form import QuadraticForm
    from sage.rings.integer import Integer
    
    # If no argument passed for c, we just pass self into qfsolve().
    if c == 0:
        x = qfsolve(self.Gram_matrix())
        if isinstance(x, Integer):
            raise ArithmeticError("no solution found (local obstruction at {})".format(x))
        return x

    # If c != 0, define a new quadratic form Q = self - c*z^2
    M = self.Gram_matrix()
    q_terms = []
    for i in range(self.dim()):
        for j in range(i, self.dim()):
            if i == j:
                q_terms.append(M[i][j])
            else:
                q_terms.append(M[i][j] * 2)
        q_terms.append(0)
    q_terms.append(-c)

    Q = QuadraticForm(self.base_ring(), self.dim() + 1, q_terms)

    # Find a solution x to Q(x) = 0, using qfsolve()
    x = qfsolve(Q.Gram_matrix())
    # Raise an error if qfsolve() doesn't find a solution
    if isinstance(x, Integer):
        raise ArithmeticError("no solution found (local obstruction at {})".format(x))
    
    # Let z be the last term of x, and remove z from x
    z = x[-1]
    x = x[:-1]
    # If z != 0, then Q(x/z) = c
    if z:
        return x * (1/z)

    # Case 2: We found a solution self(x) = 0. Let e be any vector such 
    # that B(x,e) != 0, where B is the bilinear form corresponding to self.
    # To find e, just try all unit vectors (0,..0,1,0...0). 
    # Let a = (c - self(e))/2*B(x,e) and let y = e + a*x.
    # Then self(y) = B(e + a*x, e + a*x) = self(e) + 2B(e, a*x) 
    #              = self(e) + 2([c - self(e)]/[2B(x,e)]) * B(x,e) = c.
    e = vector([0] * self.dim())
    e[0] = 1
    i = 1
    while self.bilinear_map(x, e) == 0:
        e[i-1] = 0
        e[i] = 1

    a = (c - self(e)) / (2 * self.bilinear_map(x, e))
    return e + a*x
