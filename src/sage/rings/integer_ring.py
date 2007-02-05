r"""
\protect{Ring $\Z$ of Integers}

The class \\class{IntegerRing} represents
the ring $\\mathbf{Z}$ of (arbitrary precision) integers.  Each integer
is an instance of the class \\class{Integer}, which is defined
in a Pyrex extension module that wraps GMP integers
(the \\class{mpz_t} type in GMP).

    sage: Z = IntegerRing(); Z
    Integer Ring
    sage: Z.characteristic()
    0
    sage: Z.is_field()
    False

There is a unique instances of class \\class{IntegerRing}.  To create
an \\class{Integer}, coerce either a Python int, long, or a string.
Various other types will also coerce to the integers, when it makes
sense.

    sage: a = Z(1234); b = Z(5678); print a, b
    1234 5678
    sage: type(a)
    <type 'sage.rings.integer.Integer'>
    sage: a + b
    6912
    sage: Z('94803849083985934859834583945394')
    94803849083985934859834583945394
"""

#*****************************************************************************
#
#   SAGE: System for Algebra and Geometry Experimentation
#
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from integer import ZZ, is_IntegerRing

def IntegerRing():
    return ZZ

def factor(n, algorithm='pari'):
    """
    Return the factorization of the positive integer $n$ as a list of
    tuples $(p_i,e_i)$ such that $n=\prod p_i^{e_i}$.
    """
    import sage.rings.arith
    return sage.rings.arith.factor(n, algorithm=algorithm)

def crt_basis(X, xgcd=None):
    """
    Compute and return a Chinese Remainder Theorem basis for the list
    X of coprime integers.

    INPUT:
        X -- a list of Integers that are coprime in pairs
    OUTPUT:
        E -- a list of Integers such that E[i] = 1 (mod X[i])
             and E[i] = 0 (mod X[j]) for all j!=i.

    The E[i] have the property that if A is a list of objects, e.g.,
    integers, vectors, matrices, etc., where A[i] is moduli X[i], then
    a CRT lift of A is simply
                       sum E[i] * A[i].

    ALGORITHM:
    To compute E[i], compute integers s and t such that

                s * X[i] + t * (prod over i!=j of X[j]) = 1.   (*)

    Then E[i] = t * (prod over i!=j of X[j]).  Notice that equation
    (*) implies that E[i] is congruent to 1 modulo X[i] and to 0
    modulo the other X[j] for j!=i.

    COMPLEXITY: We compute len(X) extended GCD's.

    EXAMPLES:
        sage: X = [11,20,31,51]
        sage: E = crt_basis([11,20,31,51])
        sage: E[0]%X[0]; E[1]%X[0]; E[2]%X[0]; E[3]%X[0]
        1
        0
        0
        0
        sage: E[0]%X[1]; E[1]%X[1]; E[2]%X[1]; E[3]%X[1]
        0
        1
        0
        0
        sage: E[0]%X[2]; E[1]%X[2]; E[2]%X[2]; E[3]%X[2]
        0
        0
        1
        0
        sage: E[0]%X[3]; E[1]%X[3]; E[2]%X[3]; E[3]%X[3]
        0
        0
        0
        1
    """
    if not isinstance(X, list):
        raise TypeError, "X must be a list"
    if len(X) == 0:
        return []

    P = misc.mul(X)

    Y = []
    # 2. Compute extended GCD's
    ONE=X[0].parent()(1)
    for i in range(len(X)):
        p = X[i]
        prod = P//p
        g,s,t = p.xgcd(prod)
        if g != ONE:
            raise ArithmeticError, "The elements of the list X must be coprime in pairs."
        Y.append(t*prod)
    return Y



def iterator():
    yield ZZ(0)
    n = ZZ(1)
    while True:
        yield n
        yield -n
        n += 1

