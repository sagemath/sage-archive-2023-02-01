r"""
Ring $\Z$ of Integers

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
    <type 'integer.Integer'>
    sage: a + b
    6912
    sage: Z('94803849083985934859834583945394')
    94803849083985934859834583945394
"""

#*****************************************************************************
#
#   SAGE: System for Algebra and Geometry Experimentation
#
#       Copyright (C) 2005 William Stein <wstein@ucsd.edu>
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

import random

import sage.misc.misc as misc
import principal_ideal_domain
import sage.rings.integer
import sage.rings.infinity
import sage.rings.rational
import sage.rings.integer_mod
import sage.rings.rational_field
import sage.structure.factorization as factorization
import sage.libs.pari.all
import ring

_obj = None
class _uniq_int(object):
    def __new__(cls):
        global _obj
        if _obj is None:
            _obj = object.__new__(cls)
        return _obj

class IntegerRing(principal_ideal_domain.PrincipalIdealDomain, _uniq_int):
    r"""
    The ring of integers.

    In order to introduce the ring $\Z$ of integers, we illustrate
    creation, calling a few functions, and working with its
    elements.

        sage: Z = IntegerRing(); Z
        Integer Ring
        sage: Z.characteristic()
        0
        sage: Z.is_field()
        False

    We next illustrate basic arithmetic in $\Z$:

        sage: a = Z(1234); b = Z(5678); print a, b
        1234 5678
        sage: type(a)
        <type 'integer.Integer'>
        sage: a + b
        6912
        sage: b + a
        6912
        sage: a * b
        7006652
        sage: b * a
        7006652
        sage: a - b
        -4444
        sage: b - a
        4444

    When we divide to integers using /, the result is
    automatically coerced to the field of rational numbers, even
    if the result is an integer.

        sage: a / b
        617/2839
        sage: type(a/b)
        <type 'rational.Rational'>
        sage: a/a
        1
        sage: type(a/a)
        <type 'rational.Rational'>

    For floor division, instead using the // operator:
        sage: a // b
        0
        sage: type(a//b)
        <type 'integer.Integer'>

    Next we illustrate arithmetic with automatic coercion.
    The types that coerce are: str, int, long, Integer.
        sage: a + 17
        1251
        sage: a * 374
        461516
        sage: 374 * a
        461516
        sage: a/19
        1234/19
        sage: 0 + Z(-64)
        -64

    Integers can be coerced:
        sage: a = Z(-64)
        sage: int(a)
        -64
    """

    def __repr__(self):
        return "Integer Ring"

    def _latex_(self):
        return "\\mbox{\\bf{}Z}"

    def __call__(self, x):
        if isinstance(x, sage.rings.integer_mod.IntegerMod):
            x = x.lift()
        return sage.rings.integer.Integer(x)

    def __iter__(self):
        """
        Iterate over all integers.
           0 1 -1 2 -2 3 -3 ...
        """
        yield self(0)
        n = self(1)
        while True:
            yield n
            yield -n
            n += 1

    def _coerce_(self, x):
        """
        EXAMPLES:
            sage: k = GF(7)
            sage: k._coerce_(2/3)
            Traceback (most recent call last):
            ...
            TypeError: Unable to coerce 2/3 into Finite Field of size 7
            sage: k._coerce_(5)   # works since there's a natural hom ZZ --> GF(7).
            5
            sage: ZZ._coerce_(GF(7)(2))
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion of 2 to an integer
        """

        if isinstance(x, sage.rings.integer.Integer):
            return x
        elif isinstance(x, (int, long)):
            return self(x)
        elif isinstance(x, sage.libs.pari.all.pari_gen) and x.type() == 't_INT':
            return self(x)
        raise TypeError, 'no canonical coercion of %s to an integer'%x

    def _is_valid_homomorphism_(self, codomain, im_gens):
        try:
            return im_gens[0] == codomain._coerce_(self.gen(0))
        except TypeError:
            return False

    def is_atomic_repr(self):
        """
        Return True, since elements of the integers do not have
        to be printed with paranethesis around them, when they
        are coefficients, e.g., in a polynomial.
        """
        return True

    def is_field(self):
        """
        Return False.
        """
        return False

    def is_finite(self):
        return False

    def fraction_field(self):
        return sage.rings.rational_field.Q

    def gens(self):
        return (self(1), )

    def gen(self, n=0):
        if n == 0:
            return self(1)
        else:
            raise IndexError, "n must be 0"

    def ngens(self):
        return 1

    def characteristic(self):
        """
        Return 0 as a Python int.
        """
        return 0

    def krull_dimension(self):
        """
        Return the Krull dimension of the integers, which is 1.
        """
        return 1

    def random(self, bound=5):
        """
        Return a random integer between -bound and bound,
        including both endpoints.
        """
        return self(random.randint(-bound,bound))

    def order(self):
        return sage.rings.infinity.Infinity()

    def __cmp__(self, other):
        if isinstance(other, IntegerRing):
            return 0
        if ring.is_Ring(other):
            if other.characteristic() == 0:
                return -1
            else:
                return 1
        return -1

    def zeta(self, n=2):
        if n == 1:
            return sage.rings.integer.Integer(1)
        elif n == 2:
            return sage.rings.integer.Integer(-1)
        else:
            raise ValueError, "No %sth root of unity in integer ring"%n


    #################################
    ## Coercions to interfaces
    #################################
    def _gap_init_(self):
        """
        EXAMPLES:
            sage: gap(ZZ)
            Integers
        """
        return 'Integers'

    def _magma_init_(self):
        """
        EXAMPLES:
            sage: magma(ZZ)
            Integer Ring
        """
        return 'IntegerRing()'


ZZ = IntegerRing()
Z = ZZ

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


