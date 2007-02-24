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

##########################################################################################

include "../ext/cdefs.pxi"
include "../ext/gmp.pxi"
include "../ext/stdsage.pxi"

import sage.rings.infinity
import sage.rings.rational
import sage.rings.rational_field
import sage.rings.ideal
import sage.structure.factorization as factorization
import sage.libs.pari.all
import sage.rings.ideal
from sage.structure.parent_gens import ParentWithGens
from sage.structure.parent cimport Parent

cimport integer
cimport rational

import ring

###########
from random import randrange
cdef extern from "stdlib.h":
    long random()
    void srandom(unsigned int seed)
k = randrange(0,int(2)**int(32))
srandom(k)

cdef gmp_randstate_t state
gmp_randinit_mt(state)
gmp_randseed_ui(state,k)
###########

def is_IntegerRing(x):
    return PY_TYPE_CHECK(x, IntegerRing_class)

import integer_ring_python

cdef class IntegerRing_class(PrincipalIdealDomain):
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
        <type 'sage.rings.integer.Integer'>
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
        <type 'sage.rings.rational.Rational'>
        sage: a/a
        1
        sage: type(a/a)
        <type 'sage.rings.rational.Rational'>

    For floor division, instead using the // operator:
        sage: a // b
        0
        sage: type(a//b)
        <type 'sage.rings.integer.Integer'>

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

    def __init__(self):
        ParentWithGens.__init__(self, self, ('x',), normalize=False)

    def __hash__(self):
        return 554590422

    def __richcmp__(left, right, int op):
        return (<Parent>left)._richcmp(right, op)

    cdef int _cmp_c_impl(left, Parent right) except -2:
        if isinstance(right,IntegerRing_class):
            return 0
        if isinstance(right, sage.rings.rational_field.RationalField):
            return -1
        return cmp(type(left), type(right))

    def _repr_(self):
        return "Integer Ring"

    def _latex_(self):
        return "\\mathbf{Z}"

    def __len__(self):
        raise TypeError, 'len() of unsized object'

    def _div(self, integer.Integer left, integer.Integer right):
        cdef rational.Rational x
        x = PY_NEW(rational.Rational)
        if mpz_cmp_si(right.value, 0) == 0:
            raise ZeroDivisionError, 'Rational division by zero'
        mpq_set_num(x.value, left.value)
        mpq_set_den(x.value, right.value)
        mpq_canonicalize(x.value)
        return x

    def __call__(self, x, base=0):
        """
        EXAMPLES:
            sage: ZZ(17/1)
            17
            sage: ZZ(Mod(19,23))
            19
            sage: ZZ(2 + 3*5 + O(5^3))
            17
        """
        if PY_TYPE_CHECK(x, integer.Integer):
            return x
        elif PY_TYPE_CHECK(x, rational.Rational):
            return (<rational.Rational>x)._integer_c()

        cdef integer.Integer y
        y = PY_NEW(integer.Integer)
        try:
            y.__init__(x, base)
            return y
        except TypeError, msg:
            try:
                x = x.lift()
                if PY_TYPE_CHECK(x, rational.Rational):
                    return (<rational.Rational>x)._integer_c()
                y.__init__(x, base)
            except AttributeError:
                pass
        raise TypeError, msg

    def __iter__(self):
        """
        Iterate over all integers.
           0 1 -1 2 -2 3 -3 ...

        EXAMPLES:
            sage: for n in ZZ:
            ...    if n < 3: print n
            ...    else: break
            0
            1
            -1
            2
            -2
        """
        return integer_ring_python.iterator(self)

    cdef _coerce_c_impl(self, x):
        """
        Return canonical coercion of x into the integers ZZ.

        x canonically coerces to the integers ZZ over only if x is an
        int, long or already an element of ZZ.

        EXAMPLES:
            sage: k = GF(7)
            sage: k._coerce_(2/3)
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion of x
            sage: k._coerce_(5)   # works since there's a natural hom ZZ --> GF(7).
            5
            sage: ZZ._coerce_(GF(7)(2))
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion to an integer

        The rational number 3/1 = 3 does not canonically coerce into
        the integers, since there is no canonical coercion map from
        the full field of rational numbers to the integers.

            sage: a = 3/1; parent(a)
            Rational Field
            sage: ZZ(a)
            3
            sage: ZZ._coerce_(a)
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion to an integer
        """
        if isinstance(x, (int, long)):
            return self(x)
        raise TypeError, "no canonical coercion to an integer"


    def is_subring(self, other):
        """
        Return True if ZZ is a subring of other in a natural way.

        Every ring of characteristic 0 contains ZZ as a subring.

        EXAMPLES:
            sage: ZZ.is_subring(QQ)
            True
        """
        if not ring.is_Ring(other):
            raise TypeError, "other must be a ring"
        if other.characteristic() == 0:
            return True
        else:
            return False

    def random_element(self, x=None, y=None):
        """
        Return a random integer.

            ZZ.random_element() -- random integer [-2,-1,0,1,2].
            ZZ.random_element(n) -- return an integer between 0 and n-1, inclusive.
            ZZ.random_element(min, max) -- return an integer between min and max-1, inclusive.

        EXAMPLES:
        The default is integers between -2 and 2 inclusive:
            sage: [ZZ.random_element() for _ in range(10)]
            [-2, -2, 1, 1, 0, 1, 2, -2, 1, -2]

            sage: ZZ.random_element(-10,10)
            -6
            sage: ZZ.random_element(10)
            6
            sage: ZZ.random_element(10^50)
            46451269108731711203254579547654565878787536081836
            sage: [ZZ.random_element(5) for _ in range(10)]
            [3, 3, 2, 1, 0, 4, 2, 1, 1, 0]

        Notice that the right endpoint is not included:
            sage: [ZZ.random_element(-2,2) for _ in range(10)]
            [1, 0, -1, 1, 0, -2, 0, -1, 1, 0]
        """
        cdef integer.Integer z, n_max, n_min, n_width
        z = integer.Integer()
        if y is None:
            if x is None:
                mpz_set_si(z.value, random()%5 - 2)
            else:
                n_max = self(x)
                mpz_urandomm(z.value, state, n_max.value)
        else:
            n_min = self(x)
            n_width = self(y) - n_min
            mpz_urandomm(z.value, state, n_width.value)
            mpz_add(z.value, z.value, n_min.value)
        #end if
        return z

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

    def quotient(self, I, names=None):
        r"""
        Return the quotient of $\Z$ by the ideal $I$ or integer $I$.

        EXAMPLES:
            sage: ZZ.quo(6*ZZ)
            Ring of integers modulo 6
            sage: ZZ.quo(0*ZZ)
            Integer Ring
            sage: ZZ.quo(3)
            Ring of integers modulo 3
            sage: ZZ.quo(3*QQ)
            Traceback (most recent call last):
            ...
            TypeError: I must be an ideal of ZZ
        """
        if isinstance(I, sage.rings.integer.Integer):
            n = I
        elif sage.rings.ideal.is_Ideal(I):
            if not (I.ring() is self):
                raise TypeError, "I must be an ideal of ZZ"
            n = I.gens()[0]
        else:
            raise TypeError, "I must be an ideal of ZZ or an integer"
        if n == 0:
            return self
        return sage.rings.integer_mod_ring.IntegerModRing(n)

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

    def order(self):
        return sage.rings.infinity.infinity

    def zeta(self, n=2):
        if n == 1:
            return sage.rings.integer.Integer(1)
        elif n == 2:
            return sage.rings.integer.Integer(-1)
        else:
            raise ValueError, "no nth root of unity in integer ring"


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
            sage: magma(ZZ)           # optional
            Integer Ring
        """
        return 'IntegerRing()'

ZZ = IntegerRing_class()
Z = ZZ

def IntegerRing():
    return ZZ

def factor(n, algorithm='pari'):
    """
    Return the factorization of the positive integer $n$ as a list of
    tuples $(p_i,e_i)$ such that $n=\prod p_i^{e_i}$.
    """
    import sage.rings.arith
    return sage.rings.arith.factor(n, algorithm=algorithm)

import sage.misc.misc
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

    P = sage.misc.misc.prod(X)

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
