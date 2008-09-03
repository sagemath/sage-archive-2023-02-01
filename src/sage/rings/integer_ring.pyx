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
include "../ext/interrupt.pxi"  # ctrl-c interrupt block support
include "../ext/random.pxi"

include "../ext/python_int.pxi"
include "../ext/python_list.pxi"

import sage.rings.infinity
import sage.rings.rational
import sage.rings.rational_field
import sage.rings.ideal
import sage.structure.factorization as factorization
import sage.libs.pari.all
import sage.rings.ideal
from sage.structure.parent_gens import ParentWithGens
from sage.structure.parent cimport Parent

from sage.structure.sequence import Sequence

cimport integer
cimport rational

import ring

cdef int number_of_integer_rings = 0

def is_IntegerRing(x):
    """
    Internal funtion: returns true iff x is the ring ZZ of integers

    EXAMPLES:
    sage: from sage.rings.integer_ring import  is_IntegerRing
    sage: is_IntegerRing(ZZ)
    True
    sage: is_IntegerRing(QQ)
    False
    sage: is_IntegerRing(parent(3))
    True
    sage: is_IntegerRing(parent(1/3))
    False
    """
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

    We can create integers from several types of objects.
        sage: ZZ(17/1)
        17
        sage: ZZ(Mod(19,23))
        19
        sage: ZZ(2 + 3*5 + O(5^3))
        17
    """

    def __init__(self):
        ParentWithGens.__init__(self, self, ('x',), normalize=False)
        self._populate_coercion_lists_(element_constructor=integer.Integer,
                                       init_no_parent=True,
                                       convert_method_name='_integer_')

    def __cinit__(self):
        # This is here because very old pickled integers don't have unique parents.
        global number_of_integer_rings
        if type(self) is IntegerRing_class:
            if number_of_integer_rings > 0:
                self._populate_coercion_lists_(element_constructor=integer.Integer,
                                               init_no_parent=True,
                                               convert_method_name='_integer_')
            number_of_integer_rings += 1

    def __reduce__(self):
        """
        TESTS:
            sage: loads(dumps(ZZ)) is ZZ
            True
        """
        return IntegerRing, ()

    def __hash__(self):
        return 554590422

    def __richcmp__(left, right, int op):
        return (<Parent>left)._richcmp_helper(right, op)

    def _cmp_(left, right):
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

    def __getitem__(self, x):
        """
        Return the ring ZZ[...] obtained by adjoing to the integers
        a list x of several elements.

        EXAMPLES:
            sage: ZZ[ sqrt(2), sqrt(3) ]
            Relative Order in Number Field in sqrt2 with defining polynomial x^2 - 2 over its base field
            sage: ZZ[x]
            Univariate Polynomial Ring in x over Integer Ring
            sage: ZZ['x,y']
            Multivariate Polynomial Ring in x, y over Integer Ring
            sage: R = ZZ[ sqrt(5) + 1]; R
            Order in Number Field in a with defining polynomial x^2 - 2*x - 4
            sage: R.is_maximal()
            False
            sage: R = ZZ[ (1+sqrt(5))/2 ]; R
            Order in Number Field in a with defining polynomial x^2 - x - 1
            sage: R.is_maximal()
            True

        """
        if x in self:
            return self
        if isinstance(x, str):
            return PrincipalIdealDomain.__getitem__(self, x)
        from sage.calculus.all import is_SymbolicVariable

        if is_SymbolicVariable(x):
            return PrincipalIdealDomain.__getitem__(self, repr(x))

        from sage.rings.number_field.all import is_NumberFieldElement

        if is_NumberFieldElement(x):
            K, from_K = x.parent().subfield(x)
            return K.order(K.gen())

        return PrincipalIdealDomain.__getitem__(self, x)

    def range(self, start, end=None, step=None):
        """
        Optimized range function for SAGE integer.

        AUTHOR:
          -- Robert Bradshaw (2007-09-20)

        EXAMPLES:
            sage: ZZ.range(10)
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
            sage: ZZ.range(-5,5)
            [-5, -4, -3, -2, -1, 0, 1, 2, 3, 4]
            sage: ZZ.range(0,50,5)
            [0, 5, 10, 15, 20, 25, 30, 35, 40, 45]
            sage: ZZ.range(0,50,-5)
            []
            sage: ZZ.range(50,0,-5)
            [50, 45, 40, 35, 30, 25, 20, 15, 10, 5]
            sage: ZZ.range(50,0,5)
            []
            sage: ZZ.range(50,-1,-5)
            [50, 45, 40, 35, 30, 25, 20, 15, 10, 5, 0]

        It uses different code if the step doesn't fit in a long:
            sage: ZZ.range(0,2^83,2^80)
            [0, 1208925819614629174706176, 2417851639229258349412352, 3626777458843887524118528, 4835703278458516698824704, 6044629098073145873530880, 7253554917687775048237056, 8462480737302404222943232]
        """
        if end is None:
            end = start
            start = PY_NEW(Integer) # 0
        if step is None:
            step = 1
        if not PyInt_CheckExact(step):
            if not PY_TYPE_CHECK(step, integer.Integer):
                step = integer.Integer(step)
            if mpz_fits_slong_p((<Integer>step).value):
                step = int(step)
        if not PY_TYPE_CHECK(start, integer.Integer):
            start = integer.Integer(start)
        if not PY_TYPE_CHECK(end, integer.Integer):
            start = integer.Integer(end)
        cdef integer.Integer a = <Integer>start
        cdef integer.Integer b = <Integer>end

        cdef int step_sign
        cdef long istep
        cdef integer.Integer zstep, last

        L = []
        if PyInt_CheckExact(step):
            istep = PyInt_AS_LONG(step)
            step_sign = istep
        else:
            zstep = <Integer>step
            step_sign = mpz_sgn(zstep.value)

        _sig_on
        while mpz_cmp(a.value, b.value)*step_sign < 0:
            last = a
            a = PY_NEW(Integer)
            if PyInt_CheckExact(step): # count on branch prediction...
                if istep > 0:
                    mpz_add_ui(a.value, last.value, istep)
                else:
                    mpz_sub_ui(a.value, last.value, -istep)
            else:
                mpz_add(a.value, last.value, zstep.value)
            PyList_Append(L, last)
        _sig_off
        return L

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

    cdef Integer _coerce_ZZ(self, ZZ_c *z):
        cdef integer.Integer i
        i = PY_NEW(integer.Integer)
        _sig_on
        ZZ_to_mpz(&i.value, z)
        _sig_off
        return i

    cpdef _coerce_map_from_(self, S):
        """
        x canonically coerces to the integers ZZ over only if x is an
        int, long or already an element of ZZ.

        EXAMPLES:
            sage: ZZ.coerce(int(5))
            5
            sage: ZZ.coerce(GF(7)(2))
            Traceback (most recent call last):
            ...
            TypeError: no cannonical coercion from Finite Field of size 7 to Integer Ring

        The rational number 3/1 = 3 does not canonically coerce into
        the integers, since there is no canonical coercion map from
        the full field of rational numbers to the integers.

            sage: a = 3/1; parent(a)
            Rational Field
            sage: ZZ(a)
            3
            sage: ZZ.coerce(a)
            Traceback (most recent call last):
            ...
            TypeError: no cannonical coercion from Rational Field to Integer Ring


        TESTS:
            sage: f = ZZ.coerce_map_from(int); f
            Native morphism:
              From: Set of Python objects of type 'int'
              To:   Integer Ring
            sage: f(4r)
            4
            sage: f(-7r)
            -7

        Note that the input MUST be an int.
            sage: f(10000000000000000000000r) # random
            5
        """
        if S is int:
            return sage.rings.integer.int_to_Z()
        elif S is long:
            return sage.rings.integer.long_to_Z()
        else:
            None


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

    def random_element(self, x=None, y=None, distribution=None):
        r"""
        Return a random integer.

            ZZ.random_element() -- return an integer using the default
              distribution described below
            ZZ.random_element(n) -- return an integer uniformly
              distributed between 0 and n-1, inclusive.
            ZZ.random_element(min, max) -- return an integer uniformly
              destributed between min and max-1, inclusive.

        The default distribution for ZZ.random_element() is based on
        $X = \mbox{trunc}(4/(5R))$, where $R$ is a random variable
        uniformly distributed between -1 and 1.  This gives
        $\mbox{Pr}(X = 0) = 1/5$, and $\mbox{Pr}(X = n) = 2/(5|n|(|n|+1))$ for $n \neq 0$.
        Most of the samples will be small; -1, 0, and 1 occur with
        probability 1/5 each.  But we also have a small but
        non-negligible proportion of ``outliers''; $\mbox{Pr}(|X| \geq n) = 4/(5n)$,
        so for instance, we expect that $|X| \geq 1000$ on one in
        1250 samples.

        We actually use an easy-to-compute truncation of the above
        distribution; the probabilities given above hold fairly well
        up to about $|n| = 10000$, but around $|n| = 30000$ some
        values will never be returned at all, and we will never return
        anything greater than $2^{30}$.

        EXAMPLES:
        The default distribution is on average 50% $\pm 1$:
            sage: [ZZ.random_element() for _ in range(10)]
            [-8, 2, 0, 0, 1, -1, 2, 1, -95, -1]

        The default uniform distribution is integers between -2 and 2 inclusive:
            sage: [ZZ.random_element(distribution="uniform") \
                    for _ in range(10)]
            [2, -2, 2, -2, -1, 1, -1, 2, 1, 0]


        If a range is given, the distribution is uniform in that range:
            sage: ZZ.random_element(-10,10)
            -5
            sage: ZZ.random_element(10)
            7
            sage: ZZ.random_element(10^50)
            62498971546782665598023036522931234266801185891699
            sage: [ZZ.random_element(5) for _ in range(10)]
            [1, 3, 4, 0, 3, 4, 0, 3, 0, 1]

        Notice that the right endpoint is not included:
            sage: [ZZ.random_element(-2,2) for _ in range(10)]
            [-1, -2, 0, -2, 1, -1, -1, -2, -2, 1]

        We compute a histogram over 1000 samples of the default distribution:
            sage: from collections import defaultdict
            sage: d = defaultdict(lambda: 0)
            sage: for _ in range(1000):
            ...       samp = ZZ.random_element()
            ...       d[samp] = d[samp] + 1
            sage: sorted(d.items())
            [(-1026, 1), (-248, 1), (-145, 1), (-81, 1), (-80, 1), (-79, 1), (-75, 1), (-69, 1), (-68, 1), (-63, 2), (-61, 1), (-57, 1), (-50, 1), (-37, 1), (-35, 1), (-33, 1), (-29, 2), (-27, 2), (-25, 1), (-23, 2), (-22, 2), (-20, 1), (-19, 1), (-18, 1), (-16, 4), (-15, 3), (-14, 1), (-13, 2), (-12, 2), (-11, 2), (-10, 7), (-9, 3), (-8, 3), (-7, 7), (-6, 8), (-5, 13), (-4, 24), (-3, 34), (-2, 75), (-1, 207), (0, 209), (1, 189), (2, 64), (3, 35), (4, 13), (5, 11), (6, 10), (7, 4), (8, 4), (10, 1), (11, 1), (12, 1), (13, 1), (14, 1), (16, 3), (18, 1), (19, 1), (26, 2), (27, 1), (28, 1), (29, 1), (30, 1), (32, 1), (33, 2), (35, 1), (37, 1), (39, 1), (41, 1), (42, 1), (52, 1), (91, 1), (94, 1), (106, 1), (111, 1), (113, 2), (132, 1), (134, 1), (232, 1), (240, 1), (2133, 1), (3636, 1)]
        """
        cdef integer.Integer z
        z = <integer.Integer>PY_NEW(integer.Integer)
        if x is not None and y is None and x <= 0:
            raise TypeError, "x must be > 0"
        if x is not None and y is not None and x >= y:
            raise TypeError, "x must be < y"
        self._randomize_mpz(z.value, x, y, distribution)
        return z

    cdef int _randomize_mpz(self, mpz_t value, x, y, distribution) except -1:
        cdef integer.Integer n_max, n_min, n_width
        cdef randstate rstate = current_randstate()
        cdef int den = rstate.c_random()-RAND_MAX/2
        if den == 0: den = 1
        if (distribution is None and x is None) or distribution == "1/n":
            mpz_set_si(value, (RAND_MAX/5*2) / den)
        elif distribution is None or distribution == "uniform":
            if y is None:
                if x is None:
                    mpz_set_si(value, rstate.c_random()%5 - 2)
                else:
                    n_max = x if PY_TYPE_CHECK(x, integer.Integer) else self(x)
                    mpz_urandomm(value, rstate.gmp_state, n_max.value)
            else:
                n_min = x if PY_TYPE_CHECK(x, integer.Integer) else self(x)
                n_max = y if PY_TYPE_CHECK(y, integer.Integer) else self(y)
                n_width = n_max - n_min
                if mpz_sgn(n_width.value) <= 0:
                    n_min = self(-2)
                    n_width = self(5)
                mpz_urandomm(value, rstate.gmp_state, n_width.value)
                mpz_add(value, value, n_min.value)
        elif distribution == "mpz_rrandomb":
            mpz_rrandomb(value, rstate.gmp_state, int(x))
        else:
            raise ValueError, "Unknown distribution for the integers: %s"%distribution

    def _is_valid_homomorphism_(self, codomain, im_gens):
        try:
            return im_gens[0] == codomain.coerce(self.gen(0))
        except TypeError:
            return False

    def is_noetherian(self):
        """
        Return True - the integers are a Noetherian ring.

        EXAMPLES:
            sage: ZZ.is_noetherian()
            True
        """
        return True

    def is_atomic_repr(self):
        """
        Return True, since elements of the integers do not have
        to be printed with parentheses around them, when they
        are coefficients, e.g., in a polynomial.

        EXAMPLE:
            sage: ZZ.is_atomic_repr()
            True
        """
        return True

    def is_field(self):
        """
        Return False - the integers are not a field.

        EXAMPLES:
            sage: ZZ.is_field()
            False
        """
        return False

    def is_finite(self):
        """
        Return False - the integers are an infinite ring.

        EXAMPLES:
            sage: ZZ.is_finite()
            False
        """
        return False

    def fraction_field(self):
        """
        Returns the field of rational numbers - the fraction field of the
        integers.

        EXAMPLES:
            sage: ZZ.fraction_field()
            Rational Field
            sage: ZZ.fraction_field() == QQ
            True
        """
        return sage.rings.rational_field.Q

    def extension(self, poly, names=None):
        """
        Returns the order in the number field defined by poly generated
        (as a ring) by a root of poly.

        EXAMPLES:
            sage: ZZ.extension(x^2-5, 'a')
            Order in Number Field in a with defining polynomial x^2 - 5
            sage: ZZ.extension([x^2 + 1, x^2 + 2], 'a,b')
            Relative Order in Number Field in a with defining polynomial x^2 + 1 over its base field
        """
        from sage.rings.number_field.order import EquationOrder
        return EquationOrder(poly, names)

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

    def residue_field(self, prime, check = True):
        """
        Return the residue field of the integers modulo the given prime, ie $\Z/p\Z$.

        INPUT:
            prime -- a prime number
            check -- (boolean, default True) whether or not to check the primality of prime.
        OUTPUT:
            The residue field at this prime.

        EXAMPLES:
            sage: F = ZZ.residue_field(61); F
            Residue field of Integers modulo 61
            sage: pi = F.reduction_map(); pi
            Partially defined reduction map from Rational Field to Residue field of Integers modulo 61
            sage: pi(123/234)
            6
            sage: pi(1/61)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Cannot reduce rational 1/61 modulo 61: it has negative valuation
            sage: lift = F.lift_map(); lift
            Lifting map from Residue field of Integers modulo 61 to Rational Field
            sage: lift(F(12345/67890))
            33
            sage: (12345/67890) % 61
            33

        Construction can be from a prime ideal instead of a prime:
            sage: ZZ.residue_field(ZZ.ideal(97))
            Residue field of Integers modulo 97

        TESTS:
            sage: ZZ.residue_field(ZZ.ideal(96))
            Traceback (most recent call last):
            ...
            TypeError: Principal ideal (96) of Integer Ring is not prime
            sage: ZZ.residue_field(96)
            Traceback (most recent call last):
            ...
            TypeError: 96 is not prime
        """
        if isinstance(prime, sage.rings.integer.Integer):
            p = self.ideal(prime)
        elif sage.rings.ideal.is_Ideal(prime):
            if not (prime.ring() is self):
                raise TypeError, "%s is not an ideal of ZZ"%prime
            p = prime
        else:
            raise TypeError, "%s is neither an ideal of ZZ nor an integer"%prime
        if check and not p.is_prime():
            raise TypeError, "%s is not prime"%prime
        from sage.rings.residue_field import ResidueField
        return ResidueField(p, names = None, check = check)

    def gens(self):
        """

        Returns the tuple (1,) containing a single element, the
        additive generator of the integers, which is 1.

        EXAMPLES:
            sage: ZZ.gens(); ZZ.gens()[0]
            (1,)
            1
            sage: type(ZZ.gens()[0])
            <type 'sage.rings.integer.Integer'>
        """
        return (self(1), )

    def gen(self, n=0):
        """
        Returns the additive generator of the integers, which is 1.

        EXAMPLES:
            sage: ZZ.gen()
            1
            sage: type(ZZ.gen())
            <type 'sage.rings.integer.Integer'>
        """
        if n == 0:
            return self(1)
        else:
            raise IndexError, "n must be 0"

    def ngens(self):
        """
        Returns the number of additive generators of the ring, which is 1.

        EXAMPLES:
            sage: ZZ.ngens()
            1
            sage: len(ZZ.gens())
            1
        """
        return 1

    def degree(self):
        """
        Return the degree of the integers, which is 1

        EXAMPLE:
            sage: ZZ.degree()
            1
        """
        return 1

    def absolute_degree(self):
        """
        Return the absolute degree of the integers, which is 1

        EXAMPLE:
            sage: ZZ.absolute_degree()
            1
        """
        return 1

    def characteristic(self):
        """
        Return the characteristic of the integers, which is 0

        EXAMPLE:
            sage: ZZ.characteristic()
            0
        """
        return 0

    def krull_dimension(self):
        """
        Return the Krull dimension of the integers, which is 1.

        EXAMPLE:
            sage: ZZ.krull_dimension()
            1
        """
        return 1

    def order(self):
        """
        Return the order (cardinality) of the integers, which is +Infinity.

        EXAMPLE:
            sage: ZZ.order()
            +Infinity
        """
        return sage.rings.infinity.infinity

    def zeta(self, n=2):
        """
        Return a primitive n'th root of unity in the integers, or
        raise an error if none exists

        INPUT:
            n -- a positive integer (default 2)
        OUTPUT:
            an n'th root of unity in ZZ

        EXAMPLE:
            sage: ZZ.zeta()
            -1
            sage: ZZ.zeta(1)
            1
            sage: ZZ.zeta(3)
            Traceback (most recent call last):
            ...
            ValueError: no nth root of unity in integer ring
            sage: ZZ.zeta(0)
            Traceback (most recent call last):
            ...
            ValueError: n must be positive in zeta()
        """
        if n == 1:
            return sage.rings.integer.Integer(1)
        elif n == 2:
            return sage.rings.integer.Integer(-1)
        elif n < 1:
            raise ValueError, "n must be positive in zeta()"
        else:
            raise ValueError, "no nth root of unity in integer ring"

    def parameter(self):
        """
        Returns an integer of degree 1 for the Euclidean property of
        ZZ, namely 1.

        EXAMPLES:
            sage: ZZ.parameter()
            1
        """
        return self(1)




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

    def _macaulay2_init_(self):
        """
        EXAMPLES:
            sage: macaulay2(ZZ)       #optional
            ZZ
        """
        return "ZZ"

    def _sage_input_(self, sib, coerced):
        r"""
        Produce an expression which will reproduce this value when evaluated.

        EXAMPLES:
            sage: sage_input(ZZ, verify=True)
            # Verified
            ZZ
            sage: from sage.misc.sage_input import SageInputBuilder
            sage: ZZ._sage_input_(SageInputBuilder(), False)
            {atomic:ZZ}
        """
        return sib.name('ZZ')

ZZ = IntegerRing_class()
Z = ZZ

def IntegerRing():
    """
    Return the integer ring

    EXAMPLE:
        sage: IntegerRing()
        Integer Ring
        sage: ZZ==IntegerRing()
        True
    """
    return ZZ

def factor(n, algorithm='pari', proof=True):
    """
    Return the factorization of the positive integer $n$ as a sorted
    list of tuples $(p_i,e_i)$ such that $n=\prod p_i^{e_i}$.

    For further documentation see sage.rings.arith.factor()

    EXAMPLE:
        sage: sage.rings.integer_ring.factor(420)
        2^2 * 3 * 5 * 7
    """
    import sage.rings.arith
    return sage.rings.arith.factor(n, algorithm=algorithm, proof=proof)

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
