"""
NTL wrapper

AUTHORS:
   - William Stein
   - Martin Albrecht <malb@informatik.uni-bremen.de>
   - David Harvey (2007-02): speed up getting data in/out of NTL
   - Joel B. Mohler:  fast conversions to and from sage Integer type
"""

#*****************************************************************************
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

include "../../ext/interrupt.pxi"
include "../../ext/stdsage.pxi"
include 'misc.pxi'
include 'decl.pxi'

from sage.rings.integer import Integer
from sage.rings.integer_ring import IntegerRing
from sage.rings.integer cimport Integer
from sage.rings.integer_ring cimport IntegerRing_class

ZZ_sage = IntegerRing()

##############################################################################
# ZZ: Arbitrary precision integers
##############################################################################

cdef class ntl_ZZ:
    r"""
    The \class{ZZ} class is used to represent signed, arbitrary length integers.

    Routines are provided for all of the basic arithmetic operations, as
    well as for some more advanced operations such as primality testing.
    Space is automatically managed by the constructors and destructors.

    This module also provides routines for generating small primes, and
    fast routines for performing modular arithmetic on single-precision
    numbers.
    """
    # See ntl.pxd for definition of data members

    def __dealloc__(self):
        del_ZZ(self.x)

    def __repr__(self):
        _sig_on
        return string(ZZ_to_str(self.x))

    def __reduce__(self):
        raise NotImplementedError

    def __mul__(ntl_ZZ self, other):
        cdef ntl_ZZ y
        if not isinstance(other, ntl_ZZ):
            other = ntl_ZZ(other)
        y = other
        _sig_on
        return make_ZZ(ZZ_mul(self.x, y.x))

    def __sub__(ntl_ZZ self, other):
        cdef ntl_ZZ y
        if not isinstance(other, ntl_ZZ):
            other = ntl_ZZ(other)
        y = other
        _sig_on
        return make_ZZ(ZZ_sub(self.x, y.x))

    def __add__(ntl_ZZ self, other):
        cdef ntl_ZZ y
        if not isinstance(other, ntl_ZZ):
            other = ntl_ZZ(other)
        y = other
        _sig_on
        return make_ZZ(ZZ_add(self.x, y.x))

    def __neg__(ntl_ZZ self):
        _sig_on
        return make_ZZ(ZZ_neg(self.x))

    def __pow__(ntl_ZZ self, long e, ignored):
        _sig_on
        return make_ZZ(ZZ_pow(self.x, e))

    cdef set(self, void *y):  # only used internally for initialization; assumes self.x not set yet!
        self.x = <ntl_c_ZZ*> y

    cdef int get_as_int(ntl_ZZ self):
        r"""
        Returns value as C int.
        Return value is only valid if the result fits into an int.

        AUTHOR: David Harvey (2006-08-05)
        """
        return ZZ_to_int(self.x)

    def get_as_int_doctest(self):
        r"""
        This method exists solely for automated testing of get_as_int().

        sage: x = ntl.ZZ(42)
        sage: i = x.get_as_int_doctest()
        sage: print i
         42
        sage: print type(i)
         <type 'int'>
        """
        return self.get_as_int()

    def get_as_sage_int(self):
        r"""
        Gets the value as a sage int.

        sage: n=ntl.ZZ(2983)
        sage: type(n.get_as_sage_int())
        <type 'sage.rings.integer.Integer'>

        AUTHOR: Joel B. Mohler
        """
        return (<IntegerRing_class>ZZ_sage)._coerce_ZZ(self.x)

    cdef void set_from_int(ntl_ZZ self, int value):
        r"""
        Sets the value from a C int.

        AUTHOR: David Harvey (2006-08-05)
        """
        ZZ_set_from_int(self.x, value)

    def set_from_sage_int(self, Integer value):
        r"""
        Sets the value from a sage int.

        sage: n=ntl.ZZ(2983)
        sage: n
        2983
        sage: n.set_from_sage_int(1234)
        sage: n
        1234

        AUTHOR: Joel B. Mohler
        """
        value._to_ZZ(self.x)

    def set_from_int_doctest(self, value):
        r"""
        This method exists solely for automated testing of set_from_int().

        sage: x = ntl.ZZ()
        sage: x.set_from_int_doctest(42)
        sage: x
         42
        """
        self.set_from_int(int(value))

    # todo: add wrapper for int_to_ZZ in wrap.cc?


cdef public make_ZZ(ntl_c_ZZ* x):
    cdef ntl_ZZ y
    _sig_off
    y = ntl_ZZ()
    y.x = x
    return y

def make_new_ZZ(x='0'):
    s = str(x)
    cdef ntl_ZZ n
    n = ntl_ZZ()
    _sig_on
    n.x = str_to_ZZ(s)
    _sig_off
    return n

# Random-number generation
def ntl_setSeed(x=None):
    """
    Seed the NTL random number generator.

    EXAMPLE:
        sage: ntl.ntl_setSeed(10)
        sage: ntl.ZZ_random(1000)
        776
    """
    cdef ntl_ZZ seed
    if x is None:
        from random import randint
        seed = make_new_ZZ(str(randint(0,int(2)**64)))
    else:
        seed = make_new_ZZ(str(x))
    _sig_on
    setSeed(seed.x)
    _sig_off

ntl_setSeed()

def randomBnd(q):
    r"""
    Returns cryptographically-secure random number in the range [0,n)

    EXAMPLES:
        sage: [ntl.ZZ_random(99999) for i in range(5)]
        [53357, 19674, 69528, 87029, 28752]

    AUTHOR:
        -- Didier Deshommes <dfdeshom@gmail.com>
    """
    cdef ntl_ZZ w

    if not isinstance(q, ntl_ZZ):
        q = make_new_ZZ(str(q))
    w = q
    _sig_on
    return  make_ZZ(ZZ_randomBnd(w.x))

def randomBits(long n):
    r"""
    Return a pseudo-random number between 0 and $2^n-1$

    EXAMPLES:
        sage: [ntl.ZZ_random_bits(20) for i in range(3)]
        [1025619, 177635, 766262]

    AUTHOR:
        -- Didier Deshommes <dfdeshom@gmail.com>
    """
    _sig_on
    return make_ZZ(ZZ_randomBits(n))


##############################################################################
#
# ZZX: polynomials over the integers
#
##############################################################################


cdef class ntl_ZZX:
    r"""
    The class \class{ZZX} implements polynomials in $\Z[X]$, i.e.,
    univariate polynomials with integer coefficients.

    Polynomial multiplication is very fast, and is implemented using
    one of 4 different algorithms:
    \begin{enumerate}
    \item\hspace{1em} Classical
    \item\hspace{1em} Karatsuba
    \item\hspace{1em} Schoenhage and Strassen --- performs an FFT by working
          modulo a "Fermat number" of appropriate size...
          good for polynomials with huge coefficients
          and moderate degree
    \item\hspace{1em} CRT/FFT --- performs an FFT by working modulo several
         small primes.  This is good for polynomials with moderate
         coefficients and huge degree.
    \end{enumerate}

    The choice of algorithm is somewhat heuristic, and may not always be
    perfect.

    Many thanks to Juergen Gerhard {\tt
    <jngerhar@plato.uni-paderborn.de>} for pointing out the deficiency
    in the NTL-1.0 ZZX arithmetic, and for contributing the
    Schoenhage/Strassen code.

    Extensive use is made of modular algorithms to enhance performance
    (e.g., the GCD algorithm and many others).
    """

    # See ntl_ZZX.pxd for definition of data members
    def __init__(self):
        """
        EXAMPLES:
            sage: f = ntl.ZZX([1,2,5,-9])
            sage: f
            [1 2 5 -9]
            sage: g = ntl.ZZX([0,0,0]); g
            []
            sage: g[10]=5
            sage: g
            [0 0 0 0 0 0 0 0 0 0 5]
            sage: g[10]
            5
        """
        return

    def __reduce__(self):
        raise NotImplementedError

    def __dealloc__(self):
        if self.x:
            ZZX_dealloc(self.x)

    def __repr__(self):
        return str(ZZX_repr(self.x))

    def copy(self):
        return make_ZZX(ZZX_copy(self.x))

    def __copy__(self):
        return make_ZZX(ZZX_copy(self.x))

    def __setitem__(self, long i, a):
        if i < 0:
            raise IndexError, "index (i=%s) must be >= 0"%i
        a = str(int(a))
        ZZX_setitem(self.x, i, a)

    cdef void setitem_from_int(ntl_ZZX self, long i, int value):
        r"""
        Sets ith coefficient to value.

        AUTHOR: David Harvey (2006-08-05)
        """
        ZZX_setitem_from_int(self.x, i, value)

    def setitem_from_int_doctest(self, i, value):
        r"""
        This method exists solely for automated testing of setitem_from_int().

        sage: x = ntl.ZZX([2, 3, 4])
        sage: x.setitem_from_int_doctest(5, 42)
        sage: x
         [2 3 4 0 0 42]
        """
        self.setitem_from_int(int(i), int(value))

    def __getitem__(self, unsigned int i):
        r"""
        Retrieves coefficient #i as a SAGE Integer.

        sage: x = ntl.ZZX([129381729371289371237128318293718237, 2, -3, 0, 4])
        sage: x[0]
         129381729371289371237128318293718237
        sage: type(x[0])
         <type 'sage.rings.integer.Integer'>
        sage: x[1]
         2
        sage: x[2]
         -3
        sage: x[3]
         0
        sage: x[4]
         4
        sage: x[5]
         0
        """
        cdef Integer output
        output = PY_NEW(Integer)
        ZZX_getitem_as_mpz(&output.value, self.x, i)
        return output

    cdef int getitem_as_int(ntl_ZZX self, long i):
        r"""
        Returns ith coefficient as C int.
        Return value is only valid if the result fits into an int.

        AUTHOR: David Harvey (2006-08-05)
        """
        return ZZX_getitem_as_int(self.x, i)

    def getitem_as_int_doctest(self, i):
        r"""
        This method exists solely for automated testing of getitem_as_int().

        sage: x = ntl.ZZX([2, 3, 5, -7, 11])
        sage: i = x.getitem_as_int_doctest(3)
        sage: print i
         -7
        sage: print type(i)
         <type 'int'>
        sage: print x.getitem_as_int_doctest(15)
         0
        """
        return self.getitem_as_int(i)

    def list(self):
        r"""
        Retrieves coefficients as a list of SAGE Integers.

        EXAMPLES:
            sage: x = ntl.ZZX([129381729371289371237128318293718237, 2, -3, 0, 4])
            sage: L = x.list(); L
            [129381729371289371237128318293718237, 2, -3, 0, 4]
            sage: type(L[0])
            <type 'sage.rings.integer.Integer'>
            sage: x = ntl.ZZX()
            sage: L = x.list(); L
            []
        """
        cdef int i
        return [self[i] for i from 0 <= i <= ZZX_degree(self.x)]


    def __add__(ntl_ZZX self, ntl_ZZX other):
        """
        EXAMPLES:
            sage: ntl.ZZX(range(5)) + ntl.ZZX(range(6))
            [0 2 4 6 8 5]
        """
        return make_ZZX(ZZX_add(self.x, other.x))

    def __sub__(ntl_ZZX self, ntl_ZZX other):
        """
        EXAMPLES:
            sage: ntl.ZZX(range(5)) - ntl.ZZX(range(6))
            [0 0 0 0 0 -5]
        """
        return make_ZZX(ZZX_sub(self.x, other.x))

    def __mul__(ntl_ZZX self, ntl_ZZX other):
        """
        EXAMPLES:
            sage: ntl.ZZX(range(5)) * ntl.ZZX(range(6))
            [0 0 1 4 10 20 30 34 31 20]
        """
        _sig_on
        return make_ZZX(ZZX_mul(self.x, other.x))

    def __div__(ntl_ZZX self, ntl_ZZX other):
        """
        Compute quotient self / other, if the quotient is a polynomial.
        Otherwise an Exception is raised.

        EXAMPLES:
            sage: f = ntl.ZZX([1,2,3]) * ntl.ZZX([4,5])**2
            sage: g = ntl.ZZX([4,5])
            sage: f/g
            [4 13 22 15]
            sage: ntl.ZZX([1,2,3]) * ntl.ZZX([4,5])
            [4 13 22 15]

            sage: f = ntl.ZZX(range(10)); g = ntl.ZZX([-1,0,1])
            sage: f/g
            Traceback (most recent call last):
            ...
            ArithmeticError: self (=[0 1 2 3 4 5 6 7 8 9]) is not divisible by other (=[-1 0 1])
        """
        _sig_on
        cdef int divisible
        cdef ntl_c_ZZX* q
        q = ZZX_div(self.x, other.x, &divisible)
        if not divisible:
            raise ArithmeticError, "self (=%s) is not divisible by other (=%s)"%(self, other)
        return make_ZZX(q)

    def __mod__(ntl_ZZX self, ntl_ZZX other):
        """
        Given polynomials a, b in ZZ[X], there exist polynomials q, r
        in QQ[X] such that a = b*q + r, deg(r) < deg(b).  This
        function returns q if q lies in ZZ[X], and otherwise raises an
        Exception.

        EXAMPLES:
            sage: f = ntl.ZZX([2,4,6]); g = ntl.ZZX([2])
            sage: f % g   # 0
            []

            sage: f = ntl.ZZX(range(10)); g = ntl.ZZX([-1,0,1])
            sage: f % g
            [20 25]
        """
        _sig_on
        return make_ZZX(ZZX_mod(self.x, other.x))

    def quo_rem(self, ntl_ZZX other):
        """
        Returns the unique integral q and r such that self = q*other +
        r, if they exist.  Otherwise raises an Exception.

        EXAMPLES:
           sage: f = ntl.ZZX(range(10)); g = ntl.ZZX([-1,0,1])
           sage: q, r = f.quo_rem(g)
           sage: q, r
           ([20 24 18 21 14 16 8 9], [20 25])
           sage: q*g + r == f
           True
        """
        cdef ntl_c_ZZX *r, *q
        _sig_on
        ZZX_quo_rem(self.x, other.x, &r, &q)
        return (make_ZZX(q), make_ZZX(r))

    def square(self):
        """
        Return f*f.

        EXAMPLES:
            sage: f = ntl.ZZX([-1,0,1])
            sage: f*f
            [1 0 -2 0 1]
        """
        _sig_on
        return make_ZZX(ZZX_square(self.x))

    def __pow__(ntl_ZZX self, long n, ignored):
        """
        Return the n-th nonnegative power of self.

        EXAMPLES:
            sage: g = ntl.ZZX([-1,0,1])
            sage: g**10
            [1 0 -10 0 45 0 -120 0 210 0 -252 0 210 0 -120 0 45 0 -10 0 1]
        """
        if n < 0:
            raise NotImplementedError
        import sage.rings.arith
        return sage.rings.arith.generic_power(self, n, make_new_ZZX([1]))

    def __cmp__(ntl_ZZX self, ntl_ZZX other):
        """
        Decide whether or not self and other are equal.

        EXAMPLES:
            sage: f = ntl.ZZX([1,2,3])
            sage: g = ntl.ZZX([1,2,3,0])
            sage: f == g
            True
            sage: g = ntl.ZZX([0,1,2,3])
            sage: f == g
            False
        """
        if ZZX_equal(self.x, other.x):
            return 0
        return -1

    def is_zero(self):
        """
        Return True exactly if this polynomial is 0.

        EXAMPLES:
            sage: f = ntl.ZZX([0,0,0,0])
            sage: f.is_zero()
            True
            sage: f = ntl.ZZX([0,0,1])
            sage: f
            [0 0 1]
            sage: f.is_zero()
            False
        """
        return bool(ZZX_is_zero(self.x))

    def is_one(self):
        """
        Return True exactly if this polynomial is 1.

        EXAMPLES:
            sage: f = ntl.ZZX([1,1])
            sage: f.is_one()
            False
            sage: f = ntl.ZZX([1])
            sage: f.is_one()
            True
        """
        return bool(ZZX_is_one(self.x))

    def is_monic(self):
        """
        Return True exactly if this polynomial is monic.

        EXAMPLES:
            sage: f = ntl.ZZX([2,0,0,1])
            sage: f.is_monic()
            True
            sage: g = f.reverse()
            sage: g.is_monic()
            False
            sage: g
            [1 0 0 2]
        """
        return bool(ZZX_is_monic(self.x))

    def __neg__(self):
        """
        Return the negative of self.
        EXAMPLES:
            sage: f = ntl.ZZX([2,0,0,1])
            sage: -f
            [-2 0 0 -1]
        """
        return make_ZZX(ZZX_neg(self.x))

    def left_shift(self, long n):
        """
        Return the polynomial obtained by shifting all coefficients of
        this polynomial to the left n positions.

        EXAMPLES:
            sage: f = ntl.ZZX([2,0,0,1])
            sage: f
            [2 0 0 1]
            sage: f.left_shift(2)
            [0 0 2 0 0 1]
            sage: f.left_shift(5)
            [0 0 0 0 0 2 0 0 1]

        A negative left shift is a right shift.
            sage: f.left_shift(-2)
            [0 1]
        """
        return make_ZZX(ZZX_left_shift(self.x, n))

    def right_shift(self, long n):
        """
        Return the polynomial obtained by shifting all coefficients of
        this polynomial to the right n positions.

        EXAMPLES:
            sage: f = ntl.ZZX([2,0,0,1])
            sage: f
            [2 0 0 1]
            sage: f.right_shift(2)
            [0 1]
            sage: f.right_shift(5)
            []
            sage: f.right_shift(-2)
            [0 0 2 0 0 1]
        """
        return make_ZZX(ZZX_right_shift(self.x, n))

    def content(self):
        """
        Return the content of f, which has sign the same as the
        leading coefficient of f.  Also, our convention is that the
        content of 0 is 0.

        EXAMPLES:
            sage: f = ntl.ZZX([2,0,0,2])
            sage: f.content()
            2
            sage: f = ntl.ZZX([2,0,0,-2])
            sage: f.content()
            -2
            sage: f = ntl.ZZX([6,12,3,9])
            sage: f.content()
            3
            sage: f = ntl.ZZX([])
            sage: f.content()
            0
        """
        cdef char* t
        t = ZZX_content(self.x)
        return int(string(t))

    def primitive_part(self):
        """
        Return the primitive part of f.  Our convention is that the leading
        coefficient of the primitive part is nonnegegative, and the primitive
        part of 0 is 0.

        EXAMPLES:
            sage: f = ntl.ZZX([6,12,3,9])
            sage: f.primitive_part()
            [2 4 1 3]
            sage: f
            [6 12 3 9]
            sage: f = ntl.ZZX([6,12,3,-9])
            sage: f
            [6 12 3 -9]
            sage: f.primitive_part()
            [-2 -4 -1 3]
            sage: f = ntl.ZZX()
            sage: f.primitive_part()
            []
        """
        return make_ZZX(ZZX_primitive_part(self.x))

    def pseudo_quo_rem(self, ntl_ZZX other):
        r"""
        Performs pseudo-division: computes q and r with deg(r) <
        deg(b), and \code{LeadCoeff(b)\^(deg(a)-deg(b)+1) a = b q + r}.
        Only the classical algorithm is used.

        EXAMPLES:
            sage: f = ntl.ZZX([0,1])
            sage: g = ntl.ZZX([1,2,3])
            sage: g.pseudo_quo_rem(f)
            ([2 3], [1])
            sage: f = ntl.ZZX([1,1])
            sage: g.pseudo_quo_rem(f)
            ([-1 3], [2])
        """
        cdef ntl_c_ZZX *r, *q
        _sig_on
        ZZX_pseudo_quo_rem(self.x, other.x, &r, &q)
        return (make_ZZX(q), make_ZZX(r))

    def gcd(self, ntl_ZZX other):
        """
        Return the gcd d = gcd(a, b), where by convention the leading coefficient
        of d is >= 0.  We use a multi-modular algorithm.

        EXAMPLES:
            sage: f = ntl.ZZX([1,2,3]) * ntl.ZZX([4,5])**2
            sage: g = ntl.ZZX([1,1,1])**3 * ntl.ZZX([1,2,3])
            sage: f.gcd(g)
            [1 2 3]
            sage: g.gcd(f)
            [1 2 3]
        """
        _sig_on
        return make_ZZX(ZZX_gcd(self.x, other.x))

    def lcm(self, ntl_ZZX other):
        """
        Return the least common multiple of self and other.
        """
        g = self.gcd(other)
        return (self*other).quo_rem(g)[0]

    def xgcd(self, ntl_ZZX other, proof=True):
        """
        If self and other are coprime over the rationals, return r, s,
        t such that r = s*self + t*other.  Otherwise return 0.  This
        is \emph{not} the same as the \sage function on polynomials
        over the integers, since here the return value r is always an
        integer.

        Here r is the resultant of a and b; if r != 0, then this
        function computes s and t such that: a*s + b*t = r; otherwise
        s and t are both 0.  If proof = False (*not* the default),
        then resultant computation may use a randomized strategy that
        errs with probability no more than $2^{-80}$.

        EXAMPLES:
            sage: f = ntl.ZZX([1,2,3]) * ntl.ZZX([4,5])**2
            sage: g = ntl.ZZX([1,1,1])**3 * ntl.ZZX([1,2,3])
            sage: f.xgcd(g)   # nothing since they are not coprime
            (0, [], [])

        In this example the input quadratic polynomials have a common root modulo 13.
            sage: f = ntl.ZZX([5,0,1])
            sage: g = ntl.ZZX([18,0,1])
            sage: f.xgcd(g)
            (169, [-13], [13])
        """
        cdef ntl_c_ZZX *s, *t
        cdef ntl_c_ZZ *r
        _sig_on
        ZZX_xgcd(self.x, other.x, &r, &s, &t, proof)
        return (make_ZZ(r), make_ZZX(s), make_ZZX(t))

    def degree(self):
        """
        Return the degree of this polynomial.  The degree of the 0
        polynomial is -1.

        EXAMPLES:
            sage: f = ntl.ZZX([5,0,1])
            sage: f.degree()
            2
            sage: f = ntl.ZZX(range(100))
            sage: f.degree()
            99
            sage: f = ntl.ZZX()
            sage: f.degree()
            -1
            sage: f = ntl.ZZX([1])
            sage: f.degree()
            0
        """
        return ZZX_degree(self.x)

    def leading_coefficient(self):
        """
        Return the leading coefficient of this polynomial.

        EXAMPLES:
            sage: f = ntl.ZZX([3,6,9])
            sage: f.leading_coefficient()
            9
            sage: f = ntl.ZZX()
            sage: f.leading_coefficient()
            0
        """
        return make_ZZ(ZZX_leading_coefficient(self.x))

    def constant_term(self):
        """
        Return the constant coefficient of this polynomial.

        EXAMPLES:
            sage: f = ntl.ZZX([3,6,9])
            sage: f.constant_term()
            3
            sage: f = ntl.ZZX()
            sage: f.constant_term()
            0
        """
        cdef char* t
        t = ZZX_constant_term(self.x)
        return int(string(t))

    def set_x(self):
        """
        Set this polynomial to the monomial "x".

        EXAMPLES:
            sage: f = ntl.ZZX()
            sage: f.set_x()
            sage: f
            [0 1]
            sage: g = ntl.ZZX([0,1])
            sage: f == g
            True

        Though f and g are equal, they are not the same objects in memory:
            sage: f is g
            False

        """
        ZZX_set_x(self.x)

    def is_x(self):
        """
        True if this is the polynomial "x".

        EXAMPLES:
            sage: f = ntl.ZZX()
            sage: f.set_x()
            sage: f.is_x()
            True
            sage: f = ntl.ZZX([0,1])
            sage: f.is_x()
            True
            sage: f = ntl.ZZX([1])
            sage: f.is_x()
            False
        """
        return bool(ZZX_is_x(self.x))

    def derivative(self):
        """
        Return the derivative of this polynomial.

        EXAMPLES:
            sage: f = ntl.ZZX([1,7,0,13])
            sage: f.derivative()
            [7 0 39]
        """
        return make_ZZX(ZZX_derivative(self.x))

    def reverse(self, hi=None):
        """
        Return the polynomial obtained by reversing the coefficients
        of this polynomial.  If hi is set then this function behaves
        as if this polynomial has degree hi.

        EXAMPLES:
            sage: f = ntl.ZZX([1,2,3,4,5])
            sage: f.reverse()
            [5 4 3 2 1]
            sage: f.reverse(hi=10)
            [0 0 0 0 0 0 5 4 3 2 1]
            sage: f.reverse(hi=2)
            [3 2 1]
            sage: f.reverse(hi=-2)
            []
        """
        if not (hi is None):
            return make_ZZX(ZZX_reverse_hi(self.x, int(hi)))
        else:
            return make_ZZX(ZZX_reverse(self.x))

    def truncate(self, long m):
        """
        Return the truncation of this polynomial obtained by
        removing all terms of degree >= m.

        EXAMPLES:
            sage: f = ntl.ZZX([1,2,3,4,5])
            sage: f.truncate(3)
            [1 2 3]
            sage: f.truncate(8)
            [1 2 3 4 5]
            sage: f.truncate(1)
            [1]
            sage: f.truncate(0)
            []
            sage: f.truncate(-1)
            []
            sage: f.truncate(-5)
            []
        """
        if m <= 0:
            return make_new_ZZX()
        return make_ZZX(ZZX_truncate(self.x, m))

    def multiply_and_truncate(self, ntl_ZZX other, long m):
        """
        Return self*other but with terms of degree >= m removed.

        EXAMPLES:
            sage: f = ntl.ZZX([1,2,3,4,5])
            sage: g = ntl.ZZX([10])
            sage: f.multiply_and_truncate(g, 2)
            [10 20]
            sage: g.multiply_and_truncate(f, 2)
            [10 20]
        """
        if m <= 0:
            return make_new_ZZX()
        return make_ZZX(ZZX_multiply_and_truncate(self.x, other.x, m))

    def square_and_truncate(self, long m):
        """
        Return self*self but with terms of degree >= m removed.

        EXAMPLES:
            sage: f = ntl.ZZX([1,2,3,4,5])
            sage: f.square_and_truncate(4)
            [1 4 10 20]
            sage: (f*f).truncate(4)
            [1 4 10 20]
        """
        if m < 0:
            return make_new_ZZX()
        return make_ZZX(ZZX_square_and_truncate(self.x, m))

    def invert_and_truncate(self, long m):
        """
        Compute and return the inverse of self modulo $x^m$.
        The constant term of self must be 1 or -1.

        EXAMPLES:
            sage: f = ntl.ZZX([1,2,3,4,5,6,7])
            sage: f.invert_and_truncate(20)
            [1 -2 1 0 0 0 0 8 -23 22 -7 0 0 0 64 -240 337 -210 49]
            sage: g = f.invert_and_truncate(20)
            sage: g * f
            [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -512 1344 -1176 343]
        """
        if m < 0:
            raise ArithmeticError, "m (=%s) must be positive"%m
        n = self.constant_term()
        if n != 1 and n != -1:
            raise ArithmeticError, \
                  "The constant term of self must be 1 or -1."
        _sig_on
        return make_ZZX(ZZX_invert_and_truncate(self.x, m))

    def multiply_mod(self, ntl_ZZX other, ntl_ZZX modulus):
        """
        Return self*other % modulus.  The modulus must be monic with
        deg(modulus) > 0, and self and other must have smaller degree.

        EXAMPLES:
            sage: modulus = ntl.ZZX([1,2,0,1])    # must be monic
            sage: g = ntl.ZZX([-1,0,1])
            sage: h = ntl.ZZX([3,7,13])
            sage: h.multiply_mod(g, modulus)
            [-10 -34 -36]
        """
        _sig_on
        return make_ZZX(ZZX_multiply_mod(self.x, other.x, modulus.x))

    def trace_mod(self, ntl_ZZX modulus):
        """
        Return the trace of this polynomial modulus the modulus.
        The modulus must be monic, and of positive degree degree bigger
        than the the degree of self.

        EXAMPLES:
            sage: f = ntl.ZZX([1,2,0,3])
            sage: mod = ntl.ZZX([5,3,-1,1,1])
            sage: f.trace_mod(mod)
            -37
        """
        _sig_on
        return make_ZZ(ZZX_trace_mod(self.x, modulus.x))

    def trace_list(self):
        """
        Return the list of traces of the powers $x^i$ of the
        monomial x modulo this polynomial for i = 0, ..., deg(f)-1.
        This polynomial must be monic.

        EXAMPLES:
            sage: f = ntl.ZZX([1,2,0,3,0,1])
            sage: f.trace_list()
            [5, 0, -6, 0, 10]

        The input polynomial must be monic or a ValueError is raised:
            sage: f = ntl.ZZX([1,2,0,3,0,2])
            sage: f.trace_list()
            Traceback (most recent call last):
            ...
            ValueError: polynomial must be monic.
        """
        if not self.is_monic():
            raise ValueError, "polynomial must be monic."
        _sig_on
        cdef char* t
        t = ZZX_trace_list(self.x)
        return eval(string(t).replace(' ', ','))

    def resultant(self, ntl_ZZX other, proof=True):
        """
        Return the resultant of self and other.  If proof = False (the
        default is proof=True), then this function may use a
        randomized strategy that errors with probability no more than
        $2^{-80}$.

        EXAMPLES:
            sage: f = ntl.ZZX([17,0,1,1])
            sage: g = ntl.ZZX([34,-17,18,2])
            sage: f.resultant(g)
            1345873
            sage: f.resultant(g, proof=False)
            1345873
        """
        # NOTES: Within a factor of 2 in speed compared to MAGMA.
        _sig_on
        return make_ZZ(ZZX_resultant(self.x, other.x, proof))

    def norm_mod(self, ntl_ZZX modulus, proof=True):
        """
        Return the norm of this polynomial modulo the modulus.  The
        modulus must be monic, and of positive degree strictly greater
        than the degree of self.  If proof=False (the default is
        proof=True) then it may use a randomized strategy that errors
        with probability no more than $2^{-80}$.


        EXAMPLE:
            sage: f = ntl.ZZX([1,2,0,3])
            sage: mod = ntl.ZZX([-5,2,0,0,1])
            sage: f.norm_mod(mod)
            -8846

        The norm is the constant term of the characteristic polynomial.
            sage: f.charpoly_mod(mod)
            [-8846 -594 -60 14 1]
        """
        _sig_on
        return make_ZZ(ZZX_norm_mod(self.x, modulus.x, proof))

    def discriminant(self, proof=True):
        r"""
        Return the discriminant of self, which is by definition
        $$
                (-1)^{m(m-1)/2} {\mbox{\tt resultant}}(a, a')/lc(a),
        $$
        where m = deg(a), and lc(a) is the leading coefficient of a.
        If proof is False (the default is True), then this function
        may use a randomized strategy that errors with probability no
        more than $2^{-80}$.
        EXAMPLES:
            sage: f = ntl.ZZX([1,2,0,3])
            sage: f.discriminant()
            -339
            sage: f.discriminant(proof=False)
            -339
        """
        _sig_on
        return make_ZZ(ZZX_discriminant(self.x, proof))

    #def __call__(self, ntl_ZZ a):
    #    _sig_on
    #    return make_ZZ(ZZX_polyeval(self.x, a.x))

    def charpoly_mod(self, ntl_ZZX modulus, proof=True):
        """
        Return the characteristic polynomial of this polynomial modulo
        the modulus.  The modulus must be monic of degree bigger than
        self. If proof=False (the default is True), then this function
        may use a randomized strategy that errors with probability no
        more than $2^{-80}$.

        EXAMPLES:
            sage: f = ntl.ZZX([1,2,0,3])
            sage: mod = ntl.ZZX([-5,2,0,0,1])
            sage: f.charpoly_mod(mod)
            [-8846 -594 -60 14 1]

        """
        _sig_on
        return make_ZZX(ZZX_charpoly_mod(self.x, modulus.x, proof))

    def minpoly_mod_noproof(self, ntl_ZZX modulus):
        """
        Return the minimal polynomial of this polynomial modulo the
        modulus.  The modulus must be monic of degree bigger than
        self.  In all cases, this function may use a randomized
        strategy that errors with probability no more than $2^{-80}$.

        EXAMPLES:
            sage: f = ntl.ZZX([0,0,1])
            sage: g = f*f
            sage: f.charpoly_mod(g)
            [0 0 0 0 1]

        However, since $f^2 = 0$ moduluo $g$, its minimal polynomial
        is of degree $2$.
            sage: f.minpoly_mod_noproof(g)
            [0 0 1]
        """
        _sig_on
        return make_ZZX(ZZX_minpoly_mod(self.x, modulus.x))

    def clear(self):
        """
        Reset this polynomial to 0.  Changes this polynomial in place.

        EXAMPLES:
            sage: f = ntl.ZZX([1,2,3])
            sage: f
            [1 2 3]
            sage: f.clear()
            sage: f
            []
        """
        ZZX_clear(self.x)

    def preallocate_space(self, long n):
        """
        Pre-allocate spaces for n coefficients.  The polynomial that f
        represents is unchanged.  This is useful if you know you'll be
        setting coefficients up to n, so memory isn't re-allocated as
        the polynomial grows.  (You might save a millisecond with this
        function.)

        EXAMPLES:
            sage: f = ntl.ZZX([1,2,3])
            sage: f.preallocate_space(20)
            sage: f
            [1 2 3]
            sage: f[10]=5  # no new memory is allocated
            sage: f
            [1 2 3 0 0 0 0 0 0 0 5]
        """
        _sig_on
        ZZX_preallocate_space(self.x, n)
        _sig_off


    cdef set(self, void* x):  # only used internally for initialization; assumes self.x not set yet!
        self.x = <ntl_c_ZZX*>x

cdef make_ZZX(ntl_c_ZZX* x):
    cdef ntl_ZZX y
    _sig_off
    y = ntl_ZZX()
    y.x = x
    return y

def make_new_ZZX(v=[]):
    s = str(v).replace(',',' ')
    cdef ntl_ZZX z
    z = ntl_ZZX()
    _sig_on
    z.x = str_to_ZZX(s)
    _sig_off
    return z


##############################################################################
#
# ZZ_p: integers modulo p
#
##############################################################################
cdef class ntl_ZZ_p:
    r"""
    The \class{ZZ_p} class is used to represent integers modulo $p$.
    The modulus $p$ may be any positive integer, not necessarily prime.

    Objects of the class \class{ZZ_p} are represented as a \code{ZZ} in the
    range $0, \ldots, p-1$.

    An executing program maintains a "current modulus", which is set to p
    with ntl_ZZ_p.init(p).  The current modulus should be initialized before
    any ZZ_p objects are created.

    The modulus may be changed, and a mechanism is provided for saving and
    restoring a modulus (see classes ZZ_pBak and ZZ_pContext below).

    TODO: This documentation is wrong
    """
    # See ntl.pxd for definition of data members

    def __reduce__(self):
        raise NotImplementedError

    def __dealloc__(self):
        del_ZZ_p(self.x)

    def __repr__(self):
        _sig_on
        return string(ZZ_p_to_str(self.x))

    def __cmp__(ntl_ZZ_p self, ntl_ZZ_p other):
        cdef int t
        _sig_on
        t = ZZ_p_eq(self.x, other.x)
        _sig_off
        if t:
            return 0
        return 1

    def __invert__(ntl_ZZ_p self):
        _sig_on
        return make_ZZ_p(ZZ_p_inv(self.x))


    def __mul__(ntl_ZZ_p self, other):
        cdef ntl_ZZ_p y
        if not isinstance(other, ntl_ZZ_p):
            other = ntl_ZZ_p(other)
        y = other
        _sig_on
        return make_ZZ_p(ZZ_p_mul(self.x, y.x))

    def __sub__(ntl_ZZ_p self, other):
        cdef ntl_ZZ_p y
        if not isinstance(other, ntl_ZZ_p):
            other = ntl_ZZ_p(other)
        y = other
        _sig_on
        return make_ZZ_p(ZZ_p_sub(self.x, y.x))

    def __add__(ntl_ZZ_p self, other):
        cdef ntl_ZZ_p y
        if not isinstance(other, ntl_ZZ_p):
            other = ntl_ZZ_p(other)
        y = other
        _sig_on
        return make_ZZ_p(ZZ_p_add(self.x, y.x))

    def __neg__(ntl_ZZ_p self):
        _sig_on
        return make_ZZ_p(ZZ_p_neg(self.x))

    def __pow__(ntl_ZZ_p self, long e, ignored):
        _sig_on
        return make_ZZ_p(ZZ_p_pow(self.x, e))

    cdef set(self, void *y):  # only used internally for initialization; assumes self.x not set yet!
        self.x = <ZZ_p*> y

    cdef int get_as_int(ntl_ZZ_p self):
        r"""
        Returns value as C int.
        Return value is only valid if the result fits into an int.

        AUTHOR: David Harvey (2006-08-05)
        """
        return ZZ_p_to_int(self.x)

    def get_as_int_doctest(self):
        r"""
        This method exists solely for automated testing of get_as_int().

        sage: ntl.set_modulus(ntl.ZZ(20))
        sage: x = ntl.ZZ_p(42)
        sage: i = x.get_as_int_doctest()
        sage: print i
         2
        sage: print type(i)
         <type 'int'>
        """
        return self.get_as_int()

    cdef void set_from_int(ntl_ZZ_p self, int value):
        r"""
        Sets the value from a C int.

        AUTHOR: David Harvey (2006-08-05)
        """
        ZZ_p_set_from_int(self.x, value)

    def set_from_int_doctest(self, value):
        r"""
        This method exists solely for automated testing of set_from_int().

        sage: ntl.set_modulus(ntl.ZZ(20))
        sage: x = ntl.ZZ_p()
        sage: x.set_from_int_doctest(42)
        sage: x
         2
        """
        self.set_from_int(int(value))

    # todo: add wrapper for int_to_ZZ_p in wrap.cc?


cdef public make_ZZ_p(ZZ_p* x):
    cdef ntl_ZZ_p y
    _sig_off
    y = ntl_ZZ_p()
    y.x = x
    return y

def make_new_ZZ_p(x='0'):
    s = str(x)
    cdef ntl_ZZ_p n
    n = ntl_ZZ_p()
    _sig_on
    n.x = str_to_ZZ_p(s)
    _sig_off
    return n

def set_ZZ_p_modulus(ntl_ZZ p):
    ntl_ZZ_set_modulus(<ntl_c_ZZ*>p.x)


def ntl_ZZ_p_random():
    """
    Return a random number modulo p.
    """
    _sig_on
    return make_ZZ_p(ZZ_p_random())

##############################################################################
#
# ZZ_pX  -- polynomials over the integers modulo p
#
##############################################################################

cdef class ntl_ZZ_pX:
    r"""
    The class \class{ZZ_pX} implements polynomial arithmetic modulo $p$.

    Polynomial arithmetic is implemented using the FFT, combined with
    the Chinese Remainder Theorem.  A more detailed description of the
    techniques used here can be found in [Shoup, J. Symbolic
    Comp. 20:363-397, 1995].

    Small degree polynomials are multiplied either with classical
    or Karatsuba algorithms.
    """
    # See ntl_ZZ_pX.pxd for definition of data members
    def __init__(self):
        """
        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(20))
            sage: f = ntl.ZZ_pX([1,2,5,-9])
            sage: f
            [1 2 5 11]
            sage: g = ntl.ZZ_pX([0,0,0]); g
            []
            sage: g[10]=5
            sage: g
            [0 0 0 0 0 0 0 0 0 0 5]
            sage: g[10]
            5
        """
        return

    def __reduce__(self):
        raise NotImplementedError

    def __dealloc__(self):
        if self.x:
            ZZ_pX_dealloc(self.x)


    def __repr__(self):
        return str(ZZ_pX_repr(self.x))

    def __copy__(self):
        return make_ZZ_pX(ZZ_pX_copy(self.x))

    def copy(self):
        return make_ZZ_pX(ZZ_pX_copy(self.x))

    cdef set(self, void* x):  # only used internally for initialization; assumes self.x not set yet!
        self.x = <ZZ_pX*>x

    def __setitem__(self, long i, a):
        if i < 0:
            raise IndexError, "index (i=%s) must be >= 0"%i
        a = str(int(a))
        ZZ_pX_setitem(self.x, i, a)

    cdef void setitem_from_int(ntl_ZZ_pX self, long i, int value):
        r"""
        Sets ith coefficient to value.

        AUTHOR: David Harvey (2006-08-05)
        """
        ZZ_pX_setitem_from_int(self.x, i, value)

    def setitem_from_int_doctest(self, i, value):
        r"""
        This method exists solely for automated testing of setitem_from_int().

        sage: ntl.set_modulus(ntl.ZZ(20))
        sage: x = ntl.ZZ_pX([2, 3, 4])
        sage: x.setitem_from_int_doctest(5, 42)
        sage: x
         [2 3 4 0 0 2]
        """
        self.setitem_from_int(int(i), int(value))

    def __getitem__(self, unsigned int i):
        cdef char* t
        t = ZZ_pX_getitem(self.x,i)
        return int(string(t))

    cdef int getitem_as_int(ntl_ZZ_pX self, long i):
        r"""
        Returns ith coefficient as C int.
        Return value is only valid if the result fits into an int.

        AUTHOR: David Harvey (2006-08-05)
        """
        return ZZ_pX_getitem_as_int(self.x, i)

    def getitem_as_int_doctest(self, i):
        r"""
        This method exists solely for automated testing of getitem_as_int().

        sage: ntl.set_modulus(ntl.ZZ(20))
        sage: x = ntl.ZZ_pX([2, 3, 5, -7, 11])
        sage: i = x.getitem_as_int_doctest(3)
        sage: print i
         13
        sage: print type(i)
         <type 'int'>
        sage: print x.getitem_as_int_doctest(15)
         0
        """
        return self.getitem_as_int(i)

    def list(self):
        """
        Return list of entries as a list of Python int's.
        """
        return eval(str(self).replace(' ',','))

    def __add__(ntl_ZZ_pX self, ntl_ZZ_pX other):
        """
        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(20))
            sage: ntl.ZZ_pX(range(5)) + ntl.ZZ_pX(range(6))
            [0 2 4 6 8 5]
        """
        return make_ZZ_pX(ZZ_pX_add(self.x, other.x))

    def __sub__(ntl_ZZ_pX self, ntl_ZZ_pX other):
        """
        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(20))
            sage: ntl.ZZ_pX(range(5)) - ntl.ZZ_pX(range(6))
            [0 0 0 0 0 15]
        """
        return make_ZZ_pX(ZZ_pX_sub(self.x, other.x))

    def __mul__(ntl_ZZ_pX self, ntl_ZZ_pX other):
        """
        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(20))
            sage: ntl.ZZ_pX(range(5)) * ntl.ZZ_pX(range(6))
            [0 0 1 4 10 0 10 14 11]
        """
        _sig_on
        return make_ZZ_pX(ZZ_pX_mul(self.x, other.x))

    def __div__(ntl_ZZ_pX self, ntl_ZZ_pX other):
        """
        Compute quotient self / other, if the quotient is a polynomial.
        Otherwise an Exception is raised.

        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(17))
            sage: f = ntl.ZZ_pX([1,2,3]) * ntl.ZZ_pX([4,5])**2
            sage: g = ntl.ZZ_pX([4,5])
            sage: f/g
            [4 13 5 15]
            sage: ntl.ZZ_pX([1,2,3]) * ntl.ZZ_pX([4,5])
            [4 13 5 15]

            sage: f = ntl.ZZ_pX(range(10)); g = ntl.ZZ_pX([-1,0,1])
            sage: f/g
            Traceback (most recent call last):
            ...
            ArithmeticError: self (=[0 1 2 3 4 5 6 7 8 9]) is not divisible by other (=[16 0 1])
        """
        _sig_on
        cdef int divisible
        cdef ZZ_pX* q
        q = ZZ_pX_div(self.x, other.x, &divisible)
        if not divisible:
            raise ArithmeticError, "self (=%s) is not divisible by other (=%s)"%(self, other)
        return make_ZZ_pX(q)

    def __mod__(ntl_ZZ_pX self, ntl_ZZ_pX other):
        """
        Given polynomials a, b in ZZ[X], there exist polynomials q, r
        in QQ[X] such that a = b*q + r, deg(r) < deg(b).  This
        function returns q if q lies in ZZ[X], and otherwise raises an
        Exception.

        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(17))
            sage: f = ntl.ZZ_pX([2,4,6]); g = ntl.ZZ_pX([2])
            sage: f % g   # 0
            []

            sage: f = ntl.ZZ_pX(range(10)); g = ntl.ZZ_pX([-1,0,1])
            sage: f % g
            [3 8]
        """
        _sig_on
        return make_ZZ_pX(ZZ_pX_mod(self.x, other.x))

    def quo_rem(self, ntl_ZZ_pX other):
        """
        Returns the unique integral q and r such that self = q*other +
        r, if they exist.  Otherwise raises an Exception.

        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(17))
            sage: f = ntl.ZZ_pX(range(10)); g = ntl.ZZ_pX([-1,0,1])
            sage: q, r = f.quo_rem(g)
            sage: q, r
            ([3 7 1 4 14 16 8 9], [3 8])
            sage: q*g + r == f
            True
        """
        cdef ZZ_pX *r, *q
        _sig_on
        ZZ_pX_quo_rem(self.x, other.x, &r, &q)
        return (make_ZZ_pX(q), make_ZZ_pX(r))

    def square(self):
        """
        Return f*f.

        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(17))
            sage: f = ntl.ZZ_pX([-1,0,1])
            sage: f*f
            [1 0 15 0 1]
        """
        _sig_on
        return make_ZZ_pX(ZZ_pX_square(self.x))

    def __pow__(ntl_ZZ_pX self, long n, ignored):
        """
        Return the n-th nonnegative power of self.

        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(20))
            sage: g = ntl.ZZ_pX([-1,0,1])
            sage: g**10
            [1 0 10 0 5 0 0 0 10 0 8 0 10 0 0 0 5 0 10 0 1]
        """
        if n < 0:
            raise NotImplementedError
        import sage.rings.arith
        return sage.rings.arith.generic_power(self, n, make_new_ZZ_pX([1]))

    def __cmp__(ntl_ZZ_pX self, ntl_ZZ_pX other):
        """
        Decide whether or not self and other are equal.

        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(20))
            sage: f = ntl.ZZ_pX([1,2,3])
            sage: g = ntl.ZZ_pX([1,2,3,0])
            sage: f == g
            True
            sage: g = ntl.ZZ_pX([0,1,2,3])
            sage: f == g
            False
        """
        if ZZ_pX_equal(self.x, other.x):
            return 0
        return -1

    def is_zero(self):
        """
        Return True exactly if this polynomial is 0.

        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(20))
            sage: f = ntl.ZZ_pX([0,0,0,20])
            sage: f.is_zero()
            True
            sage: f = ntl.ZZ_pX([0,0,1])
            sage: f
            [0 0 1]
            sage: f.is_zero()
            False
        """
        return bool(ZZ_pX_is_zero(self.x))

    def is_one(self):
        """
        Return True exactly if this polynomial is 1.

        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(20))
            sage: f = ntl.ZZ_pX([1,1])
            sage: f.is_one()
            False
            sage: f = ntl.ZZ_pX([1])
            sage: f.is_one()
            True
        """
        return bool(ZZ_pX_is_one(self.x))

    def is_monic(self):
        """
        Return True exactly if this polynomial is monic.

        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(20))
            sage: f = ntl.ZZ_pX([2,0,0,1])
            sage: f.is_monic()
            True
            sage: g = f.reverse()
            sage: g.is_monic()
            False
            sage: g
            [1 0 0 2]
        """
        return bool(ZZ_pX_is_monic(self.x))

    def __neg__(self):
        """
        Return the negative of self.
        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(20))
            sage: f = ntl.ZZ_pX([2,0,0,1])
            sage: -f
            [18 0 0 19]
        """
        return make_ZZ_pX(ZZ_pX_neg(self.x))

    def left_shift(self, long n):
        """
        Return the polynomial obtained by shifting all coefficients of
        this polynomial to the left n positions.

        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(20))
            sage: f = ntl.ZZ_pX([2,0,0,1])
            sage: f
            [2 0 0 1]
            sage: f.left_shift(2)
            [0 0 2 0 0 1]
            sage: f.left_shift(5)
            [0 0 0 0 0 2 0 0 1]

        A negative left shift is a right shift.
            sage: f.left_shift(-2)
            [0 1]
        """
        return make_ZZ_pX(ZZ_pX_left_shift(self.x, n))

    def right_shift(self, long n):
        """
        Return the polynomial obtained by shifting all coefficients of
        this polynomial to the right n positions.

        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(20))
            sage: f = ntl.ZZ_pX([2,0,0,1])
            sage: f
            [2 0 0 1]
            sage: f.right_shift(2)
            [0 1]
            sage: f.right_shift(5)
            []
            sage: f.right_shift(-2)
            [0 0 2 0 0 1]
        """
        return make_ZZ_pX(ZZ_pX_right_shift(self.x, n))

    def gcd(self, ntl_ZZ_pX other):
        """
        Return the gcd d = gcd(a, b), where by convention the leading coefficient
        of d is >= 0.  We uses a multi-modular algorithm.

        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(17))
            sage: f = ntl.ZZ_pX([1,2,3]) * ntl.ZZ_pX([4,5])**2
            sage: g = ntl.ZZ_pX([1,1,1])**3 * ntl.ZZ_pX([1,2,3])
            sage: f.gcd(g)
            [6 12 1]
            sage: g.gcd(f)
            [6 12 1]
        """
        _sig_on
        return make_ZZ_pX(ZZ_pX_gcd(self.x, other.x))

    def xgcd(self, ntl_ZZ_pX other, plain=True):
        """
        Returns r,s,t such that r = s*self + t*other.

        Here r is the resultant of a and b; if r != 0, then this
        function computes s and t such that: a*s + b*t = r; otherwise
        s and t are both 0.

        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(17))
            sage: f = ntl.ZZ_pX([1,2,3]) * ntl.ZZ_pX([4,5])**2
            sage: g = ntl.ZZ_pX([1,1,1])**3 * ntl.ZZ_pX([1,2,3])
            sage: f.xgcd(g)   # nothing since they are not coprime
            ([6 12 1], [15 13 6 8 7 9], [4 13])

        In this example the input quadratic polynomials have a common root modulo 13.
            sage: f = ntl.ZZ_pX([5,0,1])
            sage: g = ntl.ZZ_pX([18,0,1])
            sage: f.xgcd(g)
            ([1], [13], [4])
            """
        cdef ZZ_pX *s, *t, *r
        _sig_on
        if plain:
            ZZ_pX_plain_xgcd(&r, &s, &t, self.x, other.x)
        else:
            ZZ_pX_xgcd(&r, &s, &t, self.x, other.x)
        return (make_ZZ_pX(r), make_ZZ_pX(s), make_ZZ_pX(t))

    def degree(self):
        """
        Return the degree of this polynomial.  The degree of the 0
        polynomial is -1.

        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(20))
            sage: f = ntl.ZZ_pX([5,0,1])
            sage: f.degree()
            2
            sage: f = ntl.ZZ_pX(range(100))
            sage: f.degree()
            99
            sage: f = ntl.ZZ_pX()
            sage: f.degree()
            -1
            sage: f = ntl.ZZ_pX([1])
            sage: f.degree()
            0
        """
        return ZZ_pX_degree(self.x)

    def leading_coefficient(self):
        """
        Return the leading coefficient of this polynomial.

        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(20))
            sage: f = ntl.ZZ_pX([3,6,9])
            sage: f.leading_coefficient()
            9
            sage: f = ntl.ZZ_pX()
            sage: f.leading_coefficient()
            0
        """
        return make_ZZ_p(ZZ_pX_leading_coefficient(self.x))

    def constant_term(self):
        """
        Return the constant coefficient of this polynomial.

        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(20))
            sage: f = ntl.ZZ_pX([3,6,9])
            sage: f.constant_term()
            3
            sage: f = ntl.ZZ_pX()
            sage: f.constant_term()
            0
        """
        cdef char* t
        t = ZZ_pX_constant_term(self.x)
        return int(string(t))

    def set_x(self):
        """
        Set this polynomial to the monomial "x".

        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(20))
            sage: f = ntl.ZZ_pX()
            sage: f.set_x()
            sage: f
            [0 1]
            sage: g = ntl.ZZ_pX([0,1])
            sage: f == g
            True

        Though f and g are equal, they are not the same objects in memory:
            sage: f is g
            False

        """
        ZZ_pX_set_x(self.x)

    def is_x(self):
        """
        True if this is the polynomial "x".

        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(20))
            sage: f = ntl.ZZ_pX()
            sage: f.set_x()
            sage: f.is_x()
            True
            sage: f = ntl.ZZ_pX([0,1])
            sage: f.is_x()
            True
            sage: f = ntl.ZZ_pX([1])
            sage: f.is_x()
            False
        """
        return bool(ZZ_pX_is_x(self.x))

    def derivative(self):
        """
        Return the derivative of this polynomial.

        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(20))
            sage: f = ntl.ZZ_pX([1,7,0,13])
            sage: f.derivative()
            [7 0 19]
        """
        return make_ZZ_pX(ZZ_pX_derivative(self.x))

    def factor(self, verbose=False):
        cdef ZZ_pX** v
        cdef long* e
        cdef long i, n
        _sig_on
        ZZ_pX_factor(&v, &e, &n, self.x, verbose)
        _sig_off
        F = []
        for i from 0 <= i < n:
            F.append((make_ZZ_pX(v[i]), e[i]))
        free(v)
        free(e)
        return F

    def linear_roots(self):
        """
        Assumes that input is monic, and has deg(f) distinct roots.
        Returns the list of roots.
        """
        cdef ZZ_p** v
        cdef long i, n
        _sig_on
        ZZ_pX_linear_roots(&v, &n, self.x)
        _sig_off
        F = []
        for i from 0 <= i < n:
            F.append(make_ZZ_p(v[i]))
        free(v)
        return F

    def reverse(self, hi=None):
        """
        Return the polynomial obtained by reversing the coefficients
        of this polynomial.  If hi is set then this function behaves
        as if this polynomial has degree hi.

        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(20))
            sage: f = ntl.ZZ_pX([1,2,3,4,5])
            sage: f.reverse()
            [5 4 3 2 1]
            sage: f.reverse(hi=10)
            [0 0 0 0 0 0 5 4 3 2 1]
            sage: f.reverse(hi=2)
            [3 2 1]
            sage: f.reverse(hi=-2)
            []
        """
        if not (hi is None):
            return make_ZZ_pX(ZZ_pX_reverse_hi(self.x, int(hi)))
        else:
            return make_ZZ_pX(ZZ_pX_reverse(self.x))

    def truncate(self, long m):
        """
        Return the truncation of this polynomial obtained by
        removing all terms of degree >= m.

        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(20))
            sage: f = ntl.ZZ_pX([1,2,3,4,5])
            sage: f.truncate(3)
            [1 2 3]
            sage: f.truncate(8)
            [1 2 3 4 5]
            sage: f.truncate(1)
            [1]
            sage: f.truncate(0)
            []
            sage: f.truncate(-1)
            []
            sage: f.truncate(-5)
            []
        """
        if m <= 0:
            return make_new_ZZ_pX()
        return make_ZZ_pX(ZZ_pX_truncate(self.x, m))

    def multiply_and_truncate(self, ntl_ZZ_pX other, long m):
        """
        Return self*other but with terms of degree >= m removed.

        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(20))
            sage: f = ntl.ZZ_pX([1,2,3,4,5])
            sage: g = ntl.ZZ_pX([10])
            sage: f.multiply_and_truncate(g, 2)
            [10]
            sage: g.multiply_and_truncate(f, 2)
            [10]
        """
        if m <= 0:
            return make_new_ZZ_pX()
        return make_ZZ_pX(ZZ_pX_multiply_and_truncate(self.x, other.x, m))

    def square_and_truncate(self, long m):
        """
        Return self*self but with terms of degree >= m removed.

        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(20))
            sage: f = ntl.ZZ_pX([1,2,3,4,5])
            sage: f.square_and_truncate(4)
            [1 4 10]
            sage: (f*f).truncate(4)
            [1 4 10]
        """
        if m < 0:
            return make_new_ZZ_pX()
        return make_ZZ_pX(ZZ_pX_square_and_truncate(self.x, m))

    def invert_and_truncate(self, long m):
        """
        Compute and return the inverse of self modulo $x^m$.
        The constant term of self must be 1 or -1.

        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(20))
            sage: f = ntl.ZZ_pX([1,2,3,4,5,6,7])
            sage: f.invert_and_truncate(20)
            [1 18 1 0 0 0 0 8 17 2 13 0 0 0 4 0 17 10 9]
            sage: g = f.invert_and_truncate(20)
            sage: g * f
            [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 8 4 4 3]
        """
        if m < 0:
            raise ArithmeticError, "m (=%s) must be positive"%m
        n = self.constant_term()
        if n != 1 and n != -1:
            raise ArithmeticError, \
                  "The constant term of self must be 1 or -1."
        _sig_on
        return make_ZZ_pX(ZZ_pX_invert_and_truncate(self.x, m))

    def multiply_mod(self, ntl_ZZ_pX other, ntl_ZZ_pX modulus):
        """
        Return self*other % modulus.  The modulus must be monic with
        deg(modulus) > 0, and self and other must have smaller degree.

        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(20))
            sage: modulus = ntl.ZZ_pX([1,2,0,1])    # must be monic
            sage: g = ntl.ZZ_pX([-1,0,1])
            sage: h = ntl.ZZ_pX([3,7,13])
            sage: h.multiply_mod(g, modulus)
            [10 6 4]
        """
        _sig_on
        return make_ZZ_pX(ZZ_pX_multiply_mod(self.x, other.x, modulus.x))

    def trace_mod(self, ntl_ZZ_pX modulus):
        """
        Return the trace of this polynomial modulus the modulus.
        The modulus must be monic, and of positive degree degree bigger
        than the the degree of self.

        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(20))
            sage: f = ntl.ZZ_pX([1,2,0,3])
            sage: mod = ntl.ZZ_pX([5,3,-1,1,1])
            sage: f.trace_mod(mod)
            3
        """
        _sig_on
        return make_ZZ_p(ZZ_pX_trace_mod(self.x, modulus.x))

    def trace_list(self):
        """
        Return the list of traces of the powers $x^i$ of the
        monomial x modulo this polynomial for i = 0, ..., deg(f)-1.
        This polynomial must be monic.

        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(20))
            sage: f = ntl.ZZ_pX([1,2,0,3,0,1])
            sage: f.trace_list()
            [5, 0, 14, 0, 10]

        The input polynomial must be monic or a ValueError is raised:
            sage: ntl.set_modulus(ntl.ZZ(20))
            sage: f = ntl.ZZ_pX([1,2,0,3,0,2])
            sage: f.trace_list()
            Traceback (most recent call last):
            ...
            ValueError: polynomial must be monic.
        """
        if not self.is_monic():
            raise ValueError, "polynomial must be monic."
        _sig_on
        cdef char* t
        t = ZZ_pX_trace_list(self.x)
        return eval(string(t).replace(' ', ','))

    def resultant(self, ntl_ZZ_pX other):
        """
        Return the resultant of self and other.

        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(17))
            sage: f = ntl.ZZ_pX([17,0,1,1])
            sage: g = ntl.ZZ_pX([34,-17,18,2])
            sage: f.resultant(g)
            0
        """
        _sig_on
        return make_ZZ_p(ZZ_pX_resultant(self.x, other.x))

    def norm_mod(self, ntl_ZZ_pX modulus):
        """
        Return the norm of this polynomial modulo the modulus.  The
        modulus must be monic, and of positive degree strictly greater
        than the degree of self.


        EXAMPLE:
            sage: ntl.set_modulus(ntl.ZZ(17))
            sage: f = ntl.ZZ_pX([1,2,0,3])
            sage: mod = ntl.ZZ_pX([-5,2,0,0,1])
            sage: f.norm_mod(mod)
            11

        The norm is the constant term of the characteristic polynomial.
            sage: f.charpoly_mod(mod)
            [11 1 8 14 1]
        """
        _sig_on
        return make_ZZ_p(ZZ_pX_norm_mod(self.x, modulus.x))

    def discriminant(self):
        r"""
        Return the discriminant of a=self, which is by definition
        $$
                (-1)^{m(m-1)/2} {\mbox{\tt resultant}}(a, a')/lc(a),
        $$
        where m = deg(a), and lc(a) is the leading coefficient of a.

        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(17))
            sage: f = ntl.ZZ_pX([1,2,0,3])
            sage: f.discriminant()
            1
        """
        cdef long m

        c = ~self.leading_coefficient()
        m = self.degree()
        if (m*(m-1)/2) % 2:
            c = -c
        return c*self.resultant(self.derivative())

    def charpoly_mod(self, ntl_ZZ_pX modulus):
        """
        Return the characteristic polynomial of this polynomial modulo
        the modulus.  The modulus must be monic of degree bigger than
        self.

        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(17))
            sage: f = ntl.ZZ_pX([1,2,0,3])
            sage: mod = ntl.ZZ_pX([-5,2,0,0,1])
            sage: f.charpoly_mod(mod)
            [11 1 8 14 1]
        """
        _sig_on
        return make_ZZ_pX(ZZ_pX_charpoly_mod(self.x, modulus.x))

    def minpoly_mod(self, ntl_ZZ_pX modulus):
        """
        Return the minimal polynomial of this polynomial modulo the
        modulus.  The modulus must be monic of degree bigger than
        self.

        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(17))
            sage: f = ntl.ZZ_pX([0,0,1])
            sage: g = f*f
            sage: f.charpoly_mod(g)
            [0 0 0 0 1]

        However, since $f^2 = 0$ moduluo $g$, its minimal polynomial
        is of degree $2$.
            sage: f.minpoly_mod(g)
            [0 0 1]
        """
        _sig_on
        return make_ZZ_pX(ZZ_pX_minpoly_mod(self.x, modulus.x))

    def clear(self):
        """
        Reset this polynomial to 0.  Changes this polynomial in place.

        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(17))
            sage: f = ntl.ZZ_pX([1,2,3])
            sage: f
            [1 2 3]
            sage: f.clear()
            sage: f
            []
        """
        ZZ_pX_clear(self.x)

    def preallocate_space(self, long n):
        """
        Pre-allocate spaces for n coefficients.  The polynomial that f
        represents is unchanged.  This is useful if you know you'll be
        setting coefficients up to n, so memory isn't re-allocated as
        the polynomial grows.  (You might save a millisecond with this
        function.)

        EXAMPLES:
            sage: ntl.set_modulus(ntl.ZZ(17))
            sage: f = ntl.ZZ_pX([1,2,3])
            sage: f.preallocate_space(20)
            sage: f
            [1 2 3]
            sage: f[10]=5  # no new memory is allocated
            sage: f
            [1 2 3 0 0 0 0 0 0 0 5]
        """
        _sig_on
        ZZ_pX_preallocate_space(self.x, n)
        _sig_off


## TODO: NTL's ZZ_pX has minpolys of linear recurrence sequences!!!


cdef make_ZZ_pX(ZZ_pX* x):
    cdef ntl_ZZ_pX y
    _sig_off
    y = ntl_ZZ_pX()
    y.x = x
    return y

def make_new_ZZ_pX(v=[]):
    s = str(v).replace(',',' ').replace('L','')
    cdef ntl_ZZ_pX z
    z = ntl_ZZ_pX()
    _sig_on
    z.x = str_to_ZZ_pX(s)
    _sig_off
    return z


##############################################################################
#
# ntl_mat_ZZ: Matrices over the integers via NTL
#
##############################################################################

cdef class ntl_mat_ZZ:
    # see ntl.pxd for data members
    r"""
    The \class{mat_ZZ} class implements arithmetic with matrices over $\Z$.
    """
    def __init__(self, nrows=0,  ncols=0, v=None):
        if nrows == _INIT:
            return
        cdef unsigned long i, j
        cdef ntl_ZZ tmp
        if nrows == 0 and ncols == 0:
            return
        nrows = int(nrows)
        ncols = int(ncols)
        self.x = new_mat_ZZ(nrows, ncols)
        self.__nrows = nrows
        self.__ncols = ncols
        if v != None:
            for i from 0 <= i < nrows:
                for j from 0 <= j < ncols:
                    tmp = make_new_ZZ(v[i*ncols+j])
                    mat_ZZ_setitem(self.x, i, j, <ntl_c_ZZ*> tmp.x)


    def __reduce__(self):
        raise NotImplementedError

    def __dealloc__(self):
        del_mat_ZZ(self.x)

    def __repr__(self):
        _sig_on
        return string(mat_ZZ_to_str(self.x))

    def __mul__(ntl_mat_ZZ self, other):
        cdef ntl_mat_ZZ y
        if not isinstance(other, ntl_mat_ZZ):
            other = ntl_mat_ZZ(other)
        y = other
        _sig_on
        return make_mat_ZZ(mat_ZZ_mul(self.x, y.x))

    def __sub__(ntl_mat_ZZ self, other):
        cdef ntl_mat_ZZ y
        if not isinstance(other, ntl_mat_ZZ):
            other = ntl_mat_ZZ(other)
        y = other
        _sig_on
        return make_mat_ZZ(mat_ZZ_sub(self.x, y.x))

    def __add__(ntl_mat_ZZ self, other):
        cdef ntl_mat_ZZ y
        if not isinstance(other, ntl_mat_ZZ):
            other = ntl_mat_ZZ(other)
        y = other
        _sig_on
        return make_mat_ZZ(mat_ZZ_add(self.x, y.x))

    def __pow__(ntl_mat_ZZ self, long e, ignored):
        _sig_on
        return make_mat_ZZ(mat_ZZ_pow(self.x, e))

    def nrows(self):
        return self.__nrows

    def ncols(self):
        return self.__ncols

    def __setitem__(self, ij, x):
        cdef ntl_ZZ y
        cdef int i, j
        if not isinstance(x, ntl_ZZ):
            y = make_new_ZZ(x)
        else:
            y = x
        if not isinstance(ij, tuple) or len(ij) != 2:
            raise TypeError, 'ij must be a 2-tuple'
        i, j = int(ij[0]),int(ij[1])
        if i < 0 or i >= self.__nrows or j < 0 or j >= self.__ncols:
            raise IndexError, "array index out of range"
        _sig_on
        mat_ZZ_setitem(self.x, i, j, <ntl_c_ZZ*> y.x)
        _sig_off

    def __getitem__(self, ij):
        cdef int i, j
        if not isinstance(ij, tuple) or len(ij) != 2:
            raise TypeError, 'ij must be a 2-tuple'
        i, j = ij
        if i < 0 or i >= self.__nrows or j < 0 or j >= self.__ncols:
            raise IndexError, "array index out of range"
        return make_ZZ(mat_ZZ_getitem(self.x, i+1, j+1))

    def determinant(self, deterministic=True):
        _sig_on
        return make_ZZ(mat_ZZ_determinant(self.x, deterministic))

    def HNF(self, D=None):
        r"""
        The input matrix A=self is an n x m matrix of rank m (so n >=
        m), and D is a multiple of the determinant of the lattice L
        spanned by the rows of A.  The output W is the Hermite Normal
        Form of A; that is, W is the unique m x m matrix whose rows
        span L, such that

        - W is lower triangular,
        - the diagonal entries are positive,
        - any entry below the diagonal is a non-negative number
          strictly less than the diagonal entry in its column.

        This is implemented using the algorithm of [P. Domich,
        R. Kannan and L. Trotter, Math. Oper. Research 12:50-59,
        1987].

        TIMINGS:
        NTL isn't very good compared to MAGMA, unfortunately:

            sage.: import ntl
            sage.: a=MatrixSpace(Q,200).random_element()    # -2 to 2
            sage.: A=ntl.mat_ZZ(200,200)
            sage.: for i in xrange(a.nrows()):
               ....:     for j in xrange(a.ncols()):
               ....:         A[i,j] = a[i,j]
               ....:
            sage.: time d=A.determinant()
            Time.: 3.89 seconds
            sage.: time B=A.HNF(d)
            Time.: 27.59 seconds

        In comparison, MAGMA does this much more quickly:
        \begin{verbatim}
            > A := MatrixAlgebra(Z,200)![Random(-2,2) : i in [1..200^2]];
            > time d := Determinant(A);
            Time: 0.710
            > time H := HermiteForm(A);
            Time: 3.080
        \end{verbatim}

        Also, PARI is also faster than NTL if one uses the flag 1 to
        the mathnf routine.  The above takes 16 seconds in PARI.
        """
        cdef ntl_ZZ _D
        if D == None:
            _D = self.determinant()
        else:
            _D = ntl_ZZ(D)
        _sig_on
        return make_mat_ZZ(mat_ZZ_HNF(self.x, <ntl_c_ZZ*>_D.x))

    def charpoly(self):
        return make_ZZX(mat_ZZ_charpoly(self.x));

    def LLL(self, a=3, b=4, return_U=False, verbose=False):
        r"""
        Performs LLL reduction of self (puts self in an LLL form).

        self is an m x n matrix, viewed as m rows of n-vectors.  m may
        be less than, equal to, or greater than n, and the rows need
        not be linearly independent. self is transformed into an
        LLL-reduced basis, and the return value is the rank r of self
        so as det2 (see below).  The first m-r rows of self are zero.

        More specifically, elementary row transformations are
        performed on self so that the non-zero rows of new-self form
        an LLL-reduced basis for the lattice spanned by the rows of
        old-self.  The default reduction parameter is $\delta=3/4$,
        which means that the squared length of the first non-zero
        basis vector is no more than $2^{r-1}$ times that of the
        shortest vector in the lattice.

        det2 is calculated as the \emph{square} of the determinant of
        the lattice---note that sqrt(det2) is in general an integer
        only when r = n.

        If return_U is True, a value U is returned which is the
        transformation matrix, so that U is a unimodular m x m matrix
        with U * old-self = new-self.  Note that the first m-r rows of
        U form a basis (as a lattice) for the kernel of old-B.

        The parameters a and b allow an arbitrary reduction parameter
        $\delta=a/b$, where $1/4 < a/b \leq 1$, where a and b are positive
        integers.  For a basis reduced with parameter delta, the
        squared length of the first non-zero basis vector is no more
        than $1/(delta-1/4)^{r-1}$ times that of the shortest vector in
        the lattice.

        The algorithm employed here is essentially the one in Cohen's
        book: [H. Cohen, A Course in Computational Algebraic Number
        Theory, Springer, 1993]

        INPUT:
           a        -- parameter a as described above (default: 3)
           b        -- parameter b as described above (default: 4)
           return_U -- return U as described above
           verbose  -- if True NTL will produce some verbatim messages on
                       what's going on internally (default: False)

        OUTPUT:
            (rank,det2,[U]) where rank,det2, and U are as described
            above and U is an optional return value if return_U is
            True.

        EXAMPLE:
            sage: M=ntl.mat_ZZ(3,3,[1,2,3,4,5,6,7,8,9])
            sage: M.LLL()
            (2, 54)
            sage: M
            [[0 0 0]
            [2 1 0]
            [-1 1 3]
            ]
            sage: M=ntl.mat_ZZ(4,4,[-6,9,-15,-18,4,-6,10,12,10,-16,18,35,-24,36,-46,-82]); M
            [[-6 9 -15 -18]
            [4 -6 10 12]
            [10 -16 18 35]
            [-24 36 -46 -82]
            ]
            sage: M.LLL()
            (3, 19140)
            sage: M
            [[0 0 0 0]
            [0 -2 0 0]
            [-2 1 -5 -6]
            [0 -1 -7 5]
            ]

        WARNING: This method modifies self. So after applying this method your matrix
        will be a vector of vectors.
        """
        cdef ntl_c_ZZ *det2
        cdef mat_ZZ *U
        if return_U:
            _sig_on
            U = new_mat_ZZ(self.nrows(),self.ncols())
            rank = int(mat_ZZ_LLL_U(&det2, self.x, U, int(a), int(b), int(verbose)))
            _sig_off
            return rank, make_ZZ(det2), make_mat_ZZ(U)
        else:
            _sig_on
            rank = int(mat_ZZ_LLL(&det2,self.x,int(a),int(b),int(verbose)))
            _sig_off
            return rank,make_ZZ(det2)

cdef make_mat_ZZ(mat_ZZ* x):
    cdef ntl_mat_ZZ y
    _sig_off
    y = ntl_mat_ZZ(_INIT)
    y.x = x
    y.__nrows = mat_ZZ_nrows(y.x);
    y.__ncols = mat_ZZ_ncols(y.x);
    return y


##############################################################################
#
# ntl_GF2X: Polynomials over GF(2) via NTL
#
# AUTHORS:
#  - Martin Albrecht <malb@informatik.uni-bremen.de>
#    2006-01: initial version (based on code by William Stein)
#
##############################################################################

__have_GF2X_hex_repr = False # hex representation of GF2X


cdef class ntl_GF2X:
    """
    Polynomials over GF(2) via NTL
    """
    # See ntl.pxd for definition of data members

    def __reduce__(self):
        raise NotImplementedError

    def __dealloc__(self):
        del_GF2X(self.gf2x_x)

    def __repr__(self):
        _sig_on
        return string(GF2X_to_str(self.gf2x_x))


    def __mul__(ntl_GF2X self, other):
        cdef ntl_GF2X y
        if not isinstance(other, ntl_GF2X):
            other = ntl_GF2X(other)
        y = other
        _sig_on
        return make_GF2X(GF2X_mul(self.gf2x_x, y.gf2x_x))

    def __sub__(ntl_GF2X self, other):
        cdef ntl_GF2X y
        if not isinstance(other, ntl_GF2X):
            other = ntl_GF2X(other)
        y = other
        _sig_on
        return make_GF2X(GF2X_sub(self.gf2x_x, y.gf2x_x))

    def __add__(ntl_GF2X self, other):
        cdef ntl_GF2X y
        if not isinstance(other, ntl_GF2X):
            other = ntl_GF2X(other)
        y = other
        _sig_on
        return make_GF2X(GF2X_add(self.gf2x_x, y.gf2x_x))

    def __neg__(ntl_GF2X self):
        _sig_on
        return make_GF2X(GF2X_neg(self.gf2x_x))

    def __pow__(ntl_GF2X self, long e, ignored):
        _sig_on
        return make_GF2X(GF2X_pow(self.gf2x_x, e))


    def __cmp__(ntl_GF2X self, ntl_GF2X other):
        cdef int t
        _sig_on
        t = GF2X_eq(self.gf2x_x, other.gf2x_x)
        _sig_off
        if t:
            return 0
        return 1

    def degree(ntl_GF2X self):
        """
        Returns the degree of this polynomial
        """
        return GF2X_deg(self.gf2x_x)

    def list(ntl_GF2X self):
        """
        Represents this element as a list of binary digits.

        EXAMPLES:
             sage: e=ntl.GF2X([0,1,1])
             sage: e.list()
             [0, 1, 1]
             sage: e=ntl.GF2X('0xff')
             sage: e.list()
             [1, 1, 1, 1, 1, 1, 1, 1]

        OUTPUT:
             a list of digits representing the coefficients in this element's
             polynomial representation
        """
        #yields e.g. "[1 1 0 0 1 1 0 1]"
        _sig_on
        s = string(GF2X_to_bin(self.gf2x_x))

        #yields e.g. [1,1,0,0,1,1,0,1]
        return map(int,list(s[1:][:len(s)-2].replace(" ","")))


    def bin(ntl_GF2X self):
        """
        Returns binary representation of this element. It is
        the same as setting \code{ntl.hex_output(False)} and
        representing this element afterwards. However it should be
        faster and preserves the HexOutput state as opposed to
        the above code.

        EXAMPLES:
             sage: e=ntl.GF2X([1,1,0,1,1,1,0,0,1])
             sage: e.bin()
             '[1 1 0 1 1 1 0 0 1]'

        OUTPUT:
            string representing this element in binary digits

        """
        _sig_on
        return string(GF2X_to_bin(self.gf2x_x))

    def hex(ntl_GF2X self):
        """
        Returns hexadecimal representation of this element. It is
        the same as setting \code{ntl.hex_output(True)} and
        representing this element afterwards. However it should be
        faster and preserves the HexOutput state as opposed to
        the above code.

        EXAMPLES:
             sage: e=ntl.GF2X([1,1,0,1,1,1,0,0,1])
             sage: e.hex()
             '0xb31'

        OUTPUT:
            string representing this element in hexadecimal

        """
        _sig_on
        return string(GF2X_to_hex(self.gf2x_x))

    def _sage_(ntl_GF2X self,R=None,cache=None):
        """
        Returns a SAGE polynomial over GF(2) equivalent to
        this element. If a ring R is provided it is used
        to construct the polynomial in, otherwise
        an appropriate ring is generated.

        INPUT:
            self  -- GF2X element
            R     -- PolynomialRing over GF(2)
            cache -- optional NTL to SAGE cache (dict)

        OUTPUT:
            polynomial in R
        """
        if R==None:
            from sage.rings.polynomial_ring import PolynomialRing
            from sage.rings.finite_field import FiniteField
            R = PolynomialRing(FiniteField(2), 'x')

        if cache != None:
            try:
                return cache[self.hex()]
            except KeyError:
                cache[self.hex()] = R(self.list())
                return cache[self.hex()]

        return R(self.list())

    cdef set(self, void *y):  # only used internally for initialization; assumes self.gf2x_x not set yet!
        self.gf2x_x = <GF2X*> y

cdef public make_GF2X(GF2X* x):
    cdef ntl_GF2X y
    _sig_off
    y = ntl_GF2X()
    y.gf2x_x = x
    return y

def make_new_GF2X(x=[]):
    """
    Constructs a new polynomial over GF(2).

    A value may be passed to this constructor. If you pass a string
    to the constructor please note that byte sequences and the hexadecimal
    notation are little endian.  So e.g. '[0 1]' == '0x2' == x.

    Input types are ntl.ZZ_px, strings, lists of digits, FiniteFieldElements
    from extension fields over GF(2), Polynomials over GF(2), Integers, and finite
    extension fields over GF(2) (uses modulus).

    INPUT:
        x -- value to be assigned to this element. See examples.

    OUTPUT:
        a new ntl.GF2X element

    EXAMPLES:
        sage: ntl.GF2X(ntl.ZZ_pX([1,1,3]))
        [1 1 1]
        sage: ntl.GF2X('0x1c')
        [1 0 0 0 0 0 1 1]
        sage: ntl.GF2X('[1 0 1 0]')
        [1 0 1]
        sage: ntl.GF2X([1,0,1,0])
        [1 0 1]
        sage: ntl.GF2X(GF(2**8,'a').gen()**20)
        [0 0 1 0 1 1 0 1]
        sage: ntl.GF2X(GF(2**8,'a'))
        [1 0 1 1 1 0 0 0 1]
        sage: ntl.GF2X(2)
        [0 1]
    """
    from sage.rings.finite_field_element import FiniteField_ext_pariElement
    from sage.rings.finite_field import FiniteField_ext_pari
    from sage.rings.finite_field_givaro import FiniteField_givaro,FiniteField_givaroElement
    from sage.rings.polynomial_element_generic import Polynomial_dense_mod_p
    from sage.rings.integer import Integer

    if isinstance(x, Integer):
        #binary repr, reversed, and "["..."]" added
        x="["+x.binary()[::-1].replace(""," ")+"]"
    elif type(x) == int:
        #hex repr, "0x" stripped, reversed (!)
        x="0x"+hex(x)[2:][::-1]
    elif isinstance(x, Polynomial_dense_mod_p):
        if x.base_ring().characteristic():
            x=x._Polynomial_dense_mod_n__poly
    elif isinstance(x, (FiniteField_ext_pari,FiniteField_givaro)):
        if x.characteristic() == 2:
            x= list(x.modulus())
    elif isinstance(x, FiniteField_ext_pariElement):
        if x.parent().characteristic() == 2:
            x=x._pari_().centerlift().centerlift().subst('a',2).int_unsafe()
            x="0x"+hex(x)[2:][::-1]
    elif isinstance(x, FiniteField_givaroElement):
        x = "0x"+hex(int(x))[2:][::-1]
    s = str(x).replace(","," ")
    cdef ntl_GF2X n
    n = ntl_GF2X()
    _sig_on
    n.gf2x_x = str_to_GF2X(s)
    _sig_off
    return n



##############################################################################
#
# ntl_GF2E: GF(2**x) via NTL
#
# AUTHORS:
#  - Martin Albrecht <malb@informatik.uni-bremen.de>
#    2006-01: initial version (based on cody by William Stein)
#
##############################################################################

def GF2X_hex_repr(have_hex=None):
    """
    Represent GF2X and GF2E elements in the more compact
    hexadecimal form to the user.

    If no parameter is provided the currently set value will be
    returned.

    INPUT:
        have_hex if True hex representation will be used
    """
    global __have_GF2X_hex_repr

    if have_hex==None:
        return __have_GF2X_hex_repr

    if have_hex==True:
        GF2X_hex(1)
    else:
        GF2X_hex(0)
    __have_GF2X_hex_repr=have_hex

def ntl_GF2E_modulus(p=None):
    """
    Initializes the current modulus to P; required: deg(P) >= 1

    The input is either ntl.GF2X or is tried to be converted to a
    ntl.GF2X element.

    If no parameter p is given: Yields copy of the current GF2E
    modulus.

    INPUT:
        p -- modulus

    EXAMPLES:
        sage: ntl.GF2E_modulus([1,1,0,1,1,0,0,0,1])
        sage: ntl.GF2E_modulus().hex()
        '0xb11'
    """
    global __have_GF2E_modulus
    cdef ntl_GF2X elem

    if p==None:
        if __have_GF2E_modulus == True:
            return make_GF2X(GF2E_modulus())
        else:
            raise "NoModulus"

    if not isinstance(p,ntl_GF2X):
        elem = make_new_GF2X(p)
    else:
        elem = p

    if(elem.degree()<1):
        raise "DegreeToSmall"

    ntl_GF2E_set_modulus(<GF2X*>elem.gf2x_x)
    __have_GF2E_modulus=True

def ntl_GF2E_modulus_degree():
    """
    Returns deg(modulus) for GF2E elements
    """
    if __have_GF2E_modulus:
        return GF2E_degree()
    else:
        raise "NoModulus"

def ntl_GF2E_sage(names='a'):
    """
    Returns a SAGE FiniteField element matching the current modulus.

    EXAMPLES:
        sage: ntl.GF2E_modulus([1,1,0,1,1,0,0,0,1])
        sage: ntl.GF2E_sage()
        Finite Field in a of size 2^8
    """
    from sage.rings.finite_field import FiniteField_ext_pari
    f = ntl_GF2E_modulus()._sage_()
    return FiniteField_ext_pari(int(2)**GF2E_degree(),modulus=f,name=names)


def ntl_GF2E_random():
    """
    Returns a random element from GF2E modulo the current modulus.
    """
    _sig_on
    return make_GF2E(GF2E_random());


# make sure not to segfault
__have_GF2E_modulus = False


cdef class ntl_GF2E(ntl_GF2X):
    r"""
    The \\class{GF2E} represents a finite extension field over GF(2) using NTL.
    Elements are represented as polynomials over GF(2) modulo \\code{ntl.GF2E_modulus()}.

    This modulus must be set using \\code{ ntl.GF2E_modulus(p) } and is unique for
    alle elements in ntl.GF2E. So it is not possible at the moment e.g. to have elements
    in GF(2**4) and GF(2**8) at the same time. You might however be lucky and get away
    with not touch the elements in GF(2**4) while having the GF(2**8) modulus set and vice
    versa.
    """

    # See ntl.pxd for definition of data members
    def __reduce__(self):
        raise NotImplementedError

    def __dealloc__(self):
        del_GF2E(self.gf2e_x)

    def __repr__(self):
        _sig_on
        return string(GF2E_to_str(self.gf2e_x))


    def __mul__(ntl_GF2E self, other):
        cdef ntl_GF2E y
        if not isinstance(other, ntl_GF2E):
            other = ntl_GF2E(other)
        y = other
        _sig_on
        return make_GF2E(GF2E_mul(self.gf2e_x, y.gf2e_x))

    def __sub__(ntl_GF2E self, other):
        cdef ntl_GF2E y
        if not isinstance(other, ntl_GF2E):
            other = ntl_GF2E(other)
        y = other
        _sig_on
        return make_GF2E(GF2E_sub(self.gf2e_x, y.gf2e_x))

    def __add__(ntl_GF2E self, other):
        cdef ntl_GF2E y
        if not isinstance(other, ntl_GF2E):
            other = ntl_GF2E(other)
        y = other
        _sig_on
        return make_GF2E(GF2E_add(self.gf2e_x, y.gf2e_x))

    def __neg__(ntl_GF2E self):
        _sig_on
        return make_GF2E(GF2E_neg(self.gf2e_x))

    def __pow__(ntl_GF2E self, long e, ignored):
        _sig_on
        return make_GF2E(GF2E_pow(self.gf2e_x, e))

    def __cmp__(ntl_GF2E self, ntl_GF2E other):
        cdef int t
        _sig_on
        t = GF2E_eq(self.gf2e_x, other.gf2e_x)
        _sig_off
        if t:
            return 0
        return 1

    def is_zero(ntl_GF2E self):
        """
        Returns True if this element equals zero, False otherwise.
        """
        return bool(GF2E_is_zero(self.gf2e_x))

    def is_one(ntl_GF2E self):
        """
        Returns True if this element equals one, False otherwise.
        """
        return bool(GF2E_is_one(self.gf2e_x))

    def __copy__(self):
        return make_GF2E(GF2E_copy(self.gf2e_x))

    def copy(ntl_GF2E self):
        """
        Returns a copy of this element.
        """
        return make_GF2E(GF2E_copy(self.gf2e_x))

    def trace(ntl_GF2E self):
        """
        Returns the trace of this element.
        """
        return GF2E_trace(self.gf2e_x)

    def ntl_GF2X(ntl_GF2E self):
        """
        Returns a ntl.GF2X copy of this element.
        """
        return make_GF2X(self.gf2x_x)


    def _sage_(ntl_GF2E self, k=None, cache=None):
        """
        Returns a \class{FiniteFieldElement} representation
        of this element. If a \class{FiniteField} k is provided
        it is constructed in this field if possible. A \class{FiniteField}
        will be constructed if none is provided.

        INPUT:
            self  -- \class{GF2E} element
            k     -- optional GF(2**deg)
            cache -- optional NTL to SAGE conversion dictionary

        OUTPUT:
            FiniteFieldElement over k
        """
        cdef int i
        cdef int length
        deg= GF2E_degree()

        if k==None:
            from sage.rings.finite_field import FiniteField_ext_pari
            f = ntl_GF2E_modulus()._sage_()
            k = FiniteField_ext_pari(2**deg,modulus=f)

        if cache != None:
            try:
                return cache[self.hex()]
            except KeyError:
                pass



        a=k.gen()
        l = self.list()

        length = len(l)
        ret = 0

        for i from 0 <= i < length:
            if l[i]==1:
                ret = ret + a**i

        if cache != None:
            cache[self.hex()] = ret

        return ret


    cdef set(self, void *y):  # only used internally for initialization; assumes self.gf2e_x not set yet!
        self.gf2e_x = <GF2E*> y

cdef public make_GF2E(GF2E* x):
    cdef ntl_GF2E y
    _sig_off
    y = ntl_GF2E()
    y.gf2e_x = x
    y.gf2x_x = GF2E_ntl_GF2X(y.gf2e_x)
    return y

def make_new_GF2E(x=[]):
    """
    Constructs a new  finite field element in GF(2**x).

    A modulus \emph{must} have been set with \code{ntl.GF2E_modulus(ntl.GF2X(<value>))}
    when calling this constructor.  A value may be passed to this constructor. If you pass a string
    to the constructor please note that byte sequences and the hexadecimal notation are Little Endian in NTL.
    So e.g. '[0 1]' == '0x2' == x.

    INPUT:
        x -- value to be assigned to this element. Same types as ntl.GF2X() are accepted.

    OUTPUT:
        a new ntl.GF2E element

    EXAMPLES:
        sage: m=ntl.GF2E_modulus(ntl.GF2X([1,1,0,1,1,0,0,0,1]))
        sage: ntl.GF2E(ntl.ZZ_pX([1,1,3]))
        [1 1 1]
        sage: ntl.GF2E('0x1c')
        [1 0 0 0 0 0 1 1]
        sage: ntl.GF2E('[1 0 1 0]')
        [1 0 1]
        sage: ntl.GF2E([1,0,1,0])
        [1 0 1]
        sage: ntl.GF2E(GF(2**8,'a').gen()**20)
        [0 0 1 0 1 1 0 1]
    """
    if not __have_GF2E_modulus:
        raise "NoModulus"

    s = str(make_new_GF2X(x))
    cdef ntl_GF2E n
    n = ntl_GF2E()
    _sig_on
    n.gf2e_x = str_to_GF2E(s)
    n.gf2x_x = GF2E_ntl_GF2X(n.gf2e_x)
    _sig_off
    return n


##############################################################################
#
# ntl_GF2EX: Polynomials over GF(2) via NTL
#
# AUTHORS:
#  - Martin Albrecht <malb@informatik.uni-bremen.de> 2006-01: initial version
#
##############################################################################


cdef class ntl_GF2EX:
    r"""
    """
    # See ntl.pxd for definition of data members

    def __reduce__(self):
        raise NotImplementedError

    def __dealloc__(self):
        del_GF2EX(self.x)

    def __repr__(self):
        _sig_on
        return string(GF2EX_to_str(self.x))


    def __mul__(ntl_GF2EX self, other):
        cdef ntl_GF2EX y
        if not isinstance(other, ntl_GF2EX):
            other = ntl_GF2EX(other)
        y = other
        _sig_on
        return make_GF2EX(GF2EX_mul(self.x, y.x))

    def __sub__(ntl_GF2EX self, other):
        cdef ntl_GF2EX y
        if not isinstance(other, ntl_GF2EX):
            other = ntl_GF2EX(other)
        y = other
        _sig_on
        return make_GF2EX(GF2EX_sub(self.x, y.x))

    def __add__(ntl_GF2EX self, other):
        cdef ntl_GF2EX y
        if not isinstance(other, ntl_GF2EX):
            other = ntl_GF2EX(other)
        y = other
        _sig_on
        return make_GF2EX(GF2EX_add(self.x, y.x))

    def __neg__(ntl_GF2EX self):
        _sig_on
        return make_GF2EX(GF2EX_neg(self.x))

    def __pow__(ntl_GF2EX self, long e, ignored):
        _sig_on
        return make_GF2EX(GF2EX_pow(self.x, e))

    cdef set(self, void *y):  # only used internally for initialization; assumes self.x not set yet!
        self.x = <GF2EX*> y

cdef public make_GF2EX(GF2EX* x):
    cdef ntl_GF2EX y
    _sig_off
    y = ntl_GF2EX()
    y.x = x
    return y

def make_new_GF2EX(x='[]'):
    s = str(x)
    cdef ntl_GF2EX n
    n = ntl_GF2EX()
    _sig_on
    n.x = str_to_GF2EX(s)
    _sig_off
    return n


##############################################################################
#
# ntl_mat_GF2E: Matrices over the GF(2**x) via NTL
#
# AUTHORS:
#  - Martin Albrecht <malb@informatik.uni-bremen.de>
#    2006-01: initial version (based on cody by William Stein)
#
##############################################################################


cdef class ntl_mat_GF2E:
    r"""
    The \class{mat_GF2E} class implements arithmetic with matrices over $GF(2**x)$.
    """
    def __init__(self, nrows=0, ncols=0, v=None):
        """
        Constructs a matrix over ntl.GF2E.

        INPUT:
            nrows -- number of rows
            ncols -- nomber of columns
            v     -- either list or Matrix over GF(2**x)

        EXAMPLES:
            sage: ntl.GF2E_modulus([1,1,0,1,1,0,0,0,1])
            sage: m=ntl.mat_GF2E(10,10)
            sage: m=ntl.mat_GF2E(Matrix(GF(2**8, 'a'),10,10))
            sage: m=ntl.mat_GF2E(10,10,[ntl.GF2E_random() for x in xrange(10*10)])
        """

        if nrows is _INIT:
            return

        cdef unsigned long _nrows, _ncols

        import sage.matrix.matrix
        if sage.matrix.matrix.is_Matrix(nrows):
            _nrows = nrows.nrows()
            _ncols = nrows.ncols()
            v     = nrows.list()
        else:
            _nrows = nrows
            _ncols = ncols

        if not __have_GF2E_modulus:
            raise "NoModulus"
        cdef unsigned long i, j
        cdef ntl_GF2E tmp


        self.x = new_mat_GF2E(_nrows, _ncols)
        self.__nrows = _nrows
        self.__ncols = _ncols
        if v != None:
            _sig_on
            for i from 0 <= i < _nrows:
                for j from 0 <= j < _ncols:
                    elem = v[i*_ncols+j]
                    if not isinstance(elem, ntl_GF2E):
                        tmp=make_new_GF2E(elem)
                    else:
                        tmp=elem
                    mat_GF2E_setitem(self.x, i, j, <GF2E*> tmp.gf2e_x)
            _sig_off

    def __reduce__(self):
        raise NotImplementedError


    def __dealloc__(self):
        del_mat_GF2E(self.x)

    def __repr__(self):
        _sig_on
        return string(mat_GF2E_to_str(self.x))

    def __mul__(ntl_mat_GF2E self, other):
        cdef ntl_mat_GF2E y
        if not isinstance(other, ntl_mat_GF2E):
            other = ntl_mat_GF2E(other)
        y = other
        _sig_on
        return make_mat_GF2E(mat_GF2E_mul(self.x, y.x))

    def __sub__(ntl_mat_GF2E self, other):
        cdef ntl_mat_GF2E y
        if not isinstance(other, ntl_mat_GF2E):
            other = ntl_mat_GF2E(other)
        y = other
        _sig_on
        return make_mat_GF2E(mat_GF2E_sub(self.x, y.x))

    def __add__(ntl_mat_GF2E self, other):
        cdef ntl_mat_GF2E y
        if not isinstance(other, ntl_mat_GF2E):
            other = ntl_mat_GF2E(other)
        y = other
        _sig_on
        return make_mat_GF2E(mat_GF2E_add(self.x, y.x))

    def __pow__(ntl_mat_GF2E self, long e, ignored):
        _sig_on
        return make_mat_GF2E(mat_GF2E_pow(self.x, e))

    def nrows(self):
        """
        Number of rows
        """
        return self.__nrows

    def ncols(self):
        """
        Number of columns
        """
        return self.__ncols

    def __setitem__(self, ij, x):
        cdef ntl_GF2E y
        cdef int i, j
        if not isinstance(x, ntl_GF2E):
            x = make_new_GF2E(x)
        y = x

        from sage.rings.integer import Integer
        if isinstance(ij, tuple) and len(ij) == 2:
            i, j = ij
        elif self.ncols()==1 and (isinstance(ij, Integer) or type(ij)==int):
            i, j = ij,0
        elif self.nrows()==1 and (isinstance(ij, Integer) or type(ij)==int):
            i, j = 0,ij
        else:
            raise TypeError, 'ij is not a matrix index'

        if i < 0 or i >= self.__nrows or j < 0 or j >= self.__ncols:
            raise IndexError, "array index out of range"
        mat_GF2E_setitem(self.x, i, j, <GF2E*> y.gf2e_x)

    def __getitem__(self, ij):
        cdef int i, j
        from sage.rings.integer import Integer
        if isinstance(ij, tuple) and len(ij) == 2:
            i, j = ij
        elif self.ncols()==1 and (isinstance(ij, Integer) or type(ij)==int):
            i, j = ij,0
        elif self.nrows()==1 and (isinstance(ij, Integer) or type(ij)==int):
            i, j = 0,ij
        else:
            raise TypeError, 'ij is not a matrix index'
        if i < 0 or i >= self.__nrows or j < 0 or j >= self.__ncols:
            raise IndexError, "array index out of range"
        return make_GF2E(mat_GF2E_getitem(self.x, i+1, j+1))

    def determinant(self):
        """
        Returns the determinant.
        """
        _sig_on
        return make_GF2E(mat_GF2E_determinant(self.x))

    def echelon_form(self,ncols=0):
        """
        Performs unitary row operations so as to bring this matrix into row echelon
        form.  If the optional argument \code{ncols} is supplied, stops when first ncols
        columns are in echelon form.  The return value is the rank (or the
        rank of the first ncols columns).

        INPUT:
           ncols -- number of columns to process

        TODO: what is the output; does it return the rank?
        """
        cdef long w
        w = ncols
        _sig_on
        return mat_GF2E_gauss(self.x, w)

    def list(self):
        """
        Returns a list of the entries in this matrix

        EXAMPLES:
            sage: ntl.GF2E_modulus([1,1,0,1,1,0,0,0,1])
            sage: m = ntl.mat_GF2E(2,2,[ntl.GF2E_random() for x in xrange(2*2)])
            sage: m.list()                   # random output
            [[1 0 1 0 1 0 1 1], [0 1 0 0 0 0 1], [0 0 0 1 0 0 1], [1 1 1 0 0 0 0 1]]
        """
        cdef unsigned long i, j
        v = []
        for i from 0 <= i < self.__nrows:
            for j from 0 <= j < self.__ncols:
                v.append(self[i,j])
        return v

    def is_zero(self):
        cdef long isZero
        _sig_on
        isZero = mat_GF2E_is_zero(self.x)
        _sig_off
        if isZero==0:
            return False
        else:
            return True

    def _sage_(ntl_mat_GF2E self, k=None, cache=None):
        """
        Returns a \class{Matrix} over a FiniteField representation
        of this element.

        INPUT:
            self  -- \class{mat_GF2E} element
            k     -- optional GF(2**deg)
            cache -- optional NTL to SAGE conversion dictionary

        OUTPUT:
            Matrix over k
        """
        if k==None:
            k = ntl_GF2E_sage()
        l = []
        for elem in self.list():
            l.append(elem._sage_(k, cache))

        from sage.matrix.constructor import Matrix
        return Matrix(k,self.nrows(),self.ncols(),l)

    def transpose(ntl_mat_GF2E self):
        """
        Returns the transposed matrix of self.

        OUTPUT:
            transposed Matrix
        """
        _sig_on
        return make_mat_GF2E(mat_GF2E_transpose(self.x))

cdef make_mat_GF2E(mat_GF2E* x):
    cdef ntl_mat_GF2E y
    _sig_off
    y = ntl_mat_GF2E(_INIT)
    y.x = x
    y.__nrows = mat_GF2E_nrows(y.x);
    y.__ncols = mat_GF2E_ncols(y.x);
    return y
