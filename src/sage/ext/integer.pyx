r"""
Elements of the ring $\Z$ of integers
"""

#*****************************************************************************
#       Copyright (C) 2004 William Stein <wstein@ucsd.edu>
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

doc="""
Integers
"""

import operator

include "gmp.pxi"
include "interrupt.pxi"  # ctrl-c interrupt block support

cdef class Integer(element.EuclideanDomainElement)

import sage.rings.integer_ring
import sage.rings.coerce
import sage.rings.infinity
import sage.rings.complex_field
import rational as rational
import sage.libs.pari.all
import mpfr

cdef public int set_mpz(Integer self, mpz_t value):
    mpz_set(self.value, value)

cdef set_from_Integer(Integer self, Integer other):
    mpz_set(self.value, other.value)

cdef set_from_int(Integer self, int other):
    mpz_set_si(self.value, other)

cdef public mpz_t* get_value(Integer self):
    return &self.value

# This crashes SAGE:
#  s = 2003^100300000
# The problem is related to realloc moving all the memory
# and returning a pointer to the new block of memory, I think.

cdef extern from "stdlib.h":
    void abort()

cdef void* pymem_realloc(void *ptr, size_t old_size, size_t new_size):
    return PyMem_Realloc(ptr, new_size)
    #print "hello-realloc: ", <long> ptr, old_size, new_size
    #cdef void* p
    #p = PyMem_Realloc(ptr, new_size)
    #print "p = ", <long> p
    #if <void*> p != <void*> ptr:
    #    abort()
    #return p


cdef void pymem_free(void *ptr, size_t size):
    PyMem_Free(ptr)

cdef void* pymem_malloc(size_t size):
    return PyMem_Malloc(size)

cdef extern from "gmp.h":
    void mp_set_memory_functions (void *(*alloc_func_ptr) (size_t), void *(*realloc_func_ptr) (void *, size_t, size_t), void (*free_func_ptr) (void *, size_t))


def pmem_malloc():
    """
    Use the Python heap for GMP memory management.
    """
    mp_set_memory_functions(PyMem_Malloc, pymem_realloc, pymem_free)
    #mp_set_memory_functions(pymem_malloc, pymem_realloc, pymem_free)

pmem_malloc()

cdef class Integer(element.EuclideanDomainElement):
    """
    The \\class{Integer} class represents arbitrary precision
    integers.  It derives from the \\class{Element} class, so
    integers can be used as ring elements anywhere in SAGE.

    \\begin{notice}
    The class \\class{Integer} is implemented in Pyrex, as a wrapper
    of the GMP \\code{mpz_t} integer type.
    \\end{notice}
    """
    def __new__(self, x=None):
        mpz_init(self.value)

    def __pyxdoc__init__(self):
        """
        You can create an integer from an int, long, string literal, or
        integer modulo N.

        EXAMPLES:
            sage: Integer(495)
            495
            sage: Integer('495949209809328523')
            495949209809328523
            sage: Integer(Mod(3,7))
            3

        Integers also support the standard arithmetic operations, such
        as +,-,*,/,^, \code{abs}, \code{mod}, \code{float}:
            sage: 2^3
            8
        """
    def __init__(self, x=None):
        if not (x is None):
            if isinstance(x, Integer):
                set_from_Integer(self, x)

            elif hasattr(x, "_integer_"):
                set_from_Integer(self, x._integer_())

            elif isinstance(x, int):
                mpz_set_si(self.value, x)

            elif isinstance(x, long):
                s = "%x"%x
                mpz_set_str(self.value, s, 16)

            elif isinstance(x, str):
                if mpz_set_str(self.value, x, 0) != 0:
                    raise TypeError, "unable to convert x (=%s) to an integer"%x

            elif isinstance(x, rational.Rational):
                if x.denominator() != 1:
                    raise TypeError, "Unable to coerce rational (=%s) to an Integer."%x
                set_from_Integer(self, x.numer())


            elif isinstance(x, sage.libs.pari.all.gen):
                if x.type() == 't_INTMOD':
                    x = x.lift()
                # TODO: figure out how to convert to pari integer in base 16 ?
                s = str(x)
                if mpz_set_str(self.value, s, 0) != 0:
                    raise TypeError, "Unable to coerce PARI %s to an Integer."%x

            else:
                raise TypeError, "Unable to coerce %s (of type %s) to an Integer."%(x,type(x))


    def __reduce__(self):
        # This single line below took me HOURS to figure out.
        # It is the *trick* needed to pickle pyrex extension types.
        # The trick is that you must put a pure Python function
        # as the first argument, and that function must return
        # the result of unpickling with the argument in the second
        # tuple as input. All kinds of problems happen
        # if we don't do this.
        return sage.rings.integer.make_integer, (self.str(32),)

    def _reduce_set(self, s):
        mpz_set_str(self.value, s, 32)

    def _im_gens_(self, codomain, im_gens):
        return codomain._coerce_(self)

    #def __reduce__(self):
    #    return sage.rings.integer_ring.Z, (long(self),)

    cdef int cmp(self, Integer x):
        cdef int i
        i = mpz_cmp(self.value, x.value)
        if i < 0:
            return -1
        elif i == 0:
            return 0
        else:
            return 1

    def __richcmp__(Integer self, right, int op):
        cdef int n
        if not isinstance(right, Integer):
            try:
                n = sage.rings.coerce.cmp(self, right)
            except TypeError:
                n = -1
        else:
            n = self.cmp(right)
        return self._rich_to_bool(op, n)

    def __cmp__(self, x):
        return self.cmp(x)


    def copy(self):
        """
        Return a copy of the integer.
        """
        cdef Integer z
        z = Integer()
        set_mpz(z,self.value)
        return z

    def  __dealloc__(self):
        mpz_clear(self.value)

    def __repr__(self):
        return self.str()

    def __str_malloc(self, int base=10):
        """
        Return the string representation of \\code{self} in the given
        base.  (Uses malloc then PyMem.  This is actually slightly
        faster than self.str() below, but it is unpythonic to use
        malloc.)  However, self.str() below is nice because we know
        the size of the string ahead of time, and can work around a
        bug in GMP nicely.  There seems to be a bug in GMP, where
        non-2-power base conversion for very large integers > 10
        million digits (?) crashes GMP.
        """
        _sig_on
        cdef char *s
        s = mpz_get_str(NULL, base, self.value)
        t = str(s)
        free(s)
        _sig_off
        return t

    def str(self, int base=10):
        r"""
        Return the string representation of \code{self} in the given
        base.

        \begin{notice}
        String representation of integers with more than about {\em four
        million digits} is not available in bases other than a power of
        2. This is because of a bug in GMP.  To obtain the base-$b$
        expansion of an integer $n$ use \code{n.str(b)}.
        \end{notice}
        """
        if base < 2 or base > 36:
            raise ValueError, "base (=%s) must be between 2 and 36"%base
        cdef size_t n
        cdef char *s
        n = mpz_sizeinbase(self.value, base) + 2
        if n > 4200000 and base != 2 and base != 4 and base != 8 and base != 16 and base != 32:
            raise RuntimeError, "String representation of integers with more than 4200000 digits is not available in bases other than a power of 2. This is because of a bug in GMP.  To obtain the base-b expansion of an integer n use n.str(b)."
        s = <char *>PyMem_Malloc(n)
        if s == NULL:
            raise MemoryError, "Unable to allocate enough memory for the string representation of an integer."
        _sig_on
        mpz_get_str(s, base, self.value)
        _sig_off
        k = PyString_FromString(s)
        PyMem_Free(s)
        return k

    def __hex__(self):
        r"""
        Return the hexadecimal digits of self in lower case.

        \note{'0x' is \emph{not} prepended to the result like is done
        by the corresponding Python function on int or long.  This is
        for efficiency sake---adding and stripping the string wastes
        time; since this function is used for conversions from
        integers to other C-library structures, it is important that
        it be fast.}

        EXAMPLES:
            sage: print hex(Integer(15))
            f
            sage: print hex(Integer(16))
            10
            sage: print hex(Integer(16938402384092843092843098243))
            36bb1e3929d1a8fe2802f083
            sage: print hex(long(16938402384092843092843098243))
            0x36BB1E3929D1A8FE2802F083L
        """
        return self.str(16)

    def binary(self):
        """
        Return the binary digits of self as a string.

        EXAMPLES:
            sage: print Integer(15).binary()
            1111
            sage: print Integer(16).binary()
            10000
            sage: print Integer(16938402384092843092843098243).binary()
            1101101011101100011110001110010010100111010001101010001111111000101000000000101111000010000011
        """
        return self.str(2)



    def set_si(self, signed long int n):
        """
        Coerces $n$ to a C signed integer if possible, and sets self
        equal to $n$.

        """
        mpz_set_si(self.value, n)

    def set_str(self, s, base=10):
        """
        Set self equal to the number defined by the string s in the
        given base.
        """
        valid = mpz_set_str(self.value, s, base)
        if valid != 0:
            raise TypeError, "unable to convert x (=%s) to an integer"%s

    cdef void set_from_mpz(Integer self, mpz_t value):
        mpz_set(self.value, value)

    cdef mpz_t* get_value(Integer self):
        return &self.value

    def __add_(Integer self, Integer other):
        cdef Integer x
        x = Integer()
        mpz_add(x.value, self.value, other.value)
        return x

    def __add__(x, y):
        if isinstance(x, Integer) and isinstance(y, Integer):
            return x.__add_(y)
        return sage.rings.coerce.bin_op(x, y, operator.add)

    def __sub_(Integer self, Integer other):
        cdef Integer x
        x = Integer()
        mpz_sub(x.value, self.value, other.value)
        return x

    def __sub__(x, y):
        """
        EXAMPLES:
            sage: a = Integer(5); b = Integer(3)
            sage: b - a
            -2
        """
        if isinstance(x, Integer) and isinstance(y, Integer):
            return x.__sub_(y)
        return sage.rings.coerce.bin_op(x, y, operator.sub)

    def __mul_(Integer self, Integer other):
        cdef Integer x
        x = Integer()

        _sig_on
        mpz_mul(x.value, self.value, other.value)
        _sig_off

        return x

    def __mul__(x, y):
        if isinstance(x, Integer) and isinstance(y, Integer):
            return x.__mul_(y)
        if isinstance(x, list):
            return x * int(y)
        return sage.rings.coerce.bin_op(x, y, operator.mul)

    def __div_(Integer self, Integer other):
        cdef Integer x
        #if mpz_divisible_p(self.value, other.value):
        #    x = Integer()
        #    mpz_divexact(x.value, self.value, other.value)
        #    return x
        #else:
        return rational.Rational(self)/rational.Rational(other)
        #raise ArithmeticError, "Exact division impossible."

    def __div__(x, y):
        if isinstance(x, Integer) and isinstance(y, Integer):
            return x.__div_(y)
        return sage.rings.coerce.bin_op(x, y, operator.div)

    def __floordiv(Integer self, Integer other):
        cdef Integer x
        x = Integer()


        _sig_on
        mpz_fdiv_q(x.value, self.value, other.value)
        _sig_off


        return x


    def __floordiv__(x, y):
        if isinstance(x, Integer) and isinstance(y, Integer):
            return x.__floordiv(y)
        return sage.rings.coerce.bin_op(x, y, operator.floordiv)


    def __pow__(self, n, dummy):
        cdef Integer _self, _n
        cdef unsigned int _nval
        if not isinstance(self, Integer):
            return self.__pow__(int(n))
        _n = Integer(n)
        if _n < 0:
            return Integer(1)/(self**(-_n))
        _self = integer(self)
        cdef Integer x
        x = Integer()
        _nval = _n

        _sig_on
        mpz_pow_ui(x.value, _self.value, _nval)
        _sig_off

        return x

    def __pos__(self):
        return self

    def __neg__(self):
        cdef Integer x
        x = Integer()
        mpz_neg(x.value, self.value)
        return x

    def __abs__(self):
        """
        """
        cdef Integer x
        x = Integer()
        mpz_abs(x.value, self.value)
        return x

    def __mod__(self, modulus):
        cdef Integer _modulus, _self
        _modulus = integer(modulus)
        if not _modulus:
            raise ZeroDivisionError, "Integer modulo by zero"
        _self = integer(self)

        cdef Integer x
        x = Integer()

        _sig_on
        mpz_mod(x.value, _self.value, _modulus.value)
        _sig_off

        return x


    def quo_rem(self, other):
        cdef Integer _other, _self
        _other = integer(other)
        if not _other:
            raise ZeroDivisionError, "other (=%s) must be nonzero"%other
        _self = integer(self)

        cdef Integer q, r
        q = Integer()
        r = Integer()

        _sig_on
        mpz_tdiv_qr(q.value, r.value, _self.value, _other.value)
        _sig_off

        return q, r



    def powermod(self, exp, mod):
        """
        powermod(self, Integer exp, Integer mod):
        Compute self**exp modulo mod.
        """
        cdef Integer x, _exp, _mod
        _exp = Integer(exp); _mod = Integer(mod)
        x = Integer()

        _sig_on
        mpz_powm(x.value, self.value, _exp.value, _mod.value)
        _sig_off

        return x

    def powermodm_ui(self, unsigned long int exp, mod):
        cdef Integer x, _mod
        _mod = Integer(mod)
        x = Integer()

        _sig_on
        mpz_powm_ui(x.value, self.value, exp, _mod.value)
        _sig_off

        return x

    def __int__(self):
        cdef char *s
        s = mpz_get_str(NULL, 16, self.value)
        n = int(s,16)
        free(s)
        return n

    def __nonzero__(self):
        return not self.is_zero()

    def __long__(self):
        cdef char *s
        s = mpz_get_str(NULL, 16, self.value)
        n = long(s,16)
        free(s)
        return n

    def __float__(self):
        return mpz_get_d(self.value)

    def __hash__(self):
        return hash(self.str(16))

    def factor(self, algorithm='pari'):
        """
        Return the prime factorization of the integer as a list of
        pairs $(p,e)$, where $p$ is prime and $e$ is a positive integer.

        INPUT:
            algorithm -- string
                 * 'pari' -- (default)  use the PARI c library
                 * 'kash' -- use KASH computer algebra system (requires
                             the optional kash package be installed)
        """
        return sage.rings.integer_ring.factor(self, algorithm=algorithm)

    def coprime_integers(self, m):
        """
        Return the positive integers $< m$ that are coprime to self.

        EXAMPLES:
            sage: 8.coprime_integers(8)
            [1, 3, 5, 7]
            sage: 8.coprime_integers(11)
            [1, 3, 5, 7, 9]
            sage: 5.coprime_integers(10)
            [1, 2, 3, 4, 6, 7, 8, 9]
            sage: 5.coprime_integers(5)
            [1, 2, 3, 4]
            sage: 99.coprime_integers(99)
            [1, 2, 4, 5, 7, 8, 10, 13, 14, 16, 17, 19, 20, 23, 25, 26, 28, 29, 31, 32, 34, 35, 37, 38, 40, 41, 43, 46, 47, 49, 50, 52, 53, 56, 58, 59, 61, 62, 64, 65, 67, 68, 70, 71, 73, 74, 76, 79, 80, 82, 83, 85, 86, 89, 91, 92, 94, 95, 97, 98]

        AUTHORS:
            -- Naqi Jaffery (2006-01-24): examples
        """
        # TODO -- make VASTLY faster
        v = []
        for n in range(1,m):
            if self.gcd(n) == 1:
                v.append(Integer(n))
        return v

    def divides(self, n):
        """
        Return True if self divides n.

        EXAMPLES:
            sage: Z = IntegerRing()
            sage: Z(5).divides(Z(10))
            True
            sage: Z(0).divides(Z(5))
            False
            sage: Z(10).divides(Z(5))
            False
        """
        cdef int t
        cdef Integer _n
        _n = Integer(n)
        if mpz_cmp_si(self.value, 0) == 0:
            return bool(mpz_cmp_si(_n.value, 0) == 0)
        _sig_on
        t = mpz_divisible_p(_n.value, self.value)
        _sig_off
        return bool(t)


    def valuation(self, p):
        if self == 0:
            return sage.rings.infinity.infinity
        cdef int k
        k = 0
        while self % p == 0:
            k = k + 1
            self = self.__floordiv__(p)
        return Integer(k)

    def _lcm(self, Integer n):
        """
        Returns the least common multiple of self and $n$.
        """
        cdef mpz_t x

        mpz_init(x)

        _sig_on
        mpz_lcm(x, self.value, n.value)
        _sig_off


        cdef Integer z
        z = Integer()
        mpz_set(z.value,x)
        mpz_clear(x)
        return z

    def denominator(self):
        """
        Return the denominator of this integer.

        EXAMPLES:
            sage: x = 5
            sage: x.denominator()
            1
            sage: x = 0
            sage: x.denominator()
            1
        """
        return ONE

    def numerator(self):
        """
        Return the numerator of this integer.

        EXAMPLE:
            sage: x = 5
            sage: x.numerator()
            5

            sage: x = 0
            sage: x.numerator()
            0
        """
        return self

    def is_one(self):
        """
        Returns \\code{True} if the integers is $1$, otherwise \\code{False}.
        """
        return mpz_cmp_si(self.value, 1) == 0

    def is_zero(self):
        """
        Returns \\code{True} if the integers is $0$, otherwise \\code{False}.
        """
        return mpz_cmp_si(self.value, 0) == 0

    def is_unit(self):
        """
        Returns \\code{true} if this integer is a unit, i.e., 1 or $-1$.
        """
        return bool(mpz_cmp_si(self.value, -1) == 0 or mpz_cmp_si(self.value, 1) == 0)

    def is_square(self):
        return bool(self._pari_().issquare())

    def is_prime(self):
        return bool(self._pari_().isprime())

    def square_free_part(self):
        """
        Return the square free part of $x$, i.e., a divisor z such that $x = z y^2$,
        for a perfect square $y^2$.

        EXAMPLES:
            sage: square_free_part(100)
            1
            sage: square_free_part(12)
            3
            sage: square_free_part(17*37*37)
            17
            sage: square_free_part(-17*32)
            -34
            sage: square_free_part(1)
            1
            sage: square_free_part(-1)
            -1
            sage: square_free_part(-2)
            -2
            sage: square_free_part(-4)
            -1
        """
        if self.is_zero():
            return self
        F = self.factor()
        n = Integer(1)
        for p, e in F:
            if e % 2 != 0:
                n = n * p
        return n * F.unit()

    def next_prime(self):
        return Integer( (self._pari_()+1).nextprime())

    def additive_order(self):
        """
        Return the additive order of self.

        EXAMPLES:
            sage: ZZ(0).additive_order()
            1
            sage: ZZ(1).additive_order()
            Infinity
        """
        import sage.rings.infinity
        if self.is_zero():
            return Integer(1)
        else:
            return sage.rings.infinity.infinity

    def multiplicative_order(self):
        r"""
        Return the multiplicative order of self, if self is a unit, or raise
        \code{ArithmeticError} otherwise.

        EXAMPLES:
            sage: ZZ(1).multiplicative_order()
            1
            sage: ZZ(-1).multiplicative_order()
            2
            sage: ZZ(0).multiplicative_order()
            Traceback (most recent call last):
            ...
            ArithmeticError: no power of 0 is a unit
            sage: ZZ(2).multiplicative_order()
            Traceback (most recent call last):
            ...
            ArithmeticError: no power of 2 is a unit
        """
        if mpz_cmp_si(self.value, 1) == 0:
            return Integer(1)
        elif mpz_cmp_si(self.value, -1) == 0:
            return Integer(2)
        else:
            raise ArithmeticError, "no power of %s is a unit"%self

    def is_square_free(self):
        """
        Returns True if this integer is not divisible by the square of
        any prime and False otherwise.
        """
        return self._pari_().issquarefree()

    def _pari_(self):
        if self._pari is None:
            # better to do in hex, but I can't figure out
            # how to input/output a number in hex in PARI!!
            # TODO: (I could just think carefully about raw bytes and make this all much faster...)
            self._pari = sage.libs.pari.all.pari(str(self))
        return self._pari

    def parent(self):
        """
        Return the ring $\\Z$ of integers.
        """
        return sage.rings.integer_ring.Z

    def isqrt(self):
        """
        Returns the integer floor of the square root of self, or raises
        an \\exception{ValueError} if self is negative.

        EXAMPLE:
            sage: a = Integer(5)
            sage: a.isqrt()
            2
        """
        if self < 0:
            raise ValueError, "square root of negative number not defined."
        cdef Integer x
        x = Integer()

        _sig_on
        mpz_sqrt(x.value, self.value)
        _sig_off

        return x

    def _mpfr_(self, R):
        return R(str(self))  # TODO: use base 16 or something (!)


    def sqrt(self, bits=None):
        r"""
        Returns the positive square root of self, possibly as a
        \emph{a real number} if self is not a perfect integer
        square.

        INPUT:
            bits -- number of bits of precision.
                    If bits is not specified, the number of
                    bits of precision is at least twice the
                    number of bits of self (the precision
                    is always at least 53 bits if not specified).
        OUTPUT:
            integer, real number, or complex number.

        For the guaranteed integer square root of a perfect square
        (with error checking), use \code{self.square_root()}.

        EXAMPLE:
            sage: Z = IntegerRing()
            sage: Z(2).sqrt(53)
            1.4142135623730951
            sage: Z(2).sqrt(100)
            1.4142135623730950488016887242092
            sage: 39188072418583779289.square_root()
            6260037733
            sage: (100^100).sqrt()
            10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
            sage: (-1).sqrt()
            1.0000000000000000*I
            sage: sqrt(-2)
            1.4142135623730949*I
            sage: sqrt(97)
            9.8488578017961039
            sage: 97.sqrt(200)
            9.8488578017961047217462114149176244816961362874427641717231516
        """
        if bits is None:
            bits = max(53, 2*(mpz_sizeinbase(self.value, 2)+2))
        if self < 0:
            x = sage.rings.complex_field.ComplexField(bits)(self)
            return x.sqrt()
        else:
            try:
                return self.square_root()
            except ValueError:
                pass
            R = mpfr.RealField(bits)
            return self._mpfr_(R).sqrt()

    def square_root(self):
        """
        Return the positive integer square root of self, or raises a ValueError
        if self is not a perfect square.
        """
        n = self.isqrt()
        if n * n == self:
            return n
        raise ValueError, "self (=%s) is not a perfect square"%self


    def _xgcd(self, Integer n):
        """
        Return a triple $g, s, t \\in\\Z$ such that
        $$
           g = s \\cdot \\mbox{\\rm self} + t \\cdot n.
        $$
        """
        cdef mpz_t g, s, t
        cdef object g0, s0, t0

        mpz_init(g)
        mpz_init(s)
        mpz_init(t)

        _sig_on
        mpz_gcdext(g, s, t, self.value, n.value)
        _sig_off

        g0 = Integer()
        s0 = Integer()
        t0 = Integer()
        set_mpz(g0,g)
        set_mpz(s0,s)
        set_mpz(t0,t)
        mpz_clear(g)
        mpz_clear(s)
        mpz_clear(t)
        return g0, s0, t0

    def __lshift__(Integer self, n):
        """
        """
        n = int(n)
        return self*(Integer(2)**n)   # todo: optimize

    def __rshift__(Integer self, n):
        """
        """
        n = int(n)
        return self.__floordiv(Integer(2)**n)    # todo: optimize

    def _and(Integer self, Integer other):
        cdef Integer x
        x = Integer()
        mpz_and(x.value, self.value, other.value)
        return x

    def __and__(x, y):
        if isinstance(x, Integer) and isinstance(y, Integer):
            return x._and(y)
        return sage.rings.coerce.bin_op(x, y, operator.and_)


    def _or(Integer self, Integer other):
        cdef Integer x
        x = Integer()
        mpz_ior(x.value, self.value, other.value)
        return x

    def __or__(x, y):
        if isinstance(x, Integer) and isinstance(y, Integer):
            return x._or(y)
        return sage.rings.coerce.bin_op(x, y, operator.or_)


    def __invert__(self):
        """
        """
        return Integer(1)/self    # todo: optimize


    def inverse_mod(self, n):
        """
        Returns the inverse of self modulo $n$, if this inverse exists.
        Otherwise, raises a \exception{ZeroDivisionError} exception.

        INPUT:
           self -- Integer
           n -- Integer
        OUTPUT:
           x -- Integer such that x*self = 1 (mod m), or
                raises ZeroDivisionError.
        IMPLEMENTATION:
           Call the mpz_invert GMP library function.
        EXAMPLES:
            sage: a = Integer(189)
            sage: a.inverse_mod(10000)
            4709
            sage: a.inverse_mod(-10000)
            4709
            sage: a.inverse_mod(1890)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Inverse does not exist.
            sage: a = Integer(19)**100000
            sage: b = a*a
            sage: c=a.inverse_mod(b)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Inverse does not exist.
        """
        cdef mpz_t x
        cdef object ans
        cdef int r
        cdef Integer m
        m = Integer(n)

        if m == 1:
            return Integer(0)

        mpz_init(x)

        _sig_on
        r = mpz_invert(x, self.value, m.value)
        _sig_off

        if r == 0:
            raise ZeroDivisionError, "Inverse does not exist."
        ans = Integer()
        set_mpz(ans,x)
        mpz_clear(x)
        return ans

    def _gcd(self, Integer n):
        """
        Return the greatest common divisor of self and $n$.
        """
        cdef mpz_t g
        cdef object g0

        mpz_init(g)


        _sig_on
        mpz_gcd(g, self.value, n.value)
        _sig_off

        g0 = Integer()
        set_mpz(g0,g)
        mpz_clear(g)
        return g0

    def crt(self, y, m, n):
        """
        Return the unique integer between $0$ and $mn$ that is
        congruent to the integer modulo $m$ and to $y$ modulo $n$.  We
        assume that~$m$ and~$n$ are coprime.
        """
        cdef object g, s, t
        cdef Integer _y, _m, _n
        _y = Integer(y); _m = Integer(m); _n = Integer(n)
        g, s, t = _m.xgcd(_n)
        if not g.is_one():
            raise ArithmeticError, "CRT requires that gcd of moduli is 1."
        # Now s*m + t*n = 1, so the answer is x + (y-x)*s*m, where x=self.
        return (self + (_y-self)*s*_m) % (_m*_n)



ONE = Integer(1)

def integer(x):
    if isinstance(x, Integer):
        return x
    return Integer(x)

def factorial(unsigned long int n):
    """
    Return the factorial $n!=1 \\cdot 2 \\cdot 3 \\cdots n$.
    The input integer $n$ must fit in an \\code{unsigned long int}.
    """
    cdef mpz_t x
    cdef Integer z

    mpz_init(x)

    _sig_on
    mpz_fac_ui(x, n)
    _sig_off

    z = Integer()
    set_mpz(z, x)
    mpz_clear(x)
    return z

