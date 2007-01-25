r"""
Elements of the ring $\Z$ of integers

AUTHORS:
    -- William Stein (2005): initial version
    -- Gonzalo Tornaria (2006-03-02): vastly improved python/GMP conversion; hashing
    -- Didier Deshommes <dfdeshom@gmail.com> (2006-03-06): numerous examples and docstrings
    -- William Stein (2006-03-31): changes to reflect GMP bug fixes
    -- William Stein (2006-04-14): added GMP factorial method (since it's
                                   now very fast).
    -- David Harvey (2006-09-15): added nth_root, exact_log
    -- David Harvey (2006-09-16): attempt to optimise Integer constructor

EXAMPLES:
   Add 2 integers:
       sage: a = Integer(3) ; b = Integer(4)
       sage: a + b == 7
       True

   Add an integer and a real number:
       sage: a + 4.0
       7.00000000000000

   Add an integer and a rational number:
       sage: a + Rational(2)/5
       17/5

   Add an integer and a complex number:
       sage: b = ComplexField().0 + 1.5
       sage: loads((a+b).dumps()) == a+b
       True

   sage: z = 32
   sage: -z
   -32
   sage: z = 0; -z
   0
   sage: z = -0; -z
   0
   sage: z = -1; -z
   1

Multiplication:
    sage: a = Integer(3) ; b = Integer(4)
    sage: a * b == 12
    True
    sage: loads((a * 4.0).dumps()) == a*b
    True
    sage: a * Rational(2)/5
    6/5

    sage: list([2,3]) * 4
    [2, 3, 2, 3, 2, 3, 2, 3]

    sage: 'sage'*Integer(3)
    'sagesagesage'

Coercions:
    Returns version of this integer in the multi-precision floating
    real field R.

        sage: n = 9390823
        sage: RR = RealField(200)
        sage: RR(n)
        9390823.0000000000000000000000000000000000000000000000000000

"""

#*****************************************************************************
#       Copyright (C) 2004 William Stein <wstein@gmail.com>
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

import sys

include "../ext/gmp.pxi"
include "../ext/interrupt.pxi"  # ctrl-c interrupt block support
include "../ext/stdsage.pxi"

cdef extern from "../ext/mpz_pylong.h":
    cdef mpz_get_pylong(mpz_t src)
    cdef int mpz_set_pylong(mpz_t dst, src) except -1
    cdef long mpz_pythonhash(mpz_t src)


from sage.libs.pari.gen cimport gen as pari_gen

cdef class Integer(sage.structure.element.EuclideanDomainElement)

import sage.rings.integer_ring
import sage.rings.infinity
import sage.rings.complex_field
import rational
import sage.libs.pari.all
import real_mpfr


cdef mpz_t mpz_tmp
mpz_init(mpz_tmp)

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
    return sage_realloc(ptr, new_size)

cdef void pymem_free(void *ptr, size_t size):
    sage_free(ptr)

cdef void* pymem_malloc(size_t size):
    return sage_malloc(size)

cdef extern from "gmp.h":
    void mp_set_memory_functions (void *(*alloc_func_ptr) (size_t),  \
                                  void *(*realloc_func_ptr) (void *, size_t, size_t),    \
                                  void (*free_func_ptr) (void *, size_t))


def pmem_malloc():
    """
    Use our own memory manager for for GMP memory management.
    """
    mp_set_memory_functions(sage_malloc, pymem_realloc, pymem_free)
    #mp_set_memory_functions(PyMem_Malloc, pymem_realloc, pymem_free)
    #mp_set_memory_functions(pymem_malloc, pymem_realloc, pymem_free)

pmem_malloc()

cdef object the_integer_ring
the_integer_ring = sage.rings.integer_ring.Z

from sage.structure.sage_object cimport SageObject
from sage.structure.element cimport EuclideanDomainElement, ModuleElement
from sage.structure.element import  bin_op

cdef class Integer(sage.structure.element.EuclideanDomainElement):
    r"""
    The \class{Integer} class represents arbitrary precision
    integers.  It derives from the \class{Element} class, so
    integers can be used as ring elements anywhere in SAGE.

    \begin{notice}
    The class \class{Integer} is implemented in Pyrex, as a wrapper
    of the GMP \code{mpz_t} integer type.
    \end{notice}
    """

    # todo: It would be really nice if we could avoid the __new__ call.
    # It has python calling conventions, and our timing tests indicate the
    # overhead can be significant. The difficulty is that then we can't
    # guarantee that the initialised will be performed exactly once.

    def __new__(self, x=None, unsigned int base=0):
        mpz_init(self.value)
        self._parent = <SageObject>the_integer_ring


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
            sage: 2^3
            8
        """
    def __init__(self, x=None, unsigned int base=0):
        """
        EXAMPLES:
            sage: a = long(-901824309821093821093812093810928309183091832091)
            sage: b = ZZ(a); b
            -901824309821093821093812093810928309183091832091
            sage: ZZ(b)
            -901824309821093821093812093810928309183091832091
            sage: ZZ('-901824309821093821093812093810928309183091832091')
            -901824309821093821093812093810928309183091832091
            sage: ZZ(int(-93820984323))
            -93820984323
            sage: ZZ(ZZ(-901824309821093821093812093810928309183091832091))
            -901824309821093821093812093810928309183091832091
            sage: ZZ(QQ(-901824309821093821093812093810928309183091832091))
            -901824309821093821093812093810928309183091832091
            sage: ZZ(pari('Mod(-3,7)'))
            4
            sage: ZZ('sage')
            Traceback (most recent call last):
            ...
            TypeError: unable to convert x (=sage) to an integer
            sage: Integer('zz',36).str(36)
            'zz'
            sage: ZZ('0x3b').str(16)
            '3b'
        """

        # TODO: All the code below should somehow be in an external
        # cdef'd function.  Then e.g., if a matrix or vector or
        # polynomial is getting filled by mpz_t's, it can use the
        # rules below to do the fill construction of mpz_t's, but
        # without the overhead of creating any Python objects at all.
        # The cdef's function should be of the form
        #     mpz_init_set_sage(mpz_t y, object x)
        # Then this function becomes the one liner:
        #     mpz_init_set_sage(self.value, x)


        if not (x is None):
            # First do all the type-check versions; these are fast.

            if PY_TYPE_CHECK(x, Integer):
                set_from_Integer(self, <Integer>x)

            elif PyInt_Check(x):
                mpz_set_si(self.value, x)

            elif PyLong_Check(x):
                mpz_set_pylong(self.value, x)

            elif PyString_Check(x):
                if base < 0 or base > 36:
                    raise ValueError, "base (=%s) must be between 2 and 36"%base
                if mpz_set_str(self.value, x, base) != 0:
                    raise TypeError, "unable to convert x (=%s) to an integer"%x

            # todo: I want to skip the name lookup here (rational.Rational).
            # I tried importing/cimporting Rational in various ways, but every
            # way was broken for some mysterious reason. Perhaps a circular
            # include somewhere?  -- David Harvey
            # Sagex does not allow circular imports of cdef'd types in any situation.
            # Since Rational cimports integer, it is not possible for integer to cimport rational.
            # This might not be fixable.  -- William Stein
            elif PY_TYPE_CHECK(x, rational.Rational):
                if x.denominator() != 1:
                    raise TypeError, "Unable to coerce rational (=%s) to an Integer."%x
                set_from_Integer(self, x.numer())

            # Similarly for "sage.libs.pari.all.pari_gen"
            elif PY_TYPE_CHECK(x, pari_gen):
                if x.type() == 't_INTMOD':
                    x = x.lift()
                # TODO: figure out how to convert to pari integer in base 16 ?

                # todo: having this "s" variable around here is causing
                # pyrex to play games with refcount for the None object, which
                # seems really stupid.

                s = hex(x)
                if mpz_set_str(self.value, s, 16) != 0:
                    raise TypeError, "Unable to coerce PARI %s to an Integer."%x
            elif PyObject_HasAttrString(x, "_integer_"):
                # todo: Note that PyObject_GetAttrString returns NULL if
                # the attribute was not found. If we could test for this,
                # we could skip the double lookup. Unfortunately pyrex doesn't
                # seem to let us do this; it flags an error if the function
                # returns NULL, because it can't construct an "object" object
                # out of the NULL pointer. This really sucks. Perhaps we could
                # make the function prototype have return type void*, but
                # then how do we make Pyrex handle the reference counting?
                set_from_Integer(self, (<object> PyObject_GetAttrString(x, "_integer_"))())

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

    def __index__(self):
        """
        Needed so integers can be used as list indices.

        EXAMPLES:
            sage: v = [1,2,3,4,5]
            sage: v[Integer(3)]
            4
            sage: v[Integer(2):Integer(4)]
            [3, 4]
        """
        return int(mpz_get_pylong(self.value))

    def _im_gens_(self, codomain, im_gens):
        return codomain._coerce_(self)

    def _xor(Integer self, Integer other):
        cdef Integer x
        x = Integer()
        mpz_xor(x.value, self.value, other.value)
        return x

    def __xor__(x, y):
        """
        Compute the exclusive or of x and y.

        EXAMPLES:
            sage: n = ZZ(2); m = ZZ(3)
            sage: n.__xor__(m)
            1
        """
        if PY_TYPE_CHECK(x, Integer) and PY_TYPE_CHECK(y, Integer):
            return x._xor(y)
        return bin_op(x, y, operator.xor)


    def __richcmp__(left, right, int op):
        return (<sage.structure.element.Element>left)._richcmp(right, op)

    cdef int _cmp_c_impl(left, sage.structure.element.Element right) except -2:
        cdef int i
        i = mpz_cmp((<Integer>left).value, (<Integer>right).value)
        if i < 0: return -1
        elif i == 0: return 0
        else: return 1

    def __copy__(self):
        """
        Return a copy of the integer.

        EXAMPLES:
            sage: n = 2
            sage: copy(n)
            2
            sage: copy(n) is n
            False
        """
        cdef Integer z
        z = Integer()
        set_mpz(z,self.value)
        return z


    def list(self):
        """
        Return a list with this integer in it, to be
        compatible with the method for number fields.

        EXAMPLES:
            sage: m = 5
            sage: m.list()
            [5]
        """
        return [ self ]


    def  __dealloc__(self):
        mpz_clear(self.value)

    def __repr__(self):
        return self.str()

    def _latex_(self):
        return self.str()

    def _mathml_(self):
        return '<mn>%s</mn>'%self

    def __str_malloc(self, int base=10):
        r"""
        Return the string representation of \code{self} in the given
        base.

        However, self.str() below is nice because we know the size of
        the string ahead of time, and can work around a bug in GMP
        nicely.  There seems to be a bug in GMP, where non-2-power
        base conversion for very large integers > 10 million digits
        (?) crashes GMP.
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


        EXAMPLES:
            sage: Integer(2^10).str(2)
            '10000000000'
            sage: Integer(2^10).str(17)
            '394'

            sage: two=Integer(2)
            sage: two.str(1)
            Traceback (most recent call last):
            ...
            ValueError: base (=1) must be between 2 and 36

            sage: two.str(37)
            Traceback (most recent call last):
            ...
            ValueError: base (=37) must be between 2 and 36

            sage: big = 10^5000000
            sage: s = big.str()                 # long time (> 20 seconds)
            sage: len(s)                        # long time (depends on above defn of s)
            5000001
            sage: s[:10]                        # long time (depends on above defn of s)
            '1000000000'
        """
        if base < 2 or base > 36:
            raise ValueError, "base (=%s) must be between 2 and 36"%base
        cdef size_t n
        cdef char *s
        n = mpz_sizeinbase(self.value, base) + 2
        s = <char *>PyMem_Malloc(n)
        if s == NULL:
            raise MemoryError, "Unable to allocate enough memory for the string representation of an integer."
        _sig_on
        mpz_get_str(s, base, self.value)
        _sig_off
        k = <object> PyString_FromString(s)
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
            0x36bb1e3929d1a8fe2802f083L
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

        EXAMPLES:
            sage: n= ZZ(54)
            sage: n.set_si(-43344);n
            -43344
            sage: n.set_si(43344);n
            43344

        Note that an error occurs when we are not dealing with
        integers anymore
            sage: n.set_si(2^32);n
            Traceback (most recent call last):      # 32-bit
            ...                                     # 32-bit
            OverflowError: long int too large to convert to int   # 32-bit
            4294967296       # 64-bit
            sage: n.set_si(-2^32);n
            Traceback (most recent call last):      # 32-bit
            ...                                     # 32-bit
            OverflowError: long int too large to convert to int     # 32-bit
            -4294967296      # 64-bit
        """
        mpz_set_si(self.value, n)

    def set_str(self, s, base=10):
        """
        Set self equal to the number defined by the string $s$ in the
        given base.

        EXAMPLES:
            sage: n=100
            sage: n.set_str('100000',2)
            sage: n
            32

        If the number begins with '0X' or '0x', it is converted
        to an hex number:
            sage: n.set_str('0x13',0)
            sage: n
            19
            sage: n.set_str('0X13',0)
            sage: n
            19

        If the number begins with a '0', it is converted to an octal
        number:
            sage: n.set_str('013',0)
            sage: n
            11

        '13' is not a valid binary number so the following raises
        an exception:
            sage: n.set_str('13',2)
            Traceback (most recent call last):
            ...
            TypeError: unable to convert x (=13) to an integer in base 2
        """
        valid = mpz_set_str(self.value, s, base)
        if valid != 0:
            raise TypeError, "unable to convert x (=%s) to an integer in base %s"%(s, base)

    cdef void set_from_mpz(Integer self, mpz_t value):
        mpz_set(self.value, value)

    cdef mpz_t* get_value(Integer self):
        return &self.value

    cdef ModuleElement _add_c_impl(self, ModuleElement right):
        # self and right are guaranteed to be Integers
        cdef Integer x
        x = <Integer> PY_NEW(Integer)
        mpz_add(x.value, self.value, (<Integer>right).value)
        return x

##     def _unsafe_add_in_place(self,  ModuleElement right):
##         """
##         Do *not* use this...  unless you really know what you
##         are doing.
##         """
##         if not (right._parent is self._parent):
##             raise TypeError
##         mpz_add(self.value, self.value, (<Integer>right).value)
##     cdef _unsafe_add_in_place_c(self,  ModuleElement right):
##         """
##         Do *not* use this...  unless you really know what you
##         are doing.
##         """
##         if not (right._parent is self._parent):
##             raise TypeError
##         mpz_add(self.value, self.value, (<Integer>right).value)

    cdef ModuleElement _sub_c_impl(self, ModuleElement right):
        # self and right are guaranteed to be Integers
        cdef Integer x
        x = <Integer> PY_NEW(Integer)
        mpz_sub(x.value, self.value, (<Integer>right).value)
        return x

    cdef ModuleElement _neg_c_impl(self):
        cdef Integer x
        x = Integer()
        mpz_neg(x.value, self.value)
        return x

    def _r_action(self, s):
        if isinstance(s, (str, list)):
            return s*int(self)
        raise TypeError

    cdef RingElement _mul_c_impl(self, RingElement right):
        # self and right are guaranteed to be Integers
        cdef Integer x
        x = <Integer> PY_NEW(Integer)
        if  mpz_sizeinbase(self.value, 2) > 1000000:  # some lack of symmetry
            # We only use the signal handler (to enable ctrl-c out) in case
            # self is huge, so the product might actually take a while to compute.
            _sig_on
            mpz_mul(x.value, self.value, (<Integer>right).value)
            _sig_off
        else:
            mpz_mul(x.value, self.value, (<Integer>right).value)
        return x

    cdef RingElement _div_c_impl(self, RingElement right):
        r"""
        Computes a \over{b}

        EXAMPLES:
            sage: a = Integer(3) ; b = Integer(4)
            sage: a / b == Rational(3) / 4
            True
            sage: Integer(32) / Integer(32)
            1
        """
        # todo -- optimize -- this is really slow.
        return rational.Rational(self)/rational.Rational(right)


    def __floordiv(Integer self, Integer other):
        cdef Integer x
        x = Integer()


        _sig_on
        mpz_fdiv_q(x.value, self.value, other.value)
        _sig_off

        return x


    def __floordiv__(x, y):
        r"""
        Computes the whole part of self \over{other}

        EXAMPLES:
            sage: a = Integer(321) ; b = Integer(10)
            sage: a // b
            32
        """
        if PY_TYPE_CHECK(x, Integer) and PY_TYPE_CHECK(y, Integer):
            return x.__floordiv(y)
        return bin_op(x, y, operator.floordiv)


    def __pow__(self, n, dummy):
        r"""
        Computes $\text{self}^n$

        EXAMPLES:
            sage: 2^-6
            1/64
            sage: 2^6
            64
            sage: 2^0
            1
            sage: 2^-0
            1
            sage: (-1)^(1/3)
            Traceback (most recent call last):
            ...
            TypeError: exponent (=1/3) must be an integer.
            Coerce your numbers to real or complex numbers first.
        """
        cdef Integer _self, _n
        cdef unsigned int _nval
        if not PY_TYPE_CHECK(self, Integer):
            return self.__pow__(int(n))
        try:
            _n = Integer(n)
        except TypeError:
            raise TypeError, "exponent (=%s) must be an integer.\nCoerce your numbers to real or complex numbers first."%n
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

    def nth_root(self, int n, int report_exact=0):
        r"""
        Returns the truncated nth root of self.

        INPUT:
            n -- integer >= 1 (must fit in C int type)
            report_exact -- boolean, whether to report if the root extraction
                          was exact

        OUTPUT:
           If report_exact is 0 (default), then returns the truncation of the
           nth root of self (i.e. rounded towards zero).

           If report_exact is 1, then returns the nth root and a boolean
           indicating whether the root extraction was exact.

        AUTHOR:
           -- David Harvey (2006-09-15)

        EXAMPLES:
          sage: Integer(125).nth_root(3)
          5
          sage: Integer(124).nth_root(3)
          4
          sage: Integer(126).nth_root(3)
          5

          sage: Integer(-125).nth_root(3)
          -5
          sage: Integer(-124).nth_root(3)
          -4
          sage: Integer(-126).nth_root(3)
          -5

          sage: Integer(125).nth_root(2, True)
          (11, False)
          sage: Integer(125).nth_root(3, True)
          (5, True)

          sage: Integer(125).nth_root(-5)
          Traceback (most recent call last):
          ...
          ValueError: n (=-5) must be positive

          sage: Integer(-25).nth_root(2)
          Traceback (most recent call last):
          ...
          ValueError: cannot take even root of negative number

        """
        if n < 1:
            raise ValueError, "n (=%s) must be positive" % n
        if (self < 0) and not (n & 1):
            raise ValueError, "cannot take even root of negative number"
        cdef Integer x
        cdef int is_exact
        x = Integer()
        _sig_on
        is_exact = mpz_root(x.value, self.value, n)
        _sig_off

        if report_exact:
            return x, bool(is_exact)
        else:
            return x

    def exact_log(self, m):
        r"""
        Returns the largest integer $k$ such that $m^k \leq \text{self}$, i.e.,
        the floor of $\log_m(\text{self})$.

        This is guaranteed to return the correct answer even when the usual
        log function doesn't have sufficient precision.

        INPUT:
            m -- integer >= 2

        AUTHOR:
           -- David Harvey (2006-09-15)

        TODO:
           -- Currently this is extremely stupid code (although it should
           always work). Someone needs to think about doing this properly
           by estimating errors in the log function etc.

        EXAMPLES:
           sage: Integer(125).exact_log(5)
           3
           sage: Integer(124).exact_log(5)
           2
           sage: Integer(126).exact_log(5)
           3
           sage: Integer(3).exact_log(5)
           0
           sage: Integer(1).exact_log(5)
           0


           sage: x = 3^100000
           sage: log(RR(x), 3)
           99999.9999999999
           sage: log(RR(x + 100000), 3)
           99999.9999999999

           sage: x.exact_log(3)
           100000
           sage: (x+1).exact_log(3)
           100000
           sage: (x-1).exact_log(3)
           99999
        """
        if self <= 0:
            raise ValueError, "self must be positive"
        if m < 2:
            raise ValueError, "m must be at least 2"
        R = real_mpfr.RealField(53)
        guess = R(self).log(base = m).floor()
        power = m ** guess

        while power > self:
            power = power / m
            guess = guess - 1

        if power == self:
            return guess

        while power < self:
            power = power * m
            guess = guess + 1

        if power == self:
            return guess
        else:
            return guess - 1


    def __pos__(self):
        """
        EXAMPLES:
            sage: z=43434
            sage: z.__pos__()
            43434
        """
        return self

    def __abs__(self):
        """
        Computes $|self|$

        EXAMPLES:
            sage: z = -1
            sage: abs(z)
            1
            sage: abs(z) == abs(1)
            True
        """
        cdef Integer x
        x = Integer()
        mpz_abs(x.value, self.value)
        return x

    def __mod__(self, modulus):
        r"""
        Returns \code{self % modulus}.

        EXAMPLES:
            sage: z = 43
            sage: z % 2
            1
            sage: z % 0
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Integer modulo by zero
        """
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
        """
        Returns the quotient and the remainder of
        self divided by other.

        INPUT:
            other -- the integer the divisor

        OUTPUT:
            q   -- the quotient of self/other
            r   -- the remainder of self/other

        EXAMPLES:
            sage: z = Integer(231)
            sage: z.quo_rem(2)
            (115, 1)
            sage: z.quo_rem(-2)
            (-115, 1)
            sage: z.quo_rem(0)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: other (=0) must be nonzero
        """
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

    def div(self, other):
        """
        Returns the quotient of self divided by other.

        INPUT:
            other -- the integer the divisor

        OUTPUT:
            q   -- the quotient of self/other

        EXAMPLES:
            sage: z = Integer(231)
            sage: z.div(2)
            115
            sage: z.div(-2)
            -115
            sage: z.div(0)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: other (=0) must be nonzero
        """
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

        return q


    def powermod(self, exp, mod):
        """
        Compute self**exp modulo mod.

        EXAMPLES:
            sage: z = 2
            sage: z.powermod(31,31)
            2
            sage: z.powermod(0,31)
            1
            sage: z.powermod(-31,31) == 2^-31 % 31
            True

            As expected, the following is invalid:
            sage: z.powermod(31,0)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: cannot raise to a power modulo 0
        """
        cdef Integer x, _exp, _mod
        _exp = Integer(exp); _mod = Integer(mod)
        if mpz_cmp_si(_mod.value,0) == 0:
            raise ZeroDivisionError, "cannot raise to a power modulo 0"

        x = Integer()

        _sig_on
        mpz_powm(x.value, self.value, _exp.value, _mod.value)
        _sig_off

        return x

    def rational_reconstruction(self, Integer m):
        return rational.pyrex_rational_reconstruction(self, m)

    def powermodm_ui(self, exp, mod):
        r"""
        Computes self**exp modulo mod, where exp is an unsigned
        long integer.

        EXAMPLES:
            sage: z = 32
            sage: z.powermodm_ui(2, 4)
            0
            sage: z.powermodm_ui(2, 14)
            2
            sage: z.powermodm_ui(2^31-1, 14)
            4
            sage: z.powermodm_ui(2^31, 14)
            Traceback (most recent call last):                              # 32-bit
            ...                                                             # 32-bit
            OverflowError: exp (=2147483648) must be <= 2147483647   # 32-bit
            2              # 64-bit
            sage: z.powermodm_ui(2^63, 14)
            Traceback (most recent call last):
            ...
            OverflowError: exp (=9223372036854775808) must be <= 2147483647           # 32-bit
            OverflowError: exp (=9223372036854775808) must be <= 9223372036854775807  # 64-bit
        """
        if exp < 0:
            raise ValueError, "exp (=%s) must be nonnegative"%exp
        elif exp > sys.maxint:
            raise OverflowError, "exp (=%s) must be <= %s"%(exp, sys.maxint)
        cdef Integer x, _mod
        _mod = Integer(mod)
        x = Integer()

        _sig_on
        mpz_powm_ui(x.value, self.value, exp, _mod.value)
        _sig_off

        return x

    def __int__(self):
        return int(mpz_get_pylong(self.value))

    def __long__(self):
        return mpz_get_pylong(self.value)

    def __nonzero__(self):
        return not self.is_zero()

    def __float__(self):
        return mpz_get_d(self.value)

    def __hash__(self):
        return mpz_pythonhash(self.value)

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
            sage: n = 8
            sage: n.coprime_integers(8)
            [1, 3, 5, 7]
            sage: n.coprime_integers(11)
            [1, 3, 5, 7, 9]
            sage: n = 5; n.coprime_integers(10)
            [1, 2, 3, 4, 6, 7, 8, 9]
            sage: n.coprime_integers(5)
            [1, 2, 3, 4]
            sage: n = 99; n.coprime_integers(99)
            [1, 2, 4, 5, 7, 8, 10, 13, 14, 16, 17, 19, 20, 23, 25, 26, 28, 29, 31, 32, 34, 35, 37, 38, 40, 41, 43, 46, 47, 49, 50, 52, 53, 56, 58, 59, 61, 62, 64, 65, 67, 68, 70, 71, 73, 74, 76, 79, 80, 82, 83, 85, 86, 89, 91, 92, 94, 95, 97, 98]

        AUTHORS:
            -- Naqi Jaffery (2006-01-24): examples

        ALGORITHM: Naive -- compute lots of GCD's.  If this isn't good
        enough for you, please code something better and submit a
        patch.
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
        """
        Return the p-adic valuation of self.

        EXAMPLE:
            sage: n = 60
            sage: n.valuation(2)
            2
            sage: n.valuation(3)
            1
            sage: n.valuation(7)
            0
        """
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

        EXAMPLES:
            sage: n = 60
            sage: n._lcm(150)
            300
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

    def factorial(self):
        """
        Return the factorial $n!=1 \\cdot 2 \\cdot 3 \\cdots n$.
        Self must fit in an \\code{unsigned long int}.

        EXAMPLES:
            sage: for n in srange(7):
            ...    print n, n.factorial()
            0 1
            1 1
            2 2
            3 6
            4 24
            5 120
            6 720
        """
        if self < 0:
            raise ValueError, "factorial -- self = (%s) must be nonnegative"%self

        if mpz_cmp_ui(self.value,4294967295) > 0:
            raise ValueError, "factorial not implemented for n >= 2^32.\nThis is probably OK, since the answer would have billions of digits."

        cdef unsigned int n
        n = self

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

    def floor(self):
        """
        Return the floor of self, which is just self since self is an integer.

        EXAMPLES:
            sage: n = 6
            sage: n.floor()
            6
        """
        return self

    def ceil(self):
        """
        Return the ceiling of self, which is self since self is an integer.

        EXAMPLES:
            sage: n = 6
            sage: n.ceil()
            6
        """
        return self

    def is_one(self):
        r"""
        Returns \code{True} if the integers is $1$, otherwise \code{False}.

        EXAMPLES:
            sage: Integer(1).is_one()
            True
            sage: Integer(0).is_one()
            False
        """
        return bool(mpz_cmp_si(self.value, 1) == 0)

    def is_zero(self):
        r"""
        Returns \code{True} if the integers is $0$, otherwise \code{False}.

        EXAMPLES:
            sage: Integer(1).is_zero()
            False
            sage: Integer(0).is_zero()
            True
        """
        return bool(mpz_cmp_si(self.value, 0) == 0)

    def is_unit(self):
        r"""
        Returns \code{true} if this integer is a unit, i.e., 1 or $-1$.

        EXAMPLES:
            sage: for n in srange(-2,3):
            ...    print n, n.is_unit()
            -2 False
            -1 True
            0 False
            1 True
            2 False
        """
        return bool(mpz_cmp_si(self.value, -1) == 0 or mpz_cmp_si(self.value, 1) == 0)

    def is_square(self):
        r"""
        Returns \code{True} if self is a perfect square

        EXAMPLES:
            sage: Integer(4).is_square()
            True
            sage: Integer(41).is_square()
            False
        """
        return bool(self._pari_().issquare())

    def is_prime(self):
        r"""
        Retuns \code{True} if self is prime

        EXAMPLES:
            sage: z = 2^31 - 1
            sage: z.is_prime()
            True
            sage: z = 2^31
            sage: z.is_prime()
            False
        """
        return bool(self._pari_().isprime())

    def is_pseudoprime(self):
        r"""
        Retuns \code{True} if self is a pseudoprime

        EXAMPLES:
            sage: z = 2^31 - 1
            sage: z.is_pseudoprime()
            True
            sage: z = 2^31
            sage: z.is_pseudoprime()
            False
        """
        return bool(self._pari_().ispseudoprime())

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
        r"""
        Returns the next prime after self

        EXAMPLES:
            sage: Integer(100).next_prime()
            101
            sage: Integer(0).next_prime()
            2
            sage: Integer(1001).next_prime()
            1009
        """
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

    def is_squarefree(self):
        """
        Returns True if this integer is not divisible by the square of
        any prime and False otherwise.

        EXAMPLES:
            sage: Integer(100).is_squarefree()
            False
            sage: Integer(102).is_squarefree()
            True
        """
        return self._pari_().issquarefree()

    def _pari_(self):
        """
        Returns the PARI version of this integer.

        EXAMPLES:
            sage: n = 9390823
            sage: m = n._pari_(); m
            9390823
            sage: type(m)
            <type 'sage.libs.pari.gen.gen'>

        ALGORITHM: Use base 10 Python string conversion, hence very
        very slow for large integers. If you can figure out how to
        input a number into PARI in hex, or otherwise optimize this,
        please implement it and send me a patch.
        """
        if self._pari is None:
            # better to do in hex, but I can't figure out
            # how to input/output a number in hex in PARI!!
            # TODO: (I could just think carefully about raw bytes and make this all much faster...)
            self._pari = sage.libs.pari.all.pari(str(self))
        return self._pari

    def _interface_init_(self):
        """
        Return canonical string to coerce this integer to any other math
        software, i.e., just the string representation of this integer
        in base 10.

        EXAMPLES:
            sage: n = 9390823
            sage: n._interface_init_()
            '9390823'
        """
        return str(self)

    def isqrt(self):
        r"""
        Returns the integer floor of the square root of self, or raises
        an \exception{ValueError} if self is negative.

        EXAMPLE:
            sage: a = Integer(5)
            sage: a.isqrt()
            2

            sage: Integer(-102).isqrt()
            Traceback (most recent call last):
            ...
            ValueError: square root of negative number not defined.
        """
        if self < 0:
            raise ValueError, "square root of negative number not defined."
        cdef Integer x
        x = Integer()

        _sig_on
        mpz_sqrt(x.value, self.value)
        _sig_off

        return x


    def sqrt(self, bits=None):
        r"""
        Returns the positive square root of self, possibly as a
        \emph{a real or complex number} if self is not a perfect
        integer square.

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
            sage: Z(4).sqrt()
            2
            sage: Z(4).sqrt(53)
            2.00000000000000
            sage: Z(2).sqrt(53)
            1.41421356237309
            sage: Z(2).sqrt(100)
            1.4142135623730950488016887242
            sage: n = 39188072418583779289; n.square_root()
            6260037733
            sage: (100^100).sqrt()
            10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
            sage: (-1).sqrt()
            1.00000000000000*I
            sage: sqrt(-2)
            1.41421356237309*I
            sage: sqrt(97)
            9.84885780179610
            sage: n = 97; n.sqrt(200)
            9.8488578017961047217462114149176244816961362874427641717231
        """
        if bits is None:
            try:
                return self.square_root()
            except ValueError:
                pass
            bits = max(53, 2*(mpz_sizeinbase(self.value, 2)+2))

        if self < 0:
            x = sage.rings.complex_field.ComplexField(bits)(self)
            return x.sqrt()
        else:
            R = real_mpfr.RealField(bits)
            return R(self).sqrt()

    def square_root(self):
        """
        Return the positive integer square root of self, or raises a ValueError
        if self is not a perfect square.

        EXAMPLES:
            sage: Integer(144).square_root()
            12
            sage: Integer(102).square_root()
            Traceback (most recent call last):
            ...
            ValueError: self (=102) is not a perfect square
        """
        n = self.isqrt()
        if n * n == self:
            return n
        raise ValueError, "self (=%s) is not a perfect square"%self


    def _xgcd(self, Integer n):
        r"""
        Return a triple $g, s, t \in\Z$ such that
        $$
           g = s \cdot \mbox{\rm self} + t \cdot n.
        $$

        EXAMPLES:
            sage: n = 6
            sage: g, s, t = n._xgcd(8)
            sage: s*6 + 8*t
            2
            sage: g
            2
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

    cdef _lshift(self, unsigned long int n):
        """
        Shift self n bits to the left, i.e., quickly multiply by $2^n$.
        """
        cdef Integer x
        x = Integer()

        _sig_on
        mpz_mul_2exp(x.value, self.value, n)
        _sig_off
        return x

    def __lshift__(x,y):
        """
        Shift x y bits to the left.

        EXAMPLES:
            sage: 32 << 2
            128
            sage: 32 << int(2)
            128
            sage: int(32) << 2
            128
        """
        if PY_TYPE_CHECK(x, Integer) and isinstance(y, (Integer, int, long)):
            return (<Integer>x)._lshift(long(y))
        return bin_op(x, y, operator.lshift)

    cdef _rshift(Integer self, unsigned long int n):
        cdef Integer x
        x = Integer()
        _sig_on
        mpz_fdiv_q_2exp(x.value, self.value, n)
        _sig_off
        return x

    def __rshift__(x, y):
        """
        EXAMPLES:
            sage: 32 >> 2
            8
            sage: 32 >> int(2)
            8
            sage: int(32) >> 2
            8
        """
        if PY_TYPE_CHECK(x, Integer) and isinstance(y, (Integer, int, long)):
            return (<Integer>x)._rshift(long(y))
        return bin_op(x, y, operator.rshift)

    cdef _and(Integer self, Integer other):
        cdef Integer x
        x = Integer()
        mpz_and(x.value, self.value, other.value)
        return x

    def __and__(x, y):
        if PY_TYPE_CHECK(x, Integer) and PY_TYPE_CHECK(y, Integer):
            return (<Integer>x)._and(y)
        return bin_op(x, y, operator.and_)


    cdef _or(Integer self, Integer other):
        cdef Integer x
        x = Integer()
        mpz_ior(x.value, self.value, other.value)
        return x

    def __or__(x, y):
        """
        Return the bitwise or of the integers x and y.

        EXAMPLES:
            sage: n = 8; m = 4
            sage: n.__or__(m)
            12

        On the command line use eval to evaluate the or using the
        caret notation (which is normally exponentiation because
        of the preprocessor).
            sage: eval('n ^ m')
            12
            sage: eval('Integer(8) ^ Integer(4)')
            12
        """
        if PY_TYPE_CHECK(x, Integer) and PY_TYPE_CHECK(y, Integer):
            return (<Integer>x)._or(y)
        return bin_op(x, y, operator.or_)


    def __invert__(self):
        """
        Return the multiplicative interse of self, as a rational number.

        EXAMPLE:
            sage: n = 10
            sage: 1/n
            1/10
            sage: n.__invert__()
            1/10
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
            sage: c = a.inverse_mod(b)
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

        EXAMPLE:
            sage: gcd(-1,1)
            1
            sage: gcd(0,1)
            1
            sage: gcd(0,0)
            0
            sage: gcd(2,2^6)
            2
            sage: gcd(21,2^6)
            1
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

    def test_bit(self, index):
        r"""
        Return the bit at \code{index}.

        EXAMPLES:
            sage: w = 6
            sage: w.str(2)
            '110'
            sage: w.test_bit(2)
            1
            sage: w.test_bit(-1)
            0
        """
        cdef unsigned long int i
        i = index
        cdef Integer x
        x = Integer(self)
        return mpz_tstbit(x.value, i)


ONE = Integer(1)

def integer(x):
    if PY_TYPE_CHECK(x, Integer):
        return x
    return Integer(x)


def LCM_list(v):
    cdef int i, n
    cdef mpz_t z
    cdef Integer w

    n = len(v)

    if n == 0:
        return Integer(1)

    try:
        w = v[0]
        mpz_init_set(z, w.value)

        _sig_on
        for i from 1 <= i < n:
            w = v[i]
            mpz_lcm(z, z, w.value)
        _sig_off
    except TypeError:
        w = Integer(v[0])
        mpz_init_set(z, w.value)

        _sig_on
        for i from 1 <= i < n:
            w = Integer(v[i])
            mpz_lcm(z, z, w.value)
        _sig_off


    w = Integer()
    mpz_set(w.value, z)
    mpz_clear(z)
    return w



def GCD_list(v):
    cdef int i, n
    cdef mpz_t z
    cdef Integer w

    n = len(v)

    if n == 0:
        return Integer(1)

    try:
        w = v[0]
        mpz_init_set(z, w.value)

        _sig_on
        for i from 1 <= i < n:
            w = v[i]
            mpz_gcd(z, z, w.value)
            if mpz_cmp_si(z, 1) == 0:
                _sig_off
                return Integer(1)
        _sig_off
    except TypeError:
        w = Integer(v[0])
        mpz_init_set(z, w.value)

        _sig_on
        for i from 1 <= i < n:
            w = Integer(v[i])
            mpz_gcd(z, z, w.value)
            if mpz_cmp_si(z, 1) == 0:
                _sig_off
                return Integer(1)
        _sig_off


    w = Integer()
    mpz_set(w.value, z)
    mpz_clear(z)
    return w

def make_integer(s):
    r = Integer()
    r._reduce_set(s)
    return r
