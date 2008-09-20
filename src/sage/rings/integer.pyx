r"""
Elements of the ring $\Z$ of integers

AUTHORS:
    -- William Stein (2005): initial version

    -- Gonzalo Tornaria (2006-03-02): vastly improved python/GMP
                                   conversion; hashing
    -- Didier Deshommes <dfdeshom@gmail.com> (2006-03-06): numerous
                                   examples and docstrings
    -- William Stein (2006-03-31): changes to reflect GMP bug fixes
    -- William Stein (2006-04-14): added GMP factorial method (since it's
                                   now very fast).
    -- David Harvey (2006-09-15): added nth_root, exact_log
    -- David Harvey (2006-09-16): attempt to optimise Integer constructor
    -- Rishikesh (2007-02-25): changed quo_rem so that the rem is positive
    -- David Harvey, Martin Albrecht, Robert Bradshaw (2007-03-01):
                                   optimized Integer constructor and
                                   pool
    -- Pablo De Napoli (2007-04-01): multiplicative_order should
                                   return +infinity for non zero
                                   numbers
    -- Robert Bradshaw (2007-04-12): is_perfect_power, Jacobi symbol
                                   (with Kronecker extension). Convert
                                   some methods to use GMP directly
                                   rather than pari, Integer() ->
                                   PY_NEW(Integer)
    -- David Roe (2007-03-21): sped up valuation and is_square, added
                                   val_unit, is_power, is_power_of and
                                   divide_knowing_divisible_by
    -- Robert Bradshaw (2008-03-26): gamma function, multifactorials

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

COERCIONS:
    Returns version of this integer in the multi-precision floating
    real field R.

    sage: n = 9390823
    sage: RR = RealField(200)
    sage: RR(n)
    9.3908230000000000000000000000000000000000000000000000000000e6

"""
#*****************************************************************************
#       Copyright (C) 2004,2006 William Stein <wstein@gmail.com>
#       Copyright (C) 2006 Gonzalo Tornaria <tornaria@math.utexas.edu>
#       Copyright (C) 2006 Didier Deshommes <dfdeshom@gmail.com>
#       Copyright (C) 2007 David Harvey <dmharvey@math.harvard.edu>
#       Copyright (C) 2007 Martin Albrecht <malb@informatik.uni-bremen.de>
#       Copyright (C) 2007,2008 Robert Bradshaw <robertwb@math.washington.edu>
#       Copyright (C) 2007 David Roe <roed314@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
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
include "../ext/python_list.pxi"
include "../ext/python_number.pxi"
include "../ext/python_int.pxi"
include "../structure/coerce.pxi"   # for parent_c
include "../libs/pari/decl.pxi"

cdef extern from "mpz_pylong.h":
    cdef mpz_get_pylong(mpz_t src)
    cdef mpz_get_pyintlong(mpz_t src)
    cdef int mpz_set_pylong(mpz_t dst, src) except -1
    cdef long mpz_pythonhash(mpz_t src)

cdef extern from "convert.h":
    cdef void t_INT_to_ZZ( mpz_t value, long *g )

cdef extern from *:
    cdef void mpz_swap(mpz_t a, mpz_t b)

cdef extern from "math.h":
    cdef double log_c "log" (double)
    cdef double ceil_c "ceil" (double)

from sage.libs.pari.gen cimport gen as pari_gen, PariInstance

cdef class Integer(sage.structure.element.EuclideanDomainElement)


import sage.rings.infinity
import sage.libs.pari.all

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

cdef _digits_internal(mpz_t v,l,int offset,int power_index,power_list,digits):
    """
        Parameters:
            v -- the value whose digits we want to put into the list
            l -- the list to file
            offset -- offset from the beginning of the list that we want to fill at
            power_index -- a measure of size to fill and index to power_list
                            we're filling 1<<(power_index+1) digits
            power_list -- a list of powers of the base, precomputed in method digits
            digits -- a python sequence type with objects to use for digits
                        note that python negative index semantics are relied upon
    """
    cdef mpz_t mpz_res
    cdef mpz_t mpz_quot
    cdef Integer temp
    cdef int v_int
    if power_index == -1:
        if digits is None:
            temp = PY_NEW(Integer)
            mpz_set(temp.value, v)
            l[offset] = temp
        else:
            v_int = mpz_get_si(v)
            l[offset] = digits[v_int]
    else:
        mpz_init(mpz_quot)
        mpz_init(mpz_res)
        temp = power_list[power_index]
        mpz_tdiv_qr(mpz_quot, mpz_res, v, temp.value)
        if mpz_sgn(mpz_res) != 0:
            _digits_internal(mpz_res,l,offset,power_index-1,power_list,digits)
        if mpz_sgn(mpz_quot) != 0:
            _digits_internal(mpz_quot,l,offset+(1<<power_index),power_index-1,power_list,digits)
        mpz_clear(mpz_quot)
        mpz_clear(mpz_res)

arith = None
cdef void late_import():
    global arith
    if arith is None:
        import sage.rings.arith
        arith = sage.rings.arith

MAX_UNSIGNED_LONG = 2 * sys.maxint

# This crashes SAGE:
#  s = 2003^100300000
# The problem is related to realloc moving all the memory
# and returning a pointer to the new block of memory, I think.

from sage.structure.sage_object cimport SageObject
from sage.structure.element cimport EuclideanDomainElement, ModuleElement
from sage.structure.element import  bin_op

import integer_ring
the_integer_ring = integer_ring.ZZ

initialized = False
cdef set_zero_one_elements():
    global the_integer_ring, initialized
    if initialized: return
    the_integer_ring._zero_element = Integer(0)
    the_integer_ring._one_element = Integer(1)
    init_mpz_globals()
    initialized = True
set_zero_one_elements()

cdef zero = the_integer_ring._zero_element
cdef one = the_integer_ring._one_element

def is_Integer(x):
    """
    Return true if x is of the SAGE integer type.

    EXAMPLES:
        sage: is_Integer(2)
        True
        sage: is_Integer(2/1)
        False
        sage: is_Integer(int(2))
        False
        sage: is_Integer(long(2))
        False
        sage: is_Integer('5')
        False
    """
    return PY_TYPE_CHECK(x, Integer)

cdef class IntegerWrapper(Integer):
    """Python classes have problems inheriting from Integer directly, but
    they don't have issues with inheriting from IntegerWrapper."""
    pass

cdef class Integer(sage.structure.element.EuclideanDomainElement):
    r"""
    The \class{Integer} class represents arbitrary precision
    integers.  It derives from the \class{Element} class, so
    integers can be used as ring elements anywhere in SAGE.

    Integer() interprets numbers and strings that begin with 0 as octal
    numbers, and numbers and strings that begin with 0x as hexadecimal
    numbers.

    \begin{notice}
    The class \class{Integer} is implemented in Pyrex, as a wrapper
    of the GMP \code{mpz_t} integer type.
    \end{notice}

    EXAMPLES:
        sage: Integer(010)
        8
        sage: Integer(0x10)
        16
        sage: Integer(10)
        10
        sage: Integer('0x12')
        18
        sage: Integer('012')
        10
    """

    # todo: It would be really nice if we could avoid the __new__ call.
    # It has python calling conventions, and our timing tests indicate the
    # overhead can be significant. The difficulty is that then we can't
    # guarantee that the initialization will be performed exactly once.

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
            sage: ZZ(RR(2.0)^80)
            1208925819614629174706176
            sage: ZZ(QQbar(sqrt(28-10*sqrt(3)) + sqrt(3)))
            5
            sage: ZZ(AA(32).nth_root(5))
            2
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
            sage: ZZ( ZZ(5).digits(3) , 3)
            5
            sage: import numpy
            sage: ZZ(numpy.int64(7^7))
            823543
            sage: ZZ(numpy.ubyte(-7))
            249
            sage: ZZ(True)
            1
            sage: ZZ(False)
            0
            sage: ZZ(1==0)
            0
            sage: ZZ('+10')
            10

            sage: k = GF(2)
            sage: ZZ( (k(0),k(1)), 2)
            2

            sage: t = pari(0*ZZ[x].0 + 3)
            sage: t.type()
            't_POL'
            sage: ZZ(t)
            3
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

        cdef Integer tmp
        cdef char* xs

        if x is None:
            if mpz_sgn(self.value) != 0:
                mpz_set_si(self.value, 0)

        else:
            # First do all the type-check versions; these are fast.

            if PY_TYPE_CHECK(x, Integer):
                set_from_Integer(self, <Integer>x)

            elif PY_TYPE_CHECK(x, bool):
                mpz_set_si(self.value, PyInt_AS_LONG(x))

            elif PyInt_CheckExact(x):
                mpz_set_si(self.value, PyInt_AS_LONG(x))

            elif PyLong_CheckExact(x):
                mpz_set_pylong(self.value, x)

            elif PY_TYPE_CHECK(x, pari_gen):

                if typ((<pari_gen>x).g) == t_INT:
                    t_INT_to_ZZ(self.value, (<pari_gen>x).g)

                else:
                    if typ((<pari_gen>x).g) == t_INTMOD:
                        x = x.lift()
                    # TODO: figure out how to convert to pari integer in base 16 ?

                    # todo: having this "s" variable around here is causing
                    # pyrex to play games with refcount for the None object, which
                    # seems really stupid.

                    try:
                        s = hex(x)
                        base = 16
                    except:
                        s = str(x)
                        base = 10

                    if mpz_set_str(self.value, s, base) != 0:
                        raise TypeError, "Unable to coerce PARI %s to an Integer."%x

            elif PyString_Check(x):
                if base < 0 or base > 36:
                    raise ValueError, "base (=%s) must be between 2 and 36"%base

                xs = x
                if xs[0] == c'+':
                    xs += 1
                if mpz_set_str(self.value, xs, base) != 0:
                    raise TypeError, "unable to convert x (=%s) to an integer"%x

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

            elif (PY_TYPE_CHECK(x, list) or PY_TYPE_CHECK(x, tuple)) and base > 1:
                b = the_integer_ring(base)
                tmp = the_integer_ring(0)
                for i in range(len(x)):
                    tmp += the_integer_ring(x[i])*b**i
                mpz_set(self.value, tmp.value)

            else:
                import numpy
                if isinstance(x, numpy.integer):
                    mpz_set_pylong(self.value, x.__long__())
                else:
                    raise TypeError, "unable to coerce element to an integer"


    def __reduce__(self):
        """
        This is used when pickling integers.

        EXAMPLES:
            sage: n = 5
            sage: t = n.__reduce__(); t
            (<built-in function make_integer>, ('5',))
            sage: t[0](*t[1])
            5
            sage: loads(dumps(n)) == n
            True
        """
        # This single line below took me HOURS to figure out.
        # It is the *trick* needed to pickle pyrex extension types.
        # The trick is that you must put a pure Python function
        # as the first argument, and that function must return
        # the result of unpickling with the argument in the second
        # tuple as input. All kinds of problems happen
        # if we don't do this.
        return sage.rings.integer.make_integer, (self.str(32),)

    cdef _reduce_set(self, s):
        """
        Set this integer from a string in base 32.

        NOTE: Integers are supposed to be immutable, so you
        should not use this function.
        """
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
        return mpz_get_pyintlong(self.value)

    def _im_gens_(self, codomain, im_gens):
        """
        Return the image of self under the map that sends
        the generators of the parent to im_gens.  Since ZZ
        maps canonically in the category of rings, this is just
        the natural coercion.

        EXAMPLES:
            sage: n = -10
            sage: R = GF(17)
            sage: n._im_gens_(R, [R(1)])
            7
        """
        return codomain._coerce_(self)

    cdef _xor(Integer self, Integer other):
        cdef Integer x
        x = PY_NEW(Integer)
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
            return (<Integer>x)._xor(y)
        return bin_op(x, y, operator.xor)


    def __richcmp__(left, right, int op):
        """
        cmp for integers

        EXAMPLES:
            sage: 2 < 3
            True
            sage: 2 > 3
            False
            sage: 2 == 3
            False
            sage: 3 > 2
            True
            sage: 3 < 2
            False

        Canonical coercions are used but non-canonical ones are not.
            sage: 4 == 4/1
            True
            sage: 4 == '4'
            False
        """
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
        z = PY_NEW(Integer)
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
        """
        Return string representation of this integer.

        EXAMPLES:
            sage: n = -5; n.__repr__()
            '-5'
        """
        return self.str()

    def _latex_(self):
        """
        Return latex representation of this integer.  This is just the
        underlying string representation and nothing more.  This
        is called by the latex function.

        EXAMPLES:
            sage: n = -5; n._latex_()
            '-5'
            sage: latex(n)
            -5
        """
        return self.str()

    def _sympy_(self):
        """
        Convert Sage Integer() to SymPy Integer.

        EXAMPLES:
            sage: n = 5; n._sympy_()
            5
            sage: n = -5; n._sympy_()
            -5
        """
        import sympy
        return sympy.sympify(int(self))

    def _mathml_(self):
        """
        Return mathml representation of this integer.

        EXAMPLES:
            sage: mathml(-45)
            <mn>-45</mn>
            sage: (-45)._mathml_()
            '<mn>-45</mn>'
        """
        return '<mn>%s</mn>'%self

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

    def __oct__(self):
        r"""
        Return the digits of self in base 8.

        \note{'0' is \emph{not} prepended to the result like is done
        by the corresponding Python function on int or long. This is
        for efficiency sake---adding and stripping the string wastes
        time; since this function is used for conversions from
        integers to other C-library structures, it is important that
        it be fast.}

        EXAMPLES:
            sage: print oct(Integer(800))
            1440
            sage: print oct(Integer(8))
            10
            sage: print oct(Integer(-50))
            -62
            sage: print oct(Integer(-899))
            -1603
            sage: print oct(Integer(16938402384092843092843098243))
            15535436162247215217705000570203

            Behavior of Sage integers vs. Python integers:
            sage: oct(Integer(10))
            '12'
            sage: oct(int(10))
            '012'
            sage: oct(Integer(-23))
            '-27'
            sage: oct(int(-23))
            '-027'
        """
        return self.str(8)

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

    def bits(self):
        """
        Return the number of bits in self.

        EXAMPLES:
            sage: 500.bits()
            9
            sage: 5.bits()
            3
            sage: 12345.bits() == len(12345.binary())
            True
        """
        return int(mpz_sizeinbase(self.value, 2))

    cdef _digits_naive(self,Integer base,digits):
        """
        This function should have identical semantics to the \code{digits} method of \code{self}.

        Although, the name may imply the opposite, it can actually be quite fast for small values!

        EXAMPLE:
            sage: 9346184691964619824.digits(base=1000) #indirect doc-test
            [824, 619, 964, 691, 184, 346, 9]
        """
        cdef mpz_t mpz_value
        cdef mpz_t mpz_res  # used on one side of the 'if'
        cdef Integer z      # used on the other side of the 'if'
        cdef object l = []

        mpz_init(mpz_value)
        mpz_set(mpz_value, self.value)

        l = []
        # we aim to avoid sage Integer creation if possible
        if digits is None:
            while mpz_cmp_si(mpz_value,0):
                z = PY_NEW(Integer)
                mpz_tdiv_qr(mpz_value, z.value, mpz_value, base.value)
                PyList_Append(l, z)
        else:
            mpz_init(mpz_res)
            while mpz_cmp_si(mpz_value,0):
                mpz_tdiv_qr(mpz_value, mpz_res, mpz_value, base.value)
                PyList_Append(l, digits[mpz_get_si(mpz_res)])
            mpz_clear(mpz_res)

        mpz_clear(mpz_value)

        return l

    def digits(self, base=2, digits=None, padto=0):
        r"""
        Return a list of digits for self in the given base in little
        endian order.

        The return is unspecified if self is a negative number and the digits are given.

        INPUT:
            base -- integer (default: 2)
            digits -- optional indexable object as source for the digits
            padto -- the minimal length of the returned list, sufficient number of
                     zeros are added to make the list minimum that length  (default: 0)

        EXAMPLE:
            sage: 5.digits()
            [1, 0, 1]
            sage: 5.digits(digits=["zero","one"])
            ['one', 'zero', 'one']
            sage: 5.digits(3)
            [2, 1]
            sage: 0.digits(base=10)  # 0 has 0 digits
            []
            sage: 0.digits(base=2)  # 0 has 0 digits
            []
            sage: 10.digits(16,'0123456789abcdef')
            ['a']
            sage: 0.digits(16,'0123456789abcdef')
            []
            sage: 0.digits(16,'0123456789abcdef',padto=1)
            ['0']
            sage: 123.digits(base=10,padto=5)
            [3, 2, 1, 0, 0]
            sage: 123.digits(base=2,padto=3)       # padto is the minimal length
            [1, 1, 0, 1, 1, 1, 1]
            sage: 123.digits(base=2,padto=10,digits=(1,-1))
            [-1, -1, 1, -1, -1, -1, -1, 1, 1, 1]
            sage: a=9939082340; a.digits(10)
            [0, 4, 3, 2, 8, 0, 9, 3, 9, 9]
            sage: a.digits(512)
            [100, 302, 26, 74]
            sage: (-12).digits(10)
            [-2, -1]
            sage: (-12).digits(2)
            [0, 0, -1, -1]

        We support large bases
            sage: n=2^6000
            sage: n.digits(2^3000)
            [0, 0, 1]

            sage: base=3; n=25
            sage: l=n.digits(base)
            sage: # the next relationship should hold for all n,base
            sage: sum(base^i*l[i] for i in range(len(l)))==n
            True
            sage: base=3; n=-30; l=n.digits(base); sum(base^i*l[i] for i in range(len(l)))==n
            True

        Note:  In some cases it is faster to give a digits collection.
        This would be particularly true for computing the digits of a series of small numbers.
        In these cases, the code is careful to allocate as few python objects as reasonably possible.
            sage: digits = range(15)
            sage: l=[ZZ(i).digits(15,digits) for i in range(100)]
            sage: l[16]
            [1, 1]

        This function is comparable to \code{str} for speed.
            sage: n=3^100000
            sage: n.digits(base=10)[-1]  # slightly slower than str
            1
            sage: n=10^10000
            sage: n.digits(base=10)[-1]  # slightly faster than str
            1

        AUTHOR:
            -- Joel B. Mohler significantly rewrote this entire function (2008-03-02)
        """
        if base < 2:
            raise ValueError, "base must be >= 2"

        cdef Integer _base
        cdef Integer self_abs = self
        cdef int power_index = 0
        cdef int i
        cdef size_t s

        if mpz_sgn(self.value) < 0:
            self_abs = -self

        if mpz_sgn(self.value) == 0:
                l = []
        elif base == 2:
            s = mpz_sizeinbase(self.value, 2)
            if digits:
                o = digits[1]
                z = digits[0]
            else:
                o = the_integer_ring._one_element if self > 0 else -the_integer_ring._one_element
                z = the_integer_ring._zero_element
            l = [z]*s
            for i from 0<= i < s:
                if mpz_tstbit(self_abs.value,i):  # mpz_tstbit seems to return 0 for the high-order bit of negative numbers?!
                    l[i] = o
        else:
            s = mpz_sizeinbase(self.value, 2)
            if s > 256:
                _sig_on

            if PY_TYPE_CHECK(base, Integer):
                _base = <Integer>base
            else:
                _base = Integer(base)

            #    It turns out that it is much faster to use simple divison for small numbers.  Some
            # benchmarking helps to define what "small" means in this case.  Here's a very simplistic
            # interpretation of what the benchmark told me, but I think it actually should cut off
            # depending on some ratio of bit-size of self to bit-size of base.
            if s > 2500:
                #   We use a divide and conquer approach (suggested
                # by the prior author, malb?, of the digits method)
                # here: for base b, compute b^2, b^4, b^8,
                # ... (repeated squaring) until you get larger
                # than your number; then compute (n // b^256,
                # n % b ^ 256) (if b^512 > number) to split
                # the number in half and recurse

                # set up digits for optimal access once we get inside the worker functions
                if not digits is None:
                    # list objects have fastest access in the innermost loop
                    if not PyList_CheckExact(digits):
                        digits = [digits[i] for i in range(_base)]
                elif base < 100:  # 100 is arbitrary -- we get a speed boost by pre-allocating digit values in big cases
                    if mpz_sgn(self.value) > 0:
                        digits = [Integer(i) for i in range(_base)]
                    else:
                        # All the digits will be negated in the recursive function.
                        # we'll just compensate for python index semantics
                        digits = [Integer(i) for i in range(-_base,0)]
                        digits[0] = the_integer_ring._zero_element

                power_list = []
                bp = _base
                while bp <= self_abs:
                    power_index += 1
                    #print power_list
                    power_list.append(bp)
                    bp = bp**2

                #   Note that it may appear that the recursive calls to _digit_internal would be assigning
                # list elements i in l for anywhere from 0<=i<(1<<power_index).  However, this is not
                # the case due to the optimization of skipping assigns assigning zero.
                #   Pre-computing the exact number of digits up-front is actually faster (for large values of self)
                # than trimming off trailing zeros after the fact.  It also seems that it would avoid duplicating
                # the list in memory with a list-slice.
                z = the_integer_ring._zero_element if digits is None else digits[0]
                l = [z]*self.ndigits(_base)

                _digits_internal(self.value,l,0,power_index-1,power_list,digits)
            else:
                l = self._digits_naive(_base,digits)

            if s > 256:
                _sig_off

        if digits is None:
            zero = the_integer_ring._zero_element     # value for padding
        else:
            zero = digits[0]                          # value for padding

        return l+[zero]*(padto-len(l))                # padding with zero

    def ndigits(self, base=10):
        """
        Return the number of digits of self expressed in the given base.

        INPUT:
            base -- integer (default: 10)

        EXAMPLES:
            sage: n = 52
            sage: n.ndigits()
            2
            sage: n = -10003
            sage: n.ndigits()
            5
            sage: n = 15
            sage: n.ndigits(2)
            4
            sage: n=1000**1000000+1
            sage: n.ndigits()
            3000001
            sage: n=1000**1000000-1
            sage: n.ndigits()
            3000000
            sage: n=10**10000000-10**9999990
            sage: n.ndigits()
            10000000
        """
        if self == 0:
            return 0
        return self.abs().exact_log(base) + 1

    cdef void set_from_mpz(Integer self, mpz_t value):
        mpz_set(self.value, value)

    cdef mpz_t* get_value(Integer self):
        return &self.value

    cdef void _to_ZZ(self, ZZ_c *z):
        _sig_on
        mpz_to_ZZ(z, &self.value)
        _sig_off

    cdef ModuleElement _add_c_impl(self, ModuleElement right):
        # self and right are guaranteed to be Integers
        cdef Integer x
        x = PY_NEW(Integer)
        mpz_add(x.value, self.value, (<Integer>right).value)
        return x

    cdef ModuleElement _iadd_c_impl(self, ModuleElement right):
        # self and right are guaranteed to be Integers, self safe to mutate
        mpz_add(self.value, self.value, (<Integer>right).value)
        return self

    cdef ModuleElement _sub_c_impl(self, ModuleElement right):
        # self and right are guaranteed to be Integers
        cdef Integer x
        x = PY_NEW(Integer)
        mpz_sub(x.value, self.value, (<Integer>right).value)
        return x

    cdef ModuleElement _isub_c_impl(self, ModuleElement right):
        mpz_sub(self.value, self.value, (<Integer>right).value)
        return self

    cdef ModuleElement _neg_c_impl(self):
        cdef Integer x
        x = PY_NEW(Integer)
        mpz_neg(x.value, self.value)
        return x

    def _r_action(self, s):
        """
        EXAMPLES:
            sage: 8 * [0] #indirect doctest
            [0, 0, 0, 0, 0, 0, 0, 0]
            sage: 8 * 'hi'
            'hihihihihihihihi'
        """
        if isinstance(s, (str, list, tuple)):
            return s*int(self)
        raise TypeError

    def _l_action(self, s):
        """
        EXAMPLES:
            sage: [0] * 8 #indirect doctest
            [0, 0, 0, 0, 0, 0, 0, 0]
            sage: 'hi' * 8
            'hihihihihihihihi'
        """
        if isinstance(s, (str, list, tuple)):
            return int(self)*s
        raise TypeError

    cdef RingElement _mul_c_impl(self, RingElement right):
        # self and right are guaranteed to be Integers
        cdef Integer x
        x = PY_NEW(Integer)
        if mpz_size(self.value) + mpz_size((<Integer>right).value) > 100000:
            # We only use the signal handler (to enable ctrl-c out) when the
            # product might take a while to compute
            _sig_on
            mpz_mul(x.value, self.value, (<Integer>right).value)
            _sig_off
        else:
            mpz_mul(x.value, self.value, (<Integer>right).value)
        return x

    cdef RingElement _imul_c_impl(self, RingElement right):
        if mpz_size(self.value) + mpz_size((<Integer>right).value) > 100000:
            # We only use the signal handler (to enable ctrl-c out) when the
            # product might take a while to compute
            _sig_on
            mpz_mul(self.value, self.value, (<Integer>right).value)
            _sig_off
        else:
            mpz_mul(self.value, self.value, (<Integer>right).value)
        return self

    cdef RingElement _div_c_impl(self, RingElement right):
        r"""
        Computes \frac{a}{b}

        EXAMPLES:
            sage: a = Integer(3) ; b = Integer(4)
            sage: a / b == Rational(3) / 4
            True
            sage: Integer(32) / Integer(32)
            1
        """
        # this is vastly faster than doing it here, since here
        # we can't cimport rationals.
        return the_integer_ring._div(self, right)

    cdef _floordiv(Integer self, Integer other):
        cdef Integer x
        x = PY_NEW(Integer)

        _sig_on
        mpz_fdiv_q(x.value, self.value, other.value)
        _sig_off

        return x


    def __floordiv__(x, y):
        r"""
        Computes the whole part of $\frac{self}{other}$.

        EXAMPLES:
            sage: a = Integer(321) ; b = Integer(10)
            sage: a // b
            32
            sage: z = Integer(-231)
            sage: z // 2
            -116
            sage: z = Integer(231)
            sage: z // 2
            115
            sage: z // -2
            -116
            sage: z // 0
            Traceback (most recent call last):
            ...
            ZeroDivisionError: other must be nonzero
        """
        if PY_TYPE_CHECK(x, Integer) and PY_TYPE_CHECK(y, Integer):
            if not mpz_sgn((<Integer>y).value):
                raise ZeroDivisionError, "other must be nonzero"
            return (<Integer>x)._floordiv(y)
        return bin_op(x, y, operator.floordiv)


    def __pow__(self, n, modulus):
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
            (-1)^(1/3)

        For consistency with Python and MPFR, 0^0 is defined to be 1
        in Sage:
            sage: 0^0
            1


        The base need not be an integer (it can be a builtin
        Python type).
            sage: int(2)^10
            1024
            sage: float(2.5)^10
            9536.7431640625
            sage: 'sage'^3
            'sagesagesage'

        The exponent must fit in a long unless the base is -1, 0, or 1.
            sage: x = 2^100000000000000000000000
            Traceback (most recent call last):
            ...
            RuntimeError: exponent must be at most 2147483647  # 32-bit
            RuntimeError: exponent must be at most 9223372036854775807 # 64-bit
            sage: (-1)^100000000000000000000000
            1

        We raise 2 to various interesting exponents:
            sage: 2^x                # symbolic x
            2^x
            sage: 2^1.5              # real number
            2.82842712474619
            sage: 2^float(1.5)       # python float
            2.8284271247461903
            sage: 2^I                # complex number
            2^I
            sage: f = 2^(sin(x)-cos(x)); f
            2^(sin(x) - cos(x))
            sage: f(3)
            2^(sin(3) - cos(3))

        A symbolic sum:
            sage: x,y,z = var('x,y,z')
            sage: 2^(x+y+z)
            2^(z + y + x)
            sage: 2^(1/2)
            sqrt(2)
            sage: 2^(-1/2)
            1/sqrt(2)

        TESTS:
            sage: complex(0,1)^2
            (-1+0j)
        """
        if modulus is not None:
            #raise RuntimeError, "__pow__ dummy argument ignored"
            from sage.rings.integer_mod import Mod
            return Mod(self, modulus) ** n


        if not PY_TYPE_CHECK(self, Integer):
            if isinstance(self, str):
                return self * n
            else:
                return self ** int(n)

        cdef Integer _self = <Integer>self
        cdef long nn

        try:
            nn = PyNumber_Index(n)
        except TypeError:
            try:
                s = parent_c(n)(self)
                return s**n
            except AttributeError:
                raise TypeError, "exponent (=%s) must be an integer.\nCoerce your numbers to real or complex numbers first."%n

        except OverflowError:
            if mpz_cmp_si(_self.value, 1) == 0:
                return self
            elif mpz_cmp_si(_self.value, 0) == 0:
                return self
            elif mpz_cmp_si(_self.value, -1) == 0:
                return self if n % 2 else -self
            raise RuntimeError, "exponent must be at most %s" % sys.maxint

        if nn == 0:
            return one

        cdef Integer x = PY_NEW(Integer)

        _sig_on
        mpz_pow_ui(x.value, (<Integer>self).value, nn if nn > 0 else -nn)
        _sig_off

        if nn < 0:
            return ~x
        else:
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
        if (mpz_sgn(self.value) < 0) and not (n & 1):
            raise ValueError, "cannot take even root of negative number"
        cdef Integer x
        cdef bint is_exact
        x = PY_NEW(Integer)
        _sig_on
        is_exact = mpz_root(x.value, self.value, n)
        _sig_off

        if report_exact:
            return x, is_exact
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
           sage: RR(log(RR(x), 3))
           100000.000000000
           sage: RR(log(RR(x + 100000), 3))
           100000.000000000

           sage: x.exact_log(3)
           100000
           sage: (x+1).exact_log(3)
           100000
           sage: (x-1).exact_log(3)
           99999

           sage: x.exact_log(2.5)
           Traceback (most recent call last):
           ...
           TypeError: Attempt to coerce non-integral RealNumber to Integer
        """
        m = Integer(m)
        if mpz_sgn(self.value) <= 0:
            raise ValueError, "self must be positive"
        if m < 2:
            raise ValueError, "m must be at least 2"

        if m <= 256:
            # if the base m is at most 256, we can use mpz_sizeinbase
            # to get the following guess which is either the exact
            # log, or 1+ the exact log
            guess = Integer(mpz_sizeinbase(self.value, m)) - 1
            # if the base is 2, the guess is always correct
            if m == 2:
                return guess

            # otherwise, we need to compare self and m^guess
            # this is time-consuming for large integers, so we start by
            # doing a rough comparison using interval arithmetic
            # (suggested by David Harvey and Carl Witty)
            # "for randomly distributed integers, the chance of this
            # interval-based comparison failing is absurdly low"
            import real_mpfi
            approx_compare = real_mpfi.RIF(m)**guess
            if self > approx_compare:
                return guess
            elif self < approx_compare:
                return guess - one
            # if we reach this point, we're in an absurdly low-probability
            # case; we "manually" compare self and m^guess
            compare = Integer(m)**guess
            if self >= compare:
                return guess
            else:
                return guess - one

        # if we are here, then the base m is bigger than 256
        # TODO: optimize this
        import real_mpfr
        R = real_mpfr.RealField(53)
        guess = R(self).log(base = m).floor()
        power = Integer(m) ** guess

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

    def prime_to_m_part(self, m):
        """
        Returns the prime-to-m part of self, i.e., the largest divisor
        of self that is coprime to m.

        INPUT:
            m -- Integer
        OUTPUT:
            Integer

        EXAMPLES:
            sage: z = 43434
            sage: z.prime_to_m_part(20)
            21717
        """
        late_import()
        return sage.rings.arith.prime_to_m_part(self, m)

    def prime_divisors(self):
        """
        The prime divisors of self, sorted in increasing order.  If n
        is negative, we do *not* include -1 among the prime divisors, since -1 is
        not a prime number.

        EXAMPLES:
            sage: a = 1; a.prime_divisors()
            []
            sage: a = 100; a.prime_divisors()
            [2, 5]
            sage: a = -100; a.prime_divisors()
            [2, 5]
            sage: a = 2004; a.prime_divisors()
            [2, 3, 167]
        """
        late_import()
        return sage.rings.arith.prime_divisors(self)

    prime_factors = prime_divisors

    def divisors(self):
        """
        Returns a list of all positive integer divisors
        of the integer self.

        EXAMPLES:
            sage: a = -3; a.divisors()
            [1, 3]
            sage: a = 6; a.divisors()
            [1, 2, 3, 6]
            sage: a = 28; a.divisors()
            [1, 2, 4, 7, 14, 28]
            sage: a = 2^5; a.divisors()
            [1, 2, 4, 8, 16, 32]
            sage: a = 100; a.divisors()
            [1, 2, 4, 5, 10, 20, 25, 50, 100]
            sage: a = 1; a.divisors()
            [1]
            sage: a = 0; a.divisors()
            Traceback (most recent call last):
            ...
            ValueError: n must be nonzero
            sage: a = 2^3 * 3^2 * 17; a.divisors()
            [1, 2, 3, 4, 6, 8, 9, 12, 17, 18, 24, 34, 36, 51, 68, 72, 102, 136, 153, 204, 306, 408, 612, 1224]
        """
        late_import()
        return sage.rings.arith.divisors(self)

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
        x = PY_NEW(Integer)
        mpz_abs(x.value, self.value)
        return x

    def __mod__(self, modulus):
        r"""
        Returns self modulo the modulus.

        EXAMPLES:
            sage: z = 43
            sage: z % 2
            1
            sage: z % 0
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Integer modulo by zero
            sage: -5 % 7
            2
            sage: -5 % -7
            -5
            sage: 5 % -7
            -2
         """
        cdef Integer _modulus, _self
        _modulus = integer(modulus)
        if not _modulus:
            raise ZeroDivisionError, "Integer modulo by zero"
        _self = integer(self)

        cdef Integer x
        x = PY_NEW(Integer)

        _sig_on
        mpz_mod(x.value, _self.value, _modulus.value)
        # GMP always returns a non-negative integer, so we
        # subtract if we got something positive.
        if mpz_sgn(_modulus.value) == -1 and mpz_sgn(x.value):
            mpz_add(x.value, x.value, _modulus.value)
        _sig_off

        return x

    def quo_rem(self, other):
        """
        Returns the quotient and the remainder of self divided by
        other. Note that the remainder returned is always either
        zero or of the same sign as other.

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
            (-116, -1)
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
        q = PY_NEW(Integer)
        r = PY_NEW(Integer)

        _sig_on
        mpz_fdiv_qr(q.value, r.value, _self.value, _other.value)
        _sig_off

        return q, r

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

        x = PY_NEW(Integer)

        _sig_on
        mpz_powm(x.value, self.value, _exp.value, _mod.value)
        _sig_off

        return x

    def rational_reconstruction(self, Integer m):
        """
        Return the rational reconstruction of this integer modulo m, i.e.,
        the unique (if it exists) rational number that reduces to self
        modulo m and whose numerator and denominator is bounded by
        sqrt(m/2).

        EXAMPLES:
            sage: (3/7)%100
            29
            sage: (29).rational_reconstruction(100)
            3/7
        """
        import rational
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
            sage: z.powermodm_ui(2^32-2, 14)
            2
            sage: z.powermodm_ui(2^32-1, 14)
            Traceback (most recent call last):                              # 32-bit
            ...                                                             # 32-bit
            OverflowError: exp (=4294967295) must be <= 4294967294          # 32-bit
            8              # 64-bit
            sage: z.powermodm_ui(2^65, 14)
            Traceback (most recent call last):
            ...
            OverflowError: exp (=36893488147419103232) must be <= 4294967294  # 32-bit
            OverflowError: exp (=36893488147419103232) must be <= 18446744073709551614     # 64-bit
        """
        if exp < 0:
            raise ValueError, "exp (=%s) must be nonnegative"%exp
        elif exp > MAX_UNSIGNED_LONG:
            raise OverflowError, "exp (=%s) must be <= %s"%(exp, MAX_UNSIGNED_LONG)
        cdef Integer x, _mod
        _mod = Integer(mod)
        x = PY_NEW(Integer)

        _sig_on
        mpz_powm_ui(x.value, self.value, exp, _mod.value)
        _sig_off

        return x

    def __int__(self):
        """
        Return the Python int (or long) corresponding to this SAGE integer.

        EXAMPLES:
            sage: n = 920938; n
            920938
            sage: int(n)
            920938
            sage: type(n.__int__())
            <type 'int'>
            sage: n = 99028390823409823904823098490238409823490820938; n
            99028390823409823904823098490238409823490820938
            sage: type(n.__int__())
            <type 'long'>
        """
        # TODO -- this crashes on sage.math, since it is evidently written incorrectly.
        return mpz_get_pyintlong(self.value)
        #return int(mpz_get_pylong(self.value))

    def __long__(self):
        """
        Return long integer corresponding to this SAGE integer.

        EXAMPLES:
            sage: n = 9023408290348092849023849820934820938490234290; n
            9023408290348092849023849820934820938490234290
            sage: long(n)
            9023408290348092849023849820934820938490234290L
            sage: n = 920938; n
            920938
            sage: long(n)
            920938L
            sage: n.__long__()
            920938L
        """
        return mpz_get_pylong(self.value)

    def __float__(self):
        """
        Return double precision floating point representation of this integer.

        EXAMPLES:
            sage: n = Integer(17); float(n)
            17.0
            sage: n = Integer(902834098234908209348209834092834098); float(n)
            9.0283409823490813e+35
            sage: n = Integer(-57); float(n)
            -57.0
            sage: n.__float__()
            -57.0
            sage: type(n.__float__())
            <type 'float'>
        """
        return mpz_get_d(self.value)

    def _rpy_(self):
        """
        Returns int(self) so that rpy can convert self into an object it
        knows how to work with.

        EXAMPLES:
            sage: n = 100
            sage: n._rpy_()
            100
            sage: type(n._rpy_())
            <type 'int'>
        """
        return self.__int__()

    def __hash__(self):
        """
        Return the hash of this integer.

        This agrees with the Python hash of the corresponding
        Python int or long.

        EXAMPLES:
            sage: n = -920384; n.__hash__()
            -920384
            sage: hash(int(n))
            -920384
            sage: n = -920390823904823094890238490238484; n.__hash__()
            -873977844         # 32-bit
            6874330978542788722   # 64-bit
            sage: hash(long(n))
            -873977844            # 32-bit
            6874330978542788722   # 64-bit
        """
        return mpz_pythonhash(self.value)

    cdef hash_c(self):
        """
        A C version of the __hash__ function.
        """
        return mpz_pythonhash(self.value)

    def _factor_trial_division(self, limit):
        """
        Return partial factorization of self obtained using trial
        division for all primes up to limit, where limit must fit
        in a signed int.

        INPUT:
            limit -- integer that fits in a signed int

        ALGORITHM: Uses Pari.

        EXAMPLES:
            sage: n = 920384092842390423848290348203948092384082349082
            sage: n._factor_trial_division(1000)
            2 * 11 * 41835640583745019265831379463815822381094652231
            sage: n._factor_trial_division(2^30)
            2 * 11 * 1531 * 27325696005058797691594630609938486205809701
        """
        import sage.structure.factorization as factorization
        F = self._pari_().factor(limit)
        B, e = F
        return factorization.Factorization([(the_integer_ring(B[i]), the_integer_ring(e[i])) for i in range(len(B))])

    def factor(self, algorithm='pari', proof=True, limit=None):
        """
        Return the prime factorization of the integer as a list of
        pairs $(p,e)$, where $p$ is prime and $e$ is a positive integer.

        INPUT:
            algorithm -- string
                 * 'pari' -- (default)  use the PARI c library
                 * 'kash' -- use KASH computer algebra system (requires
                             the optional kash package be installed)
            proof -- bool (default: True) whether or not to prove primality
                    of each factor (only applicable for PARI).
            limit -- int or None (default: None) if limit is given it
                     must fit in a signed int, and the factorization
                     is done using trial divsion and primits up to limit.

        EXAMPLES:
            sage: n = 2^100 - 1; n.factor()
            3 * 5^3 * 11 * 31 * 41 * 101 * 251 * 601 * 1801 * 4051 * 8101 * 268501

        We using proof=False, which doesn't prove correctness of the
        primes that appear in the factorization:
            sage: n = 920384092842390423848290348203948092384082349082
            sage: n.factor(proof=False)
            2 * 11 * 1531 * 4402903 * 10023679 * 619162955472170540533894518173
            sage: n.factor(proof=True)
            2 * 11 * 1531 * 4402903 * 10023679 * 619162955472170540533894518173

        We factor using trial division only:
            sage: n.factor(limit=1000)
            2 * 11 * 41835640583745019265831379463815822381094652231
        """
        if limit is not None:
            return self._factor_trial_division(limit)
        import sage.rings.integer_ring
        return sage.rings.integer_ring.factor(self, algorithm=algorithm, proof=proof)

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
        cdef bint t
        cdef Integer _n
        _n = Integer(n)
        if mpz_sgn(self.value) == 0:
            return mpz_sgn(_n.value) == 0
        _sig_on
        t = mpz_divisible_p(_n.value, self.value)
        _sig_off
        return t

    cdef RingElement _valuation(Integer self, Integer p):
        r"""
        Return the p-adic valuation of self.

        We do not require that p be prime, but it must be at least 2.
        For more documentation see \code{valuation}

        AUTHOR:
            -- David Roe (3/31/07)
        """
        if mpz_sgn(self.value) == 0:
            return sage.rings.infinity.infinity
        if mpz_cmp_ui(p.value, 2) < 0:
            raise ValueError, "You can only compute the valuation with respect to a integer larger than 1."
        cdef Integer v
        v = PY_NEW(Integer)
        cdef mpz_t u
        mpz_init(u)
        _sig_on
        mpz_set_ui(v.value, mpz_remove(u, self.value, p.value))
        _sig_off
        mpz_clear(u)
        return v

    cdef object _val_unit(Integer self, Integer p):
        r"""
        Returns a pair: the p-adic valuation of self, and the p-adic unit of self.

        We do not require the p be prime, but it must be at least 2.
        For more documentation see \code{val_unit}

        AUTHOR:
            -- David Roe (3/31/07)
        """
        cdef Integer v, u
        if mpz_cmp_ui(p.value, 2) < 0:
            raise ValueError, "You can only compute the valuation with respect to a integer larger than 1."
        if self == 0:
            u = ONE
            Py_INCREF(ONE)
            return (sage.rings.infinity.infinity, u)
        v = PY_NEW(Integer)
        u = PY_NEW(Integer)
        _sig_on
        mpz_set_ui(v.value, mpz_remove(u.value, self.value, p.value))
        _sig_off
        return (v, u)

    def valuation(self, p):
        """
        Return the p-adic valuation of self.

        INPUT:
            p -- an integer at least 2.

        EXAMPLE:
        sage: n = 60
        sage: n.valuation(2)
        2
        sage: n.valuation(3)
        1
        sage: n.valuation(7)
        0
        sage: n.valuation(1)
        Traceback (most recent call last):
        ...
        ValueError: You can only compute the valuation with respect to a integer larger than 1.

        We do not require that p is a prime:
        sage: (2^11).valuation(4)
        5
        """
        return self._valuation(Integer(p))

    def val_unit(self, p):
        r"""
        Returns a pair: the p-adic valuation of self, and the p-adic unit of self.

        INPUT:
            p -- an integer at least 2.

        OUTPUT:
            v_p(self) -- the p-adic valuation of \code{self}
            u_p(self) -- \code{self} / $p^{v_p(\code{self})}$

        EXAMPLE:
        sage: n = 60
        sage: n.val_unit(2)
        (2, 15)
        sage: n.val_unit(3)
        (1, 20)
        sage: n.val_unit(7)
        (0, 60)
        sage: (2^11).val_unit(4)
        (5, 2)
        sage: 0.val_unit(2)
        (+Infinity, 1)
        """
        return self._val_unit(Integer(p))

    def ord(self, p=None):
        """Synonym for valuation

        EXAMPLES:
        sage: n=12
        sage: n.ord(3)
        1
        """
        return self.valuation(p)

    cdef Integer _divide_knowing_divisible_by(Integer self, Integer right):
        r"""
        Returns the integer self / right when self is divisible by right.

        If self is not divisible by right, the return value is undefined, and may not even be close to self/right.
        For more documentation see \code{divide_knowing_divisible_by}

        AUTHOR:
            -- David Roe (3/31/07)
        """
        if mpz_cmp_ui(right.value, 0) == 0:
            raise ZeroDivisionError
        cdef Integer x
        x = PY_NEW(Integer)
        if mpz_size(self.value) + mpz_size((<Integer>right).value) > 100000:
            # Only use the signal handler (to enable ctrl-c out) when the
            # quotient might take a while to compute
            _sig_on
            mpz_divexact(x.value, self.value, right.value)
            _sig_off
        else:
            mpz_divexact(x.value, self.value, right.value)
        return x

    def divide_knowing_divisible_by(self, right):
        r"""
        Returns the integer self / right when self is divisible by right.

        If self is not divisible by right, the return value is undefined,
        and may not even be close to self/right for multi-word integers.

        EXAMPLES:
            sage: a = 8; b = 4
            sage: a.divide_knowing_divisible_by(b)
            2
            sage: (100000).divide_knowing_divisible_by(25)
            4000
            sage: (100000).divide_knowing_divisible_by(26) # close (random)
            3846

      However, often it's way off.

            sage: a = 2^70; a
            1180591620717411303424
            sage: a // 11  # floor divide
            107326510974310118493
            sage: a.divide_knowing_divisible_by(11) # way off and possibly random
            43215361478743422388970455040

        """
        return self._divide_knowing_divisible_by(right)

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
        z = PY_NEW(Integer)
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
        if mpz_sgn(self.value) < 0:
            raise ValueError, "factorial -- self = (%s) must be nonnegative"%self

        if not mpz_fits_uint_p(self.value):
            raise ValueError, "factorial not implemented for n >= 2^32.\nThis is probably OK, since the answer would have billions of digits."

        cdef Integer z = PY_NEW(Integer)

        _sig_on
        mpz_fac_ui(z.value, mpz_get_ui(self.value))
        _sig_off

        return z

    def multifactorial(self, int k):
        r"""
        Computes the k-th factorial $n!^{(k)}$ of self. For k=1 this is the
        standard factorial, and for k greater than one it is the product of
        every k-th terms down from self to k. The recursive definition is used
        to extend this function to the negative integers.

        EXAMPLES:
            sage: 5.multifactorial(1)
            120
            sage: 5.multifactorial(2)
            15
            sage: 23.multifactorial(2)
            316234143225
            sage: prod([1..23, step=2])
            316234143225
            sage: (-29).multifactorial(7)
            1/2640
        """
        if k <= 0:
            raise ValueError, "multifactorial only defined for positive values of k"

        if not mpz_fits_sint_p(self.value):
            raise ValueError, "multifactorial not implemented for n >= 2^32.\nThis is probably OK, since the answer would have billions of digits."

        cdef int n = mpz_get_si(self.value)

        # base case
        if 0 < n < k:
            return ONE

        # easy to calculate
        elif n % k == 0:
            factorial = Integer(n/k).factorial()
            if k == 2:
                return factorial << (n/k)
            else:
                return factorial * Integer(k)**(n/k)

        # negative base case
        elif -k < n < 0:
            return ONE / (self+k)

        # reflection case
        elif n < -k:
            if (n/k) % 2:
                sign = -ONE
            else:
                sign = ONE
            return sign / Integer(-k-n).multifactorial(k)

        # compute the actual product, optimizing the number of large multiplications
        cdef int i,j

        # we need (at most) log_2(#factors) concurrent sub-products
        cdef int prod_count = <int>ceil_c(log_c(n/k+1)/log_c(2))
        cdef mpz_t* sub_prods = <mpz_t*>malloc(prod_count * sizeof(mpz_t))
        if sub_prods == NULL:
            raise MemoryError
        for i from 0 <= i < prod_count:
            mpz_init(sub_prods[i])

        _sig_on

        cdef residue = n % k
        cdef int tip = 0
        for i from 1 <= i <= n//k:
            mpz_set_ui(sub_prods[tip], k*i + residue)
            # for the i-th terms we use the bits of i to calculate how many
            # times we need to multiply "up" the stack of sub-products
            for j from 0 <= j < 32:
                if i & (1 << j):
                    break
                tip -= 1
                mpz_mul(sub_prods[tip], sub_prods[tip], sub_prods[tip+1])
            tip += 1
        cdef int last = tip-1
        for tip from last > tip >= 0:
            mpz_mul(sub_prods[tip], sub_prods[tip], sub_prods[tip+1])

        _sig_off

        cdef Integer z = PY_NEW(Integer)
        mpz_swap(z.value, sub_prods[0])

        for i from 0 <= i < prod_count:
            mpz_clear(sub_prods[i])
        free(sub_prods)

        return z


    def gamma(self):
        r"""
        The gamma function on integers is the factorial function (shifted by
        one) on positive integers, and $\pm \infty$ on non-positive integers.

        EXAMPLES:
            sage: gamma(5)
            24
            sage: gamma(0)
            -Infinity
            sage: gamma(-1)
            +Infinity
            sage: gamma(-2^150)
            -Infinity
        """
        if mpz_sgn(self.value) > 0:
            return (self-ONE).factorial()
        else:
            from sage.rings.infinity import infinity
            if mpz_even_p(self.value):
                return -infinity
            else:
                return infinity

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
        return mpz_cmp_si(self.value, 1) == 0

    def __nonzero__(self):
        r"""
        Returns \code{True} if the integers is not $0$, otherwise \code{False}.

        EXAMPLES:
            sage: Integer(1).is_zero()
            False
            sage: Integer(0).is_zero()
            True
        """
        return mpz_sgn(self.value) != 0

    def is_integral(self):
        """
        Return \code{True} since integers are integral, i.e., satisfy
        a monic polynomial with integer coefficients.

        EXAMPLES:
            sage: Integer(3).is_integral()
            True
        """
        return True

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
        return mpz_cmp_si(self.value, -1) == 0 or mpz_cmp_si(self.value, 1) == 0

    def is_square(self):
        r"""
        Returns \code{True} if self is a perfect square

        EXAMPLES:
        sage: Integer(4).is_square()
        True
        sage: Integer(41).is_square()
        False
        """
        return mpz_perfect_square_p(self.value)

    def is_power(self):
        r"""
        Returns \code{True} if self is a perfect power, ie if there
        exist integers a and b, $b > 1$ with $self = a^b$.

        EXAMPLES:
        sage: Integer(-27).is_power()
        True
        sage: Integer(12).is_power()
        False
        """
        return mpz_perfect_power_p(self.value)

    cdef bint _is_power_of(Integer self, Integer n):
        r"""
        Returns a non-zero int if there is an integer b with $self = n^b$.

        For more documentation see \code{is_power_of}

        AUTHOR:
            -- David Roe (3/31/07)
        """
        cdef int a
        cdef unsigned long b, c
        cdef mpz_t u, sabs, nabs
        a = mpz_cmp_ui(n.value, 2)
        if a <= 0: # n <= 2
            if a == 0: # n == 2
                if mpz_popcount(self.value) == 1: #number of bits set in self == 1
                    return 1
                else:
                    return 0
            a = mpz_cmp_si(n.value, -2)
            if a >= 0: # -2 <= n < 2:
                a = mpz_get_si(n.value)
                if a == 1: # n == 1
                    if mpz_cmp_ui(self.value, 1) == 0: # Only 1 is a power of 1
                        return 1
                    else:
                        return 0
                elif a == 0: # n == 0
                    if mpz_cmp_ui(self.value, 0) == 0 or mpz_cmp_ui(self.value, 1) == 0: # 0^0 = 1, 0^x = 0
                        return 1
                    else:
                        return 0
                elif a == -1: # n == -1
                    if mpz_cmp_ui(self.value, 1) == 0 or mpz_cmp_si(self.value, -1) == 0: # 1 and -1 are powers of -1
                        return 1
                    else:
                        return 0
                elif a == -2: # n == -2
                    mpz_init(sabs)
                    mpz_abs(sabs, self.value)
                    if mpz_popcount(sabs) == 1: # number of bits set in |self| == 1
                        b = mpz_scan1(sabs, 0) % 2 # b == 1 if |self| is an odd power of 2, 0 if |self| is an even power
                        mpz_clear(sabs)
                        if (b == 1 and mpz_cmp_ui(self.value, 0) < 0) or (b == 0 and mpz_cmp_ui(self.value, 0) > 0):
                            # An odd power of -2 is negative, an even power must be positive.
                            return 1
                        else: # number of bits set in |self| is not 1, so self cannot be a power of -2
                            return 0
                    else: # |self| is not a power of 2, so self cannot be a power of -2
                        return 0
            else: # n < -2
                mpz_init(nabs)
                mpz_neg(nabs, n.value)
                if mpz_popcount(nabs) == 1: # |n| = 2^k for k >= 2.  We special case this for speed
                    mpz_init(sabs)
                    mpz_abs(sabs, self.value)
                    if mpz_popcount(sabs) == 1: # |self| = 2^L for some L >= 0.
                        b = mpz_scan1(sabs, 0) # the bit that self is set at
                        c = mpz_scan1(nabs, 0) # the bit that n is set at
                        # Having obtained b and c, we're done with nabs and sabs (on this branch anyway)
                        mpz_clear(nabs)
                        mpz_clear(sabs)
                        if b % c == 0: # Now we know that |self| is a power of |n|
                            b = (b // c) % 2 # Whether b // c is even or odd determines whether (-2^c)^(b // c) is positive or negative
                            a = mpz_cmp_ui(self.value, 0)
                            if b == 0 and a > 0 or b == 1 and a < 0:
                                # These two cases are that b // c is even and self positive, or b // c is odd and self negative
                                return 1
                            else: # The sign of self is wrong
                                return 0
                        else: # Since |self| is not a power of |n|, self cannot be a power of n
                            return 0
                    else: # self is not a power of 2, and thus cannot be a power of n, which is a power of 2.
                        mpz_clear(nabs)
                        mpz_clear(sabs)
                        return 0
                else: # |n| is not a power of 2, so we use mpz_remove
                    mpz_init(u)
                    _sig_on
                    b = mpz_remove(u, self.value, nabs)
                    _sig_off
                    # Having obtained b and u, we're done with nabs
                    mpz_clear(nabs)
                    if mpz_cmp_ui(u, 1) == 0: # self is a power of |n|
                        mpz_clear(u)
                        if b % 2 == 0: # an even power of |n|, and since self > 0, this means that self is a power of n
                            return 1
                        else:
                            return 0
                    elif mpz_cmp_si(u, -1) == 0: # -self is a power of |n|
                        mpz_clear(u)
                        if b % 2 == 1: # an odd power of |n|, and thus self is a power of n
                            return 1
                        else:
                            return 1
                    else: # |self| is not a power of |n|, so self cannot be a power of n
                        mpz_clear(u)
                        return 0
        elif mpz_popcount(n.value) == 1: # n > 2 and in fact n = 2^k for k >= 2
            if mpz_popcount(self.value) == 1: # since n is a power of 2, so must self be.
                if mpz_scan1(self.value, 0) % mpz_scan1(n.value, 0) == 0: # log_2(self) is divisible by log_2(n)
                    return 1
                else:
                    return 0
            else: # self is not a power of 2, and thus not a power of n
                return 0
        else: # n > 2, but not a power of 2, so we use mpz_remove
            mpz_init(u)
            _sig_on
            mpz_remove(u, self.value, n.value)
            _sig_off
            a = mpz_cmp_ui(u, 1)
            mpz_clear(u)
            if a == 0:
                return 1
            else:
                return 0

    def is_power_of(Integer self, n):
        r"""
        Returns \code{True} if there is an integer b with $self = n^b$.

        EXAMPLES:
            sage: Integer(64).is_power_of(4)
            True
            sage: Integer(64).is_power_of(16)
            False

        TESTS:

            sage: Integer(-64).is_power_of(-4)
            True
            sage: Integer(-32).is_power_of(-2)
            True
            sage: Integer(1).is_power_of(1)
            True
            sage: Integer(-1).is_power_of(-1)
            True
            sage: Integer(0).is_power_of(1)
            False
            sage: Integer(0).is_power_of(0)
            True
            sage: Integer(1).is_power_of(0)
            True
            sage: Integer(1).is_power_of(8)
            True
            sage: Integer(-8).is_power_of(2)
            False

        NOTES:

        For large integers self, is_power_of() is faster than is_power().
        The following examples gives some indication of how much faster.

            sage: b = lcm(range(1,10000))
            sage: b.exact_log(2)
            14446
            sage: t=cputime()
            sage: for a in range(2, 1000): k = b.is_power()
            sage: cputime(t)      # random
            0.53203299999999976
            sage: t=cputime()
            sage: for a in range(2, 1000): k = b.is_power_of(2)
            sage: cputime(t)      # random
            0.0
            sage: t=cputime()
            sage: for a in range(2, 1000): k = b.is_power_of(3)
            sage: cputime(t)      # random
            0.032002000000000308

            sage: b = lcm(range(1, 1000))
            sage: b.exact_log(2)
            1437
            sage: t=cputime()
            sage: for a in range(2, 10000): k = b.is_power() # note that we change the range from the example above
            sage: cputime(t)      # random
            0.17201100000000036
            sage: t=cputime(); TWO=int(2)
            sage: for a in range(2, 10000): k = b.is_power_of(TWO)
            sage: cputime(t)      # random
            0.0040000000000000036
            sage: t=cputime()
            sage: for a in range(2, 10000): k = b.is_power_of(3)
            sage: cputime(t)      # random
            0.040003000000000011
            sage: t=cputime()
            sage: for a in range(2, 10000): k = b.is_power_of(a)
            sage: cputime(t)      # random
            0.02800199999999986
        """
        if not PY_TYPE_CHECK(n, Integer):
            n = Integer(n)
        return self._is_power_of(n)

    def is_prime_power(self, flag=0):
        r"""
        Returns True if $x$ is a prime power, and False otherwise.

        INPUT:
            flag (for primality testing) -- int
                    0 (default): use a combination of algorithms.
                    1: certify primality using the Pocklington-Lehmer Test.
                    2: certify primality using the APRCL test.

        EXAMPLES:
            sage: (-10).is_prime_power()
            False
            sage: (10).is_prime_power()
            False
            sage: (64).is_prime_power()
            True
            sage: (3^10000).is_prime_power()
            True
            sage: (10000).is_prime_power(flag=1)
            False
        """
        if self.is_zero():
            return False
        elif self.is_one():
            return True
        elif mpz_sgn(self.value) < 0:
            return False
        if self.is_prime():
            return True
        if not self.is_perfect_power():
            return False
        k, g = self._pari_().ispower()
        if not k:
            raise RuntimeError, "inconsistent results between GMP and pari"
        return g.isprime(flag=flag)

    def is_prime(self):
        r"""
        Returns \code{True} if self is prime.

        NOTE: Integer primes are by definition \emph{positive}!
        This is different than Magma, but the same as in Pari.
        See also the \code{is_irreducible()} method.

        EXAMPLES:
            sage: z = 2^31 - 1
            sage: z.is_prime()
            True
            sage: z = 2^31
            sage: z.is_prime()
            False
            sage: z = 7
            sage: z.is_prime()
            True
            sage: z = -7
            sage: z.is_prime()
            False
            sage: z.is_irreducible()
            True
        """
        return bool(self._pari_().isprime())

    def is_irreducible(self):
        r"""
        Returns \code{True} if self is irreducible, i.e. +/- prime

        EXAMPLES:
            sage: z = 2^31 - 1
            sage: z.is_irreducible()
            True
            sage: z = 2^31
            sage: z.is_irreducible()
            False
            sage: z = 7
            sage: z.is_irreducible()
            True
            sage: z = -7
            sage: z.is_irreducible()
            True
        """
        n = self if self>=0 else -self
        return bool(n._pari_().isprime())

    def is_pseudoprime(self):
        r"""
        Returns \code{True} if self is a pseudoprime

        EXAMPLES:
            sage: z = 2^31 - 1
            sage: z.is_pseudoprime()
            True
            sage: z = 2^31
            sage: z.is_pseudoprime()
            False
        """
        return bool(self._pari_().ispseudoprime())

    def is_perfect_power(self):
        r"""
        Returns \code{True} if self is a perfect power.

        EXAMPLES:
            sage: z = 8
            sage: z.is_perfect_power()
            True
            sage: z = 144
            sage: z.is_perfect_power()
            True
            sage: z = 10
            sage: z.is_perfect_power()
            False
        """
        return mpz_perfect_power_p(self.value)

    def jacobi(self, b):
        r"""
        Calculate the Jacobi symbol $\left(\frac{self}{b}\right)$.

        EXAMPLES:
            sage: z = -1
            sage: z.jacobi(17)
            1
            sage: z.jacobi(19)
            -1
            sage: z.jacobi(17*19)
            -1
            sage: (2).jacobi(17)
            1
            sage: (3).jacobi(19)
            -1
            sage: (6).jacobi(17*19)
            -1
            sage: (6).jacobi(33)
            0
            sage: a = 3; b = 7
            sage: a.jacobi(b) == -b.jacobi(a)
            True
        """
        cdef long tmp
        if PY_TYPE_CHECK(b, int):
            tmp = b
            if (tmp & 1) == 0:
                raise ValueError, "Jacobi symbol not defined for even b."
            return mpz_kronecker_si(self.value, tmp)
        if not PY_TYPE_CHECK(b, Integer):
            b = Integer(b)
        if mpz_even_p((<Integer>b).value):
            raise ValueError, "Jacobi symbol not defined for even b."
        return mpz_jacobi(self.value, (<Integer>b).value)

    def kronecker(self, b):
        r"""
        Calculate the Kronecker symbol
        $\left(\frac{self}{b}\right)$ with the Kronecker extension
        $(self/2)=(2/self)$ when self odd, or $(self/2)=0$ when $self$ even.

        EXAMPLES:
            sage: z = 5
            sage: z.kronecker(41)
            1
            sage: z.kronecker(43)
            -1
            sage: z.kronecker(8)
            -1
            sage: z.kronecker(15)
            0
            sage: a = 2; b = 5
            sage: a.kronecker(b) == b.kronecker(a)
            True
        """
        if PY_TYPE_CHECK(b, int):
            return mpz_kronecker_si(self.value, b)
        if not PY_TYPE_CHECK(b, Integer):
            b = Integer(b)
        return mpz_kronecker(self.value, (<Integer>b).value)

    def radical(self):
        r"""
        Return the product of the prime divisors of self.

        If self is 0, returns 1.

        EXAMPLES:
            sage: Integer(10).radical()
            10
            sage: Integer(20).radical()
            10
            sage: Integer(-20).radical()
            -10
            sage: Integer(0).radical()
            1
            sage: Integer(36).radical()
            6
        """
        if self.is_zero():
            return ONE
        F = self.factor()
        n = one
        for p, e in F:
            n *= p
        return n * F.unit()

    def squarefree_part(self):
        r"""
        Return the square free part of $x$ (=self), i.e., the unique
        integer $z$ that $x = z y^2$, with $y^2$ a perfect square and
        $z$ square-free.

        Use \code{self.radical()} for the product of the primes that
        divide self.

        If self is 0, just returns 0.

        EXAMPLES:
            sage: squarefree_part(100)
            1
            sage: squarefree_part(12)
            3
            sage: squarefree_part(17*37*37)
            17
            sage: squarefree_part(-17*32)
            -34
            sage: squarefree_part(1)
            1
            sage: squarefree_part(-1)
            -1
            sage: squarefree_part(-2)
            -2
            sage: squarefree_part(-4)
            -1
        """
        if self.is_zero():
            return self
        F = self.factor()
        n = one
        for p, e in F:
            if e % 2 != 0:
                n = n * p
        return n * F.unit()

    def next_probable_prime(self):
        """
        Returns the next probable prime after self, as determined by PARI.

        EXAMPLES:
            sage: (-37).next_probable_prime()
            2
            sage: (100).next_probable_prime()
            101
            sage: (2^512).next_probable_prime()
            13407807929942597099574024998205846127479365820592393377723561443721764030073546976801874298166903427690031858186486050853753882811946569946433649006084171
            sage: 0.next_probable_prime()
            2
            sage: 126.next_probable_prime()
            127
            sage: 144168.next_probable_prime()
            144169
        """
        return Integer( (self._pari_()+1).nextprime())

    def next_prime(self, proof=None):
        r"""
        Returns the next prime after self.

        INPUT:
            proof -- bool or None (default: None, see proof.arithmetic or
                            sage.structure.proof)
                        Note that the global Sage default is proof=True

        EXAMPLES:
            sage: Integer(100).next_prime()
            101

        Use Proof = False, which is way faster:
            sage: b = (2^1024).next_prime(proof=False)

            sage: Integer(0).next_prime()
            2
            sage: Integer(1001).next_prime()
            1009
        """
        if proof is None:
            from sage.structure.proof.proof import get_flag
            proof = get_flag(proof, "arithmetic")
        if self < 2:   # negatives are not prime.
            return integer_ring.ZZ(2)
        if self == 2:
            return integer_ring.ZZ(3)
        if not proof:
            return Integer( (self._pari_()+1).nextprime())
        n = self
        if n % 2 == 0:
            n += 1
        else:
            n += 2
        while not n.is_prime():  # pari isprime is provably correct
            n += 2
        return integer_ring.ZZ(n)


    def additive_order(self):
        """
        Return the additive order of self.

        EXAMPLES:
            sage: ZZ(0).additive_order()
            1
            sage: ZZ(1).additive_order()
            +Infinity
        """
        import sage.rings.infinity
        if self.is_zero():
            return one
        else:
            return sage.rings.infinity.infinity

    def multiplicative_order(self):
        r"""
        Return the multiplicative order of self.

        EXAMPLES:
            sage: ZZ(1).multiplicative_order()
            1
            sage: ZZ(-1).multiplicative_order()
            2
            sage: ZZ(0).multiplicative_order()
            +Infinity
            sage: ZZ(2).multiplicative_order()
            +Infinity
        """
        import sage.rings.infinity
        if  mpz_cmp_si(self.value, 1) == 0:
                return one
        elif mpz_cmp_si(self.value, -1) == 0:
                return Integer(2)
        else:
                return sage.rings.infinity.infinity

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

        TESTS:
            sage: n = 10^10000000
            sage: m = n._pari_() ## crash from trac 875
            sage: m % 1234567
            1041334
        """
        return self._pari_c()

    cdef _pari_c(self):

        cdef PariInstance P
        P = sage.libs.pari.gen.pari
        return P.new_gen_from_mpz_t(self.value)

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

    def _sage_input_(self, sib, coerced):
        r"""
        Produce an expression which will reproduce this value when evaluated.

        EXAMPLES:
            sage: sage_input(1, verify=True)
            # Verified
            1
            sage: sage_input(1, preparse=False)
            ZZ(1)
            sage: sage_input(-12435, verify=True)
            # Verified
            -12435
            sage: sage_input(0, verify=True)
            # Verified
            0
            sage: sage_input(-3^70, verify=True)
            # Verified
            -2503155504993241601315571986085849
            sage: sage_input(-37, preparse=False)
            -ZZ(37)
            sage: sage_input(-37 * polygen(ZZ), preparse=False)
            R = ZZ['x']
            x = R.gen()
            -37*x
            sage: from sage.misc.sage_input import SageInputBuilder
            sage: (-314159)._sage_input_(SageInputBuilder(preparse=False), False)
            {unop:- {call: {atomic:ZZ}({atomic:314159})}}
            sage: (314159)._sage_input_(SageInputBuilder(preparse=False), True)
            {atomic:314159}
        """
        if coerced or sib.preparse():
            return sib.int(self)
        else:
            if self < 0:
                return -sib.name('ZZ')(sib.int(-self))
            else:
                return sib.name('ZZ')(sib.int(self))

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
            ValueError: square root of negative integer not defined.
        """
        if mpz_sgn(self.value) < 0:
            raise ValueError, "square root of negative integer not defined."
        cdef Integer x
        x = PY_NEW(Integer)

        _sig_on
        mpz_sqrt(x.value, self.value)
        _sig_off

        return x

    def sqrt_approx(self, prec=None, all=False):
        """
        EXAMPLES:
            sage: 5.sqrt_approx(prec=200)
            2.2360679774997896964091736687312762354406183596115257242709
            sage: 5.sqrt_approx()
            2.23606797749979
            sage: 4.sqrt_approx()
            2
        """
        try:
            return self.sqrt(extend=False,all=all)
        except ValueError:
            pass
        if prec is None:
            prec = max(53, 2*(mpz_sizeinbase(self.value, 2)+2))
        return self.sqrt(prec=prec, all=all)

    def sqrt(self, prec=None, extend=True, all=False):
        """
        The square root function.

        INPUT:
            prec -- integer (default: None): if None, returns an exact
                 square root; otherwise returns a numerical square
                 root if necessary, to the given bits of precision.
            extend -- bool (default: True); if True, return a square
                 root in an extension ring, if necessary. Otherwise,
                 raise a ValueError if the square is not in the base
                 ring.
            all -- bool (default: False); if True, return all square
                 roots of self, instead of just one.

        EXAMPLES:
            sage: Integer(144).sqrt()
            12
            sage: Integer(102).sqrt()
            sqrt(102)

            sage: n = 2
            sage: n.sqrt(all=True)
            [sqrt(2), -sqrt(2)]
            sage: n.sqrt(prec=10)
            1.4
            sage: n.sqrt(prec=100)
            1.4142135623730950488016887242
            sage: n.sqrt(prec=100,all=True)
            [1.4142135623730950488016887242, -1.4142135623730950488016887242]
            sage: n.sqrt(extend=False)
            Traceback (most recent call last):
            ...
            ValueError: square root of 2 not an integer
            sage: Integer(144).sqrt(all=True)
            [12, -12]
            sage: Integer(0).sqrt(all=True)
            [0]
        """
        if mpz_sgn(self.value) == 0:
            return [self] if all else self

        if mpz_sgn(self.value) < 0:
            if not extend:
                raise ValueError, "square root of negative number not an integer"
            if prec:
                from sage.rings.complex_field import ComplexField
                K = ComplexField(prec)
                return K(self).sqrt(all=all)
            from sage.calculus.calculus import sqrt
            return sqrt(self, all=all)


        cdef int non_square
        cdef Integer z = PY_NEW(Integer)
        cdef mpz_t tmp
        _sig_on
        mpz_init(tmp)
        mpz_sqrtrem(z.value, tmp, self.value)
        non_square = mpz_sgn(tmp) != 0
        mpz_clear(tmp)
        _sig_off

        if non_square:
            if not extend:
                raise ValueError, "square root of %s not an integer"%self
            if prec:
                from sage.rings.real_mpfr import RealField
                K = RealField(prec)
                return K(self).sqrt(all=all)
            from sage.calculus.calculus import sqrt
            return sqrt(self, all=all)

        if all:
           return [z, -z]
        return z

    def _xgcd(self, Integer n, bint minimal=0):
        r"""
        Computes extended gcd of self and $n$.

        INPUT:
            self -- integer
            n -- integer
            minimal -- boolean (default false), whether to compute minimal
                       cofactors (see below)

        OUTPUT:
            A triple (g, s, t), where $g$ is the non-negative gcd of self
            and $n$, and $s$ and $t$ are cofactors satisfying the Bezout
            identity
              $$ g = s \cdot \mbox{\rm self} + t \cdot n. $$

        NOTE:
            If minimal is False, then there is no guarantee that the returned
            cofactors will be minimal in any sense; the only guarantee is that
            the Bezout identity will be satisfied (see examples below).

            If minimal is True, the cofactors will satisfy the following
            conditions. If either self or $n$ are zero, the trivial solution
            is returned. If both self and $n$ are nonzero, the function returns
            the unique solution such that $0 \leq s < |n|/g$ (which then must
            also satisfy $0 \leq |t| \leq |\mbox{\rm self}|/g$).

        EXAMPLES:
                sage: Integer(5)._xgcd(7)
                (1, -4, 3)
                sage: 5*-4 + 7*3
                1
                sage: g,s,t = Integer(58526524056)._xgcd(101294172798)
                sage: g
                22544886
                sage: 58526524056 * s + 101294172798 * t
                22544886

            Without minimality guarantees, weird things can happen:
                sage: Integer(3)._xgcd(21)
                (3, -20, 3)
                sage: Integer(3)._xgcd(24)
                (3, -15, 2)
                sage: Integer(3)._xgcd(48)
                (3, -15, 1)

            Weirdness disappears with minimal option:
                sage: Integer(3)._xgcd(21, minimal=True)
                (3, 1, 0)
                sage: Integer(3)._xgcd(24, minimal=True)
                (3, 1, 0)
                sage: Integer(3)._xgcd(48, minimal=True)
                (3, 1, 0)
                sage: Integer(21)._xgcd(3, minimal=True)
                (3, 0, 1)

            Try minimal option with various edge cases:
                sage: Integer(5)._xgcd(0, minimal=True)
                (5, 1, 0)
                sage: Integer(-5)._xgcd(0, minimal=True)
                (5, -1, 0)
                sage: Integer(0)._xgcd(5, minimal=True)
                (5, 0, 1)
                sage: Integer(0)._xgcd(-5, minimal=True)
                (5, 0, -1)
                sage: Integer(0)._xgcd(0, minimal=True)
                (0, 1, 0)

            Exhaustive tests, checking minimality conditions:
                sage: for a in srange(-20, 20):
                ...     for b in srange(-20, 20):
                ...       if a == 0 or b == 0: continue
                ...       g, s, t = a._xgcd(b)
                ...       assert g > 0
                ...       assert a % g == 0 and b % g == 0
                ...       assert a*s + b*t == g
                ...       g, s, t = a._xgcd(b, minimal=True)
                ...       assert g > 0
                ...       assert a % g == 0 and b % g == 0
                ...       assert a*s + b*t == g
                ...       assert s >= 0 and s < abs(b)/g
                ...       assert abs(t) <= abs(a)/g

        AUTHORS:
            -- David Harvey (2007-12-26): added minimality option

        """
        cdef Integer g = PY_NEW(Integer)
        cdef Integer s = PY_NEW(Integer)
        cdef Integer t = PY_NEW(Integer)

        _sig_on
        mpz_gcdext(g.value, s.value, t.value, self.value, n.value)
        _sig_off

        # Note: the GMP documentation for mpz_gcdext (or mpn_gcdext for that
        # matter) makes absolutely no claims about any minimality conditions
        # satisfied by the returned cofactors. They guarantee a non-negative
        # gcd, but that's it. So we have to do some work ourselves.

        if not minimal:
            return g, s, t

        # handle degenerate cases n == 0 and self == 0

        if not mpz_sgn(n.value):
            mpz_set_ui(t.value, 0)
            mpz_abs(g.value, self.value)
            mpz_set_si(s.value, 1 if mpz_sgn(self.value) >= 0 else -1)
            return g, s, t

        if not mpz_sgn(self.value):
            mpz_set_ui(s.value, 0)
            mpz_abs(g.value, n.value)
            mpz_set_si(t.value, 1 if mpz_sgn(n.value) >= 0 else -1)
            return g, s, t

        # both n and self are nonzero, so we need to do a division and
        # make the appropriate adjustment

        cdef mpz_t u1, u2
        mpz_init(u1)
        mpz_init(u2)
        mpz_divexact(u1, n.value, g.value)
        mpz_divexact(u2, self.value, g.value)
        if mpz_sgn(u1) > 0:
            mpz_fdiv_qr(u1, s.value, s.value, u1)
        else:
            mpz_cdiv_qr(u1, s.value, s.value, u1)
        mpz_addmul(t.value, u1, u2)
        mpz_clear(u2)
        mpz_clear(u1)

        return g, s, t


    cdef _lshift(self, long int n):
        """
        Shift self n bits to the left, i.e., quickly multiply by $2^n$.
        """
        cdef Integer x
        x = <Integer> PY_NEW(Integer)

        _sig_on
        if n < 0:
            mpz_fdiv_q_2exp(x.value, self.value, -n)
        else:
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
            sage: 1 << 2.5
            Traceback (most recent call last):
            ...
            TypeError: unsupported operands for <<
        """
        try:
            if not PY_TYPE_CHECK(x, Integer):
                x = Integer(x)
            elif not PY_TYPE_CHECK(y, Integer):
                y = Integer(y)
            return (<Integer>x)._lshift(long(y))
        except TypeError:
            raise TypeError, "unsupported operands for <<"
        if PY_TYPE_CHECK(x, Integer) and isinstance(y, (Integer, int, long)):
            return (<Integer>x)._lshift(long(y))
        return bin_op(x, y, operator.lshift)

    cdef _rshift(Integer self, long int n):
        cdef Integer x
        x = <Integer> PY_NEW(Integer)

        _sig_on
        if n < 0:
            mpz_mul_2exp(x.value, self.value, -n)
        else:
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
            sage: 1 >> 2.5
            Traceback (most recent call last):
            ...
            TypeError: unsupported operands for >>
        """
        try:
            if not PY_TYPE_CHECK(x, Integer):
                x = Integer(x)
            elif not PY_TYPE_CHECK(y, Integer):
                y = Integer(y)
            return (<Integer>x)._rshift(long(y))
        except TypeError:
            raise TypeError, "unsupported operands for >>"

        #if PY_TYPE_CHECK(x, Integer) and isinstance(y, (Integer, int, long)):
        #    return (<Integer>x)._rshift(long(y))
        #return bin_op(x, y, operator.rshift)

    cdef _and(Integer self, Integer other):
        cdef Integer x
        x = PY_NEW(Integer)
        mpz_and(x.value, self.value, other.value)
        return x

    def __and__(x, y):
        """
        Return the bitwise and two integers.

        EXAMPLES:
            sage: n = Integer(6);  m = Integer(2)
            sage: n & m
            2
            sage: n.__and__(m)
            2
        """
        if PY_TYPE_CHECK(x, Integer) and PY_TYPE_CHECK(y, Integer):
            return (<Integer>x)._and(y)
        return bin_op(x, y, operator.and_)


    cdef _or(Integer self, Integer other):
        cdef Integer x
        x = PY_NEW(Integer)
        mpz_ior(x.value, self.value, other.value)
        return x

    def __or__(x, y):
        """
        Return the bitwise or of the integers x and y.

        EXAMPLES:
            sage: n = 8; m = 4
            sage: n.__or__(m)
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
        return one/self    # todo: optimize

    def inverse_of_unit(self):
        """
        Return inverse of self if self is a unit in the integers, i.e.,
        self is -1 or 1.  Otherwise, raise a ZeroDivisionError.

        EXAMPLES:
            sage: (1).inverse_of_unit()
            1
            sage: (-1).inverse_of_unit()
            -1
            sage: 5.inverse_of_unit()
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Inverse does not exist.
            sage: 0.inverse_of_unit()
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Inverse does not exist.
        """
        if mpz_cmp_si(self.value, 1) == 0 or mpz_cmp_si(self.value, -1) == 0:
            return self
        else:
            raise ZeroDivisionError, "Inverse does not exist."

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

        if m == one:
            return zero

        mpz_init(x)

        _sig_on
        r = mpz_invert(x, self.value, m.value)
        _sig_off

        if r == 0:
            raise ZeroDivisionError, "Inverse does not exist."
        ans = PY_NEW(Integer)
        set_mpz(ans,x)
        mpz_clear(x)
        return ans

    def gcd(self, n):
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
        cdef Integer _n = Integer(n)

        mpz_init(g)


        _sig_on
        mpz_gcd(g, self.value, _n.value)
        _sig_off

        g0 = PY_NEW(Integer)
        set_mpz(g0,g)
        mpz_clear(g)
        return g0

    def crt(self, y, m, n):
        """
        Return the unique integer between $0$ and $mn$ that is
        congruent to the integer modulo $m$ and to $y$ modulo $n$.  We
        assume that~$m$ and~$n$ are coprime.

        EXAMPLES:
            sage: n = 17
            sage: m = n.crt(5, 23, 11); m
            247
            sage: m%23
            17
            sage: m%11
            5
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
Py_INCREF(ONE)

def LCM_list(v):
    """
    Return the LCM of a list v of integers.  Element of v is
    converted to a SAGE integer if it isn't one already.

    This function is used, e.g., by rings/arith.py

    INPUT:
        v -- list or tuple

    OUTPUT:
        integer

    EXAMPLES:
        sage: from sage.rings.integer import LCM_list
        sage: w = LCM_list([3,9,30]); w
        90
        sage: type(w)
        <type 'sage.rings.integer.Integer'>

    The inputs are converted to SAGE integers.
        sage: w = LCM_list([int(3), int(9), '30']); w
        90
        sage: type(w)
        <type 'sage.rings.integer.Integer'>
    """
    cdef int i, n
    cdef mpz_t z
    cdef Integer w

    n = len(v)

    if n == 0:
        return one

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


    w = PY_NEW(Integer)
    mpz_set(w.value, z)
    mpz_clear(z)
    return w



def GCD_list(v):
    """
    Return the GCD of a list v of integers.  Element of v is
    converted to a SAGE integer if it isn't one already.

    This function is used, e.g., by rings/arith.py

    INPUT:
        v -- list or tuple

    OUTPUT:
        integer

    EXAMPLES:
        sage: from sage.rings.integer import GCD_list
        sage: w = GCD_list([3,9,30]); w
        3
        sage: type(w)
        <type 'sage.rings.integer.Integer'>

    The inputs are converted to SAGE integers.
        sage: w = GCD_list([int(3), int(9), '30']); w
        3
        sage: type(w)
        <type 'sage.rings.integer.Integer'>
    """
    cdef int i, n
    cdef mpz_t z
    cdef Integer w

    n = len(v)

    if n == 0:
        return one

    try:
        w = v[0]
        mpz_init_set(z, w.value)

        _sig_on
        for i from 1 <= i < n:
            w = v[i]
            mpz_gcd(z, z, w.value)
            if mpz_cmp_si(z, 1) == 0:
                _sig_off
                return one
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
                return one
        _sig_off


    w = PY_NEW(Integer)
    mpz_set(w.value, z)
    mpz_clear(z)
    return w

def make_integer(s):
    """
    Create a SAGE integer from the base-32 Python *string* s.
    This is used in unpickling integers.

    EXAMPLES:
        sage: from sage.rings.integer import make_integer
        sage: make_integer('-29')
        -73
        sage: make_integer(29)
        Traceback (most recent call last):
        ...
        TypeError: expected string or Unicode object, sage.rings.integer.Integer found
    """
    cdef Integer r
    r = PY_NEW(Integer)
    r._reduce_set(s)
    return r

cdef class int_to_Z(Morphism):
    def __init__(self):
        import integer_ring
        import sage.categories.homset
        from sage.structure.parent import Set_PythonType
        Morphism.__init__(self, sage.categories.homset.Hom(Set_PythonType(int), integer_ring.ZZ))
    cpdef Element _call_(self, a):
        cdef Integer r
        r = <Integer>PY_NEW(Integer)
        mpz_set_si(r.value, PyInt_AS_LONG(a))
        return r
    def _repr_type(self):
        return "Native"

cdef class long_to_Z(Morphism):
    """
    EXAMPLES:
        sage: f = ZZ.coerce_map_from(long); f
        Native morphism:
          From: Set of Python objects of type 'long'
          To:   Integer Ring
        sage: f(1rL)
        1
        sage: f(-10000000000000000000001r)
        -10000000000000000000001
    """
    def __init__(self):
        import integer_ring
        import sage.categories.homset
        from sage.structure.parent import Set_PythonType
        Morphism.__init__(self, sage.categories.homset.Hom(Set_PythonType(long), integer_ring.ZZ))
    cpdef Element _call_(self, a):
        cdef Integer r
        r = <Integer>PY_NEW(Integer)
        mpz_set_pylong(r.value, a)
        return r
    def _repr_type(self):
        return "Native"


############### INTEGER CREATION CODE #####################

# We need a couple of internal GMP datatypes.

# This may be potentialy very dangerous as it reaches
# deeply into the internal structure of GMP which may not
# be consistant across future versions of GMP.
# See extensive note in the fast_tp_new() function below.

include "../ext/python_rich_object.pxi"
cdef extern from "gmp.h":
    ctypedef void* mp_ptr #"mp_ptr"

    # We allocate _mp_d directly (mpz_t is typedef of this in GMP)
    ctypedef struct __mpz_struct "__mpz_struct":
        int _mp_alloc
        int _mp_size
        mp_ptr _mp_d

    # sets the three free, alloc, and realloc function pointers to the
    # memory management functions set in GMP. Accepts NULL pointer.
    # Potentially dangerous if changed by calling
    # mp_set_memory_functions again after we initialized this module.
    void mp_get_memory_functions (void *(**alloc) (size_t), void *(**realloc)(void *, size_t, size_t), void (**free) (void *, size_t))

    # GMP's configuration of how many Bits are stuffed into a limb
    cdef int __GMP_BITS_PER_MP_LIMB

# This variable holds the size of any Integer object in bytes.
cdef int sizeof_Integer

# We use a global Integer element to steal all the references
# from.  DO NOT INITIALIZE IT AGAIN and DO NOT REFERENCE IT!
cdef Integer global_dummy_Integer
global_dummy_Integer = Integer()

# Accessing the .value attribute of an Integer object causes Pyrex to
# refcount it. This is problematic, because that causes overhead and
# more importantly an infinite loop in the destructor. If you refcount
# in the destructor and the refcount reaches zero (which is true
# everytime) the destructor is called.
#
# To avoid this we calculate the byte offset of the value member and
# remember it in this variable.
#
# Eventually this may be rendered obsolete by a change in SageX allowing
# non-reference counted extension types.
cdef long mpz_t_offset
mpz_t_offset_python = None


# stores the GMP alloc function
cdef void * (* mpz_alloc)(size_t)

# stores the GMP free function
cdef void (* mpz_free)(void *, size_t)

# A global  pool for performance when integers are rapidly created and destroyed.
# It operates on the following principles:
#
# - The pool starts out empty.
# - When an new integer is needed, one from the pool is returned
#   if available, otherwise a new Integer object is created
# - When an integer is collected, it will add it to the pool
#   if there is room, otherwise it will be deallocated.


cdef int integer_pool_size = 100

cdef PyObject** integer_pool
cdef int integer_pool_count = 0

# used for profiling the pool
cdef int total_alloc = 0
cdef int use_pool = 0

# The signature of tp_new is
# PyObject* tp_new(RichPyTypeObject *t, PyObject *a, PyObject *k).
# However we only use t in this implementation.
#
# t in this case is the Integer TypeObject.

cdef PyObject* fast_tp_new(RichPyTypeObject *t, PyObject *a, PyObject *k):

    global integer_pool, integer_pool_count, total_alloc, use_pool

    cdef RichPyObject* new

    # for profiling pool usage
    # total_alloc += 1

    # If there is a ready integer in the pool, we will
    # decrement the counter and return that.

    if integer_pool_count > 0:

        # for profiling pool usage
        # use_pool += 1

        integer_pool_count -= 1
        new = <RichPyObject *> integer_pool[integer_pool_count]

    # Otherwise, we have to create one.

    else:

        # allocate enough room for the Integer, sizeof_Integer is
        # sizeof(Integer). The use of PyObject_MALLOC directly
        # assumes that Integers are not garbage collected, i.e.
        # they do not pocess references to other Python
        # objects (Aas indicated by the Py_TPFLAGS_HAVE_GC flag).
        # See below for a more detailed description.
        new = PyObject_MALLOC( sizeof_Integer )

        # Now set every member as set in z, the global dummy Integer
        # created before this tp_new started to operate.

        memcpy(new, (<void*>global_dummy_Integer), sizeof_Integer )

        # This line is only needed if Python is compiled in debugging
        # mode './configure --with-pydebug'. If that is the case a Python
        # object has a bunch of debugging fields which are initialized
        # with this macro. For speed reasons, we don't call it if Python
        # is not compiled in debug mode. So uncomment the following line
        # if you are debugging Python.

        #PyObject_INIT(new, (<RichPyObject*>global_dummy_Integer).ob_type)

        # We take the address 'new' and move mpz_t_offset bytes (chars)
        # to the address of 'value'. We treat that address as a pointer
        # to a mpz_t struct and allocate memory for the _mp_d element of
        # that struct. We allocate one limb.
        #
        # What is done here is potentialy very dangerous as it reaches
        # deeply into the internal structure of GMP. Consequently things
        # may break if a new release of GMP changes some internals. To
        # emphazise this, this is what the GMP manual has to say about
        # the documentation for the struct we are using:
        #
        #  "This chapter is provided only for informational purposes and the
        #  various internals described here may change in future GMP releases.
        #  Applications expecting to be compatible with future releases should use
        #  only the documented interfaces described in previous chapters."
        #
        # If this line is used SAGE is not such an application.
        #
        # The clean version of the following line is:
        #
        #  mpz_init( <mpz_t>(<char *>new + mpz_t_offset) )
        #
        # We save time both by avoiding an extra function call and
        # because the rest of the mpz struct was already initalized
        # fully using the memcpy above.

        (<__mpz_struct *>( <char *>new + mpz_t_offset) )._mp_d = <mp_ptr>mpz_alloc(__GMP_BITS_PER_MP_LIMB >> 3)

    # The global_dummy_Integer may have a reference count larger than
    # one, but it is expected that newly created objects have a
    # reference count of one. This is potentially unneeded if
    # everybody plays nice, because the gobal_dummy_Integer has only
    # one reference in that case.

    # Objects from the pool have reference count zero, so this
    # needs to be set in this case.

    new.ob_refcnt = 1

    return new

cdef void fast_tp_dealloc(PyObject* o):

    # If there is room in the pool for a used integer object,
    # then put it in rather than deallocating it.

    global integer_pool, integer_pool_count

    if integer_pool_count < integer_pool_size:

        # Here we free any extra memory used by the mpz_t by
        # setting it to a single limb.
        if (<__mpz_struct *>( <char *>o + mpz_t_offset))._mp_alloc > 1:
            _mpz_realloc(<mpz_t *>( <char *>o + mpz_t_offset), 1)

        # It's cheap to zero out an integer, so do it here.
        (<__mpz_struct *>( <char *>o + mpz_t_offset))._mp_size = 0

        # And add it to the pool.
        integer_pool[integer_pool_count] = o
        integer_pool_count += 1
        return

    # Again, we move to the mpz_t and clear it. See above, why this is evil.
    # The clean version of this line would be:
    #   mpz_clear(<mpz_t>(<char *>o + mpz_t_offset))

    mpz_free((<__mpz_struct *>( <char *>o + mpz_t_offset) )._mp_d, 0)

    # Free the object. This assumes that Py_TPFLAGS_HAVE_GC is not
    # set. If it was set another free function would need to be
    # called.

    PyObject_FREE(o)


hook_fast_tp_functions()
from sage.misc.allocator cimport hook_tp_functions
cdef hook_fast_tp_functions():
    """
    Initialize the fast integer creation functions.
    """
    global global_dummy_Integer, mpz_t_offset, sizeof_Integer, integer_pool

    integer_pool = <PyObject**>sage_malloc(integer_pool_size * sizeof(PyObject*))

    cdef RichPyObject* o
    o = <RichPyObject*>global_dummy_Integer

    # calculate the offset of the GMP mpz_t to avoid casting to/from
    # an Integer which includes reference counting. Reference counting
    # is bad in constructors and destructors as it potentially calls
    # the destructor.
    # Eventually this may be rendered obsolete by a change in SageX allowing
    # non-reference counted extension types.
    mpz_t_offset = <char *>(&global_dummy_Integer.value) - <char *>o
    global mpz_t_offset_python
    mpz_t_offset_python = mpz_t_offset

    # store how much memory needs to be allocated for an Integer.
    sizeof_Integer = o.ob_type.tp_basicsize

    # get the functions to do memory management for the GMP elements
    # WARNING: if the memory management functions are changed after
    # this initialisation, we are/you are doomed.

    mp_get_memory_functions(&mpz_alloc, NULL, &mpz_free)

    # Finally replace the functions called when an Integer needs
    # to be constructed/destructed.
    hook_tp_functions(global_dummy_Integer, NULL, <newfunc>&fast_tp_new, NULL, &fast_tp_dealloc, False)

cdef integer(x):
    if PY_TYPE_CHECK(x, Integer):
        return x
    return Integer(x)


def free_integer_pool():
    cdef int i
    cdef PyObject *o

    global integer_pool_count, integer_pool_size

    for i from 0 <= i < integer_pool_count:
        o = integer_pool[i]
        mpz_clear(<__mpz_struct *>(<char*>o + mpz_t_offset))
        # Free the object. This assumes that Py_TPFLAGS_HAVE_GC is not
        # set. If it was set another free function would need to be
        # called.
        PyObject_FREE(o)

    integer_pool_size = 0
    integer_pool_count = 0
    sage_free(integer_pool)
