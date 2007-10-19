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
include "../../ext/cdefs.pxi"
include 'misc.pxi'
include 'decl.pxi'

from sage.rings.integer import Integer
from sage.rings.integer_ring import IntegerRing
from sage.rings.integer cimport Integer
from sage.rings.integer_ring cimport IntegerRing_class

ZZ_sage = IntegerRing()

cdef make_ZZ(ZZ_c* x):
    cdef ntl_ZZ y
    y = ntl_ZZ()
    y.x = x[0]
    ZZ_delete(x)
    _sig_off
    return y


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
    def __init__(self, v=None):
        r"""
        Initializes and NTL integer.

        EXAMPLES:
            sage: ntl.ZZ(12r)
            12
            sage: ntl.ZZ(Integer(95413094))
            95413094
            sage: ntl.ZZ(long(223895239852389582983))
            223895239852389582983
            sage: ntl.ZZ('-1')
            -1
            sage: ntl.ZZ('1L')
            1
            sage: ntl.ZZ('-1r')
            -1

        AUTHOR: Joel B. Mohler (2007-06-14)
        """
        if PY_TYPE_CHECK(v, ntl_ZZ):
            self.x = (<ntl_ZZ>v).x
        elif PyInt_Check(v):
            ZZ_conv_int(self.x, v)
        elif PyLong_Check(v):
            ZZ_set_pylong(self.x, v)
        elif PY_TYPE_CHECK(v, Integer):
            self.set_from_sage_int(v)
        elif v is not None:
            v = str(v)
            if len(v) == 0:
                v = '0'
            if not ((v[0].isdigit() or v[0] == '-') and \
                    (v[1:-1].isdigit() or (len(v) <= 2)) and \
                    (v[-1].isdigit() or (v[-1].lower() in ['l','r']))):
               raise ValueError, "invalid integer: %s"%v
            _sig_on
            ZZ_from_str(&self.x, v)
            _sig_off

    def __new__(self, v=None):
        ZZ_construct(&self.x)

    def __dealloc__(self):
        ZZ_destruct(&self.x)

    def __repr__(self):
        """
        Return the string representation of self.

        EXAMPLES:
            sage: ntl.ZZ(5).__repr__()
            '5'
        """
        return ZZ_to_PyString(&self.x)

    def __reduce__(self):
        """
        sage: from sage.libs.ntl.ntl_ZZ import ntl_ZZ
        sage: a = ntl_ZZ(-7)
        sage: loads(dumps(a))
        -7
        """
        return unpickle_class_value, (ntl_ZZ, self.get_as_sage_int())

    def __cmp__(self, other):
        """
        Compare self to other.

        EXAMPLES:
            sage: f = ntl.ZZ(1)
            sage: g = ntl.ZZ(2)
            sage: h = ntl.ZZ(2)
            sage: w = ntl.ZZ(7)
            sage: h == g
            True
            sage: f == g
            False
            sage: h > w ## indirect doctest
            False
            sage: h < w
            True
        """
        if (type(self) != type(other)):
            return cmp(type(self), type(other))
        diff = self.__sub__(other)
        if ZZ_IsZero( (<ntl_ZZ>diff).x ):
            return 0
        elif ZZ_sign( (<ntl_ZZ>diff).x ) == 1:
            return 1
        else:
            return -1

    def __mul__(self, other):
        """
            sage: n=ntl.ZZ(2983)*ntl.ZZ(2)
            sage: n
            5966
        """
        cdef ntl_ZZ r = PY_NEW(ntl_ZZ)
        if not PY_TYPE_CHECK(self, ntl_ZZ):
            self = ntl_ZZ(self)
        if not PY_TYPE_CHECK(other, ntl_ZZ):
            other = ntl_ZZ(other)
        mul_ZZ(r.x, (<ntl_ZZ>self).x, (<ntl_ZZ>other).x)
        return r

    def __sub__(self, other):
        """
            sage: n=ntl.ZZ(2983)-ntl.ZZ(2)
            sage: n
            2981
            sage: ntl.ZZ(2983)-2
            2981
        """
        cdef ntl_ZZ r = PY_NEW(ntl_ZZ)
        if not PY_TYPE_CHECK(self, ntl_ZZ):
            self = ntl_ZZ(self)
        if not PY_TYPE_CHECK(other, ntl_ZZ):
            other = ntl_ZZ(other)
        sub_ZZ(r.x, (<ntl_ZZ>self).x, (<ntl_ZZ>other).x)
        return r

    def __add__(self, other):
        """
            sage: n=ntl.ZZ(2983)+ntl.ZZ(2)
            sage: n
            2985
            sage: ntl.ZZ(23)+2
            25
        """
        cdef ntl_ZZ r = PY_NEW(ntl_ZZ)
        if not PY_TYPE_CHECK(self, ntl_ZZ):
            self = ntl_ZZ(self)
        if not PY_TYPE_CHECK(other, ntl_ZZ):
            other = ntl_ZZ(other)
        add_ZZ(r.x, (<ntl_ZZ>self).x, (<ntl_ZZ>other).x)
        return r

    def __neg__(ntl_ZZ self):
        """
            sage: x = ntl.ZZ(38)
            sage: -x
            -38
            sage: x.__neg__()
            -38
        """
        cdef ntl_ZZ r = PY_NEW(ntl_ZZ)
        ZZ_negate(r.x, self.x)
        return r

    def __pow__(ntl_ZZ self, long e, ignored):
        """
            sage: ntl.ZZ(23)^50
            122008981252869411022491112993141891091036959856659100591281395343249
        """
        cdef ntl_ZZ r = ntl_ZZ()
        power_ZZ(r.x, self.x, e)
        return r

    def __int__(self):
        """
        Return self as an int.

        EXAMPLES:
            sage: ntl.ZZ(22).__int__()
            22
            sage: type(ntl.ZZ(22).__int__())
            <type 'int'>

            sage: ntl.ZZ(10^30).__int__()
            1000000000000000000000000000000L
            sage: type(ntl.ZZ(10^30).__int__())
            <type 'long'>
        """
        return int(self.get_as_sage_int())

    cdef int get_as_int(ntl_ZZ self):
        r"""
        Returns value as C int.
        Return value is only valid if the result fits into an int.

        AUTHOR: David Harvey (2006-08-05)
        """
        return ZZ_to_int(&self.x)

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
        return (<IntegerRing_class>ZZ_sage)._coerce_ZZ(&self.x)

    def _integer_(self):
        r"""
        Gets the value as a sage int.

        sage: n=ntl.ZZ(2983)
        sage: type(n._integer_())
        <type 'sage.rings.integer.Integer'>

        Alias for get_as_sage_int
        """
        return (<IntegerRing_class>ZZ_sage)._coerce_ZZ(&self.x)

    cdef void set_from_int(ntl_ZZ self, int value):
        r"""
        Sets the value from a C int.

        AUTHOR: David Harvey (2006-08-05)
        """
        ZZ_set_from_int(&self.x, value)

    def set_from_sage_int(self, Integer value):
        r"""
        Sets the value from a sage int.

        EXAMPLES:
            sage: n=ntl.ZZ(2983)
            sage: n
            2983
            sage: n.set_from_sage_int(1234)
            sage: n
            1234

        AUTHOR: Joel B. Mohler
        """
        _sig_on
        value._to_ZZ(&self.x)
        _sig_off

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

def unpickle_class_value(cls, x):
    """
    Here for unpickling.

    EXAMPLES:
        sage: sage.libs.ntl.ntl_ZZ.unpickle_class_value(ntl.ZZ, 3)
        3
        sage: type(sage.libs.ntl.ntl_ZZ.unpickle_class_value(ntl.ZZ, 3))
        <type 'sage.libs.ntl.ntl_ZZ.ntl_ZZ'>
    """
    return cls(x)

def unpickle_class_args(cls, x):
    """
    Here for unpickling.

    EXAMPLES:
        sage: sage.libs.ntl.ntl_ZZ.unpickle_class_args(ntl.ZZ, [3])
        3
        sage: type(sage.libs.ntl.ntl_ZZ.unpickle_class_args(ntl.ZZ, [3]))
        <type 'sage.libs.ntl.ntl_ZZ.ntl_ZZ'>
    """
    return cls(*x)

# Random-number generation
def ntl_setSeed(x=None):
    """
    Seed the NTL random number generator.

    EXAMPLE:
        sage: ntl.ntl_setSeed(10)
        sage: ntl.ZZ_random(1000)
        776
    """
    cdef ntl_ZZ seed = ntl_ZZ(1)
    if x is None:
        from random import randint
        seed = ntl_ZZ(str(randint(0,int(2)**64)))
    else:
        seed = ntl_ZZ(str(x))
    _sig_on
    setSeed(&seed.x)
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

    if not PY_TYPE_CHECK(q, ntl_ZZ):
        q = ntl_ZZ(str(q))
    w = q
    _sig_on
    return  make_ZZ(ZZ_randomBnd(&w.x))

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
