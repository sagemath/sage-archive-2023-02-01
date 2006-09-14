r"""
Elements of $\Z/n\Z$

An element of the integers modulo $n$.
"""

include "../ext/interrupt.pxi"  # ctrl-c interrupt block support
#cdef extern from "../ext/mpz_pylong.h":
#    cdef long mpz_pythonhash(mpz_t src)

import operator

import integer_mod_ring
import arith
import rational
from sage.libs.all import pari, PariError
import integer_ring
import integer
import commutative_ring_element
import sage.interfaces.all

import sage.ext.element
cimport sage.ext.element
cimport sage.ext.integer

from sage.rings.coerce import cmp as coerce_cmp

def Mod(n, m):
    """
    Return the equivalence class of n modulo m
    as an element of $\Z/m\Z$.

    EXAMPLES:
        sage: x = Mod(12345678, 32098203845329048)
        sage: x
        12345678
        sage: x**100
        1017322209155072

    You can also use the lowercase version:
        sage: mod(12,5)
        2
    """


    return IntegerMod(integer_mod_ring.IntegerModRing(m),n)

mod = Mod


def IntegerMod(parent, value):
    cdef sage.ext.integer.Integer modulus
    modulus = parent.order()
    if mpz_cmp_si(modulus.value, INTEGER_MOD_INT32_LIMIT) < 0:
        return IntegerMod_int(parent, value)
    elif mpz_cmp_si(modulus.value, INTEGER_MOD_INT64_LIMIT) < 0:
        return IntegerMod_int64(parent, value)
    else:
        return IntegerMod_gmp(parent, value)


def is_IntegerMod(x):
    return isinstance(x, IntegerMod_abstract)


def makeNativeIntStruct(sage.ext.integer.Integer z):
    return NativeIntStruct(z)


cdef class NativeIntStruct:

    def __init__(NativeIntStruct self, sage.ext.integer.Integer z):
        self.sageInteger = z
        if mpz_cmp_si(z.value, INTEGER_MOD_INT64_LIMIT) < 0:
            self.int64 = mpz_get_si(z.value)
            if self.int64 < INTEGER_MOD_INT32_LIMIT:
                self.int32 = self.int64

    def __reduce__(NativeIntStruct self):
        return sage.rings.integer_mod_pyx.makeNativeIntStruct, (self.sageInteger, )



cdef class IntegerMod_abstract(sage.ext.element.CommutativeRingElement):

    def __init__(self, parent, value, empty=False):
        """
        EXAMPLES:
            sage: a = Mod(10,30^10); a
            10
            sage: loads(a.dumps()) == a
            True
        """
        if self.__class__ is IntegerMod:
            raise NotImplementedError, "Can't instantiate abstract IntegerMod"
        commutative_ring_element.CommutativeRingElement.__init__(self, parent)
        self.__modulus = parent._pyx_order

    def __reduce__(IntegerMod_abstract self):
        """
        EXAMPLES:
            sage: a = Mod(4,5); a
            4
            sage: loads(a.dumps()) == a
            True
            sage: a = Mod(-1,5^30)^25;
            sage: loads(a.dumps()) == a
            True
        """
        return sage.rings.integer_mod_pyx.mod, (self.lift(), self.modulus())





    #################################################################
    # Interfaces
    #################################################################
    def _pari_init_(self):
        return 'Mod(%s,%s)'%(str(self), self.__modulus.sageInteger)

    def pari(self):
        return pari(self._pari_init_()) # TODO: is this called implicitly anywhere?

    def _gap_init_(self):
        r"""
        Return string representation of corresponding GAP object.

        This can be slow since non-prime GAP finite field elements are
        represented as powers of a generator for the multiplicative
        group, so the discrete log problem must be solved.

        \note{This function will create a meaningless GAP object if the
        modulus is not a power of a prime.  Also, the modulus must
        be $\leq 65536$.}

        EXAMPLES:
            sage: a = Mod(2,19)
            sage: gap(a)
            Z(19)
            sage: a._gap_(gap)
            Z(19)
            sage: gap(a).Int()
            2
            sage: b = Mod(0,25)
            sage: gap(b)
            0*Z(5)
        """
        R = self.parent()
        m = self.__modulus.sageInteger

        if m > 65536:
            raise ValueError, "order must be at most 65536."

        if self == 0:
            return '0*Z(%s)'%m

        # I couldn't find a guarentee in the GAP docs that the
        # root of unity they use must be the smallest.   This
        # was *not* the case in MAGMA once, so who knows, especially
        # when the order of the ring is not prime.  So we make
        # no such dangerous assumptions (for now).

        # Find the root of unity used by Gap.
        from sage.interfaces.all import gap        # here to reduce dependencies
        g = int(gap.eval('Int(Z(%s))'%m))
        n = R(g).log(self)
        return 'Z(%s)^%s'%(m, n)

    def _magma_init_(self):
        """
        Coercion to Magma.

        EXAMPLES:
            sage: a = Integers(15)(4)
            sage: b = magma(a)                # optional
            sage: b.Type()                    # optional
            RngIntResElt
            sage: b^2                         # optional
            1
        """
        return '%s!%s'%(self.parent()._magma_init_(), self)

    def log(self, a):
        """
        Return $x$ such that $b^x = a$, where $b$ is self.

        INPUT:
            self, a are units in the integers modulo $N$.

        OUTPUT:
            Integer $x$ such that $a^x = b$, if it exists.
            Raises a ValueError exception if no such $x$ exists.

        EXAMPLES:
            sage: R = Integers(500)
            sage: b = R(17); a = b^19
            sage: b.log(a)
            19

        AUTHOR: David Joyner and William Stein (2005-11)
        """
        n = self._parent.unit_group_order()
        return arith.discrete_log_generic(self, a, n) # TODO update this function


    def modulus(IntegerMod_abstract self):
        return self.__modulus.sageInteger


    def is_square(self):
        return bool(self.pari().issquare()) # TODO implement directly

    def charpoly(self):
        """
        Returns the characteristic polynomial of this element.

        EXAMPLES:
            sage: k = GF(3)
            sage: a = k.gen()
            sage: a.charpoly()
            x + 2
            sage: a.charpoly()(a)
            0

        AUTHOR:
         -- Craig Citro
        """
        import polynomial_ring
        R = polynomial_ring.PolynomialRing(self._parent)
        return R([-self,1])

    def norm(self):
        """
        Returns the norm of this element, which is itself. (This
        is here for compatibility with higher order finite fields.)

        EXAMPLES:
            sage: k = GF(691)
            sage: a = k(389)
            sage: a.norm()
            389

        AUTHOR:
         -- Craig Citro
        """
        return self

    def trace(self):
        """
        Returns the trace of this element, which is itself. (This
        is here for compatibility with higher order finite fields.)

        EXAMPLES:
            sage: k = GF(691)
            sage: a = k(389)
            sage: a.trace()
            389

        AUTHOR:
         -- Craig Citro
        """
        return self

    def sqrt(self):
        """
        Same as self.square_root().
        """
        return self.square_root()

    def square_root(self):
        """
        EXAMPLES:
            sage: mod(-1, 17).square_root()
            4
            sage: mod(5, 389).square_root()
            86
        """
        try:
            return self.parent()(self.pari().sqrt())  # TODO: implement directly
        except PariError:
            raise ValueError, "self must be a square."

    def rational_reconstruction(self):
        """
        EXAMPLES:
            sage: R = IntegerModRing(97)
            sage: a = R(2) / R(3)
            sage: a
            33
            sage: a.rational_reconstruction()
            2/3
        """
        return self.lift().rational_reconstruction(self.modulus())

    def crt(IntegerMod_abstract self, IntegerMod_abstract other):
        """
        Use the Chinese Remainder Theorem to find an element of the
        integers modulo the product of the moduli that reduces to self
        and to other.  The modulus of other must be coprime to the
        modulus of self.
        EXAMPLES:
            sage: a = mod(3,5)
            sage: b = mod(2,7)
            sage: a.crt(b)
            23

            sage: a = mod(37,10^8)
            sage: b = mod(9,3^8)
            sage: a.crt(b)
            125900000037

        AUTHOR:
            -- Robert Bradshaw
        """
        cdef int_fast64_t new_modulus
        if not isinstance(self, IntegerMod_gmp) and not isinstance(other, IntegerMod_gmp):

            if other.__modulus.int64 == 1: return self
            new_modulus = self.__modulus.int64 * other.__modulus.int64
            if new_modulus < INTEGER_MOD_INT32_LIMIT:
                return self.__crt(other)

            elif new_modulus < INTEGER_MOD_INT64_LIMIT:
                if not isinstance(self, IntegerMod_int64):
                    self = IntegerMod_int64(self._parent, self.lift())
                if not isinstance(other, IntegerMod_int64):
                    other = IntegerMod_int64(other._parent, other.lift())
                return self.__crt(other)

        if not isinstance(self, IntegerMod_gmp):
            self = IntegerMod_gmp(self._parent, self.lift())

        if not isinstance(other, IntegerMod_gmp):
            other = IntegerMod_gmp(other._parent, other.lift())

        return self.__crt(other)

    def __floordiv__(self, right):
        """
        EXAMPLES:
            sage: R = Integers(25)
            sage: R(15)//5
            3
        """

        try:
            return self / right
        except ZeroDivisionError:
            return IntegerMod(self._parent, self.lift() / sage.rings.integer_ring.Z(right) )

    def order(self):
        """
        Returns the additive order of self.
        """
        return integer.Integer(self.__modulus.sageInteger.__FLOORDIV__(self.lift().gcd(n)))

    def multiplicative_order(self):
        """
        Returns the additive order of self.
        """
        if not self.is_unit():
            raise ArithmeticError, "self must be a unit"
        return integer.Integer(self.pari().order())  # pari's "order" is by default multiplicative  # TODO: implement this directly

    def _repr_(self):
        return str(self.lift())

    def _latex_(self):
        return str(self)

    def _integer_(self):
        return self.lift()

    def _rational_(self):
        return rational.Rational(self.lift())




######################################################################
#      class IntegerMod_gmp
######################################################################

cdef class IntegerMod_gmp(IntegerMod_abstract):
    """
    Elements of $\Z/n\Z$ for n not small enough to be operated on in word size
    AUTHORS:
        -- Robert Bradshaw (2006-08-24)
    """

    def __init__(IntegerMod_gmp self, parent, value, empty=False):
        mpz_init(self.value)
        IntegerMod_abstract.__init__(self, parent, value)
        if empty:
            return
        cdef sage.ext.integer.Integer z
        if isinstance(value, sage.ext.integer.Integer):
            z = value
        elif isinstance(value, rational.Rational):
            z = value % self.__modulus.sageInteger
        else:
            z = sage.rings.integer_ring.Z(value)
        self.set_from_mpz(z.value)


    def __dealloc__(self):
        mpz_clear(self.value)


    cdef void set_from_mpz(IntegerMod_gmp self, mpz_t value):
        cdef sage.ext.integer.Integer modulus
        modulus = self.__modulus.sageInteger
        if mpz_sgn(value) == -1 or mpz_cmp(value, modulus.value) >= 0:
            mpz_mod(self.value, value, modulus.value)
        else:
            mpz_set(self.value, value)

    cdef mpz_t* get_value(IntegerMod_gmp self):
        return &self.value


    def __cmp__(IntegerMod_gmp self, right):
        if not isinstance(right, IntegerMod_gmp):
            try:
                return coerce_cmp(self, right)
            except TypeError:
                return -1
        return self.cmp(right)

    def cmp(IntegerMod_gmp self, IntegerMod_gmp right):
        """
        EXAMPLES:
            sage: mod(5,13^20) == mod(5,13^20)
            True
            sage: mod(5,13^20) == mod(-5,13^20)
            False
            sage: mod(5,13^20) == mod(-5,13)
            False
        """
        if right._parent is not self._parent:
            return -1
        cdef int i
        i = mpz_cmp(self.value, right.value)
        if i < 0:
            return -1
        elif i == 0:
            return 0
        else:
            return 1

    def __richcmp__(self, right, int op):
        cdef int n
        if not isinstance(right, IntegerMod_gmp):
            try:
                n = coerce_cmp(self, right)
            except TypeError:
                n = -1
        else:
            n = self.cmp(right)
        return self._rich_to_bool(op, n)

    def is_one(IntegerMod_gmp self):
        """
        Returns \\code{True} if this is $1$, otherwise \\code{False}.

        EXAMPLES:
            sage: mod(6,5).is_one()
            True
            sage: mod(0,5).is_one()
            False
        """
        return bool(mpz_cmp_si(self.value, 1) == 0)

    def is_zero(IntegerMod_gmp self):
        """
        Returns \\code{True} if this is $0$, otherwise \\code{False}.

        EXAMPLES:
            sage: mod(13,5).is_zero()
            False
            sage: mod(25,5).is_zero()
            True
        """
        return bool(mpz_cmp_si(self.value, 0) == 0)

    def is_unit(self):
        return bool(self.lift().gcd(self.modulus()) == 1)

    def __crt(IntegerMod_gmp self, IntegerMod_gmp other):
        cdef IntegerMod_gmp lift, x
        cdef sage.ext.integer.Integer modulus, other_modulus

        modulus = self.__modulus.sageInteger
        other_modulus = other.__modulus.sageInteger
        lift = IntegerMod_gmp(integer_mod_ring.IntegerModRing(modulus*other_modulus, check_prime=False), None, empty=True)
        try:
            if mpz_cmp(self.value, other.value) > 0:
                x = (other - IntegerMod_gmp(other._parent, self.lift())) / IntegerMod_gmp(other._parent, modulus)
                mpz_mul(lift.value, x.value, modulus.value)
                mpz_add(lift.value, lift.value, self.value)
            else:
                x = (self - IntegerMod_gmp(self._parent, other.lift())) / IntegerMod_gmp(self._parent, other_modulus)
                mpz_mul(lift.value, x.value, other_modulus.value)
                mpz_add(lift.value, lift.value, other.value)
            return lift
        except ZeroDivisionError:
            raise ZeroDivisionError, "moduli must be coprime"


    def copy(IntegerMod_gmp self):
        cdef IntegerMod_gmp copy
        copy = IntegerMod_gmp(self._parent, None, empty=True)
        mpz_set(copy.value, self.value)
        return copy


    def _add_(IntegerMod_gmp self, IntegerMod_gmp right):
        """
        EXAMPLES:
            sage: R = Integers(10)
            sage: R(7) + R(8)
            5
        """
        cdef sage.ext.integer.Integer modulus
        modulus = self.__modulus.sageInteger
        cdef IntegerMod_gmp x
        x = IntegerMod_gmp(self._parent, None, empty=True)
        mpz_add(x.value, self.value, right.value)
        if mpz_cmp(x.value, modulus.value)  >= 0:
            mpz_sub(x.value, x.value, modulus.value)
        return x;

    def _sub_(IntegerMod_gmp self, IntegerMod_gmp right):
        """
        EXAMPLES:
            sage: R = Integers(10)
            sage: R(7) - R(8)
            9
        """
        cdef IntegerMod_gmp x
        x = IntegerMod_gmp(self._parent, None, empty=True)
        mpz_sub(x.value, self.value, right.value)
        if mpz_sgn(x.value) == -1:
            mpz_add(x.value, x.value, self.__modulus.sageInteger.value)
        return x;

    def _mul_(IntegerMod_gmp self, IntegerMod_gmp right):
        """
        EXAMPLES:
            sage: R = Integers(10)
            sage: R(7) * R(8)
            6
        """
        cdef sage.ext.integer.Integer modulus
        modulus = self.__modulus.sageInteger
        cdef IntegerMod_gmp x
        x = IntegerMod_gmp(self._parent, None, empty=True)
        mpz_mul(x.value, self.value, right.value)
        mpz_fdiv_r(x.value, x.value, modulus.value)
        return x;

    def _div_(IntegerMod_gmp self, IntegerMod_gmp right):
        return self._mul_(~right)

    def __int__(self):
        return int(self.lift())

    def __long__(self):
        return long(self.lift())

    def __mod__(self, right):
        right = int(right)
        if self.modulus() % right != 0:
            raise ZeroDivisionError, "Error - reduction modulo right not defined."
        return IntegerMod(integer_mod_ring.IntegerModRing(right, self))

    def __pow__(IntegerMod_gmp self, right, m): # NOTE: m ignored, always use modulus of parent ring
        """
        EXAMPLES:
            sage: R = Integers(10)
            sage: R(2)^10
            4
            sage: R = Integers(389)
            sage: R(7)^388
            1
        """
        cdef sage.ext.integer.Integer exp
        exp = sage.rings.integer_ring.Z(right)
        cdef IntegerMod_gmp x
        x = IntegerMod_gmp(self._parent, None, empty=True)
        _sig_on
        mpz_powm(x.value, self.value, exp.value, self.__modulus.sageInteger.value)
        _sig_off
        return x

    def __neg__(IntegerMod_gmp self):
        """
        EXAMPLES:
            sage: -mod(7,10)
            3
        """
        cdef IntegerMod_gmp x
        x = IntegerMod_gmp(self._parent, None, empty=True)
        mpz_sub(x.value, self.__modulus.sageInteger.value, self.value)
        return x

    def __invert__(IntegerMod_gmp self):
        """
        EXAMPLES:
            sage: ~mod(7,100)
            43
        """
        cdef IntegerMod_gmp x
        x = IntegerMod_gmp(self._parent, None, empty=True)
        _sig_on
        if (mpz_invert(x.value, self.value, self.__modulus.sageInteger.value)):
            _sig_off
            return x
        else:
            _sig_off
            raise ZeroDivisionError, "Inverse does not exist."

    def lift(IntegerMod_gmp self):
        cdef sage.ext.integer.Integer z
        z = sage.ext.integer.Integer()
        z.set_from_mpz(self.value)
        return z

    def __float__(self):
        return float(self.lift())

    def __hash__(self):
#        return mpz_pythonhash(self.value)
        return hash(self.lift())



######################################################################
#      class IntegerMod_int
######################################################################

cdef class IntegerMod_int(IntegerMod_abstract):
    """
    Elements of $\Z/n\Z$ for n small enough to be operated on in 32 bits
    AUTHORS:
        -- Robert Bradshaw (2006-08-24)
    """

    def __init__(self, parent, value, empty=False):
        """
        EXAMPLES:
            sage: a = Mod(10,30); a
            10
            sage: loads(a.dumps()) == a
            True
        """
        IntegerMod_abstract.__init__(self, parent, value)
        if empty:
            return
        cdef sage.ext.integer.Integer z
        if isinstance(value, sage.ext.integer.Integer):
            z = value
        elif isinstance(value, rational.Rational):
            z = value % self.__modulus.sageInteger
        else:
            z = sage.rings.integer_ring.Z(value)
        self.set_from_mpz(z.value)


    cdef void set_from_mpz(IntegerMod_int self, mpz_t value):
        if mpz_sgn(value) == -1 or mpz_cmp_si(value, self.__modulus.int32) >= 0:
            self.ivalue = mpz_fdiv_ui(value, self.__modulus.int32)
        else:
            self.ivalue = mpz_get_si(value)

    cdef void set_from_int(IntegerMod_int self, int_fast32_t ivalue):
        if ivalue < 0:
            self.ivalue = self.__modulus.int32 + (ivalue % self.__modulus.int32)
        elif ivalue >= self.__modulus.int32:
            self.ivalue = ivalue % self.__modulus.int32
        else:
            self.ivalue = ivalue

    cdef int_fast32_t get_int_value(IntegerMod_int self):
        return self.ivalue


    def __cmp__(IntegerMod_int self, right):
        if not isinstance(right, IntegerMod_int):
            try:
                return coerce_cmp(self, right)
            except TypeError:
                return -1
        return self.cmp(right)

    def cmp(IntegerMod_int self, IntegerMod_int right):
        """
        EXAMPLES:
            sage: mod(5,13) == mod(-8,13)
            True
            sage: mod(5,13) == mod(8,13)
            False
            sage: mod(5,13) == mod(5,24)
            False
            sage: mod(0, 13) == 0
            True
            sage: mod(0, 13) == int(0)
            True
        """
        if right._parent != self._parent:
            return -1
        if self.ivalue == right.ivalue: return 0
        elif self.ivalue < right.ivalue: return -1
        else: return 1

    def __richcmp__(self, right, int op):
        cdef int n
        if not isinstance(right, IntegerMod_int):
            try:
                n = coerce_cmp(self, right)
            except TypeError:
                n = -1
        else:
            n = self.cmp(right)
        return self._rich_to_bool(op, n)


    def is_one(IntegerMod_int self):
        """
        Returns \\code{True} if this is $1$, otherwise \\code{False}.

        EXAMPLES:
            sage: mod(6,5).is_one()
            True
            sage: mod(0,5).is_one()
            False
        """
        return bool(self.ivalue == 1)

    def is_zero(IntegerMod_int self):
        """
        Returns \\code{True} if this is $0$, otherwise \\code{False}.

        EXAMPLES:
            sage: mod(13,5).is_zero()
            False
            sage: mod(25,5).is_zero()
            True
        """
        return bool(self.ivalue == 0)

    def is_unit(IntegerMod_int self):
        return bool(gcd_int(self.ivalue, self.__modulus.int32) == 1)

    def __crt(IntegerMod_int self, IntegerMod_int other):
        """
        Use the Chinese Remainder Theorem to find an element of the
        integers modulo the product of the moduli that reduces to self
        and to other.  The modulus of other must be coprime to the
        modulus of self.
        EXAMPLES:
            sage: a = mod(3,5)
            sage: b = mod(2,7)
            sage: a.crt(b)
            23

        AUTHOR:
            -- Robert Bradshaw
        """
        cdef IntegerMod_int lift
        cdef int_fast32_t x

        lift = IntegerMod_int(integer_mod_ring.IntegerModRing(self.__modulus.int32 * other.__modulus.int32, check_prime=False), None, empty=True)

        try:
            x = (other.ivalue - self.ivalue % other.__modulus.int32) * mod_inverse_int(self.__modulus.int32, other.__modulus.int32)
            lift.set_from_int( x * self.__modulus.int32 + self.ivalue )
            return lift
        except ZeroDivisionError:
            raise ZeroDivisionError, "moduli must be coprime"


    def copy(IntegerMod_int self):
        cdef IntegerMod_int copy
        copy = IntegerMod_int(self._parent, None, empty=True)
        copy.ivalue = self.ivalue
        return copy

    def _add_(IntegerMod_int self, IntegerMod_int right):
        """
        EXAMPLES:
            sage: R = Integers(10)
            sage: R(7) + R(8)
            5
        """
        cdef IntegerMod_int x
        x = IntegerMod_int(self._parent, None, empty=True)
        x.ivalue = self.ivalue + right.ivalue
        if x.ivalue >= self.__modulus.int32:
            x.ivalue = x.ivalue - self.__modulus.int32
        return x;

    def _sub_(IntegerMod_int self, IntegerMod_int right):
        """
        EXAMPLES:
            sage: R = Integers(10)
            sage: R(7) - R(8)
            9
        """
        cdef IntegerMod_int x
        x = IntegerMod_int(self._parent, None, empty=True)
        x.ivalue = self.ivalue - right.ivalue
        if x.ivalue < 0:
            x.ivalue = x.ivalue + self.__modulus.int32
        return x;

    def _mul_(IntegerMod_int self, IntegerMod_int right):
        """
        EXAMPLES:
            sage: R = Integers(10)
            sage: R(7) * R(8)
            6
        """
        cdef IntegerMod_int x
        x = IntegerMod_int(self._parent, None, empty=True)
        x.ivalue = (self.ivalue * right.ivalue) % self.__modulus.int32
        return x;

    def _div_(IntegerMod_int self, IntegerMod_int right):
        """
        EXAMPLES:
            sage: R = Integers(10)
            sage: R(2)/3
            4
        """
        cdef IntegerMod_int x
        try:
            x = IntegerMod_int(self.parent(), None, empty=True)
            x.imodulus = self.imodulus
            x.ivalue = (self.ivalue * mod_inverse_int(right.ivalue, self.imodulus) ) % self.imodulus
            return x;
        except ZeroDivisionError:
            return IntegerMod_int(self.parent(), self.lift() / right.lift() )

    def __int__(IntegerMod_int self):
        return self.ivalue

    def __long__(IntegerMod_int self):
        return self.ivalue

    def __mod__(IntegerMod_int self, right):
        right = int(right)
        if self.__modulus.int32 % right != 0:
            raise ZeroDivisionError, "Error - reduction modulo right not defined."
        return integer_mod_ring.IntegerModRing(right)(self)

    def __pow__(IntegerMod_int self, right, m): # NOTE: m ignored, always use modulus of parent ring
        """
        EXAMPLES:
            sage: R = Integers(10)
            sage: R(2)^10
            4
            sage: R = Integers(389)
            sage: R(7)^388
            1
        """
        cdef sage.ext.integer.Integer exp, base
        exp = sage.rings.integer_ring.Z(right)
        cdef IntegerMod_int x
        cdef mpz_t x_mpz
        x = IntegerMod_int(self._parent, None, empty=True)
        if mpz_sgn(exp.value) >= 0 and mpz_cmp_si(exp.value, 100000) < 0:  # TODO: test to find a good threshold
            x.ivalue = mod_pow_int(self.ivalue, mpz_get_si(exp.value), self.__modulus.int32)
        else:
            mpz_init(x_mpz)
            _sig_on
            base = self.lift()
            mpz_powm(x_mpz, base.value, exp.value, self.__modulus.sageInteger.value)
            _sig_off
            x.ivalue = mpz_get_si(x_mpz)
            mpz_clear(x_mpz)
        return x


    def __neg__(IntegerMod_int self):
        """
        EXAMPLES:
            sage: -mod(7,10)
            3
        """
        cdef IntegerMod_int x
        x = IntegerMod_int(self._parent, None, empty=True)
        x.ivalue = self.__modulus.int32 - self.ivalue
        return x

    def __invert__(IntegerMod_int self):
        """
        EXAMPLES:
            sage: ~mod(7,100)
            43
        """
        cdef IntegerMod_int x
        x = IntegerMod_int(self._parent, None, empty=True)
        x.ivalue = mod_inverse_int(self.ivalue, self.__modulus.int32)
        return x

    def lift(IntegerMod_int self):
        cdef sage.ext.integer.Integer z
        z = sage.ext.integer.Integer()
        mpz_set_si(z.value, self.ivalue)
        return z

    def __float__(IntegerMod_int self):
        return float(self.ivalue)

    def __hash__(self):
        return hash(self.ivalue)

### End of class


cdef int_fast32_t gcd_int(int_fast32_t a, int_fast32_t b):
    """
    Returns the gcd of a and b
    For use with IntegerMod_int
    AUTHOR:
      -- Robert Bradshaw
    """
    cdef int_fast32_t tmp
    if a < b:
        tmp = b
        b = a
        a = tmp
    while b:
        tmp = b
        b = a % b
        a = tmp
    return a


cdef int_fast32_t mod_inverse_int(int_fast32_t x, int_fast32_t n) except 0:
    """
    Returns y such that xy=1 mod n
    For use in IntegerMod_int
    AUTHOR:
      -- Robert Bradshaw
    """
    cdef int_fast32_t tmp, a, b, last_t, t, next_t, q
    a = n
    b = x
    t = 0
    next_t = 1
    while b:
        # a = s * n + t * x
        if b == 1:
            next_t = next_t % n
            if next_t < 0:
                next_t = next_t + n
            return next_t
        q = a / b
        tmp = b
        b = a % b
        a = tmp
        last_t = t
        t = next_t
        next_t = last_t - q * t
    raise ZeroDivisionError, "Inverse does not exist."


cdef int_fast32_t mod_pow_int(int_fast32_t base, int_fast32_t exp, int_fast32_t n):
    """
    Returns base^exp mod n
    For use in IntegerMod_int
    AUTHOR:
      -- Robert Bradshaw
    """
    cdef int_fast32_t prod, pow2
    if exp <= 5:
        if exp == 0: return 1
        if exp == 1: return base
        prod = base * base % n
        if exp == 2: return prod
        if exp == 3: return (prod * base) % n
        if exp == 4: return (prod * prod) % n

    pow2 = base
    if exp % 2: prod = base
    else: prod = 1
    exp = exp >> 1
    while(exp != 0):
        pow2 = pow2 * pow2
        if pow2 >= INTEGER_MOD_INT32_LIMIT: pow2 = pow2 % n
        if exp % 2:
            prod = prod * pow2
            if prod >= INTEGER_MOD_INT32_LIMIT: prod = prod % n
        exp = exp >> 1

    if prod > n:
        prod = prod % n
    return prod


def test_gcd(a, b):
    return gcd_int(int(a), int(b))

def test_mod_inverse(a, b):
    return mod_inverse_int(int(a), int(b))



