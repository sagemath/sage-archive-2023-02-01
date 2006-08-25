r"""
Elements of $\Z/n\Z$

An element of the integers modulo $n$.
"""

include "../ext/interrupt.pxi"  # ctrl-c interrupt block support


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

import sage.rings.coerce as coerce

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


def is_IntegerMod(x):
    return isinstance(x, (integer_mod.IntegerMod, integer_mod_int.IntegerMod_int))


cdef public int set_from_mpz

cdef class IntegerMod(sage.ext.element.CommutativeRingElement):

    def __new__(self, parent, x=None, empty=False):
        mpz_init(self.value)

    def __init__(self, parent, value, construct=False, empty=False):
        """
        EXAMPLES:
            sage: a = Mod(10,30^10); a
            10
            sage: loads(a.dumps()) == a
            True
        """
        commutative_ring_element.CommutativeRingElement.__init__(self, parent)
        if construct: # TODO: can I ever use this?
            self.__value = value
            return
        if empty:
            return
        if isinstance(value, rational.Rational):
            value = value % parent.order()
        elif not isinstance(value, sage.ext.integer.Integer):
            value = sage.rings.integer_ring.Z(value)
        cdef sage.ext.integer.Integer z
        z = value
        self.set_from_mpz(z.value)

    def __dealloc__(self):
        mpz_clear(self.value)

    def __cmp__(self, right):
        """
        EXAMPLES:
            sage: mod(5,13) == mod(-8,13)
            True
        """
        if not isinstance(self, IntegerMod) or not isinstance(right, IntegerMod) \
               or right.parent() != self.parent():
            return coerce_cmp(self, right)
        cdef IntegerMod x
        x = right
        return mpz_cmp(self.value, x.value)

    cdef void set_from_mpz(IntegerMod self, mpz_t value):
        cdef sage.ext.integer.Integer modulus
        modulus = self.modulus()
        if mpz_sgn(value) == -1 or mpz_cmp(value, modulus.value) >= 0:
            mpz_mod(self.value, value, modulus.value)
        else:
            mpz_set(self.value, value)

    cdef mpz_t* get_value(IntegerMod self):
        return &self.value

    cdef mpz_t* mpz_modulus(IntegerMod self):
        cdef sage.ext.integer.Integer modulus
        modulus = self.modulus()
        return &modulus.value



    #################################################################
    # Interfaces
    #################################################################
    def _pari_init_(self):
        return 'Mod(%s,%s)'%(str(self), self.parent().order())

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

        if R.order() > 65536:
            raise ValueError, "order must be at most 65536."

        if self == 0:
            return '0*Z(%s)'%R.order()

        # I couldn't find a guarentee in the GAP docs that the
        # root of unity they use must be the smallest.   This
        # was *not* the case in MAGMA once, so who knows, especially
        # when the order of the ring is not prime.  So we make
        # no such dangerous assumptions (for now).

        # Find the root of unity used by Gap.
        m = R.order()
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
        n = (self.parent()).unit_group_order()
        return arith.discrete_log_generic(self, a, n) # TODO update this function


    def modulus(self):
        return self.parent().order()


    def is_one(IntegerMod self):
        """
        Returns \\code{True} if this is $1$, otherwise \\code{False}.

        EXAMPLES:
            sage: mod(6,5).is_one()
            True
            sage: mod(0,5).is_one()
            False
        """
        return bool(mpz_cmp_si(self.value, 1) == 0)

    def is_zero(IntegerMod self):
        """
        Returns \\code{True} if this is $0$, otherwise \\code{False}.

        EXAMPLES:
            sage: mod(13,5).is_zero()
            False
            sage: mod(25,5).is_zero()
            True
        """
        return bool(mpz_cmp_si(self.value, 0) == 0)

    def is_square(self):
        return bool(self.pari().issquare()) # TODO implement directly

    def is_unit(self):
        return bool(self.lift().gcd(self.modulus()) == 1)

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
        R = polynomial_ring.PolynomialRing(self.parent())
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
        return lift().rational_reconstruction(self.modulus())


    def crt(IntegerMod self, py_other):
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
            sage: a = m.mod(37,10^8)
            sage: b = m.mod(9,3^8)
            sage: a.crt(b)
            125900000037
        AUTHOR:
            -- Robert Bradshaw
        """
        if isinstance(py_other, IntegerMod_int):
            py_other = IntegerMod(py_other.parent(), py_other.lift())

        if not isinstance(py_other, IntegerMod):
            raise TypeError, "other must be an integer mod"

        cdef IntegerMod lift, other, x
        cdef sage.ext.integer.Integer modulus, other_modulus

        other = py_other
        modulus = self.modulus()
        other_modulus = other.modulus()
        lift = IntegerMod(integer_mod_ring.IntegerModRing(modulus*other_modulus, check_prime=False), None, empty=True)
        try:
          if mpz_cmp(self.value, other.value) > 0:
              x = (other - IntegerMod(other.parent(), self.lift())) / IntegerMod(other.parent(), modulus)
              mpz_mul(lift.value, x.value, modulus.value)
              mpz_add(lift.value, lift.value, self.value)
          else:
              x = (self - IntegerMod(self.parent(), other.lift())) / IntegerMod(self.parent(), other_modulus)
              mpz_mul(lift.value, x.value, other_modulus.value)
              mpz_add(lift.value, lift.value, other.value)
          return lift
        except ZeroDivisionError:
            raise ZeroDivisionError, "moduli must be coprime"


    def order(self):
        """
        Returns the additive order of self.
        """
        n = self.modulus()
        return integer.Integer(n.__FLOORDIV__(self.lift().gcd(n)))

    def multiplicative_order(self):
        """
        Returns the additive order of self.
        """
        if not self.is_unit():
            raise ArithmeticError, "self must be a unit"
        return integer.Integer(self.pari().order())  # pari's "order" is by default multiplicative  # TODO: implement this directly

    def copy(IntegerMod self):
        cdef IntegerMod copy
        copy = IntegerMod(self.parent(), None, empty=True)
        mpz_set(copy.value, self.value)
        return copy

    def _repr_(self):
        return str(self.lift())

    def _latex_(self):
        return str(self)

    def _add_(IntegerMod self, IntegerMod right):
        """
        EXAMPLES:
            sage: R = Integers(10)
            sage: R(7) + R(8)
            5
        """
        cdef sage.ext.integer.Integer modulus
        modulus = self.modulus()
        cdef IntegerMod x
        x = IntegerMod(self.parent(), None, empty=True)
        mpz_add(x.value, self.value, right.value)
        if mpz_cmp(x.value, modulus.value)  >= 0:
            mpz_sub(x.value, x.value, modulus.value)
        return x;

    def _sub_(IntegerMod self, IntegerMod right):
        """
        EXAMPLES:
            sage: R = Integers(10)
            sage: R(7) - R(8)
            9
        """
        cdef sage.ext.integer.Integer modulus
        modulus = self.modulus()
        cdef IntegerMod x
        x = IntegerMod(self.parent(), None, empty=True)
        mpz_sub(x.value, self.value, right.value)
        if mpz_sgn(x.value) == -1:
            mpz_add(x.value, x.value, modulus.value)
        return x;

    def _mul_(IntegerMod self, IntegerMod right):
        """
        EXAMPLES:
            sage: R = Integers(10)
            sage: R(7) * R(8)
            6
        """
        cdef sage.ext.integer.Integer modulus
        modulus = self.modulus()
        cdef IntegerMod x
        x = IntegerMod(self.parent(), None, empty=True)
        mpz_mul(x.value, self.value, right.value)
        mpz_fdiv_r(x.value, x.value, modulus.value)
        return x;

    def _div_(IntegerMod self, IntegerMod right):
        """
        EXAMPLES:
            sage: R = Integers(25)
            sage: R(15)/5
            3
        """
        try:
            return self._mul_(~right)
        except ZeroDivisionError:
            return IntegerMod(self.parent(), self.lift() / right.lift() )

    def __floordiv__(self, right):
        return self._div_(right)

    def __int__(self):
        return int(self.lift())

    def _integer_(self):
        return self.lift()

    def _rational_(self):
        return rational.Rational(self.lift())

    def __long__(self):
        return long(self.lift())

    def __mod__(self, right):
        right = int(right)
        if self.modulus() % right != 0:
            raise ZeroDivisionError, "Error - reduction modulo right not defined."
        return integer_mod_ring.IntegerModRing(right)(self)

    def __pow__(IntegerMod self, right, m): # NOTE: m ignored, always use modulus of parent ring
        """
        EXAMPLES:
            sage: R = Integers(10)
            sage: R(2)^10
            4
            sage: R = Integers(389)
            sage: R(7)^388
            1
        """
        cdef sage.ext.integer.Integer modulus, exp
        modulus = self.modulus()
        exp = sage.rings.integer_ring.Z(right)
        cdef IntegerMod x
        x = IntegerMod(self.parent(), None, empty=True)
        _sig_on
        mpz_powm(x.value, self.value, exp.value, modulus.value)
        _sig_off
        return x

    def __neg__(IntegerMod self):
        """
        EXAMPLES:
            sage: -mod(7,10)
            3
        """
        cdef sage.ext.integer.Integer modulus
        modulus = self.modulus()
        cdef IntegerMod x
        x = IntegerMod(self.parent(), None, empty=True)
        mpz_sub(x.value, modulus.value, self.value)
        return x

    def __invert__(IntegerMod self):
        """
        EXAMPLES:
            sage: ~mod(7,100)
            43
        """
        cdef sage.ext.integer.Integer modulus
        modulus = self.modulus()
        cdef IntegerMod x
        x = IntegerMod(self.parent(), None, empty=True)
        _sig_on
        if (mpz_invert(x.value, self.value, modulus.value)):
            _sig_off
            return x
        else:
            _sig_off
            raise ZeroDivisionError, "Inverse does not exist."

    def lift(IntegerMod self):
        cdef sage.ext.integer.Integer z
        z = sage.ext.integer.Integer()
        z.set_from_mpz(self.value)
        return z

    def __float__(self):
        return float(self.lift())


######################################################################
#      class IntegerMod_int
######################################################################
"""
Elements of $\Z/n\Z$ for n small enough to be operated on in word size
TODO: modify constructors/ring to call this class for n<INTEGER_MOD_INT_LIMIT
AUTHORS:
    -- Robert Bradshaw (2006-08-24)
"""

# NOTE: I've gone back and forth between making IntegerMod_int
#       a descendant or sibling class of IntegerMod
#   Pros of descending:
#     - lots of duplicate functionality that does not depend on representation
#     - easier isinstance checking
#   Cons sibling:
#     - mpz_t value cannot be assumed to be initalized/used
#     - must be sure to override all added functionality thad depends on representation

cdef class IntegerMod_int(IntegerMod):

    def __init__(self, parent, value, construct=False, empty=False):
        """
        EXAMPLES:
            sage: a = Mod(10,30); a
            10
            sage: loads(a.dumps()) == a
            True
        """
        commutative_ring_element.CommutativeRingElement.__init__(self, parent)
        if empty:
            return
        if isinstance(value, rational.Rational):
            value = value % parent.order()
        elif not isinstance(value, sage.ext.integer.Integer):
            value = sage.rings.integer_ring.Z(value)
        cdef sage.ext.integer.Integer modulus
        modulus = self.modulus()
        self.imodulus = (int)(mpz_get_si(modulus.value))
        cdef sage.ext.integer.Integer z
        z = value
        self.set_from_mpz(z.value)

    def __cmp__(IntegerMod_int self, right):
        """
        EXAMPLES:
            sage: mod(5,13) == mod(-8,13)
            True
        """
        if not isinstance(self, IntegerMod_int) or not isinstance(right, IntegerMod_int) \
               or right.parent() != self.parent():
            return coerce_cmp(self, right)
        cdef IntegerMod_int x
        x = right
        if self.ivalue == x.ivalue: return 0
        elif self.ivalue < x.ivalue: return -1
        else: return 1

    cdef void set_from_mpz(IntegerMod_int self, mpz_t value):
        if mpz_sgn(value) == -1 or mpz_cmp_si(value, self.imodulus) >= 0:
            self.ivalue = (int)(mpz_fdiv_ui(value, self.imodulus))
        else:
            self.ivalue = (int)(mpz_get_si(value))

    cdef void set_from_int(IntegerMod_int self, int ivalue):
        if ivalue < 0:
            self.ivalue = self.imodulus + (ivalue % self.imodulus)
        elif ivalue >= self.imodulus:
            self.ivalue = ivalue % self.imodulus
        else:
            self.ivalue = ivalue

    cdef int get_int_value(IntegerMod_int self):
        return self.ivalue


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
        return bool(gcd_int(self.ivalue, self.imodulus) == 1)

    def crt(IntegerMod_int self, py_other):
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
            sage: a = m.mod(37,10^8)
            sage: b = m.mod(9,3^8)
            sage: a.crt(b)
            125900000037
        AUTHOR:
            -- Robert Bradshaw
        """
        if not isinstance(py_other, IntegerMod):
            raise TypeError, "other must be an integer mod"

        if not isinstance(py_other, IntegerMod_int):
            return py_other.crt(self)

        cdef IntegerMod_int lift, other
        other = py_other
        if self.imodulus * other.imodulus > INTEGER_MOD_INT_LIMIT:
            return IntegerMod(self.parent(), lift()).crt(py_other)

        lift = IntegerMod(integer_mod_ring.IntegerModRing(self.imodulus * other.imodulus, check_prime=False), None, empty=True)
        lift.imodulus = self.imodulus * other.imodulus

        try:
          lift.set_from_int( (other.ivalue - self.ivalue % other.imodulus) * mod_inverse_int(other.imodulus, self.imodulus) )
          return lift
        except ZeroDivisionError:
            raise ZeroDivisionError, "moduli must be coprime"


    def copy(IntegerMod_int self):
        cdef IntegerMod_int copy
        copy = IntegerMod_int(self.parent(), None, empty=True)
        copy.imodulus = self.imodulus
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
        x = IntegerMod_int(self.parent(), None, empty=True)
        x.imodulus = self.imodulus
        x.ivalue = self.ivalue + right.ivalue
        if x.ivalue >= self.imodulus:
            x.ivalue = x.ivalue - self.imodulus
        return x;

    def _sub_(IntegerMod_int self, IntegerMod_int right):
        """
        EXAMPLES:
            sage: R = Integers(10)
            sage: R(7) - R(8)
            9
        """
        cdef IntegerMod_int x
        x = IntegerMod_int(self.parent(), None, empty=True)
        x.imodulus = self.imodulus
        x.ivalue = self.ivalue - right.ivalue
        if x.ivalue < 0:
            x.ivalue = x.ivalue + self.imodulus
        return x;

    def _mul_(IntegerMod_int self, IntegerMod_int right):
        """
        EXAMPLES:
            sage: R = Integers(10)
            sage: R(7) * R(8)
            6
        """
        cdef IntegerMod_int x
        x = IntegerMod_int(self.parent(), None, empty=True)
        x.imodulus = self.imodulus
        x.ivalue = (self.ivalue * right.ivalue) % self.imodulus
        return x;

    def _div_(IntegerMod_int self, IntegerMod_int right):
        """
        EXAMPLES:
            sage: R = Integers(25)
            sage: R(15)/5
            3
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
        return int(self.ivalue)

    def __long__(IntegerMod_int self):
        return long(self.ivalue)

    def __mod__(IntegerMod_int self, right):
        right = int(right)
        if self.imodulus % right != 0:
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
        cdef sage.ext.integer.Integer modulus, exp, base
        modulus = self.modulus()
        base = self.lift()
        exp = sage.rings.integer_ring.Z(right)
        cdef IntegerMod_int x
        x = IntegerMod_int(self.parent(), None, empty=True)
        x.imodulus = self.imodulus
        cdef mpz_t x_mpz
        mpz_init(x_mpz)
        _sig_on
        mpz_powm(x_mpz, base.value, exp.value, modulus.value)
        _sig_off
        x.ivalue = (int)(mpz_get_ui(x_mpz))
        mpz_clear(x_mpz)
        return x

    def __neg__(IntegerMod_int self):
        """
        EXAMPLES:
            sage: -mod(7,10)
            3
        """
        cdef IntegerMod_int x
        x = IntegerMod_int(self.parent(), None, empty=True)
        x.imodulus = self.imodulus
        x.ivalue = self.imodulus - self.ivalue
        return x

    def __invert__(IntegerMod_int self):
        """
        EXAMPLES:
            sage: ~mod(7,100)
            43
        """
        cdef sage.ext.integer.Integer modulus
        modulus = self.modulus()
        cdef IntegerMod_int x
        x = IntegerMod_int(self.parent(), None, empty=True)
        x.imodulus = self.imodulus
        x.ivalue = mod_inverse_int(self.ivalue, self.imodulus)
        return x

    def lift(IntegerMod_int self):
        cdef sage.ext.integer.Integer z
        z = sage.ext.integer.Integer()
        mpz_set_si(z.value, self.ivalue)
        return z

    def __float__(IntegerMod_int self):
        return float(self.ivalue)



cdef int gcd_int(int a, int b):
    """
    Returns the gcd of a and b
    For use with IntegerMod_int
    AUTHOR:
      -- Robert Bradshaw
    """
    cdef int tmp
    if a < b:
        tmp = b
        b = a
        a = tmp
    while b:
        tmp = b
        b = a % b
        a = tmp
    return a


cdef int mod_inverse_int(int x, int n) except 0:
    """
    Returns y such that xy=1 mod n
    For use in IntegerMod_int
    AUTHOR:
      -- Robert Bradshaw
    """
    cdef int tmp, a, b, last_t, t, next_t, q
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


def test_gcd(a, b):
    return gcd_int(int(a), int(b))

def test_mod_inverse(a, b):
    return mod_inverse_int(int(a), int(b))



