r"""
Elements of $\Z/n\Z$

An element of the integers modulo $n$.
"""

import operator

import integer_mod_ring
import arith
import rational
from sage.libs.all import pari, PariError
import integer_ring
import integer
import commutative_ring_element
import sage.interfaces.all

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

class IntegerMod(commutative_ring_element.CommutativeRingElement):
    def __init__(self, parent, value, construct=False):
        """
        EXAMPLES:
            sage: a = Mod(10,30); a
            10
            sage: loads(a.dumps()) == a
            True
        """
        commutative_ring_element.CommutativeRingElement.__init__(self, parent)
        if construct:
            self.__value = value
            return
        if isinstance(value, rational.Rational):
            value = value % parent.order()
        self.__value = pari('Mod(%s,%s)'%(value, parent.order()))


    def __cmp__(self, right):
        if not isinstance(self, IntegerMod) or not isinstance(right, IntegerMod) \
               or right.parent() != self.parent():
            return coerce.cmp(self, right)
        return self._cmp(right)

    #################################################################
    # Interfaces
    #################################################################
    def _pari_init_(self):
        return str(self.__value)

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
        return arith.discrete_log_generic(self, a, n)


    def modulus(self):
        return self.parent().order()

    def is_square(self):
        return bool(self.__value.issquare())

    def is_unit(self):
        return self.lift().gcd(self.modulus()) == 1

    def sqrt(self):
        """
        Same as self.square_root().
        """
        return self.square_root()

    def square_root(self):
        try:
            return self.parent()(self.__value.sqrt(), construct=True)
        except PariError:
            raise ValueError, "self must be a square."

    def _pari_(self):
        return self.__value

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
        return integer.Integer(self.__value).rational_reconstruction(self.modulus())


    def crt(self, other):
        """
        Use the Chinese Remainder Theorem to find an element of the
        integers modulo the product of the moduli that reduces to self
        and to other.  The modulus of other must be coprime to the
        modulus of self.
        """
        if not isinstance(other, IntegerMod):
            raise TypeError, "other must be an integer mod"
        return self.parent()(self.__value.chinese(other.__value), construct=True)

    def order(self):
        """
        Returns the additive order of self.
        """
        n = self.modulus()
        return integer.Integer(n//int(self.__value.lift().gcd(n)))

    def multiplicative_order(self):
        """
        Returns the additive order of self.
        """
        if not self.is_unit():
            raise ArithmeticError, "self must be a unit"
        return integer.Integer(self.__value.order())  # pari's "order" is by default multiplicative

    def copy(self):
        return self.parent()(self.__value, construct=True)

    def _repr_(self):
        return str(self.__value.lift())

    def _latex_(self):
        return str(self)

    def _add_(self, right):
        return self.parent()(self.__value + right.__value, construct=True)

    def _sub_(self, right):
        return self.parent()(self.__value - right.__value, construct=True)

    def _mul_(self, right):
        return self.parent()(self.__value * right.__value, construct=True)

    def _div_(self, right):
        """
        EXAMPLES:
            sage: R = Integers(25)
            sage: R(15)/5
            3
        """
        try:
            return self.parent()(self.__value / right.__value, construct=True)
        except PariError:
            P = self.parent()
            t = (self.__value.lift() / right.__value.lift()).Mod(P.order())
            return P(t, construct=True)

    def __int__(self):
        return int(self.__value.lift())

    def _integer_(self):
        return integer.Integer(self.__value.lift())

    def _rational_(self):
        return rational.Rational(self.__value.lift())

    def __long__(self):
        return long(self.__value.lift())

    def __mod__(self, right):
        right = int(right)
        if self.modulus() % right != 0:
            raise ZeroDivisionError, "Error - reduction modulo right not defined."
        return integer_mod_ring.IntegerModRing(right)(self)

    # commented out because PARI (used for .__value) prints
    # out crazy warnings when the exponent is LARGE -- this
    # is even a problem in gp!!!
    #def __pow__(self, right):
    #    right = int(right)
    #    return self.parent()(self.__value**right, construct=True)

    def __neg__(self):
        return self.parent()(-self.__value, construct=True)

    def __invert__(self):
        return self.parent()(~self.__value, construct=True)

    def lift(self):
        return integer.Integer(self.__value.lift())

    def __float__(self):
        return float(self.__value.lift())

    def __cmp__(self, other):
        if not isinstance(other, IntegerMod) or self.parent() != other.parent():
            return coerce.cmp(self, other)
        return cmp(self.__value, other.__value)

