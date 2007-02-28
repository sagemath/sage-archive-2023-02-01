"""
Elements of p-Adic Rings with Fixed Modulus

AUTHOR:
    -- David Roe
    -- Genya Zaytman: documentation
    -- David Harvey: doctests
"""

#*****************************************************************************
#       Copyright (C) 2007 David Roe <roed@math.harvard.edu>
#                          William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

#import sage.rings.arith
#import sage.rings.infinity
#import sage.libs.all
import sage.rings.padics.padic_ring_generic_element
import sage.rings.padics.padic_field_generic_element
import sage.rings.padics.padic_lazy_element
#import sage.rings.padics.padic_field_capped_relative_element
#import sage.rings.commutative_ring_element
#import sage.rings.integer_mod
#import sage.rings.finite_field_element
#import sage.rings.integer
#import sage.rings.rational
#import sage.rings.padics.precision_error

infinity = sage.rings.infinity.infinity
PrecisionError = sage.rings.padics.precision_error.PrecisionError
Integer = sage.rings.integer.Integer
Rational = sage.rings.rational.Rational
#Zp = sage.rings.padics.padic_ring_generic.Zp
Mod = sage.rings.integer_mod.Mod
pAdicLazyElement = sage.rings.padics.padic_lazy_element.pAdicLazyElement
pAdicRingGenericElement = sage.rings.padics.padic_ring_generic_element.pAdicRingGenericElement
pAdicGenericElement = sage.rings.padics.padic_generic_element.pAdicGenericElement
#pAdicCappedRelativeElement = sage.rings.padics.padic_field_capped_relative_element.pAdicFieldCappedRelativeElement
pari = sage.libs.pari.gen.pari
pari_gen = sage.libs.pari.gen.gen
PariError = sage.libs.pari.gen.PariError

class pAdicRingFixedModElement(pAdicRingGenericElement):
    def __init__(self, parent, x, construct=False):
        r"""
        INPUT:
            parent -- a pAdicRingFixedMod object.
        Types currently supported:
            Integers
            Rationals -- denominator must be relatively prime to p
            FixedMod p-adics
        Types that should be supported:
            Finite precision p-adics
            Lazy p-adics
            Elements of local extensions of THIS p-adic ring that actually lie in Zp
            Elements of IntegerModRing(p^k) for k less than or equal to the modulus
        """
        sage.rings.commutative_ring_element.CommutativeRingElement.__init__(self,parent)
        if construct:
            self._value = x
            return

        if isinstance(x, pAdicGenericElement) and x.valuation() < 0:
            raise ValueError, "element has negative valuation"
        if isinstance(x, pAdicLazyElement):
            try:
                x.set_precision_absolute(prec)
            except PrecisionError:
                pass
        if isinstance(x, pAdicGenericElement):
            self._value = Mod(x.lift(), self.parent().prime_pow(self.parent().precision_cap()))
            return

        if isinstance(x, pari_gen) and x.type() == "t_PADIC":
            t = x.lift()
            if t.type() == 't_INT':
                x = int(t)
            else:
                x = QQ(t)

        if sage.rings.integer_mod.is_IntegerMod(x):
            k, p = pari(x.modulus()).ispower()
            if not k or p != parent.prime():
                raise TypeError, "cannot change primes in creating p-adic elements"
            x = x.lift()

        if sage.rings.finite_field_element.is_FiniteFieldElement(x):
            if x.parent().order() != parent.prime():
                raise TypeError, "can only create p-adic element out of finite field when order of field is p"
            x = x.lift()

        #Now use the code below to convert from integer or rational, so don't make the next line elif

        if isinstance(x, Integer):
            self._value = Mod(x, self.parent().prime_pow(self.parent().precision_cap()))
        elif isinstance(x, Rational):
            if x.valuation(self.parent().prime()) < 0:
                raise ValueError, "p divides the denominator"
            else:
                self._value = Mod(x, self.parent().prime_pow(self.parent().precision_cap()))
        elif isinstance(x, (int, long)):
            self._value = Mod(x, self.parent().prime_pow(self.parent().precision_cap()))
        else:
            raise TypeError, "unable to create p-adic element"

    def __invert__(self):
        r"""
        Returns multiplicative inverse of this element. Its valuation
        must be zero.

        EXAMPLES:
            sage: R = Zp(7, 4, 'fixed-mod', 'series')
            sage: ~R(2)
            4 + 3*7 + 3*7^2 + 3*7^3 + O(7^4)
            sage: ~R(0)
            Traceback (most recent call last):
            ...
            ValueError: cannot invert non-unit
            sage: ~R(7)
            Traceback (most recent call last):
            ...
            ValueError: cannot invert non-unit
        """
        if self.valuation() > 0:
            raise ValueError, "cannot invert non-unit"
        else:
            # todo: use ~(self._value) perhaps? -- dmharvey
            inverse = 1 / self._value
            return pAdicRingFixedModElement(self.parent(), inverse, construct = True)

    def __mod__(self, right):
        r"""
        Returns self modulo right.

        This doesn't make a whole lot of sense :-) but things are defined so
        that always (x // y) * y + (x % y) == x.

        TODO:
            -- write down a full and proper explanation of the exact semantics
            of floordiv and mod for all padic rings/fields; this should go
            in a single overall doc file somewhere, and a reference to it
            plus a summary should go in this docstring
            -- make this work when "right" is e.g. an Integer -- perhaps
            the mod operator needs to be brought under the SAGE arithmetic
            architecture umbrella

        """
        val = self.valuation()
        rval = right.valuation()
        quotient =  self / right._unit_part()
        return pAdicRingFixedModElement(self.parent(),
                     Mod(quotient.lift() % self.parent().prime_pow(rval),
                     self.parent().prime_pow(self.parent().precision_cap())),
                     construct = True)

    def _neg_(self):
        r"""
        Returns negative of self.

        EXAMPLES:
            sage: R = Zp(7, 4, 'fixed-mod', 'series')
            sage: -R(7)
            6*7 + 6*7^2 + 6*7^3 + O(7^4)
        """
        return pAdicRingFixedModElement(self.parent(), -self._value, construct = True)


    def __pow__(self, right):
        right = Integer(right) #Need to make sure that this works for p-adic exponents
        return pAdicRingFixedModElement(self.parent(), self._value**right, construct = True)

    def _add_(self, right):
        r"""
        Returns sum of self and right.

        EXAMPLES:
            sage: R = Zp(7, 4, 'fixed-mod', 'series')
            sage: x = R(1721); x
            6 + 5*7^3 + O(7^4)
            sage: y = R(1373); y
            1 + 4*7^3 + O(7^4)
            sage: x + y
            7 + 2*7^3 + O(7^4)
        """
        return pAdicRingFixedModElement(self.parent(),
                            self._value + right._value, construct = True)

    def _div_(self, right):
        r"""
        Returns quotient of self and right. The latter must have
        valuation zero.

        EXAMPLES:
            sage: R = Zp(7, 4, 'fixed-mod', 'series')
            sage: R(3) / R(2)
            5 + 3*7 + 3*7^2 + 3*7^3 + O(7^4)
            sage: R(5) / R(0)
            Traceback (most recent call last):
            ...
            ValueError: cannot invert non-unit
            sage: R(7) / R(49)
            Traceback (most recent call last):
            ...
            ValueError: cannot invert non-unit
        """
        if right.valuation() > 0:
            raise ValueError, "cannot invert non-unit"
        else:
            return pAdicRingFixedModElement(self.parent(),
                            self._value / right.value, construct = True)

    def __floordiv__(self, right):
        if isinstance(right, Integer):
            right = pAdicRingFixedModElement(self.parent(), right)
        ppow = self.parent().prime_pow(right.valuation())
        runit = right._unit_part()
        quotient = Mod((self._value / runit).lift() // ppow, self.parent().prime_pow(self.parent().precision_cap()))
        return pAdicRingFixedModElement(self.parent(), quotient, construct = True)

    def _integer_(self):
        r"""
        Return an integer congruent to self modulo self's precision.

        See the lift() method.
        """
        return self._value.lift()

    def _mul_(self, right):
        r"""
        Returns product of self and right.

        EXAMPLES:
            sage: R = Zp(7, 4, 'fixed-mod', 'series')
            sage: R(3) * R(2)
            7 + O(7^4)
            sage: R(1/2) * R(2)
            1 + O(7^4)
        """
        return pAdicRingFixedModElement(self.parent(), self._value * right._value, construct = True)

    def _sub_(self, right):
        r"""
        Returns difference of self and right.

        EXAMPLES:
            sage: R = Zp(7, 4, 'fixed-mod', 'series')
            sage: x = R(1721); x
            6 + 5*7^3 + O(7^4)
            sage: y = R(1373); y
            1 + 4*7^3 + O(7^4)
            sage: x - y
            5 + 7^3 + O(7^4)
        """
        return pAdicRingFixedModElement(self.parent(), self._value - right._value, construct = True)

    def add_bigoh(self, prec):
        """
        Returns a new element with precision decreased to prec
        INPUT:
            self -- a p-adic element
            prec -- an integer
        OUTPUT:
            element -- self with precision set to the minimum of  self's precision and prec

        EXAMPLES:
            sage: R = Zp(7,4,'fixed-mod','series'); a = R(8); a.add_bigoh(1)
            1 + O(7^4)
        """
        return pAdicRingFixedModElement(self.parent(), Mod(Mod(self._value, self.parent().prime_pow(prec)), self.parent().prime_pow(self.parent().precision_cap())), construct = True)

    def copy(self):
        return pAdicRingFixedModElement(self.parent(), self._value, construct = True)

    def exp_artin_hasse(self):
        raise NotImplementedError

    def gamma(self):
        raise NotImplementedError

    def is_zero(self, prec):
        r"""
        Returns whether self is zero modulo $p^{\mbox{prec}}$.

        INPUT:
            self -- a p-adic element
            prec -- an integer
        OUTPUT:
            boolean -- whether self is zero

        """
        return Mod(self._value, self.parent().prime_pow(prec)) == 0

    def is_equal_to(self, right, prec): #assumes they have the same parent
        r"""
        Returns whether self is equal to right modulo $p^{\mbox{prec}}$.

        INPUT:
            self -- a p-adic element
            right -- a p-addic element
            prec -- an integer
        OUTPUT:
            boolean -- whether self is equal to right

        """
        return Mod(self._value, self.parent().prime_pow(prec)) == Mod(right._value, self.parent().prime_pow(prec))

    def lift(self):
        r"""
        Return an integer congruent to self modulo self's precision.

        INPUT:
            self -- a p-adic element
        OUTPUT:
            integer -- a integer congruent to self mod $p^{\mbox{prec}}$
        EXAMPLES:
            sage: R = Zp(7,4,'fixed-mod'); a = R(8); a.lift()
            8
            sage: type(a.lift())
            <type 'sage.rings.integer.Integer'>
        """
        return self._value.lift()

    def list(self):
        r"""
        Returns a list of coefficients of p starting with $p^0$.

        INPUT:
            self -- a p-adic element
        OUTPUT:
            list -- the list of coefficients of self. Precisely $n$
                coefficients are returned, where the modulus of this ring
                is $p^n$.
        EXAMPLES:
            sage: R = Zp(7,4,'fixed-mod'); a = R(2*7+7**2); a.list()
            [0, 2, 1, 0]

        NOTE:
            this differs from the list method of padic_field_element
        """
        def plist(n, p):
            if n == 0:
                return []
            else:
                return [n % p] + plist(n // p, p)
        rep = plist(self._value.lift(), self.parent().prime())
        return rep + [0 for w in range(len(rep), self.parent().precision_cap())]

    def log_artin_hasse(self):
        raise NotImplementedError

    def multiplicative_order(self):
        r"""
        Returns the multiplicative order of self, where self is considered to
        be 1 if it is 1 modulo $p^{\mbox{prec}}$.

        INPUT:
            self -- a p-adic element
            prec -- an integer
        OUTPUT:
            integer -- the multiplicative order of self
        """
        if self._value % self.parent().prime() == 0:
            return infinity
        if self._value == 1:
            return Integer(1)
        if self._value == -1:
            return Integer(2)
        elif self.is_equal_to(self.parent().teichmuller(self.unit_part().residue(1)), prec): #need to improve efficiency
            return self.residue(1).multiplicative_order()
        else:
            return infinity

    def precision_absolute(self):
        """
        Returns the absolute precision of self
         INPUT:
            self -- a p-adic element
        OUTPUT:
            integer -- the absolute precision of self
        EXAMPLES:
            sage: R = Zp(7,4,'fixed-mod'); a = R(7); a.precision_absolute()
            4
        """
        return self.parent().precision_cap()

    def precision_relative(self):
        r"""
        Returns the relative precision of self
         INPUT:
            self -- a p-adic element
        OUTPUT:
            integer -- the relative precision of self
        EXAMPLES:
            sage: R = Zp(7,4,'fixed-mod'); a = R(7); a.precision_relative()
            3
        """
        return self.parent().precision_cap() - self.valuation()

    def residue(self, prec):
        r"""
        Reduces this mod $p^prec$
        INPUT:
            self -- a p-adic element
            prec - an integer
        OUTPUT:
            element of Z/(p^prec Z) -- self reduced mod p^prec
        EXAMPLES:
            sage: R = Zp(7,4,'fixed-mod'); a = R(8); a.residue(1)
            1
        """
        return Mod(self._value, self.parent().prime_pow(prec))

    def square_root(self):
        r"""
        Returns the square root of this p-adic number

        INPUT:
            self -- a p-adic element
        OUTPUT:
            p-adic element -- the square root of this p-adic number

            The square root chosen is the one whose reduction mod p is in
            the range [0, p/2).

            Note that because this is a fixed modulus ring, garbage digits
            may be introduced, if either
            (a) the valuation of the input is positive, or
            (b) p = 2.

            If no square root exists, a ValueError is raised.
            (This may be changed later to return an element of an extension
            field.)

        EXAMPLES:
            sage: R = Zp(3,20,'fixed-mod')
            sage: R(0).square_root()
                O(3^20)
            sage: R(1).square_root()
                1 + O(3^20)
            sage: R(2).square_root()
            Traceback (most recent call last):
            ...
            ValueError: element is not a square
            sage: R(4).square_root() == R(-2)
                True
            sage: R(9).square_root()
                3 * 1 + O(3^20)
            sage: R2 = Zp(2,20,'fixed-mod')
            sage: R2(0).square_root()
                O(2^20)
            sage: R2(1).square_root()
                1 + O(2^20)
            sage: R2(4).square_root()
                2 * 1 + O(2^20)
            sage: R2(9).square_root() == R2(3) or R2(9).square_root() == R2(-3)
                True
            sage: R2(17).square_root()
                206569 + O(2^20)
            sage: R3 = Zp(5,20,'fixed-mod')
            sage: R3(0).square_root()
                O(5^20)
            sage: R3(1).square_root()
                1 + O(5^20)
            sage: R3(-1).square_root() == R3.teichmuller(2) or R3(-1).square_root() == R3.teichmuller(3)
                True
        """
        try:
            # use pari
            return self.parent()(pari(self).sqrt())
        except PariError:
            # todo: should eventually change to return an element of
            # an extension field
            raise ValueError, "element is not a square"


    def _unit_part(self):
        return Mod(self._value.lift() // self.parent().prime_pow(self.valuation()), self.parent().prime_pow(self.parent().precision_cap()))

    def unit_part(self):
        r"""
        Returns the unit part of self.

        INPUT:
            self -- a p-adic element
        OUTPUT:
            p-adic element -- the unit part of self
        EXAMPLES:
            sage: R = Zp(17,4,'fixed-mod')
            sage: a = R(18*17)
            sage: a.unit_part()
            18 + O(17^4)
            sage: type(a)
            <class 'sage.rings.padics.padic_ring_fixed_mod_element.pAdicRingFixedModElement'>
        """
        if self._value == 0:
            raise PrecisionError, "Not enough precision to determine unit part."
        return pAdicRingFixedModElement(self.parent(), self._unit_part(), construct = True)

    def valuation(self):
        """
        Returns the valuation of self.

        INPUT:
            self -- a p-adic element
        OUTPUT:
            integer -- the valuation of self

        EXAMPLES:
            sage: R = Zp(17, 4,'fixed-mod')
            sage: a = R(2*17^2)
            sage: a.valuation()
            2
            sage: R = Zp(5, 4,'fixed-mod')
            sage: R(0).valuation()
            4
            sage: R(1).valuation()
            0
            sage: R(2).valuation()
            0
            sage: R(5).valuation()
            1
            sage: R(10).valuation()
            1
            sage: R(25).valuation()
            2
            sage: R(50).valuation()
            2
        """
        val = sage.rings.arith.valuation(self.lift(),self.parent().prime())
        if val == infinity:
            return self.parent().precision_cap()
        else:
            return val
