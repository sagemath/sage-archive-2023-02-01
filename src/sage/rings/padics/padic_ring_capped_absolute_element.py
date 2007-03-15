"""
Elements of p-Adic Rings with Absolute Precision Cap

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

#import sage.rings.integer
#import sage.rings.rational
import sage.rings.padics.padic_ring_fixed_mod_element
#import sage.rings.padics.padic_ring_capped_relative_element
#import sage.rings.padics.padic_ring_generic_element
#import sage.libs.all
#import sage.rings.infinity
#import sage.rings.arith
#import sage.rings.integer_mod
#import sage.rings.finite_field_element
#import sage.rings.padics.precision_error

Mod = sage.rings.integer_mod.Mod
Integer = sage.rings.integer.Integer
Rational = sage.rings.rational.Rational
infinity = sage.rings.infinity.infinity
#Zp = sage.rings.padics.padic_ring.Zp
pAdicRingFixedModElement = sage.rings.padics.padic_ring_fixed_mod_element.pAdicRingFixedModElement
pAdicGenericElement = sage.rings.padics.padic_generic_element.pAdicGenericElement
pAdicRingGenericElement = sage.rings.padics.padic_ring_generic_element.pAdicRingGenericElement
pAdicFieldGenericElement = sage.rings.padics.padic_field_generic_element.pAdicFieldGenericElement
pAdicLazyElement = sage.rings.padics.padic_lazy_element.pAdicLazyElement
#pAdicFieldCappedRelativeElement = sage.rings.padics.padic_ring_capped_relative_element
pari = sage.libs.pari.gen.pari
pari_gen = sage.libs.pari.gen.gen
PariError = sage.libs.pari.gen.PariError
PrecisionError = sage.rings.padics.precision_error.PrecisionError

class pAdicRingCappedAbsoluteElement(pAdicRingFixedModElement):

    def __init__(self, parent, x, absprec=infinity, relprec = infinity, construct=False):
        sage.rings.commutative_ring_element.CommutativeRingElement.__init__(self,parent)
        if construct:
            (self._value, self._absprec) = x
            return
        if not absprec is infinity and not relprec is infinity:
            raise ValueError, "can only specify one of absprec and relprec"
        if relprec is infinity:
            if absprec > parent.precision_cap():
                absprec = parent.precision_cap()

        if isinstance(x, pAdicGenericElement) and x.valuation() < 0:
            raise ValueError, "element valuation cannot be negative."
        if isinstance(x, pAdicLazyElement):
            if relprec is infinity:
                absprec = min(absprec, parent.precision_cap())
                try:
                    x.set_precision_absolute(absprec)
                except PrecisionError:
                    pass
                self._absprec = min(absprec, x.precision_absolute())
            else:
                try:
                    x.set_precision_relative(relprec)
                except PrecisionError:
                    pass
                self._absprec = min(parent.precision_cap(), x.precision_absolute())
            self._value = Mod(x.residue(self._absprec), self.parent().prime_pow(self.parent().precision_cap()))
            return
        if isinstance(x, pAdicGenericElement):
            if relprec is infinity:
                self._absprec = min(x.precision_absolute(), absprec)
            else:
                self._absprec = min(x.precision_absolute(), x.valuation() + relprec)
            self._value = Mod(Mod(x.lift(), self.parent().prime_pow(self._absprec)), self.parent().prime_pow(self.parent().precision_cap()))
            return

        if isinstance(x, pari_gen):
            if x.type() == "t_PADIC":
                absprec = min(x.padicprec(parent.prime()), absprec)
                x = x.lift()
            if x.type() == "t_INT":
                x = Integer(x)
            elif x.type() == "t_FRAC":
                x = Rational(x)
            else:
                raise TypeError, "unsupported coercion from pari: only p-adics, integers and rationals allowed"

        if sage.rings.integer_mod.is_IntegerMod(x):
            k, p = pari(x.modulus()).ispower()
            if not k or p != parent.prime():
                raise TypeError, "cannot change primes in creating p-adic elements"
            if absprec is infinity and relprec is infinity: #this allows you to lift integer_mod elements to higher precision than they are actually defined.
                absprec = min(absprec, k)
            x = x.lift()

        #if sage.rings.finite_field_element.is_FiniteFieldElement(x):
        #    if x.parent().order() != parent.prime():
        #        raise TypeError, "can only create p-adic element out of finite field when order of field is p"
        #    prec = min(prec, 1)
        #    x = x.lift()

        #self._absprec = prec


        #Now use the code below to convert from integer or rational, so don't make the next line elif

        if isinstance(x, (int, long)):
            val = sage.rings.arith.valuation(x, parent.prime())
        elif isinstance(x, Integer):
            val = x.valuation(parent.prime())
        elif isinstance(x, Rational):
            val = x.valuation(parent.prime())
            if val < 0:
                raise ValueError, "p divides the denominator"
        if isinstance(x, (int, long, Integer, Rational)):
            self._absprec = min(absprec, val + relprec, parent.precision_cap())
            self._value = Mod(Mod(x, parent.prime_pow(self._absprec)), parent.prime_pow(parent.precision_cap()))
        else:
            raise TypeError, "unable to create p-adic element"

    def __invert__(self):
        return self.parent().fraction_field()(self).__invert__()

    def _neg_(self):
        """
        Returns -x.

        INPUT:
        x -- a p-adic capped absolute element
        OUTPUT:
        -x

        EXAMPLES:
        sage: R = Zp(5, prec=10, type='capped-abs')
        sage: a = R(1)
        sage: -a
        4 + 4*5 + 4*5^2 + 4*5^3 + 4*5^4 + 4*5^5 + 4*5^6 + 4*5^7 + 4*5^8 + 4*5^9 + O(5^10)
        """
        return pAdicRingCappedAbsoluteElement(self.parent(), (-self._value, self._absprec), construct = True)

    def __pow__(self, right):
        new = Integer(right) #Need to make sure that this works for p-adic exponents
        val = self.valuation()
        if (val > 0) and (isinstance(right, pAdicRingGenericElement) or isinstance(right, pAdicFieldGenericElement)):
            raise ValueError, "Can only have p-adic exponent if base is a unit"
        return pAdicRingCappedAbsoluteElement(self.parent(), (self._value**new, min(self._absprec - val + new * val, self.parent().precision_cap())), construct = True)

    def _add_(self, right):
        return pAdicRingCappedAbsoluteElement(self.parent(), (self._value + right._value, min(self.precision_absolute(), right.precision_absolute())), construct = True)

    def _div_(self, right):
        return self * right.__invert__()

    #def __floordiv__(self, right):
    #    #There is still a bug in here
    #    if isinstance(right, Integer):
    #        right = pAdicRingCappedAbsoluteElement(self.parent(), right)
    #    rval = right.valuation()
    #    ppow = self.parent().prime_pow(rval)
    #    runit = right._unit_part()
    #    prec = min(self.precision_absolute(), right.precision_relative() + self.valuation())
    #    #print "self = %s, right = %s, rval = %s, runit = %s, prec = %s"%(self, right, rval, runit, prec)
    #    quotient = Mod(self._value.lift(), self.parent().prime_pow(prec)) / Mod(runit, self.parent().prime_pow(prec))
    #    return pAdicRingCappedAbsoluteElement(self.parent(), (Mod(quotient.lift() // ppow, self.parent().prime_pow(self.parent().precision_cap())), min(self.precision_relative(), right.precision_relative()) + self.valuation() - rval), construct = True)

    #def __mod__(self, right):
    #    rval = right.valuation()
    #    ppow = self.parent().prime_pow(rval)
    #    runit = right.unit_part()
    #    quotient = self / runit
    #    return pAdicRingCappedAbsoluteElement(self.parent(), (Mod(quotient.lift() % ppow, self.parent().prime_pow(self.parent().precision_cap())), self.parent().precision_cap()), construct = True)

    def __lshift__(self, shift):
        """
        EXAMPLES:
        We create a capped relative field:
            sage: R = Zp(5, 20, 'capped-rel'); a = R(1000); a
            3*5^3 + 5^4 + O(5^23)

        Shifting to the right is the same as dividing by a power of
        the uniformizer $p$ of the $p$-adic ring.
            sage: a >> 1
            3*5^2 + 5^3 + O(5^22)

        Shifting to the left is the same as multiplying by a power of $p$:
            sage: a << 2
            3*5^5 + 5^6 + O(5^25)
            sage: a*5^2
            3*5^5 + 5^6 + O(5^25)

        Shifting by a negative integer to the left is the same as right shifting
        by the absolute value:
            sage: a << -3
            3 + 5 + O(5^20)
            sage: a >> 3
            3 + 5 + O(5^20)
        """
        shift = Integer(shift)
        if shift < 0:
            return self.__rshift__(-shift)
        return pAdicRingCappedAbsoluteElement(self.parent(),
                                              (Mod(self._value.lift() *self.parent().prime_pow(shift),
                                                   self.parent().prime_pow(self.parent().precision_cap())),
                                               self.precision_absolute() + shift), construct = True)

    def __rshift__(self, shift):
        """
        EXAMPLES:
            sage: R = Zp(997, 7, 'capped-rel'); a = R(123456878908); a
            964*997 + 572*997^2 + 124*997^3 + O(997^8)

        Shifting to the right divides by a power of p, but dropping terms with
        negative valuation:
            sage: a >> 3
            124 + O(997^5)

        Shifting to the left multiplies by that power of p.
            sage: a << 3
            964*997^4 + 572*997^5 + 124*997^6 + O(997^11)
        """
        shift = Integer(shift)
        if shift < 0:
            return self.__lshift__(-shift)
        if shift > self.precision_absolute():
            return pAdicRingCappedAbsoluteElement(self.parent(), (Mod(0, self.parent().prime_pow(self.parent().precision_cap())), Integer(0)), construct = True)
        return pAdicRingCappedAbsoluteElement(self.parent(), (Mod(self._value.lift() // self.parent().prime_pow(shift), self.parent().prime_pow(self.parent().precision_cap())), self.precision_absolute() - shift), construct = True)

    def _mul_(self, right):
        """
        EXAMPLES:

        """
        return pAdicRingCappedAbsoluteElement(self.parent(), (self._value * right._value, min(self.valuation() + right.valuation() + min(self.precision_relative(), right.precision_relative()), self.parent().precision_cap())), construct = True)

    def _sub_(self, right):
        return pAdicRingCappedAbsoluteElement(self.parent(), (self._value - right._value, min(self.precision_absolute(), right.precision_absolute())), construct = True)

    def add_bigoh(self, prec):
        """
        Returns a new element with absolute precision decreased to prec

        INPUT:
            self -- a p-adic element
            prec -- an integer
        OUTPUT:
            element -- self with precision set to the minimum of self's precision and prec

        EXAMPLES:
            sage: R = Zp(7,4,'capped-abs','series'); a = R(8); a.add_bigoh(1)
            1 + O(7)
        """
        return pAdicRingCappedAbsoluteElement(self.parent(), (Mod(Mod(self._value, self.parent().prime_pow(prec)), self.parent().prime_pow(self.parent().precision_cap())), min(prec, self.parent().precision_cap())), construct = True)

    def copy(self):
        return pAdicRingCappedAbsoluteElement(self.parent(), (self._value, self._absprec), construct = True)

    def exp_artin_hasse(self):
        raise NotImplementedError

    def gamma(self):
        raise NotImplementedError

    def lift(self):
        return self._value.lift() % self.parent().prime_pow(self.precision_absolute())

    def lift_to_precision(self, absprec):
        absprec = min(absprec, self.parent().precision_cap())
        return pAdicRingCappedAbsoluteElement(self.parent(), (Mod(self._value, self.parent().prime_pow(absprec)), absprec), construct = True)

    def list(self):
        """
        Returns a list of coeficiants of p starting with $p^0$
        INPUT:
            self -- a p-adic element
        OUTPUT:
            list -- the list of coeficients of self
        EXAMPLES:
            sage: R = Zp(7,4,'capped-abs'); a = R(2*7+7**2); a.list()
                [0, 2, 1, 0]

        NOTE:
            this differs from the list method of padic_field_element
            use slice operators to get a particular range

        """
        def plist(n, p):
            if n == 0:
                return []
            else:
                return [n % p] + plist(n // p, p)
        rep = plist(self._value.lift(), self.parent().prime())
        return rep + [0 for w in range(len(rep), self.precision_absolute())]

    def log_artin_hasse(self):
        raise NotImplementedError

    def multiplicative_order(self):
        r"""
        Returns the multiplicative order of self, where self is considered to be one if it is one modulo $p^{\mbox{prec}}$.

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
        else:
            return Mod(self._value, self.parent().prime_pow(self.precision_absolute())).multiplicative_order()

    def padded_list(self, n):
        """
        Returns a list of coeficiants of p starting with $p^0$ up to $p^n$ exclusive (padded with zeros if needed)
        INPUT:

self -- a p-adic element
            n - an integer
        OUTPUT:
            list -- the list of coeficients of self
        EXAMPLES:
            sage: R = Zp(7,4,'capped-abs'); a = R(2*7+7**2); a.padded_list(5)
                [0, 2, 1, 0, 0]

        NOTE:
            this differs from the padded_list method of padic_field_element
            the slice operators throw an error if asked for a slice above the precision, while this function works
        """
        return self.list()[:n] + [0 for w in range(self.precision_absolute(), n)]

    def precision_absolute(self):
        """
        Returns the absolute precision of self
         INPUT:
            self -- a p-adic element
        OUTPUT:
            integer -- the absolute precision of self
        EXAMPLES:
            sage: R = Zp(7,4,'capped-abs'); a = R(7); a.precision_absolute()
                4
       """
        return self._absprec

    def precision_relative(self):
        """
        Returns the relative precision of self
         INPUT:
            self -- a p-adic element
        OUTPUT:
            integer -- the relative precision of self
        EXAMPLES:
            sage: R = Zp(7,4,'capped-abs'); a = R(7); a.precision_relative()
                3
       """
        return self._absprec - self.valuation()

    def square_root(self):
        r"""
        Returns the square root of this p-adic number

        IMPORTANT NOTE!!!!!!:
            doctests for square_root currently disabled until trac 268 resolved

        INPUT:
            self -- a p-adic element
        OUTPUT:
            p-adic element -- the square root of this p-adic number
        EXAMPLES:
            sage.: R = Zp(3,20,'capped-abs')
            sage.: R(0).square_root()
                0 + O(3^20)
            sage.: R(1).square_root()
                1 + O(3^20)
            sage.: R(4).square_root() == R(-2)
                True
            sage.: R(9).square_root()
                3 + O(3^20)
            sage.: R2 = Zp(2,20,'capped-abs')
            sage.: R2(0).square_root()
                0 + O(2^20)
            sage.: R2(1).square_root()
                1 + O(2^19)
            sage.: R2(4).square_root()
                2 + O(2^20)
            sage.: R2(9).square_root() == R2(3) or R2(9).square_root() == R2(-3)
                True
            sage.: R2(17).square_root()
                206569 + O(2^19)
            sage.: R3 = Zp(5,20,'capped-abs')
            sage.: R3(0).square_root()
                0
            sage.: R3(1).square_root()
                1 + O(5^20)
            sage.: R3(-1).square_root() == R3.teichmuller(2) or R3(-1).square_root() == R3.teichmuller(3)
                True
        """
        if self.is_square():
            return pAdicRingCappedAbsoluteElement(self.parent(), (self._value.square_root(), self.precision_absolute() - self.valuation() // 2), construct = True)
        else:
            raise ValueError, "element is not a square" # should eventually change to return an element of an extension field

    def unit_part(self):
        r"""
        Returns the unit part of self.

        INPUT:
            self -- a p-adic element
        OUTPUT:
            p-adic element -- the unit part of self
        EXAMPLES:
            sage: R = Zp(17,4,'capped-abs', 'val-unit')
            sage: a = R(18*17)
            sage: a.unit_part()
                18 + O(17^3)
            sage: type(a)
                <class 'sage.rings.padics.padic_ring_capped_absolute_element.pAdicRingCappedAbsoluteElement'>
        """
        return pAdicRingCappedAbsoluteElement(self.parent(), (self._unit_part(), self.precision_relative()), construct = True)
