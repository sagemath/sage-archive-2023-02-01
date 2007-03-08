"""
Elements of p-Adic Rings with Capped Relative Precision

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

#import sage.libs.all
#import sage.rings.arith
#import sage.rings.integer_mod
import sage.rings.padics.padic_ring_generic_element
import sage.rings.padics.padic_field_generic_element
import sage.rings.padics.padic_lazy_element
#import sage.rings.padics.padic_field_capped_relative_element
#import sage.rings.commutative_ring_element
#import sage.rings.finite_field_element
#import sage.rings.integer
#import sage.rings.rational
#import sage.rings.padics.precision_error
#import sage.rings.infinity

#pAdicFieldCappedRelativeElement = sage.rings.padics.padic_field_capped_relative_element.pAdicFieldCappedRelativeElement
pAdicRingGenericElement = sage.rings.padics.padic_ring_generic_element.pAdicRingGenericElement
pAdicFieldGenericElement = sage.rings.padics.padic_field_generic_element.pAdicFieldGenericElement
pAdicLazyElement = sage.rings.padics.padic_lazy_element.pAdicLazyElement
pAdicGenericElement = sage.rings.padics.padic_generic_element.pAdicGenericElement
PrecisionError = sage.rings.padics.precision_error.PrecisionError
pari = sage.libs.pari.gen.pari
pari_gen = sage.libs.pari.gen.gen
PariError = sage.libs.pari.gen.PariError
#Zp = sage.rings.padics.zp.Zp
Mod = sage.rings.integer_mod.Mod
Integer = sage.rings.integer.Integer
Rational = sage.rings.rational.Rational
infinity = sage.rings.infinity.infinity

class pAdicRingCappedRelativeElement(pAdicRingGenericElement):
    def __init__(self, parent, x, absprec=None, relprec=None, construct=False):
        """
        Constructs new element with given parent and value.

        INPUT:
            x -- value to coerce into ring
            prec -- number of digits of precision, or None to use the ring's
                precision cap
            construct -- boolean, default False. True is for internal use,
                in which case x is a triple to be assigned directly.

        EXAMPLES:
            sage: R = Zp(5, 10, 'capped-rel')

        Construct from integers:
            sage: R(3)
            3 + O(5^10)
            sage: R(75)
            3*5^2 + O(5^12)
            sage: R(0)
            0
            sage: R(-1)
            4 + 4*5 + 4*5^2 + 4*5^3 + 4*5^4 + 4*5^5 + 4*5^6 + 4*5^7 + 4*5^8 + 4*5^9 + O(5^10)
            sage: R(-5)
            4*5 + 4*5^2 + 4*5^3 + 4*5^4 + 4*5^5 + 4*5^6 + 4*5^7 + 4*5^8 + 4*5^9 + 4*5^10 + O(5^11)
            sage: R(-7*25)
            3*5^2 + 3*5^3 + 4*5^4 + 4*5^5 + 4*5^6 + 4*5^7 + 4*5^8 + 4*5^9 + 4*5^10 + 4*5^11 + O(5^12)

        Construct from rationals:
            sage: R(1/2)
            3 + 2*5 + 2*5^2 + 2*5^3 + 2*5^4 + 2*5^5 + 2*5^6 + 2*5^7 + 2*5^8 + 2*5^9 + O(5^10)
            sage: R(-7875/874)
            3*5^3 + 2*5^4 + 2*5^5 + 5^6 + 3*5^7 + 2*5^8 + 3*5^10 + 3*5^11 + 3*5^12 + O(5^13)
            sage: R(15/425)
            Traceback (most recent call last):
            ...
            ValueError: element not a p-adic integer.

        Construct from IntegerMod:
            sage: R(Integers(125)(3))
            3 + O(5^3)
            sage: R(Integers(5)(3))
            3 + O(5)
            sage: R(Integers(5^30)(3))
            3 + O(5^10)
            sage: R(Integers(5^30)(1+5^23))
            1 + O(5^10)
            sage: R(Integers(49)(3))
            Traceback (most recent call last):
            ...
            TypeError: cannot change primes in creating p-adic elements

        # todo: should the above TypeError be another type of error?

            sage: R(Integers(48)(3))
            Traceback (most recent call last):
            ...
            TypeError: cannot change primes in creating p-adic elements

        # todo: the error message for the above TypeError is not quite accurate

        Some other conversions:
            sage: R(R(5))
            5 + O(5^11)

        # todo: doctests for converting from other types of p-adic rings

        """
        sage.rings.commutative_ring_element.CommutativeRingElement.__init__(self, parent)
        if construct:
            (self._ordp, self._unit, self._relprec) = x
            return

        if not absprec is None and not relprec is None:
            raise ValueError, "can only specify one of absprec and relprec"
        if absprec is None:
            if relprec is None or relprec > parent.precision_cap():
                relprec = parent.precision_cap()

        if isinstance(x, pAdicGenericElement) and x.parent().is_field() and x.valuation() < 0:
            raise ValueError, "element has negative valuation."
        if isinstance(x, pAdicLazyElement):
            if relprec is None:
                try:
                    x.set_precision_absolute(absprec)
                except PrecisionError:
                    pass
                self._relprec = min(parent.precision_cap(), x.precision_relative())
            else:
                relprec = min(relprec, parent.precision_cap())
                try:
                    x.set_precision_relative(relprec)
                except PrecisionError:
                    pass
                self._relprec = min(relprec, x.precision_relative())
            self._ordp = x._min_valuation()
            self._unit = Mod(x._unit_part(), parent.prime_pow(self._relprec))
            return
        if isinstance(x, pAdicRingGenericElement) or isinstance(x, pAdicFieldGenericElement):
            if parent.prime() != x.parent().prime():
                raise ValueError, "Cannot coerce between p-adic rings with different primes."
            self._ordp = x.valuation()
            if relprec is None:
                relprec = absprec - self._ordp
            self._relprec = min(relprec, x.precision_relative(), parent.precision_cap())
            self._unit = Mod(x._unit_part(), parent.prime_pow(self._relprec))
            return

        if isinstance(x, pari_gen):
            if x.type() == "t_PADIC":
                if not absprec is None:
                    absprec = min(x.padicprec(parent.prime()), absprec)
                else:
                    absprec = x.padicprec(parent.prime())
                x = x.lift()
            if x.type() == "t_INT":
                x = Integer(x)
            elif x.type() == "t_FRAC":
                x = Rational(x)
            else:
                raise TypeError, "unsupported coercion from pari: only p-adics, integers and rationals allowed"

        #if sage.rings.finite_field_element.is_FiniteFieldElement(x):
        #    if x.parent().order() != parent.prime():
        #        raise TypeError, "can only create p-adic element out of finite field when order of field is p"
        #    #prec = min(prec, 1)
        #    x = x.lift()

        elif sage.rings.integer_mod.is_IntegerMod(x):
            k, p = pari(x.modulus()).ispower()
            if not k or p != parent.prime():
                raise TypeError, "cannot change primes in creating p-adic elements"
            x = x.lift()
            if absprec is None:
                absprec = k
            else:
                absprec = min(k, absprec)

            # We now use the code, below, so don't make the next line elif
        if isinstance(x, (int, long)):
            self._ordp = sage.rings.arith.valuation(x, parent.prime())
        elif isinstance(x, (Integer, Rational)):
            self._ordp = x.valuation(self.parent().prime())
        else:
            raise TypeError, "cannot create a p-adic out of %s"%(type(x))
        if self._ordp < 0:
            raise ValueError, "element not a p-adic integer."
        elif self._ordp == infinity:
            self._unit = Mod(0, 1)
            self._relprec = 0
            return
        x = x / self.parent().prime_pow(self._ordp)
        if relprec is None:
            self._relprec = min(absprec - self._ordp, parent.precision_cap())
        elif absprec is None:
            self._relprec = relprec
        else:
            self._relprec = min(relprec, absprec - self._ordp, parent.precision_cap())
        self._unit = Mod(x, self.parent().prime_pow(self._relprec))
        return

    def _repr_(self, mode = None, do_latex = False):
        return sage.rings.padics.padic_generic_element.pAdicGenericElement._repr_(self, mode, do_latex, True)

    def __mod__(self, right):
        val = self.valuation()
        rval = right.valuation()
        if rval > val:
            raise PrecisionError, "not enough precision to reduce"
        if right == self.parent().prime_pow(rval):
            return Mod(self.lift(), right)
        else:
            raise ValueError, "modulus must be a power of p"

    def _neg_(self):
        """
        EXAMPLES:
            sage: R = Zp(5, 20, 'capped-rel', 'val-unit')
            sage: -R(1)
            95367431640624 + O(5^20)
            sage: -R(5)
            5 * 95367431640624 + O(5^21)
            sage: -R(0)
            0
        """
        return pAdicRingCappedRelativeElement(self.parent(), (self._ordp, -self._unit, self._relprec), construct=True)

    def __pow__(self, right):
        """
        EXAMPLES:
            sage: R = Zp(19, 5, 'capped-rel','series')
            sage: a = R(-1); a
                18 + 18*19 + 18*19^2 + 18*19^3 + 18*19^4 + O(19^5)
            sage: a^2
                1 + O(19^5)
            sage: a^3
                18 + 18*19 + 18*19^2 + 18*19^3 + 18*19^4 + O(19^5)
            sage: R(5)^30
                11 + 14*19 + 19^2 + 7*19^3 + O(19^5)
        """
        right = Integer(right)
        if self == 0:
            if right == 0:
                raise ValueError, "0^0 not defined"
            return pAdicRingCappedRelativeElement(self.parent(), (self.valuation() * right, Mod(0, self.parent().prime()), 0), construct = True) #this isn't quite right.
        if right < 0:
            inv = 1/self
            return inv**(-right)
        if right == 0:
            return pAdicRingCappedRelativeElement(self.parent(), (0, Mod(1, self.parent().prime_pow(self.parent().precision_cap())), self.parent().precision_cap()), construct = True)
        ordp = right * self.valuation()
        prec = self.precision_relative()
        unit = Mod(sage.rings.arith.power_mod(self._unit_part().lift(), right, self.parent().prime_pow(prec)),self.parent().prime_pow(prec))
        return pAdicRingCappedRelativeElement(self.parent(), (ordp, unit, prec), construct = True)

    def _add_(self, right): #assumes both have the same parent (ie same p and same relative cap)
        if self.valuation() == infinity:
            return right
        if right.valuation() == infinity:
            return self
        if self.valuation() == right.valuation():
            if self._relprec == right._relprec:
                u = self._unit + right._unit
                rprec = self._relprec
            elif self._relprec < right._relprec:
                u = self._unit + Mod(right._unit, self.parent().prime_pow(self._relprec))
                rprec = self._relprec
            else:
                u = Mod(self._unit, self.parent().prime_pow(right._relprec)) + right._unit
                rprec = right._relprec
            ul = u.lift()
            if ul == 0: #In this case we set the valuation of the sum to be the minimum possible valuation and say that we have zero precision.
                v = min(self._relprec, right._relprec)
                rprec = 0
                u = Mod(0, 1)
            else:
                v = ul.valuation(self.parent().prime())
                if v > 0:
                    u = Mod(ul // self.parent().prime_pow(v), self.parent().prime_pow(rprec - v))
                    rprec = rprec - v
        else:
            rprec = min(min(self.valuation() + self.precision_relative(), right.valuation() + right.precision_relative()) - min(self.valuation(), right.valuation()), self.parent().precision_cap())
            uself = self._unit.lift() * self.parent().prime_pow(self.valuation() - min(self.valuation(), right.valuation()))
            uright = right._unit.lift() * self.parent().prime_pow(right.valuation() - min(self.valuation(), right.valuation()))
            u = Mod(uself + uright, self.parent().prime_pow(rprec))
            v = 0
        return pAdicRingCappedRelativeElement(self.parent(), (min(self.valuation(), right.valuation()) + v, u, rprec), construct = True)

    def __floordiv__(self, right):
        if isinstance(right, Integer):
            right = pAdicRingCappedRelativeElement(self.parent(), right)
        val = self.valuation()
        rval = right.valuation()
        if val >= rval:
            prec = min(self.precision_absolute() - rval, right.precision_relative())
            return pAdicRingCappedRelativeElement(self.parent(), (val - rval, Mod(self._unit.lift(), self.parent().prime_pow(prec)) / Mod(right._unit.lift(), self.parent().prime_pow(prec)), prec), construct = True)
        else:
            ppow = self.parent().prime_pow(rval - val)
            u = (self._unit / right._unit).lift() // ppow
            uval = u.valuation(self.parent().prime())
            prec = min(self.precision_absolute() - rval - uval, right.precision_relative())
            return pAdicRingCappedRelativeElement(self.parent(), (uval, Mod(u / self.parent().prime_pow(uval), self.parent().prime_pow(prec)), prec), construct = True)

    def __mod__(self, right):
        if isinstance(right, Integer):
            right = pAdicRingCappedRelativeElement(self.parent(), right)
        val = self.valuation()
        rval = right.valuation()
        if val >= rval:
            return pAdicRingCappedRelativeElement(self.parent(), 0)
        else:
            ppow = self.parent().prime_pow(rval - val)
            u = (self._unit / right._unit).lift() % ppow
            return pAdicRingCappedRelativeElement(self.parent(), (self.valuation(), Mod(u, self.parent().prime_pow(self.parent().precision_cap())), self.parent().precision_cap()), construct = True)

    def _integer_(self):
        return self._unit.lift() * self.parent().prime_pow(self.valuation())

    def _mul_(self, right):
        rprec = min(self._relprec, right._relprec)
        return pAdicRingCappedRelativeElement(self.parent(), (self.valuation() + right.valuation(), Mod(self._unit, self.parent().prime_pow(rprec)) * Mod(right._unit, self.parent().prime_pow(rprec)), rprec), construct = True)

    def add_bigoh(self, prec):
        """
        Returns a new element with absolute precision decreased to prec
        INPUT:
            self -- a p-adic element
            prec -- an integer
        OUTPUT:
            element -- self with precision set to the minimum of  self's precision and prec

        EXAMPLE:
            sage: R = Zp(7,4,'capped-rel','series'); a = R(8); a.add_bigoh(1)
                1 + O(7)
            sage: b = R(0); b.add_bigoh(3)
                O(7^3)
        """
        rprec = min(self._relprec, prec - self.valuation())
        if rprec <= 0:
            return pAdicRingCappedRelativeElement(self.parent(), (prec, Mod(0, 1), 0), construct = True)
        return pAdicRingCappedRelativeElement(self.parent(), (self.valuation(), Mod(self._unit, self.parent().prime_pow(rprec)), rprec), construct = True)

    def copy(self):
        return pAdicRingCappedRelativeElement(self.parent(), (self.valuation(), self._unit, self._relprec), construct = True)

    def exp_artin_hasse(self):
        raise NotImplementedError
        #E_p(x) = exp(x + x^p/p + x^(p^2)/p^2 + ...)

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
        return (self._relprec <= 0) or (self.valuation() >= prec)

    def is_equal_to(self, right, prec):
        r"""
        Returns whether self is equal to right modulo $p^{\mbox{prec}}$.

        INPUT:
            self -- a p-adic element
            right -- a p-addic element
            prec -- an integer
        OUTPUT:
            boolean -- whether self is equal to right

        """
        return (self - right).is_zero(prec)

    def lift(self):
        """
        Return an integer congruent to self modulo self's precision

        INPUT:
            self -- a p-adic element
        OUTPUT:
            integer -- a integer congruent to self mod $p^{\mbox{prec}}$
        EXAMPLE:
            sage: R = Zp(7,4,'capped-rel'); a = R(8); a.lift()
                8
        """
        if self.valuation() == infinity:
            return 0
        else:
            return self.parent().prime_pow(self.valuation()) * self._unit_part().lift()

    def list(self):
        """
        Returns a list of coeficiants of p starting with $p^0$
        INPUT:
            self -- a p-adic element
        OUTPUT:
            list -- the list of coeficients of self
        EXAMPLES:
            sage: R = Zp(7,3,'capped-rel'); a = R(2*7+7**2); a.list()
                [0, 2, 1, 0]

        NOTE:
            this differs from the list method of padic_field_element
            use slice operators to get a particular range

        """
        if (self.valuation() == infinity) or (self.precision_relative() <= 0):
            return []
        else:
            def plist(n, p, prec):
                if prec == 0:
                    return []
                else:
                    return [n % p] + plist(n // p, p, prec - 1)
            return [0 for w in range(self.valuation())] + plist(self._unit.lift(), self.parent().prime(), self.precision_relative())

    def log_artin_hasse(self):
        raise NotImplementedError

    def padded_list(self, n):
        """
        Returns a list of coeficiants of p starting with $p^0$ up to $p^n$ exclusive (padded with zeros if needed)
        INPUT:
            self -- a p-adic element
            n - an integer
        OUTPUT:
            list -- the list of coeficients of self
        EXAMPLES:
            sage: R = Zp(7,3,'capped-rel'); a = R(2*7+7**2); a.padded_list(5)
                [0, 2, 1, 0, 0]

        NOTE:
            this differs from the padded_list method of padic_field_element
            the slice operators throw an error if asked for a slice above the precision
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
            sage: R = Zp(7,3,'capped-rel'); a = R(7); a.precision_absolute()
                4
       """
        return self._ordp + self._relprec

    def precision_relative(self):
        """
        Returns the relative precision of self
         INPUT:
            self -- a p-adic element
        OUTPUT:
            integer -- the relative precision of self
        EXAMPLES:
            sage: R = Zp(7,3,'capped-rel'); a = R(7); a.precision_relative()
                3
       """
        return self._relprec

    def residue(self, prec):
        """
        Reduces this mod $p^prec$
        INPUT:
            self -- a p-adic element
            prec - an integer
        OUTPUT:
            element of Z/(p^prec Z) -- self reduced mod p^prec
        EXAMPLES:
            sage: R = Zp(7,4,'capped-rel'); a = R(8); a.residue(1)
                1
        """
        if prec <= self.precision_absolute():
            return Mod(self.parent().prime_pow(self.valuation()) * self._unit.lift(), self.parent().prime_pow(prec))
        else:
            raise PrecisionError, "not enough precision to compute residue"

    def unit_part(self):
        r"""
        Returns the unit part of self.

        INPUT:
            self -- a p-adic element
        OUTPUT:
            p-adic element -- the unit part of self
        EXAMPLES:
            sage: R = Zp(17,4,'capped-rel')
            sage: a = R(18*17)
            sage: a.unit_part()
                1 + 17 + O(17^4)
            sage: type(a)
                <class 'sage.rings.padics.padic_ring_capped_relative_element.pAdicRingCappedRelativeElement'>
        """
        if self.precision_relative() == 0:
            raise PrecisionError, "Not enough precision to determine unit part"
        return pAdicRingCappedRelativeElement(self.parent(), (0, self._unit, self._relprec), construct = True)

    def _unit_part(self):
        return self._unit

    def valuation(self):
        """
        Returns the valuation of self.

        INPUT:
            self -- a p-adic element
        OUTPUT:
            integer -- the valuation of self

        EXAMPLES:
            sage: R = Zp(17, 4,'capped-rel')
            sage: a = R(2*17^2)
            sage: a.valuation()
                2
            sage: R = Zp(5, 4,'capped-rel')
            sage: R(0).valuation()
                +Infinity
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
        return self._ordp
