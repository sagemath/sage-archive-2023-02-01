"""
Elements of p-Adic Fields with Capped Relative Precision

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
import sage.rings.padics.padic_field_generic_element
import sage.rings.padics.padic_ring_generic_element
import sage.rings.padics.padic_lazy_element
#import sage.rings.padics.padic_ring_capped_relative_element
#import sage.rings.commutative_ring_element
#import sage.rings.finite_field_element
#import sage.rings.integer
#import sage.rings.rational
#import sage.rings.padics.precision_error
#import sage.rings.infinity

#pAdicRingCappedRelativeElement = sage.rings.padics.padic_ring_capped_relative_element.pAdicRingCappedRelativeElement
pAdicLazyElement = sage.rings.padics.padic_lazy_element.pAdicLazyElement
pAdicFieldGenericElement = sage.rings.padics.padic_field_generic_element.pAdicFieldGenericElement
pAdicGenericElement = sage.rings.padics.padic_generic_element.pAdicGenericElement
PrecisionError = sage.rings.padics.precision_error.PrecisionError
pari = sage.libs.pari.gen.pari
pari_gen = sage.libs.pari.gen.gen
PariError = sage.libs.pari.gen.PariError
#Zp = sage.rings.padics.padic_ring.Zp
Mod = sage.rings.integer_mod.Mod
Integer = sage.rings.integer.Integer
Rational = sage.rings.rational.Rational
infinity = sage.rings.infinity.infinity
QQ = sage.rings.rational_field.QQ

class pAdicFieldCappedRelativeElement(sage.rings.padics.padic_field_generic_element.pAdicFieldGenericElement):
    def __init__(self, parent, x, prec=None, construct=False):
        sage.rings.commutative_ring_element.CommutativeRingElement.__init__(self,parent)
        if construct:
            (self._ordp, self._unit, self._relprec) = x
            return
        if prec == None or prec > parent.precision_cap():
            prec = parent.precision_cap()
        if isinstance(x, pAdicLazyElement):
            try:
                x.set_precision_relative(prec)
            except PrecisionError:
                pass
        if isinstance(x, pAdicGenericElement):
            if parent.prime() != x.parent().prime():
                raise ValueError, "Cannot coerce between p-adic rings with different primes."
            self._ordp = x.valuation()
            self._relprec = min(x.precision_relative(), prec)
            self._unit = Mod(x._unit_part(), x.parent().prime_pow(self._relprec))
            return

        if isinstance(x, pari_gen) and x.type() == "t_PADIC":
            t = x.lift()
            prec = min(x.padicprec(parent.prime()) - x.valuation(parent.prime()), prec)
            if t.type() == 't_INT':
                x = int(t)
            else:
                x = QQ(t)

	if isinstance(x, pari_gen) and x.type() == "t_INT":
	    x = int(x)

        if sage.rings.finite_field_element.is_FiniteFieldElement(x):
            if x.parent().order() != parent.prime():
                raise TypeError, "can only create p-adic element out of finite field when order of field is p"
            #prec = min(prec, 1)
            x = x.lift()

        elif sage.rings.integer_mod.is_IntegerMod(x):
            k, p = pari(x.modulus()).ispower()
            if not k or p != parent.prime():
                raise TypeError, "cannot change primes in creating p-adic elements"
            x = x.lift()
            k = k - x.valuation(p)
            prec = min(prec, k)

            # We now use the code, below, so don't make the next line elif
        if isinstance(x, (int, long)):
            self._ordp = sage.rings.arith.valuation(x, parent.prime())
        elif isinstance(x, (Integer, Rational)):
            self._ordp = x.valuation(self.parent().prime())
        else:
            raise TypeError, "unable to compute ordp"
        if self._ordp == infinity:
            self._unit = Mod(0, 1)
            self._relprec = prec
            return
        x = x / self.parent().prime_pow(self._ordp)
        self._unit = Mod(x, self.parent().prime_pow(prec))
        self._relprec = prec
        return

    def _repr_(self, mode = None, do_latex = False):
        return sage.rings.padics.padic_generic_element.pAdicGenericElement._repr_(self, mode, do_latex, True)

    def __invert__(self, prec=infinity):
        r"""
        Returns the multiplicative inverse of self.

        EXAMPLES:
            sage: R = Qp(7,4,'capped-rel','series'); a = R(3); a
            3 + O(7^4)
            sage: ~a
            5 + 4*7 + 4*7^2 + 4*7^3 + O(7^4)
        """
        if self.precision_relative() == 0:
            raise PrecisionError, "cannot divide by something indistinguishable from zero."
        prec = min(prec, self.precision_relative())
        ppow = self.parent().prime_pow(prec)
        unit = Mod(sage.rings.arith.inverse_mod(self._unit_part(), ppow), ppow)
        return pAdicFieldCappedRelativeElement(self.parent(), (-self.valuation(), unit, prec), construct = True)

    def _neg_(self):
        return pAdicFieldCappedRelativeElement(self.parent(), (self._ordp, -self._unit, self._relprec), construct=True)

    def __pow__(self, right):
        """
        EXAMPLES:
            sage: K = Qp(19, 5, 'capped-rel','series')
            sage: a = K(-1); a
            18 + 18*19 + 18*19^2 + 18*19^3 + 18*19^4 + O(19^5)
            sage: a^2
            1 + O(19^5)
            sage: a^3
            18 + 18*19 + 18*19^2 + 18*19^3 + 18*19^4 + O(19^5)
            sage: K(5)^30
            11 + 14*19 + 19^2 + 7*19^3 + O(19^5)
        """
        right = Integer(right)
        if self == 0:
            if right == 0:
                raise ValueError, "0^0 not defined"
            return pAdicFieldCappedRelativeElement(self.parent(), (self.valuation() * right, Mod(0, 1), 0), construct = True) #this isn't quite right.
        if right < 0:
            inv = 1/self
            return inv**(-right)
        if right == 0:
            return pAdicFieldCappedRelativeElement(self.parent(), (0, Mod(1, self.parent().prime_pow(self.parent().precision_cap())), self.parent().precision_cap()), construct = True)
        ordp = right * self.valuation()
        prec = self.precision_relative()
        unit = Mod(sage.rings.arith.power_mod(self._unit_part().lift(), right, self.parent().prime_pow(prec)),self.parent().prime_pow(prec))
        return pAdicFieldCappedRelativeElement(self.parent(), (ordp, unit, prec), construct = True)

    def _add_(self, right):
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
        return pAdicFieldCappedRelativeElement(self.parent(), (min(self.valuation(), right.valuation()) + v, u, rprec), construct = True)

    def __floordiv__(self, right):
        return self / right

    def _integer_(self):
        return self._unit.lift() * self.parent().prime_pow(self.valuation())

    def _mul_(self, right):
        rprec = min(self._relprec, right._relprec)
        return pAdicFieldCappedRelativeElement(self.parent(), (self.valuation() + right.valuation(), Mod(self._unit, self.parent().prime_pow(rprec)) * Mod(right._unit, self.parent().prime_pow(rprec)), rprec), construct = True)

    def add_bigoh(self, prec):
        """
        Returns a new element with absolute precision decreased to prec
        INPUT:
            self -- a p-adic element
            prec -- an integer
        OUTPUT:
            element -- self with precision set to the minimum of  self's precision and prec

        EXAMPLES:
            sage: R = Qp(7,4); R.set_print_mode('series'); a = R(8); a.add_bigoh(1)
                1 + O(7)
            sage: b = R(0); b.add_bigoh(3)
                O(7^3)
        """
        rprec = min(self._relprec, prec - self.valuation())
        if rprec <= 0:
            return pAdicFieldCappedRelativeElement(self.parent(), (prec, Mod(0, 1), 0), construct = True)
        return pAdicFieldCappedRelativeElement(self.parent(), (self.valuation(), Mod(self._unit, self.parent().prime_pow(rprec)), rprec), construct = True)

    def copy(self):
        return pAdicFieldCappedRelativeElement(self.parent(), (self.valuation(), self._unit, self._relprec), construct = True)

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
        Return a rational "congruent" to self modulo self's precision

        INPUT:
            self -- a p-adic element
        OUTPUT:
            rational -- a rational equal to self upto our precision
        EXAMPLES:
            sage: R = Qp(7,4); a = R(8); a.lift()
                8
            sage: R = Qp(7,4); a = R(8/7); a.lift()
                8/7
        """
        if self.valuation() == infinity:
            return 0
        else:
            return self.parent().prime_pow(self.valuation()) * self._unit_part().lift()

    def list(self):
        """
        Returns a list of coeficiants of p starting with the lowest power of p with non-zero coefficient
        INPUT:
            self -- a p-adic element
        OUTPUT:
            list -- the list of coeficients of self
        EXAMPLES:
            sage: R = Qp(7,4); a = R(2*7+7**2); a.list()
                [2, 1, 0, 0]

        NOTE:
            this differs from the list method of padic_ring_element
        """
        if (self.valuation() == infinity) or (self.precision_relative() <= 0):
            return []
        else:
            def plist(n, p, prec):
                if prec == 0:
                    return []
                else:
                    return [n % p] + plist(n // p, p, prec - 1)
            return plist(self._unit.lift(), self.parent().prime(), self.precision_relative())

    def padded_list(self, n):
        """
        Returns a list of coeficiants of p starting with the lowest power with non-zero coeficient up to $p^n$ exclusive (padded with zeros if needed)
        INPUT:
            self -- a p-adic element
            n - an integer
        OUTPUT:
            list -- the list of coeficients of self
        EXAMPLES:
            sage: R = Qp(7,3); a = R(2*7+7**2); a.padded_list(5)
                [2, 1, 0, 0]
            sage: a.padded_list(3)
                [2, 1]

        NOTE:
            this differs from the padded_list method of padic_ring_element
            the slice operators throw an error if asked for a slice above the precision

        """
        return self.list()[:(n - self.valuation())] + [0 for w in range(self.precision_relative(), (n - self.valuation()))]

    def precision_absolute(self):
        """
        Returns the absolute precision of self
         INPUT:
            self -- a p-adic element
        OUTPUT:
            integer -- the absolute precision of self
        EXAMPLES:
            sage: R = Qp(7,3); a = R(7); a.precision_absolute()
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
            sage: R = Qp(7,3); a = R(7); a.precision_relative()
                3
       """
        return self._relprec

    def residue(self, prec):
        """
        Reduces this mod $p^prec$
        INPUT:
            self -- a p-adic element that must be integral
            prec - an integer
        OUTPUT:
            element of Z/(p^prec Z) -- self reduced mod p^prec
        EXAMPLES:
            sage: R = Qp(7,4,'capped-rel'); a = R(8); a.residue(1)
                1
        """
        if prec > self.precision_absolute():
            raise PrecisionError, "Not enough precision known in order to compute residue."
        elif self.valuation() >= 0:
            return Mod(self.parent().prime_pow(self.valuation()) * self._unit.lift(), self.parent().prime_pow(prec))
        else:
            raise ValueError, "Element must have non-negative valuation in order to compute residue."

    def unit_part(self):
        r"""
        Returns the unit part of self.

        INPUT:
            self -- a p-adic element
        OUTPUT:
            p-adic element -- the unit part of self
        EXAMPLES:
            sage: R = Qp(17,4,'capped-rel')
            sage: a = R(18*17)
            sage: a.unit_part()
                1 + 17 + O(17^4)
            sage: type(a)
                <class 'sage.rings.padics.padic_field_capped_relative_element.pAdicFieldCappedRelativeElement'>
        """
        if self.precision_relative() == 0:
            raise PrecisionError, "Not enough precision to compute unit part."
        return pAdicFieldCappedRelativeElement(self.parent(), (0, self._unit, self._relprec), construct = True)

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
            sage: R = Qp(17, 4)
            sage: a = R(2*17^2)
            sage: a.valuation()
                2
            sage: R = Qp(5, 4)
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
            sage: R(1/2).valuation()
                0
            sage: R(1/5).valuation()
                -1
            sage: R(1/10).valuation()
                -1
            sage: R(1/25).valuation()
                -2
            sage: R(1/50).valuation()
                -2
        """
        return self._ordp
