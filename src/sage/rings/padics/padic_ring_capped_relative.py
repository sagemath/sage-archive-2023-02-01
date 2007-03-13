"""
p-Adic Rings with capped relative precision

AUTHOR:
    -- David Roe
    -- Genya Zaytman: documentation
    -- David Harvey: doctests

EXAMPLES:
p-Adic Rings are examples of inexact structures, as the reals are.  That means that elements cannot generally be stored exactly: to do so would take an infinite amount of storage.  Instead, we store an approximation to the elements with varying precision.

There are two types of precision for a p-adic element.  The first is relative precision, which gives the number of known p-adic digits.
    sage: R = Zp(5, 20, 'capped-rel', 'series'); a = R(675); a
        2*5^2 + 5^4 + O(5^22)
    sage: a.precision_relative()
        20

The second type of precision is absolute precision, which gives the power of p that this element is stored modulo.
    sage: a.precision_absolute()
        22

The number of times that p divides the element is called the valuation, and can be accessed with the functions valuation() and ordp()
    sage: a.valuation()
        2

The following relationship holds: self.valuation() + self.precision_relative() == self.precision_absolute().
    sage: a.valuation() + a.precision_relative() == a.precision_absolute()
        True

In the capped relative case, the relative precision of an element is restricted to be at most a certain value, specified at the creation of the field.  Individual elements also store their own precision, so the effect of various arithmetic operations on precision is tracked.  When you cast an exact element into a capped relative field, it truncates it to the precision cap of the field.
    sage: R = Zp(5, 5, 'capped-rel', 'series'); a = R(4006); a
        1 + 5 + 2*5^3 + 5^4 + O(5^5)
    sage: b = R(17/3); b
        4 + 2*5 + 3*5^2 + 5^3 + 3*5^4 + O(5^5)
    sage: c = R(4025); c
        5^2 + 2*5^3 + 5^4 + 5^5 + O(5^7)
    sage: a + b
	4*5 + 3*5^2 + 3*5^3 + 4*5^4 + O(5^5)
    sage: a + b + c
        4*5 + 4*5^2 + 5^4 + O(5^5)

p-Adic rings should be created using the creation function Zp as above.  This will ensure that there is only one instance of $\Z_p$ of a given type, p and precision.  It also saves typing very long class names.
    sage: Zp(17,10,'capped-rel')
        17-adic Ring with capped relative precision 10
    sage: Zp(7, prec = 30, type = 'lazy', print_mode = 'val-unit')
        Lazy 7-adic Ring
    sage: R = Zp(7, prec = 20, type = 'capped-rel', print_mode = 'val-unit'); S = Zp(7, prec = 20, type = 'capped-rel', print_mode = 'series'); R is S
        True
    sage: Zp(2)
        2-adic Ring with capped relative precision 20

Once one has a p-Adic ring, one can cast elements into it in the standard way.  Integers, ints, longs, Rationals, other p-Adic types, pari p-adics and elements of $\Z / p^n \Z$ can all be cast into a p-Adic ring.
    sage: R = Zp(5, 5, 'capped-rel','series'); a = R(16); a
        1 + 3*5 + O(5^5)
    sage: b = R(25/3); b
        2*5^2 + 3*5^3 + 5^4 + 3*5^5 + 5^6 + O(5^7)
    sage: S = Zp(5, 5, 'fixed-mod','val-unit'); c = S(Mod(75,125)); c
        5^2 * 3 + O(5^5)
    sage: R(c)
        3*5^2 + O(5^5)

Note that in the last example, since fixed-mod elements don't keep track of their precision, we assume that it has the full precision of the ring.  This is why you have to cast manually there.

While you can cast explicitly as above, the chains of automatic coercion are more restricted.  As always in SAGE, the following arrows are transitive and the diagram is commutative.

int -> long -> Integer -> Rational -> Zp lazy -> pari p-adic -> Zp capped-rel -> Zp capped-abs
                                       /    \                          \          /
                                      /      \                          \        /
                                     V        V                          V      V
                           Zp fixed-mod     Qp lazy ----------------> Qp capped-rel
                                     \_______________________          /
                                                             \        /
                                                              V      V
                                                             IntegerMod

In addition, there are arrows within each type from higher precision_cap to lower.
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
#import sage.rings.integer_mod
import sage.rings.padics.padic_ring_generic
import sage.rings.padics.padic_field_generic
import sage.rings.padics.padic_ring_capped_relative_element
import sage.rings.padics.padic_ring_fixed_mod
import sage.rings.padics.padic_lazy_element
import sage.rings.infinity

from sage.rings.integer_ring import ZZ

Integer = sage.rings.integer.Integer
Mod = sage.rings.integer_mod.Mod
infinity = sage.rings.infinity.infinity
pAdicRingBaseGeneric = sage.rings.padics.padic_ring_generic.pAdicRingBaseGeneric
pAdicFieldBaseGeneric = sage.rings.padics.padic_field_generic.pAdicFieldBaseGeneric
pAdicRingFixedMod = sage.rings.padics.padic_ring_fixed_mod.pAdicRingFixedMod
pAdicRingCappedRelativeElement = sage.rings.padics.padic_ring_capped_relative_element.pAdicRingCappedRelativeElement
pAdicLazyElement = sage.rings.padics.padic_lazy_element.pAdicLazyElement

class pAdicRingCappedRelative(pAdicRingBaseGeneric):
    r"""
    An implementation of the p-adic integers with capped relative precision.
    """

    def __call__(self, x, absprec = infinity, relprec = infinity):
        r"""
            Casts x into self.  Uses the constructor from pAdicRingCappedRelativeElement.
        """
        return pAdicRingCappedRelativeElement(self, x, absprec, relprec)

    def __cmp__(self, other):
        if isinstance(other, pAdicRingCappedRelative):
            if self.prime() < other.prime():
                return -1
            elif self.prime() > other.prime():
                return 1
            elif self.precision_cap() < other.precision_cap():
                return -1
            elif self.precision_cap() > other.precision_cap():
                return 1
            else:
                return 0
        elif isinstance(other, pAdicRingFixedMod) or isinstance(other, pAdicFieldBaseGeneric):
            return 1
        else:
            return -1

    def __contains__(self, x):
        if isinstance(x, (int, long, Integer)):
            return True
        if isinstance(x, pAdicRingCappedRelativeElement) and x.parent().prime() == self.prime() and x.parent.precision_cap() >= self.precision_cap():
            return True
        if isinstance(x, pAdicLazyElement) and not x.parent().is_field() and x.parent().prime() == self.prime():
            return True
        return False

    def _coerce_impl(self, x):
        if self.__contains__(x):
            return pAdicRingCappedRelativeElement(self, x)
        else:
            raise TypeError, "no canonical coercion of x"

    def _repr_(self, do_latex=False):
        return "%s-adic Ring with capped relative precision %s"%(self.prime(), self.precision_cap())

    def teichmuller(self, x, prec = None):
        r"""
        Returns the teichmuller representative of x.

        INPUT:
            self -- a p-adic ring
            x -- an integer or element of $\Z / p\Z$ that is not divisible by $p$
        OUTPUT:
            element -- the teichmuller lift of x
        EXAMPLES:
            sage: R = Zp(5, 10, 'capped-rel', 'series')
            sage: R.teichmuller(2)
                2 + 5 + 2*5^2 + 5^3 + 3*5^4 + 4*5^5 + 2*5^6 + 3*5^7 + 3*5^9 + O(5^10)
        """
        if prec is None:
            prec = self.precision_cap()
        p = self.prime()
        x = Mod(x,self.prime_pow(prec))
        xnew = x**p
        while x != xnew:
            x = xnew
            xnew = x**p
        return pAdicRingCappedRelativeElement(self, x)

    def integer_ring(self):
        r"""
        Returns the integer ring of self, i.e. self.
        """
        return self

    def fraction_field(self):
        r"""
        Returns the fraction field of self.

        INPUT:
            self -- a p-adic ring
        OUTPUT:
            the p-adic field that is the fraction field of this ring

        """
        from sage.rings.padics.qp import Qp
        return Qp(self.prime(), self.precision_cap(), 'capped-rel', self.print_mode())

    def random_element(self, algorithm='default'):
        r"""
        Returns a random element of self, optionally using the algorithm
        argument to decide how it generates the element. Algorithms currently
        implemented:

            default: Choose a_i, i >= 0, randomly between 0 and p-1 until
              a nonzero choice is made. Then continue choosing a_i randomly
              between 0 and p-1 until we reach precision_cap, and return
              \sum a_i p^i.
        """
        if (algorithm == 'default'):
            i = 0
            a_i = ZZ.random_element(self.prime())
            while a_i.is_zero():
                i += 1
                a_i = ZZ.random_element(self.prime())
            return self((self.prime()**i)*(a_i + self.prime()*ZZ.random_element(self.prime()**(self.precision_cap()-1))))
        else:
            raise NotImplementedError, "Don't know %s algorithm"%algorithm

    def unit_group(self):
        raise NotImplementedError

    def unit_group_gens(self):
        raise NotImplementedError

    def principal_unit_group(self):
        raise NotImplementedError

