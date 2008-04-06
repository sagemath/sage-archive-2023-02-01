"""
p-Adic Rings with capped absolute precision

AUTHOR:
    -- David Roe
    -- Genya Zaytman: documentation
    -- David Harvey: doctests

EXAMPLES:
p-Adic Rings are examples of inexact structures, as the reals are.
That means that elements cannot generally be stored exactly: to do so
would take an infinite amount of storage.  Instead, we store an
approximation to the elements with varying precision.

There are two types of precision for a p-adic element.  The first is
relative precision, which gives the number of known p-adic digits.
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

In the capped absolute type, instead of having a cap on the relative precision of an element there is instead a cap on the absolute precision.  Elements still store their own precisions, and as with the capped relative case, exact elements are truncated when cast into the ring.
    sage: R = Zp(5, 5, 'capped-abs', 'series'); a = R(4005); a
    5 + 2*5^3 + 5^4 + O(5^5)
    sage: b = R(4025); b
    5^2 + 2*5^3 + 5^4 + O(5^5)
    sage: a * b
    5^3 + 2*5^4 + O(5^5)
    sage: (a * b) // 5^3
    1 + 2*5 + O(5^2)
    sage: type((a * b) // 5^3)
    <type 'sage.rings.padics.padic_capped_absolute_element.pAdicCappedAbsoluteElement'>
    sage: (a * b) / 5^3
    1 + 2*5 + O(5^2)
    sage: type((a * b) / 5^3)
    <type 'sage.rings.padics.padic_capped_relative_element.pAdicCappedRelativeElement'>

p-Adic rings should be created using the creation function Zp as above.  This will ensure that there is only one instance of $\Z_p$ of a given type, p and precision.  It also saves typing very long class names.
    sage: Zp(17,10,'capped-rel')
    17-adic Ring with capped relative precision 10

    #sage: Zp(7, prec = 30, type = 'lazy', print_mode = 'val-unit')
    #Lazy 7-adic Ring
    sage: R = Zp(7, prec = 20, type = 'capped-rel', print_mode = 'val-unit'); S = Zp(7, prec = 20, type = 'capped-rel', print_mode = 'val-unit'); R is S
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

from sage.rings.integer_ring import ZZ
import sage.rings.integer
#import sage.rings.integer_mod
#import sage.rings.padics.padic_ring_generic
#import sage.rings.padics.padic_field_generic
import padic_ring_base_generic
import padic_capped_absolute_ring_generic
#import sage.rings.padics.padic_ring_fixed_mod
from sage.rings.padics.padic_capped_absolute_element import pAdicCappedAbsoluteElement
#import sage.rings.padics.padic_ring_capped_relative_element
#import sage.rings.padics.padic_lazy_element
#import sage.rings.infinity
import padic_base_generic_element

Integer = sage.rings.integer.Integer
#Mod = sage.rings.integer_mod.Mod
#infinity = sage.rings.infinity.infinity
pAdicRingBaseGeneric = padic_ring_base_generic.pAdicRingBaseGeneric
pAdicBaseGenericElement = padic_base_generic_element.pAdicBaseGenericElement
#pAdicFieldBaseGeneric = sage.rings.padics.padic_field_base_generic.pAdicFieldBaseGeneric
#pAdicRingFixedMod = sage.rings.padics.padic_ring_fixed_mod.pAdicRingFixedMod
#pAdicCappedAbsoluteElement = sage.rings.padics.padic_ring_capped_absolute_element.pAdicCappedAbsoluteElement
#pAdicRingCappedRelativeElement = sage.rings.padics.padic_ring_capped_relative_element.pAdicRingCappedRelativeElement
#pAdicLazyElement = sage.rings.padics.padic_lazy_element.pAdicLazyElement
pAdicCappedAbsoluteRingGeneric = padic_capped_absolute_ring_generic.pAdicCappedAbsoluteRingGeneric

class pAdicRingCappedAbsolute(pAdicRingBaseGeneric, pAdicCappedAbsoluteRingGeneric):
    r"""
    An implementation of the p-adic integers with capped absolute precision.
    """
    def __init__(self, p, prec, print_mode, names):
        pAdicRingBaseGeneric.__init__(self, p, prec, print_mode, names, pAdicCappedAbsoluteElement)

    def _repr_(self, do_latex = False):
        return "%s-adic Ring with capped absolute precision %s"%(self.prime(), self.precision_cap())

    def __contains__(self, x):
        if isinstance(x, (int, long, Integer)):
            return True
        if not isinstance(x, pAdicBaseGenericElement):
            return False
        par = x.parent()
        if par.is_field():
            return False
        if par.prime() != self.prime():
            return False
        if par.is_capped_absolute() or par.is_capped_relative() or par.is_lazy():
            if par.is_capped_absolute() and par.precision_cap() < self.precision_cap():
                return False
            return True
        return False

    def fraction_field(self):
        r"""
        Returns the fraction field of self.
        """
        from sage.rings.padics.factory import Qp
        return Qp(self.prime(), self.precision_cap(), 'capped-rel', self.print_mode())

    def random_element(self):
        """
        Returns a random element of self.
        """
        return self(ZZ.random_element(self.prime()**self.precision_cap()))

