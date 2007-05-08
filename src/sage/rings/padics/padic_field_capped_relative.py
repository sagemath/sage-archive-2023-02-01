"""
p-Adic Fields with Capped Relative Precision

AUTHOR:
    -- David Roe
    -- Genya Zaytman: documentation
    -- David Harvey: doctests
    -- William Stein: doctest updates

EXAMPLES:

p-Adic Fields are examples of inexact structures, as the reals are.
That means that elements cannot generally be stored exactly: to do so
would take an infinite amount of storage.  Instead, we store an
approximation to the elements with varying precision.

There are two types of precision for a p-adic element.  The first is
relative precision, which gives the number of known p-adic digits:

    sage: R = Qp(5, 20, 'capped-rel', 'series'); a = R(675); a
    2*5^2 + 5^4 + O(5^22)
    sage: a.precision_relative()
    20

The second type of precision is absolute precision, which gives the
power of p that this element is stored modulo:
    sage: a.precision_absolute()
    22

The number of times that p divides the element is called the
valuation, and can be accessed with the functions valuation() and
ordp():
    sage: a.valuation()
    2

The following relationship holds: self.valuation() + self.precision_relative() == self.precision_absolute().
    sage: a.valuation() + a.precision_relative() == a.precision_absolute()
    True

In the capped relative case, the relative precision of an element is
restricted to be at most a certain value, specified at the creation of
the field.  Individual elements also store their own precision, so the
effect of various arithmetic operations on precision is tracked.  When
you cast an exact element into a capped relative field, it truncates
it to the precision cap of the field.
    sage: R = Qp(5, 5, 'capped-rel', 'series'); a = R(4006); a
    1 + 5 + 2*5^3 + 5^4 + O(5^5)
    sage: b = R(17/3); b
    4 + 2*5 + 3*5^2 + 5^3 + 3*5^4 + O(5^5)
    sage: c = R(4025); c
    5^2 + 2*5^3 + 5^4 + 5^5 + O(5^7)
    sage: a + b
    4*5 + 3*5^2 + 3*5^3 + 4*5^4 + O(5^5)
    sage: a + b + c
    4*5 + 4*5^2 + 5^4 + O(5^5)

p-Adic fields should be created using the creation function Qp as
above.  This will ensure that there is only one instance of $\Q_p$ of
a given type, p, print mode and precision.  It also saves typing very
long class names.
    sage: Qp(17,10,'capped-rel')
    17-adic Field with capped relative precision 10
    sage: Qp(7, prec = 30, type = 'lazy', print_mode = 'val-unit')
    Lazy 7-adic Field
    sage: R = Qp(7, prec = 20, type = 'capped-rel', print_mode = 'val-unit'); S = Qp(7, prec = 20, type = 'capped-rel', print_mode = 'val-unit'); R is S
    True
    sage: Qp(2)
    2-adic Field with capped relative precision 20

Once one has a p-Adic field, one can cast elements into it in the
standard way.  Integers, ints, longs, Rationals, other p-Adic types,
pari p-adics and elements of $\Z / p^n \Z$ can all be cast into a
p-Adic field.
    sage: R = Qp(5, 5, 'capped-rel','series'); a = R(16); a
    1 + 3*5 + O(5^5)
    sage: b = R(23/15); b
    5^-1 + 3 + 3*5 + 5^2 + 3*5^3 + O(5^4)
    sage: S = Zp(5, 5, 'fixed-mod','val-unit'); c = S(Mod(75,125)); c
    5^2 * 3 + O(5^5)
    sage: R(c)
    3*5^2 + O(5^5)

In the previous example, since fixed-mod elements don't keep track of their precision, we assume that it has the full precision of the ring.  This is why you have to cast manually here.

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

#import sage.libs.all
#import sage.rings.integer_mod
import sage.rings.integer
import sage.rings.rational
#import sage.rings.finite_field
import sage.rings.padics.padic_capped_relative_element
import sage.rings.padics.padic_field_base_generic
import sage.rings.padics.padic_capped_relative_field_generic
import padic_base_generic_element
import padic_lazy_element
#import sage.rings.padics.padic_lazy_element
#import sage.rings.padics.padic_ring_capped_relative_element
#import sage.rings.padics.padic_ring_capped_absolute_element
#import sage.rings.padics.padic_ring_generic

#pari = sage.libs.pari.gen.pari
Integer = sage.rings.integer.Integer
Rational = sage.rings.rational.Rational
#Mod = sage.rings.integer_mod.Mod
#infinity = sage.rings.infinity.infinity
pAdicFieldBaseGeneric = sage.rings.padics.padic_field_base_generic.pAdicFieldBaseGeneric
pAdicCappedRelativeFieldGeneric = sage.rings.padics.padic_capped_relative_field_generic.pAdicCappedRelativeFieldGeneric
#pAdicLazyElement = sage.rings.padics.padic_lazy_element.pAdicLazyElement
pAdicBaseGenericElement = padic_base_generic_element.pAdicBaseGenericElement
pAdicLazyElement = padic_lazy_element.pAdicLazyElement
pAdicCappedRelativeElement = sage.rings.padics.padic_capped_relative_element.pAdicCappedRelativeElement
#pAdicRingCappedRelativeElement = sage.rings.padics.padic_ring_capped_relative_element.pAdicRingCappedRelativeElement
#pAdicRingCappedAbsoluteElement = sage.rings.padics.padic_ring_capped_absolute_element.pAdicRingCappedAbsoluteElement
#Zp = sage.rings.padics.zp.Zp

class pAdicFieldCappedRelative(pAdicFieldBaseGeneric, pAdicCappedRelativeFieldGeneric):
    def __init__(self, p, prec, print_mode, names):
        pAdicFieldBaseGeneric.__init__(self, p, prec, print_mode, names, pAdicCappedRelativeElement)

    r"""
    An implementation of p-adic fields with capped relative precision.
    """
    def __contains__(self, x):
        if isinstance(x, (int, long, Integer, Rational)):
            return True
        if not isinstance(x, (pAdicBaseGenericElement, pAdicLazyElement)):
            return False
        par = x.parent()
        if par.prime() != self.prime():
            return False
        if par.is_capped_relative() or par.is_capped_absolute() or par.is_lazy():
            if par.is_capped_relative() and par.is_field() and par.precision_cap() < self.precision_cap():
                return False
            return True
        return False

    def _repr_(self, do_latex=False):
        return "%s-adic Field with capped relative precision %s"%(self.prime(), self.precision_cap())

    def integer_ring(self):
        r"""
        Returns the integer ring of self, i.e. an appropriate implementation of $\Z_p$.

        INPUT:
            self -- a p-adic field
        OUTPUT:
            the p-adic ring that is the integer ring of this field
        """
        from sage.rings.padics.factory import Zp
        return Zp(self.prime(), self.precision_cap(), 'capped-rel', self.print_mode())
