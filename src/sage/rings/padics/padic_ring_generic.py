"""
p-Adic Rings

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

There are four types of p-adic ring: capped relative, capped absolute, fixed modulus and lazy.
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
    <class 'sage.rings.padics.padic_ring_capped_absolute_element.pAdicRingCappedAbsoluteElement'>
    sage: (a * b) / 5^3
    1 + 2*5 + O(5^2)
    sage: type((a * b) / 5^3)
    <class 'sage.rings.padics.padic_capped_relative_element.pAdicCappedRelativeElement'>

The fixed modulus type is the leanest of the p-adic rings: it is basically just a wrapper around $\Z / p^n \Z$ providing a unified interface with the rest of the p-adics.  This is the type you should use if your primary interest is in speed.  It does not track precision of elements.
    sage: R = Zp(5,5,'fixed-mod','series'); a = R(4005); a
    5 + 2*5^3 + 5^4 + O(5^5)
    sage: a // 5
    1 + 2*5^2 + 5^3 + O(5^5)

In the lazy case, a certain number of digits are computed and stored, and in addition a function is stored so that additional digits can be computed later.  In order to set the number of known digits you call_set_precision_absolute().
    sage: R = Zp(5, 5, 'lazy', 'series'); a = R(4006); a
    1  + 5 + 2*5^3 + 5^4 + O(5^5)
    sage: b = R(50127); b
    2 + 5^3 + O(5^5)
    sage: c = a * b; c
    2 + 2*5 + 4*5^4 + O(5^5)
    sage: c.set_precision_absolute(15)
    sage: c
    2 + 2*5 + 4*5^4 + 3*5^5 + 5^6 + 4*5^8 + 2*5^9 + 4*5^11 + O(5^15)

There is some performance penalty for carrying the function around, but it is minimized if you determine the precision you will need going into a computation and set the precision appropriately at the outset.

p-Adic rings should be created using the creation function Zp as
above.  This will ensure that there is only one instance of $\Z_p$ of
a given type, p, print mode and precision.  It also saves typing very long class
names.
    sage: Zp(17,10,'capped-rel')
    17-adic Ring with capped relative precision 10
    sage: Zp(7, prec = 30, type = 'lazy', print_mode = 'val-unit')
    Lazy 7-adic Ring
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

Note that in the last example, since fixed-mod elements don't keep
track of their precision, we assume that it has the full precision of
the ring.  This is why you have to cast manually there.

While you can cast explicitly as above, the chains of automatic
coercion are more restricted.  As always in SAGE, the following arrows
are transitive and the diagram is commutative.

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


#import weakref
import sage.rings.padics.padic_generic
import sage.rings.ring
#import sage.rings.padics.padic_ring_capped_relative
#import sage.rings.padics.padic_ring_capped_absolute
#import sage.rings.padics.padic_ring_fixed_mod
#import sage.rings.padics.padic_ring_lazy
from sage.rings.integer import Integer


def is_pAdicRing(R):
    return isinstance(R, pAdicRingGeneric)

class pAdicRingGeneric(sage.rings.padics.padic_generic.pAdicGeneric, sage.rings.ring.EuclideanDomain):
    def is_field(self):
        return False

    def _repr_(self, do_latex = False):
        return "Generic %s-adic Ring"%(self.prime())

    def integer_ring(self):
        r"""
        Returns the integer ring of self, i.e. self.
        """
        return self

    def krull_dimension(self):
        r"""
        Returns the Krull dimension of self, i.e. 1

        INPUT:
            self -- a p-adic ring
        OUTPUT:
            the Krull dimension of self.  Since self is a p-adic ring, this is 1.
        """
        return 1
