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

The fixed modulus type is the leanest of the p-adic rings: it is basically just a wrapper around $\Z / p^n \Z$ providing a unified interface with the rest of the p-adics.  This is the type you should use if your primary interest is in speed.  It does not track precision of elements.
    sage: R = Zp(5,5,'fixed-mod','series'); a = R(4005); a
        5 + 2*5^3 + 5^4 + O(5^5)
    sage: a // 5
        1 + 2*5^2 + 5^3 + O(5^5)

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

from sage.rings.integer import Integer
from sage.rings.integer_mod_ring import Integers
from sage.rings.padics.padic_ring_base_generic import pAdicRingBaseGeneric
from sage.rings.padics.padic_fixed_mod_element import pAdicFixedModElement
from sage.rings.padics.padic_fixed_mod_ring_generic import pAdicFixedModRingGeneric

class pAdicRingFixedMod(pAdicRingBaseGeneric, pAdicFixedModRingGeneric):
    r"""
    An implementation of the p-adic integers using fixed modulus.
    """
    def __init__(self, p, prec, print_mode, names):
        pAdicRingBaseGeneric.__init__(self, p, prec, print_mode, names, pAdicFixedModElement)

    def __contains__(self, x):
        if isinstance(x, (int, long, Integer)):
            return True
        if isinstance(x.parent(), pAdicRingFixedMod) and x.parent().prime() != self.prime() and x.parent().precision_cap() >= self.precision_cap():
            return True
        return False

    def _repr_(self, do_latex=False):
        if do_latex:
            return "%s-adic Ring of fixed modulus %s^{%s}"%(self.prime(), self.prime(), self.precision_cap())
        else:
            return "%s-adic Ring of fixed modulus %s^%s"%(self.prime(), self.prime(), self.precision_cap())

    def fraction_field(self, print_mode = None):
        r"""
        Would normally return $\Q_p$, but there is no implementation of $Q_p$ matching this ring so this raises an error

        If you want to be able to divide with elements of a fixed modulus p-adic ring, you must cast explicitly.
        """
        raise TypeError, "This implementation of the p-adic ring does not support fields of fractions."

    def integer_ring(self, print_mode=None):
        r"""
        Returns the integer ring of self, possibly with print_mode changed.
        """
        if print_mode is None:
            return self
        from sage.rings.padics.factory import ZpFM
        if print_mode is None:
            print_mode = {}
        elif isinstance(print_mode, str):
            print_mode = {'mode': print_mode}
        for option in ['mode', 'pos', 'ram_name', 'unram_name', 'var_name', 'max_ram_terms', 'max_unram_terms', 'max_terse_terms', 'sep', 'alphabet']:
            if not print_mode.has_key(option):
                print_mode[option] = self._printer.dict()[option]
        return ZpFM(self.prime(), self.precision_cap(), print_mode=print_mode)

    def random_element(self):
        """
        Returns a random element of self, selected from a uniform distribution.
        INPUT:
            self -- a p-adic ring
        OUTPUT:
            element -- a random element of self
        """
        return self(Integers(self.prime()**self.precision_cap()).random_element())

