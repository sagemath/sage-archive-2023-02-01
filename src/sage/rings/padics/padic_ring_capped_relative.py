r"""
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

The second type of precision is absolute precision, which gives the
power of p that this element is stored modulo.
    sage: a.precision_absolute()
    22

The number of times that p divides the element is called the valuation, and can be accessed with the functions valuation() and ordp()
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

    #sage: Zp(7, prec = 30, type = 'lazy', print_mode = 'val-unit')
    #Lazy 7-adic Ring
    sage: R = Zp(7, prec = 20, type = 'capped-rel', print_mode = 'val-unit'); S = Zp(7, prec = 20, type = 'capped-rel', print_mode = 'val-unit'); R is S
    True
    sage: Zp(2)
    2-adic Ring with capped relative precision 20

Once one has a p-Adic ring, one can cast elements into it in the
standard way.  Integers, ints, longs, Rationals, other p-Adic types,
pari p-adics and elements of $\Z / p^n \Z$ can all be cast into a
p-Adic ring.

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

(WARNING: This is wrong.)
\begin{verbatim}
int -> long -> Integer -> Rational -> Zp lazy -> pari p-adic -> Zp capped-rel -> Zp capped-abs
                                       /    \                          \          /
                                      /      \                          \        /
                                     V        V                          V      V
                           Zp fixed-mod     Qp lazy ----------------> Qp capped-rel
                                     \_______________________          /
                                                             \        /
                                                              V      V
                                                             IntegerMod
\end{verbatim}

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



import sage.rings.padics.padic_ring_base_generic
#import sage.rings.padics.padic_generic
import sage.rings.padics.padic_capped_relative_element
#import sage.rings.padics.padic_ring_fixed_mod
#import sage.rings.padics.padic_lazy_element
#import sage.rings.infinity
import padic_base_generic_element
import padic_capped_relative_ring_generic
import sage.rings.integer
import sage.rings.integer_ring

#from sage.rings.integer_ring import ZZ

ZZ = sage.rings.integer_ring.ZZ
Integer = sage.rings.integer.Integer
#Mod = sage.rings.integer_mod.Mod
#infinity = sage.rings.infinity.infinity
#Element = sage.structure.element.Element
pAdicRingBaseGeneric = sage.rings.padics.padic_ring_base_generic.pAdicRingBaseGeneric
pAdicBaseGenericElement = padic_base_generic_element.pAdicBaseGenericElement
pAdicCappedRelativeRingGeneric = sage.rings.padics.padic_capped_relative_ring_generic.pAdicCappedRelativeRingGeneric
#pAdicGeneric = sage.rings.padics.padic_generic.pAdicGeneric
#pAdicRingFixedMod = sage.rings.padics.padic_ring_fixed_mod.pAdicRingFixedMod
pAdicCappedRelativeElement = sage.rings.padics.padic_capped_relative_element.pAdicCappedRelativeElement
#pAdicLazyElement = sage.rings.padics.padic_lazy_element.pAdicLazyElement

class pAdicRingCappedRelative(pAdicRingBaseGeneric, pAdicCappedRelativeRingGeneric):
    def __init__(self, p, prec, print_mode, names):
        pAdicRingBaseGeneric.__init__(self, p, prec, print_mode, names, pAdicCappedRelativeElement)

    r"""
    An implementation of the p-adic integers with capped relative precision.
    """
    def __contains__(self, x):
        if isinstance(x, (int, long, Integer)):
            return True
        if not isinstance(x, pAdicBaseGenericElement):
            return False
        if x.parent().prime() != self.prime():
            return False
        if x.parent().is_field():
            return False
        if x.parent().is_capped_relative() and x.parent().precision_cap() >= self.precision_cap():
            return True
        if x.parent().is_lazy():
            return True
        return False

    def _repr_(self, do_latex=False):
        return "%s-adic Ring with capped relative precision %s"%(self.prime(), self.precision_cap())

    def fraction_field(self, print_mode=None):
        r"""
        Returns the fraction field of self.

        INPUT:
        print_mode - a dictionary containing print options.  Defaults to the same options as this ring.
        OUTPUT:
        the fraction field of self.

        EXAMPLES:
            sage: R = Zp(5, print_mode='digits')
            sage: K = R.fraction_field(); repr(K(1/3))
            '...31313131313131313132'
            sage: L = R.fraction_field({'max_ram_terms':4}); repr(L(1/3))
            '...3132'
        """
        from sage.rings.padics.factory import Qp
        if print_mode is None:
            print_mode = {}
        elif isinstance(print_mode, str):
            print_mode = {'mode': print_mode}
        for option in ['mode', 'pos', 'ram_name', 'unram_name', 'var_name', 'max_ram_terms', 'max_unram_terms', 'max_terse_terms', 'sep', 'alphabet']:
            if not print_mode.has_key(option):
                print_mode[option] = self._printer.dict()[option]
        return Qp(self.prime(), self.precision_cap(), print_mode=print_mode)

    def integer_ring(self, print_mode=None):
        r"""
        Returns the integer ring of self, possibly with print_mode changed.
        """
        if print_mode is None:
            return self
        from sage.rings.padics.factory import ZpCR
        if print_mode is None:
            print_mode = {}
        elif isinstance(print_mode, str):
            print_mode = {'mode': print_mode}
        for option in ['mode', 'pos', 'ram_name', 'unram_name', 'var_name', 'max_ram_terms', 'max_unram_terms', 'max_terse_terms', 'sep', 'alphabet']:
            if not print_mode.has_key(option):
                print_mode[option] = self._printer.dict()[option]
        return ZpCR(self.prime(), self.precision_cap(), print_mode=print_mode)


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

