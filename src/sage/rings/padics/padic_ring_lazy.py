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

In the lazy case, a certain number of digits are computed and stored, and in addition a function is stored so that additional digits can be computed later.  In order to set the number of known digits you call set_precision_absolute(prec) or set_precision_relative(prec).
    #sage: R = Zp(5, 5, 'lazy', 'series'); a = R(4006); a
    #    1  + 5 + 2*5^3 + 5^4 + O(5^5)
    #sage: b = R(50127); b
    #    2 + 5^3 + O(5^5)
    #sage: c = a * b; c
    #    2 + 2*5 + 4*5^4 + O(5^5)
    #sage: c.set_precision_absolute(15)
    #sage: c
    #    2 + 2*5 + 4*5^4 + 3*5^5 + 5^6 + 4*5^8 + 2*5^9 + 4*5^11 + O(5^15)

There is some performance penalty for carrying the function around, but it is minimized if you determine the precision you will need going into a computation and set the cache precision appropriately at the outset.

p-Adic rings should be created using the creation function Zp as above.  This will ensure that there is only one instance of $\Z_p$ of a given type, p and precision.  It also saves typing very long class names.
    #sage: Zp(7, prec = 30, type = 'lazy', print_mode = 'val-unit')
    #    Lazy 7-adic Ring
    #sage: R = Zp(7, prec = 20, type = 'lazy', print_mode = 'val-unit'); S = Zp(7, prec = 20, type = 'lazy', print_mode = 'val-unit'); R is S
    #    True
    sage: Zp(2)
        2-adic Ring with capped relative precision 20

Once one has a p-Adic ring, one can cast elements into it in the standard way.  Integers, ints, longs, Rationals, other p-Adic types, pari p-adics and elements of $\Z / p^n \Z$ can all be cast into a p-Adic ring.
    #sage: R = Zp(5, 5, 'lazy','series'); a = R(16); a
    #    1 + 3*5 + O(5^5)
    #sage: b = R(25/3); b
    #    2*5^2 + 3*5^3 + 5^4 + 3*5^5 + 5^6 + O(5^7)
    #sage: S = Zp(5, 5, 'fixed-mod','val-unit'); c = S(Mod(75,125)); c
    #    5^2 * 3 + O(5^5)
    #sage: R(c)
    #    3*5^2 + O(5^5)

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

import sage.rings.padics.padic_lazy_element as lazy
import padic_lazy_ring_generic
import padic_ring_base_generic
import sage.rings.infinity
import copy

pari = sage.libs.pari.gen.pari
pari_gen = sage.libs.pari.gen.gen
PariError = sage.libs.pari.gen.PariError
infinity = sage.rings.infinity.infinity
Mod = sage.rings.integer_mod.Mod
Integer = sage.rings.integer.Integer
Rational = sage.rings.rational.Rational
pAdicLazyRingGeneric = padic_lazy_ring_generic.pAdicLazyRingGeneric
pAdicGenericElement = sage.rings.padics.padic_generic_element.pAdicGenericElement
pAdicRingBaseGeneric = padic_ring_base_generic.pAdicRingBaseGeneric

class ErrorReporter:
    def __init__(self, par):
        self._obj = par

    def __call__(self, x, absprec = infinity, relprec = infinity):
        raise Exception, "There is a bug in the lazy code\n %s's generic methods are trying to create elements directly."

class pAdicRingLazy(pAdicRingBaseGeneric, pAdicLazyRingGeneric):
    r"""
    An implementation of the p-adic integers with lazily evaluated elements.
    """
    def __init__(self, p, prec, print_mode, halt, names):
        pAdicRingBaseGeneric.__init__(self, p, prec, print_mode, names, ErrorReporter(self))
        self._halt = halt

    def __call__(self, x, absprec = infinity, relprec = infinity):
        r"""
            Casts x into self.
        """
        if x == 0:
            return lazy.pAdicLazy_zero(self)
        if isinstance(x, (int, long)):
            x = Integer(x)

        if absprec < infinity and relprec < infinity:
            raise ValueError, "only one of absprec and relprec may be specified"
        if isinstance(x, Integer):
            return lazy.pAdicLazy_integer(self, x, absprec, relprec)
        if isinstance(x, Rational):
            if x.valuation(self.prime()) < 0:
                raise ValueError, "element not a p-adic integer"
            return lazy.pAdicLazy_rational(self, x, absprec, relprec)
        if isinstance(x, lazy.pAdicLazyElement):
            if x.parent().prime() != self.prime():
                raise TypeError, "cannot change primes in creating p-adic elements"
            if x.valuation() < 0:
                raise ValueError, "element not a p-adic integer"
            # In order to save code complexity, we copy x.  Instead, we could have a big if statment
            # and have pAdicLazy_**** accept a pAdicLazy_**** as an argument in addition to ****
            y = copy.copy(x)
            y._set_parent(self)
            if absprec < infinity:
                y.set_precision_absolute(absprec)
            elif relprec < infinity:
                y.set_precision_relative(relprec)
            return y
        if isinstance(x, pAdicGenericElement):
            if x.valuation() < 0:
                raise ValueError, "element not a p-adic integer"
            if x.parent().prime() == self.prime():
                return lazy.pAdicLazy_otherpadic(self, x, absprec, relprec)
            raise TypeError, "cannot change primes in creating p-adic elements"
        #if sage.rings.finite_field_element.is_FiniteFieldElement(x):
        #    if x.parent().order() == self.prime():
        #        return lazy.pAdicLazy_integer(self, x.lift(), absprec, relprec) #Note that this is cast in with full precision.  This may not be the desired behavior.
        #    raise TypeError, "cannot change primes in creating p-adic elements"
        if sage.rings.integer_mod.is_IntegerMod(x):
            k, p = pari(x.modulus()).ispower()
            if not k or p != parent.prime():
                raise TypeError, "cannot change primes in creating p-adic elements"
            return lazy.pAdicLazy_mod(self, x.lift(), min(absprec, k), relprec)
        if isinstance(x, pari_gen):
            if x.type() == "t_PADIC":
                try:
                    from sage.rings.padics.factory import Qp
                    return lazy.pAdicLazy_otherpadic(self, Qp(parent.prime(), x.padicprec(parent.prime()) - x.valuation(parent.prime()), 'capped-rel')(x), absprec, relprec)
                except PariError:
                    raise TypeError, "cannot change primes in creating p-adic elements"
            elif x.type() == "t_INT":
                return lazy.pAdicLazy_integer(self,Integer(x), absprec, relprec)
            elif x.type() == "t_FRAC":
                return lazy.pAdicLazy_rational(self, Rational(x), absprec, relprec)
            else:
                raise TypeError, "unsupported coercion from pari: only p-adics, integers and rationals allowed"
        raise TypeError, "Cannot create a p-adic out of %s"%(type(x))

    def __contains__(self, x):
        if isinstance(x, (int, long, Integer)):
            return True
        if isinstance(x, lazy.pAdicLazyElement) and x.parent().prime() == self.prime() and not x.parent().is_field():
            if x.parent().halting_parameter() > self.halting_parameter():
                return True
            if x.parent().halting_parameter() == self.halting_parameter() and x.parent().precision_cap() >= self.precision_cap():
                return True
        return False

    def _repr_(self, do_latex=False):
        return "Lazy %s-adic Ring"%(self.prime())

    def teichmuller(self, x, prec = None):
        r"""
        Returns the teichmuller representative of x.

        INPUT:
            self -- a p-adic ring
            x -- an integer or element of $\Z / p\Z$ that is not divisible by $p$
        OUTPUT:
            element -- the teichmuller lift of x
        EXAMPLES:
            #sage: R = Zp(5, 10, 'lazy', 'series')
            #sage: R.teichmuller(2)
            #    2 + 5 + 2*5^2 + 5^3 + 3*5^4 + 4*5^5 + 2*5^6 + 3*5^7 + 3*5^9 + O(5^10)
        """
        x = self.residue_class_field()(x)
        if x == 0:
            return lazy.pAdicLazy_zero(self)
        if prec is None:
            prec = self.precision_cap()
        return lazy.pAdicLazy_teichmuller(self, x, prec)

    def halting_parameter(self):
        return self._halt

    def fraction_field(self):
        r"""
        Returns the fraction field of self.
        """
        from sage.rings.padics.factory import Qp
        return Qp(self.prime(), self.precision_cap(), 'lazy', self.print_mode(), self.halting_parameter())

    def random_element(self):
        """
        Returns a random element of self.
        """
        raise NotImplementedError

    def unit_group(self):
        raise NotImplementedError

    def unit_group_gens(self):
        raise NotImplementedError

    def principal_unit_group(self):
        raise NotImplementedError




