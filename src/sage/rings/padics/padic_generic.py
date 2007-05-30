"""
p-Adics

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

While you can cast explicitly, the chains of automatic coercion are more restricted.  As always in SAGE, the following arrows are transitive and the diagram is commutative.

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

import sage.rings.infinity
import sage.rings.finite_field
import sage.rings.padics.local_generic
import sage.rings.integer_mod
import sage.rings.ring

from sage.rings.padics.pow_computer import PowComputer
from sage.rings.integer import Integer

infinity = sage.rings.infinity.infinity
Mod = sage.rings.integer_mod.Mod

class pAdicGeneric(sage.rings.ring.PrincipalIdealDomain,
                   sage.rings.padics.local_generic.LocalGeneric):
    def __init__(self, p, prec, print_mode, names, element_class):
        sage.rings.padics.local_generic.LocalGeneric.__init__(self, prec, names)
        #if prec <= 100:
        #    self.prime_pow = PowComputer(p, prec, prec, self.is_field())
        #else:
        #    self.prime_pow = PowComputer(p, 3, prec, self.is_field())
        self.prime_pow = PowComputer(p, prec, self.is_field())
        self.__set_print_mode(print_mode)
        self._element_class = element_class

    def _repr_(self, do_latex = False):
        return "Generic %s-adic Parent."%(self.prime())

    def __call__(self, x, absprec = infinity, relprec = infinity):
        r"""
            Casts x into self.  Uses the constructor of self._element_class.
        """
        return self._element_class(self, x, absprec, relprec)

    def _coerce_impl(self, x):
        if self.__contains__(x):
            return self.__call__(x)
        else:
            raise TypeError, "no canonical coercion of %s of type %s into %s"%(x, type(x), self)

    def __cmp__(self, other):
        if isinstance(other, type(self)):
            if self.prime() < other.prime():
                return -1
            elif self.prime() > other.prime():
                return 1
            try:
                if self.halting_parameter() < other.halting_parameter():
                    return -1
                elif self.halting_parameter() > other.halting_parameter():
                    return 1
            except AttributeError:
                pass
            if self.precision_cap() < other.precision_cap():
                return -1
            elif self.precision_cap() > other.precision_cap():
                return 1
            else:
                return 0
        else:
            return cmp(type(self), type(other))

    def ngens(self):
        return 1

    def gen(self, n = 0):
        if n != 0:
            raise IndexError, "only one generator"
        return self(self.prime())

    def print_mode(self):
        r"""
        Returns the current print mode as a string.

        INPUT:
            self -- a p-adic field

        OUTPUT:
            string -- self's print mode

        EXAMPLES:
            sage: R = Qp(7,5, 'capped-rel')
            sage: R.print_mode()
            'series'
        """
        return self._print_mode

    def __set_print_mode(self, print_mode):
        """
        Sets the print mode.

        WARNING: You should not use this function.

        INPUT:
            self -- a p-adic ring
            print_mode -- string (see NOTES)

        NOTES:
            The options for print_mode are:
            'val-unit' -- elements are displayed as p^k*u
            'terse' -- elements are displayed as an integer if positive valuation, as u/ppow or u/p^k if negative valuation
            'series' -- elements are displayed as series in p, where p is self.variable_name() (default, e.g., "5")
        """
        if (print_mode in ['val-unit', 'terse', 'series']):
            try:
                old = self._print_mode
                self._print_mode = print_mode
                return old
            except AttributeError:
                self._print_mode = print_mode
        else:
            raise ValueError, "print_mode=%s must be either val-unit, terse, series"%print_mode

    def _element_class(self):
        return self._element_class

    def characteristic(self):
        r"""
        Returns the characteristic of self, which is always 0.

        INPUT:
            self -- a p-adic parent

        OUTPUT:
            integer -- self's characteristic, i.e., 0

        EXAMPLES:
            sage: R = Zp(3, 10,'fixed-mod'); R.characteristic()
            0
        """
        return 0

    def prime(self):
        """
        Returns the prime, ie the characteristic of the residue field.

        INPUT:
            self -- a p-adic parent

        OUTPUT:
            integer -- the characteristic of the residue field

        EXAMPLES:
            sage: R = Zp(3,5,'fixed-mod')
            sage: R.prime()
                3
        """
        return self.prime_pow._prime()

    def uniformizer_pow(self, n):
        if n is infinity:
            return self(0)
        return self.uniformizer()**n

    def residue_characteristic(self):
        """
        Returns the prime, i.e., the characteristic of the residue field.

        INPUT:
            self -- a p-adic ring

        OUTPUT:
            integer -- the characteristic of the residue field

        EXAMPLES:
            sage: R = Zp(3,5,'fixed-mod')
            sage: R.residue_characteristic()
                3
        """
        return self.prime()

    def residue_class_field(self):
        """
        Returns the residue class field.

        INPUT:
            self -- a p-adic ring

        OUTPUT:
            the residue field

        EXAMPLES:
            sage: R = Zp(3,5,'fixed-mod')
            sage: k = R.residue_class_field()
            sage: k
                Finite Field of size 3
        """
        return sage.rings.finite_field.GF(self.prime())

    def residue_system(self):
        """
        Returns a list of elements representing all the residue classes.

        INPUT:
            self -- a p-adic ring

        OUTPUT:
            list of elementss -- a list of elements representing all the residue classes

        EXAMPLES:
            sage: R = Zp(3, 5,'fixed-mod')
            sage: R.residue_system()
                [O(3^5), 1 + O(3^5), 2 + O(3^5)]
        """
        return [self(i) for i in self.residue_class_field()]

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
            sage: R = Qp(5, 10,'capped-rel','series')
            sage: R.teichmuller(2)
            2 + 5 + 2*5^2 + 5^3 + 3*5^4 + 4*5^5 + 2*5^6 + 3*5^7 + 3*5^9 + O(5^10)
            sage: R = Zp(5, 10, 'capped-abs', 'series')
            sage: R.teichmuller(2)
            2 + 5 + 2*5^2 + 5^3 + 3*5^4 + 4*5^5 + 2*5^6 + 3*5^7 + 3*5^9 + O(5^10)
            sage: R = Zp(5, 10, 'fixed-mod', 'series')
            sage: R.teichmuller(2)
            2 + 5 + 2*5^2 + 5^3 + 3*5^4 + 4*5^5 + 2*5^6 + 3*5^7 + 3*5^9 + O(5^10)

        AUTHORS:
        Initial version: David Roe
        Quadratic time version: Kiran Kedlaya <kedlaya@math.mit.edu> (3/27/07)
        """
        if prec is None:
            prec = self.precision_cap()
        else:
            prec = min(Integer(prec), self.precision_cap())
        x = Integer(x)
        ans = self._element_class(self, None, empty = True)
        ans._teichmuller_set(x, prec)
        return ans

    def teichmuller_system(self):
        r"""
        Returns a set of teichmuller representatives for the invertible elements of $\Z / p\Z$.

        INPUT:
            self -- a p-adic ring

        OUTPUT:
            list of elements -- a list of teichmuller representatives for the invertible elements of $\Z / p\Z$

        EXAMPLES:
            sage: R = Zp(3, 5,'fixed-mod', 'terse')
            sage: R.teichmuller_system()
            [1 + O(3^5), 242 + O(3^5)]

        NOTES:
            Should this return 0 as well?
        """
        return [self.teichmuller(i) for i in self.residue_class_field() if i != 0]

    def absolute_discriminant(self):
        """
        Returns the absolute discriminant of this p-adic ring

        INPUT:
            self -- a p-adic ring

        OUTPUT:
            integer -- the absolute discriminant of this ring, i.e., 1
        """
        return 1

    def discriminant(self, K=None):
        """
        Returns the discriminant of this p-adic ring over K

        INPUT:
            self -- a p-adic ring
            K -- a sub-ring of self or None (default: None)

        OUTPUT:
            integer -- the discriminant of this ring over K (or the absolute discriminant if K is None)
        """
        if (K is None or K is self):
            return 1
        else:
            raise ValueError, "Ground Ring must be a subring of self"

    def different(self):
        raise NotImplementedError

    def automorphisms(self):
        r"""
        Returns the group of automorphisms of $\Z_p$, i.e. the trivial group.
        """
        raise NotImplementedError

    def galois_group(self):
        r"""
        Returns the Galois group of $\Z_p$, i.e. the trivial group.
        """
        raise NotImplementedError

    def is_abelian(self):
        """
        Returns whether the Galois group is abelian, i.e. True.  #should this be automorphism group?

        INPUT:
            self -- a p-adic ring

        OUTPUT:
            boolean -- whether self is abelian, i.e., True

        EXAMPLES:
            sage: R = Zp(3, 10,'fixed-mod'); R.is_abelian()
                True
        """
        return True

    def is_normal(self):
        """
        Returns whether or not this is a normal extension, i.e. True.

        INPUT:
            self -- a p-adic ring

        OUTPUT:
            boolean -- whether self is normal, i.e., True

        EXAMPLES:
            sage: R = Zp(3, 10,'fixed-mod'); R.is_normal()
                True
        """
        return True

    def uniformizer(self):
        """
        Returns a uniformizer for this ring.

        INPUT:
            self -- a p-adic ring

        OUTPUT:
            element -- self's uniformizer

        EXAMPLES:
            sage: R = Zp(3,5,'fixed-mod', 'series')
            sage: R.uniformizer()
                3 + O(3^5)
        """
        return self(self.prime_pow._prime())

    def has_pth_root(self):
        r"""
        Returns whether or not $\Z_p$ has a $p^{\mbox{th}}$ root of unity.

        INPUT:
            self -- a p-adic ring

        OUTPUT:
            boolean -- whether self has $p^{\mbox{th}}$ root of unity
        """
        return (self.prime() == 2)

    def has_root_of_unity(self, n):
        r"""
        Returns whether or not $\Z_p$ has an $n^{\mbox{th}}$ root of unity.

        INPUT:
            self -- a p-adic ring
            n -- an integer

        OUTPUT:
            boolean -- whether self has $n^{\mbox{th}}$ root of unity
        """
##        p = self.prime()
##        if (p == 2) and (n % 2 == 0):
##            return True
##        if p.divides(n):
##            return False
##        if n == 1:
##            return True
##        if gcd(n, p-1) > 1:
##            return True
##        return False
	##
	## I'm not sure why the above definition existed. I'm keeping
	## it for now until at least one other person looks and says
	## I'm not missing something.
	##
	if (self.prime() == 2):
	    return n.divides(2)
        else:
	    return n.divides(self.prime() - 1)

    def hasGNB(self):
        r"""
        Returns whether or not $\Z_p$ has a Gauss Normal Basis.
        """
        raise NotImplementedError

    def zeta(self, n=None):
        r"""
        Returns a generator of the group of roots of unity.

        INPUT:
            self -- a p-adic ring
            n -- an integer ot None (default: None)

        OUTPUT:
            element -- a generator of the $n^\mbox{th}$ roots of unity, or a generator of the full group of roots of unity if n is None
        """
        if (self.prime() == 2):
            if (n == None) or (n == 2):
                return self(-1)
            if n == 1:
                return self(1)
            else:
                raise ValueError, "No, %sth root of unity in self"%n
        else:
            return self.teichmuller(sage.rings.finite_field.GF(self.prime()).zeta(n).lift())

    def zeta_order(self):
        """
        Returns the order of the group of roots of unity.

        INPUT:
            self -- a p-adic ring

        OUTPUT:
            integer -- the order of the group of roots of unity
        """
        if (self.prime() == 2):
            return 2
        else:
            return self.prime() - 1

    def extension(self, modulus, prec = None, names = None, print_mode = None, halt = None):
        """
        Create an extension of this p-adic ring.

        WARNING -- this isn't ready for general use yet.

        EXAMPLES:
            sage: k = Qp(5)
            sage: R.<x> = k[]
            sage: l.<a> = k.extension(x^2 + 5)
            Traceback (most recent call last):
            ...
            NotImplementedError: This is not yet ready for general use.
        """
        raise NotImplementedError, "This is not yet ready for general use."
        if not self is modulus.base_ring():
            modulus = modulus.parent().change_ring(self)(modulus)
        from sage.rings.padics.factory import ExtensionFactory
        return ExtensionFactory(modulus, prec, print_mode, halt, names, check = True)

    ext = extension

class local_print_mode:
    r"""
    Context manager for safely temporarily changing the print_mode
    of a p-adic ring/field.

    EXAMPLES:
        sage: R = Zp(5)
        sage: R(45)
        4*5 + 5^2 + O(5^21)
        sage: with local_print_mode(R, 'val-unit'):
        ...       print R(45)
        ...
        5 * 9 + O(5^21)

    NOTES:  For more documentation see localvars in parent_gens.pyx
    """
    def __init__(self, obj, print_mode):
        self._obj = obj
        self._print_mode = print_mode

    def __enter__(self):
        self._orig = self._obj._pAdicGeneric__set_print_mode(self._print_mode)

    def __exit__(self, type, value, traceback):
        self._obj._pAdicGeneric__set_print_mode(self._orig)
