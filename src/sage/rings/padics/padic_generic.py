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
#from sage.rings.padics.hack_infinity import plusinfinity

infinity = sage.rings.infinity.infinity

class pAdicGeneric(sage.rings.padics.local_generic.LocalGeneric):

    def __init__(self, p, prec, print_mode, names):
        sage.rings.padics.local_generic.LocalGeneric.__init__(self, prec, names)
        self._p = p
        if (print_mode in ['val-unit', 'terse', 'series']):
            self._print_mode = print_mode
        else:
            raise ValueError, "print_mode must be either val-unit, terse, series"

    def _repr_(self, do_latex = False):
        return "Generic %s-adic Parent."%(self.prime())

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

    #def set_print_mode(self, print_mode):
    #    """
    #    Sets the print mode.
    #
    #    INPUT:
    #        self -- a p-adic ring
    #        print_mode -- string (see NOTES)
    #
    #    EXAMPLES:
    #        sage: R = Zp(3,5,'fixed-mod'); R.set_print_mode('val-unit')
    #        sage: a = R(117); a
    #        3^2 * 13 + O(3^5)
    #        sage: R.set_print_mode('terse'); a
    #        117 + O(3^5)
    #        sage: R.set_print_mode('series'); a
    #        3^2 + 3^3 + 3^4 + O(3^5)
    #
    #    NOTES:
    #        The options for print_mode are:
    #        'val-unit' -- elements are displayed as p^k*u
    #        'terse' -- elements are displayed as an integer if positive valuation, as u/ppow or u/p^k if negative valuation
    #        'series' -- elements are displayed as series in p, where p is self.variable_name() (default, e.g., "5")
    #    """


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
        return self._p

    def prime_pow(self, n):
        """
        Returns the prime raised to the nth power.

        INPUT:
            self -- a p-adic field

        OUTPUT:
            integer -- p^n

        EXAMPLES:
            sage: R = Zp(3)
            sage: R.prime_pow(5)
                243
        """
        if n is infinity:
            return 0
        else:
            return self._p ** n

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
        return [self(i) for i in range(self.prime())]

    def teichmuller(self, x):
        raise NotImplementedError

    def teichmuller_system(self):
        r"""
        Returns a set of teichmuller representatives for the invertible elements of $\Z / p\Z$.

        INPUT:
            self -- a p-adic ring

        OUTPUT:
            list of elements -- a list of teichmuller representatives for the invertible elements of $\Z / p\Z$

        EXAMPLES:
            sage: R = Zp(3, 5,'fixed-mod', 'integer')
            sage: R.teichmuller_system()
                [1 + O(3^5), 242 + O(3^5)]
        """
        return [self.teichmuller(i) for i in range(1, self.prime())]

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
        return self(self._p)

    def _uniformizer_sym(self, do_latex = False):
        return "%s"%(self._p)

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
        from sage.rings.padics.extension_factory import ExtensionFactory
        return ExtensionFactory(self, modulus, prec, names, print_mode, halt, check)

    ext = extension
