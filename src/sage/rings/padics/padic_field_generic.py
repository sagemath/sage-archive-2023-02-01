"""
p-Adic Fields

AUTHOR:
    -- David Roe
    -- Genya Zaytman: documentation
    -- David Harvey: doctests

EXAMPLES:
p-Adic Fields are examples of inexact structures, as the reals are.  That means that elements cannot generally be stored exactly: to do so would take an infinite amount of storage.  Instead, we store an approximation to the elements with varying precision.

There are two types of precision for a p-adic element.  The first is relative precision, which gives the number of known p-adic digits.
    sage: R = Qp(5, 20, 'capped-rel', 'series'); a = R(675); a
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

There are two types of p-adic field: capped relative and lazy.
In the capped relative case, the relative precision of an element is restricted to be at most a certain value, specified at the creation of the field.  Individual elements also store their own precision, so the effect of various arithmetic operations on precision is tracked.  When you cast an exact element into a capped relative field, it truncates it to the precision cap of the field.
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

In the lazy case, a certain number of digits are computed and stored, and in addition a function is stored so that additional digits can be computed later.  In order to set the number of known digits you call set_precision_absolute().
    sage: R = Qp(5, 5, 'lazy', 'series'); a = R(4006); a
        1  + 5 + 2*5^3 + 5^4 + O(5^5)
    sage: b = R(50127); b
        2 + 5^3 + O(5^5)
    sage: c = a * b; c
        2 + 2*5 + 4*5^4 + O(5^5)
    sage: c.set_precision_absolute(15)
    sage: c
        2 + 2*5 + 4*5^4 + 3*5^5 + 5^6 + 4*5^8 + 2*5^9 + 4*5^11 + O(5^15)

There is some performance penalty for carrying the function around, but it is minimized if you determine the precision you will need going into a computation and set the precision appropriately at the outset.

p-Adic fields should be created using the creation function Qp as above.  This will ensure that there is only one instance of $\Q_p$ of a given type, p and precision.  It also saves typing very long class names.
    sage: Qp(17,10,'capped-rel')
        17-adic Field with capped relative precision 10
    sage: Qp(7, prec = 30, type = 'lazy', print_mode = 'val-unit')
        Lazy 7-adic Field
    sage: R = Qp(7, prec = 20, type = 'capped-rel', print_mode = 'val-unit'); S = Qp(7, prec = 20, type = 'capped-rel', print_mode = 'series'); R is S
        True
    sage: Qp(2)
        2-adic Field with capped relative precision 20

Once one has a p-Adic field, one can cast elements into it in the standard way.  Integers, ints, longs, Rationals, other p-Adic types, pari p-adics and elements of $\Z / p^n \Z$ can all be cast into a p-Adic field.
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


#import weakref
import sage.rings.padics.padic_generic
#import sage.rings.padics.padic_field_capped_relative
#import sage.rings.padics.padic_field_lazy

#pAdicFieldCappedRelative = sage.rings.padics.padic_field_capped_relative.pAdicFieldCappedRelative
#pAdicFieldLazy = sage.rings.padics.padic_field_lazy.pAdicFieldLazy

def is_pAdicField(R):
    return isinstance(R, pAdicFieldGeneric)

class pAdicFieldGeneric(sage.rings.padics.padic_generic.pAdicGeneric):

    def __init__(self, p, prec, print_mode):
        sage.rings.padics.padic_generic.pAdicGeneric.__init__(self, p, prec)
        self.set_print_mode(print_mode)

    def _repr_(self, do_latex = False):
        return "Generic %s-adic Field"%(self.prime())

    def get_print_mode(self):
        r"""
        Returns the current print mode as a string.

        INPUT:
            self -- a p-adic field

        OUTPUT:
            string -- self's print mode

        EXAMPLES:
            sage: R = Qp(7,5, 'capped-rel')
            sage: R.get_print_mode()
                'val-unit'
        """
        return self._print_mode

    def set_print_mode(self, print_mode):
        """
        Sets the print mode.

        INPUT:
            self -- a p-adic field
            print_mode -- string (see NOTES)

        EXAMPLES:
            sage: R = Qp(3,5,'capped-rel','val-unit')
            sage: a = R(117); a
                3^2 * 13 + O(3^7)
            sage: R.set_print_mode('series'); a
                3^2 + 3^3 + 3^4 + O(3^7)
            sage: R.set_print_mode('val-unit-p'); a
                p^2 * 13 + O(p^7)
            sage: R.set_print_mode('series-p'); a
                p^2 + p^3 + p^4 + O(p^7)

        NOTES:
            The options are:
            'val-unit' -- elements are displayed as p^k*u
            'integer' -- elements are displayed as an integer
            'series' -- elements are displayed as series in p
            'val-unit-p' -- same as val-unit, except that p is written as "p"
            'integer-p' -- same as integer, except that p is written as "p"
            'series-p' -- same as series, except that p is written as "p"
        """
        if (print_mode in ['val-unit', 'series', 'val-unit-p', 'series-p']):
            self._print_mode = print_mode
        else:
            raise ValueError, "print_mode must be either val-unit, integer, series, val-unit-p, integer-p, or series-p"

    def krull_dimension(self):
        """
        Returns the Krull dimension of $\Q_p$, i.e. 0.

        INPUT:
            self -- a p-adic field

        OUTPUT:
            integer -- self's krull dimension, i.e., 0

        EXAMPLES:
            sage: K = Qp(3, 10,'capped-rel'); K.krull_dimension()
                0
        """
        return 0

    def is_field(self):
        """
        Returns whether this p-adic field is a field, i.e. True.

        INPUT:
            self -- a p-adic field

        OUTPUT:
            boolean -- whether self is a field, i.e., True

        """
        return True

    def is_isomorphic(self, ring):
        r"""
        Returns whether self and ring are isomorphic, i.e. whether ring is an implementation of $\Q_p$ for the same prime as self.

        INPUT:
            self -- a p-adic field
            ring -- a ring

        OUTPUT:
            boolean -- whether ring is an implementation of $\Q_p$ for the same prime as self.
        """
        return is_instance(ring, pAdicFieldGeneric) and self.prime() == ring.prime()

    def class_field(self, group=None, map=None, generators=None):
        raise NotImplementedError

    def composite(self, subfield1, subfield2):
        r"""
        Returns the composite of two subfields of self, i.e., the largest subfield containing both

        INPUT:
            self -- a p-adic field
            subfield1 -- a subfield
            subfield2 -- a subfield

        OUTPUT:
            subfield -- the composite of subfield1 and subfield2
        """
        if (subfield1 is self) and (subfield2 is self):
            return self
        raise ValueError, "Arguments must be subfields of self."

    def norm_equation(self):
        raise NotImplementedError

    def norm_group(self):
        raise NotImplementedError

    def norm_group_discriminant(self, group=None, map=None, generators=None):
        raise NotImplementedError

    def number_of_extensions(self, degree, discriminant=None, e=None, f=None):
        raise NotImplementedError

    def list_of_extensions(self, degree, discriminant=None, e=None, f=None):
        raise NotImplementedError

    def subfield(self, list):
        r"""
        Returns the subfield generated by the elements in list

        INPUT:
            self -- a p-adic field
            list -- a list of elements of self

        OUTPUT:
            subfield -- the subfield of self generated by the elements of list
        """
        for x in list:
            if not self.__contains__(x):
                raise TypeError, "Members of the list of generators must be elements of self."
        return self

    def subfield_lattice(self):
        raise NotImplementedError

    def subfields_of_degree(self, n):
        r"""
        Returns the number of subfields of self of degree n

        INPUT:
            self -- a p-adic field
            n -- an integer

        OUTPUT:
            integer -- the number of subfields of degree n
        """
        if n == 1:
            return 1
        else:
            return 0

class pAdicFieldBaseGeneric(pAdicFieldGeneric):
    pass
