"""
p-Adic Fields

AUTHOR:
    -- David Roe
    -- Genya Zaytman: documentation
    -- David Harvey: doctests

EXAMPLES:

p-Adic Fields are examples of inexact structures, as the reals are.
That means that elements cannot generally be stored exactly: to do so
would take an infinite amount of storage.  Instead, we store an
approximation to the elements with varying precision.

There are two types of precision for a p-adic element.  The first is
relative precision, which gives the number of known p-adic digits.

    sage: R = Qp(5, 20, 'capped-rel', 'series'); a = R(675); a
        2*5^2 + 5^4 + O(5^22)
    sage: a.precision_relative()
        20

The second type of precision is absolute precision, which gives the
power of p that this element is stored modulo.

    sage: a.precision_absolute()
        22

The number of times that p divides the element is called the
valuation, and can be accessed with the functions valuation() and
ordp()

    sage: a.valuation()
        2

The following relationship holds: self.valuation() + self.precision_relative() == self.precision_absolute().
    sage: a.valuation() + a.precision_relative() == a.precision_absolute()
        True

There are two types of p-adic field: capped relative and lazy.  In the
capped relative case, the relative precision of an element is
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

In the lazy case, a certain number of digits are computed and stored,
and in addition a function is stored so that additional digits can be
computed later.  In order to set the number of known digits you call
set_precision_absolute().
    #sage: R = Qp(5, 5, 'lazy', 'series'); a = R(4006); a
    #1  + 5 + 2*5^3 + 5^4 + O(5^5)
    #sage: b = R(50127); b
    #2 + 5^3 + O(5^5)
    #sage: c = a * b; c
    #2 + 2*5 + 4*5^4 + O(5^5)
    #sage: c.set_precision_absolute(15)
    #sage: c
    #2 + 2*5 + 4*5^4 + 3*5^5 + 5^6 + 4*5^8 + 2*5^9 + 4*5^11 + O(5^15)

There is some performance penalty for carrying the function around,
but it is minimized if you determine the precision you will need going
into a computation and set the precision appropriately at the outset.

p-Adic fields should be created using the creation function Qp as
above.  This will ensure that there is only one instance of $\Q_p$ of
a given type, p and precision.  It also saves typing very long class
names.

    sage: Qp(17,10,'capped-rel')
    17-adic Field with capped relative precision 10

    #sage: Qp(7, prec = 30, type = 'lazy', print_mode = 'val-unit')
    #Lazy 7-adic Field
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

In the previous example, since fixed-mod elements don't keep track of
their precision, we assume that it has the full precision of the ring.
This is why you have to cast manually here.

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
import padic_generic
import sage.rings.ring
#import sage.rings.padics.padic_field_capped_relative
#import sage.rings.padics.padic_field_lazy

#pAdicFieldCappedRelative = sage.rings.padics.padic_field_capped_relative.pAdicFieldCappedRelative
#pAdicFieldLazy = sage.rings.padics.padic_field_lazy.pAdicFieldLazy

def is_pAdicField(R):
    return isinstance(R, pAdicFieldGeneric)

class pAdicFieldGeneric(padic_generic.pAdicGeneric, sage.rings.ring.Field):
    def _repr_(self, do_latex = False):
        return "Generic %s-adic Field"%(self.prime())

    def _latex_(self):
        return "\\QQ_{%s}" % self.prime()

    def fraction_field(self, print_mode=None):
        r"""
        Returns the fraction field of self, with given print_mode (defaults to the print_mode of self).

        INPUT:
        print_mode - a dictionary with print options.
        """
        return self

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
        #should be overridden for extension fields
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
