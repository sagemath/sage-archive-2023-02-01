r"""
Ring $\Z_p$ of $p$-adic Integers
"""

import weakref

import field
from infinity import infinity
from integer_mod_ring import IntegerModRing
import padic
import integer
import rational
#Need to import more/less than this

class pAdicRing_generic(ring.CommutativeRing):
    r"""
    Abstract superclass for $p$-adic integer rings.
    """

#need to implement weak references to make sure that there's exactly one Z_p
class pAdicRing(pAdicRing_generic, ring.DiscreteValuationRing):
    r"""
    EXAMPLES:
        sage: K = pAdicRing(17); K
        17-adic Field

        sage: loads(K.dumps()) == K
        True
    """
    def __init__(self, p, series_print=True, print_prec=infinity):
        self.__p = p
        self.__series_print = series_print
        self.__print_prec = print_prec

    def series_print(self, n=None):
        if n is None:
            return self.__series_print
        self.__series_print = n

    def print_prec(self, n=None):
        """
        If you call print_prec(n), then printing of elements
        is truncated by only printing the unit up to
	$O(p^n)$.   Calling print_prec() with
        no arguments returns n.  This command only affects printing,
        and does not alter the actual values of elements of the ring.
        """
        if n==None:
            return self.__print_prec
        self.__print_prec = n

    def prime(self):
        """
        The prime p such that this is the ring Zp.

        EXAMPLES:
            sage: K = Zp(7)
            sage: K.prime()
            7
        """
        return self.__p

    def residue_characteristic(self):
        """
        The characteristic of the residue class field of Zp.

        EXAMPLES:
            sage: K = Zp(7)
            sage: K.residue_characteristic()
            7
        """
        return self.__p

    def residue_class_field(self):
        """
        The residue class field of the ring Zp, i.e., the field Z/pZ.

        EXAMPLES:
            sage: K = Zp(3)
            sage: K.residue_class_field()
            Finite Field of size 3
        """
        return IntegerModRing(self.__p)

    def _repr_(self):
        return "%s-adic Field"%self.__p

    def __call__(self, x, prec=infinity):
        return padic.pAdic(self, x, prec)

    def _coerce_(self, x):
        if x.parent() is self:
            return x
        if isinstance(x, (integer.Integer, rational.Rational)): #need to add lots more coersion
            return self(x)
        raise TypeError

    def __cmp__(self, other):
        if not isinstance(other, pAdicField_generic):
            return -1
        if self.__p < other.__p:
            return -1
        elif self.__p > other.__p:
            return 1
        return 0





Zp = pAdicRing

def is_pAdicRing(x):
    return isinstance(x, pAdicRing_generic)

