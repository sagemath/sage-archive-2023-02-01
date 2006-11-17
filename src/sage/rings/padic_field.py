r"""
Field $\Q_p$ of $p$-adic Numbers
"""

import weakref

import field
from infinity import infinity
from integer_mod_ring import IntegerModRing
import padic
import integer
import rational

padics = {}
def pAdicField(p, prec=20):
    key = (p, prec)
    if padics.has_key(key):
        K = padics[key]()
        if K != None:
            return K
    K = pAdicField_generic(p, prec)
    padics[key] = weakref.ref(K)
    return K

class pAdicField_generic(field.Field):
    r"""
    Field $\Q_p$ of $p$-adic numbers.

    EXAMPLES:
        sage: K = pAdicField(17); K
        17-adic Field

        sage: loads(K.dumps()) == K
        True
    """
    def __init__(self, p, prec, series_print=True, print_prec=infinity):
        self.__p = p
        self.__prec = prec
        self.__series_print = series_print
        self.__print_prec = print_prec

    def __hash__(self):
        return hash((self.__p, self.__prec))

    def prec(self):
        return self.__prec

    def series_print(self, n=None):
        if n is None:
            return self.__series_print
        self.__series_print = n

    def print_prec(self, n=None):
        """
        If you call print_prec(n), then printing of elements in this
        p-adic field is truncated at $O(p^n)$.   Calling print_prec() with
        no arguments returns n.  This command only affects printing,
        and does not alter the actual values of elements of this field.
        """
        if n==None:
            return self.__print_prec
        self.__print_prec = n

    def characteristic(self):
        """
        The characteristic of the field $\Q_p$, which is always 0.

        EXAMPLES:
            sage: K = Qp(7)
            sage: K.characteristic()
            0
        """
        return 0

    def prime(self):
        """
        The prime p such that this is the field Qp.

        EXAMPLES:
            sage: K = Qp(7)
            sage: K.prime()
            7
        """
        return self.__p

    def residue_characteristic(self):
        """
        The characteristic of the residue class field Qp.

        EXAMPLES:
            sage: K = Qp(7)
            sage: K.residue_characteristic()
            7
        """
        return self.__p

    def residue_class_field(self):
        """
        The residue class field of the ring Zp of integers of Qp, i.e., the field Z/pZ.

        EXAMPLES:
            sage: K = Qp(3)
            sage: K.residue_class_field()
            Finite Field of size 3
        """
        return IntegerModRing(self.__p)

    def _repr_(self):
        return "%s-adic Field"%self.__p

    def __call__(self, x, prec=infinity):
        return padic.pAdic(self, x, prec)

    def _coerce_impl(self, x):
        P = x.parent()
        if is_pAdicField(P) and P.prime() == self.prime() and self.prec() <= P.prec():
            return self(x)
        if isinstance(x, (integer.Integer, rational.Rational)):
            return self(x)
        raise TypeError

    def __cmp__(self, other):
        if not isinstance(other, pAdicField_generic):
            return -1
        return cmp((self.__p, self.__prec), (other.__p, other.__prec))





Qp = pAdicField

def is_pAdicField(x):
    return isinstance(x, pAdicField_generic)

