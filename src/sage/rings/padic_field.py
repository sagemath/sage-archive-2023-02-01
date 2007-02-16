r"""
Field $\Q_p$ of $p$-adic Numbers
"""

import weakref

import field
from infinity import infinity
import finite_field
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

    def is_exact(self):
        return False

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

        INPUT:
            n -- integer
                 infinity
                 None

        EXAMPLES:
            sage: K = Qp(13, 6)
            sage: K.print_prec(infinity)
            sage: a = K(1/2); a
            7 + 6*13 + 6*13^2 + 6*13^3 + 6*13^4 + 6*13^5 + O(13^6)
            sage: K.print_prec(3)
            sage: a
            7 + 6*13 + 6*13^2 + O(13^6)
        """
        if n is None:
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
        return finite_field.FiniteField(self.__p)

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
            return cmp(type(self), type(other))
        return cmp((self.__p, self.__prec), (other.__p, other.__prec))


    def teichmuller(self, a):
        """
        INPUT:
            a -- an integer between 1 and p-1 inclusive

        OUTPUT:
            the Teichmuller lift of a as an element of
            this p-adic field, computed to the default
            precision of this field.

        EXAMPLES:
            sage: K = pAdicField(13,prec=6)
            sage: a = K.teichmuller(2); a
            2 + 6*13 + 2*13^2 + 2*13^3 + 4*13^4 + 2*13^5 + O(13^6)
            sage: a^12
            1 + O(13^6)
            sage: a = K.teichmuller(0); a
            Traceback (most recent call last):
            ...
            ValueError: a must be nonzero modulo p
        """
        _a = int(a)
        if _a != a:
            raise TypeError, "a must be coercible to an integer"
        p = self.prime()
        if _a <= 0 or _a >= p:
            _a %= p
        if _a == 0:
            raise ValueError, "a must be nonzero modulo p"
        try:
            return self.__teich[_a]
        except AttributeError:
            pass
        v = [0]*p
        # compute a (p-1)st root of unity in Z_p.
        zeta = self.zeta(p-1)
        z = zeta
        for i in range(p-1):
            j = (z % p).lift()
            v[j] = z
            z *= zeta
        self.__teich = v
        return v[_a]





Qp = pAdicField

def is_pAdicField(x):
    return isinstance(x, pAdicField_generic)

