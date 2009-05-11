"""
p-Adic Base Leaves.

Implementations of Zp and Qp

AUTHORS::

    - David Roe
    - Genya Zaytman: documentation
    - David Harvey: doctests
    - William Stein: doctest updates

EXAMPLES::

    p-Adic rings and fields are examples of inexact structures, as the
    reals are.  That means that elements cannot generally be stored
    exactly: to do so would take an infinite amount of storage.
    Instead, we store an approximation to the elements with varying
    precision.

    There are two types of precision for a p-adic element.  The first
    is relative precision, which gives the number of known p-adic
    digits::

    sage: R = Qp(5, 20, 'capped-rel', 'series'); a = R(675); a
    2*5^2 + 5^4 + O(5^22)
    sage: a.precision_relative()
    20

    The second type of precision is absolute precision, which gives
    the power of p that this element is stored modulo::

    sage: a.precision_absolute()
    22

    The number of times that p divides the element is called the
    valuation, and can be accessed with the functions valuation() and
    ordp():

    sage: a.valuation()
    2

    The following relationship holds: self.valuation() +
    self.precision_relative() == self.precision_absolute().

    sage: a.valuation() + a.precision_relative() == a.precision_absolute()
    True

    In the capped relative case, the relative precision of an element
    is restricted to be at most a certain value, specified at the
    creation of the field.  Individual elements also store their own
    precision, so the effect of various arithmetic operations on
    precision is tracked.  When you cast an exact element into a
    capped relative field, it truncates it to the precision cap of the
    field.::

    sage: R = Qp(5, 5); a = R(4006); a
    1 + 5 + 2*5^3 + 5^4 + O(5^5)
    sage: b = R(17/3); b
    4 + 2*5 + 3*5^2 + 5^3 + 3*5^4 + O(5^5)
    sage: c = R(4025); c
    5^2 + 2*5^3 + 5^4 + 5^5 + O(5^7)
    sage: a + b
    4*5 + 3*5^2 + 3*5^3 + 4*5^4 + O(5^5)
    sage: a + b + c
    4*5 + 4*5^2 + 5^4 + O(5^5)

    ::

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

    In the capped absolute type, instead of having a cap on the
    relative precision of an element there is instead a cap on the
    absolute precision.  Elements still store their own precisions,
    and as with the capped relative case, exact elements are truncated
    when cast into the ring.::

    sage: R = ZpCA(5, 5); a = R(4005); a
    5 + 2*5^3 + 5^4 + O(5^5)
    sage: b = R(4025); b
    5^2 + 2*5^3 + 5^4 + O(5^5)
    sage: a * b
    5^3 + 2*5^4 + O(5^5)
    sage: (a * b) // 5^3
    1 + 2*5 + O(5^2)
    sage: type((a * b) // 5^3)
    <type 'sage.rings.padics.padic_capped_absolute_element.pAdicCappedAbsoluteElement'>
    sage: (a * b) / 5^3
    1 + 2*5 + O(5^2)
    sage: type((a * b) / 5^3)
    <type 'sage.rings.padics.padic_capped_relative_element.pAdicCappedRelativeElement'>

    The fixed modulus type is the leanest of the p-adic rings: it is
    basically just a wrapper around $\Z / p^n \Z$ providing a unified
    interface with the rest of the p-adics.  This is the type you
    should use if your primary interest is in speed (though it's not
    all that much faster than other p-adic types).  It does not track
    precision of elements.::

    sage: R = ZpFM(5, 5); a = R(4005); a
    5 + 2*5^3 + 5^4 + O(5^5)
    sage: a // 5
    1 + 2*5^2 + 5^3 + O(5^5)

    p-Adic rings and fields should be created using the creation functions Zp and Qp as
    above.  This will ensure that there is only one instance of $\Z_p$ and $\Q_p$
    of a given type, p, print mode and precision.  It also saves
    typing very long class names.::

    sage: Qp(17,10)
    17-adic Field with capped relative precision 10
    sage: R = Qp(7, prec = 20, print_mode = 'val-unit'); S = Qp(7, prec = 20, print_mode = 'val-unit'); R is S
    True
    sage: Qp(2)
    2-adic Field with capped relative precision 20

    Once one has a p-Adic ring or field, one can cast elements into it in the
    standard way.  Integers, ints, longs, Rationals, other p-Adic
    types, pari p-adics and elements of $\Z / p^n \Z$ can all be cast
    into a p-Adic field.::

    sage: R = Qp(5, 5, 'capped-rel','series'); a = R(16); a
    1 + 3*5 + O(5^5)
    sage: b = R(23/15); b
    5^-1 + 3 + 3*5 + 5^2 + 3*5^3 + O(5^4)
    sage: S = Zp(5, 5, 'fixed-mod','val-unit'); c = S(Mod(75,125)); c
    5^2 * 3 + O(5^5)
    sage: R(c)
    3*5^2 + O(5^5)

    In the previous example, since fixed-mod elements don't keep track
    of their precision, we assume that it has the full precision of
    the ring.  This is why you have to cast manually here.

    While you can cast explicitly as above, the chains of automatic
    coercion are more restricted.  As always in SAGE, the following
    arrows are transitive and the diagram is commutative.

    int -> long -> Integer -> Zp capped-rel -> Zp capped_abs -> IntegerMod
    Integer -> Zp fixed-mod -> IntegerMod
    Integer -> Zp capped-abs -> Qp capped-rel

    In addition, there are arrows within each type.  For capped
    relative and capped absolute rings and fields, these arrows go
    from lower precision cap to higher precision cap.  This works
    since elements track their own precision: choosing the parent with
    higher precision cap means that precision is less likely to be
    truncated unnecessarily.  For fixed modulus parents, the arrow
    goes from higher precision cap to lower; the fact that elements
    don't track precision necessitates this choice in order to not
    produce incorrect results.

TESTS::

    sage: R = Qp(5, 15, print_mode='bars', print_sep='&'); S = loads(dumps(R))
    sage: R == S
    True
    sage: repr(S(2777))[3:]
    '4&2&1&0&2'

    sage: R = Zp(5, 15, print_mode='bars', print_sep='&'); S = loads(dumps(R))
    sage: R == S
    True
    sage: repr(S(2777))[3:]
    '4&2&1&0&2'

    sage: R = ZpCA(5, 15, print_mode='bars', print_sep='&'); S = loads(dumps(R))
    sage: R == S
    True
    sage: repr(S(2777))[3:]
    '4&2&1&0&2'

"""

from generic_nodes import pAdicFieldBaseGeneric, \
                          pAdicCappedRelativeFieldGeneric, \
                          pAdicRingBaseGeneric, \
                          pAdicCappedRelativeRingGeneric, \
                          pAdicFixedModRingGeneric, \
                          pAdicCappedAbsoluteRingGeneric
from padic_capped_relative_element import pAdicCappedRelativeElement
from padic_capped_absolute_element import pAdicCappedAbsoluteElement
from padic_fixed_mod_element import pAdicFixedModElement
from sage.rings.integer_ring import ZZ

class pAdicRingCappedRelative(pAdicRingBaseGeneric, pAdicCappedRelativeRingGeneric):
    r"""
    An implementation of the p-adic integers with capped relative precision.
    """
    def __init__(self, p, prec, print_mode, names):
        """
        Initialization.

        INPUTS::

            - p -- prime
            - prec -- precision cap
            - print_mode -- dictionary with print options.
            - names -- how to print the prime.

        EXAMPLES::

            sage: R = ZpCR(next_prime(10^60)) #indirect doctest
        """
        pAdicRingBaseGeneric.__init__(self, p, prec, print_mode, names, pAdicCappedRelativeElement)

    def _coerce_map_from_(self, R):
        """
        Returns True if there is a coerce map from R to self.

        EXAMPLES::

            sage: K = Zp(17)
            sage: K(1) + 1 #indirect doctest
            2 + O(17^20)
            sage: K.has_coerce_map_from(ZZ)
            True
            sage: K.has_coerce_map_from(int)
            True
            sage: K.has_coerce_map_from(QQ)
            False
            sage: K.has_coerce_map_from(RR)
            False
            sage: K.has_coerce_map_from(Qp(7))
            False
            sage: K.has_coerce_map_from(Zp(17,40))
            False
            sage: K.has_coerce_map_from(Zp(17,10))
            True
            sage: K.has_coerce_map_from(ZpCA(17,40))
            False
        """
        #if isistance(R, pAdicRingLazy) and R.prime() == self.prime():
        #    return True
        if isinstance(R, pAdicRingCappedRelative) and R.prime() == self.prime():
            if R.precision_cap() < self.precision_cap():
                return True
            elif R.precision_cap() == self.precision_cap() and self._printer.cmp_modes(R._printer) <= 0:
                return True

    def _repr_(self, do_latex=False):
        r"""
        Print representation.

        EXAMPLES::

            sage: K = Zp(17); K #indirect doctest
            17-adic Ring with capped relative precision 20
            sage: latex(K)
            \mathbf{Z}_{17}
        """
        if do_latex:
            return "\\mathbf{Z}_{%s}" % self.prime()
        return "%s-adic Ring with capped relative precision %s"%(self.prime(), self.precision_cap())

class pAdicRingCappedAbsolute(pAdicRingBaseGeneric, pAdicCappedAbsoluteRingGeneric):
    r"""
    An implementation of the p-adic integers with capped absolute precision.
    """
    def __init__(self, p, prec, print_mode, names):
        """
        Initialization.

        INPUTS::

            - p -- prime
            - prec -- precision cap
            - print_mode -- dictionary with print options.
            - names -- how to print the prime.

        EXAMPLES::

            sage: R = ZpCA(next_prime(10^60)) #indirect doctest
        """
        pAdicRingBaseGeneric.__init__(self, p, prec, print_mode, names, pAdicCappedAbsoluteElement)

    def _repr_(self, do_latex = False):
        r"""
        Print representation.

        EXAMPLES::

            sage: K = ZpCA(17); K #indirect doctest
            17-adic Ring with capped absolute precision 20
            sage: latex(K)
            \mathbf{Z}_{17}
        """
        if do_latex:
            return "\\mathbf{Z}_{%s}" % self.prime()
        return "%s-adic Ring with capped absolute precision %s"%(self.prime(), self.precision_cap())

    def _coerce_map_from_(self, R):
        """
        Returns True if there is a coerce map from R to self.

        EXAMPLES::

            sage: K = ZpCA(17)
            sage: K(1) + 1 #indirect doctest
            2 + O(17^20)
            sage: K.has_coerce_map_from(ZZ)
            True
            sage: K.has_coerce_map_from(int)
            True
            sage: K.has_coerce_map_from(QQ)
            False
            sage: K.has_coerce_map_from(RR)
            False
            sage: K.has_coerce_map_from(Qp(7))
            False
            sage: K.has_coerce_map_from(ZpCA(17,40))
            False
            sage: K.has_coerce_map_from(ZpCA(17,10))
            True
            sage: K.has_coerce_map_from(Zp(17,40))
            True
        """
        #if isistance(R, pAdicRingLazy) and R.prime() == self.prime():
        #    return True
        if isinstance(R, pAdicRingCappedRelative) and R.prime() == self.prime():
            return True
        if isinstance(R, pAdicRingCappedAbsolute) and R.prime() == self.prime():
            if R.precision_cap() < self.precision_cap():
                return True
            if R.precision_cap() == self.precision_cap() and self._printer.cmp_modes(R._printer) <= 0:
                return True

class pAdicRingFixedMod(pAdicRingBaseGeneric, pAdicFixedModRingGeneric):
    r"""
    An implementation of the p-adic integers using fixed modulus.
    """
    def __init__(self, p, prec, print_mode, names):
        """
        Initialization

        INPUTS::

            - p -- prime
            - prec -- precision cap
            - print_mode -- dictionary with print options.
            - names -- how to print the prime.

        EXAMPLES::

            sage: R = ZpFM(next_prime(10^60)) #indirect doctest
        """
        pAdicRingBaseGeneric.__init__(self, p, prec, print_mode, names, pAdicFixedModElement)

    def _coerce_map_from_(self, R):
        """
        Returns True if there is a coerce map from R to self.

        EXAMPLES::

            sage: K = ZpFM(17)
            sage: K(1) + 1 #indirect doctest
            2 + O(17^20)
            sage: K.has_coerce_map_from(ZZ)
            True
            sage: K.has_coerce_map_from(int)
            True
            sage: K.has_coerce_map_from(QQ)
            False
            sage: K.has_coerce_map_from(RR)
            False
            sage: K.has_coerce_map_from(Zp(7))
            False
            sage: K.has_coerce_map_from(ZpFM(17,40))
            True
            sage: K.has_coerce_map_from(ZpFM(17,10))
            False
            sage: K.has_coerce_map_from(Zp(17,40))
            False
        """
        #if isistance(R, pAdicRingLazy) and R.prime() == self.prime():
        #    return True
        if isinstance(R, pAdicRingFixedMod) and R.prime() == self.prime():
            if R.precision_cap() > self.precision_cap():
                return True
            if R.precision_cap() == self.precision_cap() and self._printer.cmp_modes(R._printer) <= 0:
                return True

    def _repr_(self, do_latex=False):
        r"""
        Representation

        EXAMPLES::

            sage: K = ZpFM(7); K
            7-adic Ring of fixed modulus 7^20
            sage: latex(K) #indirect doctest
            \mathbf{Z}_{7}
        """
        if do_latex:
            return "\\mathbf{Z}_{%s}" % self.prime()
        return "%s-adic Ring of fixed modulus %s^%s"%(self.prime(), self.prime(), self.precision_cap())

    def fraction_field(self, print_mode = None):
        r"""
        Would normally return $\Q_p$, but there is no implementation
        of $Q_p$ matching this ring so this raises an error

        If you want to be able to divide with elements of a fixed
        modulus p-adic ring, you must cast explicitly.

        EXAMPLES::

            sage: ZpFM(5).fraction_field()
            Traceback (most recent call last):
            ...
            TypeError: This implementation of the p-adic ring does not support fields of fractions.
        """
        raise TypeError, "This implementation of the p-adic ring does not support fields of fractions."


class pAdicFieldCappedRelative(pAdicFieldBaseGeneric, pAdicCappedRelativeFieldGeneric):
    r"""
    An implementation of p-adic fields with capped relative precision.

    EXAMPLES::

        sage: K = Qp(17, 1000000)
        sage: K = Qp(next_prime(10^60))
    """

    def __init__(self, p, prec, print_mode, names):
        """
        Initialization.

        INPUTS::

            - p -- the prime
            - prec -- relative precision cap
            - print_mode -- dictionary of print options.
            - names -- how to print the prime p.

        EXAMPLES::

            sage: K = Qp(101) #indirect doctest
        """
        pAdicFieldBaseGeneric.__init__(self, p, prec, print_mode, names, pAdicCappedRelativeElement)

    def _coerce_map_from_(self, R):
        """
        Returns True if there is a coerce map from R to self.

        EXAMPLES::

            sage: K = Qp(17)
            sage: K(1) + 1 #indirect doctest
            2 + O(17^20)
            sage: K.has_coerce_map_from(ZZ)
            True
            sage: K.has_coerce_map_from(int)
            True
            sage: K.has_coerce_map_from(QQ)
            True
            sage: K.has_coerce_map_from(RR)
            False
            sage: K.has_coerce_map_from(Qp(7))
            False
            sage: K.has_coerce_map_from(Qp(17,40))
            False
            sage: K.has_coerce_map_from(Qp(17,10))
            True
            sage: K.has_coerce_map_from(Zp(17,40))
            True

        """
        #if isinstance(R, pAdicRingLazy) or isinstance(R, pAdicFieldLazy) and R.prime() == self.prime():
        #    return True
        if isinstance(R, (pAdicRingCappedRelative, pAdicRingCappedAbsolute)) and R.prime() == self.prime():
            return True
        if isinstance(R, pAdicFieldCappedRelative) and R.prime() == self.prime():
            if R.precision_cap() < self.precision_cap():
                return True
            elif R.precision_cap() == self.precision_cap() and self._printer.cmp_modes(R._printer) <= 0:
                return True

    def _repr_(self, do_latex=False):
        r"""
        Returns a string representation of self.

        EXAMPLES::

            sage: K = Qp(17); K #indirect doctest
            17-adic Field with capped relative precision 20
            sage: latex(K)
            \mathbf{Q}_{17}
        """
        if do_latex:
            return "\\mathbf{Q}_{%s}" % self.prime()
        return "%s-adic Field with capped relative precision %s"%(self.prime(), self.precision_cap())


    def random_element(self, algorithm='default'):
        r"""
        Returns a random element of self, optionally using the algorithm
        argument to decide how it generates the element. Algorithms currently
        implemented:

            default: Choose an integer k using the standard
              distribution on the integers.  Then choose an integer a
              uniformly in the range 0 <= a < p^self.precision_cap().
              Return self(p^k * a, absprec = k + self.precision_cap()).

        EXAMPLES::

            sage: Qp(17,6).random_element()
            15*17^-8 + 10*17^-7 + 3*17^-6 + 2*17^-5 + 11*17^-4 + 6*17^-3 + O(17^-2)
        """
        if (algorithm == 'default'):
            k = ZZ.random_element()
            a = ZZ.random_element(self.prime()**self.precision_cap())
            return self(self.prime()**k * a, absprec = k + self.precision_cap())
        else:
            raise NotImplementedError, "Don't know %s algorithm"%algorithm
