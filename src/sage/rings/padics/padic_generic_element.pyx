"""
p-Adic Generic Element

Elements of p-Adic Rings and Fields

AUTHOR:

- David Roe
- Genya Zaytman: documentation
- David Harvey: doctests
"""

#*****************************************************************************
#       Copyright (C) 2007 David Roe <roed@math.harvard.edu>
#                          William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include "../../ext/gmp.pxi"
include "../../ext/interrupt.pxi"
include "../../ext/stdsage.pxi"

import sys

cimport sage.rings.padics.local_generic_element
from sage.rings.padics.local_generic_element cimport LocalGenericElement
from sage.rings.rational cimport Rational
#cimport sage.structure.element
#from sage.structure.element cimport Element
#cimport pow_computer
from sage.rings.integer cimport Integer
#from sage.rings.rational import Rational
from sage.rings.infinity import infinity
from sage.libs.pari.gen import pari
from sage.libs.pari.gen import PariError
from sage.rings.rational_field import QQ
from sage.rings.real_mpfr import RR
import sage.rings.rational_field
#from sage.rings.padics.pow_computer cimport PowComputer_base


#Rational = sage.rings.rational.Rational
#infinity = sage.rings.infinity.infinity
#PariError = sage.libs.pari.gen.PariError
#pari = sage.libs.pari.gen.pari
#QQ = sage.rings.rational_field.QQ

cdef long maxordp = (1L << (sizeof(long) * 8 - 2)) - 1

cdef class pAdicGenericElement(LocalGenericElement):
    def __richcmp__(left, right, int op):
        """
        Comparison.

        EXAMPLES::

            sage: R = Zp(5); a = R(5, 6); b = R(5 + 5^6, 8)
            sage: a == b #indirect doctest
            True
        """
        return (<Element>left)._richcmp(right, op)

    cdef int _cmp_c_impl(left, Element right) except -2:
        """
        First compare valuations, then compare normalized
        residue of unit part.

        EXAMPLES::

            sage: R = Zp(19, 5,'capped-rel','series'); K = Qp(19, 5, 'capped-rel', 'series')
            sage: a = R(2); a
            2 + O(19^5)
            sage: b = R(3); b
            3 + O(19^5)
            sage: a < b
            True
            sage: a = K(2); a
            2 + O(19^5)
            sage: b = K(3); b
            3 + O(19^5)
            sage: a < b
            True
        """
        m = min(left.precision_absolute(), right.precision_absolute())
        x_ordp = left.valuation()
        if x_ordp >= m :
            x_ordp = infinity
        y_ordp = right.valuation()
        if y_ordp >= m :
            y_ordp = infinity
        if x_ordp < y_ordp:
            return -1
        elif x_ordp > y_ordp:
            return 1
        else:  # equal ordp
            if x_ordp is infinity:
                return 0 # since both are zero
            else:
                return (<pAdicGenericElement>left.unit_part())._cmp_units(right.unit_part())

    cdef int _cmp_units(left, pAdicGenericElement right) except -2:
        raise NotImplementedError

    cdef int _set_from_Integer(self, Integer x, absprec, relprec) except -1:
        raise NotImplementedError
    cdef int _set_from_mpz(self, mpz_t x) except -1:
        raise NotImplementedError
    cdef int _set_from_mpz_rel(self, mpz_t x, long relprec) except -1:
        raise NotImplementedError
    cdef int _set_from_mpz_abs(self, mpz_t value, long absprec) except -1:
        raise NotImplementedError
    cdef int _set_from_mpz_both(self, mpz_t x, long absprec, long relprec) except -1:
        raise NotImplementedError

    cdef int _set_from_Rational(self, Rational x, absprec, relprec) except -1:
        raise NotImplementedError
    cdef int _set_from_mpq(self, mpq_t x) except -1:
        raise NotImplementedError
    cdef int _set_from_mpq_rel(self, mpq_t x, long relprec) except -1:
        raise NotImplementedError
    cdef int _set_from_mpq_abs(self, mpq_t value, long absprec) except -1:
        raise NotImplementedError
    cdef int _set_from_mpq_both(self, mpq_t x, long absprec, long relprec) except -1:
        raise NotImplementedError

    cdef int _pshift_self(self, long shift) except -1:
        raise NotImplementedError

    cdef int _set_inexact_zero(self, long absprec) except -1:
        raise NotImplementedError
    cdef int _set_exact_zero(self) except -1:
        raise TypeError, "this type of p-adic does not support exact zeros"

    cpdef bint _is_exact_zero(self) except -1:
        """
        Returns True if self is exactly zero.  Since
        non-capped-relative elements cannot be exact, this function
        always returns False.

        EXAMPLES::

            sage: ZpCA(17)(0)._is_exact_zero()
            False
        """
        return False

    cpdef bint _is_inexact_zero(self) except -1:
        """
        Returns True if self is indistinguishable from zero, but not
        exactly zero.

        EXAMPLES::

            sage: Zp(5)(0,5)._is_inexact_zero()
            True
        """
        raise NotImplementedError

    cpdef bint _is_zero_rep(self) except -1:
        """
        Returns True is self is indistinguishable from zero.

        EXAMPLES::

            sage: ZpCA(17)(0,15)._is_zero_rep()
            True
        """
        return self._is_inexact_zero() or self._is_exact_zero()

    cdef bint _set_prec_abs(self, long absprec) except -1:
        self._set_prec_both(absprec, (<PowComputer_class>self.parent().prime_pow).prec_cap)

    cdef bint _set_prec_rel(self, long relprec) except -1:
        self._set_prec_both((<PowComputer_class>self.parent().prime_pow).prec_cap, relprec)

    cdef bint _set_prec_both(self, long absprec, long relprec) except -1:
        return 0

    #def _pari_(self):
    #    """
    #    Returns a pari version of this element.

    #    EXAMPLES::

    #        sage: R = Zp(5)
    #        sage: pari(R(1777))
    #        2 + 5^2 + 4*5^3 + 2*5^4 + O(5^20)
    #    """
    #    return pari(self._pari_init_())

    def __floordiv__(self, right):
        """
        Divides self by right and throws away the nonintegral part if
        self.parent() is not a field.

        There are a number of reasonable definitions for floor
        division.  Any definition should satisfy the following
        identity:

        (1) a = (a // b) * b + a % b

        If a and b lie in a field, then setting a % b = 0 and a // b =
        a / b provides a canonical way of satisfying this equation.

        However, for elements of integer rings, there are many choices
        of definitions for a // b and a % b that satisfy this
        equation.  Since p-adic rings in Sage come equipped with a
        uniformizer pi, we can use the choice of uniformizer in our
        definitions.  Here are some other criteria we might ask for:

        (2) If b = pi^k, the series expansion (in terms of pi) of a //
        b is just the series expansion of a, shifted over by k terms.

        (2') The series expansion of a % pi^k has no terms above
        pi^(k-1).

        The conditions (2) and (2') are equivalent.  But when we
        generalize these conditions to arbitrary b they diverge.

        (3) For general b, the series expansion of a // b is just the
        series expansion of a / b, truncating terms with negative
        exponents of pi.

        (4) For general b, the series expansion of a % b has no terms
        above b.valuation() - 1.

        In order to satisfy (3), one defines

        a // b = (a / b.unit_part()) >> b.valuation()
        a % b = a - (a // b) * b

        In order to satisfy (4), one defines

        a % b = a.lift() % pi.lift()^b.valuation()
        a // b = ((a - a % b) >> b.valuation()) / b.unit_part()


        In Sage we choose option (3), mainly because it is more easily
        defined in terms of shifting and thus generalizes more easily
        to extension rings.

        EXAMPLES::

            sage: R = ZpCA(5); a = R(129378); b = R(2398125)
            sage: a // b #indirect doctest
            3 + 3*5 + 4*5^2 + 2*5^4 + 2*5^6 + 4*5^7 + 5^9 + 5^10 + 5^11 + O(5^12)
            sage: a / b
            4*5^-4 + 3*5^-3 + 2*5^-2 + 5^-1 + 3 + 3*5 + 4*5^2 + 2*5^4 + 2*5^6 + 4*5^7 + 5^9 + 5^10 + 5^11 + O(5^12)
            sage: a % b
            3 + 5^4 + 3*5^5 + 2*5^6 + 4*5^7 + 5^8 + O(5^16)
            sage: (a // b) * b + a % b
            3 + 2*5^4 + 5^5 + 3*5^6 + 5^7 + O(5^16)

            The alternative definition:

            sage: a
            3 + 2*5^4 + 5^5 + 3*5^6 + 5^7 + O(5^20)
            sage: c = ((a - 3)>>4)/b.unit_part(); c
            1 + 2*5 + 2*5^3 + 4*5^4 + 5^6 + 5^7 + 5^8 + 4*5^9 + 2*5^10 + 4*5^11 + 4*5^12 + 2*5^13 + 3*5^14 + O(5^16)
            sage: c*b + 3
            3 + 2*5^4 + 5^5 + 3*5^6 + 5^7 + O(5^20)
        """
        if right == 0:
            raise ZeroDivisionError
        P = self.parent()
        if P.is_field():
            return self / right
        else:
            right = P(right)
            v, u = right.val_unit()
            return P(self / u) >> v


    def __getitem__(self, n):
        r"""
        Returns the coefficient of `p^n` in the series expansion of this
        element, as an integer in the range `0` to `p-1`.

        EXAMPLES::

            sage: R = Zp(7,4,'capped-rel','series'); a = R(1/3); a
            5 + 4*7 + 4*7^2 + 4*7^3 + O(7^4)
            sage: a[0] #indirect doctest
            5
            sage: a[1]
            4

        Negative indices do not have the special meaning they have for regular
        python lists. In the following example, ``a[-1]`` is simply the
        coefficient of `7^{-1}`::

            sage: K = Qp(7,4,'capped-rel')
            sage: b = K(1/7 + 7); b
            7^-1 + 7 + O(7^3)
            sage: b[-2]
            0
            sage: b[-1]
            1
            sage: b[0]
            0
            sage: b[1]
            1
            sage: b[2]
            0

        It is an error to access coefficients which are beyond the precision
        bound::

            sage: b[3]
            Traceback (most recent call last):
            ...
            IndexError: list index out of range
            sage: b[-2]
            0

        Slices also work::

            sage: a[0:2]
            5 + 4*7 + O(7^2)
            sage: a[-1:3:2]
            5 + 4*7^2 + O(7^3)
            sage: b[0:2]
            7 + O(7^2)
            sage: b[-1:3:2]
            7^-1 + 7 + O(7^3)

        If the slice includes coefficients which are beyond the precision
        bound, they are ignored. This is similar to the behaviour of slices of
        python lists::

            sage: a[3:7]
            4*7^3 + O(7^4)
            sage: b[3:7]
            O(7^3)

        For extension elements, "zeros" match the behavior of
        ``list``::

            sage: S.<a> = Qq(125)
            sage: a[-2]
            []

        .. SEEALSO::

            :meth:`sage.rings.padics.local_generic_element.LocalGenericElement.slice`
        """
        if isinstance(n, slice):
            return self.slice(n.start, n.stop, n.step)
        if self.parent().f() == 1:
            zero = Integer(0)
        else:
            zero = []
        if n < self.valuation():
            return zero
        if n >= self.precision_absolute():
            raise IndexError("list index out of range")

        if self.parent().is_field():
            n -= self.valuation()

        # trailing coefficients which are zero are not stored in self.list() -
        # we catch an IndexError to check for this.
        try:
            return self.list()[n]
        except IndexError:
            return zero

    def __invert__(self):
        r"""
        Returns the multiplicative inverse of self.

        EXAMPLE::

            sage: R = Zp(7,4,'capped-rel','series'); a = R(3); a
            3 + O(7^4)
            sage: ~a #indirect doctest
            5 + 4*7 + 4*7^2 + 4*7^3 + O(7^4)

        NOTE: The element returned is an element of the fraction field.
        """
        return self.parent().fraction_field()(self, relprec = self.precision_relative()).__invert__()

    def __mod__(self, right):
        """
        If self is in a field, returns 0.  If in a ring, returns a
        p-adic integer such that

        (1) a = (a // b) * b + a % b

        holds.

        WARNING: The series expansion of a % b continues above the
        valuation of b.

        The definitions of a // b and a % b are intertwined by
        equation (1).  If a and b lie in a field, then setting a % b =
        0 and a // b = a / b provides a canonical way of satisfying
        this equation.

        However, for elements of integer rings, there are many choices
        of definitions for a // b and a % b that satisfy this
        equation.  Since p-adic rings in Sage come equipped with a
        uniformizer pi, we can use the choice of uniformizer in our
        definitions.  Here are some other criteria we might ask for:

        (2) If b = pi^k, the series expansion (in terms of pi) of a //
        b is just the series expansion of a, shifted over by k terms.

        (2') The series expansion of a % pi^k has no terms above
        pi^(k-1).

        The conditions (2) and (2') are equivalent.  But when we
        generalize these conditions to arbitrary b they diverge.

        (3) For general b, the series expansion of a // b is just the
        series expansion of a / b, truncating terms with negative
        exponents of pi.

        (4) For general b, the series expansion of a % b has no terms
        above b.valuation() - 1.

        In order to satisfy (3), one defines

        a // b = (a / b.unit_part()) >> b.valuation()
        a % b = a - (a // b) * b

        In order to satisfy (4), one defines

        a % b = a.lift() % pi.lift()^b.valuation()
        a // b = ((a - a % b) >> b.valuation()) / b.unit_part()


        In Sage we choose option (3), mainly because it is more easily
        defined in terms of shifting and thus generalizes more easily
        to extension rings.

        EXAMPLES::

            sage: R = ZpCA(5); a = R(129378); b = R(2398125)
            sage: a // b
            3 + 3*5 + 4*5^2 + 2*5^4 + 2*5^6 + 4*5^7 + 5^9 + 5^10 + 5^11 + O(5^12)
            sage: a / b
            4*5^-4 + 3*5^-3 + 2*5^-2 + 5^-1 + 3 + 3*5 + 4*5^2 + 2*5^4 + 2*5^6 + 4*5^7 + 5^9 + 5^10 + 5^11 + O(5^12)
            sage: a % b #indirect doctest
            3 + 5^4 + 3*5^5 + 2*5^6 + 4*5^7 + 5^8 + O(5^16)

            The alternative definition:

            sage: a
            3 + 2*5^4 + 5^5 + 3*5^6 + 5^7 + O(5^20)
            sage: c = ((a - 3)>>4)/b.unit_part(); c
            1 + 2*5 + 2*5^3 + 4*5^4 + 5^6 + 5^7 + 5^8 + 4*5^9 + 2*5^10 + 4*5^11 + 4*5^12 + 2*5^13 + 3*5^14 + O(5^16)
            sage: c*b + 3
            3 + 2*5^4 + 5^5 + 3*5^6 + 5^7 + O(5^20)
        """
        if right == 0:
            raise ZeroDivisionError
        if self.parent().is_field():
            return self.parent()(0)
        else:
            return self - (self // right) * right

    #def _is_exact_zero(self):
    #    return False

    #def _is_inexact_zero(self):
    #    return self.is_zero() and not self._is_exact_zero()

    def str(self, mode=None):
        """
        Returns a string representation of self.

        EXAMPLES::

            sage: Zp(5,5,print_mode='bars')(1/3).str()[3:]
            '1|3|1|3|2'
        """
        return self._repr_(mode=mode)

    def _repr_(self, mode=None, do_latex=False):
        """
        Returns a string representation of self.

        INPUTS:

        - ``mode`` -- allows one to override the default print mode of
          the parent (default: ``None``).

        - ``do_latex`` -- whether to return a latex representation or
          a normal one.

        EXAMPLES::

            sage: Zp(5,5)(1/3) # indirect doctest
            2 + 3*5 + 5^2 + 3*5^3 + 5^4 + O(5^5)
        """
        return self.parent()._printer.repr_gen(self, do_latex, mode=mode)

    def additive_order(self, prec):
        r"""
        Returns the additive order of self, where self is considered
        to be zero if it is zero modulo `p^{\mbox{prec}}`.

        INPUT:

        - ``self`` -- a p-adic element
        - ``prec`` -- an integer

        OUTPUT:

        integer -- the additive order of self

        EXAMPLES::

            sage: R = Zp(7, 4, 'capped-rel', 'series'); a = R(7^3); a.additive_order(3)
            1
            sage: a.additive_order(4)
            +Infinity
            sage: R = Zp(7, 4, 'fixed-mod', 'series'); a = R(7^5); a.additive_order(6)
            1
        """
        if self.is_zero(prec):
            return Integer(1)
        else:
            return infinity

    def algdep(self, n):
        """
        Returns a polynomial of degree at most `n` which is approximately
        satisfied by this number. Note that the returned polynomial need not be
        irreducible, and indeed usually won't be if this number is a good
        approximation to an algebraic number of degree less than `n`.

        ALGORITHM: Uses the PARI C-library ``algdep`` command.

        INPUT:

        - ``self`` -- a p-adic element
        - ``n`` -- an integer

        OUTPUT:

        polynomial -- degree n polynomial approximately satisfied by self

        EXAMPLES::

            sage: K = Qp(3,20,'capped-rel','series'); R = Zp(3,20,'capped-rel','series')
            sage: a = K(7/19); a
            1 + 2*3 + 3^2 + 3^3 + 2*3^4 + 2*3^5 + 3^8 + 2*3^9 + 3^11 + 3^12 + 2*3^15 + 2*3^16 + 3^17 + 2*3^19 + O(3^20)
            sage: a.algdep(1)
            19*x - 7
            sage: K2 = Qp(7,20,'capped-rel')
            sage: b = K2.zeta(); b.algdep(2)
            x^2 - x + 1
            sage: K2 = Qp(11,20,'capped-rel')
            sage: b = K2.zeta(); b.algdep(4)
            x^4 - x^3 + x^2 - x + 1
            sage: a = R(7/19); a
            1 + 2*3 + 3^2 + 3^3 + 2*3^4 + 2*3^5 + 3^8 + 2*3^9 + 3^11 + 3^12 + 2*3^15 + 2*3^16 + 3^17 + 2*3^19 + O(3^20)
            sage: a.algdep(1)
            19*x - 7
            sage: R2 = Zp(7,20,'capped-rel')
            sage: b = R2.zeta(); b.algdep(2)
            x^2 - x + 1
            sage: R2 = Zp(11,20,'capped-rel')
            sage: b = R2.zeta(); b.algdep(4)
            x^4 - x^3 + x^2 - x + 1
        """
        # TODO: figure out if this works for extension rings.  If not, move this to padic_base_generic_element.
        return sage.rings.arith.algdep(self, n)

    def algebraic_dependency(self, n):
        """
        Returns a polynomial of degree at most `n` which is approximately
        satisfied by this number.  Note that the returned polynomial need not
        be irreducible, and indeed usually won't be if this number is a good
        approximation to an algebraic number of degree less than `n`.

        ALGORITHM: Uses the PARI C-library algdep command.

        INPUT:

        - ``self`` -- a p-adic element
        - ``n`` -- an integer

        OUTPUT:

        polynomial -- degree n polynomial approximately satisfied by self

        EXAMPLES::

            sage: K = Qp(3,20,'capped-rel','series'); R = Zp(3,20,'capped-rel','series')
            sage: a = K(7/19); a
            1 + 2*3 + 3^2 + 3^3 + 2*3^4 + 2*3^5 + 3^8 + 2*3^9 + 3^11 + 3^12 + 2*3^15 + 2*3^16 + 3^17 + 2*3^19 + O(3^20)
            sage: a.algebraic_dependency(1)
            19*x - 7
            sage: K2 = Qp(7,20,'capped-rel')
            sage: b = K2.zeta(); b.algebraic_dependency(2)
            x^2 - x + 1
            sage: K2 = Qp(11,20,'capped-rel')
            sage: b = K2.zeta(); b.algebraic_dependency(4)
            x^4 - x^3 + x^2 - x + 1
            sage: a = R(7/19); a
            1 + 2*3 + 3^2 + 3^3 + 2*3^4 + 2*3^5 + 3^8 + 2*3^9 + 3^11 + 3^12 + 2*3^15 + 2*3^16 + 3^17 + 2*3^19 + O(3^20)
            sage: a.algebraic_dependency(1)
            19*x - 7
            sage: R2 = Zp(7,20,'capped-rel')
            sage: b = R2.zeta(); b.algebraic_dependency(2)
            x^2 - x + 1
            sage: R2 = Zp(11,20,'capped-rel')
            sage: b = R2.zeta(); b.algebraic_dependency(4)
            x^4 - x^3 + x^2 - x + 1
        """
        return self.algdep(n)

    #def exp_artin_hasse(self):
    #    """
    #    Returns the Artin-Hasse exponential of self.

    #    This is defined by: E_p(x) = exp(x + x^p/p + x^(p^2)/p^2 + ...)
    #    """
    #    raise NotImplementedError

    #def gamma(self):
    #    raise NotImplementedError

    def is_square(self): #should be overridden for lazy elements
        """
        Returns whether self is a square

        INPUT:

        - ``self`` -- a p-adic element

        OUTPUT:

        boolean -- whether self is a square

        EXAMPLES::

            sage: R = Zp(3,20,'capped-rel')
            sage: R(0).is_square()
            True
            sage: R(1).is_square()
            True
            sage: R(2).is_square()
            False

        TESTS::

            sage: R(3).is_square()
            False
            sage: R(4).is_square()
            True
            sage: R(6).is_square()
            False
            sage: R(9).is_square()
            True

            sage: R2 = Zp(2,20,'capped-rel')
            sage: R2(0).is_square()
            True
            sage: R2(1).is_square()
            True
            sage: R2(2).is_square()
            False
            sage: R2(3).is_square()
            False
            sage: R2(4).is_square()
            True
            sage: R2(5).is_square()
            False
            sage: R2(6).is_square()
            False
            sage: R2(7).is_square()
            False
            sage: R2(8).is_square()
            False
            sage: R2(9).is_square()
            True

            sage: K = Qp(3,20,'capped-rel')
            sage: K(0).is_square()
            True
            sage: K(1).is_square()
            True
            sage: K(2).is_square()
            False
            sage: K(3).is_square()
            False
            sage: K(4).is_square()
            True
            sage: K(6).is_square()
            False
            sage: K(9).is_square()
            True
            sage: K(1/3).is_square()
            False
            sage: K(1/9).is_square()
            True

            sage: K2 = Qp(2,20,'capped-rel')
            sage: K2(0).is_square()
            True
            sage: K2(1).is_square()
            True
            sage: K2(2).is_square()
            False
            sage: K2(3).is_square()
            False
            sage: K2(4).is_square()
            True
            sage: K2(5).is_square()
            False
            sage: K2(6).is_square()
            False
            sage: K2(7).is_square()
            False
            sage: K2(8).is_square()
            False
            sage: K2(9).is_square()
            True
            sage: K2(1/2).is_square()
            False
            sage: K2(1/4).is_square()
            True
       """
        if self._is_exact_zero() or self._is_inexact_zero():
            return True
        elif self.parent().prime() != 2:
            return (self.valuation() % 2 == 0) and (self.unit_part().residue(1).is_square())
        else:
            #won't work for general extensions...
            return (self.valuation() % 2 == 0) and (self.unit_part().residue(3) == 1)

    #def log_artin_hasse(self):
    #    raise NotImplementedError

    def multiplicative_order(self, prec = None): #needs to be rewritten for lazy elements
        r"""
        Returns the multiplicative order of self, where self is
        considered to be one if it is one modulo `p^{\mbox{prec}}`.

        INPUT:

        - ``self`` -- a p-adic element
        - ``prec`` -- an integer

        OUTPUT:

        - integer -- the multiplicative order of self

        EXAMPLES::

            sage: K = Qp(5,20,'capped-rel')
            sage: K(-1).multiplicative_order(20)
            2
            sage: K(1).multiplicative_order(20)
            1
            sage: K(2).multiplicative_order(20)
            +Infinity
            sage: K(3).multiplicative_order(20)
            +Infinity
            sage: K(4).multiplicative_order(20)
            +Infinity
            sage: K(5).multiplicative_order(20)
            +Infinity
            sage: K(25).multiplicative_order(20)
            +Infinity
            sage: K(1/5).multiplicative_order(20)
            +Infinity
            sage: K(1/25).multiplicative_order(20)
            +Infinity
            sage: K.zeta().multiplicative_order(20)
            4

            sage: R = Zp(5,20,'capped-rel')
            sage: R(-1).multiplicative_order(20)
            2
            sage: R(1).multiplicative_order(20)
            1
            sage: R(2).multiplicative_order(20)
            +Infinity
            sage: R(3).multiplicative_order(20)
            +Infinity
            sage: R(4).multiplicative_order(20)
            +Infinity
            sage: R(5).multiplicative_order(20)
            +Infinity
            sage: R(25).multiplicative_order(20)
            +Infinity
            sage: R.zeta().multiplicative_order(20)
            4
        """
        if self.valuation() != 0:
            return infinity
        res = self.residue(1)
        if self.is_equal_to(self.parent().teichmuller(res.lift()),prec): #should this be made more efficient?
            return res.multiplicative_order()
        else:
            return infinity

    def valuation(self, p = None):
        """
        Returns the valuation of self.

        INPUT:

        - ``self`` -- a p-adic element
        - ``p`` -- a prime (default: None). If specified, will make sure that p==self.parent().prime()

        NOTE: The optional argument p is used for consistency with the valuation methods on integer and rational.

        OUTPUT:

        integer -- the valuation of self

        EXAMPLES::

            sage: R = Zp(17, 4,'capped-rel')
            sage: a = R(2*17^2)
            sage: a.valuation()
            2
            sage: R = Zp(5, 4,'capped-rel')
            sage: R(0).valuation()
            +Infinity

        TESTS::

            sage: R(1).valuation()
            0
            sage: R(2).valuation()
            0
            sage: R(5).valuation()
            1
            sage: R(10).valuation()
            1
            sage: R(25).valuation()
            2
            sage: R(50).valuation()
            2
            sage: R = Qp(17, 4)
            sage: a = R(2*17^2)
            sage: a.valuation()
            2
            sage: R = Qp(5, 4)
            sage: R(0).valuation()
            +Infinity
            sage: R(1).valuation()
            0
            sage: R(2).valuation()
            0
            sage: R(5).valuation()
            1
            sage: R(10).valuation()
            1
            sage: R(25).valuation()
            2
            sage: R(50).valuation()
            2
            sage: R(1/2).valuation()
            0
            sage: R(1/5).valuation()
            -1
            sage: R(1/10).valuation()
            -1
            sage: R(1/25).valuation()
            -2
            sage: R(1/50).valuation()
            -2

            sage: K.<a> = Qq(25)
            sage: K(0).valuation()
            +Infinity

            sage: R(1/50).valuation(5)
            -2
            sage: R(1/50).valuation(3)
            Traceback (most recent call last):
            ...
            ValueError: Ring (5-adic Field with capped relative precision 4) residue field of the wrong characteristic.
        """
        if not p is None and p != self.parent().prime():
            raise ValueError, 'Ring (%s) residue field of the wrong characteristic.'%self.parent()
        cdef long v = self.valuation_c()
        if v == maxordp:
            return infinity
        cdef Integer ans = PY_NEW(Integer)
        mpz_set_si(ans.value, v)
        return ans

    cdef long valuation_c(self):
        """
        This function is overridden in subclasses to provide an
        actual implementation of valuation.

        For exact zeros, ``maxordp`` is returned, rather than infinity.

        EXAMPLES:

        For example, the valuation function on pAdicCappedRelativeElements
        uses an overridden version of this function.

        ::

            sage: Zp(5)(5).valuation() #indirect doctest
            1
        """
        raise NotImplementedError

    cpdef val_unit(self):
        """
        Return ``(self.valuation(), self.unit_part())``. To be overridden in
        derived classes.

        EXAMPLES::

            sage: Zp(5,5)(5).val_unit()
            (1, 1 + O(5^5))
        """
        raise NotImplementedError

    def ordp(self, p = None):
        r"""
        Returns the valuation of self, normalized so that the valuation of `p` is 1

        INPUT:

        - ``self`` -- a p-adic element
        - ``p`` -- a prime (default: ``None``). If specified, will make sure that ``p == self.parent().prime()``

        NOTE: The optional argument p is used for consistency with the valuation methods on integer and rational.


        OUTPUT:

        integer -- the valuation of self, normalized so that the valuation of `p` is 1

        EXAMPLES::

            sage: R = Zp(5,20,'capped-rel')
            sage: R(0).ordp()
            +Infinity
            sage: R(1).ordp()
            0
            sage: R(2).ordp()
            0
            sage: R(5).ordp()
            1
            sage: R(10).ordp()
            1
            sage: R(25).ordp()
            2
            sage: R(50).ordp()
            2
            sage: R(1/2).ordp()
            0
        """
        return self.valuation(p) / self.parent().ramification_index()

    def rational_reconstruction(self):
        r"""
        Returns a rational approximation to this p-adic number

        INPUT:

        - ``self`` -- a p-adic element

        OUTPUT:

        rational -- an approximation to self

        EXAMPLES::

            sage: R = Zp(5,20,'capped-rel')
            sage: for i in range(11):
            ...       for j in range(1,10):
            ...           if j == 5:
            ...               continue
            ...           assert i/j == R(i/j).rational_reconstruction()
        """
        if self.is_zero(self.precision_absolute()):
            return Rational(0)
        p = self.parent().prime()
        alpha = self.unit_part().lift()
        m = Integer(p**self.precision_relative())
        r = sage.rings.arith.rational_reconstruction(alpha, m)
        return (Rational(p)**self.valuation())*r

    def _shifted_log(self, prec, check = True):
        r"""
        Returns `-\log(1-x)`.

        INPUT:

        - ``self`` -- a p-adic element
        - ``prec`` -- integer. Adds terms of the series only up to
          the indicated power of precision `p^\mbox{prec}`

        OUTPUT:

        - a p-adic element. The value of `-\log(1-x)` correct up to `p^\mbox{prec}`.

        EXAMPLE::

            sage: r = Qp(5,prec=4)(5)
            sage: r._shifted_log(2)
            5 + 3*5^2 + 2*5^3 + 2*5^4 + O(5^5)
            sage: r._shifted_log(4)
            5 + 3*5^2 + 4*5^3 + 4*5^4 + O(5^5)

        """
        alpha=self.valuation()
        if check and alpha<=0:
            raise ValueError, 'Input value (=%s) should have strictly positive valuation' % self

        if hasattr(self.parent(),'_rat_inverses'):
            inv_integers=self.parent()._rat_inverses
        else:
            inv_integers= [QQ(0)]+[1/QQ(ii) for ii in range(1,self.parent().precision_cap()+1)]
            self.parent()._rat_inverses=inv_integers
        e=self.parent().ramification_index()
        p=self.parent().prime()

        if hasattr(self,'lift'):
            x=QQ(self.lift())
        else:
            x=self

        # Compute the power series using Horner's method
        # One needs to add n2 = floor( prec/alpha ) terms for sure,
        # but some extra terms are needed depending of what n2 is.
        # At the end of the first loop we calculate the other
        # coefficients that contribute up to the requested precision.

        total=x.parent()(0)
        n1=1
        k=0
        p2k=1
        n2=p2k*((prec+k*e)/(p2k*alpha)).floor()
        # The while block gets executed at least once, and usually only
        # once. Sometimes a few more coefficients are needed to get the
        # requested precision.
        while n1 <= n2:
            new=0
            xx=x**(p2k)
            for n in range(n2,n1-1,-p2k):
                new*=xx
                try:
                    new += inv_integers[n]
                except IndexError:
                    new +1/QQ(n)
            # The integer n is not necessarily n1 in subsequent loops,
            # since n decreases by -p2k at each iteration.
            new*=x**(n-1)
            n1=n2+1
            total+=new
            k+=1
            p2k*=p
            # Compute the next range of coefficients
            n2=p2k*((prec+k*e)/(p2k*alpha)).floor()
        try:
            if hasattr(self,'add_bigoh'):
                return self.parent()(x*total).add_bigoh(self.precision_absolute())
            else:
                return self.parent()(x*total)
        except:
            return x*total



    def log(self, branch = None,  prec = None):
        r"""
        Compute the `p`-adic logarithm of any unit in `\ZZ_p`. (See below for
        normalization.)

        The usual power series for log with values in the additive group of
        `\ZZ_p` only converges for 1-units (units congruent to 1 modulo `p`).
        However, there is a unique extension of log to a homomorphism defined
        on all the units. If `u = a \cdot v` is a unit with `v \equiv 1
        \pmod{p}` and `a` a Teichmuller representative, then we define `log(u)
        = log(v)`.  This is the correct extension because the units `U` of
        `\ZZ_p` split as a product `U = V \times \langle w \rangle`, where `V`
        is the subgroup of 1-units and `w` is a `(p-1)`-st root of unity.  The
        `\langle w \rangle` factor is torsion, so must go to 0 under any
        homomorphism to the torsion free group `(\ZZ_p, +)`.

        INPUTS:

        - ``self`` -- a `p`-adic element
        - ``branch`` -- A choice of branch, i.e. a choice of logarithm of the
          uniformizer.  This choice can be made arbitrarily.
        - ``prec`` -- integer (default: None) if specified, computes only up to
          the indicated absolute precision indicated as a power of `p`. This parameter
          is at most the precision cap of the parent, and will be ignored if larger.
          The precision of ``self`` is left unchanged.


        NOTES:

        What some other systems do:

        * PARI:  Seems to define log the same way as we do.

        * MAGMA: Gives an error when unit is not a 1-unit.

        ALGORITHM:

        Input: Some p-adic unit `u` (or non-unit if ``branch`` is
        specified).

        1. Check that the input `p`-adic number is really a unit
           (i.e., valuation 0), or take the unit part and multiply by
           ``branch * self.valuation()`` if not and ``branch``
           specified.

        2. Raise `u` to `q-1` where `q` is the degree of the ring extension, to obtain a one-unit.

        3. Use the series expansion

        .. math::

            \log(1-x) = F(x) = -x - 1/2 x^2 - 1/3 x^3 - 1/4 x^4 - 1/5 x^5 - \cdots

        to compute the logarithm `log(u)`.


        EXAMPLES::

            sage: Z13 = Zp(13, 10, print_mode='series')
            sage: a = Z13(14); a
            1 + 13 + O(13^10)

        Note that the relative precision decreases when we take log -- it is
        the absolute precision that is preserved::

            sage: a.log()
            13 + 6*13^2 + 2*13^3 + 5*13^4 + 10*13^6 + 13^7 + 11*13^8 + 8*13^9 + O(13^10)
            sage: Q13 = Qp(13, 10, print_mode='series')
            sage: a = Q13(14); a
            1 + 13 + O(13^10)
            sage: a.log()
            13 + 6*13^2 + 2*13^3 + 5*13^4 + 10*13^6 + 13^7 + 11*13^8 + 8*13^9 + O(13^10)

        The next few examples illustrate precision when computing `p`-adic
        logs. First we create a field with *default* precision 10::

            sage: R = Zp(5,10, print_mode='series')
            sage: e = R(389); e
            4 + 2*5 + 3*5^3 + O(5^10)
            sage: e.log()
            2*5 + 2*5^2 + 4*5^3 + 3*5^4 + 5^5 + 3*5^7 + 2*5^8 + 4*5^9 + O(5^10)
            sage: K = Qp(5,10, print_mode='series')
            sage: e = K(389); e
            4 + 2*5 + 3*5^3 + O(5^10)
            sage: e.log()
            2*5 + 2*5^2 + 4*5^3 + 3*5^4 + 5^5 + 3*5^7 + 2*5^8 + 4*5^9 + O(5^10)

        Check that results are consistent over a range of precision::

            sage: max_prec = 40
            sage: p = 3
            sage: K = Zp(p, max_prec)
            sage: full_log = (K(1 + p)).log()
            sage: for prec in range(2, max_prec):
            ...       ll1 = (K(1 + p).add_bigoh(prec)).log()
            ...       ll2 = K(1+p).log(prec)
            ...       assert ll1 == full_log
            ...       assert ll2 == full_log
            ...       assert ll1.precision_absolute() == prec

        First, the Eisenstein case::

            sage: R = ZpCR(5,5)
            sage: S.<x> = R[]
            sage: f = x^4 + 15*x^2 + 625*x - 5
            sage: W.<w> = R.ext(f)
            sage: z = 1 + w^2 + 4*w^7; z
            1 + w^2 + 4*w^7 + O(w^20)
            sage: z.log()
            w^2 + 2*w^4 + 3*w^6 + 4*w^7 + w^9 + 4*w^10 + 4*w^11 + 4*w^12 + 3*w^14 + w^15 + w^17 + 3*w^18 + 3*w^19 + O(w^20)

        Check that log is multiplicative::

            sage: y = 1 + 3*w^4 + w^5
            sage: y.log() + z.log() - (y*z).log()
            O(w^20)

        Now an unramified example::

            sage: g = x^3 + 3*x + 3
            sage: A.<a> = R.ext(g)
            sage: b = 1 + 5*(1 + a^2) + 5^3*(3 + 2*a)
            sage: b.log()
            (a^2 + 1)*5 + (3*a^2 + 4*a + 2)*5^2 + (3*a^2 + 2*a)*5^3 + (3*a^2 + 2*a + 2)*5^4 + O(5^5)

        Check that log is multiplicative::

            sage: c = 3 + 5^2*(2 + 4*a)
            sage: b.log() + c.log() - (b*c).log()
            O(5^5)

        We illustrate the effect of the precision argument::

            sage: R = ZpCA(7,10)
            sage: x = R(123456789)/R(3); x
            5 + 3*7^2 + 4*7^3 + 3*7^4 + 5*7^5 + 6*7^6 + 7^9 + O(7^10)
            sage: x.log(prec = 5)
            7 + 3*7^2 + 4*7^3 + 3*7^4 + O(7^5)
            sage: x.log(prec = 7)
            7 + 3*7^2 + 4*7^3 + 3*7^4 + 7^5 + 3*7^6 + O(7^7)
            sage: x.log()
            7 + 3*7^2 + 4*7^3 + 3*7^4 + 7^5 + 3*7^6 + 7^7 + 3*7^8 + 4*7^9 + O(7^10)

        AUTHORS:

        - William Stein: initial version
        - David Harvey (2006-09-13): corrected subtle precision bug
          (need to take denominators into account! -- see trac \#53)
        - Genya Zaytman (2007-02-14): adapted to new p-adic class
        - Amnon Besser, Marc Masdeu (2012-02-21): complete rewrite, valid for
          generic p-adic rings.
        """
        if self.is_zero():
            return ValueError,'Log is not defined at zero'

        e=self.parent().ramification_index()
        p=self.parent().prime()
        q=p**self.parent().degree()
        if not self.is_padic_unit():
            if branch is None:
                raise ValueError,'Not a unit: please specify a branch'
            if e > 1:
                try:
                    u_log=self.parent()._log_u
                except AttributeError:
                    u_log=-((self.parent().uniformizer()**e)/p).log(prec = self.parent().precision_cap())
                    self.parent()._log_u=u_log
                unif_log=(branch-u_log)/e
            else:
                unif_log=branch
            y=self.unit_part()
            total=unif_log*self.valuation()
        else:
            y=self
            total=self.parent()(0)
        if (y-1).valuation()>0:
            denom=Integer(1)
        else:
            y=y**(q-1)
            denom=Integer(q-1)

        # if prec is None or prec > self.parent().precision_cap():
        #     prec=self.parent().precision_cap()
        if prec is None or prec > self.precision_absolute():
            prec=self.precision_absolute()


        # Need to be able to carry computations in higher precision.
        # Otherwise this will lose precision, so we deactivate it for now.
        if False:
            k=0
            y1 = y**p
            while (y1-1).valuation()>(y-1).valuation():
                y=y1
                k+=1
                denom*=p
                y1 = y**p
            return total-((1-y)._shifted_log(prec = prec+1+k, check = False))/denom

        if hasattr(self,'add_bigoh'):
            return (total-((1-y)._shifted_log(prec = prec+1, check = False))/denom).add_bigoh(prec)
        else:
            return total-((1-y)._shifted_log(prec = prec+1, check = False))/denom



    def exp(self, prec = None):
        r"""
        Compute the `p`-adic exponential of any element of `\ZZ_p` where the
        series converges.

        INPUT:

        - ``self`` -- a `p`-adic number
        - ``prec`` -- integer (default: None) if specified, computes only up to
          the indicated precision.

        ALGORITHM: If self has a ``lift`` method (which should happen for
        elements of `\QQ_p` and `\ZZ_p`), then one uses the rule:
        `\exp(x)=\exp(p)^{x/p}` modulo the precision. The value of `\exp(p)` is
        precomputed. Otherwise, use the power series expansion of `\exp`,
        evaluating a certain number of terms which does about `O(prec)`
        multiplications.

        EXAMPLES:

        Borrowed from log::

            sage: Z13 = Zp(13, 10, print_mode='series')
            sage: a = Z13(13 + 6*13**2 + 2*13**3 + 5*13**4 + 10*13**6 + 13**7 + 11*13**8 + 8*13**9).add_bigoh(10); a
            13 + 6*13^2 + 2*13^3 + 5*13^4 + 10*13^6 + 13^7 + 11*13^8 + 8*13^9 + O(13^10)
            sage: a.exp()
            1 + 13 + O(13^10)
            sage: Q13 = Qp(13, 10, print_mode='series')
            sage: a = Q13(13 + 6*13**2 + 2*13**3 + 5*13**4 + 10*13**6 + 13**7 + 11*13**8 + 8*13**9).add_bigoh(10); a
            13 + 6*13^2 + 2*13^3 + 5*13^4 + 10*13^6 + 13^7 + 11*13^8 + 8*13^9 + O(13^10)
            sage: a.exp()
            1 + 13 + O(13^10)

        The next few examples illustrate precision when computing `p`-adic
        exps. First we create a field with *default* precision 10::

            sage: R = Zp(5,10, print_mode='series')
            sage: e = R(2*5 + 2*5**2 + 4*5**3 + 3*5**4 + 5**5 + 3*5**7 + 2*5**8 + 4*5**9).add_bigoh(10); e
            2*5 + 2*5^2 + 4*5^3 + 3*5^4 + 5^5 + 3*5^7 + 2*5^8 + 4*5^9 + O(5^10)
            sage: e.exp()*R.teichmuller(4)
            4 + 2*5 + 3*5^3 + O(5^10)

        ::

            sage: K = Qp(5,10, print_mode='series')
            sage: e = K(2*5 + 2*5**2 + 4*5**3 + 3*5**4 + 5**5 + 3*5**7 + 2*5**8 + 4*5**9).add_bigoh(10); e
            2*5 + 2*5^2 + 4*5^3 + 3*5^4 + 5^5 + 3*5^7 + 2*5^8 + 4*5^9 + O(5^10)
            sage: e.exp()*K.teichmuller(4)
            4 + 2*5 + 3*5^3 + O(5^10)

        TESTS:

        Check that results are consistent over a range of precision::

            sage: max_prec = 40
            sage: p = 3
            sage: K = Zp(p, max_prec)
            sage: full_exp = (K(p)).exp()
            sage: for prec in range(2, max_prec):
            ...       ll = (K(p).add_bigoh(prec)).exp()
            ...       assert ll == full_exp
            ...       assert ll.precision_absolute() == prec
            sage: K = Qp(p, max_prec)
            sage: full_exp = (K(p)).exp()
            sage: for prec in range(2, max_prec):
            ...       ll = (K(p).add_bigoh(prec)).exp()
            ...       assert ll == full_exp
            ...       assert ll.precision_absolute() == prec

        AUTHORS:

        - Genya Zaytman (2007-02-15)
        - Amnon Besser, Marc Masdeu (2012-02-23): Complete rewrite
        """
        p=self.parent().prime()
        if self.valuation()<=QQ(1)/(p-1):
            raise ValueError, 'Exponential does not converge for that input.'
        if prec is None or prec > self.parent().precision_cap():
            prec=self.parent().precision_cap()

        if hasattr(self,'lift'):
            y=QQ(self.lift())
            if p == 2: p=4
            if hasattr(self.parent(),'_exp_p'):
                exp_p=self.parent()._exp_p
            else:
                exp_p=self.parent()(p)._exp(self.parent().precision_cap())
                self.parent()._exp_p=exp_p
            return (exp_p**Integer(y/p)).add_bigoh(self.precision_absolute())
        return self._exp(prec)

    def _exp(self, prec):
        r"""
        Does the actual computation of the exponential power series, using
        Horner's evaluation and only one division.

        INPUT:

        - ``self`` -- a `p`-adic element
        - ``prec`` -- the precison to which to compute the exponential

        EXAMPLES::

            sage: R.<w> = Zq(7^2,5)
            sage: x = R(7*w)
            sage: x._exp(5)
            1 + w*7 + (4*w + 2)*7^2 + (w + 6)*7^3 + 5*7^4 + O(7^5)

        AUTHORS:

        - Genya Zaytman (2007-02-15)
        - Amnon Besser, Marc Masdeu (2012-02-23): Complete rewrite
        """
        p=self.parent().prime()
        e=self.parent().ramification_index()
        val=self.valuation()
        # Copied the bound from a previous implementation
        cdef long max_term = ((p-1)*(prec-1))//((p-1)*val - 1) + 1
        cdef long n
        if hasattr(self,'lift'):
            x=QQ(self.lift())
        else:
            x=self
        term=x.parent()(1)
        factor=Integer(1)
        for n in range(max_term,0,-1):
            factor *= n
            term *= x
            term += factor
        if not term.parent() is QQ:
            if e ==1 :
                vfact=factor.valuation(p)
            else:
                vfact=self.parent()(factor).ordp()
            fact1=factor/(p**vfact)
            term1=(term<<(-e*vfact)).lift_to_precision(self.parent().precision_cap())
            return term1 / fact1
        else:
            fact1=factor
            term1=term
            return self.parent()(term1/fact1)

    def square_root(self, extend = True, all = False):
        r"""
        Returns the square root of this p-adic number

        INPUT:

        - ``self`` -- a p-adic element
        - ``extend`` -- bool (default: True); if True, return a square root in
          an extension if necessary; if False and no root exists in the given
          ring or field, raise a ValueError
        - ``all`` -- bool (default: False); if True, return a list of all
          square roots

        OUTPUT:

        p-adic element -- the square root of this p-adic number

        If ``all=False``, the square root chosen is the one whose
        reduction mod `p` is in the range `[0, p/2)`.

        EXAMPLES::

            sage: R = Zp(3,20,'capped-rel', 'val-unit')
            sage: R(0).square_root()
            0
            sage: R(1).square_root()
            1 + O(3^20)
            sage: R(2).square_root(extend = False)
            Traceback (most recent call last):
            ...
            ValueError: element is not a square
            sage: R(4).square_root() == R(-2)
            True
            sage: R(9).square_root()
            3 * 1 + O(3^21)

        When p = 2, the precision of the square root is one less
        than the input::

            sage: R2 = Zp(2,20,'capped-rel')
            sage: R2(0).square_root()
            0
            sage: R2(1).square_root()
            1 + O(2^19)
            sage: R2(4).square_root()
            2 + O(2^20)

            sage: R2(9).square_root() == R2(3, 19) or R2(9).square_root() == R2(-3, 19)
            True

            sage: R2(17).square_root()
            1 + 2^3 + 2^5 + 2^6 + 2^7 + 2^9 + 2^10 + 2^13 + 2^16 + 2^17 + O(2^19)

            sage: R3 = Zp(5,20,'capped-rel')
            sage: R3(0).square_root()
            0
            sage: R3(1).square_root()
            1 + O(5^20)
            sage: R3(-1).square_root() == R3.teichmuller(2) or R3(-1).square_root() == R3.teichmuller(3)
            True


        TESTS::

            sage: R = Qp(3,20,'capped-rel')
            sage: R(0).square_root()
            0
            sage: R(1).square_root()
            1 + O(3^20)
            sage: R(4).square_root() == R(-2)
            True
            sage: R(9).square_root()
            3 + O(3^21)
            sage: R(1/9).square_root()
            3^-1 + O(3^19)

            sage: R2 = Qp(2,20,'capped-rel')
            sage: R2(0).square_root()
            0
            sage: R2(1).square_root()
            1 + O(2^19)
            sage: R2(4).square_root()
            2 + O(2^20)
            sage: R2(9).square_root() == R2(3,19) or R2(9).square_root() == R2(-3,19)
            True
            sage: R2(17).square_root()
            1 + 2^3 + 2^5 + 2^6 + 2^7 + 2^9 + 2^10 + 2^13 + 2^16 + 2^17 + O(2^19)

            sage: R3 = Qp(5,20,'capped-rel')
            sage: R3(0).square_root()
            0
            sage: R3(1).square_root()
            1 + O(5^20)
            sage: R3(-1).square_root() == R3.teichmuller(2) or R3(-1).square_root() == R3.teichmuller(3)
            True
        """
        # need special case for zero since pari(self) is the *integer* zero
        # whose square root is a real number....!
        if self.valuation() is infinity:
            return self
        try:
            # use pari
            ans = self.parent()(pari(self).sqrt())
            if all:
                return [ans, -ans]
            else:
                return ans
        except PariError:
            # todo: should eventually change to return an element of
            # an extension field
            if extend:
                raise NotImplementedError, "extending using the sqrt function not yet implemented"
            elif all:
                return []
            else:
                raise ValueError, "element is not a square"

    #def _unit_part(self):
    #    raise NotImplementedError

    cpdef abs(self, prec=None):
        """
        Return the `p`-adic absolute value of ``self``.

        This is normalized so that the absolute value of `p` is `1/p`.

        INPUT:

        - ``prec`` -- Integer.  The precision of the real field in which
          the answer is returned.  If ``None``, returns a rational for
          absolutely unramified fields, or a real with 53 bits of
          precision for ramified fields.

        EXAMPLES::

            sage: a = Qp(5)(15); a.abs()
            1/5
            sage: a.abs(53)
            0.200000000000000
            sage: Qp(7)(0).abs()
            0
            sage: Qp(7)(0).abs(prec=20)
            0.00000

        An unramified extension::

            sage: R = Zp(5,5)
            sage: P.<x> = PolynomialRing(R)
            sage: Z25.<u> = R.ext(x^2 - 3)
            sage: u.abs()
            1
            sage: (u^24-1).abs()
            1/5

        A ramified extension::

            sage: W.<w> = R.ext(x^5 + 75*x^3 - 15*x^2 + 125*x - 5)
            sage: w.abs()
            0.724779663677696
            sage: W(0).abs()
            0.000000000000000
        """
        K = self.parent()
        if not prec and K.e() > 1:
            prec = 53
        if prec:
            from sage.rings.real_mpfr import RealField
            if self.is_zero():
                return RealField(prec).zero()
            return RealField(prec)(K.prime())**(-self.ordp())
        else:
            if self.is_zero():
                return Rational(0)
            return Rational(K.prime())**(-self.valuation())

