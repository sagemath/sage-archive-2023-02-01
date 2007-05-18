"""
Elements of p-Adic Rings and Fields

AUTHOR:
    -- David Roe
    -- Genya Zaytman: documentation
    -- David Harvey: doctests
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

from __future__ import with_statement
cimport sage.rings.padics.local_generic_element
from sage.rings.padics.local_generic_element cimport LocalGenericElement
cimport sage.structure.element
from sage.structure.element cimport Element
cimport pow_computer
from sage.rings.integer cimport Integer
import sage.rings.rational_field

Rational = sage.rings.rational.Rational
infinity = sage.rings.infinity.infinity
PariError = sage.libs.pari.gen.PariError
pari = sage.libs.pari.gen.pari
QQ = sage.rings.rational_field.QQ

cdef class pAdicGenericElement(LocalGenericElement):
    def __init__(self, parent):
        self.prime_pow = <PowComputer_class> parent.prime_pow
        LocalGenericElement.__init__(self, parent)

    def __richcmp__(left, right, int op):
        return (<Element>left)._richcmp(right, op)

    cdef int _cmp_c_impl(left, Element right) except -2:
        """
        First compare valuations, then compare normalized
        residue of unit part.

        EXAMPLES:
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
            if x_ordp == infinity:
                return 0 # since both are zero
            else:
                p = left.prime_pow.prime
                a = left.unit_part().lift()
                b = right.unit_part().lift()
                prec = min(left.precision_relative(), right.precision_relative())
                ppow = p**prec
                a %= ppow
                b %= ppow
                if a < b:
                    return -1
                elif a == b:
                    return 0
                else:
                    return 1

    def __floordiv__(self, right):
        if self.parent().is_field():
            return self / right
        else:
            right = self.parent()(right)
            v, u = right._val_unit()
            return self.parent()(self / u).__rshift__(v)

    def __getitem__(self, n):
        r"""
        Returns the coefficient of $p^n$ in the series expansion of self, as an integer in the range $0$ to $p-1$.

        EXAMPLE:
            sage: R = Zp(7,4,'capped-rel','series'); a = R(1/3); a
            5 + 4*7 + 4*7^2 + 4*7^3 + O(7^4)
            sage: a[0]
            5
            sage: a[1]
            4
        """
        if isinstance(n, slice):
            if n.start == 0:
                raise ValueError, "due to limitations in Python 2.5, you must call the slice() function rather than using the [:] syntax in this case"
            if n.stop == sys.maxint:
                return self.slice(n.start, None, n.step)
            return self.slice(n.start, n.stop, n.step)
        if n < self.valuation():
            return self.parent()(0)
        if self.parent().is_field():
            return self.list()[n - self.valuation()]
        return self.list()[n]

    def __invert__(self):
        r"""
        Returns the multiplicative inverse of self.

        EXAMPLE:
            sage: R = Zp(7,4,'capped-rel','series'); a = R(3); a
            3 + O(7^4)
            sage: ~a
            5 + 4*7 + 4*7^2 + 4*7^3 + O(7^4)

        NOTES:
        The element returned is an element of the fraction field.
        """
        return self.parent().fraction_field()(self, relprec = self.relprec).__invert__()

    def __mod__(self, right):
        if right == 0:
            raise ZeroDivisionError
        if self.parent().is_field():
            return self.parent()(0)
        else:
            return self - (self // right) * right

    def _integer_(self):
        return Integer(self.lift())

    def _is_exact_zero(self):
        return False

    def str(self, mode=None):
        return self._repr(mode=mode)

    def _repr(self, mode = None, do_latex = False):
        r"""
        Prints a string representation of the element.  See __init__ for more details on print modes.

        EXAMPLES:
            sage: R = Zp(7,4,'capped-rel','val-unit'); a = R(364); a
            7 * 52 + O(7^5)
            sage: print a.str('terse')
            364 + O(7^5)
            sage: print a.str('series')
            3*7 + 7^3 + O(7^5)
            sage: K = Qp(7,4,'capped-rel','val-unit'); a = K(364); a
            7 * 52 + O(7^5)
            sage: print a.str('series')
            3*7 + 7^3 + O(7^5)
        """
        if self._is_exact_zero():
            return "0"
        if mode is None:
            mode = self.parent().print_mode()
        elif not ((mode == 'val-unit') or (mode == 'series') or (mode == 'terse')):
            raise TypeError, "printing mode must be one of 'val-unit', 'series' or 'terse'"
        pprint = self.parent().variable_name()
        if self.lift() == 0:
            if mode == 'val-unit' or mode == 'series':
                if do_latex:
                    return "O(%s^{%s})"%(pprint, self.precision_absolute())
                else:
                    return "O(%s^%s)"%(pprint, self.precision_absolute())
            elif mode == 'terse':
                if do_latex:
                    return "0 + O(%s^{%s})"%(pprint, self.precision_absolute())
                else:
                    return "0 + O(%s^%s)"%(pprint, self.precision_absolute())
        if mode == 'val-unit':
            if do_latex:
                if self.valuation() == 0:
                    return "%s + O(%s^{%s})"%(self.lift(), pprint, self.precision_absolute())
                if self.valuation() == 1:
                    return "%s \\cdot %s + O(%s^{%s})"%(pprint, self.unit_part().lift(), pprint, self.precision_absolute())
                return "%s^{%s} \\cdot %s + O(%s^{%s})"%(pprint, self.valuation(), self.unit_part().lift(), pprint, self.precision_absolute())
            else:
                if self.valuation() == 0:
                    return "%s + O(%s^%s)"%(self.lift(), pprint, self.precision_absolute())
                if self.valuation() == 1:
                    return "%s * %s + O(%s^%s)"%(pprint, self.unit_part().lift(), pprint, self.precision_absolute())
                return "%s^%s * %s + O(%s^%s)"%(pprint, self.valuation(), self.unit_part().lift(), pprint, self.precision_absolute())
        elif mode == 'terse':
            if self.valuation() < 0:
                ppow1 = str(self.prime() ** (-self.valuation()))
                if do_latex:
                    ppow2 = "%s^{%s}"%(pprint, -self.valuation())
                    if len(ppow1) < len(ppow2):
                        ppow = ppow1
                    else:
                        ppow = ppow2
                    return "%s/%s + O(%s^{%s})"%(self.unit_part().lift(), ppow, pprint, self.precision_absolute())
                else:
                    ppow2 = "%s^%s"%(p, -self.valuation())
                    if len(ppow1) < len(ppow2):
                        ppow = ppow1
                    else:
                        ppow = ppow2
                    return "%s/%s + O(%s^{%s})"%(self.unit_part().lift(), ppow, pprint, self.precision_absolute())
            else:
                if do_latex:
                    return "%s + O(%s^{%s})"%(self.lift(), pprint, self.precision_absolute())
                else:
                    return "%s + O(%s^%s)"%(self.lift(), pprint, self.precision_absolute())
        else:
            slist = self.list()
            s = ""
            p = self.prime_pow.prime
            if self.parent().is_field():
                exp = self.valuation()
            else:
                exp = 0
            for a in slist:
                if a != 0:
                    if exp == 0:
                        s += "%s + "%a
                    else:
                        var = pprint
                        if exp != 1:
                            if do_latex:
                                var += "^{%s}"%exp
                            else:
                                var += "^%s"%exp
                        if a != 1:
                            if do_latex:
                                s += "%s \\cdot %s + "%(a, var)
                            else:
                                s += "%s*%s + "%(a, var)
                        else:
                            s += "%s + "%var
                exp += 1
            s += "O(%s"%(pprint)
            if self.precision_absolute() == 1:
                s += ")"
            else:
                if do_latex:
                    s += "^{%s})"%self.precision_absolute()
                else:
                    s += "^%s)"%self.precision_absolute()
            return s

    def additive_order(self, prec):
        r"""
        Returns the additive order of self, where self is considered to be zero if it is zero modulo $p^{\mbox{prec}}$.

        INPUT:
            self -- a p-adic element
            prec -- an integer
        OUTPUT:
            integer -- the additive order of self

        """
        if self.is_zero(prec):
            return Integer(1)
        else:
            return infinity

    def algdep(self, n):
        """
        Returns a polynomial of degree at most $n$ which is approximately
        satisfied by this number.  Note that the returned polynomial
        need not be irreducible, and indeed usually won't be if this number
        is a good approximation to an algebraic number of degree less than $n$.

        ALGORITHM: Uses the PARI C-library algdep command.

        INPUT:
            self -- a p-adic element
            n -- an integer
        OUTPUT:
            polynomial -- degree n polynomial aproximately satisfied by self

        EXAMPLES:
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
        return sage.rings.arith.algdep(self, n)

    def algebraic_dependency(self, n):
        return self.algdep(n)

    def exp(self):
        r"""
        Compute the p-adic exp of any element of $\Z_p$ where the
        series converges.

        EXAMPLES:
        Borrowed from log.

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

        The next few examples illustrate precision when computing $p$-adic exps.
        First we create a field with {\em default} precision 10.
            sage: R = Zp(5,10, print_mode='series')
            sage: e = R(2*5 + 2*5**2 + 4*5**3 + 3*5**4 + 5**5 + 3*5**7 + 2*5**8 + 4*5**9).add_bigoh(10); e
                2*5 + 2*5^2 + 4*5^3 + 3*5^4 + 5^5 + 3*5^7 + 2*5^8 + 4*5^9 + O(5^10)
            sage: e.exp()*R.teichmuller(4)
                4 + 2*5 + 3*5^3 + O(5^10)

            sage: K = Qp(5,10, print_mode='series')
            sage: e = K(2*5 + 2*5**2 + 4*5**3 + 3*5**4 + 5**5 + 3*5**7 + 2*5**8 + 4*5**9).add_bigoh(10); e
                2*5 + 2*5^2 + 4*5^3 + 3*5^4 + 5^5 + 3*5^7 + 2*5^8 + 4*5^9 + O(5^10)
            sage: e.exp()*K.teichmuller(4)
                4 + 2*5 + 3*5^3 + O(5^10)

        TESTS:
        Check that results are consistent over a range of precision:

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
            -- Genya Zaytman (2007-02-15)

        """

        val = self.valuation()
        p = self.parent().prime()
        if val > 1 or (val > 0 and p != 2):

            prec = self.precision_absolute()

            # I believe this works
            max_term = ((p-1)*(prec-1))//((p-1)*val - 1) + 1

            # Need extra precision to take into account powers of p
            # in the denominators of the series. (Indeed, it's a
            # not-entirely-trivial fact that if x is given mod p^n, that
            # exp(x) is well-defined mod p^n !) .
            extra_prec = max_term//(p-1)

            from sage.rings.padics.factory import Zp
            working_ring = Zp(p, prec + extra_prec, type = 'capped-abs')
            x = working_ring(self.lift())
            term = ans = working_ring(Integer(1))
            for n in range(1, max_term):
                term *=x
                term = term // working_ring(Integer(n))
                ans += term
            # Note that it is the absolute precision that is respected by exp: even when p == 2?
            return self.parent()(ans).add_bigoh(prec)
        else:
            raise ValueError, "series doesn't converge"

    def exp_artin_hasse(self):
        raise NotImplementedError

    def gamma(self):
        raise NotImplementedError

    def is_square(self): #should be overridden for lazy elements
        """
        Returns whether self is a square

        INPUT:
            self -- a p-adic element
        OUTPUT:
            boolean -- whether self is a square

        EXAMPLES:
            sage: R = Zp(3,20,'capped-rel')
            sage: R(0).is_square()
            True
            sage: R(1).is_square()
            True
            sage: R(2).is_square()
            False
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
        if self.valuation() == infinity:
            return True
        elif self.precision_relative() < 1:
            return True
        elif self.parent().prime() != 2:
            return bool((self.valuation() % 2 == 0) and (self.unit_part().residue(1).is_square()))
        else:
            return bool((self.valuation() % 2 == 0) and (self.unit_part().residue(3) == 1))

    def log(self, branch = None):
        r"""
        Compute the p-adic logarithm of any unit in $\Z_p$.
        (See below for normalization.)

        The usual power series for log with values in the additive
        group of $\Z_p$ only converges for 1-units (units congruent to
        1 modulo p).  However, there is a unique extension of log to a
        homomorphism defined on all the units.  If u = a*v is a unit
        with v = 1 (mod p), then we define log(u) = log(v).  This is
        the correct extension because the units U of Z_p splits as a
        product U = V x <w>, where V is the subgroup of 1-units and w
        is a (p-1)st root of unity.  The <w> factor is torsion, so
        must go to 0 under any homomorphism to the torsion free group
        $(\Z_p, +)$.

        Notes -- What some other systems do:
           PARI:  Seems to define log the same way as we do.
           MAGMA: Gives an error when unit is not a 1-unit.

        Algorithm:
           Input: Some p-adic unit u.
           1. Check that the input p-adic number is really a unit
              (i.e., valuation 0)
           2. Let $1-x = u^{p-1}$, which is a 1-unit.
           3. Use the series expansion
              $$
                \log(1-x) = F(x) = -x - 1/2*x^2 - 1/3*x^3 - 1/4*x^4 - 1/5*x^5 - ...
              $$
              to compute the logarithm log(u**(p-1)).  Use enough
              terms so that terms added on are zero
           4. Then $$\log(u) = log(u^{p-1})/(p-1) = F(1-u^{p-1})/(p-1).$$

        EXAMPLES:
            sage: Z13 = Zp(13, 10, print_mode='series')
            sage: a = Z13(14); a
            1 + 13 + O(13^10)

        Note that the relative precision decreases when we take log: it is the absolute
        precision that is preserved.
            sage: a.log()
            13 + 6*13^2 + 2*13^3 + 5*13^4 + 10*13^6 + 13^7 + 11*13^8 + 8*13^9 + O(13^10)
            sage: Q13 = Qp(13, 10, print_mode='series')
            sage: a = Q13(14); a
            1 + 13 + O(13^10)
            sage: a.log()
            13 + 6*13^2 + 2*13^3 + 5*13^4 + 10*13^6 + 13^7 + 11*13^8 + 8*13^9 + O(13^10)

        The next few examples illustrate precision when computing $p$-adic logs.
        First we create a field with {\em default} precision 10.
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

        Check that results are consistent over a range of precision:
            sage: max_prec = 40
            sage: p = 3
            sage: K = Zp(p, max_prec)
            sage: full_log = (K(1 + p)).log()
            sage: for prec in range(2, max_prec):
            ...       ll = (K(1 + p).add_bigoh(prec)).log()
            ...       assert ll == full_log
            ...       assert ll.precision_absolute() == prec


        AUTHORS:
            -- William Stein: initial version
            -- David Harvey (2006-09-13): corrected subtle precision bug
               (need to take denominators into account! -- see trac \#53)
            -- Genya Zaytman (2007-02-14): adapted to new p-adic class

        TODO:
            -- Currently implemented as $O(N^2)$. This can be improved to
            soft-$O(N)$ using algorithm described by Dan Bernstein:
            http://cr.yp.to/lineartime/multapps-20041007.pdf

        """

        p = self.parent().prime()
        # Step 1 -- a unit?
        if self.is_unit() and self.unit_part().residue(1) == 1:
            # It's already a 1-unit, so just use the series
            # (base case of "induction")

            prec = self.precision_absolute()

            # Need extra precision to take into account powers of p
            # in the denominators of the series. (Indeed, it's a
            # not-entirely-trivial fact that if x is given mod p^n, that
            # log(x) is well-defined mod p^n !) Specifically:
            # we are only guaranteed that $x^j/j$ is zero mod $p^n$ if
            # j >= floor(log_p(j)) + n.
            extra_prec = 0
            while extra_prec < Integer(prec + extra_prec).exact_log(p):
                extra_prec += 1

            x = Integer(1) - self
            from sage.rings.padics.factory import Zp
            working_ring = Zp(p, prec + extra_prec, type = 'capped-abs')
            x = working_ring(x.lift())
            xpow = x
            ans = working_ring(Integer(0))
            for n in range(1, prec + extra_prec):
                ans -= xpow//working_ring(Integer(n))
                xpow *= x
            # Note that it is the absolute precision that is respected by log
            return self.parent()(ans.lift()).add_bigoh(prec)
        elif self.is_unit():
            z = self.unit_part()
            return (z**Integer(p-1)).log() // Integer(p-1)
        elif not branch is None and self.parent().__contains__(branch):
            branch = self.parent()(branch)
            return self.unit_part().log() + branch*self.valuation()
        else:
            raise ValueError, "not a unit"

    def log_artin_hasse(self):
        raise NotImplementedError

    def minimal_polynomial(self, name):
        """
        Returns a minimal polynomial of this p-adic element, i.e., x - self

        INPUT:
            self -- a p-adic element
            name -- string the name of the variable

        OUTPUT:
            polynomial -- a minimal polynomial of this p-adic element, i.e., x - self
        """
        import sage.rings.polynomial.polynomial_ring
        R = sage.rings.polynomial.polynomial_ring.PolynomialRing(self.parent(), name)
        return R.gen() - R(self)

    def multiplicative_order(self, prec = None): #needs to be rewritten for lazy elements
        r"""
        Returns the multiplicative order of self, where self is considered to be one if it is one modulo $p^{\mbox{prec}}$.

        INPUT:
            self -- a p-adic element
            prec -- an integer
        OUTPUT:
            integer -- the multiplicative order of self
        EXAMPLES:
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
        if prec is None:
            if self == self.parent().teichmuller(res):
                return res.multiplicative_order()
            else:
                return infinity
        if self.is_equal_to(self.parent().teichmuller(res),prec): #should this be made more efficient?
            return res.multiplicative_order()
        else:
            return infinity

    def norm(self, ground=None):
        """
        Returns the norm of this p-adic element over the ground ring

        INPUT:
            self -- a p-adic element
            ground -- a subring of the ground ring (default: base ring)

        OUTPUT:
            element -- the norm of this p-adic element over the ground ring
        """
        if (ground != None) and (ground != self.parent()):
            raise ValueError, "Ground Field not a subfield"
        else:
            return self

    def ordp(self):
        r"""
        Returns the valuation of self, normalized so that the valuation of p is 1

        INPUT:
            self -- a p-adic element
        OUTPUT:
            integer -- the valuation of self, normalized so that the valuation of p is 1
        EXAMPLES:
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
        return self.valuation()

    def rational_reconstruction(self):
        r"""
        Returns a rational approximation to this p-adic number

        INPUT:
            self -- a p-adic element
        OUTPUT:
            rational -- an approximation to self
        EXAMPLES:
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

    def square_root(self, extend = True, all = False):
        r"""
        Returns the square root of this p-adic number

        INPUT:
            self -- a p-adic element
            extend -- bool (default: True); if True, return a square root
                in an extension if necessary; if False and no root exists
                in the given ring or field, raise a ValueError
            all -- bool (default: False); if True, return a list of all
                square roots
        OUTPUT:
            p-adic element -- the square root of this p-adic number

            If all = False, the square root chosen is the one whose reduction mod p is in
            the range [0, p/2).

        EXAMPLES:
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

        When p = 2, the precision of the square root is one less than the
        input:
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


        Tests for fields:
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
            else:
                raise ValueError, "element is not a square"

    def trace(self, ground=None):
        """
        Returns the trace of this p-adic element over the ground ring

        INPUT:
            self -- a p-adic element
            ground -- a subring of the ground ring (default: base ring)

        OUTPUT:
            element -- the trace of this p-adic element over the ground ring
        """
        if (ground != None) and (ground != self.parent()):
            raise ValueError, "Ground ring not a subring"
        else:
            return self

    #def _unit_part(self):
    #    raise NotImplementedError

    def _val_unit(self):
        return self.valuation(), self.unit_part().lift()

    cdef base_p_list(self, mpz_t value, lift_mode):
        cdef mpz_t tmp, halfp
        cdef int neg, curpower
        cdef Integer list_elt
        cdef unsigned long preccap
        preccap = self.prime_pow._cache_limit
        ans = PyList_New(0)
        mpz_init_set(tmp, value)
        if lift_mode == 'simple':
            while mpz_sgn(tmp) != 0:
                list_elt = PY_NEW(Integer)
                mpz_mod(list_elt.value, tmp, self.prime_pow.prime.value)
                mpz_sub(tmp, tmp, list_elt.value)
                mpz_divexact(tmp, tmp, self.prime_pow.prime.value)
                PyList_Append(ans, list_elt)
        elif lift_mode == 'smallest':
            neg = 0
            curpower = preccap
            mpz_init(halfp)
            mpz_fdiv_q_2exp(halfp, self.prime_pow.prime.value, 1)
            while mpz_sgn(tmp) != 0:
                curpower -= 1
                list_elt = PY_NEW(Integer)
                mpz_mod(list_elt.value, tmp, self.prime_pow.prime.value)
                if mpz_cmp(list_elt.value, halfp) >= 0:
                    mpz_sub(list_elt.value, list_elt.value, self.prime_pow.prime.value)
                    neg = 1
                else:
                    neg = 0
                mpz_sub(tmp, tmp, list_elt.value)
                mpz_divexact(tmp, tmp, self.prime_pow.prime.value)
                if neg == 1:
                    if mpz_cmp(tmp, self.prime_pow.dense_list[curpower]) >= 0:
                        mpz_sub(tmp, tmp, self.prime_pow.dense_list[curpower])
                PyList_Append(ans, list_elt)
            mpz_clear(halfp)
        else:
            mpz_clear(tmp)
            raise ValueError, "lift mode must be one of 'simple', 'smallest' or 'teichmuller'"
        mpz_clear(tmp)
        return ans

cdef public void teichmuller_set_c(mpz_t value, mpz_t p, mpz_t ppow):
    r"""
    Sets value to the integer between 0 and p^prec that is congruent to the Teichmuller lift of value to Z_p.

    INPUT:

    value -- An mpz_t currently holding an approximation to the
             Teichmuller representative (this approximation can
             be any integer).  It will be set to the actual
             Teichmuller lift
    p -- An mpz_t holding the value of the prime
    ppow -- An mpz_t holding the value p^prec, where prec is the desired precision of the Teichmuller lift
    """
    cdef mpz_t u, xnew
    if mpz_divisible_p(value, p) != 0:
        mpz_set_ui(value, 0)
        return
    if mpz_sgn(value) < 0 or mpz_cmp(value, ppow) >= 0:
        mpz_mod(value, value, ppow)
    mpz_init(u)
    mpz_init(xnew)
    # u = 1 / Mod(1 - p, ppow)
    mpz_sub(u, ppow, p)
    mpz_add_ui(u, u, 1)
    mpz_invert(u, u, ppow)
    # Consider x as Mod(self.value, ppow)
    # xnew = x + u*(x^p - x)
    mpz_powm(xnew, value, p, ppow)
    mpz_sub(xnew, xnew, value)
    mpz_mul(xnew, xnew, u)
    mpz_add(xnew, xnew, value)
    mpz_mod(xnew, xnew, ppow)
    # while x != xnew:
    #     x = xnew
    #     xnew = x + u*(x^p - x)
    while mpz_cmp(value, xnew) != 0:
        mpz_set(value, xnew)
        mpz_powm(xnew, value, p, ppow)
        mpz_sub(xnew, xnew, value)
        mpz_mul(xnew, xnew, u)
        mpz_add(xnew, xnew, value)
        mpz_mod(xnew, xnew, ppow)
    mpz_clear(u)
    mpz_clear(xnew)

