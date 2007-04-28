"""
Power Series

SAGE provides an implementation of dense and sparse power series over
any SAGE base ring.

AUTHOR:
   -- William Stein
   -- David Harvey (2006-09-11): added solve_linear_de() method
   -- Robert Bradshaw (2007-04): sqrt, rmul, lmul, shifting
   -- Robert Bradshaw (2007-04): SageX version

EXAMPLE:
    sage: R.<x> = PowerSeriesRing(ZZ)
    sage: R([1,2,3])
    1 + 2*x + 3*x^2
    sage: R([1,2,3], 10)
    1 + 2*x + 3*x^2 + O(x^10)
    sage: f = 1 + 2*x - 3*x^3 + O(x^4); f
    1 + 2*x - 3*x^3 + O(x^4)
    sage: f^10
    1 + 20*x + 180*x^2 + 930*x^3 + O(x^4)
    sage: g = 1/f; g
    1 - 2*x + 4*x^2 - 5*x^3 + O(x^4)
    sage: g * f
    1 + O(x^4)

In Python (as opposted to SAGE) create the power series ring and its
generator as follows:

    sage: R, x = objgen(PowerSeriesRing(ZZ, 'x'))
    sage: parent(x)
    Power Series Ring in x over Integer Ring

EXAMPLE: COERCION
This example illustrates that coercion for power series rings
is consistent with coercion for polynomial rings.

    sage: poly_ring1.<gen1> = PolynomialRing(QQ)
    sage: poly_ring2.<gen2> = PolynomialRing(QQ)
    sage: huge_ring.<x> = PolynomialRing(poly_ring1)

The generator of the first ring gets coerced in as itself,
since it is the base ring.
    sage: huge_ring(gen1)
    gen1

The generator of the second ring gets mapped via the
natural map sending one generator to the other.
    sage: huge_ring(gen2)
    x

With power series the behavior is the same.
    sage: power_ring1.<gen1> = PowerSeriesRing(QQ)
    sage: power_ring2.<gen2> = PowerSeriesRing(QQ)
    sage: huge_power_ring.<x> = PowerSeriesRing(power_ring1)
    sage: huge_power_ring(gen1)
    gen1
    sage: huge_power_ring(gen2)
    x
"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include "../ext/stdsage.pxi"

import operator

from infinity import infinity, is_Infinite
from polynomial_ring import PolynomialRing
import polynomial_element_generic
import polynomial_element
import power_series_ring
import sage.misc.misc
import ring_element
import arith
import sage.misc.latex
import sage.structure.coerce
import rational_field, integer_ring
from integer import Integer
from integer_mod_ring import IntegerModRing
import sage.libs.pari.all
from sage.libs.all import PariError
from sage.misc.functional import sqrt, log
from sage.rings.arith import ceil

from sage.rings.ring import is_Field

Polynomial = polynomial_element.Polynomial_generic_dense

from sage.structure.element cimport AlgebraElement, RingElement, ModuleElement, Element

def is_PowerSeries(x):
    return isinstance(x, PowerSeries)

cdef class PowerSeries(AlgebraElement):
    """
    A power series.
    """

    def __init__(self, parent, prec, is_gen=False):
        """
        Initialize a power series.
        """
        AlgebraElement.__init__(self, parent)
        self.__is_gen = is_gen
        if not (prec is infinity):
            prec = int(prec)
        self._prec = prec

    def __reduce__(self):
        return make_element_from_parent, (self._parent, self.polynomial(), self.prec())

    def is_sparse(self):
        """
        EXAMPLES:
            sage: R.<t> = PowerSeriesRing(ZZ)
            sage: t.is_sparse()
            False
            sage: R.<t> = PowerSeriesRing(ZZ, sparse=True)
            sage: t.is_sparse()
            True
        """
        return self._parent.is_sparse()

    def is_dense(self):
        """
        EXAMPLES:
            sage: R.<t> = PowerSeriesRing(ZZ)
            sage: t.is_dense()
            True
            sage: R.<t> = PowerSeriesRing(ZZ, sparse=True)
            sage: t.is_dense()
            False
        """
        return self._parent.is_dense()

    def is_gen(self):
        """
        Returns True if this the generator (the variable) of the power series ring.

        EXAMPLES:
            sage: R.<t> = QQ[[]]
            sage: t.is_gen()
            True
            sage: (1 + 2*t).is_gen()
            False

        Note that this only returns true on the actual generator, not on
        something that happens to be equal to it.
            sage: (1*t).is_gen()
            False
            sage: 1*t == t
            True
        """
        return bool(self.__is_gen)

    def _im_gens_(self, codomain, im_gens):
        """
        Returns the image of this series under the map that sends the generators
        to im_gens.   This is used internally for computing homomorphisms.

        EXAMPLES:
            sage: R.<t> = QQ[[]]
            sage: f = 1 + t + t^2
            sage: f._im_gens_(ZZ, [3])
            13
        """
        return codomain(self(im_gens[0]))

    def base_extend(self, R):
        """
        Return a copy of this power series but with coefficients in R.

        The following coercion uses base_extend implicitly:
            sage: R.<t> = ZZ[['t']]
            sage: (t - t^2) * Mod(1, 3)
            t + 2*t^2
        """
        S = self._parent.base_extend(R)
        return S(self)

    def change_ring(self, R):
        """
        Change if possible the coefficients of self to lie in R.

        EXAMPLES:
            sage: R.<T> = QQ[[]]; R
            Power Series Ring in T over Rational Field
            sage: f = 1 - 1/2*T + 1/3*T^2 + O(T^3)
            sage: f.base_extend(GF(5))
            Traceback (most recent call last):
            ...
            TypeError: no base extension defined
            sage: f.change_ring(GF(5))
            1 + 2*T + 2*T^2 + O(T^3)
            sage: f.change_ring(GF(3))
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Inverse does not exist.

        We can only change irng if there is a __call__ coercion defined.
        The following succeeds because ZZ(K(4)) is defined.

            sage: K.<a> = NumberField(cyclotomic_polynomial(3), 'a')
            sage: R.<t> = K[['t']]
            sage: (4*t).change_ring(ZZ)
            4*t

        This does not succeed because ZZ(K(a+1)) is not defined.
            sage: K.<a> = NumberField(cyclotomic_polynomial(3), 'a')
            sage: R.<t> = K[['t']]
            sage: ((a+1)*t).change_ring(ZZ)
            Traceback (most recent call last):
            ...
            TypeError: Unable to coerce a + 1 to an integer

        """
        S = self._parent.change_ring(R)
        return S(self)

    def __cmp__(left, right):
        return (<Element>left)._cmp(right)

    cdef int _cmp_c_impl(self, Element right) except -2:
        r"""
        Comparison of self and right.

        We say two approximate power series are equal, if they agree
        for all coefficients up to the *minimum* of the precisions of
        each.  Thus, e.g., $f=1+q+O(q^2)$ is equal to $g=1+O(q)$.
        This is how PARI defines equality of power series, but not how
        MAGMA defines equality.  (MAGMA would declare f and g
        unequal.)  I side with PARI, because even if $g=1+q+O(q^2)$,
        we don't really know whether f equals g, since we don't know
        the coefficients of $q^2$.

        Comparison is done in dictionary order from lowest degree to
        highest degree coefficients (this is different than
        polynomials).

        EXAMPLES:
            sage: R.<q> = ZZ[[ ]]; R
            Power Series Ring in q over Integer Ring
            sage: f=1+q+O(q^2); g = 1+O(q)
            sage: f == g
            True
            sage: 1 - 2*q + q^2 +O(q^3) == 1 - 2*q^2 + q^2 + O(q^4)
            False
        """
        # A very common case throughout code
        # TODO: change code to use is_zero() or equivalent
        if PY_TYPE_CHECK(right, int):
            return self.is_zero()

        prec = self.common_prec(right)
        x = self.list()
        y = right.list()
        if not (prec is infinity):
            x = x[:prec]
            y = y[:prec]
        return cmp(x,y)

    def __call__(self, x):   # you *MUST* overrride this in the derived class
        raise NotImplementedError

    def list(self):          # you *MUST* overrride this in the derived class
        raise NotImplementedError

    def polynomial(self):          # you *MUST* overrride this in the derived class
        raise NotImplementedError

    def __setitem__(self, n, value):   # you *MUST* overrride this in the derived class
        raise NotImplementedError

    def __copy__(self):
        """
        Return this power series.  Power series are immutable so copy
        can safely just return the same polynomial.

        EXAMPLES:
            sage: R.<q> = ZZ[[ ]]; R
            Power Series Ring in q over Integer Ring
            sage: f = 1 + 3*q + O(q^10)
            sage: copy(f) is f       # !!! ok since power series are immutable.
            True
        """
        return self

    def base_ring(self):
        """
        Return the base ring that this power series is defined over.

        EXAMPLES:
            sage: R.<t> = GF(49,'alpha')[[]]
            sage: (t^2 + O(t^3)).base_ring()
            Finite Field in alpha of size 7^2
        """
        return self._parent.base_ring()

    def padded_list(self, n):
        """
        Return list of coefficients of self up to (but not include $q^n$).

        Includes 0's in the list on the right so that the list has
        length $n$.

        EXAMPLES:
            sage: R.<q> = PowerSeriesRing(QQ)
            sage: f = 1 - 17*q + 13*q^2 + 10*q^4 + O(q^7)
            sage: f.list()
            [1, -17, 13, 0, 10]
            sage: f.padded_list(7)
            [1, -17, 13, 0, 10, 0, 0]
            sage: f.padded_list(10)
            [1, -17, 13, 0, 10, 0, 0, 0, 0, 0]
            sage: f.padded_list(3)
            [1, -17, 13]
        """
        v = self.list()
        if len(v) < n:
            z = self._parent.base_ring()(0)
            return v + [z]*(n - len(v))
        else:
            return v[:int(n)]

    def prec(self):
        """
        The precision of $...+O(x^r)$ is by definition $r$.

        EXAMPLES:
            sage: R.<t> = ZZ[[]]
            sage: (t^2 + O(t^3)).prec()
            3
            sage: (1 - t^2 + O(t^100)).prec()
            100
        """
        return self._prec

    def _repr_(self):
        """
        Return string represenation of this power series.

        EXAMPLES:
            sage: R.<t> = ZZ[[]]
            sage: (t^2 + O(t^3))._repr_()
            't^2 + O(t^3)'

            sage: R.<t> = QQ[[]]
            sage: 1 / (1+2*t +O(t^5))
            1 - 2*t + 4*t^2 - 8*t^3 + 16*t^4 + O(t^5)

            sage: R.<t> = PowerSeriesRing(QQ, sparse=True)
            sage: 1 / (1+2*t +O(t^5))
            1 - 2*t + 4*t^2 - 8*t^3 + 16*t^4 + O(t^5)
            sage: -13/2 * t^3  + 5*t^5 + O(t^10)
            -13/2*t^3 + 5*t^5 + O(t^10)

        """
        if self.is_zero():
            if self.prec() is infinity:
                return "0"
            else:
                return "O(%s^%s)"%(self._parent.variable_name(),self.prec())

        atomic_repr = self._parent.base_ring().is_atomic_repr()
        X = self._parent.variable_name()

        s = " "
        if self.is_sparse():
            f = self.polynomial()
            m = f.degree() + 1
            d = f._dict_unsafe()
            coeffs = list(d.iteritems())
            coeffs.sort()
            for (n, x) in coeffs:
                if x != 0:
                    if s != ' ':
                        s += " + "
                    x = str(x)
                    if not atomic_repr and n > 0 and (x.find("+") != -1 or x.find("-") != -1):
                        x = "(%s)"%x
                    if n > 1:
                        var = "*%s^%s"%(X,n)
                    elif n==1:
                        var = "*%s"%X
                    else:
                        var = ""
                    s += "%s%s"%(x,var)
        else:
            v = self.list()
            m = len(v)
            first = True
            for n in xrange(m):
                x = v[n]
                if x != 0:
                    if not first:
                        s += " + "
                    x = str(x)
                    if not atomic_repr and n > 0 and (x[1:].find("+") != -1 or x[1:].find("-") != -1):
                        x = "(%s)"%x
                    if n > 1:
                        var = "*%s^%s"%(X,n)
                    elif n==1:
                        var = "*%s"%X
                    else:
                        var = ""
                    s += "%s%s"%(x,var)
                    first = False
        # end

        if atomic_repr:
            s = s.replace(" + -", " - ")
        s = s.replace(" 1*"," ")
        s = s.replace(" -1*", " -")
        if not (self._prec is infinity):
            if self._prec == 0:
                bigoh = "O(1)"
            elif self._prec == 1:
                bigoh = "O(%s)"%self._parent.variable_name()
            else:
                bigoh = "O(%s^%s)"%(self._parent.variable_name(),self._prec)
            if s==" ":
                return bigoh
            s += " + %s"%bigoh
        return s[1:]

    def _latex_(self):
        r"""
        Return latex representation of this power series.

        EXAMPLES:
            sage: R.<t> = QQ[[]]
            sage: f = -1/2 * t + 2/3*t^2 + -9/7 * t^15 + O(t^20); f
            -1/2*t + 2/3*t^2 - 9/7*t^15 + O(t^20)
            sage: latex(f)
            -\frac{1}{2}t + \frac{2}{3}t^{2} - \frac{9}{7}t^{15} + O(\text{t}^{20})
        """
        if self.is_zero():
            if self.prec() is infinity:
                return "0"
            else:
                return "0 + \\cdots"
        s = " "
        v = self.list()
        m = len(v)
        X = self._parent.variable_name()
        atomic_repr = self._parent.base_ring().is_atomic_repr()
        first = True
        for n in xrange(m):
            x = v[n]
            if x != 0:
                if not first:
                    s += " + "
                x = sage.misc.latex.latex(x)
                if not atomic_repr and n > 0 and (x[1:].find("+") != -1 or x[1:].find("-") != -1):
                    x = "\\left(%s\\right)"%x
                if n > 1:
                    var = "%s^{%s}"%(X,n)
                elif n==1:
                    var = "%s"%X
                else:
                    var = ""
                if n > 0:
                    s += "%s|%s"%(x,var)
                else:
                    s += str(x)
                first = False

        if atomic_repr:
            s = s.replace(" + -", " - ")
        s = s.replace(" -1|", " -")
        s = s.replace(" 1|"," ")
        s = s.replace("|","")
        if not (self._prec is infinity):
            if self._prec == 0:
                bigoh = "O(1)"
            elif self._prec == 1:
                bigoh = "O(%s)"%sage.misc.latex.latex(self._parent.variable_name())
            else:
                bigoh = "O(%s^{%s})"%(sage.misc.latex.latex(self._parent.variable_name()),self._prec)
            if s == " ":
                return bigoh
            s += " + %s"%bigoh
        return s[1:]


    def truncate(self, prec=infinity):
        """
        The polynomial obtained from power series by truncation.

        EXAMPLES:
            sage: R.<I> = GF(2)[[]]
            sage: f = 1/(1+I+O(I^8)); f
            1 + I + I^2 + I^3 + I^4 + I^5 + I^6 + I^7 + O(I^8)
            sage: f.truncate(5)
            I^4 + I^3 + I^2 + I + 1
        """
        if prec is infinity:
            prec = self._prec
        a = self.list()
        v = [a[i] for i in range(min(prec, len(a)))]
        return self._parent._poly_ring()(v)

    def add_bigoh(self, prec):
        r"""
        Returns the power series of precision at most prec got by
        adding $O(q^\text{prec})$ to f, where q is the variable.

        EXAMPLES:
            sage: R.<A> = RDF[[]]
            sage: f = (1+A+O(A^5))^5; f
            1.0 + 5.0*A + 10.0*A^2 + 10.0*A^3 + 5.0*A^4 + O(A^5)
            sage: f.add_bigoh(3)
            1.0 + 5.0*A + 10.0*A^2 + O(A^3)
        """
        if prec is infinity or prec >= self.prec():
            return self
        a = self.list()
        v = [a[i] for i in range(min(prec, len(a)))]
        return self._parent(v, prec)

    def __getitem__(self,n):
        r"""
        Return the coefficient of $t^n$ in this power series, where
        $t$ is the indeterminate of the power series ring.

        If n is negative returns 0.  If n is beyond the precision,
        raises an IndexError.

        EXAMPLES:
            sage: R.<m> = CDF[[]]
            sage: f = pi^2 + m^3 + e*m^4 + O(m^10); f
            9.86960440109 + 1.0*m^3 + 2.71828182846*m^4 + O(m^10)
            sage: f[-5]
            0
            sage: f[0]
            9.86960440109
            sage: f[4]
            2.71828182846
            sage: f[9]
            0
            sage: f[10]
            Traceback (most recent call last):
            ...
            IndexError: coefficient not known
            sage: f[1000]
            Traceback (most recent call last):
            ...
            IndexError: coefficient not known
        """
        if n<0:
            return self.base_ring()(0)
        c = self.list()
        if n >= len(c):
            if self._prec > n:
                return self.base_ring()(0)
            else:
                raise IndexError, "coefficient not known"
        return c[n]

    def common_prec(self, f):
        r"""
        Returns minimum precision of $f$ and self.

        EXAMPLES:
            sage: R.<t> = PowerSeriesRing(QQ)

            sage: f = t + t^2 + O(t^3)
            sage: g = t + t^3 + t^4 + O(t^4)
            sage: f.common_prec(g)
            3
            sage: g.common_prec(f)
            3

            sage: f = t + t^2 + O(t^3)
            sage: g = t^2
            sage: f.common_prec(g)
            3
            sage: g.common_prec(f)
            3

            sage: f = t + t^2
            sage: f = t^2
            sage: f.common_prec(g)
            +Infinity
        """
        if self.prec() is infinity:
            return f.prec()
        elif f.prec() is infinity:
            return self.prec()
        return min(self.prec(), f.prec())

    cdef common_prec_c(self, PowerSeries f):
        if self._prec is infinity:
            return f._prec
        elif f._prec is infinity:
            return self._prec
        elif self._prec < f._prec:
            return self._prec
        else:
            return f._prec

    cdef RingElement _mul_c_impl(self, RingElement right_r):
        cdef PowerSeries right = <PowerSeries>right_r
        if self.is_zero():
            return self
        if right.is_zero():
            return right
        sp = self._prec
        rp = right_.prec
        if sp is infinity:
            if rp is infinity:
                prec = infinity
            else:
                prec = rp + self.valuation()
        else:  # sp != infinity
            if rp is infinity:
                prec = sp + right.valuation()
            else:
                prec = min(rp + self.valuation(), sp + right.valuation())
        # endif
        return self._mul_(right, prec) # ???

    def is_zero(self):
        """
        Return True if this power series equals 0.

        EXAMPLES:
            sage: R.<q> = ZZ[[ ]]; R
            Power Series Ring in q over Integer Ring
            sage: f = 1 + 3*q + O(q^10)
            sage: f.is_zero()
            False
            sage: (0 + O(q^2)).is_zero()
            True
            sage: R(0).is_zero()
            True
            sage: (0 + O(q^1000)).is_zero()
            True
        """
        return self.polynomial().is_zero()

    def is_unit(self):
        """
        Returns whether this power series is invertible, which is the
        case precisely when the constant term is invertible.

        EXAMPLES:
            sage: R.<t> = PowerSeriesRing(ZZ)
            sage: (-1 + t - t^5).is_unit()
            True
            sage: (3 + t - t^5).is_unit()
            False

        AUTHOR: David Harvey (2006-09-03)
        """
        return self[0].is_unit()

    def __invert__(self):
        """
        Inverse of the power series (i.e. a series Y such that XY = 1).
        The first nonzero coefficient must be a unit in the coefficient ring.
        If the valuation of the series is positive, this function will return
        a Laurent series.

        ALGORITHM:
            Uses Newton's method. Complexity is around $O(M(n) \log n)$,
            where $n$ is the precision and $M(n)$ is the time required to
            multiply polynomials of length $n$.

        EXAMPLES:
            sage: R.<q> = QQ[[]]
            sage: 1/(1+q + O(q**2))
            1 - q + O(q^2)
            sage: 1/(1+q)
            1 - q + q^2 - q^3 + q^4 - q^5 + q^6 - q^7 + q^8 - q^9 + q^10 - q^11 + q^12 - q^13 + q^14 - q^15 + q^16 - q^17 + q^18 - q^19 + O(q^20)
            sage: prec = R.default_prec(); prec
            20
            sage: R.set_default_prec(5)
            sage: 1/(1+q)
            1 - q + q^2 - q^3 + q^4 + O(q^5)

            sage: 1/(q + q^2)
             q^-1 - 1 + q - q^2 + q^3 + O(q^4)
            sage: g = 1/(q + q^2 + O(q^5))
            sage: g; g.parent()
            q^-1 - 1 + q - q^2 + O(q^3)
            Laurent Series Ring in q over Rational Field

            sage: 1/g
            q + q^2 + O(q^5)
            sage: (1/g).parent()
            Laurent Series Ring in q over Rational Field

            sage: 1/(2 + q)
             1/2 - 1/4*q + 1/8*q^2 - 1/16*q^3 + 1/32*q^4 + O(q^5)

            sage: R.<q> = QQ[['q']]
            sage: R.set_default_prec(5)
            sage: f = 1 + q + q^2 + O(q^50)
            sage: f/10
            1/10 + 1/10*q + 1/10*q^2 + O(q^50)
            sage: f/(10+q)
            1/10 + 9/100*q + 91/1000*q^2 - 91/10000*q^3 + 91/100000*q^4 + O(q^5)

            sage: R.<t> = PowerSeriesRing(QQ, sparse=True)
            sage: u = 17 + 3*t^2 + 19*t^10 + O(t^12)
            sage: v = ~u; v
            1/17 - 3/289*t^2 + 9/4913*t^4 - 27/83521*t^6 + 81/1419857*t^8 - 1587142/24137569*t^10 + O(t^12)
            sage: u*v
            1 + O(t^12)

        AUTHORS:
            -- David Harvey (2006-09-09): changed to use Newton's method

        """
        if self == 1:
            return self
        prec = self.prec()
        if prec is infinity and self.degree() > 0:
            prec = self._parent.default_prec()
        if self.valuation() > 0:
            u = ~self.valuation_zero_part()    # inverse of unit part
            R = self._parent.laurent_series_ring()
            return R(u, -self.valuation())

        # Use Newton's method, i.e. start with single term approximation,
        # and then iteratively compute $x' = 2x - Ax^2$, where $A$ is the
        # series we're trying to invert.

        try:
            first_coeff = ~self[0]
        except ValueError, ZeroDivisionError:
            raise ZeroDivisionError, "leading coefficient must be a unit"

        if prec is infinity:
            return self._parent(first_coeff, prec=prec)

        A = self.truncate()
        R = A.parent()     # R is the corresponding polynomial ring
        current = R(first_coeff)

        # todo: in the case that the underlying polynomial ring is
        # implemented via NTL, the truncate() method should use NTL's
        # methods. Currently it is very slow because it uses generic code
        # that has to pull all the data in and out of the polynomials.

        # todo: also, NTL has built-in series inversion. We should use
        # that when available.

        for next_prec in sage.misc.misc.newton_method_sizes(prec)[1:]:
            z = current.square() * A.truncate(next_prec)
            current = 2*current - z.truncate(next_prec)

        return self._parent(current, prec=prec)

        # Here is the old code, which uses a simple recursion, and is
        # asymptotically inferior:
        #
        #a0 = self[0]
        #try:
        #    b = [~a0]
        #except ValueError, ZeroDivisionError:
        #    raise ZeroDivisionError, "leading coefficient must be a unit"

        ##  By multiplying through we may assume that the leading coefficient
        ##  of f=self is 1.   If f = 1 + a_1*q + a_2*q^2 + ..., then we have
        ##  the following recursive formula for the coefficients b_n of the
        ##  expansion of f^(-1):
        ##        b_n = -b_0*(b_{n-1}*a_1 + b_{n-2}*a_2 + ... + b_0 a_n).
        #if self.degree() > 0:
        #    a = self.list()
        #    for n in range(1,prec):
        #        b.append(-b[0]*sum([b[n-i]*a[i] for i in range(1,n+1) if i < len(a)]))
        #return self.parent()(b, prec=prec)

    def valuation_zero_part(self):
        r"""
        Factor self as as $q^n\cdot (a_0 + a_1 q + \cdots)$ with $a_0$
        nonzero.  Then this function returns $a_0 + a_1 q + \cdots $.

        NOTE: this valuation zero part need not be a unit if, e.g.,
        $a_0$ is not invertible in the base ring.

        EXAMPLES:
            sage: R.<t> = PowerSeriesRing(QQ)
            sage: ((1/3)*t^5*(17-2/3*t^3)).valuation_zero_part()
            17/3 - 2/9*t^3

        In this example the valuation 0 part is not a unit:
            sage: R.<t> = PowerSeriesRing(ZZ, sparse=True)
            sage: u = (-2*t^5*(17-t^3)).valuation_zero_part(); u
            -34 + 2*t^3
            sage: u.is_unit()
            False
            sage: u.valuation()
            0
        """
        if self.is_zero():
            raise ValueError, "power series has no valuation 0 part"
        n = self.valuation()
        if n == 0:
            return self
        elif self.is_dense():
            v = self.list()[int(n):]
        else:
            n = int(n)
            v = {}
            for k, x in self.dict().iteritems():
                if k >= n:
                    v[k-n] = x
        return self._parent(v, self.prec()-n)

    cdef RingElement _div_c_impl(self, RingElement denom_r):
        """
        EXAMPLES:
            sage: k.<t> = QQ[[]]
            sage: t/t
            1
            sage: (t/(t^3 + 1)) * (t^3 + 1)
            t + O(t^21)
            sage: (t^5/(t^2 - 2)) * (t^2 -2 )
            t^5 + O(t^25)
        """
        denom = <PowerSeries>denom_r
        if denom.is_zero():
            raise ZeroDivisionError, "Can't divide by something indistinguishable from 0"
        u = denom.valuation_zero_part()
        inv = ~u  # inverse

        v = denom.valuation()
        if v > self.valuation():
            R = self._parent.laurent_series_ring()
            return R(self)/R(denom)

        # Algorithm: Cancel common factors of q from top and bottom,
        # then invert the denominator.  We do the cancellation first
        # because we can only invert a unit (and remain in the ring
        # of power series).

        if v > 0:
            num = self >> v
        else:
            num = self
        return num*inv

    def __mod__(self, other):
        """
        EXAMPLES:
            sage: R.<T> = Qp(7)[[]]
            sage: f = (48*67 + 46*67^2)*T + (1 + 42*67 + 5*67^3)*T^2 + O(T^3)
            sage: f % 67
            T^2 + O(T^3)
        """
        if isinstance(other,(int,Integer,long)):
            return power_series_ring.PowerSeriesRing(IntegerModRing(other), self.variable())(self)
        raise NotImplementedError, "Mod on power series ring elements not defined except modulo an integer."

    def shift(self, n):
        r"""
        Returns this power series multiplied by the power $t^n$. If $n$
        is negative, terms below $t^n$ will be discarded. Does not
        change this power series.

        NOTE:
            Despite the fact that higher order terms are printed to the
            right in a power series, right shifting decreases the powers
            of $t$, while left shifting increases them. This is to be
            consistant with polynomials, integers, etc.

        EXAMPLES:
            sage: R.<t> = PowerSeriesRing(QQ['y'], 't', 5)
            sage: f = ~(1+t); f
            1 + -t + t^2 + -t^3 + t^4 + O(t^5)
            sage: f.shift(3)
            t^3 + -t^4 + t^5 + -t^6 + t^7 + O(t^8)
            sage: f >> 2
            1 + -t + t^2 + O(t^3)
            sage: f << 10
            t^10 + -t^11 + t^12 + -t^13 + t^14 + O(t^15)
            sage: t << 29
            t^30

        AUTHOR:
            -- Robert Bradshaw (2007-04-18)
        """
        return self._parent(self.polynomial().shift(n), self._prec + n)

    def __lshift__(self, n):
        return self.parent()(self.polynomial() << n, self.prec())

    def __rshift__(self, n):
        return self.parent()(self.polynomial() >> n, self.prec())

    def is_square(self):
        """
        Returns True if this function has a square root in this ring,
        e.g. there is an element $y$ in \code{self.parent()} such that
        $y^2 = \code{self}$.

        ALGORITHM:
            If the basering is a field, this is true whenver the power
            series has even valuation and the leading coefficent is a
            perfect square.

            For an integral domain, it operates attempts the square root
            in the fraction field and tests whether or not the result
            lies in the original ring.

        EXAMPLES:
            sage: K.<t> = PowerSeriesRing(QQ, 't', 5)
            sage: (1+t).is_square()
            True
            sage: (2+t).is_square()
            False
            sage: (2+t.change_ring(RR)).is_square()
            True
            sage: t.is_square()
            False
            sage: K.<t> = PowerSeriesRing(ZZ, 't', 5)
            sage: (1+t).is_square()
            False
            sage: f = (1+t)^100
            sage: f.is_square()
            True

        """
        val = self.valuation()
        if val is not infinity and val % 2 == 1:
            return False
        elif not self[val].is_square():
            return False
        elif is_Field(self.base_ring()):
            return True
        else:
            try:
                self.parent()(self.sqrt())
                return True
            except TypeError:
                return False

    def sqrt(self):
        r"""
        Return the square root of self, up to the precision of parent.
        The leading power of x must be even.

        ALGORITHM:
            Newton's method
            $x_{i+1} = \frac{1}{2}( x_i + self/x_i )$

        EXAMPLES:
            sage: K.<t> = PowerSeriesRing(QQ, 't', 5)
            sage: sqrt(t^2)
            t
            sage: sqrt(1+t)
            1 + 1/2*t - 1/8*t^2 + 1/16*t^3 - 5/128*t^4 + O(t^5)
            sage: sqrt(4+t)
            2 + 1/4*t - 1/64*t^2 + 1/512*t^3 - 5/16384*t^4 + O(t^5)
            sage: sqrt(2+t)
            1.41421356237309 + 0.353553390593274*t - 0.0441941738241592*t^2 + 0.0110485434560399*t^3 - 0.00345266983001233*t^4 + O(t^5)

            sage: K.<t> = PowerSeriesRing(QQ, 't', 50)
            sage: sqrt(1+2*t+t^2)
            1 + t + O(t^50)
            sage: sqrt(t^2 +2*t^4 + t^6)
            t + t^3 + O(t^51)
            sage: sqrt(1 + t + t^2 + 7*t^3)^2
            1 + t + t^2 + 7*t^3 + O(t^50)
            sage: sqrt(K(0))
            0
            sage: sqrt(t^2)
            t

        AUTHOR:
            -- Robert Bradshaw
        """
        if self.is_zero():
            return self._parent(0).O(self.prec()/2)

        val = self.valuation()
        if val is not infinity and val % 2 == 1:
            raise ValueError, "Square root not defined for power series of odd valuation."

        if self.degree() == 0:
            x = self.valuation_zero_part()[0].sqrt()
            return self.parent().base_extend(x.parent())([x], val/2)

        prec = self.prec()
        if prec == infinity:
            prec = self._parent.default_prec()

        a = self.valuation_zero_part()
        x = a[0].sqrt()
        newp = self._parent.base_extend(x.parent())
        a = newp(a)
        half = ~newp.base_ring()(2)

        for i in range (ceil(log(prec, 2))):
            x = half * (x + a/x)

        return newp.gen(0)**(val/2) * x

    def square_root(self):
        """
        Return the square root of self in this ring. If this cannot be done
        then an error will be raised.

        This function succeeds if and only if \code{self.is_square()}

        EXAMPLES:
            sage: K.<t> = PowerSeriesRing(QQ, 't', 5)
            sage: (1+t).square_root()
            1 + 1/2*t - 1/8*t^2 + 1/16*t^3 - 5/128*t^4 + O(t^5)
            sage: (2+t).square_root()
            Traceback (most recent call last):
            ...
            ValueError: Square root does not live in this ring.
            sage: (2+t.change_ring(RR)).square_root()
            1.41421356237309 + 0.353553390593274*t - 0.0441941738241592*t^2 + 0.0110485434560399*t^3 - 0.00345266983001233*t^4 + O(t^5)
            sage: t.square_root()
            Traceback (most recent call last):
            ...
            ValueError: Square root not defined for power series of odd valuation.
            sage: K.<t> = PowerSeriesRing(ZZ, 't', 5)
            sage: f = (1+t)^20
            sage: f.square_root()
            1 + 10*t + 45*t^2 + 120*t^3 + 210*t^4 + O(t^5)
            sage: f = 1+t
            sage: f.square_root()
            Traceback (most recent call last):
            ...
            ValueError: Square root does not live in this ring.

        AUTHOR:
            -- Robert Bradshaw
        """
        val = self.valuation()
        if val is not infinity and val % 2 == 1:
            raise ValueError, "Square root not defined for power series of odd valuation."
        elif not self[val].is_square():
            raise ValueError, "Square root does not live in this ring."
        elif is_Field(self.base_ring()):
            return self.sqrt()
        else:
            try:
                return self.parent()(self.sqrt())
            except TypeError:
                raise ValueError, "Square root does not live in this ring."

    def O(self, prec):
        r"""
        Return this series plus $O(x^\text{prec})$.  Does not change
        self.
        """
        if prec is infinity or prec >= self.prec():
            return self
        coeffs = self[:prec]
        return self._parent(coeffs, prec)


    def solve_linear_de(self, prec = infinity, b = None, f0 = None):
        r"""
        Obtains a power series solution to an inhomogeneous linear
        differential equation of the form:
           $$  f'(t) = a(t) f(t) + b(t). $$

        INPUT:
            self -- the power series $a(t)$
            b -- the power series $b(t)$ (default is zero)
            f0 -- the constant term of $f$ (``initial condition'')
                 (default is 1)
            prec -- desired precision of result (this will be reduced if
                    either a or b have less precision available)

        OUTPUT:
            the power series f, to indicated precision

        ALGORITHM:
            A divide-and-conquer strategy; see the source code. Running time
            is approximately $M(n) \log n$, where $M(n)$ is the time required
            for a polynomial multiplication of length $n$ over the coefficient
            ring. (If you're working over something like RationalField(),
            running time analysis can be a little complicated because the
            coefficients tend to explode.)

        NOTES:
            -- If the coefficient ring is a field of characteristic zero,
               then the solution will exist and is unique.
            -- For other coefficient rings, things are more complicated.
               A solution may not exist, and if it does it may not be unique.
               Generally, by the time the nth term has been computed, the
               algorithm will have attempted divisions by $n!$ in the
               coefficient ring. So if your coefficient ring has enough
               ``precision'', and if your coefficient ring can perform
               divisions even when the answer is not unique, and if you know
               in advance that a solution exists, then this function will find
               a solution (otherwise it will probably crash).

        AUTHORS:
          -- David Harvey (2006-09-11): factored functionality out from
             exp() function, cleaned up precision tests a bit

        EXAMPLES:
           sage: R.<t> = PowerSeriesRing(QQ, default_prec=10)

           sage: a = 2 - 3*t + 4*t^2 + O(t^10)
           sage: b = 3 - 4*t^2 + O(t^7)
           sage: f = a.solve_linear_de(prec=5, b=b, f0=3/5)
           sage: f
            3/5 + 21/5*t + 33/10*t^2 - 38/15*t^3 + 11/24*t^4 + O(t^5)
           sage: f.derivative() - a*f - b
            O(t^4)

           sage: a = 2 - 3*t + 4*t^2
           sage: b = b = 3 - 4*t^2
           sage: f = a.solve_linear_de(b=b, f0=3/5)
           Traceback (most recent call last):
           ...
           ValueError: cannot solve differential equation to infinite precision

           sage: a.solve_linear_de(prec=5, b=b, f0=3/5)
            3/5 + 21/5*t + 33/10*t^2 - 38/15*t^3 + 11/24*t^4 + O(t^5)
        """
        if b is None:
            b = self._parent(0)

        if f0 is None:
            f0 = 1
        f0 = self.base_ring()(f0)

        # reduce precision to whatever is available from a and b
        prec = min(prec, self.prec() + 1, b.prec() + 1)
        if prec is infinity:
            raise ValueError, \
                 "cannot solve differential equation to infinite precision"
        prec = int(prec)
        if prec == 0:
            return self._parent(0, 0)
        if prec < 0:
            raise ValueError, \
                  "prec (=%s) must be a non-negative integer" % prec

        base_ring = self._parent.base_ring()
        R = PolynomialRing(base_ring, self._parent.variable_name())

        a_list = self.list()
        b_list = [base_ring(0)]
        b_list.extend(b.list())
        while len(b_list) < prec:
            b_list.append(base_ring(0))
        f = _solve_linear_de(R, 0, prec, a_list, b_list, f0)
        return self._parent(f, prec)


    def exp(self, prec=None):
        r"""
        Returns exp of this power series to the indicated precision.

        INPUT:
            prec -- integer; default is self.parent().default_prec

        ALGORITHM:
            See PowerSeries.solve_linear_de().

        NOTES:
            -- Screwy things can happen if the coefficient ring is not
               a field of characteristic zero.
               See PowerSeries.solve_linear_de().

        AUTHORS:
          -- David Harvey (2006-09-08): rewrote to use simplest possible
             ``lazy'' algorithm.
          -- David Harvey (2006-09-10): rewrote to use divide-and-conquer
             strategy.
          -- David Harvey (2006-09-11): factored functionality out to
             solve_linear_de().

        EXAMPLES:
           sage: R.<t> = PowerSeriesRing(QQ, default_prec=10)

         Check that $\exp(t)$ is, well, $\exp(t)$:
           sage: (t + O(t^10)).exp()
           1 + t + 1/2*t^2 + 1/6*t^3 + 1/24*t^4 + 1/120*t^5 + 1/720*t^6 + 1/5040*t^7 + 1/40320*t^8 + 1/362880*t^9 + O(t^10)

         Check that $\exp(\log(1+t))$ is $1+t$:
           sage: (sum([-(-t)^n/n for n in range(1, 10)]) + O(t^10)).exp()
           1 + t + O(t^10)

         Check that $\exp(2t + t^2 - t^5)$ is whatever it is:
           sage: (2*t + t^2 - t^5 + O(t^10)).exp()
           1 + 2*t + 3*t^2 + 10/3*t^3 + 19/6*t^4 + 8/5*t^5 - 7/90*t^6 - 538/315*t^7 - 425/168*t^8 - 30629/11340*t^9 + O(t^10)

         Check requesting lower precision:
           sage: (t + t^2 - t^5 + O(t^10)).exp(5)
           1 + t + 3/2*t^2 + 7/6*t^3 + 25/24*t^4 + O(t^5)

         Can't get more precision than the input:
           sage: (t + t^2 + O(t^3)).exp(10)
           1 + t + 3/2*t^2 + O(t^3)

         Check some boundary cases:
           sage: (t + O(t^2)).exp(1)
           1 + O(t)
           sage: (t + O(t^2)).exp(0)
           O(t^0)

        """
        if prec is None:
            prec = self._parent.default_prec()
        if not self[0].is_zero():
            raise ValueError, "constant term must to zero"
        return self.derivative().solve_linear_de(prec)


    def V(self, n):
        """
        If $f = \sum a_m x^m$, then this function returns $\sum a_m x^{nm}$.
        """
        v = self.list()
        m = 0
        w = []
        zero = self.base_ring()(0)
        for i in range(len(v)*n):
            if i%n != 0:
                w.append(zero)
            else:
                w.append(v[m])
                m += 1
        return self._parent(w, self.prec()*n)

    def valuation(self):
        """
        Return the valuation of this power series.

        This is equal to the valuation of the underlying polynomial.

        EXAMPLES:
        Sparse examples:
            sage: R.<t> = PowerSeriesRing(QQ, sparse=True)
            sage: f = t^100000 + O(t^10000000)
            sage: f.valuation()
            100000
            sage: R(0).valuation()
            +Infinity

        Dense examples:
            sage: R.<t> = PowerSeriesRing(ZZ)
            sage: f = 17*t^100 +O(t^110)
            sage: f.valuation()
            100
            sage: t.valuation()
            1
        """
        return self.polynomial().valuation()

    def variable(self):
        """
        EXAMPLES:
            sage: R.<x> = PowerSeriesRing(Rationals())
            sage: f = x^2 + 3*x^4 + O(x^7)
            sage: f.variable()
            'x'

        AUTHOR:
            -- David Harvey (2006-08-08)
        """
        return self._parent.variable_name()

    def degree(self):
        """
        Return the degree of this power series, which is by definition
        the degree of the underlying polynomial.

        EXAMPLES:
            sage: R.<t> = PowerSeriesRing(QQ, sparse=True)
            sage: f = t^100000 + O(t^10000000)
            sage: f.degree()
            100000
        """
        return self.polynomial().degree()

    def derivative(self):
        raise NotImplementedError

#def _normalize(v, prec):
#    n = len(v)-1
#    while n>=0 and (v[n] == 0 or (prec != infinity and n >= prec)):
#        del v[n]
#        n -= 1

cdef class PowerSeries_poly(PowerSeries):

    def __init__(self, parent, f=0, prec=infinity, int check=1, is_gen=0):
        """
        EXAMPLES:
            sage: R, q = PowerSeriesRing(CC, 'q').objgen()
            sage: R
            Power Series Ring in q over Complex Field with 53 bits of precision
            sage: loads(q.dumps()) == q
            True
        """
        R = parent._poly_ring()
        if PY_TYPE_CHECK(f, Element):
            if (<Element>f)._parent is R:
                pass
            elif (<Element>f)._parent is R.base_ring():
                f = R([f])
            elif PY_TYPE_CHECK(f, PowerSeries_poly):
                prec = (<PowerSeries_poly>f)._prec
                f = R((<PowerSeries_poly>f).__f)
            else:
                f = R(f, check=check)
        else:
            f = R(f, check=check)

        self.__f = f
        if check and not (prec is infinity):
            self.__f = self.__f.truncate(prec)
        PowerSeries.__init__(self, parent, prec, is_gen)

    def __reduce__(self):
        # TODO: find out what's going on with ring pickling
        params = (self._parent.base_ring(), self._parent.variable_name(), self._parent.default_prec(), self._parent.variable_names(), self._parent.is_sparse())
        return make_powerseries_poly, (params, self.__f, self._prec, 0, self.__is_gen)

    def __richcmp__(left, right, int op):
        return (<Element>left)._richcmp(right, op)

    def __pow__(self_t, r, dummy):
        cdef PowerSeries_poly self = self_t
        cdef int right = r
        if right != r:
            raise ValueError, "exponent must be an integer"
        if right < 0:
            return (~self)**(-right)
        if right == 0:
            return self._parent(1)
        if self.__is_gen:
            return PowerSeries_poly(self._parent, self.__f**right, check=False)
        if self.is_zero():
            return self
        return arith.generic_power(self, right, self._parent(1))

    def polynomial(self):
        """
        EXAMPLE:
            sage: R.<t> = GF(7)[[]]
            sage: f = 3 - t^3 + O(t^5)
            sage: f.polynomial()
            6*t^3 + 3
        """
        return self.__f

    def valuation(self):
        return self.__f.valuation()

    def degree(self):
        return self.__f.degree()

    def __call__(self, *xs):
        """
        EXAMPLE:
            sage: R.<t> = GF(7)[[]]
            sage: f = 3 - t^3 + O(t^5)
            sage: f(1)
            2
        """
        if isinstance(xs[0], tuple):
            xs = xs[0]
        x = xs[0]
        try:
            if x.parent() is self._parent:
                if not (self.prec() is infinity):
                    x = x.add_bigoh(self.prec()*x.valuation())
                    xs = list(xs); xs[0] = x; xs = tuple(xs) # tuples are immutable
        except AttributeError:
            pass
        return self.__f(xs)

    def __setitem__(self, n, value):
        """
        EXAMPLES:
            sage: R.<t> = ZZ[[]]
            sage: f = 3 - t^3 + O(t^5)
            sage: f[1] = 5
            Traceback (most recent call last):
            ...
            IndexError: power series are immutable
        """
        raise IndexError, "power series are immutable"

    def __getslice__(self, i, j):
        r"""
        Return slice of coefficient of this power series.

        This calls slice on the underlying polynomial, and makes a power
        series out of the result, with precision the precision of self.

        EXAMPLES:
            sage: R.<t> = ZZ[[]]
            sage: f = (2-t)^5 + O(t^7); f
            32 - 80*t + 80*t^2 - 40*t^3 + 10*t^4 - t^5 + O(t^7)
            sage: f[2:4]
            80*t^2 - 40*t^3 + O(t^7)
        """
        return PowerSeries_poly(self._parent, self.__f[i:j], prec=self.prec(), check=False)

    def _unsafe_mutate(self, i, value):
        """
        SAGE assumes throughout that commutative ring elements are immutable.
        This is relevant for caching, etc.  But sometimes you need to change
        a power series and you really know what you're doing.  That's
        when this function is for you.

        ** DO NOT USE THIS ** unless you know what you're doing.

        EXAMPLES:
            sage: R.<t> = GF(7)[[]]
            sage: f = 3 + 6*t^3 + O(t^5)
            sage: f._unsafe_mutate(0, 5)
            sage: f
            5 + 6*t^3 + O(t^5)

        Mutating can even bump up the precision.
            sage: f._unsafe_mutate(7,2)
            sage: f
            5 + 6*t^3 + 2*t^7 + O(t^8)
        """
        self.__f._unsafe_mutate(i, value)
        self._prec = max(self._prec, i+1)

    def __getitem__(self, n):
        """
        Return the n-th coefficient.

        Returns 0 for negative coefficients.  Raises an IndexError if
        try to access beyond known coefficients.

        EXAMPLES:
            sage: R.<t> = QQ[[]]
            sage: f = 3/2 - 17/5*t^3 + O(t^5)
            sage: f[3]
            -17/5
            sage: f[-2]
            0
            sage: f[4]
            0
            sage: f[5]
            Traceback (most recent call last):
            ...
            IndexError: coefficient not known
            sage: f[1:4]
            -17/5*t^3 + O(t^5)
        """
        if n<0:
            return self.base_ring()(0)
        if n > self.__f.degree():
            if self._prec > n:
                return self.base_ring()(0)
            #elif isinstance(n, slice):
                # It makes no sense that this is needed and that
                # __getslice__ isn't just called by default...
            #    return self.__getslice__(slice[0],slice[1])
            else:
                raise IndexError, "coefficient not known"
        return self.__f[n]

    def __iter__(self):
        """
        Return an interator over the coefficients of this power series.

        EXAMPLES:
            sage: R.<t> = QQ[[]]
            sage: f = t + 17/5*t^3 + 2*t^4 + O(t^5)
            sage: for a in f: print a,
            0 1 0 17/5 2
        """
        return iter(self.__f)

    def __neg__(self):
        """
        Return the negative of this power series.

        EXAMPLES:
            sage: R.<t> = QQ[[]]
            sage: f = t + 17/5*t^3 + 2*t^4 + O(t^5)
            sage: -f
            -t - 17/5*t^3 - 2*t^4 + O(t^5)
        """
        return PowerSeries_poly(self._parent, -self.__f,
                                         self._prec, check=False)

    cdef ModuleElement _add_c_impl(self, ModuleElement right_m):
        """
        EXAMPLES:
            sage: R.<x> = PowerSeriesRing(ZZ)
            sage: f = x^4 + O(x^5); f
            x^4 + O(x^5)
            sage: g = x^2 + O(x^3); g
            x^2 + O(x^3)
            sage: f+g
            x^2 + O(x^3)
        """
        cdef PowerSeries_poly right = <PowerSeries_poly>right_m
        return PowerSeries_poly(self._parent, self.__f + right.__f, \
                                         self.common_prec_c(right), check=True)

    cdef ModuleElement _sub_c_impl(self, ModuleElement right_m):
        """
        Return difference of two power series.

        EXAMPLES:
            sage: k.<w> = ZZ[]
            sage: R.<t> = k[[]]
            sage: w*t^2 -w*t +13 - (w*t^2 + w*t)
            13 + -2*w*t
        """
        cdef PowerSeries_poly right = <PowerSeries_poly>right_m
        return PowerSeries_poly(self._parent, self.__f - right.__f, \
                                         self.common_prec_c(right), check=True)

    cdef RingElement _mul_c_impl(self, RingElement right_r):
        """
        Return the product of two power series.

        EXAMPLES:
            sage: k.<w> = ZZ[[]]
            sage: (1+17*w+15*w^3+O(w^5))*(19*w^10+O(w^12))
            19*w^10 + 323*w^11 + O(w^12)
        """
        cdef PowerSeries_poly right = <PowerSeries_poly>right_r
        sp = self._prec
        rp = right._prec
        if is_Infinite(sp):
            if is_Infinite(rp):
                prec = infinity
            else:
                prec = rp + self.valuation()
        else:  # sp != infinity
            if is_Infinite(rp):
                prec = sp + right.valuation()
            else:
                prec = min(rp + self.valuation(), sp + right.valuation())
        return PowerSeries_poly(self._parent,
                                         self.__f * right.__f,
                                         prec,
                                         check=True)  # check, since truncation may be needed

    def _rmul_(self, c):
        return PowerSeries_poly(self._parent, c * self.__f, self._prec, check=False)

    def _lmul_(self, c):
        return PowerSeries_poly(self._parent, self.__f * c, self._prec, check=False)


    def __floordiv__(self, denom):
        try:
            return PowerSeries.__div__(self, denom)
        except (PariError, ZeroDivisionError), e: # PariError to general?
            if is_PowerSeries(denom) and denom.degree() == 0 and denom[0] in self._parent.base_ring():
                denom = denom[0]
            elif not denom in self._parent.base_ring():
                raise ZeroDivisionError, e
            return PowerSeries_poly(self._parent,
                                             self.__f // denom, self._prec)

    def __lshift__(_self, n):
        cdef PowerSeries_poly self = _self
        return PowerSeries_poly(self._parent, self.__f << n, self._prec + n)

    def __rshift__(_self, n):
        cdef PowerSeries_poly self = _self
        return PowerSeries_poly(self._parent, self.__f >> n, self._prec - n)

    def copy(self):
        return PowerSeries_poly(self._parent, self.__f, self.prec(), check=False)

    def list(self):
        return self.__f.list()

    def dict(self):
        return self.__f.dict()

    def derivative(self):
        """
        Return the derivative of this power series.

        EXAMPLES:
            sage: R.<t> = PowerSeriesRing(QQ, sparse=True)
            sage: f = 2 + 3*t^2 + t^100000 + O(t^10000000); f
            2 + 3*t^2 + t^100000 + O(t^10000000)
            sage: f.derivative()
            6*t + 100000*t^99999 + O(t^9999999)
        """
        return PowerSeries_poly(self._parent, self.__f.derivative(),
                                         self.prec()-1, check=False)

    def integral(self):
        """
        The integral of this power series with 0 constant term.

        EXAMPLES:
            sage: k.<w> = QQ[[]]
            sage: (1+17*w+15*w^3+O(w^5)).integral()
            w + 17/2*w^2 + 15/4*w^4 + O(w^6)
        """
        return PowerSeries_poly(self._parent, self.__f.integral(),
                                         self.prec()+1, check=False)

    def laurent_series(self):
        """
        Return the Laurent series associated to this power series, i.e., this
        series considered as a Laurent series.

        EXAMPLES:
            sage: k.<w> = QQ[[]]
            sage: f = 1+17*w+15*w^3+O(w^5)
            sage: parent(f)
            Power Series Ring in w over Rational Field
            sage: g = f.laurent_series(); g
            1 + 17*w + 15*w^3 + O(w^5)
        """
        return self._parent.laurent_series_ring()(self)

    def ogf(self):
        r"""
        Returns the ordinary generating function associated to self.

        This function is known as serlaplace in GP/PARI.

        EXAMPLES:
            sage: R.<t> = PowerSeriesRing(QQ)
            sage: f = t + t^2/factorial(2) + 2*t^3/factorial(3)
            sage: f.ogf()
            t + t^2 + 2*t^3
        """
        return self.parent()([self[i] * arith.factorial(i) for i in range(self.degree()+1)])

    def egf(self):
        r"""
        Returns the exponential generating function associated to self.

        This function is known as serlaplace in GP/PARI.

        EXAMPLES:
            sage: R.<t> = PowerSeriesRing(QQ)
            sage: f = t + t^2 + 2*t^3
            sage: f.egf()
            t + 1/2*t^2 + 1/3*t^3
        """
        return self.parent()([self[i] / arith.factorial(i) for i in range(self.degree()+1)])

    def reversion(self):
        """
        Return the reversion of f, i.e., the series g such that
        g(f(x)) = x.

        EXAMPLES:
            sage: R.<x> = PowerSeriesRing(QQ)
            sage: f = 2*x + 3*x**2 - x**4 + O(x**5)
            sage: g = f.reversion()
            sage: g
            1/2*x - 3/8*x^2 + 9/16*x^3 - 131/128*x^4 + O(x^5)
            sage: f(g)
            x + O(x^5)
            sage: g(f)
            x + O(x^5)
        """
        if not isinstance(self.parent().base_ring(), rational_field.RationalField):
            raise NotImplementedError
        if self.prec() is infinity:
            raise RuntimeError, "series must have finite precision for reversion."
        f = self._pari_()
        g = f.serreverse()
        return PowerSeries_poly(self.parent(),g.Vecrev(),self.prec())

    def _pari_(self):
        """
        Return PARI power series corresponding to this series.

        This is currently only implemented over QQ and ZZ.

        EXAMPLES:
            sage: k.<w> = QQ[[]]
            sage: f = 1+17*w+15*w^3+O(w^5)
            sage: pari(f)
            1 + 17*w + 15*w^3 + O(w^5)
            sage: pari(1 - 19*w + w^5)
            Traceback (most recent call last):
            ...
            RuntimeError: series precision must be finite for conversion to pari object.
        """
        if not isinstance(self.parent().base_ring(),
                          (rational_field.RationalField, integer_ring.IntegerRing)):
            raise NotImplementedError
        if self.prec() is infinity:
            raise RuntimeError, "series precision must be finite for conversion to pari object."
        return sage.libs.pari.all.pari(str(self))



def _solve_linear_de(R, N, L, a, b, f0):
    r"""
    Internal function used by PowerSeries.solve_linear_de().

    INPUT:
       R -- a PolynomialRing
       N -- integer >= 0
       L -- integer >= 1
       a -- list of coefficients of $a$, any length, all coefficients should
            belong to base ring of R.
       b -- list of coefficients of $b$, length at least $L$
            (only first $L$ coefficients are used), all coefficients
            should belong to base ring of R.
       f0 -- constant term of $f$ (only used if $N == 0$), should belong
            to base ring of R.

    OUTPUT:
       List of coefficients of $f$ (length exactly $L$), where $f$ is a
       solution to the linear inhomogeneous differential equation:
          $$ (t^N f)'  =  a t^N f  +  t^{N-1} b  +  O(t^{N+L-1}).$$
       The constant term of $f$ is taken to be f0 if $N == 0$, and otherwise
       is determined by $N f_0 = b_0$.

    ALGORITHM:
        The algorithm implemented here is inspired by the ``fast lazy
        multiplication algorithm'' described in the paper ``Lazy multiplication
        of formal power series'' by J van der Hoeven (1997 ISSAC conference).

        Our situation is a bit simpler than the one described there,
        because only one of the series being multiplied needs the lazy
        treatment; the other one is known fully in advance. I have
        reformulated the algorithm as an explicit divide-and-conquer
        recursion. Possibly it is slower than the one described by van der
        Hoeven, by a constant factor, but it seems easier to implement.
        Also, because we repeatedly split in half starting at the top level,
        instead of repeatedly doubling in size from the bottom level, it's
        easier to avoid problems with ``overshooting'' in the last iteration.

        The idea is to split the problem into two instances with $L$
        about half the size. Take $L'$ to be the ceiling of
        $(L/2)$. First recursively find $g$ modulo $t^{L'}$ such that
          $$ (t^N g)'  =  a t^N g  +  t^{N-1} b  +  O(t^{N+L'-1}).$$

        Next we want to find $h$ modulo $t^{L-L'}$ such that
        $f = g + t^{L'} h$ is a solution of the original equation.
        We can find a suitable $h$ by recursively solving another
        differential equation of the same form, namely
          $$ (t^{N+L'} h)'  =  a t^{N+L'} h  +  c t^{N+L'-1} + O(t^{N+L'-1}),$$
        where $c$ is determined by
          $$ (t^N g)' - a t^N g - t^{N-1} b  =  -c t^{N+L'-1} + O(t^{N+L-1}).$$
        Once $g$ is known, $c$ can be recovered from this relation by computing
        the coefficients of $t^j$ of the product $a g$ for $j$ in the range
        $L'-1 \leq j < L-1$.

        At the bottom of the recursion we have $L = 1$, in which case
        the solution is simply given by $f_0 = b_0/N$ (or by the supplied
        initial condition $f_0$ when $N == 0$).

        The algorithm has to do one multiplication of length $N$, two of
        length $N/2$, four of length $N/4$, etc, hence giving runtime
        $O(M(N) \log N)$.

    AUTHOR:
        -- David Harvey (2006-09-11): factored out of PowerSeries.exp().

    """
    if L == 1:
        # base case
        if N == 0:
            return [f0]
        else:
            return [b[0] / N]

    L2 = (L + 1) >> 1    # ceil(L/2)

    g = _solve_linear_de(R, N, L2, a, b, f0)

    term1 = R(g, check=False)
    term2 = R(a[:L], check=False)
    product = (term1 * term2).list()

    # todo: perhaps next loop could be made more efficient
    c = b[L2 : L]
    for j in range(L2 - 1, min(L-1, len(product))):
        c[j - (L2-1)] = c[j - (L2-1)] + product[j]

    h = _solve_linear_de(R, N + L2, L - L2, a, c, f0)

    g.extend(h)
    return g


# TODO: fix pickling of PowerSeriesRing
def make_powerseries_poly(params, *args):
    return PowerSeries_poly(power_series_ring.PowerSeriesRing(*params), *args)

def make_element_from_parent(parent, *args):
    return parent(*args)