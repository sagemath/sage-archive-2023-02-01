"""
Power Series

AUTHOR:
   -- William Stein
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

import operator

from infinity import infinity
import polynomial_element as polynomial
import power_series_ring
import sage.misc.misc as misc
import ring_element
import arith
import sage.misc.latex as latex
from coerce import bin_op
import sage.ext.coerce
import rational_field, integer_ring
import sage.libs.pari.all as pari
import sage.misc.latex as latex
from sage.structure.element import Element_cmp_
from sage.libs.all import PariError


Polynomial = polynomial.Polynomial_generic_dense

def is_PowerSeries(x):
    return isinstance(x, PowerSeries)

class PowerSeries(Element_cmp_, ring_element.RingElement):
    def __init__(self, parent, prec, is_gen=False):
        ring_element.RingElement.__init__(self, parent)
        self.__is_gen = is_gen
        if not (prec is infinity):
            prec = int(prec)
        self._prec = prec

    def is_gen(self):
        return self.__is_gen

    def _im_gens_(self, codomain, im_gens):
        return codomain(self(im_gens[0]))

    def _cmp_(self, right):
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
        """
        prec = self.common_prec(right)
        x = self.list()
        y = right.list()
        if prec != infinity:
            x = x[:prec]
            y = y[:prec]
        return misc.generic_cmp(x,y)

    def __call__(self, x):
        raise NotImplementedError

    def base_ring(self):
        return self.parent().base_ring()

    def list(self):
        raise NotImplementedError

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
            z = self.parent().base_ring()(0)
            return v + [z]*(n - len(v))
        else:
            return v[:int(n)]

    def prec(self):
        """
        The precision of $...+O(x^r)$ is by definition $r$.
        """
        return self._prec

    def _repr_(self):
        if self.is_zero():
            if self.prec() == infinity:
                return "0"
            else:
                return "O(%s^%s)"%(self.parent().variable_name(),self.prec())
        s = " "
        v = self.list()
        m = len(v)
        X = self.parent().variable_name()
        atomic_repr = self.parent().base_ring().is_atomic_repr()
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

        if atomic_repr:
            s = s.replace(" + -", " - ")
        s = s.replace(" 1*"," ")
        s = s.replace(" -1*", " -")
        if self._prec == 0:
            bigoh = "O(1)"
        elif self._prec == 1:
            bigoh = "O(%s)"%self.parent().variable_name()
        else:
            bigoh = "O(%s^%s)"%(self.parent().variable_name(),self._prec)
        if self._prec != infinity:
            if s==" ":
                return bigoh
            s += " + %s"%bigoh
        return s[1:]

    def _latex_(self):
        if self.is_zero():
            if self.prec() == infinity:
                return "0"
            else:
                return "0 + \\cdots"
        s = " "
        v = self.list()
        m = len(v)
        X = self.parent().variable_name()
        atomic_repr = self.parent().base_ring().is_atomic_repr()
        first = True
        for n in xrange(m):
            x = v[n]
            if x != 0:
                if not first:
                    s += " + "
                x = latex.latex(x)
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
        if self._prec != infinity:
            if s==" ":
                return "0 + \\cdots"
            s += " + \\cdots"
        return s[1:]


    def truncate(self, prec=infinity):
        """
        The polynomial obtained from power series by truncation.
        """
        if prec == infinity:
            prec = self._prec
        a = self.list()
        v = [a[i] for i in range(min(prec, len(a)))]
        return self.parent()._poly_ring()(v)

    def add_bigoh(self, prec):
        r"""
        Returns the power series of precision at most prec got by
        adding $O(q^\text{prec})$ to f, where q is the variable.
        """
        if prec == infinity or prec >= self.prec():
            return self
        a = self.list()
        v = [a[i] for i in range(min(prec, len(a)))]
        return self.parent()(v, prec)

    def __getitem__(self,n):
        if n<0:
            return self.base_ring()(0)
        c = self.list()
        if n >= len(c):
            if self._prec > n:
                return self.base_ring()(0)
            else:
                raise IndexError, "coefficient not known"
        return c[n]

    def __getslice__(self, i, j):
        if j>self.prec():
            j = self.prec()
        return self.list()[int(i):int(j)]

    def __setitem__(self, n, value):
        raise NotImplementedError

    def common_prec(self, f):
        if self.prec() == infinity:
            return f.prec()
        elif f.prec() == infinity:
            return self.prec()
        return min(self.prec(), f.prec())

    def __add__(self):
        raise NotImplementedError

    def __radd__(self, left):
        return self.parent()(left) + self

    def __sub__(self, right):
        raise NotImplementedError

    def __mul__(self, right):
        if not isinstance(right,PowerSeries) or self.parent() != right.parent():
            return bin_op(self, right, operator.mul)
        if self.is_zero():
            return self
        if right.is_zero():
            return right
        sp = self.prec()
        rp = right.prec()
        if sp == infinity:
            if rp == infinity:
                prec = infinity
            else:
                prec = rp + self.valuation()
        else:  # sp != infinity
            if rp == infinity:
                prec = sp + right.valuation()
            else:
                prec = min(rp + self.valuation(), sp + right.valuation())
        # endif
        return self._mul_(right, prec)

    def _mul_(self, right, prec):
        raise NotImplementedError

    def is_unit(self):
        """ Returns whether this power series is invertible, which is
        the case precisely when the constant term is invertible.

        AUTHOR: David Harvey (2006-09-03)
        """
        return self[0].is_unit()

    def __invert__(self):
        """
        Inverse of the power series, which we assume to have nonzero
        constant term so that the inverse is again a power series.

        EXAMPLES:
            sage: R = PowerSeriesRing(RationalField(), 'q')
            sage: q = R.gen()
            sage: 1/(1+q + O(q**2))
            1 - q + O(q^2)
            sage: 1/(1+q)
            1 - q + q^2 - q^3 + q^4 - q^5 + q^6 - q^7 + q^8 - q^9 + q^10 - q^11 + q^12 - q^13 + q^14 - q^15 + q^16 - q^17 + q^18 - q^19 + O(q^20)
            sage: prec = R.default_prec(); prec
            20
            sage: R.set_default_prec(5)
            sage: 1/(1+q)
            1 - q + q^2 - q^3 + q^4 + O(q^5)

            sage: R.<q> = QQ[['q']]
            sage: R.set_default_prec(5)
            sage: f = 1 + q + q^2 + O(q^50)
            sage: f/10
            1/10 + 1/10*q + 1/10*q^2 + O(q^50)
            sage: f/(10+q)
            1/10 + 9/100*q + 91/1000*q^2 - 91/10000*q^3 + 91/100000*q^4 + O(q^5)
        """
        if self == 1:
            return self
        prec = self.prec()
        if prec is infinity and self.degree() > 0:
            prec = self.parent().default_prec()
        if self.valuation() > 0:
            u = ~self.unit_part()    # inverse of unit part
            R = self.parent().laurent_series_ring()
            return R(u, -self.valuation())
        a0 = self[0]
        try:
            b = [~a0]
        except ValueError, ZeroDivisionError:
            raise ZeroDivisionError, "leading coefficient must be a unit"

        #  By multiplying through we may assume that the leading coefficient
        #  of f=self is 1.   If f = 1 + a_1*q + a_2*q^2 + ..., then we have
        #  the following recursive formula for the coefficients b_n of the
        #  expansion of f^(-1):
        #        b_n = -b_0*(b_{n-1}*a_1 + b_{n-2}*a_2 + ... + b_0 a_n).
        if self.degree() > 0:
            a = self.list()
            for n in range(1,prec):
                b.append(-b[0]*sum([b[n-i]*a[i] for i in range(1,n+1) if i < len(a)]))
        return self.parent()(b, prec=prec)

    def unit_part(self):
        """
        Suppose self factors as $q^n\cdot (a_0 + a_1 q + \cdots)$
        with $a_0$ nonzero.  Then this function returns $a_0 + a_1 q +
        \cdots $.
        """
        #assert not self.is_zero(), "Argument must be nonzero to have a unit part."
        if self.is_zero():
            return self.parent()(0, self.prec())
        n = self.valuation()
        v = self[int(n):]
        return self.parent()(v, self.prec()-n)

    def __div__(self, denom):
        if not isinstance(denom, PowerSeries) or \
               self.parent() != denom.parent():

            return bin_op(self, denom, operator.div)

        u = denom.unit_part()
        inv = ~u  # inverse

        v = denom.valuation()
        if v > self.valuation():
            R = self.parent().laurent_series_ring()
            return R(self)/R(denom)

        # Algorithm: Cancel common factors of q from top and bottom,
        # then invert the denominator.  We do the cancellation first
        # because we can only invert a unit (and remain in the ring
        # of power series).

        u = denom.unit_part()
        if v > 0:
            coeffs = self[int(v):]
            num = self.parent()(coeffs, self.prec()-v)
        else:
            num = self
        return num*inv

    def __rdiv__(self, left):
        return self.parent()(left) / self

    def __pow__(self, right):
        right=int(right)
        if right < 0:
            return (~self)**(-right)
        if right == 0:
            return self.parent()(1)
        if self.is_zero():
            return self
        if self.valuation() == 1 and self.degree() == 1:   # alpha*q
            v = [0 for _ in range(right)]
            v.append(self[1]**right)
            return self.parent()(v, prec=self.prec()+right-1)
        return arith.generic_power(self, right)

    def __mod__(self, other):
        if isinstance(other,(int,Integer,long)):
            try:
                return power_series_ring.PowerSeriesRing(Zmod(other))(self)
            except TypeError:
                raise ZeroDivisionError
        raise NotImplementedError, "Mod on power series ring not defined."

    def O(self, prec):
        r"""
        Return this series plus $O(x^\text{prec})$.  Does not change
        self.
        """
        if prec == infinity or prec >= self.prec():
            return self
        coeffs = self[:prec]
        return self.parent()(coeffs, prec)

    def exp(self, prec = infinity):
        if prec == infinity:
            if self.prec() == infinity:
                raise RuntimeError, "prec must be finite for exponentiation " + \
                      "of infinite-precision power series"
            prec = self.prec()
        if self.valuation() <= 0:
            raise ArithmeticError, "valuation must be positive"
        x = self.parent()(self.list(), prec)
        ans = self.parent()(1)
        k = self.parent().base_ring()
        for n in range(1,prec):
            ans += x * (k(1)/k(arith.factorial(n)))
            x *= self
        return self.parent()(ans, prec)

    def copy(self):
        raise NotImplementedError

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
        return self.parent()(w, self.prec()*n)

    def valuation(self):
        prec = self.prec()
        v = self.list()
        if prec == infinity:
            if len(v) == 0:
                return infinity
        n = 0
        while n < len(v) and v[n] == 0:
            n += 1
        return n

    def variable(self):
        """
        EXAMPLES:
        sage: R.<x> = PowerSeriesRing(Rationals())
        sage: f = x^2 + 3*x^4 + O(x^7)
        sage: f.variable()
         'x'

        AUTHOR:
            -- David Harvey (2006-08-08); copied from LaurentSeriesRingElement
        """
        return self.parent().variable_name()

    def degree(self):
        return len(self.list())-1

    def derivative(self):
        raise NotImplementedError

#def _normalize(v, prec):
#    n = len(v)-1
#    while n>=0 and (v[n] == 0 or (prec != infinity and n >= prec)):
#        del v[n]
#        n -= 1

class PowerSeries_generic_dense(PowerSeries):
    def __init__(self, parent, f=0, prec=infinity, check=True, is_gen=False):
        """
        EXAMPLES:
            sage: R, q = PowerSeriesRing(CC, 'q').objgen()
            sage: R
            Power Series Ring in q over Complex Field with 53 bits of precision
            sage: loads(q.dumps()) == q
            True
        """
        R = parent._poly_ring()
        if not (isinstance(f, polynomial.Polynomial) and f.parent() == R):
            if isinstance(f, PowerSeries_generic_dense):
                prec = f.prec()
                f = R(f.__f)
            else:
                f = R(f, check=check)
        self.__f = f
        if check and prec != infinity:
            self.__f = self.__f.truncate(prec)
        PowerSeries.__init__(self, parent, prec, is_gen)

    def __call__(self, x):
        try:
            if x.parent() is self.parent():
                x = x.add_bigoh(self.prec()*x.valuation())
        except AttributeError:
            pass
        return self.__f(x)

    def __setitem__(self, n, value):
        self.__f[n] = value
        self._prec = max(self._prec, n+1)

    def __getitem__(self, n):
        if n<0:
            return self.base_ring()(0)
        if n > self.__f.degree():
            if self._prec > n:
                return self.base_ring()(0)
            elif isinstance(n, slice):
                # It makes no sense that this is needed and that
                # __getslice__ isn't just called by default...
                return self.__getslice__(slice[0],slice[1])
            else:
                raise IndexError, "coefficient not known"
        return self.__f[n]

    def __getslice__(self, i, j):
        return self.__f[i:j]

    def __iter__(self):
        for i in range(self.__f.degree()+1):
            yield self.__f[i]

    def __neg__(self):
        return PowerSeries_generic_dense(self.parent(), -self.__f,
                                         self._prec, check=False)

    def __add__(self, right):
        """
        EXAMPLES:
            sage: x = PowerSeriesRing(ZZ).gen()
            sage: f = x^4 + O(x^5); f
            x^4 + O(x^5)
            sage: g = x^2 + O(x^3); g
            x^2 + O(x^3)
            sage: f+g
            x^2 + O(x^3)
        """
        if not isinstance(right,PowerSeries_generic_dense) or self.parent() != right.parent():
            return bin_op(self, right, operator.add)
        return PowerSeries_generic_dense(self.parent(), self.__f + right.__f, \
                                         self.common_prec(right), check=True)

    def __sub__(self, right):
        if not isinstance(right,PowerSeries_generic_dense) or self.parent() != right.parent():
            return bin_op(self, right, operator.sub)
        return PowerSeries_generic_dense(self.parent(), self.__f - right.__f, \
                                         self.common_prec(right), check=True)

    def _mul_(self, right, prec):
        return PowerSeries_generic_dense(self.parent(),
                                         self.__f * right.__f, prec)

    def __floordiv__(self, denom):
        try:
            return PowerSeries.__div__(self, denom)
        except (PariError, ZeroDivisionError), e: # PariError to general?
            if is_PowerSeries(denom) and denom.degree() == 0 and denom[0] in self.parent().base_ring():
                denom = denom[0]
            elif not denom in self.parent().base_ring():
                raise ZeroDivisionError, e
            return PowerSeries_generic_dense(self.parent(),
                                             self.__f // denom, self._prec)



    def copy(self):
        return PowerSeries_generic_dense(self.parent(), self.__f, self.prec(),
                                         check=False)

    def list(self):
        return self.__f.list()

    def derivative(self):
        return PowerSeries_generic_dense(self.parent(), self.__f.derivative(),
                                         self.prec()-1, check=False)

    def integral(self):
        """
        The integral of this power series with 0 constant term.
        """
        return PowerSeries_generic_dense(self.parent(), self.__f.integral(),
                                         self.prec()+1, check=False)

    def laurent_series(self):
        return self.parent().laurent_series_ring()(self)

    def reversion(self):
        """
        Return the reversion of f, i.e., the series g such that
        g(f(x)) = x.

        EXAMPLES:
            sage: x = PowerSeriesRing(RationalField()).gen()
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
        if self.prec() == infinity:
            raise RuntimeError, "series must have finite precision for reversion."
        f = self._pari_()
        g = f.serreverse()
        return PowerSeries_generic_dense(self.parent(),g.Vecrev(),self.prec())

    def _pari_(self):
        if not isinstance(self.parent().base_ring(),
                          (rational_field.RationalField, integer_ring.IntegerRing)):
            raise NotImplementedError
        if self.prec() == infinity:
            raise RuntimeError, "series precision must be finite for conversion to pari object."
        return pari.pari(str(self))


