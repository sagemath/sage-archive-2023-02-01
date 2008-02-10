"""
Laurent Series

EXAMPLES:

    sage: R.<t> = LaurentSeriesRing(GF(7), 't'); R
    Laurent Series Ring in t over Finite Field of size 7
    sage: f = 1/(1-t+O(t^10)); f
    1 + t + t^2 + t^3 + t^4 + t^5 + t^6 + t^7 + t^8 + t^9 + O(t^10)

Laurent series are immutable:
    sage: f[2]
    1
    sage: f[2] = 5
    Traceback (most recent call last):
    ...
    IndexError: Laurent series are immutable

We compute with a Laurent series over the complex mpfr numbers.
    sage: K.<q> = Frac(CC[['q']])
    sage: K
    Laurent Series Ring in q over Complex Field with 53 bits of precision
    sage: q
    1.00000000000000*q

Saving and loading.
    sage: loads(q.dumps()) == q
    True
    sage: loads(K.dumps()) == K
    True

IMPLEMENTATION: Laurent series in SAGE are represented internally as a
    power of the variable times the unit part (which need not be a
    unit -- it's a polynomial with nonzero constant term).  The zero
    Laurent series has unit part 0.


AUTHORS:
    -- William Stein: original version
    -- David Joyner: added examples 2006-01-22
    -- Robert Bradshaw: optimizations, shifting 2007-04
    -- Robert Bradshaw: SageX version
"""

import operator

from infinity import infinity

import laurent_series_ring
import power_series_ring_element
import power_series_ring
import sage.rings.polynomial.polynomial_element as polynomial
import sage.misc.latex
import sage.rings.ring_element as ring_element
from sage.rings.integer import Integer

from sage.structure.element cimport Element, ModuleElement, RingElement, AlgebraElement

include "../ext/stdsage.pxi"

def is_LaurentSeries(x):
    return isinstance(x, LaurentSeries)


cdef class LaurentSeries(AlgebraElement):
    """
    A Laurent Series.
    """

    def __init__(self, parent, f, n=0):
        r"""
        Create the Laurent series $t^n \cdot f$.  The default is n=0.

        INPUT:
            parent -- a Laurent series ring
            f -- a power series (or something can be coerced to one); note that
                 f does *not* have to be a unit.
            n -- integer (default 0)

        OUTPUT:
            a Laurent series

        EXAMPLES:
            sage: R.<q> = LaurentSeriesRing(ZZ)
            sage: R([1,2,3])
            1 + 2*q + 3*q^2
            sage: R([1,2,3],-5)
            q^-5 + 2*q^-4 + 3*q^-3

            sage: S.<s> = LaurentSeriesRing(GF(5))
            sage: T.<t> = PowerSeriesRing(pAdicRing(5))
            sage: S(t)
            s
            sage: parent(S(t))
            Laurent Series Ring in s over Finite Field of size 5
            sage: parent(S(t)[1])
            Finite Field of size 5
        """
        AlgebraElement.__init__(self, parent)

        if PY_TYPE_CHECK(f, LaurentSeries):
            n += (<LaurentSeries>f).__n
            if (<LaurentSeries>f).__u._parent is parent.power_series_ring():
                f = (<LaurentSeries>f).__u
            else:
                f = parent.power_series_ring()((<LaurentSeries>f).__u)
        elif not PY_TYPE_CHECK(f, PowerSeries):
            f = parent.power_series_ring()(f)
        ## now this is a power series, over a different ring ...
        ## requires that power series rings with same vars over the
        ## same parent are unique.
        elif parent is not f.parent():
            f = parent.power_series_ring()(f)


        # self is that t^n * u:
        cdef long val
        if not f:
            if n == infinity:
                self.__n = 0
                self.__u = parent.power_series_ring()(0)
            else:
                self.__n = n
                self.__u = f
        else:
            val = f.valuation()
            if val == 0:
                self.__n = n    # power of the variable
                self.__u = f    # unit part
            else:
                self.__n = n + val
                self.__u = f >> val

    def __reduce__(self):
        return make_element_from_parent, (self._parent, self.__u, self.__n)

    def change_ring(self, R):
        return self.parent().change_ring(R)(self)

    def is_unit(self):
        """
        Returns True if this is Laurent series is a unit in this ring.

        EXAMPLES:
            sage: R.<t> = LaurentSeriesRing(QQ)
            sage: (2+t).is_unit()
            True
            sage: f = 2+t^2+O(t^10); f.is_unit()
            True
            sage: 1/f
            1/2 - 1/4*t^2 + 1/8*t^4 - 1/16*t^6 + 1/32*t^8 + O(t^10)
            sage: R(0).is_unit()
            False
            sage: R.<s> = LaurentSeriesRing(ZZ)
            sage: f = 2 + s^2 + O(s^10)
            sage: f.is_unit()
            False
            sage: 1/f
            Traceback (most recent call last):
            ...
            ArithmeticError: division not defined

        ALGORITHM: A Laurent series is a unit if and only if
        its "unit part" is a unit.
        """
        return self.__u.is_unit()

    def is_zero(self):
        """
        EXAMPLES:
            sage: x = Frac(QQ[['x']]).0
            sage: f = 1/x + x + x^2 + 3*x^4 + O(x^7)
            sage: f.is_zero()
            0
            sage: z = 0*f
            sage: z.is_zero()
            1
        """
        return self.__u.is_zero()

    def __nonzero__(self):
        return not not self.__u

    def _im_gens_(self, codomain, im_gens):
        return codomain(self(im_gens[0]))

    def __normalize(self):
        r"""
        A Laurent series is a pair (u(t), n), where either u=0 (to
        some precision) or u is a unit.  This pair corresponds to
        $t^n\cdot u(t)$.
        """
        if self.is_zero():
            return
        v = self.__u.valuation()
        if v == 0:
            return
        self.__n += v
        self.__u = self.__u.valuation_zero_part()

    def _repr_(self):
        """
        EXAMPLES:
            sage: R.<t> = LaurentSeriesRing(QQ)
            sage: (2 + (2/3)*t^3).__repr__()
            '2 + 2/3*t^3'
        """
        if self.is_zero():
            if self.prec() == infinity:
                return "0"
            else:
                return "O(%s^%s)"%(self._parent.variable_name(),self.prec())
        s = " "
        v = self.__u.list()
        valuation = self.__n
        m = len(v)
        X = self._parent.variable_name()
        atomic_repr = self._parent.base_ring().is_atomic_repr()
        first = True
        for n in xrange(m):
            x = v[n]
            e = n + valuation
            x = str(x)
            if x != '0':
                if not first:
                    s += " + "
                if not atomic_repr and (x[1:].find("+") != -1 or x[1:].find("-") != -1):
                    x = "(%s)"%x
                if e == 1:
                    var = "*%s"%X
                elif e == 0:
                    var = ""
                else:
                    var = "*%s^%s"%(X,e)
                s += "%s%s"%(x,var)
                first = False
        s = s.replace(" + -", " - ")
        s = s.replace(" 1*"," ")
        s = s.replace(" -1*", " -")
        if self.prec() == 0:
            bigoh = "O(1)"
        elif self.prec() == 1:
            bigoh = "O(%s)"%self._parent.variable_name()
        else:
            bigoh = "O(%s^%s)"%(self._parent.variable_name(),self.prec())
        if self.prec() != infinity:
            if s == " ":
                return bigoh
            s += " + %s"%bigoh
        return s[1:]

    def _latex_(self):
        r"""
        EXAMPLES:
            sage: x = Frac(QQ[['x']]).0
            sage: f = (17/2)*x^-2 + x + x^2 + 3*x^4 + O(x^7)
            sage: latex(f)
            \frac{\frac{17}{2}}{x^{2}} + x + x^{2} + 3x^{4} + O(x^{7})
        """
        if self.is_zero():
            if self.prec() == infinity:
                return "0"
            else:
                return "0 + \\cdots"
        s = " "
        v = self.__u.list()
        valuation = self.__n
        m = len(v)
        X = self._parent.latex_variable_names()[0]
        atomic_repr = self._parent.base_ring().is_atomic_repr()
        first = True
        for n in xrange(m):
            x = v[n]
            e = n + valuation
            x = sage.misc.latex.latex(x)
            if x != '0':
                if not first:
                    s += " + "
                if not atomic_repr and n > 0 and (x[1:].find("+") != -1 or x[1:].find("-") != -1):
                    x = "\\left(%s\\right)"%x
                if e == 1:
                    var = "|%s"%X
                elif e == 0:
                    var = ""
                elif e > 0:
                    var = "|%s^{%s}"%(X,e)
                if e >= 0:
                    s += "%s%s"%(x,var)
                else: # negative e
                    if e == -1:
                        s += "\\frac{%s}{%s}"%(x, X)
                    else:
                        s += "\\frac{%s}{%s^{%s}}"%(x, X,-e)
                first = False
        s = s.replace(" + -", " - ")
        s = s.replace(" 1|"," ")
        s = s.replace(" -1|", " -")
        s = s.replace("|","")
        pr = self.prec()
        if pr != infinity:
            if pr == 0:
                bigoh = "O(1)"
            elif pr == 1:
                bigoh = "O(%s)"%(X,)
            else:
                bigoh = "O(%s^{%s})"%(X,pr)
            if s == " ":
                return bigoh
            s += " + %s"%bigoh
        return s[1:]

    def __hash__(self):
        return hash(self.__u) ^ self.__n

    def __getitem__(self, i):
        """
        EXAMPLES:
            sage: R.<t> = LaurentSeriesRing(QQ)
            sage: f = -5/t^(10) + t + t^2 - 10/3*t^3; f
            -5*t^-10 + t + t^2 - 10/3*t^3
            sage: f[-10]
            -5
            sage: f[1]
            1
            sage: f[3]
            -10/3
            sage: f[-9]
            0
        """
        return self.__u[i-self.__n]

    def __getslice__(self, i, j):
        """
        EXAMPLES:
            sage: R.<t> = LaurentSeriesRing(QQ)
            sage: f = -5/t^(10) + 1/3 + t + t^2 - 10/3*t^3 + O(t^5); f
            -5*t^-10 + 1/3 + t + t^2 - 10/3*t^3 + O(t^5)
            sage: f[-10:2]
            -5*t^-10 + 1/3 + t + O(t^5)
            sage: f[0:]
            1/3 + t + t^2 - 10/3*t^3 + O(t^5)
        """
        if j > self.__u.degree():
            j = self.__u.degree()
        f = self.__u[i-self.__n:j-self.__n]
        return LaurentSeries(self._parent, f, self.__n)

    def __iter__(self):
        """
        Iterate through the coefficients from the first nonzero one to
        the last nonzero one.

        EXAMPLES:
            sage: R.<t> = LaurentSeriesRing(QQ)
            sage: f = -5/t^(2) + t + t^2 - 10/3*t^3; f
            -5*t^-2 + t + t^2 - 10/3*t^3
            sage: for a in f: print a
            -5
            0
            0
            1
            1
            -10/3
        """
        return iter(self.__u)

    def list(self):
        return self.__u.list()

    def __setitem__(self, n, value):
        """
        EXAMPLES:
            sage: R.<t> = LaurentSeriesRing(QQ)
            sage: f = t^2 + t^3 + O(t^10)
            sage: f[2] = 5
            Traceback (most recent call last):
            ...
            IndexError: Laurent series are immutable
        """
        raise IndexError, "Laurent series are immutable"

    def _unsafe_mutate(self, i, value):
        """
        SAGE assumes throughout that commutative ring elements are immutable.
        This is relevant for caching, etc.  But sometimes you need to change
        a Laurent series and you really know what you're doing.  That's
        when this function is for you.

        EXAMPLES:

        """
        j = i - self.__n
        if j >= 0:
            self.__u._unsafe_mutate(j, value)
        else: # off to the left
            if value != 0:
                self.__n = self.__n + j
                R = self._parent.base_ring()
                coeffs = [value] + [R(0) for _ in range(1,-j)] + self.__u.list()
                self.__u = self.__u._parent(coeffs)
        self.__normalize()

    cdef ModuleElement _add_c_impl(self, ModuleElement right_m):
        """
        Add two power series with the same parent.

        EXAMPLES:
            sage: R.<t> = LaurentSeriesRing(QQ)
            sage: t + t
            2*t
            sage: f = 1/t + t^2 + t^3 - 17/3 * t^4 + O(t^5)
            sage: g = 1/(1-t + O(t^7)); g
            1 + t + t^2 + t^3 + t^4 + t^5 + t^6 + O(t^7)
            sage: f + g
            t^-1 + 1 + t + 2*t^2 + 2*t^3 - 14/3*t^4 + O(t^5)
            sage: f + 0
            t^-1 + t^2 + t^3 - 17/3*t^4 + O(t^5)
            sage: 0 + f
            t^-1 + t^2 + t^3 - 17/3*t^4 + O(t^5)
            sage: R(0) + R(0)
            0
            sage: (t^3 + O(t^10)) + (t^-3 +O(t^9))
            t^-3 + t^3 + O(t^9)

        ALGORITHM:
            Shift the unit parts to align them, then add.
        """
        cdef LaurentSeries right = <LaurentSeries>right_m
        cdef long m

        # 1. Special case when one or the other is 0.
        if not right:
            return self.add_bigoh(right.prec())
        if not self:
            return right.add_bigoh(self.prec())

        # 2. Align the unit parts.
        if self.__n < right.__n:
            m = self.__n
            f1 = self.__u
            f2 = right.__u << right.__n - m
        elif self.__n > right.__n:
            m = right.__n
            f1 = self.__u << self.__n - m
            f2 = right.__u
        else:
            m = self.__n
            f1 = self.__u
            f2 = right.__u
        # 3. Add
        return LaurentSeries(self._parent, f1 + f2, m)

    cdef ModuleElement _iadd_c_impl(self, ModuleElement right_m):
        """
        EXAMPLES:
            sage: R.<t> = LaurentSeriesRing(QQ)
            sage: f = t+t
            sage: f += t; f
            3*t
            sage: f += O(t^5); f
            3*t + O(t^5)
        """
        cdef LaurentSeries right = <LaurentSeries>right_m
        if self.__n == right.__n:
            self.__u += right.__u
            return self
        else:
            return self._add_c_impl(right)

    cdef ModuleElement _sub_c_impl(self, ModuleElement right_m):
        """
        Subtract two power series with the same parent.

        EXAMPLES:
            sage: R.<t> = LaurentSeriesRing(QQ)
            sage: t - t
            0
            sage: t^5 + 2 * t^-5
            2*t^-5 + t^5

        ALGORITHM:
            Shift the unit parts to align them, then subtract.
        """
        cdef LaurentSeries right = <LaurentSeries>right_m
        cdef long m

        # 1. Special case when one or the other is 0.
        if not right:
            return self.add_bigoh(right.prec())
        if not self:
            return -right.add_bigoh(self.prec())

        # 2. Align the unit parts.
        if self.__n < right.__n:
            m = self.__n
            f1 = self.__u
            f2 = right.__u << right.__n - m
        else:
            m = right.__n
            f1 = self.__u << self.__n - m
            f2 = right.__u
        # 3. Subtract
        return LaurentSeries(self._parent, f1 - f2, m)


    def add_bigoh(self, prec):
        """
        EXAMPLES:
            sage: R.<t> = LaurentSeriesRing(QQ)
            sage: f = t^2 + t^3 + O(t^10); f
            t^2 + t^3 + O(t^10)
            sage: f.add_bigoh(5)
            t^2 + t^3 + O(t^5)
        """
        if prec == infinity or prec >= self.prec():
            return self
        u = self.__u.add_bigoh(prec - self.__n)
        return LaurentSeries(self._parent, u, self.__n)

    def degree(self):
        """
        Return the degree of a polynomial equivalent to this power
        series modulo big oh of the precision.

        EXAMPLES:
            sage: x = Frac(QQ[['x']]).0
            sage: g = x^2 - x^4 + O(x^8)
            sage: g.degree()
            4
            sage: g = -10/x^5 + x^2 - x^4 + O(x^8)
            sage: g.degree()
            4
        """
        return self.__u.degree() + self.__n


    def __neg__(self):
        """
        sage: R.<t> = LaurentSeriesRing(QQ)
        sage: -(1+t^5)
        -1 - t^5
        sage: -(1/(1+t+O(t^5)))
        -1 + t - t^2 + t^3 - t^4 + O(t^5)
        """
        return LaurentSeries(self._parent, -self.__u, self.__n)

    cdef RingElement _mul_c_impl(self, RingElement right_r):
        """
        EXAMPLES:
            sage: x = Frac(QQ[['x']]).0
            sage: f = 1/x^3 + x + x^2 + 3*x^4 + O(x^7)
            sage: g = 1 - x + x^2 - x^4 + O(x^8)
            sage: f*g
            x^-3 - x^-2 + x^-1 + 4*x^4 + O(x^5)
        """
        cdef LaurentSeries right = <LaurentSeries>right_r
        return LaurentSeries(self._parent,
                             self.__u * right.__u,
                             self.__n + right.__n)

    cdef RingElement _imul_c_impl(self, RingElement right_r):
        """
        EXAMPLES:
            sage: x = Frac(QQ[['x']]).0
            sage: f = 1/x^3 + x + x^2 + 3*x^4 + O(x^7)
            sage: g = 1 - x + x^2 - x^4 + O(x^8)
            sage: f *= g; f
            x^-3 - x^-2 + x^-1 + 4*x^4 + O(x^5)
        """
        cdef LaurentSeries right = <LaurentSeries>right_r
        self.__u *= right.__u
        self.__n += right.__n
        return self

    cdef ModuleElement _rmul_c_impl(self, RingElement c):
        return LaurentSeries(self._parent, self.__u._rmul_c(c), self.__n)

    cdef ModuleElement _lmul_c_impl(self, RingElement c):
        return LaurentSeries(self._parent, self.__u._lmul_c(c), self.__n)

    cdef ModuleElement _ilmul_c_impl(self, RingElement c):
        self.__u *= c
        return self

    def __pow__(_self, r, dummy):
        """
        EXAMPLES:
            sage: x = Frac(QQ[['x']]).0
            sage: f = x + x^2 + 3*x^4 + O(x^7)
            sage: g = 1/x^10 - x + x^2 - x^4 + O(x^8)
            sage: f^7
            x^7 + 7*x^8 + 21*x^9 + 56*x^10 + 161*x^11 + 336*x^12 + O(x^13)
            sage: g^7
            x^-70 - 7*x^-59 + 7*x^-58 - 7*x^-56 + O(x^-52)
        """
        cdef LaurentSeries self = _self
        right=int(r)
        if right != r:
            raise ValueError, "exponent must be an integer"
        return LaurentSeries(self._parent, self.__u**right, self.__n*right)

    def shift(self, k):
        r"""
        Returns this laurent series multiplied by the power $t^n$. Does not
        change this series.

        NOTE:
            Despite the fact that higher order terms are printed to the
            right in a power series, right shifting decreases the powers
            of $t$, while left shifting increases them. This is to be
            consistant with polynomials, integers, etc.

        EXAMPLES:
            sage: R.<t> = LaurentSeriesRing(QQ['y'])
            sage: f = (t+t^-1)^4; f
            t^-4 + 4*t^-2 + 6 + 4*t^2 + t^4
            sage: f.shift(10)
            t^6 + 4*t^8 + 6*t^10 + 4*t^12 + t^14
            sage: f >> 10
            t^-14 + 4*t^-12 + 6*t^-10 + 4*t^-8 + t^-6
            sage: t << 4
            t^5
            sage: t + O(t^3) >> 4
            t^-3 + O(t^-1)

        AUTHOR:
            -- Robert Bradshaw (2007-04-18)
        """
        return LaurentSeries(self._parent, self.__u, self.__n + k)

    def __lshift__(LaurentSeries self, k):
        return LaurentSeries(self._parent, self.__u, self.__n + k)

    def __rshift__(LaurentSeries self, k):
        return LaurentSeries(self._parent, self.__u, self.__n - k)

    def truncate(self, long n):
        r"""
        Returns the laurent series of degree $ < n$ which is equivalent to self
        modulo $x^n$.
        """
        if n <= self.__n:
            return LaurentSeries(self._parent, 0)
        else:
            return LaurentSeries(self._parent, self.__u.truncate_powerseries(n - self.__n), self.__n)

    def truncate_neg(self, long n):
        r"""
        Returns the laurent series equivalent to self except without any degree < n terms.

        This is equivalent to $\code{self - self.truncate(n)}$.
        """
        return LaurentSeries(self._parent, self.__u >> (n - self.__n), n)

    cdef RingElement _div_c_impl(self, RingElement right_r):
        """
        EXAMPLES:
            sage: x = Frac(QQ[['x']]).0
            sage: f = x + x^2 + 3*x^4 + O(x^7)
            sage: g = 1/x^7 - x + x^2 - x^4 + O(x^8)
            sage: f/x
            1 + x + 3*x^3 + O(x^6)
            sage: f/g
            x^8 + x^9 + 3*x^11 + O(x^14)
        """
        cdef LaurentSeries right = <LaurentSeries>right_r
        if right.__u.is_zero():
            raise ZeroDivisionError
        try:
            return LaurentSeries(self._parent,
                             self.__u / right.__u,
                             self.__n - right.__n)
        except TypeError, msg:
            # todo: this could also make something in the formal fraction field.
            raise ArithmeticError, "division not defined"


    def common_prec(self, f):
        r"""
        Returns minimum precision of $f$ and self.

        EXAMPLES:
            sage: R.<t> = LaurentSeriesRing(QQ)

            sage: f = t^(-1) + t + t^2 + O(t^3)
            sage: g = t + t^3 + t^4 + O(t^4)
            sage: f.common_prec(g)
            3
            sage: g.common_prec(f)
            3

            sage: f = t + t^2 + O(t^3)
            sage: g = t^(-3) + t^2
            sage: f.common_prec(g)
            3
            sage: g.common_prec(f)
            3

            sage: f = t + t^2
            sage: f = t^2
            sage: f.common_prec(g)
            +Infinity

            sage: f = t^(-3) + O(t^(-2))
            sage: g = t^(-5) + O(t^(-1))
            sage: f.common_prec(g)
            -2
        """
        if self.prec() is infinity:
            return f.prec()
        elif f.prec() is infinity:
            return self.prec()
        return min(self.prec(), f.prec())

    def __richcmp__(left, right, int op):
        return (<Element>left)._richcmp(right, op)

    cdef int _cmp_c_impl(self, Element right_r) except -2:
        r"""
        Comparison of self and right.

        We say two approximate laurent series are equal, if they agree
        for all coefficients up to the *minimum* of the precisions of
        each. Comparison is done in dictionary order from lowest degree to
        highest degree coefficients (this is different than
        polynomials).

        See power_series_ring_element.__cmp__() for more information.

        EXAMPLES:
            sage: R.<x> = LaurentSeriesRing(QQ)
            sage: f = x^(-1) + 1 + x + O(x^2)
            sage: g = x^(-1) + 1 + O(x)
            sage: f == g
            True

            sage: f = x^(-1) + 1 + x + O(x^2)
            sage: g = x^(-1) + 2 + O(x)
            sage: f == g
            False
            sage: f < g
            True
            sage: f > g
            False

            sage: f = x^(-2) + 1 + x + O(x^2)
            sage: g = x^(-1) + 2 + O(x)
            sage: f == g
            False
            sage: f < g
            False
            sage: f > g
            True
        """
        cdef LaurentSeries right = <LaurentSeries>right_r

        prec = self.common_prec(right)
        if not prec:
            return 0
        zero = self.base_ring()(0)

        if not self and not right:
            if self.__n < right.__n:
                return cmp(self.__u[0], zero)
            elif self.__n > right.__n:
                return cmp(zero, right.__u[0])

        # zero pad coefficients on the left, to line them up for comparison
        cdef long n = min(self.__n, right.__n)
        x = [zero] * (self.__n - n) + self.__u.list()
        y = [zero] * (right.__n - n) + right.__u.list()

        # zero pad on right to make the lists the same length
        # (this is necessary since the power series list() function just
        # returns the coefficients of the underlying polynomial, which may
        # have zeroes in the high coefficients)
        if len(x) < len(y):
            x.extend([zero] * (len(y) - len(x)))
        elif len(y) < len(x):
            y.extend([zero] * (len(x) - len(y)))

        if not (prec is infinity):
            x = x[:(prec-n)]
            y = y[:(prec-n)]

        return cmp(x,y)

    def valuation_zero_part(self):
        """
        EXAMPLES:
            sage: x = Frac(QQ[['x']]).0
            sage: f = x + x^2 + 3*x^4 + O(x^7)
            sage: f/x
            1 + x + 3*x^3 + O(x^6)
            sage: f.valuation_zero_part()
            1 + x + 3*x^3 + O(x^6)
            sage: g = 1/x^7 - x + x^2 - x^4 + O(x^8)
            sage: g.valuation_zero_part()
            1 - x^8 + x^9 - x^11 + O(x^15)
        """
        return self.__u

    def valuation(self):
        """
        EXAMPLES:
            sage: x = Frac(QQ[['x']]).0
            sage: f = 1/x + x^2 + 3*x^4 + O(x^7)
            sage: g = 1 - x + x^2 - x^4 + O(x^8)
            sage: f.valuation()
            -1
            sage: g.valuation()
            0
        """
        return self.__n

    def variable(self):
        """
        EXAMPLES:
            sage: x = Frac(QQ[['x']]).0
            sage: f = 1/x + x^2 + 3*x^4 + O(x^7)
            sage: f.variable()
            'x'
        """
        return self._parent.variable_name()

    def prec(self):
        """
        This function returns the n so that the Laurent series is
        of the form (stuff) + $O(t^n)$.  It doesn't matter how many
        negative powers appear in the expansion.  In particular,
        prec could be negative.

        EXAMPLES:
            sage: x = Frac(QQ[['x']]).0
            sage: f = x^2 + 3*x^4 + O(x^7)
            sage: f.prec()
            7
            sage: g = 1/x^10 - x + x^2 - x^4 + O(x^8)
            sage: g.prec()
            8
        """
        return self.__u.prec() + self.__n

    def copy(self):
        return LaurentSeries(self._parent, self.__u.copy(), self.__n)

    def derivative(self):
        """
        The formal derivative of this Laurent series.

        EXAMPLES:
            sage: x = Frac(QQ[['x']]).0
            sage: f = x^2 + 3*x^4 + O(x^7)
            sage: f.derivative()
            2*x + 12*x^3 + O(x^6)
            sage: g = 1/x^10 - x + x^2 - x^4 + O(x^8)
            sage: g.derivative()
            -10*x^-11 - 1 + 2*x - 4*x^3 + O(x^7)
        """
        if self.is_zero():
            return LaurentSeries(self._parent, 0, self.__u.prec() - 1)
        cdef long m, n = self.__n
        a = self.__u.list()
        v = [(n+m)*a[m] for m from 0 <= m < len(a)]
        u = self._parent.power_series_ring()(v, self.__u.prec())
        return LaurentSeries(self._parent, u, n-1)

    def integral(self):
        r"""
        The formal integral of this Laurent series with 0 constant term.

        EXAMPLES:
        The integral may or may not be defined if the base ring
        is not a field.
            sage: t = LaurentSeriesRing(ZZ, 't').0
            sage: f = 2*t^-3 + 3*t^2 + O(t^4)
            sage: f.integral()
            -t^-2 + t^3 + O(t^5)

            sage: f = t^3
            sage: f.integral()
            Traceback (most recent call last):
            ...
            ArithmeticError: Coefficients of integral cannot be coerced into the base ring


        The integral of 1/t is $\log(t)$, which is not given by a Laurent series:

            sage: t = Frac(QQ[['t']]).0
            sage: f = -1/t^3 - 31/t + O(t^3)
            sage: f.integral()
            Traceback (most recent call last):
            ...
            ArithmeticError: The integral of is not a Laurent series, since t^-1 has nonzero coefficient.

        Another example with just one negative coefficient:
            sage: A.<t> = QQ[[]]
	    sage: f = -2*t^(-4) + O(t^8)
	    sage: f.integral()
	    2/3*t^-3 + O(t^9)
	    sage: f.integral().derivative() == f
	    True

        """
        cdef long i, n = self.__n
        a = self.__u.list()
        if self[-1] != 0:
            raise ArithmeticError, \
                  "The integral of is not a Laurent series, since t^-1 has nonzero coefficient."

        if n < 0:
            v = [a[i]/(n+i+1) for i in range(min(-1-n,len(a)))] + [0]
        else:
            v = []
        v += [a[i]/(n+i+1) for i in range(max(-n,0), len(a))]
        try:
            u = self._parent.power_series_ring()(v, self.__u.prec())
        except TypeError:
            raise ArithmeticError, "Coefficients of integral cannot be coerced into the base ring"
        return LaurentSeries(self._parent, u, n+1)


    def power_series(self):
        """
        EXAMPLES:
            sage: R.<t> = LaurentSeriesRing(ZZ)
            sage: f = 1/(1-t+O(t^10)); f.parent()
            Laurent Series Ring in t over Integer Ring
            sage: g = f.power_series(); g
            1 + t + t^2 + t^3 + t^4 + t^5 + t^6 + t^7 + t^8 + t^9 + O(t^10)
            sage: parent(g)
            Power Series Ring in t over Integer Ring
            sage: f = 3/t^2 +  t^2 + t^3 + O(t^10)
            sage: f.power_series()
            Traceback (most recent call last):
            ...
            ArithmeticError: self is a not a power series
        """
        if self.__n < 0:
            raise ArithmeticError, "self is a not a power series"
        u = self.__u
        t = u.parent().gen()
        return t**(self.__n) * u

    def __call__(self, *x):
        """
        Compute value of this Laurent series at x.

        EXAMPLES:
            sage: R.<t> = LaurentSeriesRing(ZZ)
            sage: f = t^(-2) + t^2 + O(t^8)
            sage: f(2)
            17/4
            sage: f(-1)
            2
            sage: f(1/3)
            82/9
        """
        if isinstance(x[0], tuple):
            x = x[0]
        return self.__u(x) * (x[0]**self.__n)


def make_element_from_parent(parent, *args):
    return parent(*args)
