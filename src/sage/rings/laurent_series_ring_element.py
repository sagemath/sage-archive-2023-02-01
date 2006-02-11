"""
Laurent Series

AUTHORS:
    -- William Stein: original version
    -- David Joyner: added examples 2006-01-22
"""

import operator

from coerce import bin_op
from infinity import infinity

import laurent_series_ring
import power_series_ring_element
import power_series_ring
import polynomial_element as polynomial
import sage.misc.latex as latex
import sage.rings.ring_element as ring_element
import sage.ext.coerce
from sage.structure.element import Element_cmp_

class LaurentSeries(Element_cmp_, ring_element.RingElement):
    def __init__(self, parent, f, n=0):
        r"""
        Create the Laurent series $t^n \cdot f$.  The default is n=0.

        INPUT:
            parent -- a Laurent series ring
            f -- a power series (or something can be coerced to one)
            n -- integer (default 0)

        OUTPUT:
            a Laurent series

        EXAMPLES:
            sage: K, q = Frac(C[['q']]).objgen()
            sage: K
            Laurent Series Ring in q over Complex Field with 53 bits of precision
            sage: q
            1.0000000000000000*q

        Saving and loading.
            sage: loads(q.dumps()) == q
            True
            sage: loads(K.dumps()) == K
            True
        """
        if not isinstance(parent, laurent_series_ring.LaurentSeriesRing_generic):
            raise TypeError, "parent must be a LaurentSeriesRing"

        ring_element.RingElement.__init__(self, parent)

        if not isinstance(f, power_series_ring_element.PowerSeries):
            if isinstance(f, LaurentSeries):
                n += f.__u
                f = f.unit_part()
            else:
                f = parent.power_series_ring()(f)

        self.__n = n + f.valuation()
        self.__u = f.unit_part()

    def is_zero(self):
        """
        EXAMPLES:
            sage: x = Frac(Q[['x']]).0
            sage: f = 1/x + x + x^2 + 3*x^4 + O(x^7)
            sage: f.is_zero()
            0
            sage: z = 0*f
            sage: z.is_zero()
            1
        """
        return self.__u.is_zero()

    def _im_gens_(self, codomain, im_gens):
        return codomain(self(im_gens[0]))

    def __normalize(self):
        """
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
        self.__u = self.__u.unit_part()

    def __repr__(self):
        if self.is_zero():
            if self.prec() == infinity:
                return "0"
            else:
                return "O(%s^%s)"%(self.parent().variable_name(),self.prec())
        s = " "
        v = self.__u.list()
        valuation = self.__n
        m = len(v)
        X = self.parent().variable_name()
        atomic_repr = self.parent().base_ring().is_atomic_repr()
        first = True
        for n in xrange(m):
            x = v[n]
            e = n + valuation
            if x != 0:
                if not first:
                    s += " + "
                x = str(x)
                if not atomic_repr and n > 0 and (x[1:].find("+") != -1 or x[1:].find("-") != -1):
                    x = "(%s)"%x
                if e == 1:
                    var = "*%s"%X
                elif e == 0:
                    var = ""
                else:
                    var = "*%s^%s"%(X,e)
                s += "%s%s"%(x,var)
                first = False
        if atomic_repr:
            s = s.replace(" + -", " - ")
        s = s.replace(" 1*"," ")
        s = s.replace(" -1*", " -")
        if self.prec() == 0:
            bigoh = "O(1)"
        elif self.prec() == 1:
            bigoh = "O(%s)"%self.parent().variable_name()
        else:
            bigoh = "O(%s^%s)"%(self.parent().variable_name(),self.prec())
        if self.prec() != infinity:
            if s == " ":
                return bigoh
            s += " + %s"%bigoh
        return s[1:]

    def _latex_(self):
        """
        EXAMPLES:
            sage: x = Frac(Q[['x']]).0
            sage: f = (17/2)*x^-2 + x + x^2 + 3*x^4 + O(x^7)
            sage: f._latex_()
            '\\frac{17}{2}x^{-2} + x + x^{2} + 3x^{4} + \\cdots'
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
        X = self.parent().variable_name()
        atomic_repr = self.parent().base_ring().is_atomic_repr()
        first = True
        for n in xrange(m):
            x = v[n]
            e = n + valuation
            if x != 0:
                if not first:
                    s += " + "
                x = latex.latex(x)
                if not atomic_repr and n > 0 and (x[1:].find("+") != -1 or x[1:].find("-") != -1):
                    x = "\\left(%s\\right)"%x
                if e == 1:
                    var = "|%s"%X
                elif e == 0:
                    var = ""
                else:
                    var = "|%s^{%s}"%(X,e)
                s += "%s%s"%(x,var)
                first = False
        if atomic_repr:
            s = s.replace(" + -", " - ")
        s = s.replace(" 1|"," ")
        s = s.replace(" -1|", " -")
        s = s.replace("|","")
        if self.prec() != infinity:
            if s == " ":
                return "0 + \\cdots"
            s += " + \\cdots"
        return s[1:]

    def __getitem__(self, i):
        return self.__u[i-self.__n]

    def __setitem__(self, i, value):
        j = i - self.__n
        if j >= 0:
            self.__u[j] = value
        else: # off to the left
            if value != 0:
                self.__n = self.__n + j
                R = self.parent().base_ring()
                coeffs = [value] + [R(0) for _ in range(1,-j)] + self.__u.list()
                self.__u = self.__u.parent()(coeffs)
        self.__normalize()

    def __add__(self, right):
        """
        EXAMPLES:
            sage: x = Frac(Q[['x']]).0
            sage: f = 1/x^10 + x + x^2 + 3*x^4 + O(x^7)
            sage: g = 1 - x + x^2 - x^4 + O(x^8)
            sage: f*g
            x^-10 - x^-9 + x^-8 - x^-6 + O(x^-2)
        """
        if not isinstance(right, LaurentSeries):
            return bin_op(self, right, operator.add)
        # Algorithm: Make a copy of self then modify coefficients
        # using those of right.
        if right.is_zero():
            return self.add_bigoh(right.prec())
        if self.is_zero():
            return right.add_bigoh(self.prec())
        s = self.copy()
        prec = min(self.prec(),right.prec())
        if prec != infinity:
            m = prec
        else:
            m = max(self.degree(), right.degree())+1
        for i in range(right.valuation(), m):
            s[i] += right[i]
        return s.add_bigoh(prec)

    def add_bigoh(self, prec):
        if prec == infinity or prec >= self.prec():
            return self
        u = self.__u.add_bigoh(prec - self.__n)
        return LaurentSeries(self.parent(), u, self.__n)

    def degree(self):
        """
        Return the degree of a polynomial equivalent to this power
        series modulo big oh of the precision.

        EXAMPLES:
            sage: x = Frac(Q[['x']]).0
            sage: g = x^2 - x^4 + O(x^8)
            sage: g.degree()
            4
            sage: g = -10/x^5 + x^2 - x^4 + O(x^8)
            sage: g.degree()
            4
        """
        return self.__u.degree() + self.valuation()

    def __sub__(self, right):
        return self + (-1)*right

    def __neg__(self):
        return (-1)*self

    def __mul__(self, right):
        """
        EXAMPLES:
            sage: x = Frac(Q[['x']]).0
            sage: f = 1/x^3 + x + x^2 + 3*x^4 + O(x^7)
            sage: g = 1 - x + x^2 - x^4 + O(x^8)
            sage: f*g
            x^-3 - x^-2 + x^-1 + 4*x^4 + O(x^5)
        """
        if not isinstance(right, LaurentSeries):
            return bin_op(self, right, operator.mul)
        return LaurentSeries(self.parent(),
                             self.__u * right.__u,
                             self.__n + right.__n)
    def __pow__(self, right):
        """
        EXAMPLES:
            sage: x = Frac(Q[['x']]).0
            sage: f = x + x^2 + 3*x^4 + O(x^7)
            sage: g = 1/x^10 - x + x^2 - x^4 + O(x^8)
            sage: f^7
            x^7 + 7*x^8 + 21*x^9 + 56*x^10 + 161*x^11 + 336*x^12 + O(x^13)
            sage: g^7
            x^-70 - 7*x^-59 + 7*x^-58 - 7*x^-56 + O(x^-52)
        """
        right=int(right)
        return LaurentSeries(self.parent(), self.__u**right, self.__n*right)


    def __div__(self, right):
        """
        EXAMPLES:
            sage: x = Frac(Q[['x']]).0
            sage: f = x + x^2 + 3*x^4 + O(x^7)
            sage: g = 1/x^7 - x + x^2 - x^4 + O(x^8)
            sage: f/x
            1 + x + 3*x^3 + O(x^6)
            sage: f/g
            x^8 + x^9 + 3*x^11 + O(x^14)
        """
        if not isinstance(right, LaurentSeries):
            return bin_op(self, right, operator.div)
        return LaurentSeries(self.parent(),
                             self.__u / right.__u,
                             self.__n - right.__n)


    def _cmp_(self, right):
        """
        EXAMPLES:
            sage: x = Frac(Q[['x']]).0
            sage: f = x + x^2 + 3*x^4 + O(x^7)
            sage: g = 1/x^7 - x + x^2 - x^4 + O(x^8)
            sage: f<g
            False
            sage: f>g
            True
        """
        if self.__n < right.__n:
            return -1
        elif self.__n > right.__n:
            return 1
        return cmp(self.__u, right.__u)

    def unit_part(self):
        """
        EXAMPLES:
            sage: x = Frac(Q[['x']]).0
            sage: f = x + x^2 + 3*x^4 + O(x^7)
            sage: f/x
            1 + x + 3*x^3 + O(x^6)
            sage: f.unit_part()
            1 + x + 3*x^3 + O(x^6)
            sage: g = 1/x^7 - x + x^2 - x^4 + O(x^8)
            sage: g.unit_part()
            1 - x^8 + x^9 - x^11 + O(x^15)
        """
        return self.__u

    def valuation(self):
        """
        EXAMPLES:
            sage: x = Frac(Q[['x']]).0
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
            sage: x = Frac(Q[['x']]).0
            sage: f = 1/x + x^2 + 3*x^4 + O(x^7)
            sage: f.variable()
            'x'
        """
        return self.parent().variable_name()

    def prec(self):
        """
        This function returns the n so that the Laurent series is
        of the form (stuff) + $O(t^n)$.  It doesn't matter how many
        negative powers appear in the expansion.  In particular,
        prec could be negative.

        EXAMPLES:
            sage: x = Frac(Q[['x']]).0
            sage: f = x^2 + 3*x^4 + O(x^7)
            sage: f.prec()
            7
            sage: g = 1/x^10 - x + x^2 - x^4 + O(x^8)
            sage: g.prec()
            8
        """
        return self.__u.prec() + self.valuation()

    def copy(self):
        return LaurentSeries(self.parent(), self.__u.copy(), self.__n)

    def derivative(self):
        """
        The formal derivative of this Laurent series.

        EXAMPLES:
            sage: x = Frac(Q[['x']]).0
            sage: f = x^2 + 3*x^4 + O(x^7)
            sage: f.derivative()
            2*x + 12*x^3 + O(x^6)
            sage: g = 1/x^10 - x + x^2 - x^4 + O(x^8)
            sage: g.derivative()
            -10*x^-11 - 1 + 2*x - 4*x^3 + O(x^7)
        """
        if self.is_zero():
            return LaurentSeries(self.parent(), 0, self.__u.prec() - 1)
        n = self.__n
        a = self.__u.list()
        v = [(n+m)*a[m] for m in range(len(a))]
        u = self.parent().power_series_ring()(v, self.__u.prec())
        return LaurentSeries(self.parent(), u, n-1)

    def integral(self):
        r"""
        The formal integral of this Laurent series with 0 constant term.

        EXAMPLES:
        The integral may or may not be defined if the base ring
        is not a field.
            sage: t = LaurentSeriesRing(Z, 't').0
            sage: f = 2*t^-3 + 3*t^2 + O(t^4)
            sage: f.integral()
            -t^-2 + t^3 + O(t^5)

            sage: f = t^3
            sage: f.integral()
            Traceback (most recent call last):
            ...
            ArithmeticError: Coefficients of integral of t^3 cannot be coerced into the base ring


        The integral of 1/t is $\log(t)$, which is not given by a Laurent series:

            sage: t = Frac(Q[['t']]).gen()
            sage: f = -1/t^3 - 31/t + O(t^3)
            sage: f.integral()
            Traceback (most recent call last):
            ...
            ArithmeticError: The integral of -t^-3 - 31*t^-1 + O(t^3) is not a Laurent series, since t^-1 has nonzero coefficient -31.
        """
        n = self.__n
        a = self.__u.list()
        if self[-1] != 0:
            raise ArithmeticError, \
                  "The integral of %s is not a Laurent series, since t^-1 has nonzero coefficient %s."%(self,self[-1])

        if n < 0:
            v = [a[i]/(n+i+1) for i in range(-1-n)] + [0]
        else:
            v = []
        v += [a[i]/(n+i+1) for i in range(max(-n,0), len(a))]
        try:
            u = self.parent().power_series_ring()(v, self.__u.prec())
        except TypeError:
            raise ArithmeticError, "Coefficients of integral of %s cannot be coerced into the base ring"%self
        return LaurentSeries(self.parent(), u, n+1)


    def power_series(self):
        if self.__n < 0:
            raise ArithmeticError, "self (=%s) is a not a power series"%self
        u = self.__u
        t = u.parent().gen()
        return t**(self.__n) * u

    def __call__(self, x):
        """
        Compute value of this Laurent series at x.

        EXAMPLES:
            sage: t = LaurentSeriesRing(Z, 't').0
            sage: f = t^(-2) + t^2 + O(t^8)
            sage: f(2)
            17/4
            sage: f(-1)
            2
            sage: f(1/3)
            82/9
        """
        if x.parent() == self.parent():
            x = x.add_bigoh(self.prec()*x.valuation())
        return self.__u(x) * (x**self.__n)


