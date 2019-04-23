r"""
Puiseux Series Ring Element
===========================

Behavior of the elements of the :class:`PuiseuxSeriesRing`.

This defines how to construct an element given a Laurent series,
center, and ramification index. The class is designed to use and be
compatible with Sage's coercion model.

A Puiseux series is a series of the form .. MATH::

    p(x) = \sum_{n=N}^\infty a_n (x-a)^{n/e}

where the integer :math:`e` is called the *ramification index* of the series
and the number :math:`a` is the *center*. A Puiseux series is essentially a
Laurent series but with fractional exponents.

EXAMPLES:

We begin by constructing the ring of Puiseux series in `x` with coefficients
in the rationals. ::

    sage: R.<x> = PuiseuxSeriesRing(QQ)

This command also defines `x` as the generator of this ring.

When constructing a Puiseux series, the ramification index is automatically
determined from the greatest common divisor of the exponents. ::

    sage: p = x^(1/2); p
    x^(1/2)
    sage: p.ramification_index()
    2
    sage: q = x^(1/2) + x**(1/3); q
    x^(1/3) + x^(1/2)
    sage: q.ramification_index()
    6

Other arithmetic can be performed with Puiseux Series. ::

    sage: p + q
    x^(1/3) + 2*x^(1/2)
    sage: p - q
    -x^(1/3)
    sage: p * q
    x^(5/6) + x
    sage: (p / q).add_bigoh(4/3)
    x^(1/6) - x^(1/3) + x^(1/2) - x^(2/3) + x^(5/6) - x + x^(7/6) + O(x^(4/3))

Mind the base ring. However, the base ring can be changed. ::

    sage: I*q
    Traceback (most recent call last):
    ...
    TypeError: unsupported operand parent(s) for '*': 'Symbolic Ring' and 'Puiseux Series Ring in x over Rational Field'
    sage: I*q.change_ring(SR)
    I*x^(1/3) + I*x^(1/2)

Other properties of the Puiseux series can be easily obtained::

    sage: r = (3*x^(-1/5) + 7*x^(2/5) + (1/2)*x).add_bigoh(6/5); r
    3*x^(-1/5) + 7*x^(2/5) + 1/2*x + O(x^(6/5))
    sage: r.valuation()
    -1/5
    sage: r.prec()
    6/5
    sage: r.precision_absolute()
    6/5
    sage: r.precision_relative()
    7/5
    sage: r.exponents()
    [-1/5, 2/5, 1]
    sage: r.coefficients()
    [3, 7, 1/2]

Finally, Puiseux series are compatible with other objects in Sage. For example,
you can perform arithmetic with Laurent series. ::

    sage: L.<x> = LaurentSeriesRing(ZZ)
    sage: l = 3*x^(-2) + x^(-1) + 2 + x**3
    sage: r + l
    3*x^-2 + x^-1 + 3*x^(-1/5) + 2 + 7*x^(2/5) + 1/2*x + O(x^(6/5))
"""
from sage.arith.misc import lcm
from sage.ext.fast_callable import fast_callable
from sage.rings.big_oh import O
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.complex_field import ComplexField
from sage.rings.infinity import infinity
from sage.rings.laurent_series_ring_element cimport LaurentSeries
from sage.rings.laurent_series_ring import LaurentSeriesRing
from sage.rings.power_series_ring_element cimport PowerSeries
from sage.structure.element cimport (Element, ModuleElement,
                                     RingElement, AlgebraElement)


cdef class PuiseuxSeries(AlgebraElement):
    r"""
    We store a Puiseux series

    .. MATH::

        \sum_{n=-N}^\infty a_n x^{n/e}

    as a Laurent series

    .. MATH::

        \sum_{n=-N}^\infty a_n t^n

    where `t = x^{1/e}`.
    """
    @property
    def laurent_part(self):
        """
        Return the underlying Laurent series.

        EXAMPLES::

            sage: R.<x> = PuiseuxSeriesRing(QQ)
            sage: p = x^(1/2) + 3/4 * x^(2/3)
            sage: p.laurent_part()
            x^3 + 3/4*x^4
        """
        return self.__l

    def ramification_index(self):
        """
        Return the ramification index.

        EXAMPLES::

            sage: R.<x> = PuiseuxSeriesRing(QQ)
            sage: p = x^(1/2) + 3/4 * x^(2/3)
            sage: p.ramification_index()
            6
        """
        return self.__e

    def __init__(self, parent, f, e=1):
        r"""
        INPUT:

        parent : Ring, the target parent.

        f : one of the following types of inputs:

            - `PuiseuxSeries`
            - `LaurentSeries`

        EXAMPLES::

            sage: R.<x> = PuiseuxSeriesRing(QQ)
            sage: p = x^(1/2) + x**3; p
            x^(1/2) + x^3
        """
        AlgebraElement.__init__(self, parent)

        if isinstance(f, PuiseuxSeries):
            if (<PuiseuxSeries>f).__l._parent is parent.laurent_series_ring():
                l = (<PuiseuxSeries>f).__l
                e = (<PuiseuxSeries>f).__e
            else:
                l = parent.laurent_series_ring()((<PuiseuxSeries>f).__l)
                e = parent.laurent_series_ring()((<PuiseuxSeries>f).__e)
        else:
            l = parent.laurent_series_ring()(f)

        self.__l = l
        self.__e = long(abs(e))

    def __reduce__(self):
        return self.parent()(self.__l, self.__e)

    def _im_gens_(self, codomain, im_gens):
        return codomain(self(im_gens[0]))

    def _repr_(self):
        """
        Return a string representation.

        EXAMPLES::

            sage: R.<x> = PuiseuxSeriesRing(QQ)
            sage: p = x^(1/2) + x**3-x**(-1/4); p
            -x^(-1/4) + x^(1/2) + x^3
        """
        X = self.parent().variable_name()

        # extract coefficients and exponents of the laurent part.
        #
        # NOTE: self.__l.coefficients() is bugged when the coefficients are in
        # QQbar but coerced into SR. Therefore, we use self.__l.list() instead
        # (which works) and manually extract the coefficients and exponents
        lst = self.__l.list()
        val = self.valuation()
        coeff = []
        exp = []
        for n in range(len(lst)):
            c = lst[n]
            if not c.is_zero():
                coeff.append(c)
                exp.append(QQ(n) / self.__e + val)

        # print each term
        s = ''
        first = True
        for coeff, exp in zip(coeff, exp):
            # omit ' +' in the first term of the expression
            if first:
                s += str(coeff)
                first = False
            else:
                # if the coefficient itself is a sum (e.g. complex number or
                # expression) then wrap with parens
                coeff = str(coeff)
                if coeff[1:].find("+") != -1 or coeff[1:].find("-") != -1:
                    s += ' + (%s)' % coeff
                else:
                    s += ' + %s' % coeff

            # don't print (x-a)^0
            if exp:
                # don't print (x-a)^1
                if exp == 1:
                    s += '*%s' % X
                else:
                    # place parentheses around exponent if rational
                    s += '*%s^' % X
                    if exp.denominator() == 1:
                        s += str(exp)
                    else:
                        s += '(%s)' % exp

        # big oh
        prec = self.prec()
        if prec != infinity:
            prec = QQ(prec)
            if prec == 0:
                bigoh = 'O(1)'
            elif prec == 1:
                bigoh = 'O(%s)' % X
            elif prec.denominator() == 1:
                bigoh = 'O(%s^%s)' % (X, prec)
            else:
                bigoh = 'O(%s^(%s))' % (X, prec)

            if not s:
                return bigoh
            s += ' + %s' % bigoh

        # cleanup
        s = s.replace(' + -', ' - ')
        s = s.replace(' - -', ' + ')
        s = s.replace('1*','')
        s = s.replace('-1*', '-')
        return s

    def __call__(self, x):
        r"""
        Evaluate this Puiseux series.

        INPUT:

        x -- element of a ring

        EXAMPLES::

            sage: R.<x> = PuiseuxSeriesRing(QQ)
            sage: p = x^(1/2) + x**3-x**(-1/4)
            sage: p(16)
            8199/2
            sage: p(pi.n())
            32.0276049867404
        """
        # use x.nth_root since x**(1/self.__e) returns oo when x = 0
        if isinstance(x, int):
            x = ZZ(x)
        elif isinstance(x, float):
            x = ComplexField()(x)
        t = x.nth_root(self.__e)
        return self.__l.laurent_polynomial()(t)

    def _common_ramification_index(self, PuiseuxSeries right):
        r"""
        Return a ramification index common to self and right.

        In order to perform arithmetic on Puiseux series it is useful to find a
        common ramification index between two operands. That is, given Puiseux
        series :math:`p` and :math:`q` of ramification indices :math:`e` and
        :math:`f` we write both as series :math:`\tilde{f}` and
        :math:`\tilde{g}` in :math:`(x-a)^(1/g)` such that,

        .. MATH::

            f = \tilde{f}((x-a)^M), g = \tilde{g}((x-a)^N).

        INPUT:

        right : PuiseuxSeries

        OUTPUT:

        - g : int -- a ramification index common to self and right.
        - M, N : int -- scaling factors on self and right, respectively.

        """
        m = self.__e
        n = right.__e
        g = lcm(m, n)
        m = g / m
        n = g / n
        return g, m, n

    cpdef ModuleElement _add_(self, ModuleElement right_m):
        """
        Return the sum.

        EXAMPLES::

            sage: R.<x> = PuiseuxSeriesRing(QQ)
            sage: p = x^(1/2) + 3/4 * x^(2/3)
            sage: q = 2*x^(1/3) + 3/4 * x^(2/5)
            sage: p + q
            2*x^(1/3) + 3/4*x^(2/5) + x^(1/2) + 3/4*x^(2/3)
        """
        cdef PuiseuxSeries right = <PuiseuxSeries>right_m
        cdef LaurentSeries l, l1, l2
        cdef long g, m, n

        g, m, n = self._common_ramification_index(right)
        l1 = self.__l.V(m)
        l2 = right.__l.V(n)
        l = l1 + l2
        return PuiseuxSeries(self._parent, l, g)

    cpdef ModuleElement _sub_(self, ModuleElement right_m):
        """
        Return the difference.

        EXAMPLES::

            sage: R.<x> = PuiseuxSeriesRing(QQ)
            sage: p = x^(1/2) + 3/4 * x^(2/3)
            sage: q = 2*x^(1/3) + 3/4 * x^(2/5)
            sage: p - q
            -2*x^(1/3) - 3/4*x^(2/5) + x^(1/2) + 3/4*x^(2/3)
        """
        cdef PuiseuxSeries right = <PuiseuxSeries>right_m
        cdef LaurentSeries l, l1, l2
        cdef long g, m, n

        g, m, n = self._common_ramification_index(right)
        l1 = self.__l.V(m)
        l2 = right.__l.V(n)
        l = l1 - l2
        return PuiseuxSeries(self._parent, l, g)

    cpdef RingElement _mul_(self, RingElement right_r):
        """
        Return the product.

        EXAMPLES::

            sage: R.<x> = PuiseuxSeriesRing(QQ)
            sage: p = x^(1/2) + 3/4 * x^(2/3)
            sage: q = 2*x^(1/3) + 3/4 * x^(2/5)
            sage: p * q
            2*x^(5/6) + 3/4*x^(9/10) + 3/2*x + 9/16*x^(16/15)
        """
        cdef PuiseuxSeries right = <PuiseuxSeries>right_r
        cdef LaurentSeries l, l1, l2
        cdef long g, m, n

        g, m, n = self._common_ramification_index(right)
        l1 = self.__l.V(m)
        l2 = right.__l.V(n)
        l = l1 * l2
        return PuiseuxSeries(self._parent, l, g)

    cpdef ModuleElement _rmul_(self, RingElement c):
        return PuiseuxSeries(self._parent, self.__l._rmul_(c), self.__e)

    cpdef ModuleElement _lmul_(self, RingElement c):
        return PuiseuxSeries(self._parent, self.__l._lmul_(c), self.__e)

    cpdef RingElement _div_(self, RingElement right_r):
        """
        Return the quotient.

        EXAMPLES::

            sage: R.<x> = PuiseuxSeriesRing(QQ)
            sage: p = x^(1/2) + 3/4 * x^(2/3)
            sage: q = 2*x^(1/3) + 3/4 * x^(2/5)
            sage: p / q
            1/2*x^(1/6) - 3/16*x^(7/30) + 9/128*x^(3/10) + 3/8*x^(1/3) - 27/1024*x^(11/30) - 9/64*x^(2/5) + 81/8192*x^(13/30) + 27/512*x^(7/15) - 243/65536*x^(1/2) - 81/4096*x^(8/15) + 729/524288*x^(17/30) + 243/32768*x^(3/5) - 2187/4194304*x^(19/30) - 729/262144*x^(2/3) + 6561/33554432*x^(7/10) + 2187/2097152*x^(11/15) - 19683/268435456*x^(23/30) - 6561/16777216*x^(4/5) + O(x^(5/6))
        """
        cdef PuiseuxSeries right = <PuiseuxSeries>right_r
        cdef LaurentSeries l, l1, l2
        cdef long g, m, n

        g, m, n = self._common_ramification_index(right)
        l1 = self.__l.V(m)
        l2 = right.__l.V(n)
        l = l1 / l2
        return PuiseuxSeries(self._parent, l, g)

    def __pow__(_self, r, dummy):
        """
        Return the power.

        EXAMPLES::

            sage: R.<x> = PuiseuxSeriesRing(QQ)
            sage: p = x^(1/2) + 3/4 * x^(2/3)
            sage: p ** 3
            x^(3/2) + 9/4*x^(5/3) + 27/16*x^(11/6) + 27/64*x^2
        """
        cdef PuiseuxSeries self = _self
        cdef LaurentSeries l
        cdef long e

        r = QQ(r)
        numer = r.numerator()
        denom = r.denominator()

        # if the exponent is integral then do normal exponentiation
        if denom == 1:
            l = self.__l ** int(numer)
            e = self.__e
        # otherwise, we only exponentiate by a rational number if there is a
        # single term in the Puiseux series
        #
        # (I suppose we could use Taylor series expansions in the general case)
        else:
            if not self.is_monomial():
                raise ValueError('Can only exponentiate single '
                                 'term by rational')
            l = self.__l.V(numer)
            e = self.__e * int(denom)
        return PuiseuxSeries(self._parent, l, e)

    cpdef int _cmp_(self, Element right_r) except -2:
        r"""
        Comparison of self and right.

        As with Laurent series, two Puiseux series are equal if they agree for
        all coefficients up to the minimum of the precisions of each.

        EXAMPLES::

            sage: R.<x> = PuiseuxSeriesRing(QQ)
            sage: p = x^(1/2) + 3/4 * x^(2/3)
            sage: q = 2*x^(1/3) + 3/4 * x^(2/5)
            sage: p < q
            False
        """
        # scale each laurent series by their ramification indices and compare
        # the laurent series.
        cdef PuiseuxSeries right = <PuiseuxSeries>right_r
        exponents_l = self.exponents()
        exponents_r = right.exponents()
        coefficients_l = self.coefficients()
        coefficients_r = right.coefficients()
        d = min(len(exponents_l), len(exponents_r))

        # first compare each exponent and then each coefficient
        for i in range(d):
            c = cmp(exponents_l[i], exponents_r[i])
            if c:
                return c
            c = cmp(coefficients_l[i], coefficients_r[i])
            if c:
                return c
        return 0

    def __lshift__(self, r):
        return self.shift(r)

    def __rshift__(self, r):
        return self.shift(-r)

    def __nonzero__(self):
        """
        Return whether this is not zero.

        EXAMPLES::

            sage: R.<x> = PuiseuxSeriesRing(QQ)
            sage: p = x^(1/2) + 3/4 * x^(2/3)
            sage: p != 0
            True
            sage: R.zero() != 0
            False
        """
        return not not self.__l

    def __hash__(self):
        """
        Return a hash of self.

        EXAMPLES::

            sage: R.<x> = PuiseuxSeriesRing(QQ)
            sage: p = x^(-7/2) + 3 + 5*x^(1/2) - 7*x**3
            sage: hash(p)  # indirect doctest
            -15360174648385722
        """
        return hash(self.__l) ^ self.__e

    def __getitem__(self, r):
        r"""
        Return the coefficient with exponent `r`.

        EXAMPLES::

            sage: R.<x> = PuiseuxSeriesRing(QQ)
            sage: p = x^(-7/2) + 3 + 5*x^(1/2) - 7*x**3
            sage: p[-7/2]
            1
            sage: p[0]
            3
            sage: p[1/2]
            5
            sage: p[3]
            -7
            sage: p[100]
            0
        """
        if isinstance(r, slice):
            start, stop, step = r.start, r.stop, r.step
            n = slice(start * self.__e, stop * self.__e, step * self.__e)
            return PuiseuxSeries(self._parent, self.__l[start:stop:step], self.__e)
        else:
            n = int(r * self.__e)
            return self.__l[n]

    def __iter__(self):
        """
        Return an iterator over the coefficients.

        EXAMPLES::

            sage: R.<x> = PuiseuxSeriesRing(QQ)
            sage: p = x^(-7/2) + 3 + 5*x^(1/2) - 7*x**3
            sage: list(p)
            [1, 0, 0, 0, 0, 0, 0, 3, 5, 0, 0, 0, 0, -7]
        """
        return iter(self.__l)

    def __copy__(self):
        """
        Return a copy of self.
        """
        return PuiseuxSeries(self._parent, self.__l.copy(), self.__e)

    def valuation(self):
        """
        Return the valuation of self.

        EXAMPLES::

            sage: R.<x> = PuiseuxSeriesRing(QQ)
            sage: p = x^(-7/2) + 3 + 5*x^(1/2) - 7*x**3
            sage: p.valuation()
            -7/2
        """
        val = self.__l.valuation() / QQ(self.__e)
        if val == infinity:
            return 0
        return val

    def add_bigoh(self, prec):
        """
        Return the truncated series at chosen precision.

        EXAMPLES::

            sage: R.<x> = PuiseuxSeriesRing(QQ)
            sage: p = x^(-7/2) + 3 + 5*x^(1/2) - 7*x**3
            sage: p.add_bigoh(2)
            x^(-7/2) + 3 + 5*x^(1/2) + O(x^2)
            sage: p.add_bigoh(0)
            x^(-7/2) + O(1)
            sage: p.add_bigoh(-1)
            x^(-7/2) + O(x^-1)
        """
        if prec == infinity or prec >= self.prec():
            return self

        # the following is here due to a bug in Sage: adding a bigoh of order
        # less than the series precision raises an error due to attempting to
        # build an underlying power series with negative precision. this fix
        # makes sure that if the requested precision is less than that if the
        # Laurent series it will just return a bigoh
        l_prec = int(prec * self.__e)
        try:
            l = self.__l.add_bigoh(l_prec)
        except ValueError:
            x = self.__l.parent().gen()
            l = O(x ** l_prec)
        return PuiseuxSeries(self._parent, l, self.__e)

    def change_ring(self, R):
        """
        Return self over a the new ring R.

        EXAMPLES::

            sage: R.<x> = PuiseuxSeriesRing(ZZ)
            sage: p = x^(-7/2) + 3 + 5*x^(1/2) - 7*x**3
            sage: q = p.change_ring(QQ); q
            x^(-7/2) + 3 + 5*x^(1/2) - 7*x^3
            sage: q.parent()
            Puiseux Series Ring in x over Rational Field
        """
        return self._parent.change_ring(R)(self)

    def is_unit(self):
        """
        Return whether self is a unit.

        A Puiseux series is a unit if and only if its leading coefficient is.

        EXAMPLES::

            sage: R.<x> = PuiseuxSeriesRing(ZZ)
            sage: p = x^(-7/2) + 3 + 5*x^(1/2) - 7*x**3
            sage: p.is_unit()
            True
            sage: q = 4 * x^(-7/2) + 3 * x**4
            sage: q.is_unit()
            False
        """
        return self.__l.is_unit()

    def is_zero(self):
        """
        Return whether self is zero.

        EXAMPLES::

            sage: R.<x> = PuiseuxSeriesRing(QQ)
            sage: p = x^(1/2) + 3/4 * x^(2/3)
            sage: p.is_zero()
            False
            sage: R.zero().is_zero()
            True
        """
        return self.__l.is_zero()

    def is_monomial(self):
        """
        Return whether self is a monomial.

        This is True if and only if self is `x^p` for some rational `p`.

        EXAMPLES::

            sage: R.<x> = PuiseuxSeriesRing(QQ)
            sage: p = x^(1/2) + 3/4 * x^(2/3)
            sage: p.is_monomial()
            False
            sage: q = x**(11/13)
            sage: q.is_monomial()
            True
            sage: q = 4*x**(11/13)
            sage: q.is_monomial()
            False
        """
        return self.__l.is_monomial()

    def list(self):
        """
        Return the list of coefficients ?
        """
        return self.__l.list()

    def coefficients(self):
        """
        Return coefficients and exponents of the Laurent part.

        .. NOTE::

            self.__l.coefficients() is bugged when the coefficients
            are in QQbar but coerced into SR. Therefore, we use
            self.__l.list() instead (which works) and manually extract
            the coefficients and exponents
        """
        lst = self.__l.list()
        val = self.valuation()
        coeff = []
        for n in range(len(lst)):
            c = lst[n]
            if not c.is_zero():
                coeff.append(c)
        return coeff

    def exponents(self):
        """
        Return exponents of the Laurent part.

        NOTE: self.__l.coefficients() is bugged when the coefficients are in
        QQbar but coerced into SR. Therefore, we use self.__l.list() instead
        (which works) and manually extract the coefficients and exponents
        """
        lst = self.__l.list()
        val = self.valuation()
        exp = []
        for n in range(len(lst)):
            c = lst[n]
            if not c.is_zero():
                exp.append(QQ(n) / self.__e + val)
        return exp

    def __setitem__(self, n, value):
        raise IndexError('Puiseux series are immutable')

    def degree(self):
        """
        Return the degree of self.
        """
        return self.__l.degree() / self.__e

    def shift(self, r):
        r"""
        Return this Puiseux series multiplied by `x^r`.
        """
        cdef LaurentSeries l = self.__l.shift(r * self.__e)
        return PuiseuxSeries(self._parent, l, self.__e)

    def truncate(self, r):
        r"""
        Return the Puiseux series of degree ` < r`.

        This is equivalent to self modulo `x^r`.
        """
        l = self.__l.truncate(r * self.__e)
        return PuiseuxSeries(self._parent, l, self.__e)

    def prec(self):
        """
        Return the precision of self.
        """
        if self.__l.prec() == infinity:
            return infinity
        return self.__l.prec() / self.__e

    precision_absolute = prec

    def precision_relative(self):
        r"""
        Return the relative precision of the series.

        The relative precision of the Puiseux series is the difference between
        its absolute precision and its valuation.
        """
        if self.is_zero():
            return 0
        return self.prec() - self.valuation()

    def common_prec(self, p):
        r"""
        Return the minimum precision of `p` and self.
        """
        if self.prec() is infinity:
            return p.prec()
        elif p.prec() is infinity:
            return self.prec()
        return min(self.prec(), p.prec())

    def variable(self):
        r"""
        Return the variable of self.

        EXAMPLES::

            sage: R.<x> = PuiseuxSeriesRing(QQ)
            sage: p = x^(-7/2) + 3 + 5*x^(1/2) - 7*x**3
            sage: p.variable()
            'x'
        """
        return self._parent.variable_name()

    def laurent_series(self):
        r"""
        If self is a Laurent series, return it as a Laurent series.
        """
        if self.__e != 1:
            raise ArithmeticError('self is not a Laurent series')
        return self.__l

    def power_series(self):
        r"""
        If self is a power series, return it as a power series.
        """
        try:
            l = self.laurent_series()
            return l.power_series()
        except:
            raise ArithmeticError('self is not a power series')

    def inverse(self):
        r"""
        Return the inverse of self.

        EXAMPLES::

            sage: R.<x> = PuiseuxSeriesRing(QQ)
            sage: p = x^(-7/2) + 3 + 5*x^(1/2) - 7*x**3
            sage: 1/p
            x^(7/2) - 3*x^7 - 5*x^(15/2) + 7*x^10 + 9*x^(21/2) + 30*x^11 +
            25*x^(23/2) + O(x^(27/2))
        """
        return self.__invert__()
