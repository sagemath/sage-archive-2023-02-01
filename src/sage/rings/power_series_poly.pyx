include "../ext/stdsage.pxi"

from power_series_ring_element cimport PowerSeries
from sage.structure.element cimport Element, ModuleElement, RingElement
from infinity import infinity, is_Infinite
import arith
from sage.libs.all import PariError
from power_series_ring_element import is_PowerSeries
import rational_field

cdef class PowerSeries_poly(PowerSeries):

    def __init__(self, parent, f=0, prec=infinity, int check=1, is_gen=0):
        """
        EXAMPLES:
            sage: R, q = PowerSeriesRing(CC, 'q').objgen()
            sage: R
            Power Series Ring in q over Complex Field with 53 bits of precision
            sage: loads(q.dumps()) == q
            True

            sage: R.<t> = QQ[[]]
            sage: f = 3 - t^3 + O(t^5)
            sage: a = f^3; a
            27 - 27*t^3 + O(t^5)
            sage: b = f^-3; b
            1/27 + 1/27*t^3 + O(t^5)
            sage: a*b
            1 + O(t^5)
        """
        R = parent._poly_ring()
        if PY_TYPE_CHECK(f, Element):
            if (<Element>f)._parent is R:
                pass
            elif (<Element>f)._parent == R.base_ring():
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
            self.__f = self.__f.truncate_c(prec)
        PowerSeries.__init__(self, parent, prec, is_gen)

    def __hash__(self):
        return hash(self.__f)

    def __reduce__(self):
        # do *not* delete old versions.
        return make_powerseries_poly_v0, (self._parent, self.__f, self._prec, self.__is_gen)

    def __richcmp__(left, right, int op):
        return (<Element>left)._richcmp(right, op)

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

    def __nonzero__(self):
        return not not self.__f


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

    cdef ModuleElement _iadd_c_impl(self, ModuleElement right_m):
        """
        EXAMPLES:
            sage: R.<x> = PowerSeriesRing(ZZ)
            sage: f = x^4
            sage: f += x; f
            x + x^4
            sage: g = x^2 + O(x^3); g
            x^2 + O(x^3)
            sage: f += g; f
            x + x^2 + O(x^3)
        """
        cdef PowerSeries_poly right = <PowerSeries_poly>right_m
        self.__f += right.__f
        if self._prec is not infinity:
            if self._prec < right._prec:
                self.__f = self.__f._inplace_truncate(self._prec)
            else:
                self.__f = self.__f._inplace_truncate(right._prec)
                self._prec = right._prec
        elif right._prec is not infinity:
            self.__f = self.__f._inplace_truncate(right._prec)
            self._prec = right._prec
        return self

    cdef ModuleElement _sub_c_impl(self, ModuleElement right_m):
        """
        Return difference of two power series.

        EXAMPLES:
            sage: k.<w> = ZZ[]
            sage: R.<t> = k[[]]
            sage: w*t^2 -w*t +13 - (w*t^2 + w*t)
            13 - 2*w*t
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
        prec = self._mul_prec(right_r)
        return PowerSeries_poly(self._parent,
                                self.__f * (<PowerSeries_poly>right_r).__f,
                                prec = prec,
                                check = True)  # check, since truncation may be needed

    cdef RingElement _imul_c_impl(self, RingElement right_r):
        """
        Return the product of two power series.

        EXAMPLES:
            sage: k.<w> = ZZ[[]]
            sage: (1+17*w+15*w^3+O(w^5))*(19*w^10+O(w^12))
            19*w^10 + 323*w^11 + O(w^12)
        """
        prec = self._mul_prec(right_r)
        self.__f *= (<PowerSeries_poly>right_r).__f
        if prec is not infinity:
            self.__f = self.__f._inplace_truncate(prec)
            self._prec = prec
        return self

    cdef ModuleElement _rmul_c_impl(self, RingElement c):
        return PowerSeries_poly(self._parent, self.__f._rmul_c(c), self._prec, check=False)

    cdef ModuleElement _lmul_c_impl(self, RingElement c):
        return PowerSeries_poly(self._parent, self.__f._lmul_c(c), self._prec, check=False)

    cdef ModuleElement _ilmul_c_impl(self, RingElement c):
#        print "f", type(self.__f), self.__f
#        print "c", type(c), c
#        print "f*c", type(self.__f*c), self.__f*c
#        ff = self.__f
        self.__f *= c
#        ff *= c
#        print "ff", type(ff), ff
#        self.__f = ff
#        self.__f = self.__f * c
        return self


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

    def __lshift__(PowerSeries_poly self, n):
        if n:
            return PowerSeries_poly(self._parent, self.__f << n, self._prec + n)
        else:
            return self

    def __rshift__(PowerSeries_poly self, n):
        if n:
            return PowerSeries_poly(self._parent, self.__f >> n, self._prec - n)
        else:
            return self

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
            return self.__f
        else:
            return self.__f.truncate_c(prec)

    cdef _inplace_truncate(self, long prec):
        self.__f = self.__f._inplace_truncate(prec)
        self.prec = prec
        return self

    def truncate_powerseries(self, long prec):
        r"""
        Returns the power series of degree $ < n$ which is equivalent to self
        modulo $x^n$.
        """
        return PowerSeries_poly(self._parent, self.__f.truncate_c(prec), self._prec if self._prec < prec else infinity, check=False)

    def copy(self):
        return PowerSeries_poly(self._parent, self.__f, self._prec, check=False)

    def list(self):
        return self.__f.list()

    def dict(self):
        return self.__f.dict()

    def _derivative(self, var=None):
        """
        Return the derivative of this power series with respect
        to the variable var.

        If var is None or is the generator of this ring, we take the derivative
        with respect to the generator.

        Otherwise, we call _derivative(var) on each coefficient of
        the series.

        SEE ALSO:
            self.derivative()

        EXAMPLES:
            sage: R.<t> = PowerSeriesRing(QQ, sparse=True)
            sage: f = 2 + 3*t^2 + t^100000 + O(t^10000000); f
            2 + 3*t^2 + t^100000 + O(t^10000000)
            sage: f._derivative()
            6*t + 100000*t^99999 + O(t^9999999)
            sage: f._derivative(t)
            6*t + 100000*t^99999 + O(t^9999999)

            sage: R.<x> = PolynomialRing(ZZ)
            sage: S.<y> = PowerSeriesRing(R, sparse=True)
            sage: f = x^3*y^4 + O(y^5)
            sage: f._derivative()
            4*x^3*y^3 + O(y^4)
            sage: f._derivative(y)
            4*x^3*y^3 + O(y^4)
            sage: f._derivative(x)
            3*x^2*y^4 + O(y^5)
        """
        if var is not None and var is not self._parent.gen():
            # call _derivative() recursively on coefficients
            return PowerSeries_poly(self._parent, self.__f._derivative(var),
                                    self.prec(), check=False)

        # compute formal derivative with respect to generator
        return PowerSeries_poly(self._parent, self.__f._derivative(),
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

def make_powerseries_poly_v0(parent,  f, prec, is_gen):
    return PowerSeries_poly(parent, f, prec, 0, is_gen)
