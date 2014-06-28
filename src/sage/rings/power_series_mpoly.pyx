# NOT ready to be used -- possibly should be deleted.

include "sage/ext/stdsage.pxi"

from power_series_ring_element cimport PowerSeries
from sage.structure.element cimport Element, ModuleElement, RingElement
from infinity import infinity, is_Infinite
import arith
from sage.libs.all import PariError
from power_series_ring_element import is_PowerSeries
import rational_field
from polynomial.multi_polynomial_ring_generic import is_MPolynomialRing
import power_series_poly

cdef class PowerSeries_mpoly(PowerSeries):

    def __init__(self, parent, f=0, prec=infinity, int check=1, is_gen=0):
        """
        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: R.<y> = S[[]]
            sage: f = x + 2*y + x*y
            sage: loads(f.dumps()) == f
            True
        """
        S = parent._mpoly_ring()
        if PY_TYPE_CHECK(f, Element) and (<Element>f)._parent is S:
            #if check and not (prec is infinity):
            #    self.__f = f.truncate(S.gens()[-1], prec)
            #    self._truncated = 1
            #else:
            self.__f = f
        else:
            # We use the generic code, since the coercion rules can be
            # very complicated.  This is non-optimal, but much easier
            # to maintain.
            g = power_series_poly.PowerSeries_poly(parent, f=f,
                                      prec=prec, check=check).polynomial()

            # Now g is a polynomial in the indeterminate over the base
            # ring.  We have to construct a multivariate polynomial
            # from g in S efficiently.

            # Let d be the dictionary that will represent this object
            # that we're creating.  We compute d explicitly below.

            v = g.list()
            # Take each of the coefficients of g, make into a polydict,
            # and then create corresponding entries of d.
            B = parent.base_ring()
            i = S.ngens() - 1

            # We divide the computation of d into 2 cases in order to
            # avoid having an if statement in the inner loop of a
            # doubly-nested for loop.
            d = {}
            if is_MPolynomialRing(B):
                for i in range(len(v)):
                    for n, c in v[i].dict().iteritems():
                        d[tuple(n) + (i,)] = c
            else:
                for i in range(len(v)):
                    for n, c in v[i].dict().iteritems():
                        d[(n,i)] = c

            self.__f = S(d)

        PowerSeries.__init__(self, parent, prec, is_gen)

    def __reduce__(self):
        # do *not* delete old versions.
        return make_powerseries_mpoly_v0, (self._parent, self.__f, self._prec, self.__is_gen)

    def __call__(self, *args, **kwds):
        if len(kwds) == 0 and len(args) == 1:
            R = self.parent()._mpoly_ring()
            return self.__f.substitute({R.gen(0):args[0]})
        else:
            return self.__f(*args, **kwds)

    def do_truncation(self):
        if self._truncated:
            return
        S = self.parent()._mpoly_ring()
        self.__f = self.__f.truncate(S.gens()[-1], self._prec)
        self._truncated = 1

    def _repr_(self):
        if not self._truncated:
            self.do_truncation()
        return PowerSeries._repr_(self)

    def list(self):
        if self.__list is None:
            self.__list = self.polynomial().list()
        return self.__list

    def polynomial(self):
        if self.__poly is None:
            S = self.parent()._mpoly_ring()
            self.__poly = self.__f.polynomial(S.gens()[-1])
        return self.__poly

    def _mpoly(self):
        return self.__f

    cpdef RingElement _mul_(self, RingElement right_r):
        """
        Return the product of two power series.
        """
        prec = self._mul_prec(right_r)
        return PowerSeries_mpoly(self._parent,
                                 self.__f * (<PowerSeries_mpoly>right_r).__f,
                                 prec = prec,
                                 check =True)


    def __iter__(self):
        """
        Return an interator over the coefficients of this power series.
        """
        return iter(self.__f)

    def __neg__(self):
        """
        Return the negative of this power series.
        """
        return PowerSeries_mpoly(self._parent, -self.__f,
                                         self._prec, check=False)

    cpdef ModuleElement _add_(self, ModuleElement right_m):
        """
        EXAMPLES:
        """
        cdef PowerSeries_mpoly right = <PowerSeries_mpoly>right_m
        return PowerSeries_mpoly(self._parent, self.__f + right.__f, \
                                         self.common_prec_c(right), check=True)

    cpdef ModuleElement _sub_(self, ModuleElement right_m):
        """
        Return difference of two power series.

        EXAMPLES:
        """
        cdef PowerSeries_mpoly right = <PowerSeries_mpoly>right_m
        return PowerSeries_mpoly(self._parent, self.__f - right.__f, \
                                         self.common_prec_c(right), check=True)

    cpdef ModuleElement _rmul_(self, RingElement c):
        return PowerSeries_mpoly(self._parent, self.__f._rmul_(c), self._prec, check=False)

    cpdef ModuleElement _lmul_(self, RingElement c):
        return PowerSeries_mpoly(self._parent, self.__f._lmul_(c), self._prec, check=False)


def make_powerseries_mpoly_v0(parent,  f, prec, is_gen):
    return PowerSeries_mpoly(parent, f, prec, 0, is_gen)
