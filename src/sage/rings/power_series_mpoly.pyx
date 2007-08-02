include "../ext/stdsage.pxi"

from power_series_ring_element cimport PowerSeries
from sage.structure.element cimport Element, ModuleElement, RingElement
from infinity import infinity, is_Infinite
import arith
from sage.libs.all import PariError
from power_series_ring_element import is_PowerSeries
import rational_field

cdef class PowerSeries_mpoly(PowerSeries):

    def __init__(self, parent, f=0, prec=infinity, int check=1, is_gen=0):
        """
        EXAMPLES:
            sage: S.<x> = QQ
            sage: R.<y> = S[[]]
            sage: f = x + 2*y + x*y
            sage: loads(f.dumps()) == f
            True
        """
        R = parent._mpoly_ring()
        if PY_TYPE_CHECK(f, Element):
            if (<Element>f)._parent is R:
                pass

    def __call__(self, *args, **kwds):
        """
        EXAMPLE:
            sage: S.<x> = GF(7)
            sage: R.<t> = S[[]]
            sage: f = 3 - x*t^3 + O(t^5)
            sage: f(1)
            2
        """
        if len(kwds) == 0 and len(args) == 1:
            R = self.parent()._mpoly_ring()
            return self.__f.substitute({R.gen(0):args[0]})
        else:
            return self.__f(*args, **kwds)

    def list(self):
        if self.__list is None:
            R = self.parent().base_ring()
            self.__list = self.__f.list(dense=True, ring=R)
        return self.__list

##         x = parent._mpoly_ring().gen(0)
##         d = self.__f.degree(x)
##         v = [R(0)]*(d+1)
##         for i from 0 <= i <= d:
##             v[i] =

##         return self.__list

    def polynomial(self):
        if self.__poly is None:
            R = self.parent()._poly_ring()
            self.__poly = R(self.list())
        return self.__poly



def make_powerseries_poly_v0(parent,  f, prec, is_gen):
    return PowerSeries_mpoly(parent, f, prec, 0, is_gen)
