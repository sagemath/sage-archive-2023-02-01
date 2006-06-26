"""
Univariate Power Series Rings

EXAMPLES:
    sage: R.<t> = PowerSeriesRing(RationalField())
    sage: R.random_element(6)
    -t - t^2 - t^3 - t^4 + O(t^6)

    sage: S = R([1, 3, 5, 7], 10); S
    1 + 3*t + 5*t^2 + 7*t^3 + O(t^10)

    sage: S.truncate(3)
    5*t^2 + 3*t + 1

AUTHOR:
    -- William Stein: the code
    -- Jeremy Cho (2006-05-17): some examples (above)
"""

import weakref
import power_series_ring_element
import polynomial_ring
import laurent_series_ring
import commutative_ring
import integral_domain
import field
import integer
import sage.structure.gens as gens
from infinity import infinity
import sage.misc.latex as latex
from sage.structure.nonexact import Nonexact

from sage.interfaces.magma import MagmaElement
from sage.misc.sage_eval import sage_eval

#_objsPowerSeriesRing = {}
def PowerSeriesRing(base_ring, name=None, default_prec=20):
    #global _objsPowerSeriesRing
    #key = (base_ring, name, default_prec)
    #if _objsPowerSeriesRing.has_key(key):
    #    x = _objsPowerSeriesRing[key]()
    #    if x != None: return x
    if isinstance(base_ring, field.Field):
        R = PowerSeriesRing_over_field(base_ring, name, default_prec)
    elif isinstance(base_ring, integral_domain.IntegralDomain):
        R = PowerSeriesRing_domain(base_ring, name, default_prec)
    elif isinstance(base_ring, commutative_ring.CommutativeRing):
        R = PowerSeriesRing_generic(base_ring, name, default_prec)
    else:
        raise TypeError, "base_ring must be a commutative ring"
    #_objsPowerSeriesRing[key] = weakref.ref(R)
    return R

def is_PowerSeriesRing(x):
    return isinstance(x, PowerSeriesRing_generic)

class PowerSeriesRing_generic(commutative_ring.CommutativeRing, Nonexact):
    def __init__(self, base_ring, name=None, default_prec=20):
        Nonexact.__init__(self, default_prec)
        self.__base_ring = base_ring
        self.__poly_ring = polynomial_ring.PolynomialRing(base_ring, name)
        self.__power_series_class = power_series_ring_element.PowerSeries_generic_dense
        self.__generator = self.__power_series_class(self, [0,1], check=True, is_gen=True)
        self.assign_names(name)

    def assign_names(self, name):
        commutative_ring.CommutativeRing.assign_names(self, name)
        self.__poly_ring.assign_names(name)

    def __repr__(self):
        return "Power Series Ring in %s over %s"%(self.variable_name(), self.base_ring())

    def _latex_(self):
        return "%s[[%s]]"%(latex.latex(self.base_ring()), self.variable_name())

    def __call__(self, f, prec=infinity, check=True):
        if isinstance(f, power_series_ring_element.PowerSeries) and f.parent() == self:
            if prec >= f.prec():
                return f
            f = f.truncate(prec)
        elif isinstance(f, MagmaElement) and str(f.Type()) == 'RngSerPowElt':
            v = sage_eval(f.Eltseq())
            return self(v) * (self.gen(0)**f.Valuation())
        return self.__power_series_class(self, f, prec, check=check)

    def _is_valid_homomorphism_(self, codomain, im_gens):
        ## NOTE: There are no ring homomorphisms from the ring of
        ## all formal power series to most rings, e.g, the p-adic
        ## field, since you can always (mathematically!) construct
        ## some power series that doesn't converge.
        ## Note that 0 is not a *ring* homomorphism.
        from laurent_series_ring import is_LaurentSeriesRing
        if is_PowerSeriesRing(codomain) or is_LaurentSeriesRing(codomain):
            return im_gens[0].valuation() > 0
        return False

    def base_ring(self):
        return self.__base_ring

    def _poly_ring(self):
        return self.__poly_ring

    def gen(self, n=0):
        if n != 0:
            raise IndexError, "Generator %s not defined."%n
        return self.__generator

    def ngens(self):
        return 1

    def random_element(self, prec, bound=0):
        """
        Return a random power series.

        INPUT:
            prec -- an int
            bound -- an int (default: 0, which tries to spread choice across ring, if implemented)

        OUTPUT:
            Polynomial -- A polynomial such that the coefficient of x^i,
            for i up to degree, are coercisions to the base ring of
            random integers between -bound and bound.
        """
        return self(self.__poly_ring.random_element(prec, bound), prec)

    def __cmp__(self, other):
        if not isinstance(other, PowerSeriesRing_generic):
            return -1
        if self.base_ring() != other.base_ring():
            return -1
        if self.variable_name() != other.variable_name():
            return -1
        return 0

    def __contains__(self, x):
        if not isinstance(x, PowerSeries):
            return False
        return x.parent() == self

    def is_atomic_repr(self):
        return False

    def is_field(self):
        return False

    def is_finite(self):
        return False

    def characteristic(self):
        return self.base_ring().characteristic()

    def laurent_series_ring(self):
        """
        If this is the power series ring R[[t]], return the Laurent
        series ring R((t)).
        """
        try:
            return self.__laurent_series_ring
        except AttributeError:
            self.__laurent_series_ring = laurent_series_ring.LaurentSeriesRing(
                                                 self.base_ring(), self.variable_name())
            return self.__laurent_series_ring

class PowerSeriesRing_domain(integral_domain.IntegralDomain, PowerSeriesRing_generic):
    def __init__(self, base_ring, name=None, default_prec=20):
        PowerSeriesRing_generic.__init__(self, base_ring, name, default_prec)


class PowerSeriesRing_over_field(PowerSeriesRing_domain):
    def __init__(self, base_ring, name=None, default_prec=20):
        PowerSeriesRing_generic.__init__(self, base_ring, name, default_prec)

    def fraction_field(self):
        return self.laurent_series_ring()

