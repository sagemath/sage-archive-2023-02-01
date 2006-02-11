"""
Laurent Series Rings
"""

import weakref

from rational_field import Q
import laurent_series_ring_element
import power_series_ring
import commutative_ring
import integral_domain
import field


_objsLaurentSeriesRing = {}
def LaurentSeriesRing(base_ring, name=None):
    """
    EXAMPLES:
        sage: R = LaurentSeriesRing(Q, 'x'); R
        Laurent Series Ring in x over Rational Field
        sage: x = R.0
        sage: g = 1 - x + x^2 - x^4 +O(x^8); g
        1 - x + x^2 - x^4 + O(x^8)
        sage: g = 10*x^(-3) + 2006 - 19*x + x^2 - x^4 +O(x^8); g
        10*x^-3 + 2006 - 19*x + x^2 - x^4 + O(x^8)

    You can also use more mathematical notation when the base is a field:
        sage: Frac(Q[['x']])
        Laurent Series Ring in x over Rational Field
        sage: Frac(GF(5)['y'])
        Fraction Field of Univariate Polynomial Ring in y over Finite Field of size 5

    Here the fraction field is not just the Laurent series ring, so you can't
    use the \code{Frac} notation to make the Laurent series ring.
        sage: Frac(Z[['t']])
        Fraction Field of Power Series Ring in t over Integer Ring
    """
    global _objsLaurentSeriesRing
    key = (base_ring, name)
    if _objsLaurentSeriesRing.has_key(key):
        x = _objsLaurentSeriesRing[key]()
        if x != None: return x

    if isinstance(base_ring, field.Field):
        R = LaurentSeriesRing_field(base_ring, name)
    elif isinstance(base_ring, integral_domain.IntegralDomain):
        R = LaurentSeriesRing_domain(base_ring, name)
    elif isinstance(base_ring, commutative_ring.CommutativeRing):
        R = LaurentSeriesRing_generic(base_ring, name)
    else:
        raise TypeError, "base_ring must be a commutative ring"
    _objsLaurentSeriesRing[key] = weakref.ref(R)
    return R

def is_LaurentSeriesRing(x):
    return isinstance(x, LaurentSeriesRing_generic)

class LaurentSeriesRing_generic(commutative_ring.CommutativeRing):
    """
    Univariate Laurent Series Ring
    EXAMPLES:
        sage: K, q = LaurentSeriesRing(C, 'q').objgen(); K
        Laurent Series Ring in q over Complex Field with 53 bits of precision
        sage: loads(K.dumps()) == K
        True
    """

    def __init__(self, base_ring, name=None):
        self.__base_ring = base_ring
        self.assign_names(name)

    def __reduce__(self):
        return self.__class__, (self.__base_ring, self.variable_name())

    def __repr__(self):
        return "Laurent Series Ring in %s over %s"%(self.variable_name(), self.base_ring())

    def __call__(self, x, n=0):
        if isinstance(x, laurent_series_ring_element.LaurentSeries) and n==0 and self == x.parent():
            return x
        return laurent_series_ring_element.LaurentSeries(self, x, n)

    def __cmp__(self, other):
        if not isinstance(other, LaurentSeriesRing_generic):
            return -1
        if self.base_ring() != other.base_ring():
            return -1
        if self.variable_name() != other.variable_name():
            return -1
        return 0


    def _is_valid_homomorphism_(self, codomain, im_gens):
        ## NOTE: There are no ring homomorphisms from the ring of
        ## all formal power series to most rings, e.g, the p-adic
        ## field, since you can always (mathematically!) construct
        ## some power series that doesn't converge.
        ## Note that 0 is not a *ring* homomorphism.
        from power_series_ring import is_PowerSeriesRing
        if is_PowerSeriesRing(codomain) or is_LaurentSeriesRing(codomain):
            return im_gens[0].valuation() > 0
        return False

    def characteristic(self):
        return self.base_ring().characteristic()

    def set_default_prec(self, n):
        self.power_series_ring().set_default_prec(n)

    def default_prec(self):
        return self.power_series_ring().default_prec()

    def base_ring(self):
        """
        EXAMPLES:
            sage: R = LaurentSeriesRing(QQ, "x")
            sage: R.base_ring()
            Rational Field
            sage: S = LaurentSeriesRing(GF(17)['x'], 'y')
            sage: S
            Laurent Series Ring in y over Univariate Polynomial Ring in x over
            Finite Field of size 17
            sage: S.base_ring()
            Univariate Polynomial Ring in x over Finite Field of size 17
        """
        return self.__base_ring

    def gen(self, n=0):
        if n != 0:
            raise IndexError, "Generator %s not defined."%n
        try:
            return self.__generator
        except AttributeError:
            self.__generator = laurent_series_ring_element.LaurentSeries(self, [0,1])
            return self.__generator

    def ngens(self):
        return 1

    def power_series_ring(self):
        r"""
        If this is the Laurent series ring $R((t))$, return the power
        series ring $R[[t]]$.

        EXAMPLES:
            sage: R = LaurentSeriesRing(QQ, "x")
            sage: R.power_series_ring()
            Power Series Ring in x over Rational Field
        """
        try:
            return self.__power_series_ring
        except AttributeError:
            self.__power_series_ring = power_series_ring.PowerSeriesRing(
                                         self.base_ring(), self.variable_name())
            return self.__power_series_ring

class LaurentSeriesRing_domain(integral_domain.IntegralDomain, LaurentSeriesRing_generic):
    def __init__(self, base_ring, name=None):
        LaurentSeriesRing_generic.__init__(self, base_ring, name)

    def fraction_field(self):
        try:
            return self.__fraction_field
        except AttributeError:
            self.__fraction_field = LaurentSeriesRing(self.base_ring().fraction_field(), name)
            return self.__fraction_field

class LaurentSeriesRing_field(field.Field, LaurentSeriesRing_generic):
    def __init__(self, base_ring, name=None):
        LaurentSeriesRing_generic.__init__(self, base_ring, name)

