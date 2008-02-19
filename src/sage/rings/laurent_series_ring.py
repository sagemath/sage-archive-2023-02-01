"""
Laurent Series Rings

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

import weakref

import laurent_series_ring_element
import power_series_ring
import commutative_ring
import integral_domain
import field

from sage.structure.parent_gens import ParentWithGens


laurent_series = {}
def LaurentSeriesRing(base_ring, name=None, names=None, sparse=False):
    """
    EXAMPLES:
        sage: R = LaurentSeriesRing(QQ, 'x'); R
        Laurent Series Ring in x over Rational Field
        sage: x = R.0
        sage: g = 1 - x + x^2 - x^4 +O(x^8); g
        1 - x + x^2 - x^4 + O(x^8)
        sage: g = 10*x^(-3) + 2006 - 19*x + x^2 - x^4 +O(x^8); g
        10*x^-3 + 2006 - 19*x + x^2 - x^4 + O(x^8)

    You can also use more mathematical notation when the base is a field:
        sage: Frac(QQ[['x']])
        Laurent Series Ring in x over Rational Field
        sage: Frac(GF(5)['y'])
        Fraction Field of Univariate Polynomial Ring in y over Finite Field of size 5

    Here the fraction field is not just the Laurent series ring, so you can't
    use the \code{Frac} notation to make the Laurent series ring.
        sage: Frac(ZZ[['t']])
        Fraction Field of Power Series Ring in t over Integer Ring

    Laurent series rings are determined by their variable and the base ring,
    and are globally unique.
        sage: K = Qp(5, prec = 5)
        sage: L = Qp(5, prec = 200)
        sage: R.<x> = LaurentSeriesRing(K)
        sage: S.<y> = LaurentSeriesRing(L)
        sage: R is S
        False
        sage: T.<y> = LaurentSeriesRing(Qp(5,prec=200))
        sage: S is T
        True
        sage: W.<y> = LaurentSeriesRing(Qp(5,prec=199))
        sage: W is T
        False
    """
    if not names is None: name = names
    if name is None:
        raise TypeError, "You must specify the name of the indeterminate of the Laurent series ring."

    global laurent_series
    key = (base_ring, name, sparse)
    if laurent_series.has_key(key):
        x = laurent_series[key]()
        if x != None: return x

    if isinstance(base_ring, field.Field):
        R = LaurentSeriesRing_field(base_ring, name, sparse)
    elif isinstance(base_ring, integral_domain.IntegralDomain):
        R = LaurentSeriesRing_domain(base_ring, name, sparse)
    elif isinstance(base_ring, commutative_ring.CommutativeRing):
        R = LaurentSeriesRing_generic(base_ring, name, sparse)
    else:
        raise TypeError, "base_ring must be a commutative ring"
    laurent_series[key] = weakref.ref(R)
    return R

def is_LaurentSeriesRing(x):
    return isinstance(x, LaurentSeriesRing_generic)

class LaurentSeriesRing_generic(commutative_ring.CommutativeRing):
    """
    Univariate Laurent Series Ring
    EXAMPLES:
        sage: K, q = LaurentSeriesRing(CC, 'q').objgen(); K
        Laurent Series Ring in q over Complex Field with 53 bits of precision
        sage: loads(K.dumps()) == K
        True
    """

    def __init__(self, base_ring, name=None, sparse=False):
        ParentWithGens.__init__(self, base_ring, name)
        self.__sparse = sparse

    def base_extend(self, R):
        """
        Returns the laurent series ring over R in the same variable as
        self, assuming there is a canonical coerce map from the base
        ring of self to R.
        """
        if R.has_coerce_map_from(self.base_ring()):
            return self.change_ring(R)
        else:
            raise TypeError, "no valid base extension defined"

    def change_ring(self, R):
        return LaurentSeriesRing(R, self.variable_name(), sparse=self.__sparse)

    def is_sparse(self):
        return self.__sparse

    def is_field(self):
        return True

    def is_dense(self):
        return not self.__sparse

    def __reduce__(self):
        return self.__class__, (self.base_ring(), self.variable_name())

    def __repr__(self):
        s = "Laurent Series Ring in %s over %s"%(self.variable_name(), self.base_ring())
        if self.is_sparse():
            s = 'Sparse ' + s
        return s

    def __call__(self, x, n=0):
        r"""
        Coerces the element x into this Laurent series ring.

        INPUT:
            x -- the element to coerce
            n -- the result of the coercion will be multiplied by $t^n$ (default: 0)

        EXAMPLES:
            sage: R.<u> = LaurentSeriesRing(Qp(5, 10))
            sage: S.<t> = LaurentSeriesRing(RationalField())
            sage: print R(t + t^2 + O(t^3))
            (1 + O(5^10))*u + (1 + O(5^10))*u^2 + O(u^3)

        Note that coercing an element into its own parent just produces
        that element again (since Laurent series are immutable):

            sage: u is R(u)
            True

        Rational functions are accepted:

            sage: I = sqrt(-1)
            sage: K.<I> = QQ[I]
            sage: P.<t> = PolynomialRing(K)
            sage: L.<u> = LaurentSeriesRing(QQ[I])
            sage: L((t*I)/(t^3+I*2*t))
            1/2 + 1/4*I*u^2 - 1/8*u^4 - 1/16*I*u^6 + 1/32*u^8 +
            1/64*I*u^10 - 1/128*u^12 - 1/256*I*u^14 + 1/512*u^16 +
            1/1024*I*u^18 + O(u^20)

            sage: L(t*I) / L(t^3+I*2*t)
            1/2 + 1/4*I*u^2 - 1/8*u^4 - 1/16*I*u^6 + 1/32*u^8 +
            1/64*I*u^10 - 1/128*u^12 - 1/256*I*u^14 + 1/512*u^16 +
            1/1024*I*u^18 + O(u^20)

        """
        from sage.rings.fraction_field_element import is_FractionFieldElement
        from sage.rings.polynomial.polynomial_element import is_Polynomial
        from sage.rings.polynomial.multi_polynomial_element import is_MPolynomial

        if isinstance(x, laurent_series_ring_element.LaurentSeries) and n==0 and self is x.parent():
            return x  # ok, since Laurent series are immutable (no need to make a copy)
        elif is_FractionFieldElement(x) and \
             (x.base_ring() is self.base_ring() or x.base_ring() == self.base_ring()) and \
             (is_Polynomial(x.numerator()) or is_MPolynomial(x.numerator())):
            x = self(x.numerator())/self(x.denominator())
            return self.gen()**n * x
        else:
            return laurent_series_ring_element.LaurentSeries(self, x, n)

    def _coerce_impl(self, x):
        """
        Return canonical coercion of x into self.

        Rings that canonically coerce to this power series ring R:

           * R itself
           * Any laurent series ring in the same variable whose base ring canonically coerces to
             the base ring of R.
           * Any ring that canonically coerces to the power series ring over the base ring of R.
           * Any ring that canonically coerces to the base ring of R
        """
        try:
            P = x.parent()
            if is_LaurentSeriesRing(P):
                if P.variable_name() == self.variable_name():
                    if self.has_coerce_map_from(P.base_ring()):
                        return self(x)
                    else:
                        raise TypeError, "no natural map between bases of power series rings"
        except AttributeError:
            pass

        return self._coerce_try(x, [self.power_series_ring(), self.base_ring()])

    def __cmp__(self, other):
        if not isinstance(other, LaurentSeriesRing_generic):
            return cmp(type(self),type(other))
        c = cmp(self.base_ring(), other.base_ring())
        if c: return c
        c = cmp(self.variable_name(), other.variable_name())
        if c: return c
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

    def is_exact(self):
        return False

    def gen(self, n=0):
        if n != 0:
            raise IndexError, "Generator n not defined."
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
                                         self.base_ring(), self.variable_name(), sparse=self.is_sparse())
            return self.__power_series_ring

class LaurentSeriesRing_domain(LaurentSeriesRing_generic, integral_domain.IntegralDomain):
    def __init__(self, base_ring, name=None, sparse=False):
        LaurentSeriesRing_generic.__init__(self, base_ring, name, sparse)

class LaurentSeriesRing_field(LaurentSeriesRing_generic, field.Field):
    def __init__(self, base_ring, name=None, sparse=False):
        LaurentSeriesRing_generic.__init__(self, base_ring, name, sparse)

