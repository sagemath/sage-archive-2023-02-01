from sage.structure.element cimport CommutativeAlgebraElement
from sage.rings.infinity import Infinity

from sage.rings.integer_ring import ZZ

cdef class TateAlgebraElement(CommutativeAlgebraElement):
    def __init__(self, parent, x, prec=None, reduce=True):
        CommutativeAlgebraElement.__init__(self, parent)
        self._prec = Infinity
        self._start = Infinity
        S = parent._polynomial_ring
        if isinstance(x, TateAlgebraElement):
            if x.parent() is parent:
                self._poly = x._poly
                self._prec = x._prec
            else:
                self._poly = S(x._poly)
                self._prec = x._prec * S.absolute_e() // x.base_ring().absolute_e()
        else:
            self._poly = S(x)
        if prec is not None:
            self._prec = min(self._prec, prec)
        if reduce:
            self._poly = self._poly.map_coefficients(lambda c: c.add_bigoh(self._prec))

    def polynomial(self):
        return self._poly

    def _add_(self, other):
        return self.__class__()(self._poly + other.polynomial(), min(self._prec, other.precision_absolute()))

    def _neg_(self, other):
        return self.__class__()(-self._poly, self._prec, reduce=False)

    def _sub_(self, other):
        return self.__class__()(self._poly - other.polynomial(), min(self._prec, other.precision_absolute()))

    def _mul_(self, other):
        prec = min(self._prec + other.valuation(), other.precision_absolute() + self.valuation())
        return self.__class__()(self._poly * other.polynomial(), prec)

    def is_zero(self):
        # This assumes that self._poly is reduced
        return self._poly == 0

    def inverse_of_unit(self):
        pass

    @cached_method
    def _lt(self):
        if self.is_zero():
            return ((ZZ(0),) * self.ngens(), self.base_ring()(0))
        terms = list(self._poly.dict().iteritems())
        T = self._parent._order
        sortkey = lambda c: (-c[1].valuation(), T.sortkey(c[0]))
        return max(terms, key=sortkey)
  
    def leading_term(self):
        lt = self._lt()
        poly = self._parent._polynomial_ring.monomial(*lt[0])
        return lt[1] * self.__class__(poly)

    def leading_coefficient(self):
        lt = self._lt()
        return lt[1]
        
    def leading_monomial(self):
        lt = self._lt()
        return self._parent._polynomial_ring.monomial(*lt[0])

    def precision_absolute(self):
        return self._prec

    def valuation(self):
        return self.leading_coefficient().valuation()

    def precision_relative(self):
        return self._prec - self.valuation()

    def weierstrass_degree(self):
        v = self.valuation()
        return self.residue(v+1).degree()

    def degree(self):
        return self.weierstrass_degree()

    def weierstrass_degrees(self):
        v = self.valuation()
        return self.residue(v+1).degrees()

    def degrees(self):
        return self.weierstrass_degrees()

    def residue(self, n=1):
        try:
            Rn = self.base_ring().residue_ring(n)
        except (AttributeError, NotImplementedError):
            Rn = self.base_ring().change(type="fixed-mod", prec=n)
        return self._poly.change_ring(Rn)

    def quo_rem(self, *divisors):
        parent = self._parent
        nvars = parent.ngens()
        q = [ parent(0) ] * len(divisors)
        r = parent(0)
        ltds = [ d._lt() for d in divisors ]
        # TODO: do everything on self._poly
        while not self.is_zero():
            lt = self._lt()
            found = False
            while not found and i < len(divisors):
                ltd = ltds[i]
                is_divisible = ltd[1].valuation() <= lt[1].valuation()
                j = 0
                while is_divisible and j < nvars:
                    is_divisible = ltd[0][j] <= lt[0][j]
                    j += 1
                if is_divisible:
                    found = True
                    break
                i += 1
            if found:
                exponent = ( lt[0][j] - ltd[0][j] for j in range(nvars) )
                factor = (lt[1]//ltd[1]) * parent._polynomial_ring.monomial(*exponent)
                self -= factor * divisors[i]
                q[i] += factor
            else:
                self -= self.leading_term()
                r += self.leading_term()
        return q, r

    def __mod__(self, divisors):
        if not isinstance(divisors, (list, tuple)):
            divisors = [ divisors ]
        _, r = self.quo_rem(*divisors)
        return r

    def __floordiv__(self, divisors):
        if not isinstance(divisors, (list, tuple)):
            q, _ = self.quo_rem(divisors)
            return q[0]
        else:
            q, _ = self.quo_rem(*divisors)
            return q

    def reduce(self, I):
        g = I.groebner_basis()
        return self.quo_rem(g)
