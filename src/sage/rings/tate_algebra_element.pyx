from sage.structure.element cimport Element
from sage.structure.element cimport MonoidElement
from sage.structure.element cimport CommutativeAlgebraElement

from sage.rings.infinity import Infinity
from sage.rings.integer_ring import ZZ

from sage.structure.element import coerce_binop

from sage.rings.polynomial.polydict cimport ETuple


cdef class TateAlgebraTerm(MonoidElement):
    def __init__(self, parent, coeff, exponent):
        MonoidElement.__init__(self, parent)
        self._field = parent.base_ring().fraction_field()
        self._coeff = self._field(coeff)
        if isinstance(exponent, ETuple):
            self._exponent = exponent
        else:
            self._exponent = ETuple(exponent)
        if len(self._exponent) != parent.ngens():
            raise ValueError("The length of the exponent does not match the number of variables")

    def _repr_(self):
        parent = self._parent
        s = "(%s)" % self._coeff
        for i in range(parent._ngens):
            if self._exponent[i] == 1:
                s += "*%s" % parent._names[i]
            elif self._exponent[i] > 1:
                s += "*%s^%s" % (parent._names[i], self._exponent[i])
        return s

    def coefficient(self):
        return self._coeff

    def exponent(self):
        return self._exponent

    cpdef _mul_(self, other):
        exponent = self._exponent.eadd(other.exponent())
        return self._parent(self._coeff * other.coefficient(), exponent)

    cpdef int _cmp_(self, other) except -2:
        c = cmp(-self.valuation(), -other.valuation())
        if c: return c
        T = self._parent.term_order()
        return cmp(T.sortkey(self._exponent), T.sortkey(other.exponent()))

    def valuation(self):
        parent = self._parent
        return self._coeff.valuation() - sum(self._exponent[i]*parent._log_radii[i] for i in range(parent._ngens))

    @coerce_binop
    def gcd(self, other):
        exponent = self._exponent.emin(other.exponent())
        val = min(self._coeff.valuation(), other.coefficient().valuation())
        return self._parent(self._field.uniformizer_pow(val), exponent)

    @coerce_binop
    def lcm(self, other):
        exponent = self._exponent.emax(other.exponent())
        val = max(self._coeff.valuation(), other.coefficient().valuation())
        return self._parent(self._field.uniformizer_pow(val), exponent)

    @coerce_binop
    def is_divisible_by(self, other, integral=False):
        parent = self._parent
        if (integral or not parent.base_ring().is_field()) and self.valuation() < other.valuation():
            return False
        for i in range(parent._ngens):
            if self._exponent[i] < other.exponent()[i]:
                return False
        return True

    def divides(self, other):
        return other.is_divisible_by(self)

    @coerce_binop
    def __floordiv__(self, other):
        parent = self.parent()
        if not parent.base_ring().is_field() and self.valuation() < other.valuation():
            raise ValueError("The division is not exact")
        exponent = [ ]
        for i in range(parent._ngens):
            if self.exponent()[i] < other.exponent()[i]:
                raise ValueError("The division is not exact")
        exponent = self.exponent().esub(other.exponent())
        return parent(self.coefficient() / other.coefficient(), ETuple(exponent))


cdef class TateAlgebraElement(CommutativeAlgebraElement):
    def __init__(self, parent, x, prec=None, reduce=True):
        cdef TateAlgebraElement xc
        CommutativeAlgebraElement.__init__(self, parent)
        self._prec = Infinity
        S = parent._polynomial_ring
        if isinstance(x, TateAlgebraElement):
            xc = <TateAlgebraElement>x
            xparent = x.parent()
            if xparent is parent:
                self._poly = xc._poly
                self._prec = xc._prec
            elif xparent.variable_names() == parent.variable_names():
                ratio = parent._base.absolute_e() / xparent.base_ring().absolute_e()
                for i in range(parent.ngens()):
                    if parent.log_radii()[i] > xparent.log_radii()[i] * ratio:
                        raise ValueError("Cannot restrict to a bigger domain")
                self._poly = S(xc._poly)
                if xc._prec is not Infinity:
                    self._prec = (xc._prec * ratio).ceil()
            else:
                raise TypeError("Variable names do not match")
        elif isinstance(x, TateAlgebraTerm):
            xparent = x.parent()
            if xparent.variable_names() == parent.variable_names():
                ratio = parent._base.absolute_e() / xparent.base_ring().absolute_e()
                for i in range(parent.ngens()):
                    if parent.log_radii()[i] > xparent.log_radii()[i] * ratio:
                        raise ValueError("Cannot restrict to a bigger domain")
                self._poly = x.coefficient() * S.monomial(*x.exponent())
        else:
            self._poly = S(x)
        if prec is not None:
            self._prec = min(self._prec, prec)
        if reduce:
            self._poly = self._poly.map_coefficients(lambda c: c.add_bigoh(self._prec))
        self._terms = None

    def _repr_(self):
        base = self._parent.base_ring()
        nvars = self._parent.ngens()
        vars = self._parent.variable_names()
        s = ""
        for t in self.terms():
            if t.valuation() >= self._prec:
                continue
            s += " + %s" % t
        if self._prec is not Infinity:
            # Only works with padics...
            s += " + %s" % base.change(show_prec="bigoh")(0, self._prec)
        if s == "":
            return "0"
        return s[3:]

    def polynomial(self):
        return self._poly

    cpdef _add_(self, other):
        return self._parent(self._poly + other.polynomial(), min(self._prec, other.precision_absolute()))

    cpdef _neg_(self):
        return self._parent(-self._poly, self._prec, reduce=False)

    cpdef _sub_(self, other):
        return self._parent(self._poly - other.polynomial(), min(self._prec, other.precision_absolute()))

    cpdef _mul_(self, other):
        prec = min(self._prec + other.valuation(), other.precision_absolute() + self.valuation())
        return self._parent(self._poly * other.polynomial(), prec)

    cpdef _lmul_(self, Element right):
        other = self._parent._base(right)
        prec = self._prec + other.valuation()
        return self._parent(self._poly * other, prec)

    def _lshift(self, n):
        parent = self._parent
        base = parent.base_ring()
        if base.is_field():
            poly = self._poly.map_coefficients(lambda c: c << n)
            return parent(poly, self._prec + n, reduce=False)
        else:
            coeffs = { }
            field = base.fraction_field()
            ngens = parent.ngens()
            for exp, c in self._poly.dict().iteritems():
                minval = sum([ exp[i] * parent.log_radii()[i] for i in range(parent.ngens()) ]).ceil()
                coeffs[exp] = field(base(c) >> (minval-n)) << minval
            return parent(coeffs, self._prec + n, reduce=False)

    def __lshift__(self, n):
        return self._lshift(n)

    def __rshift__(self, n):
        return self._lshift(-n)

    def is_zero(self, prec=None):
        if prec is None:
            # This assumes that self._poly is reduced
            return self._poly == 0
        else:
            return self.valuation() >= prec

    def inverse_of_unit(self):
        if not self.is_unit():
            raise ValueError("This series in not invertible")
        t = self.leading_term()
        c = t.coefficient()
        parent = self._parent
        v, c = c.val_unit()
        cap = parent.precision_cap()
        inv = parent(c.inverse_of_unit(), prec=cap)
        x = self.add_bigoh(cap)
        prec = 1
        while prec < cap:
            prec *= 2
            inv = 2*inv - self*inv*inv
        return inv << v

    def is_unit(self):
        if self.is_zero():
            return False
        t = self.leading_term()
        if max(t.exponent()) != 0:
            return False
        base = self.base_ring()
        if not base.is_field() and t.valuation() > 0:
            return False
        return True

    def restriction(self, log_radii):
        parent = self._parent
        from sage.rings.tate_algebra import TateAlgebra
        ring = TateAlgebra(self.base_ring(), parent.variable_names(), log_radii, parent.precision_cap(), parent.term_order())
        return ring(self)

    def dict(self):
        return self._poly.dict()

    def terms(self):
        if self._terms is not None:
            return self._terms
        if self.is_zero():
            self._terms = []
        else:
            parent = self._parent
            parent_terms = parent.monoid_of_terms()
            self._terms = [ parent_terms(c,e) for (e,c) in self.dict().iteritems() ]
            self._terms.sort(reverse=True)
        return self._terms

    def coefficients(self):
        return [ t.coefficient() for t in self.terms() ]

    def add_bigoh(self, n):
        return self._parent(self, prec=n)

    def precision_absolute(self):
        return self._prec

    def valuation(self):
        cap = self._parent.precision_cap()
        if self.is_zero():
            return cap
        else:
            return min(self.leading_term().valuation(), cap)

    def precision_relative(self):
        return self._prec - self.valuation()

    def leading_term(self):
        if self.is_zero():
            return None
        return self.terms()[0]

    def leading_coefficient(self):
        if self.is_zero():
            return self.base_ring(0)
        return self.leading_term().coefficient()

    def monic(self):
        c = self.leading_coefficient()
        poly = ~c * self._poly
        return self._parent(poly, self._prec - c.valuation())
        
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
        for r in self._parent.log_radii():
            if r != 0:
                raise NotImplementedError("Residues are only implemented for radius 1")
        try:
            Rn = self.base_ring().residue_ring(n)
        except (AttributeError, NotImplementedError):
            Rn = self.base_ring().change(field=False, type="fixed-mod", prec=n)
        return self._poly.change_ring(Rn)

    def quo_rem(self, divisors, integral=False):
        parent = self._parent
        nvars = parent.ngens()
        if not isinstance(divisors, (list, tuple)):
            divisors = [ divisors ]
        ltds = [ d.leading_term() for d in divisors ]
        # TODO: do everything on self._poly
        cap = parent.precision_cap()
        f = self.add_bigoh(cap)
        q = [ parent(0, cap - d.valuation()) for d in divisors ]
        r = parent(0, cap)
        while not f.is_zero():
            lt = f.leading_term()
            for i in range(len(divisors)):
                if lt.is_divisible_by(ltds[i], integral=integral):
                    factor = lt // ltds[i]
                    f -= factor * divisors[i]
                    q[i] += factor
                    break
            else:
                f -= lt; r += lt
        return q, r

    def __mod__(self, divisors):
        if not isinstance(divisors, (list, tuple)):
            divisors = [ divisors ]
        _, r = self.quo_rem(divisors)
        return r

    def __floordiv__(self, divisors):
        q, _ = self.quo_rem(divisors)
        if not isinstance(divisors, (list, tuple)):
            return q[0]
        else:
            q, _ = self.quo_rem(divisors)
            return q

    def reduce(self, I):
        g = I.groebner_basis()
        return self.quo_rem(g)

    @coerce_binop
    def Spoly(self, other):
        if self.is_zero() or other.is_zero():
            raise ValueError("Cannot compute the S-polynomial of zero")
        st = self.leading_term()
        ot = other.leading_term()
        t = st.lcm(ot)
        return (t//st)*self - (t//ot)*other
