from sage.structure.element cimport Element
from sage.structure.element cimport MonoidElement
from sage.structure.element cimport CommutativeAlgebraElement

from sage.rings.infinity import Infinity
from sage.rings.integer_ring import ZZ

from sage.structure.element import coerce_binop

from sage.rings.polynomial.polydict cimport ETuple
from sage.rings.padics.padic_generic_element cimport pAdicGenericElement


cdef class TateAlgebraTerm(MonoidElement):
    def __init__(self, parent, coeff, exponent):
        MonoidElement.__init__(self, parent)
        field = parent.base_ring().fraction_field()
        self._coeff = field(coeff)
        if isinstance(exponent, ETuple):
            self._exponent = exponent
        else:
            self._exponent = ETuple(exponent)
        if len(self._exponent) != parent.ngens():
            raise ValueError("The length of the exponent does not match the number of variables")

    cdef TateAlgebraTerm _new_c(self):
        cdef TateAlgebraTerm ans = TateAlgebraTerm.__new__(TateAlgebraTerm)
        ans._parent = self._parent
        return ans

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
        cdef TateAlgebraTerm ans = self._new_c()
        ans._exponent = self._exponent.eadd((<TateAlgebraTerm>other)._exponent)
        ans._coeff = self._coeff * (<TateAlgebraTerm>other)._coeff
        return ans

    cpdef int _cmp_(self, other) except -2:
        cdef int c
        c = cmp(-(<TateAlgebraTerm>self)._valuation_c(), 
                -(<TateAlgebraTerm>other)._valuation_c())
        if c: return c
        T = (<TateAlgebraTerm>self)._parent.term_order()
        c = cmp(T.sortkey((<TateAlgebraTerm>self)._exponent), 
                T.sortkey((<TateAlgebraTerm>other)._exponent))
        if c: return c
        return cmp((<TateAlgebraTerm>self)._coeff, 
                   (<TateAlgebraTerm>other)._coeff)

    def valuation(self):
        return ZZ(self._valuation_c())

    cdef long _valuation_c(self):
        return (<pAdicGenericElement>self._coeff).valuation_c() - <long>self._exponent.dotprod(self._parent._log_radii)

    @coerce_binop
    def is_coprime_with(self, other):
        for i in range(self._parent.ngens()):
            if self._exponent[i] > 0 and other.exponent()[i] > 0:
                return False
        if self._parent.base_ring().is_field():
            return True
        else:
            return self.valuation() == 0 or other.valuation() == 0

    @coerce_binop
    def gcd(self, other):
        return self._gcd_c(other)

    cdef TateAlgebraTerm _gcd_c(self, TateAlgebraTerm other):
        cdef TateAlgebraTerm ans = self._new_c()
        ans._exponent = self._exponent.emin(other._exponent)
        if self._coeff.valuation() < other._coeff.valuation():
            ans._coeff = self._coeff
        else:
            ans._coeff = other._coeff
        return ans

    @coerce_binop
    def lcm(self, other):
        return self._lcm_c(other)

    cdef TateAlgebraTerm _lcm_c(self, TateAlgebraTerm other):
        cdef TateAlgebraTerm ans = self._new_c()
        ans._exponent = self._exponent.emax(other._exponent)
        if self._coeff.valuation() < other._coeff.valuation():
            ans._coeff = other._coeff
        else:
            ans._coeff = self._coeff
        return ans

    @coerce_binop
    def is_divisible_by(self, other, integral=False):
        return (<TateAlgebraTerm?>other)._divides_c(self, integral)
    @coerce_binop
    def divides(self, other, integral=False):
        return self._divides_c(other, integral)

    cdef bint _divides_c(self, TateAlgebraTerm other, bint integral):
        parent = self._parent
        if (integral or not parent.base_ring().is_field()) and self.valuation() > other.valuation():
            return False
        for i in range(parent._ngens):
            if self._exponent[i] > other._exponent[i]:
                return False
        return True

    @coerce_binop
    def __floordiv__(self, other):
        parent = self.parent()
        if not parent.base_ring().is_field() and self.valuation() < other.valuation():
            raise ValueError("The division is not exact")
        for i in range(parent._ngens):
            if self.exponent()[i] < other.exponent()[i]:
                raise ValueError("The division is not exact")
        return (<TateAlgebraTerm>self)._floordiv_c(<TateAlgebraTerm>other)

    cdef TateAlgebraTerm _floordiv_c(self, TateAlgebraTerm other):
        cdef TateAlgebraTerm ans = self._new_c()
        ans._exponent = self.exponent().esub(other.exponent())
        ans._coeff = self._coeff / other._coeff
        return ans


cdef class TateAlgebraElement(CommutativeAlgebraElement):
    def __init__(self, parent, x, prec=None):
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
        #self._normalize()
        self._terms = None

    cdef TateAlgebraElement _new_c(self):
        cdef TateAlgebraElement ans = TateAlgebraElement.__new__(TateAlgebraElement)
        ans._parent = self._parent
        return ans

    cdef _normalize(self):
        if self._prec is Infinity: return
        cdef TateAlgebraTerm t
        cdef dict coeffs = self._poly.dict()
        cdef int v
        for (e,c) in coeffs.iteritems():
            v = (<ETuple>self._parent._log_radii).dotprod(<ETuple>e)
            coeffs[e] = coeffs[e].add_bigoh((self._prec - v).ceil())
        self._poly = self._parent._polynomial_ring(coeffs)

    def _repr_(self):
        self._normalize()
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
        cdef TateAlgebraElement ans = self._new_c()
        ans._poly = self._poly + (<TateAlgebraElement>other)._poly
        ans._prec = min(self._prec, (<TateAlgebraElement>other)._prec)
        if self._prec != (<TateAlgebraElement>other)._prec:
            ans._normalize()
        return ans

    cpdef _neg_(self):
        cdef TateAlgebraElement ans = self._new_c()
        ans._poly = -self._poly
        ans._prec = self._prec
        return ans

    cpdef _sub_(self, other):
        cdef TateAlgebraElement ans = self._new_c()
        ans._poly = self._poly - (<TateAlgebraElement>other)._poly
        ans._prec = min(self._prec, (<TateAlgebraElement>other)._prec)
        if self._prec != (<TateAlgebraElement>other)._prec:
            ans._normalize()
        return ans

    cpdef _mul_(self, other):
        cdef TateAlgebraElement ans = self._new_c()
        a = self._prec + (<TateAlgebraElement>other).valuation()
        b = self.valuation() + (<TateAlgebraElement>other)._prec
        ans._poly = self._poly * (<TateAlgebraElement>other)._poly
        ans._prec = min(a, b)
        #ans._normalize()
        return ans

    cpdef _lmul_(self, Element right):
        cdef TateAlgebraElement ans = self._new_c()
        ans._poly = self._poly * right
        ans._prec = self._prec + (<pAdicGenericElement>self._parent._base(right)).valuation_c()
        return ans

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
            prec = self._prec
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

    def terms(self):
        self._normalize()
        self._terms = None
        return self._terms_c()

    cdef list _terms_c(self):
        if self._terms is not None:
            return self._terms
        cdef TateAlgebraTerm oneterm = self._parent._oneterm
        cdef TateAlgebraTerm term
        self._terms = []
        for (e,c) in self._poly.dict().iteritems():
            term = oneterm._new_c()
            term._coeff = c
            term._exponent = e
            if term.valuation() < self._prec:
                self._terms.append(term)
        self._terms.sort(reverse=True)
        return self._terms

    def coefficients(self):
        return [ t.coefficient() for t in self.terms() ]

    def add_bigoh(self, n):
        return self._parent(self, prec=n)

    def precision_absolute(self):
        return self._prec

    cpdef valuation(self):
        cdef TateAlgebraTerm t
        cdef list terms = self._terms_c()
        if terms:
            return min(terms[0].valuation(), self._prec)
        else:
            return self._prec

    def precision_relative(self):
        return self._prec - self.valuation()

    def leading_term(self):
        terms = self.terms()
        if terms:
            return terms[0]

    def leading_coefficient(self):
        terms = self.terms()
        if terms:
            return terms[0].coefficient()
        else:
            return self.base_ring()(0)

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
        try:
            return self._Spoly_c(other)
        except IndexError:
            raise ValueError("Cannot compute the S-polynomial of zero")

    cdef TateAlgebraElement _Spoly_c(self, TateAlgebraElement other):
        cdef TateAlgebraTerm st = self._terms_c()[0]
        cdef TateAlgebraTerm ot = other._terms_c()[0]
        cdef TateAlgebraTerm t = st._lcm_c(ot)
        return (t._floordiv_c(st))*self - (t._floordiv_c(ot))*other
