from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.ring import CommutativeAlgebra
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ

from sage.categories.commutative_algebras import CommutativeAlgebras
from sage.categories.complete_discrete_valuation import CompleteDiscreteValuationRings
from sage.categories.complete_discrete_valuation import CompleteDiscreteValuationFields

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.tate_algebra_element import TateAlgebraElement

from sage.categories.pushout import pushout

DEFAULT_CAP = 20

class TateAlgebra(CommutativeAlgebra, UniqueRepresentation):
    def __init__(self, base, names, log_radii=QQ(0), prec=None, order='degrevlex'):
        if base not in CompleteDiscreteValuationRings() and base not in CompleteDiscreteValuationFields():
            raise TypeError("The base ring must be a complete discrete valuation ring or field")
        if isinstance(names, (list, tuple)):
            names = [ str(var) for var in names ]
        else:
            names = [ str(names) ]
        self._ngens = len(names)
        self.element_class = TateAlgebraElement
        CommutativeAlgebra.__init__(self, base, names, category=CommutativeAlgebras(base))
        if not isinstance(log_radii, (list, tuple)):
            self._log_radii = [ QQ(log_radii) ] * self._ngens
        elif len(log_radii) != self._ngens:
            raise ValueError("The number of radii does not match the number of variables")
        else:
            self._log_radii = [ QQ(r) for r in log_radii ]
        field = base.fraction_field()
        self._polynomial_ring = PolynomialRing(field, names, order=order)
        self._names = self._polynomial_ring.variable_names()
        self._order = self._polynomial_ring.term_order()
        self._gens = tuple([ self((field(1) << self._log_radii[i].ceil()) * self._polynomial_ring.gen(i)) for i in range(self._ngens) ])
        if prec is None:
            try:
                self._cap = base.precision_cap()
            except AttributeError:
                try:
                    self._cap = base.default_prec()
                except AttributeError:
                    self._cap = DEFAULT_CAP
        else:
            self._cap = ZZ(prec)

    def _an_element_(self):
        return self.element_class(0)

    def _coerce_map_from_(self, R):
        base = self._base
        if base.has_coerce_map_from(R):
            return True
        if isinstance(R, TateAlgebra):
            Rbase = R.base_ring()
            if base.has_coerce_map_from(Rbase) and self._names == R.variable_names() and self._order == R.term_order():
                ratio = base.absolute_e() // Rbase.absolute_e()
                for i in range(self._ngens):
                    if self._log_radii[i] != R.log_radii()[i] * ratio:
                        return False
                return True

    #def _pushout_(self, R):
    #    if not isinstance(R, TateAlgebra):  # should we allow PolynomialRing as well?
    #        return None
    #    if self._names != R.variable_names():
    #        return None
    #    if self._order != R.term_order():
    #        return None
    #    Sbase = self._base
    #    Rbase = R.base_ring()
    #    base = pushout(Sbase, Rbase)
    #    Sratio = base.absolute_e() // Sbase.absolute_e()
    #    Rratio = base.absolute_e() // Rbase.absolute_e()
    #    log_radii = tuple([ min(self._log_radii[i] * Sratio, R.log_radii()[i] * Rratio) for i in range(self._ngens) ])
    #    cap = min(self._cap * Sratio, R.precision_cap() * Rratio)
    #    return TateAlgebra(base, self._names, log_radii, cap, self._order)

    def gen(self, n=0):
        return self._gens[n]

    def gens(self):
        return self._gens

    def ngens(self):
        return self._ngens

    def _repr_(self):
        vars = ", ".join(self._names)
        return "Tate Algebra in %s over %s" % (vars, self.base_ring())

    def variable_names(self):
        return self._names

    def log_radii(self):
        return self._log_radii

    def term_order(self):
        return self._order

    def precision_cap(self):
        return self._cap

    def characteristic(self):
        return self.base_ring().characteristic()

    #def ideal(self, gens):
    #    pass
