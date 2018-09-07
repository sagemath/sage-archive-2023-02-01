from sage.structure.unique_representation import UniqueRepresentation
from sage.monoids.monoid import Monoid_class
from sage.rings.ring import CommutativeAlgebra
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ

from sage.categories.monoids import Monoids
from sage.categories.commutative_algebras import CommutativeAlgebras
from sage.categories.complete_discrete_valuation import CompleteDiscreteValuationRings
from sage.categories.complete_discrete_valuation import CompleteDiscreteValuationFields

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.tate_algebra_element import TateAlgebraTerm
from sage.rings.tate_algebra_element import TateAlgebraElement

from sage.categories.pushout import pushout

DEFAULT_CAP = 20

class TateTermMonoid(Monoid_class, UniqueRepresentation):
    def __init__(self, base, names, log_radii, order):
        # This function is not exposed to the user
        # so we do not check the inputs
        self.element_class = TateAlgebraTerm
        Monoid_class.__init__(self, names)
        self._base = base
        self._names = names
        self._ngens = len(self._names)
        self._log_radii = log_radii
        self._order = order

    def _repr_(self):
        if self._ngens == 0:
            return "Monoid of terms over %s" % self._base
        vars = ""
        for i in range(self._ngens):
            vars += ", %s (val >= %s)" % (self._names[i], -self._log_radii[i])
        return "Monoid of terms in %s over %s" % (vars[2:], self._base)

    def base_ring(self):
        return self._base

    def variable_names(self):
        return self._names

    def log_radii(self):
        return self._log_radii

    def term_order(self):
        return self._order

    def ngens(self):
        return self._ngens



class TateAlgebra(CommutativeAlgebra, UniqueRepresentation):
    def __init__(self, base, names, log_radii=QQ(0), prec=None, order='degrevlex'):
        r"""
        Create a Tate series ring over a given complete discrete valuation
        ring.

        Given a complete discrete valuation ring R, variables
        X1,...,Xk and weights r1,..., rn in RR_{>0}, the Tate algebra
        R{X1,...,Xk} is the algebra of power series with coefficients
        a_{i1,...,in} in R and such that
        |a_{i1,...,in}|*r1^i1*...*rn^in tends to 0 as i1,...,in go
        towards infinity.

        
        INPUT:

        - ``base`` - a complete discrete valuation ring or field

        - ``names`` - names of the indeterminates

        - ``log-radii`` - (default: ``0``) the value(s) -log(ri). If only
          one number l is given, all ri's are defined with -log(ri)=l.

        - ``prec`` - the default precision used if an exact object
          must be changed to an approximate object in order to do an
          arithmetic operation. If left as ``None``, it will be set to
          the default precision of the base ring, if any. Otherwise,
          it will be set to the default cap in precision of the base
          ring, if any. Otherwise, it will be set to the global
          default (20).

        - ``order`` - (default: ``degrevlex``) the monomial ordering 
          used to break ties when comparing terms with the same 
          coefficient valuation

        EXAMPLES::

            sage: R = Zp(2, 10, print_mode='digits'); R
            2-adic Ring with capped relative precision 10

        ::
        
            sage: A.<x,y> = TateAlgebra(R, order='lex'); A
            Tate Algebra in x, y over 2-adic Ring with capped relative precision 10

        The term ordering is used to determine how series are
        displayed. Terms are compared first according to the valuation
        of their coefficient, and ties are broken using the monomial
        ordering.

            sage: A.term_order()
            Lexicographic term order
            sage: f = 2 + y^5 + x^2; f
            (...0000000001)*x^2 + (...0000000001)*y^5 + (...00000000010)
            sage: B.<x,y> = TateAlgebra(R); B
            Tate Algebra in x, y over 2-adic Ring with capped relative precision 1
            sage: B.term_order()
            Degree reverse lexicographic term order
            sage: B(f)
            (...0000000001)*y^5 + (...0000000001)*x^2 + (...00000000010)

        """
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
            self._log_radii = (QQ(log_radii),) * self._ngens
        elif len(log_radii) != self._ngens:
            raise ValueError("The number of radii does not match the number of variables")
        else:
            self._log_radii = ( QQ(r) for r in log_radii )
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
        self._parent_terms = TateTermMonoid(self._base, self._names, self._log_radii, self._order)

    def _an_element_(self):
        return self.element_class(0)

    def _coerce_map_from_(self, R):
        base = self._base
        if base.has_coerce_map_from(R):
            return True
        if isinstance(R, (TateTermMonoid, TateAlgebra)):
            Rbase = R.base_ring()
            if base.has_coerce_map_from(Rbase) and self._names == R.variable_names() and self._order == R.term_order():
                ratio = base.absolute_e() // Rbase.absolute_e()
                for i in range(self._ngens):
                    if self._log_radii[i] != R.log_radii()[i] * ratio:
                        return False
                return True

    def _ideal_class_(self, n):
        from sage.rings.tate_algebra_ideal import TateAlgebraIdeal
        return TateAlgebraIdeal

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
        if self._ngens == 0:
            return "Tate Algebra over %s" % self._base
        vars = ""
        for i in range(self._ngens):
            vars += ", %s (val >= %s)" % (self._names[i], -self._log_radii[i])
        return "Tate Algebra in %s over %s" % (vars[2:], self._base)

    def variable_names(self):
        return self._names

    def log_radii(self):
        return self._log_radii

    def monoid_of_terms(self):
        return self._parent_terms

    def term_order(self):
        return self._order

    def precision_cap(self):
        return self._cap

    def characteristic(self):
        return self.base_ring().characteristic()

    #def ideal(self, gens):
    #    pass
