from sage.structure.element import CommutativeAlgebraElement

cdef class TateAlgebraElement(CommutativeAlgebraElement):
    def __init__(self, parent, x, prec, start=0):
        pass

    def precision_absolute(self):
        pass

    def _add_(self, other):
        pass

    def _sub_(self, other):
        pass

    def _mul_(self, other):
        pass

    def inverse_of_unit(self):
        pass

    def valuation(self):
        pass

    def leading_term(self):
        pass

    def leading_coefficient(self):
        pass

    def leading_monomial(self):
        pass

    def weierstrass_degree(self):
        pass

    def degree(self):
        return weierstrass_degree(self)

    def weierstrass_degrees(self):
        pass

    def degrees(self):
        return weierstrass_degrees(self)

    def residue(self, n=1):
        pass

    def newton_polytope(self):
        pass

    def quo_rem(self, divisors):
        pass

    def reduce(self, I):
        pass

