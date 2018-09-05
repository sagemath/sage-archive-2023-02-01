from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.ring import CommutativeAlgebra
from sage.rings.rational_ring import QQ

from sage.categories.commutative_algebras import CommutativeAlgebras

from sage.rings.tate_algebra_element import TateAlgebraElement


class TateAlgebra(CommutativeAlgebra, UniqueRepresentation):
    def __init__(self, base, names, radii=QQ(1), order=None):
        if base not in CompleteDiscreteValuationRings() and base not in CompleteDiscreteValuationFields():
            raise TypeError("The base ring must be a complete discrete valuation ring or field")
        CommutativeAlgebra.__init__(self, base, names, category=CommutativeAlgebras(base))
        if isinstance(names, (list, tuple)):
            names = [ str(var) for var in names ]
        else:
            names = [ str(names) ]
        nvars = len(names)
        if not isinstance(radii, (list, tuple)):
            self._radii = [ QQ(radii) ] * nvars
        elif len(radii) != nvars:
            raise ValueError("The length of the radii does not match the number of variables")
        else:
            self._radii = [ QQ(r) for r in radii ]
        self._polynomial_ring = PolynomialRing(base, names)
        self._names = self._polynomial_ring.variable_names()
        self._order = order
        self.element_class = TateAlgebraElement

    def _coerce_map_from_(self, R):
        if self.base_ring().has_coerce_map_from(R):
            return True
        # Define other coercion maps (from polynomial rings, from other Tate algebras...)

    def _repr_(self):
        vars = ", ".join(self._names)
        return "Tate Algebra in %s over %s" % (vars, self.base_ring())

    def ngens(self):
        return len(self._names)

    def variable_names(self):
        return self._names

    def radii(self):
        return self._radii

    def term_order(self):
        raise NotImplementedError
        return self._order

    def characteristic(self):
        return self.base_ring().characteristic()

    #def ideal(self, gens):
    #    pass
