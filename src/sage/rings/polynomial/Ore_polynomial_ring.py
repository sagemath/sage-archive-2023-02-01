from __future__ import print_function, absolute_import, division
    
from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.ring import Algebra
from sage.structure.element import Element
from sage.categories.rings import Rings
from sage.rings.polynomial.Ore_polynomial_element import OrePolynomial1

class OrePolynomialRing_general(Algebra, UniqueRepresentation):

    def __init__(self, base_ring, twist_map, derivation, name, sparse = False):
        r"""

        INPUT :
        - ''base_ring'' -- a commutative ring
        - ''twist_map'' -- an automorphism of the base ring
        - ''derivation'' -- a derivation of the base ring
        - ''name'' --string or list od strings representing the name of
        the variables of ring
        - ''sparse'' -- boolean (defaut : ''False'')
        - ''element_class'' -- class representing the type of element to be used in ring
        """
        if sparse:
            raise NotImplementedError()
        self._is_sparse = sparse
        self.Element = OrePolynomial1
        self._map = twist_map
        self._derivation = derivation
        Algebra.__init__(self, base_ring, names = name, normalize = True, category = Rings())


    def _coerce_map_from_(self, R):
        if self.base_ring().has_coerce_map_from(R):
            return True
        
    
    def _repr_(self):
        s = "Ore Polynomial Ring in %s over %s twisted by %s and derivate by %s" % (self.variable_name(),
                                                                                    self.base_ring(),
                                                                                    self._map._repr_short(),
                                                                                    self._derivation._repr_())
        if self.is_sparse():
            s = "Sparse " + s
        return s

    def characteristic(self):
        return self.base_ring().characteristic()

    def twist_map(self, n=1):
        try:
            return self._map ** n
        except TypeError as e:
            if n<0:
                raise NotImplementedError("inversion of the twist map %s" % self._map)
            else:
                raise ValueError("Unexpected error in iterating the twist map : %s", e)

    def gen(self, n=0):
        if n != 0:
            raise IndexError("generator %s not defined" %n)
        return self([0,1])

    parameter = gen

    def gens_dict(self):
        return dict(zip(self.variable_names(), self.gens()))

    def is_finite(self):
        R = self.base_ring()
        return R.is_finite() and R.order() == 1

    def is_exact(self):
        return self.base_ring().is_exact()
    
    def is_sparse(self):
        return self._is_sparse

    def ngens(self):
        return 1

    def is_commutative(self):
        return self.twist_map().is_identity() and self.derivation.is_zero()
        
    def zero(self):
        return OrePolynomial1(self, [0])
    def base_ring_element(self):
        return OrePolynomial1(self, [0, 1])
    def base_multa(self, a):
        return OrePolynomial1(self, [self._derivation(a), self._map(a)])

    def random_element(self, k):
        return OrePolynomial1(self, [self.base_ring().random_element() for _ in range(k)])
