from discrete_valuation import DiscreteValuation

from sage.rings.all import ZZ, infinity

class RationalFunctionFieldValuation(DiscreteValuation):
    def __init__(self, domain, base_valuation):
        if base_valuation.domain() is not domain._ring:
            raise ValueError

        self._base_valuation = base_valuation
        DiscreteValuation.__init__(self, domain)

    def _call_(self, x):
        return self._base_valuation(x.numerator()) - self._base_valuation(x.denominator())

    def __hash__(self):
        return hash(self._base_valuation) + hash(self.domain())

    def _cache_key(self):
        return self._base_valuation, self.domain()

    def __cmp__(self, other):
        if type(self) != type(other):
            return cmp(type(self),type(other))
        return cmp(self._base_valuation,other._base_valuation)

    def shift(self, f, v):
        if f.parent() is not self.domain():
            raise ValueError

        if v == 0:
            return f
        elif v < 0:
            return f/self._base_valuation.element_with_valuation(-v)
        else:
            return f*self._base_valuation.element_with_valuation(v)

    def value_group(self):
        return self._base_valuation.value_group()

    def residue_field(self):
        from sage.rings.all import FunctionField
        return FunctionField(self._base_valuation.residue_field(), names=self.domain().variable_names())

    def reduce(self, f):
        if f.parent() is not self.domain():
            raise ValueError
        if self(f) > 0:
            return self.residue_field().zero()
        if self(f) < 0:
            raise ValueError

        base = self._base_valuation

        num = f.numerator()
        den = f.denominator()

        assert base(num) == base(den)
        shift = base.element_with_valuation(-base(num))
        num *= shift
        den *= shift
        ret = base.reduce(num) / base.reduce(den)
        assert not ret.is_zero()
        return self.residue_field()(ret)

    def lift(self, F):
        if F.parent() is not self.residue_field():
            raise ValueError

        return self.domain()(self._base_valuation.lift(F.numerator()) / self._base_valuation.lift(F.denominator()))

    def mac_lane_approximants(self, G, precision_cap=None, assume_squarefree=False):
        """
        Some difficult cases from Mark van Hoeij::

        sage: from sage.rings.padics.function_field_valuation import RationalFunctionFieldValuation, TrivialValuation
        sage: k = GF(2)
        sage: K.<x> = FunctionField(k)
        sage: R.<y> = K[]
        sage: F = y^21 + x*y^20 + (x^3 + x + 1)*y^18 + (x^3 + 1)*y^17 + (x^4 + x)*y^16 + (x^7 + x^6 + x^3 + x + 1)*y^15 + x^7*y^14 + (x^8 + x^7 + x^6 + x^4 + x^3 + 1)*y^13 + (x^9 + x^8 + x^4 + 1)*y^12 + (x^11 + x^9 + x^8 + x^5 + x^4 + x^3 + x^2)*y^11 + (x^12 + x^9 + x^8 + x^7 + x^5 + x^3 + x + 1)*y^10 + (x^14 + x^13 + x^10 + x^9 + x^8 + x^7 + x^6 + x^3 + x^2 + 1)*y^9 + (x^13 + x^9 + x^8 + x^6 + x^4 + x^3 + x)*y^8 + (x^16 + x^15 + x^13 + x^12 + x^11 + x^7 + x^3 + x)*y^7 + (x^17 + x^16 + x^13 + x^9 + x^8 + x)*y^6 + (x^17 + x^16 + x^12 + x^7 + x^5 + x^2 + x + 1)*y^5 + (x^19 + x^16 + x^15 + x^12 + x^6 + x^5 + x^3 + 1)*y^4 + (x^18 + x^15 + x^12 + x^10 + x^9 + x^7 + x^4 + x)*y^3 + (x^22 + x^21 + x^20 + x^18 + x^13 + x^12 + x^9 + x^8 + x^7 + x^5 + x^4 + x^3)*y^2 + (x^23 + x^22 + x^20 + x^17 + x^15 + x^14 + x^12 + x^9)*y + x^25 + x^23 + x^19 + x^17 + x^15 + x^13 + x^11 + x^5
        sage: x = K._ring.gen()
        sage: v0 = RationalFunctionFieldValuation(K, GaussValuation(K._ring, TrivialValuation(k)).extension(x,1))
        sage: v0.mac_lane_approximants(F, assume_squarefree=True) # assume_squarefree for speed, not tested - factorization over function fields over finite fields is missing
        [[ Gauss valuation induced by Valuation on rational function field induced by [ Gauss valuation induced by Trivial valuation, v(x) = 1 ], v(y + x + 1) = 3/2 ],
         [ Gauss valuation induced by Valuation on rational function field induced by [ Gauss valuation induced by Trivial valuation, v(x) = 1 ], v(y) = 4/3, v(y^3 + x^4) = 13/3 ],
         [ Gauss valuation induced by Valuation on rational function field induced by [ Gauss valuation induced by Trivial valuation, v(x) = 1 ], v(y + x) = 2 ],
         [ Gauss valuation induced by Valuation on rational function field induced by [ Gauss valuation induced by Trivial valuation, v(x) = 1 ], v(y^15 + y^13 + (x + 1)*y^12 + x*y^11 + (x + 1)*y^10 + y^9 + y^8 + x*y^6 + x*y^5 + y^4 + y^3 + y^2 + (x + 1)*y + x + 1) = 2 ]]
        sage: v0 = RationalFunctionFieldValuation(K, GaussValuation(K._ring, TrivialValuation(k)).extension(x+1,1))
        sage: v0.mac_lane_approximants(F, assume_squarefree=True) # not tested, factorization over function fields over finite fields is missing
        [[ Gauss valuation induced by Valuation on rational function field induced by [ Gauss valuation induced by Trivial valuation, v(x + 1) = 1 ], v(y) = 7/2, v(y^2 + x^7 + x^6 + x^5 + x^4 + x^3 + x^2 + x + 1) = 15/2 ],
         [ Gauss valuation induced by Valuation on rational function field induced by [ Gauss valuation induced by Trivial valuation, v(x + 1) = 1 ], v(y + x^2 + 1) = 7/2 ],
         [ Gauss valuation induced by Valuation on rational function field induced by [ Gauss valuation induced by Trivial valuation, v(x + 1) = 1 ], v(y) = 3/4, v(y^4 + x^3 + x^2 + x + 1) = 15/4 ],
         [ Gauss valuation induced by Valuation on rational function field induced by [ Gauss valuation induced by Trivial valuation, v(x + 1) = 1 ], v(y^13 + x*y^12 + y^10 + (x + 1)*y^9 + (x + 1)*y^8 + x*y^7 + x*y^6 + (x + 1)*y^4 + y^3 + (x + 1)*y^2 + 1) = 2 ]]
        sage: v0 = RationalFunctionFieldValuation(K, GaussValuation(K._ring, TrivialValuation(k)).extension(x^3+x^2+1,1))
        sage: v0.mac_lane_approximants(F, assume_squarefree=True) # not tested, factorization over function fields over finite fields is missing
        [[ Gauss valuation induced by Valuation on rational function field induced by [ Gauss valuation induced by Trivial valuation, v(x^3 + x^2 + 1) = 1 ], v(y + x^3 + x^2 + x) = 2, v(y^2 + (x^6 + x^4 + 1)*y + x^14 + x^10 + x^9 + x^8 + x^5 + x^4 + x^3 + x^2 + x) = 5 ],
         [ Gauss valuation induced by Valuation on rational function field induced by [ Gauss valuation induced by Trivial valuation, v(x^3 + x^2 + 1) = 1 ], v(y^2 + (x^7 + x^5 + x^4 + x^3 + x^2 + x)*y + x^7 + x^5 + x + 1) = 3 ],
         [ Gauss valuation induced by Valuation on rational function field induced by [ Gauss valuation induced by Trivial valuation, v(x^3 + x^2 + 1) = 1 ], v(y^3 + (x^8 + x^5 + x^4 + x^3 + x + 1)*y^2 + (x^7 + x^6 + x^5)*y + x^8 + x^5 + x^4 + x^3 + 1) = 3 ],
         [ Gauss valuation induced by Valuation on rational function field induced by [ Gauss valuation induced by Trivial valuation, v(x^3 + x^2 + 1) = 1 ], v(y^3 + (x^8 + x^4 + x^3 + x + 1)*y^2 + (x^4 + x^3 + 1)*y + x^8 + x^7 + x^4 + x + 1) = 3 ],
         [ Gauss valuation induced by Valuation on rational function field induced by [ Gauss valuation induced by Trivial valuation, v(x^3 + x^2 + 1) = 1 ], v(y^4 + (x^8 + x^7 + x^6 + x^5 + x^4 + x^3 + x^2 + x + 1)*y^3 + (x^8 + x^5 + x^4 + x^3 + x^2 + x + 1)*y^2 + (x^8 + x^7 + x^6 + x^5 + x^3 + x^2 + 1)*y + x^8 + x^7 + x^6 + x^5 + x^3 + 1) = 3 ],
         [ Gauss valuation induced by Valuation on rational function field induced by [ Gauss valuation induced by Trivial valuation, v(x^3 + x^2 + 1) = 1 ], v(y^7 + (x^8 + x^5 + x^4 + x)*y^6 + (x^7 + 1)*y^5 + (x^4 + x^2)*y^4 + (x^8 + x^3 + x + 1)*y^3 + (x^7 + x^6 + x^4 + x^2 + x + 1)*y^2 + (x^8 + x^7 + x^5 + x^3 + 1)*y + x^7 + x^6 + x^5 + x^4 + x^3 + x^2) = 3 ]]

        Some cases with trivial residue field extensions::

            sage: K.<x> = FunctionField(QQ)
            sage: S.<y> = K[]
            sage: F = y^2 - x^2 - x^3 - 3
            sage: v0 = GaussValuation(K._ring,pAdicValuation(QQ,3))
            sage: v1 = v0.extension(K._ring.gen(),1/3)
            sage: mu0 = RationalFunctionFieldValuation(K, v1)
            sage: mu0.mac_lane_approximants(F)
            [[ Gauss valuation induced by Valuation on rational function field induced by [ Gauss valuation induced by 3-adic valuation, v(x) = 1/3 ], v(y + x) = 2/3 ],
             [ Gauss valuation induced by Valuation on rational function field induced by [ Gauss valuation induced by 3-adic valuation, v(x) = 1/3 ], v(y + 2*x) = 2/3 ]]

        """
        # TODO: move this to discrete valuation
        R = G.parent()
        if R.base_ring() is not self.domain():
            raise ValueError("G must be defined over the domain of this valuation")
        if not assume_squarefree and not G.is_squarefree():
            raise ValueError("G must be squarefree")

        from sage.rings.all import infinity
        from gauss_valuation import GaussValuation

        leaves = [ GaussValuation(R, self)]
        while True:
            ef = [ v.E()*v.F() for v in leaves]
            if sum(ef) == G.degree():
                if precision_cap is None or all([v(v.phi())>precision_cap for v in leaves]):
                    return leaves

            expandables = []
            new_leaves = []
            for v in leaves:
                if v(G) is infinity:
                    new_leaves.append(v)
                else:
                    expandables.append(v)
            leaves = new_leaves

            if not expandables:
                return leaves

            for v in expandables:
                leaves.extend(v.mac_lane_step(G))

    def _repr_(self):
        return "Valuation on rational function field induced by %s"%self._base_valuation

    def extension(self, L, algorithm="mac_lane"):
        K = self.domain()
        if L is K:
            return self
        if L.base() is not K:
            raise ValueError("L must be a simple finite extension of %s"%K)

        if algorithm == "mac_lane":
            W = self.mac_lane_approximants(L.polynomial(),precision_cap=infinity)
            if len(W) > 1:
                raise ValueError("extension to %s is not unique"%L)
            return FunctionFieldPolymodValuation(L, W[0])
        else: raise ValueError()

    def uniformizer(self):
        return self.domain()(self._base_valuation.uniformizer())

class FunctionFieldPolymodValuation(DiscreteValuation):
    def __init__(self, domain, base_valuation):
        from sage.rings.all import infinity
        if base_valuation.domain() is not domain._ring:
            raise ValueError
        if base_valuation(domain.polynomial()) is not infinity:
            raise ValueError

        self._base_valuation = base_valuation
        DiscreteValuation.__init__(self, domain)

    def _call_(self, x):
        return self._base_valuation(x.element())

    def __hash__(self):
        return hash(self._base_valuation) + hash(self.domain())

    def _cache_key(self):
        return self._base_valuation, self.domain()

    def __cmp__(self, other):
        if type(self) != type(other):
            return cmp(type(self),type(other))
        return cmp(self._base_valuation,other._base_valuation)

    def shift(self, f, v):
        if f.parent() is not self.domain():
            raise ValueError

        if v == 0:
            return f
        elif v < 0:
            return f/self._base_valuation.element_with_valuation(-v)
        else:
            return f*self._base_valuation.element_with_valuation(v)

    def value_group(self):
        return self._base_valuation.value_group()

    def residue_field(self):
        return self._base_valuation.residue_field()

    def reduce(self, f):
        return self.residue_field()(self._base_valuation._base_valuation.reduce(f.element()))

    def lift(self, F):
        if F.parent() is not self.residue_field():
            raise ValueError

        return self.domain()(self._base_valuation._base_valuation.lift(F.element()))

    def mac_lane_approximants(self, G, precision_cap=None, assume_squarefree=False):
        # TODO: move this to discrete valuation
        R = G.parent()
        if R.base_ring() is not self.domain():
            raise ValueError("G must be defined over the domain of this valuation")
        if not assume_squarefree and not G.is_squarefree():
            raise ValueError("G must be squarefree")

        from sage.rings.all import infinity
        from gauss_valuation import GaussValuation

        leaves = [ GaussValuation(R, self)]
        while True:
            ef = [ v.E()*v.F() for v in leaves]
            if sum(ef) == G.degree():
                if precision_cap is None or all([v(v.phi())>precision_cap for v in leaves]):
                    return leaves

            expandables = []
            new_leaves = []
            for v in leaves:
                if v(G) is infinity:
                    new_leaves.append(v)
                else:
                    expandables.append(v)
            leaves = new_leaves

            if not expandables:
                return leaves

            for v in expandables:
                leaves.extend(v.mac_lane_step(G))

    def _repr_(self):
        return "Valuation on rational function field induced by %s"%self._base_valuation

from sage.structure.unique_representation import UniqueRepresentation
class TrivialValuation(UniqueRepresentation, DiscreteValuation):
    from gauss_valuation import classmaker
    __metaclass__ = classmaker()

    def __init__(self, domain):
        DiscreteValuation.__init__(self, domain)

    def _call_(self, x):
        if x.is_zero():
            return infinity
        else:
            return ZZ.zero()

    def shift(self, x, v):
        if x.parent() is not self.domain():
            raise ValueError

        if v == 0:
            return x
        else:
            raise ValueError

    def value_group(self):
        return self._value_group(0)

    def residue_field(self):
        if self.domain().is_field():
            return self.domain()
        else:
            raise NotImplementedError("residue ring is not a field")

    def reduce(self, x):
        if x.parent() is not self.domain():
            raise ValueError

        return x

    lift = reduce

    def _repr_(self):
        return "Trivial valuation"
