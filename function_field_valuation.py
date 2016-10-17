# -*- coding: utf-8 -*-
r"""
Discrete valuations on function fields

AUTHORS:

- Julian Rueth (2016-10-16): initial version

"""
#*****************************************************************************
#       Copyright (C) 2016 Julian Rüth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

# Fix doctests so they work in standalone mode (when invoked with sage -t, they run within the mac_lane/ directory)
import sys, os
if hasattr(sys.modules['__main__'], 'DC') and 'standalone' in sys.modules['__main__'].DC.options.optional:
    sys.path.append(os.getcwd())
    sys.path.append(os.path.dirname(os.getcwd()))

from sage.structure.factory import UniqueFactory

from valuation import DiscretePseudoValuation
from trivial_valuation import TrivialValuation

class FunctionFieldValuationFactory(UniqueFactory):
    r"""
    Create a valuation on ``domain`` corresponding to ``prime``.

    INPUT:

    - ``domain`` -- a function field

    - ``prime`` -- a place of the function field or a valuation on a subring

    EXAMPLES::
    
        sage: from mac_lane import * # optional: standalone
        sage: K.<x> = FunctionField(QQ)
    
    We create a valuation that correspond to a finite rational place of a function
    field::

        sage: v = FunctionFieldValuation(K, 1); v
        (x - 1)-adic valuation
        sage: v(x)
        0
        sage: v(x - 1)
        1

    A place can also be specified with an irreducible polynomial::

        sage: v = FunctionFieldValuation(K, x - 1); v
        (x - 1)-adic valuation
    
    Similarly, for a finite non-rational place::
        
        sage: v = FunctionFieldValuation(K, x^2 + 1); v
        (x^2 + 1)-adic valuation
        sage: v(x^2 + 1)
        1
        sage: v(x)
        0

    Or for the infinite place::
    
        sage: v = FunctionFieldValuation(K, 1/x); v
        sage: v(x)
        -1
    
    Instead of specifying a generator of a place, we can define a valuation on a
    rational function field by giving a discrete valuation on the underlying
    polynomial ring::
    
        sage: R.<x> = QQ[]
        sage: w = GaussValuation(R, TrivialValuation(QQ)).extend(x - 1, 1)
        sage: v = FunctionFieldValuation(K, w); v
    
    Note that this allows us to specify valuations which do not correspond to a
    place of the function field::
    
        sage: w = GaussValuation(R, pAdicValuation(QQ, 2))
        sage: v = FunctionFieldValuation(K, w); v

    """
    def create_key(self, domain, prime):
        r"""
        Create a unique key which identifies the valuation given by ``prime``
        on ``domain``.

        TESTS:

        We specify a valuation on a function field by two different means and
        get the same object::

            sage: from mac_lane import * # optional: standalone
            sage: K.<x> = FunctionField(QQ)
            sage: v = FunctionFieldValuation(K, x - 1) # indirect doctest

            sage: R.<x> = QQ[]
            sage: w = GaussValuation(R, TrivialValuation(QQ)).extend(x - 1, 1)
            sage: FunctionFieldValuation(K, w) is v
            True
        
        """
        from sage.categories.function_fields import FunctionFields
        from sage.rings.valuation.valuation_space import DiscretePseudoValuationSpace
        if domain not in FunctionFields():
            raise ValueError("Domain must be a function field.")

        if prime in domain:
            # prime defines a place
            return self.create_key_from_place(domain, prime)
        elif prime in DiscretePseudoValuationSpace(domain._ring):
            # prime is a discrete (pseudo-)valuation on the polynomial ring
            # that the domain is constructed from
            return self.create_key_from_valuation(domain, prime)
        elif prime in DiscretePseudoValuationSpace(domain.base_field()):
            # prime is a discrete (pseudo-)valuation on the domain or a subfield
            # that has a unique extension
            return self.create_key_from_valuation(domain, prime)
        elif domain.base_field() is not domain:
            # prime might define a valuation on a subfield of domain that has a
            # unique extension
            base_valuation = FunctionFieldValuation(domain.base_field(), prime)
            return create_key(domain, base_valuation)

        raise NotImplementedError("argument must be a place or a pseudo-valuation on an underlying polynomial ring")

    def create_key_from_place(self, domain, generator):
        r"""
        Create a unique key which identifies the valuation at the place
        specified by ``generator``.

        TESTS:

            sage: from mac_lane import * # optional: standalone
            sage: K.<x> = FunctionField(QQ)
            sage: v = FunctionFieldValuation(K, 1/x) # indirect doctest

        """
        if generator not in domain.base_field():
            raise NotImplementedError("a place must be defined over a rational function field")

        if domain.base_field() is not domain:
            # if this is an extension field, construct the unique place over
            # the place on the subfield
            return self.create_key(domain.base_field(), FunctionFieldValuation(domain.base_field(), generator))

        if generator in domain.constant_base_field():
            # generator is a constant, we associate to it the place which
            # corresponds to the polynomial (x - generator)
            return self.create_key(domain, domain.gen() - generator)

        if generator in domain._ring:
            # generator is a polynomial
            generator = domain._ring(generator)
            if not generator.is_monic() or not generator.is_irreducible():
                raise ValueError("place must be defined by a monic irreducible polynomial")
            # we construct the corresponding valuation on the polynomial ring
            # with v(generator) = 1
            from sage.rings.valuation.gauss_valuation import GaussValuation
            valuation = GaussValuation(domain._ring, TrivialValuation(domain.constant_base_field())).extension(generator, 1)
            return self.create_key(domain, valuation)
        elif generator == ~domain.gen():
            # generator is 1/x, the infinite place
            return (domain, ~domain.gen())
        else:
            raise ValueError("a place must be given by an irreducible polynomial or the inverse of the generator")

    def create_key_from_valuation(self, domain, valuation):
        r"""
        Create a unique key which identifies the valuation which extends
        ``valuation``.

        TESTS:

            sage: from mac_lane import * # optional: standalone
            sage: K.<x> = FunctionField(QQ)
            sage: R.<x> = QQ[]
            sage: w = GaussValuation(R, TrivialValuation(QQ)).extend(x - 1, 1)
            sage: v = FunctionFieldValuation(K, w) # indirect doctest

        """
        if valuation.domain() is domain._ring:
            # TODO: check validity
            return domain, valuation
        raise NotImplementedError()

    def create_object(self, version, key, **extra_args):
        r"""
        Create the valuation specified by ``key``.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: K.<x> = FunctionField(QQ)
            sage: w = GaussValuation(R, pAdicValuation(QQ, 2))
            sage: v = FunctionFieldValuation(K, w); v # indirect doctest

        """
        domain, valuation = key
        if valuation in domain:
            assert(valuation == ~domain.gen())
            raise NotImplementedError
        from sage.rings.valuation.valuation_space import DiscreteValuationSpace
        parent = DiscreteValuationSpace(domain)
        return parent.__make_element_class__(RationalFunctionFieldValuation)(parent, valuation)

FunctionFieldValuation = FunctionFieldValuationFactory("FunctionFieldValuation")

from sage.rings.all import QQ, ZZ, infinity

class RationalFunctionFieldValuation(DiscretePseudoValuation):
    def __init__(self, parent, base_valuation):
        domain = parent.domain()
        if base_valuation.domain() is not domain._ring:
            raise ValueError

        self._base_valuation = base_valuation
        DiscretePseudoValuation.__init__(self, parent)

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
        # This is incorrect. For classical valuations on function fields, this is not correct.
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
        from sage.rings.valuation.augmented_valuation import AugmentedValuation
        if isinstance(self._base_valuation, AugmentedValuation):
            from sage.rings.valuation.gauss_valuation import GaussValuation
            if self._base_valuation._base_valuation == GaussValuation(self.domain()._ring, TrivialValuation(self.domain().constant_field())):
                if self._base_valuation._mu == 1:
                    return "(%r)-adic valuation"%(self._base_valuation.phi())
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

class FunctionFieldPolymodValuation(DiscretePseudoValuation):
    def __init__(self, domain, base_valuation):
        from sage.rings.all import infinity
        if base_valuation.domain() is not domain._ring:
            raise ValueError
        if base_valuation(domain.polynomial()) is not infinity:
            raise ValueError

        self._base_valuation = base_valuation
        DiscretePseudoValuation.__init__(self, domain)

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
