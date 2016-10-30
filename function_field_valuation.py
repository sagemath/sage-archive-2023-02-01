# -*- coding: utf-8 -*-
r"""
Discrete valuations on function fields

EXAMPLES:

We can create classical valuations that correspond to finite and infinite
places on a rational function field::

    sage: from mac_lane import * # optional: standalone
    sage: K.<x> = FunctionField(QQ)
    sage: v = FunctionFieldValuation(K, 1); v
    (x - 1)-adic valuation
    sage: v = FunctionFieldValuation(K, x^2 + 1); v
    (x^2 + 1)-adic valuation
    sage: v = FunctionFieldValuation(K, 1/x); v
    Valuation at the infinite place

Note that we can also specify valuations which do not correspond to a place of
the function field::

    sage: R.<x> = QQ[]
    sage: w = GaussValuation(R, pAdicValuation(QQ, 2))
    sage: v = FunctionFieldValuation(K, w); v
    2-adic valuation

Valuations on a rational function field can then be extended to finite
extensions::

    sage: v = FunctionFieldValuation(K, x - 1); v
    (x - 1)-adic valuation
    sage: R.<y> = K[]
    sage: L.<y> = K.extension(y^2 - x)
    sage: v.extensions(L)
    [[ Gauss valuation induced by (x - 1)-adic valuation, v(y - 1) = 1 , … ],
     [ Gauss valuation induced by (x - 1)-adic valuation, v(y + 1) = 1 , … ]]

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
from sage.rings.all import QQ, ZZ, infinity

from valuation import DiscreteValuation
from trivial_valuation import TrivialValuation
from limit_valuation import FiniteExtensionFromLimitValuation

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
        Valuation at the infinite place
        sage: v(x)
        -1
    
    Instead of specifying a generator of a place, we can define a valuation on a
    rational function field by giving a discrete valuation on the underlying
    polynomial ring::
    
        sage: R.<x> = QQ[]
        sage: w = GaussValuation(R, TrivialValuation(QQ)).augmentation(x - 1, 1)
        sage: v = FunctionFieldValuation(K, w); v
        (x - 1)-adic valuation
    
    Note that this allows us to specify valuations which do not correspond to a
    place of the function field::
    
        sage: w = GaussValuation(R, pAdicValuation(QQ, 2))
        sage: v = FunctionFieldValuation(K, w); v
        2-adic valuation

    Note that classical valuations at finite places or the infinite place are
    always normalized such that the uniformizing element has valuation 1::

        sage: K.<t> = FunctionField(GF(3))
        sage: M.<x> = FunctionField(K)
        sage: v = FunctionFieldValuation(M, x^3 - t)
        sage: v(x^3 - t)
        1

    However, if such a valuation comes out of a base change of the ground
    field, this is not the case anymore. In the example below, the unique
    extension of ``v`` to ``L`` still has valuation 1 on ``x^3 - t`` but it has
    valuation ``1/3`` on its uniformizing element  ``x - w``::

        sage: R.<w> = K[]
        sage: L.<w> = K.extension(w^3 - t)
        sage: N.<x> = FunctionField(L)
        sage: w = v.extension(N) # optional: integrated
        sage: w(x^3 - t) # optional: integrated
        1
        sage: w(x - w) # optional: integrated
        1/3

    There are several ways to create valuations on extensions of rational
    function fields::

        sage: K.<x> = FunctionField(QQ)
        sage: R.<y> = K[]
        sage: L.<y> = K.extension(y^2 - x); L
        Function field in y defined by y^2 - x

    A place that has a unique extension can just be defined downstairs::

        sage: v = FunctionFieldValuation(L, x); v
        [ Gauss valuation induced by (x)-adic valuation, v(y) = 1/2 , … ]

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
            sage: w = GaussValuation(R, TrivialValuation(QQ)).augmentation(x - 1, 1)
            sage: FunctionFieldValuation(K, w) is v
            True
        
        """
        from sage.categories.function_fields import FunctionFields
        from sage.rings.valuation.valuation_space import DiscretePseudoValuationSpace
        if domain not in FunctionFields():
            raise ValueError("Domain must be a function field.")

        if prime in DiscretePseudoValuationSpace(domain):
            # prime is already a valuation of the requested domain
            # if we returned (domain, prime), we would break caching
            # because this element has been created from a different key
            # Instead, we return the key that was used to create prime
            # so the caller gets back a correctly cached version of prime
            if not hasattr(prime, "_factory_data"):
               raise NotImplementedError("Valuations on function fields must be unique and come out of the FunctionFieldValuatino factory but %r has been created by other means"%(prime,))
            return prime._factory_data[2]

        if prime in domain:
            # prime defines a place
            return self.create_key_from_place(domain, prime)
        if prime in DiscretePseudoValuationSpace(domain._ring):
            # prime is a discrete (pseudo-)valuation on the polynomial ring
            # that the domain is constructed from
            return self.create_key_from_valuation(domain, prime)
        if domain.base_field() is not domain:
            # prime might define a valuation on a subring of domain and have a
            # unique extension to domain
            base_valuation = FunctionFieldValuation(domain.base_field(), prime)
            return self.create_key_from_valuation(domain, base_valuation)

        raise NotImplementedError("argument must be a place or a pseudo-valuation on a supported subring but %r does not satisfy this for the domain %r"%(prime, domain))

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
            return self.create_key(domain, FunctionFieldValuation(domain.base_field(), generator))

        if generator in domain.constant_base_field():
            # generator is a constant, we associate to it the place which
            # corresponds to the polynomial (x - generator)
            return self.create_key(domain, domain.gen() - generator)

        if generator in domain._ring:
            # generator is a polynomial
            generator = domain._ring(generator)
            if not generator.is_monic():
                raise ValueError("place must be defined by a monic polynomiala but %r is not monic"%(generator,))
            if not generator.is_irreducible():
                raise ValueError("place must be defined by an irreducible polynomial but %r factors over %r"%(generator, domain._ring))
            # we construct the corresponding valuation on the polynomial ring
            # with v(generator) = 1
            from sage.rings.valuation.gauss_valuation import GaussValuation
            valuation = GaussValuation(domain._ring, TrivialValuation(domain.constant_base_field())).augmentation(generator, 1)
            return self.create_key(domain, valuation)
        elif generator == ~domain.gen():
            # generator is 1/x, the infinite place
            return domain, ~domain.gen()
        else:
            raise ValueError("a place must be given by an irreducible polynomial or the inverse of the generator; %r does not define a place over %r"%(generator, domain))

    def create_key_from_valuation(self, domain, valuation):
        r"""
        Create a unique key which identifies the valuation which extends
        ``valuation``.

        TESTS:

            sage: from mac_lane import * # optional: standalone
            sage: K.<x> = FunctionField(QQ)
            sage: R.<x> = QQ[]
            sage: w = GaussValuation(R, TrivialValuation(QQ)).augmentation(x - 1, 1)
            sage: v = FunctionFieldValuation(K, w) # indirect doctest

        """
        # this should have been handled by create_key already
        assert valuation.domain() is not domain

        if valuation.domain() is domain._ring:
            if domain.base_field() is not domain:
                if valuation.constant_valuation().domain() is not domain.base_field():
                    raise ValueError("valuation must extend a valuation on the base field but %r extends %r whose domain is not %r"%(valuation, valuation.constant_valuation(), domain.base_field()))
                # Valuation is an approximant that describes a single valuation
                # on domain.
                # For uniqueness of valuations (which provides better caching
                # and easier pickling) we need to find a normal form of
                # valuation, i.e., the smallest approximant that describes this
                # valuation
                return domain, valuation.constant_valuation().mac_lane_approximant(domain.polynomial(), valuation)
            else:
                # on a rational function field K(x), any valuation on K[x] that
                # is not only a pseudo-valuation extends to a valuation on K(x)
                if not valuation.is_discrete_valuation():
                    raise ValueError("valuation must be a discrete valuation but %r is not."%(valuation,))
                return domain, valuation

        if valuation.domain().is_subring(domain.base_field()):
            # valuation is defined on a subring of this function field, try to lift it
            return self.create_key(domain, valuation.extension(domain))

        raise NotImplementedError("extension of valuation from %r to %r not implemented yet"%(valuation.domain(), domain))

    def create_object(self, version, key, **extra_args):
        r"""
        Create the valuation specified by ``key``.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: K.<x> = FunctionField(QQ)
            sage: R.<x> = QQ[]
            sage: w = GaussValuation(R, pAdicValuation(QQ, 2))
            sage: v = FunctionFieldValuation(K, w); v # indirect doctest
            2-adic valuation

        """
        domain, valuation = key
        from sage.rings.valuation.valuation_space import DiscreteValuationSpace
        parent = DiscreteValuationSpace(domain)

        if valuation in domain:
            # only the infinite place is handled with an element instead of a base valuation
            # any other element should have been replace with a valuation by _create_key_from_place
            assert(valuation == ~domain.gen())
            return parent.__make_element_class__(InfiniteRationalFunctionFieldValuation)(parent)

        if domain is valuation.domain():
            # we can not just return valuation in this case
            # as this would break uniqueness and pickling
            raise ValueError("valuation must not be a valuation on domain yet but %r is a valuation on %r"%(valuation, domain))

        if domain.base_field() is domain:
            # valuation is a base valuation on K[x] that induces a valuation on K(x)
            return parent.__make_element_class__(FiniteRationalFunctionFieldValuation)(parent, valuation)
        else:
            # valuation is a limit valuation that singles out an extension
            return parent.__make_element_class__(FunctionFieldFromLimitValuation)(parent, valuation, domain.polynomial())

        raise NotImplementedError("valuation on %r from %r on %r"%(domain, valuation, valuation.domain()))

FunctionFieldValuation = FunctionFieldValuationFactory("FunctionFieldValuation")

class FunctionFieldValuation_base(DiscreteValuation):
    r"""
    Base class for valuations on function fields.

    TESTS::

        sage: from mac_lane import * # optional: standalone
        sage: K.<x> = FunctionField(QQ)
        sage: v = FunctionFieldValuation(K, x) # indirect doctest
        sage: isinstance(v, FunctionFieldValuation_base)
        True
        sage: TestSuite(v).run()

    """
    def _extensions_from_constant_field(self, L):
        raise NotImplementedError

    def extensions(self, L):
        K = self.domain()
        from sage.categories.function_fields import FunctionFields
        if L is self.domain():
            return [self]
        if L in FunctionFields():
            if K.is_subring(K):
                if L.base() is K:
                    # L is a simple extension of the domain of this valuation
                    return [FunctionFieldValuation(L, w) for w in self.mac_lane_approximants(L.polynomial())]
                elif L.base() is not L and K.is_subring(L):
                    # recursively call this method for the tower of fields
                    from operator import add
                    return reduce(add, A, [])
                elif L.constant_field() is not K.constant_field() and K.constant_field().is_subring(L):
                    return self._extensions_from_constant_field(L)
        raise NotImplementedError("extension of %r from %r to %r not implemented"%(self, K, L))



class RationalFunctionFieldValuation_base(FunctionFieldValuation_base):
    r"""
    Base class for discrete valuations on rational function fields.

    TESTS::

        sage: from mac_lane import * # optional: standalone
        sage: K.<x> = FunctionField(GF(2))
        sage: v = FunctionFieldValuation(K, x) # indirect doctest
        sage: isinstance(v, RationalFunctionFieldValuation_base)
        True
        sage: TestSuite(v).run()

    """

class ClassicalFunctionFieldValuation_base(FunctionFieldValuation_base):
    r"""
    Base class for discrete valuations on rational function fields that come
    from points on the projective line.

    TESTS::

        sage: from mac_lane import * # optional: standalone
        sage: K.<x> = FunctionField(GF(5))
        sage: v = FunctionFieldValuation(K, x) # indirect doctest
        sage: isinstance(v, ClassicalFunctionFieldValuation_base)
        True
        sage: TestSuite(v).run()

    """
    def _test_classical_residue_field(self, **options):
        r"""
        Check correctness of the residue field of a discrete valuation at a
        classical point.

        TESTS::

            sage: from mac_lane import * # optional: standalone
            sage: K.<x> = FunctionField(QQ)
            sage: v = FunctionFieldValuation(K, x^2 + 1)
            sage: v._test_classical_residue_field()

        """
        tester = self._tester(**options)

        try:
            tester.assertTrue(self.domain().constant_field().is_subring(self.residue_field()))
        except NotImplementedError:
            # is_injective is often not implemented by morphisms
            pass

class InducedFunctionFieldValuation_base(FunctionFieldValuation_base):
    r"""
    Base class for function field valuation induced by a valuation on the
    underlying polynomial ring.

    TESTS::

        sage: from mac_lane import * # optional: standalone
        sage: K.<x> = FunctionField(QQ)
        sage: v = FunctionFieldValuation(K, x^2 + 1) # indirect doctest
        sage: TestSuite(v).run()

    """
    def __init__(self, parent, base_valuation):
        r"""
        TESTS::

            sage: from mac_lane import * # optional: standalone
            sage: K.<x> = FunctionField(QQ)
            sage: v = FunctionFieldValuation(K, x) # indirect doctest
            sage: isinstance(v, InducedFunctionFieldValuation_base)
            True
            
        """
        domain = parent.domain()
        if base_valuation.domain() is not domain._ring:
            raise ValueError("base valuation must be defined on %r but %r is defined on %r"%(domain._ring, base_valuation, base_valuation.domain()))

        self._base_valuation = base_valuation
        RationalFunctionFieldValuation_base.__init__(self, parent)

    def uniformizer(self):
        r"""
        Return a uniformizing element for this valuation.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: K.<x> = FunctionField(QQ)
            sage: FunctionFieldValuation(K, x).uniformizer()
            x
            
        """
        return self.domain()(self._base_valuation.uniformizer())

    def lift(self, F):
        r"""
        Return a lift of ``F`` to the :meth:`domain` of this valuation such
        that :meth:`reduce` returns the original element.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: K.<x> = FunctionField(QQ)
            sage: v = FunctionFieldValuation(K, x)
            sage: v.lift(0)
            0
            sage: v.lift(1)
            1

        """
        F = self.residue_ring().coerce(F)
        if F in self._base_valuation.residue_ring():
            f = self._base_valuation.lift(F)
            return self.domain()(f)
        raise NotImplementedError("lifting not implemented for this valuation")

    def value_group(self):
        r"""
        Return the value group of this valuation.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: K.<x> = FunctionField(QQ)
            sage: FunctionFieldValuation(K, x).value_group()
            Additive Abelian Group generated by 1

        """
        return self._base_valuation.value_group()

    def reduce(self, f):
        r"""
        Return the reduction of ``f`` in :meth:`residue_ring`.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: K.<x> = FunctionField(QQ)
            sage: v = FunctionFieldValuation(K, x^2 + 1)
            sage: v.reduce(x)
            u1

        """
        f = self.domain().coerce(f)

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

    def mac_lane_approximants(self, G, precision_cap=None, assume_squarefree=False):
        """
        Some difficult cases from Mark van Hoeij::

        sage: from mac_lane import * # optional: standalone
        sage: k = GF(2)
        sage: K.<x> = FunctionField(k)
        sage: R.<y> = K[]
        sage: F = y^21 + x*y^20 + (x^3 + x + 1)*y^18 + (x^3 + 1)*y^17 + (x^4 + x)*y^16 + (x^7 + x^6 + x^3 + x + 1)*y^15 + x^7*y^14 + (x^8 + x^7 + x^6 + x^4 + x^3 + 1)*y^13 + (x^9 + x^8 + x^4 + 1)*y^12 + (x^11 + x^9 + x^8 + x^5 + x^4 + x^3 + x^2)*y^11 + (x^12 + x^9 + x^8 + x^7 + x^5 + x^3 + x + 1)*y^10 + (x^14 + x^13 + x^10 + x^9 + x^8 + x^7 + x^6 + x^3 + x^2 + 1)*y^9 + (x^13 + x^9 + x^8 + x^6 + x^4 + x^3 + x)*y^8 + (x^16 + x^15 + x^13 + x^12 + x^11 + x^7 + x^3 + x)*y^7 + (x^17 + x^16 + x^13 + x^9 + x^8 + x)*y^6 + (x^17 + x^16 + x^12 + x^7 + x^5 + x^2 + x + 1)*y^5 + (x^19 + x^16 + x^15 + x^12 + x^6 + x^5 + x^3 + 1)*y^4 + (x^18 + x^15 + x^12 + x^10 + x^9 + x^7 + x^4 + x)*y^3 + (x^22 + x^21 + x^20 + x^18 + x^13 + x^12 + x^9 + x^8 + x^7 + x^5 + x^4 + x^3)*y^2 + (x^23 + x^22 + x^20 + x^17 + x^15 + x^14 + x^12 + x^9)*y + x^25 + x^23 + x^19 + x^17 + x^15 + x^13 + x^11 + x^5
        sage: x = K._ring.gen()
        sage: v0 = FunctionFieldValuation(K, GaussValuation(K._ring, TrivialValuation(k)).augmentation(x,1))
        sage: v0.mac_lane_approximants(F, assume_squarefree=True) # optional: integrated; assumes squarefree for speed
        [[ Gauss valuation induced by Valuation on rational function field induced by [ Gauss valuation induced by Trivial valuation, v(x) = 1 ], v(y + x + 1) = 3/2 ],
         [ Gauss valuation induced by Valuation on rational function field induced by [ Gauss valuation induced by Trivial valuation, v(x) = 1 ], v(y) = 4/3, v(y^3 + x^4) = 13/3 ],
         [ Gauss valuation induced by Valuation on rational function field induced by [ Gauss valuation induced by Trivial valuation, v(x) = 1 ], v(y + x) = 2 ],
         [ Gauss valuation induced by Valuation on rational function field induced by [ Gauss valuation induced by Trivial valuation, v(x) = 1 ], v(y^15 + y^13 + (x + 1)*y^12 + x*y^11 + (x + 1)*y^10 + y^9 + y^8 + x*y^6 + x*y^5 + y^4 + y^3 + y^2 + (x + 1)*y + x + 1) = 2 ]]
        sage: v0 = FunctionFieldValuation(K, GaussValuation(K._ring, TrivialValuation(k)).augmentation(x+1,1))
        sage: v0.mac_lane_approximants(F, assume_squarefree=True) # optional: integrated; assumes squarefree for speed
        [[ Gauss valuation induced by Valuation on rational function field induced by [ Gauss valuation induced by Trivial valuation, v(x + 1) = 1 ], v(y) = 7/2, v(y^2 + x^7 + x^6 + x^5 + x^4 + x^3 + x^2 + x + 1) = 15/2 ],
         [ Gauss valuation induced by Valuation on rational function field induced by [ Gauss valuation induced by Trivial valuation, v(x + 1) = 1 ], v(y + x^2 + 1) = 7/2 ],
         [ Gauss valuation induced by Valuation on rational function field induced by [ Gauss valuation induced by Trivial valuation, v(x + 1) = 1 ], v(y) = 3/4, v(y^4 + x^3 + x^2 + x + 1) = 15/4 ],
         [ Gauss valuation induced by Valuation on rational function field induced by [ Gauss valuation induced by Trivial valuation, v(x + 1) = 1 ], v(y^13 + x*y^12 + y^10 + (x + 1)*y^9 + (x + 1)*y^8 + x*y^7 + x*y^6 + (x + 1)*y^4 + y^3 + (x + 1)*y^2 + 1) = 2 ]]
        sage: v0 = FunctionFieldValuation(K, GaussValuation(K._ring, TrivialValuation(k)).augmentation(x^3+x^2+1,1))
        sage: v0.mac_lane_approximants(F, assume_squarefree=True) # optional: integrated; assumes squarefree for speed
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
            sage: v1 = v0.augmentation(K._ring.gen(),1/3)
            sage: mu0 = FunctionFieldValuation(K, v1)
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
        from sage.rings.valuation.gauss_valuation import GaussValuation
        if isinstance(self._base_valuation, AugmentedValuation):
            if self._base_valuation._base_valuation == GaussValuation(self.domain()._ring, TrivialValuation(self.domain().constant_field())):
                if self._base_valuation._mu == 1:
                    return "(%r)-adic valuation"%(self._base_valuation.phi())
        if self._base_valuation == GaussValuation(self.domain()._ring, self._base_valuation.constant_valuation()):
            return repr(self._base_valuation.constant_valuation())
        return "Valuation on rational function field induced by %s"%self._base_valuation

    def _extensions_from_constant_field(self, L):
        # extend the underlying valuation on the polynomial ring
        W = self._base_valuation.extensions(L._ring)
        return [FunctionFieldValuation(L, w) for w in W]

class InfiniteRationalFunctionFieldValuation(RationalFunctionFieldValuation_base, ClassicalFunctionFieldValuation_base):
    r"""
    Valuation of the infinite place of a function field.

    EXAMPLES::

        sage: from mac_lane import * # optional: standalone
        sage: K.<x> = FunctionField(QQ)
        sage: v = FunctionFieldValuation(K, 1/x) # indirect doctest

    TESTS::

        sage: TestSuite(v).run()

    """
    def _call_(self, f):
        r"""
        Evaluate this valuation at the rational function ``f``.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: K.<x> = FunctionField(QQ)
            sage: v = FunctionFieldValuation(K, 1/x)
            sage: v(x/(x^2 + 1))
            1

        """
        if f.is_zero():
            return infinity
        return f.denominator().degree() - f.numerator().degree()

    def _repr_(self):
        r"""
        Return a printable representation of this valuation.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: K.<x> = FunctionField(QQ)
            sage: FunctionFieldValuation(K, 1/x) # indirect doctest
            Valuation at the infinite place

        """
        return "Valuation at the infinite place"

    def residue_ring(self):
        r"""
        Return the residue field of this valuation.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: K.<x> = FunctionField(QQ)
            sage: FunctionFieldValuation(K, 1/x).residue_ring()
            Rational Field

        """
        return self.domain().constant_field()

    def uniformizer(self):
        r"""
        Return a uniformizer of this valuation.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: K.<x> = FunctionField(QQ)
            sage: FunctionFieldValuation(K, 1/x).uniformizer()
            1/x

        """
        return ~self.domain().gen()

    def reduce(self, f):
        r"""
        Return the reduction of ``f`` with respect to this valuation.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: K.<x> = FunctionField(GF(2))
            sage: v = FunctionFieldValuation(K, 1/x)
            sage: v.reduce(1/x)
            0
            sage: v.reduce( x*(x+1) / (x^2+x+1) )
            1
            sage: v.reduce(x)
            Traceback (most recent call last):
            ...
            ValueError: reduction is only defined for elements of non-negative valuation

        """
        f = self.domain().coerce(f)

        if self(f) < 0:
            raise ValueError("reduction is only defined for elements of non-negative valuation")
        if self(f) > 0:
            return self.residue_ring().zero()

        from sage.rings.function_field.function_field_valuation import FunctionFieldValuation
        return FunctionFieldValuation(self.domain(), self.domain().gen()).reduce(f.element()(~self.domain().gen()))

    def lift(self, F):
        r"""
        Lift ``F`` from :meth:`residue_ring` to :meth:`domain` such that
        :meth:`reduce` of the lift produces ``F``.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: K.<x> = FunctionField(GF(2))
            sage: v = FunctionFieldValuation(K, 1/x)

            sage: v.lift(v.residue_ring().one())
            1

        """
        F = self.residue_ring().coerce(F)
        # the infinite place has a trivial residue field extension
        assert F in self.domain()
        return self.domain()(F)

class FiniteRationalFunctionFieldValuation(ClassicalFunctionFieldValuation_base, InducedFunctionFieldValuation_base, RationalFunctionFieldValuation_base):
    r"""
    Valuation of a finite place of a function field.

    EXAMPLES::

        sage: from mac_lane import * # optional: standalone
        sage: K.<x> = FunctionField(QQ)
        sage: v = FunctionFieldValuation(K, x + 1); v # indirect doctest
        (x + 1)-adic valuation

    A finite place with residual degree::

        sage: w = FunctionFieldValuation(K, x^2 + 1); w
        (x^2 + 1)-adic valuation

    A finite place with ramification::

        sage: K.<t> = FunctionField(GF(3))
        sage: L.<x> = FunctionField(K)
        sage: u = FunctionFieldValuation(L, x^3 - t); u
        (x^3 + 2*t)-adic valuation

    A finite place with residual degree and ramification::

        sage: q = FunctionFieldValuation(L, x^6 - t); q
        (x^6 + 2*t)-adic valuation

    TESTS::

        sage: TestSuite(v).run()
        sage: TestSuite(w).run()
        sage: TestSuite(u).run() # long time
        sage: TestSuite(q).run() # long time

    """
    def _call_(self, f):
        r"""
        Evaluate this valuation at the function ``f``.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: K.<x> = FunctionField(QQ)
            sage: v = FunctionFieldValuation(K, x) # indirect doctest
            sage: v((x+1)/x^2)
            -2
            
        """
        return self._base_valuation(f.numerator()) - self._base_valuation(f.denominator())

    def residue_ring(self):
        r"""
        Return the residue field of this valuation.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: K.<x> = FunctionField(QQ)
            sage: FunctionFieldValuation(K, 1/x).residue_ring()
            Rational Field

        """
        return self._base_valuation.residue_field()

class FunctionFieldFromLimitValuation(FiniteExtensionFromLimitValuation):
    def _to_base_domain(self, f):
        return f.element()
