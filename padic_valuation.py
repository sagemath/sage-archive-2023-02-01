"""
The discrete valuation of a `p`-adic ring

This file makes it possible to use `p`-adics, integers, and rationals in the
general discrete valuation framework.

AUTHORS:

- Julian Rueth (2013-03-16): initial version

"""
#*****************************************************************************
#       Copyright (C) 2013 Julian Rueth <julian.rueth@fsfe.org>
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

from valuation import DiscretePseudoValuation
from limit_valuation import FiniteExtensionFromLimitValuation
from sage.structure.factory import UniqueFactory
from sage.misc.cachefunc import cached_method
from sage.misc.fast_methods import WithEqualityById

class PadicValuationFactory(UniqueFactory):
    """
    Create a ``prime``-adic valuation on ``R``.

    INPUT:

    - ``R`` -- a subring of a number field or a subring of a local field in
      characteristic zero.

    - ``prime`` -- a prime that does not split, a discrete (pseudo-)valuation,
      a fractional ideal, or ``None`` (default: ``None``)

    EXAMPLES:

    For integers and rational numbers, ``prime`` is just a prime of the
    integers::

        sage: from mac_lane import * # optional: standalone
        sage: pAdicValuation(ZZ, 3)
        3-adic valuation

        sage: pAdicValuation(QQ, 3)
        3-adic valuation

    ``prime`` may be ``None`` for local rings::

        sage: pAdicValuation(Qp(2))
        2-adic valuation

        sage: pAdicValuation(Zp(2))
        2-adic valuation

    But it must be specified in all other cases::

        sage: pAdicValuation(ZZ)
        Traceback (most recent call last):
        ...
        ValueError: prime must be specified for this ring

    For number fields, ``prime`` can be an integer that is completely ramified
    in ``R``::

        sage: pAdicValuation(GaussianIntegers().fraction_field(), 2)
        2-adic valuation

    For number fields, ``prime`` can be an integer that is unramified in ``R``:

        sage: pAdicValuation(GaussianIntegers().fraction_field(), 3)
        3-adic valuation

    The same applies if ``R`` is a subring of a number field::
    
        sage: pAdicValuation(GaussianIntegers(), 3)
        3-adic valuation

    However, this is only supported if ``prime`` does not factor into
    pairwise distinct factors::

        sage: pAdicValuation(GaussianIntegers(), 5)
        Traceback (most recent call last):
        ...
        ValueError: 5 does not single out a unique extension of 5-adic valuation to Number Field in I with defining polynomial x^2 + 1

    When ``R`` is an absolute or relative number field, or a subring thereof,
    ``prime`` can also be specified by providing a valuation on the base ring
    that has a unique extension::

        sage: pAdicValuation(GaussianIntegers(), pAdicValuation(ZZ, 3))
        3-adic valuation

    When the extension is not unique, this does not work::

        sage: pAdicValuation(GaussianIntegers(), pAdicValuation(ZZ, 5))
        Traceback (most recent call last):
        ...
        ValueError: 5-adic valuation does not single out a unique extension of 5-adic valuation to Number Field in I with defining polynomial x^2 + 1

    For a number field which is of the form `K[x]/(G)`, you can specify a
    valuation by providing a discrete pseudo-valuation on `K[x]` which sends
    `G` to `\infty`. This lets us specify which extension of the 5-adic
    valuation we care about in the above example::

        sage: R.<x> = QQ[]
        sage: v = pAdicValuation(GaussianIntegers(), GaussValuation(R, pAdicValuation(QQ, 5)).extension(x + 2, infinity))
        sage: w = pAdicValuation(GaussianIntegers(), GaussValuation(R, pAdicValuation(QQ, 5)).extension(x + 1/2, infinity))
        sage: v == w
        False

    Note that you get the same valuation, even if you write down the
    pseudo-valuation differently::

        sage: ww = pAdicValuation(GaussianIntegers(), GaussValuation(R, pAdicValuation(QQ, 5)).extension(x + 3, infinity))
        sage: w is ww
        True

    The valuation ``prime`` does not need to send the defining polynomial `G`
    to `\infty`. It is sufficient if it singles out one of the valuations on
    the number field.  This is important if the prime only factors over the
    completion, i.e., if it is not possible to write down one of the factors
    within the number field::

        sage: TODO

    Finally, ``prime`` can also be a fractional ideal of a number field if it
    singles out an extension of a `p`-adic valuation of the base field::

        sage: TODO

    TESTS:

    Check that we produce different valuations for distinguishable number
    fields:: 

        sage: from mac_lane import * # optional: standalone
        sage: R.<x> = QQ[]
        sage: L.<s> = QQ.extension(x^2 - 2)
        sage: R.<y> = L[]
        sage: N.<t> = L.extension(y^2 - s)
        sage: N = N.absolute_field('u')

        sage: M.<u> = QQ.extension(x^4 - 2)
        sage: M is N
        False
        sage: M == N
        True

        sage: pAdicValuation(QQ, 2).extension(M) is pAdicValuation(QQ, 2).extension(N)
        False

    A case that took very long due to the hashing of number fields::

        sage: R.<x> = QQ[]
        sage: Delta1= x^12 + 20*x^11 + 154*x^10 + 664*x^9 + 1873*x^8 + 3808*x^7 + 5980*x^6 + 7560*x^5 + 7799*x^4 + 6508*x^3 + 4290*x^2 + 2224*x + 887 
        sage: K.<theta>=NumberField(x^9-2)
        sage: vK=pAdicValuation(QQ, 2).extension(K)
        sage: G=Delta1.change_ring(K)
        sage: v0=GaussValuation(G.parent(), vK)
        sage: V=v0.mac_lane_step(G)
        sage: while len(V) == 1: V=V[0].mac_lane_step(G)
        sage: v=V[0]
        sage: while v(v.phi()) < 50: v=v.mac_lane_step(G)[0]
        sage: F=v.phi()
        sage: L.<alpha>=K.extension(F)
        sage: vK.mac_lane_approximants(F)
        [[ Gauss valuation induced by theta-adic valuation, v(x + 1) = 9/4 ]]
        sage: vL = vK.extension(L)

    """
    def create_key_and_extra_args(self, R, prime=None):
        r"""
        Create a unique key identifying the valuation of ``R`` with respect to
        ``prime``.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: pAdicValuation(QQ, 2) # indirect doctest
            2-adic valuation

        """
        from sage.rings.all import ZZ, QQ
        from sage.rings.padics.padic_generic import pAdicGeneric
        from sage.rings.number_field.number_field import is_NumberField

        if R.characteristic() != 0:
            # We do not support equal characteristic yet
            raise ValueError("R must be a ring of characteristic zero.")

        if R is ZZ or R is QQ:
            return self.create_key_for_integers(R, prime), {}
        elif isinstance(R, pAdicGeneric):
            return self.create_key_for_local_ring(R, prime), {}
        elif is_NumberField(R.fraction_field()):
            return self.create_key_and_extra_args_for_number_field(R, prime)
        else:
            raise NotImplementedError("p-adic valuations not implemented for %r"%(R,))

    def create_key_for_integers(self, R, prime):
        from sage.rings.all import ZZ
        if prime is None:
            raise ValueError("prime must be specified for this ring")
        if isinstance(prime, DiscretePseudoValuation):
            prime = prime.uniformizer()
        if prime not in ZZ or not ZZ(prime).is_prime():
            raise ValueError("prime must be a prime in the integers but %s is not"%(prime,))
        return R, prime

    def create_key_for_local_ring(self, R, prime):
        # We do not care much about the value of prime since there is only one
        # reasonable p-adic valuation here
        if prime is not None and R(prime).valuation() <= 0:
            raise ValueError("prime must be an element of positive valuation")

        return (R,)

    def create_key_and_extra_args_for_number_field(self, R, prime):
        from sage.rings.number_field.number_field_ideal import NumberFieldFractionalIdeal
        # To make our lives easier, we move prime to the fraction field of R
        # which we denote in the following as L = K[x]/(G), do all computations
        # there and then come back the original ring
        L = R.fraction_field()
        G = L.relative_polynomial()
        K = L.base_ring()

        if isinstance(prime, DiscretePseudoValuation):
            return self.create_key_and_extra_args_for_number_field_from_valuation(R, prime, prime)
        elif prime in K:
            return self.create_key_and_extra_args_for_number_field_from_valuation(R, pAdicValuation(K, prime), prime)
        elif prime in L or isinstance(prime, NumberFieldFractionalIdeal):
            return self.create_key_and_extra_args_for_number_field_from_ideal(R, L.fractional_ideal(prime), prime)
        else:
            raise ValueError("prime must be a discrete pseudo-valuation, a prime in the base ring, or a fractional ideal")

    def create_key_and_extra_args_for_number_field_from_valuation(self, R, v, prime):
        r"""
        prime is only used to provide more meaningful error messages
        """
        # To make our lives easier, we move prime to the fraction field of R
        # which we denote in the following as L = K[x]/(G), do all computations
        # there and then come back the original ring
        L = R.fraction_field()
        G = L.relative_polynomial()
        K = L.base_ring()

        # Lift v to a pseudo-valuation on K[x]
        # First, we lift valuations defined on subrings of K to valuations on K[x]
        if v.domain().fraction_field() is not G.parent().fraction_field():
            if v.domain() is not K:
                v = pAdicValuation(K, v)
            from gauss_valuation import GaussValuation
            v = GaussValuation(G.parent(), v)
        # Then, we lift valuations defined on polynmial rings which are subrings of K[x] to K[x]
        if v.domain() != G.parent():
            v = v.extension(G.parent())
        assert(v.domain() is G.parent())

        # To obtain uniqueness of p-adic valuations, we need a canonical
        # description of v. We consider all extensions of vK to L and select
        # the one approximated by v.
        vK = v.constant_valuation()
        candidates = vK.extensions(L)

        # We make use of the partial order on discrete pseudo-valuations to
        # single out our candidate

        # First, we see if v extends an approximant (necessarily there can be only one)
        below_v = [c for c in candidates if v >= c]
        assert(len(below_v) <= 1)
        if len(below_v) == 1:
            ret = below_v[0]
            # equality of number fields has it quirks since it says that two
            # fields are == even if they are distinguishable (because they come
            # from different constructions.)
            # Including construction() into the key seems to be a way to distinguish such cases properly.
            return (R, ret, L.construction()), {'approximants':  candidates}

        # Even if it does not extend an approximant, v could be an
        # approximation of a single approximant
        over_v = [c for c in candidates if v <= c]
        if len(over_v) > 1:
            raise ValueError("%s does not single out a unique extension of %s to %s"%(prime, vK, L))
        elif len(over_v) == 0:
            raise ValueError("%s is unrelated to extensions of %s to %s"%(prime, vK, L))
        else:
            # equality of number fields has it quirks since it says that two
            # fields are == even if they are distinguishable (because they come
            # from different constructions.)
            # Including structure() into the key seems to be a way to distinguish such cases properly.
            return (R, over_v[0], L.construction()), {'approximants': candidates}

    def create_key_and_extra_args_for_number_field_from_ideal(self, R, I, prime):
        r"""
        prime only for error messages
        """
        # To make our lives easier, we move prime to the fraction field of R
        # which we denote in the following as L = K[x]/(G), do all computations
        # there and then come back the original ring
        L = R.fraction_field()
        G = L.relative_polynomial()
        K = L.base_ring()

        # To obtain uniqueness of p-adic valuations, we need a canonical
        # description of v. We consider all extensions of vK to L and select
        # the one approximated by v.
        # Of course, this only works if I comes from a single prime downstairs.
        vK = pAdicValuation(K, I.relative_norm())
        candidates = vK.mac_lane_approximants(G)

        candidates_for_I = [c for c in candidates if all(c(g) > 0 for g in I.gens())]
        assert(len(candidates_for_I > 0)) # This should not be possible, unless I contains a unit
        if len(candidates_for_I) > 1:
            raise ValueError("%s does not single out a unique extension of %s to %s"%(prime, vK, L))
        else:
            # equality of number fields has it quirks since it says that two
            # fields are == even if they are distinguishable (because they come
            # from different constructions.)
            # Including structure() into the key seems to be a way to distinguish such cases properly.
            return (R, candidates_for_I[0], L.construction()), {'approximants': candidates}

    def create_object(self, version, key, **extra_args):
        from sage.rings.all import ZZ, QQ
        from sage.rings.padics.padic_generic import pAdicGeneric
        from valuation_space import DiscreteValuationSpace
        R = key[0]
        parent = DiscreteValuationSpace(R)
        if isinstance(R, pAdicGeneric):
            assert(len(key)==1)
            return parent.__make_element_class__(pAdicValuation_padic)(R)
        elif R is ZZ or R is QQ:
            prime = key[1]
            assert(len(key)==2)
            return parent.__make_element_class__(pAdicValuation_int)(R, prime)
        else: # Number field case
            v = key[1]
            _ = key[2] # ignored
            approximants = extra_args['approximants']
            parent = DiscreteValuationSpace(R)
            return parent.__make_element_class__(pAdicFromLimitValuation)(parent, v, R.relative_polynomial())
            # return parent.__make_element_class__(pAdicValuation_number_field)(R, v, approximants)

pAdicValuation = PadicValuationFactory("pAdicValuation")

class pAdicValuation_base(DiscretePseudoValuation):
    """
    Common base class for `p`-adic valuations on integral domains.

    INPUT:

    - ``ring`` -- an integral domain whose elements provide a method
      ``valuation`` which takes ``prime`` as a parameter

    - ``prime`` -- a prime

    EXAMPLES::

        sage: pAdicValuation(ZZ, 3)
        3-adic valuation

        sage: pAdicValuation(QQ, 5)
        5-adic valuation

     For `p`-adic rings, ``prime`` has to match the prime of the ring.

        sage: v = pAdicValuation(Zp(3), 2); v
        Traceback (most recent call last):
        ...
        ValueError: prime must be an element of positive valuation

    TESTS::

        sage: TestSuite(pAdicValuation(ZZ, 3)).run()
        sage: TestSuite(pAdicValuation(QQ, 5)).run()
        sage: TestSuite(pAdicValuation(Zp(5), 5)).run()

    """
    def __init__(self, ring, prime):
        """
        Initialization.

        TESTS::

            sage: from sage.rings.padics.padic_valuation import pAdicValuation_int # optional: integrated
            sage: isinstance(pAdicValuation(ZZ, 3), pAdicValuation_int)
            True

        """
        from rings_with_valuation import RingsWithDiscreteValuation
        from valuation_space import DiscreteValuationSpace
        DiscretePseudoValuation.__init__(self, DiscreteValuationSpace(ring))
        self._prime = prime

    def prime(self):
        """
        Return the prime `p` of this `p`-adic valuation.

        EXAMPLES::

            sage: pAdicValuation(ZZ, 3).prime()
            3

        """
        return self._prime

    def value_group(self):
        r"""
        Return the value group of this valuation.

        EXAMPLES::

            sage: pAdicValuation(ZZ, 3).value_group()
            DiscreteValueGroup(1)

        """
        from value_group import DiscreteValueGroup
        return DiscreteValueGroup(1)

    def _repr_(self):
        """
        Return a printable representation of this valuation.

        EXAMPLES::

            sage: pAdicValuation(ZZ, 3)._repr_()
            '3-adic valuation'

        """
        return "%s-adic valuation"%(self.prime())
        if isinstance(self.prime(), DiscretePseudoValuation):
            vals = self.prime()._augmentations()
            if vals[0].is_gauss_valuation():
                return "[ %r over %r ]-adic valuation"%(", ".join("v(%r) = %r"%(v.phi().change_var(self.variable_name()), v(v.phi())) for v in vals[1:]), vals[0].constant_valuation())
            else:
                return vals[0]
                return "%s-adic valuation"%(self.prime())
        else:
            return "%s-adic valuation"%(self.prime())

    def _call_(self, x):
        """
        Evaluate this valuation at ``x``.

        INPUT::

        - ``x`` --  an element in the domain of this valuation

        EXAMPLES::

            sage: pAdicValuation(ZZ, 3)(9)
            2

        """
        return x.valuation(self.prime())

    def reduce(self, x):
        """
        Reduce ``x`` modulo the ideal of elements of positive valuation.

        INPUT:

        - ``x`` -- an element in the domain of this valuation

        OUTPUT:

        An element of the :meth:`residue_field`.

        EXAMPLES::

            sage: pAdicValuation(ZZ, 3).reduce(4)
            1

        """
        x = self.domain().coerce(x)
        if self(x) < 0:
            raise ValueError("reduction is only defined for elements of non-negative valuation")
        return self.residue_field()(x)

    def lift(self, x):
        """
        Lift ``x`` from the residue field to the domain of this valuation.

        INPUT:

        - ``x`` -- an element of the :meth:`residue_field`

        EXAMPLES::

            sage: v = pAdicValuation(ZZ, 3)
            sage: xbar = v.reduce(4)
            sage: v.lift(xbar)
            1

        """
        if x.parent() is not self.residue_field():
            raise ValueError("x must be in the residue field of the valuation")

        return self.domain()(x)

    def uniformizer(self):
        """
        Return a uniformizer of this `p`-adic valuation, i.e., `p` as an
        element of the domain.

        EXAMPLES::

            sage: v = pAdicValuation(ZZ, 3)
            sage: v.uniformizer()
            3

        """
        return self.domain()(self._prime)

    def residue_ring(self):
        """
        Return the residue field of this valuation.

        EXAMPLES::

            sage: v = pAdicValuation(ZZ, 3)
            sage: v.residue_ring()
            Finite Field of size 3

        """
        from sage.rings.all import GF
        return GF(self._prime)

    def is_unramified(self, G, include_steps=False, assume_squarefree=False):
        """
        Return whether ``G`` defines an unramified extension of the `p`-adic
        completion of the domain this valuation.

        INPUT:

        - ``G`` -- a monic polynomial over the domain of this valuation

        - ``include_steps`` -- a boolean (default: ``False``); whether or not
          to include the approximate valuations that were used to determine
          the result in the return value.

        - ``assume_squarefree`` -- a boolean (default: ``False``); whether or
          not to assume that ``G`` is square-free over the `p`-adic completion
          of the domain of this valuation. Setting this to ``True`` can
          significantly improve the performance.

        EXAMPLES:

        We consider an extension as unramified if its ramification index is 1.
        Hence, a trivial extension is unramified::

            sage: pass


        NOTE: It is not enough to check whether ``G`` is irreducible in
        reduction for this.

        """
        R = G.parent()
        if R.base_ring() is not self.domain():
            raise ValueError("G must be defined over the domain of this valuation")
        if not assume_squarefree and not G.is_squarefree():
            raise ValueError("G must be squarefree")

        from sage.rings.all import infinity
        from gauss_valuation import GaussValuation

        steps = [ GaussValuation(R, self) ]
        while True:
            v = steps[-1]
            if v.E() > 1:
                ret = False
                break
            if v.F() == G.degree():
                ret = True
                break

            assert v(G) is not infinity
            if v.is_key(G):
                ret = True
                break

            next = v.mac_lane_step(G, assume_squarefree=True)
            if len(next)>1:
                ret = False
                break
            steps.append(next[0])

        if include_steps:
            return ret, steps
        else:
            return ret

    def is_totally_ramified(self, G, include_steps=False, assume_squarefree=False):
        """
        Return whether ``G`` defines a totally ramified extension.

        INPUT:

        - ``G`` -- a monic polynomial over the domain of this valuation

        - ``include_steps`` -- a boolean (default: ``False``), include the the
          valuations produced during the process of checking whether ``G`` is
          totally ramified in the return value

        - ``assume_squarefree`` -- a boolean (default: ``False``), whether to
          assume the input to be squarefree

        ALGORITHM:

        This is simplified version of :meth:`mac_lane_approximants`.

        EXAMPLES::

            sage: k=Qp(5,4)
            sage: v = pAdicValuation(k)
            sage: R.<x>=k[]
            sage: G = x^2 + 1
            sage: v.is_totally_ramified(G)
            False
            sage: G = x + 1
            sage: v.is_totally_ramified(G)
            True
            sage: G = x^2 + 2
            sage: v.is_totally_ramified(G)
            False
            sage: G = x^2 + 5
            sage: v.is_totally_ramified(G)
            True
            sage: v.is_totally_ramified(G, include_steps=True)
            (True, [Gauss valuation induced by 5-adic valuation, [ Gauss valuation induced by 5-adic valuation, v((1 + O(5^4))*x) = 1/2 ]])

        """
        R = G.parent()
        if R.base_ring() is not self.domain():
            raise ValueError("G must be defined over the domain of this valuation")
        if not assume_squarefree and not G.is_squarefree():
            raise ValueError("G must be squarefree")

        from sage.rings.all import infinity
        from gauss_valuation import GaussValuation

        steps = [ GaussValuation(R, self) ]
        while True:
            v = steps[-1]
            if v.F() > 1:
                ret = False
                break
            if v.E() == G.degree():
                ret = True
                break

            assert v(G) is not infinity
            if v.is_key(G):
                ret = False
                break

            next = v.mac_lane_step(G, assume_squarefree=True)
            if len(next)>1:
                ret = False
                break
            steps.append(next[0])

        if include_steps:
            return ret, steps
        else:
            return ret

    def montes_factorization(self, G, assume_squarefree=False):
        """
        Factor ``G`` by computing :meth:`mac_lane_approximants` to arbitrary
        precision.

        INPUT:

        - ``G`` -- a monic polynomial over the domain of this valuation

        - ``assume_squarefree`` -- a boolean (default: ``False``), whether to
          assume the input to be squarefree

        EXAMPLES::

            sage: k=Qp(5,4)
            sage: v = pAdicValuation(k)
            sage: R.<x>=k[]
            sage: G = x^2 + 1
            sage: v.montes_factorization(G) # long time
            ((1 + O(5^4))*x + 2 + 5 + 2*5^2 + 5^3 + O(5^4)) * ((1 + O(5^4))*x + 3 + 3*5 + 2*5^2 + 3*5^3 + O(5^4))

        REFERENCES:

        .. [GMN2008] Jordi Guardia, Jesus Montes, Enric Nart (2008). Newton
        polygons of higher order in algebraic number theory. arXiv:0807.2620
        [math.NT]

        """
        R = G.parent()
        if R.base_ring() is not self.domain():
            raise ValueError("G must be defined over the domain of this valuation")
        if not G.is_monic():
            raise ValueError("G must be monic")
        if not all([self(c)>=0 for c in G.coefficients()]):
            raise ValueError("G must be integral")

        from sage.rings.all import infinity
        # W contains approximate factors of G
        W = self.mac_lane_approximants(G, precision_cap=infinity,assume_squarefree=assume_squarefree)
        ret = [w.phi() for w in W]

        # the polynomials in ret give "a" factorization; there is no guarantee
        # that every factorization is of this form... TODO
        from sage.structure.factorization import Factorization
        return Factorization([ (g,1) for g in ret ], simplify=False)

    def mac_lane_approximants(self, G, precision_cap=None, assume_squarefree=False):
        """
        Compute the extensions `w` of this valuation to a polynomial ring over
        the domain which have `w(G)=\infty`.

        INPUT:

        - ``G`` -- a monic polynomial over the domain of this valuation

        - ``precision_cap`` -- a rational, ``infinity`` or ``None`` (default:
          ``None``), stop computation as soon as no new branches can show up in
          the tree of valuations and the valuation of the key polynomials is at
          least ``precision_cap`` (if not ``None``).

        - ``assume_squarefree`` -- a boolean (default: ``False``), whether to
          assume the input to be squarefree

        EXAMPLES::

            sage: k=Qp(2,10)
            sage: v = pAdicValuation(k)

            sage: R.<x>=k[]
            sage: G = x
            sage: v.mac_lane_approximants(G)
            [Gauss valuation induced by 2-adic valuation]
            sage: v.mac_lane_approximants(G, precision_cap = infinity)
            [[ Gauss valuation induced by 2-adic valuation, v((1 + O(2^10))*x) = +Infinity ]]

            sage: G = x^2 + 1
            sage: v.mac_lane_approximants(G)
            [[ Gauss valuation induced by 2-adic valuation, v((1 + O(2^10))*x + (1 + O(2^10))) = 1/2 ]]
            sage: v.mac_lane_approximants(G, precision_cap = infinity)
            [[ Gauss valuation induced by 2-adic valuation, v((1 + O(2^10))*x + (1 + O(2^10))) = 1/2, v((1 + O(2^10))*x^2 + (1 + O(2^10))) = +Infinity ]]

            sage: G = x^4 + 2*x^3 + 2*x^2 - 2*x + 2
            sage: v.mac_lane_approximants(G)
            [[ Gauss valuation induced by 2-adic valuation, v((1 + O(2^10))*x) = 1/4 ]]
            sage: v.mac_lane_approximants(G,infinity)
            [[ Gauss valuation induced by 2-adic valuation, v((1 + O(2^10))*x) = 1/4, v((1 + O(2^10))*x^4 + (2 + O(2^11))*x^3 + (2 + O(2^11))*x^2 + (2 + 2^2 + 2^3 + 2^4 + 2^5 + 2^6 + 2^7 + 2^8 + 2^9 + 2^10 + O(2^11))*x + (2 + O(2^11))) = +Infinity ]]

        The factorization of primes in the Gaussian integers can be read off
        the Mac Lane approximants::

            sage: v0 = pAdicValuation(QQ, 2)
            sage: R.<x> = QQ[]
            sage: G = x^2 + 1
            sage: v0.mac_lane_approximants(G)
            [[ Gauss valuation induced by 2-adic valuation, v(x + 1) = 1/2 ]]

            sage: v0 = pAdicValuation(QQ, 3)
            sage: v0.mac_lane_approximants(G)
            [[ Gauss valuation induced by 3-adic valuation, v(x^2 + 1) = +Infinity ]]

            sage: v0 = pAdicValuation(QQ, 5)
            sage: v0.mac_lane_approximants(G)
            [[ Gauss valuation induced by 5-adic valuation, v(x + 2) = 1 ], [ Gauss valuation induced by 5-adic valuation, v(x + 3) = 1 ]]
            sage: v0.mac_lane_approximants(G, precision_cap = 10) # long time
            [[ Gauss valuation induced by 5-adic valuation, v(x + 25670807) = 11 ], [ Gauss valuation induced by 5-adic valuation, v(x + 23157318) = 11 ]]

        The same example over the 5-adic numbers. In the quadratic extension
        `\QQ[x]/(x^2+1)`, 5 factors `-(x - 2)(x + 2)`, this behaviour can be
        read off the Mac Lane approximants::

            sage: k=Qp(5,4)
            sage: v = pAdicValuation(k)
            sage: R.<x>=k[]
            sage: G = x^2 + 1
            sage: v1,v2 = v.mac_lane_approximants(G); v1,v2
            ([ Gauss valuation induced by 5-adic valuation, v((1 + O(5^4))*x + (2 + O(5^4))) = 1 ], [ Gauss valuation induced by 5-adic valuation, v((1 + O(5^4))*x + (3 + O(5^4))) = 1 ])
            sage: w1, w2 = v.mac_lane_approximants(G,precision_cap=2); w1,w2
            ([ Gauss valuation induced by 5-adic valuation, v((1 + O(5^4))*x + (2 + 5 + 2*5^2 + O(5^4))) = 3 ], [ Gauss valuation induced by 5-adic valuation, v((1 + O(5^4))*x + (3 + 3*5 + 2*5^2 + O(5^4))) = 3 ])

        Note how the latter give a better approximation to the factors of `x^2 + 1`::

            sage: v1.phi() * v2.phi() - G
            (5 + O(5^4))*x + 5 + O(5^4)
            sage: w1.phi() * w2.phi() - G
            (5^3 + O(5^4))*x + 5^3 + O(5^4)

        In this example, the process stops with a factorization of `x^2 + 1`::

            sage: v.mac_lane_approximants(G, precision_cap=infinity)
            [[ Gauss valuation induced by 5-adic valuation, v((1 + O(5^4))*x + (2 + 5 + 2*5^2 + 5^3 + O(5^4))) = +Infinity ],
             [ Gauss valuation induced by 5-adic valuation, v((1 + O(5^4))*x + (3 + 3*5 + 2*5^2 + 3*5^3 + O(5^4))) = +Infinity ]]

        This obviously cannot happen over the rationals where we only get an
        approximate factorization::

            sage: v = pAdicValuation(QQ, 5)
            sage: R.<x>=QQ[]
            sage: G = x^2 + 1
            sage: v.mac_lane_approximants(G)
            [[ Gauss valuation induced by 5-adic valuation, v(x + 2) = 1 ], [ Gauss valuation induced by 5-adic valuation, v(x + 3) = 1 ]]
            sage: v.mac_lane_approximants(G, precision_cap=5)
            [[ Gauss valuation induced by 5-adic valuation, v(x + 14557) = 6 ], [ Gauss valuation induced by 5-adic valuation, v(x + 32318) = 7 ]]

        TESTS:

        Initial versions ran into problems with the trivial residue field
        extensions in this case::

            sage: K = Qp(3,20)
            sage: R.<T> = K[]

            sage: alpha = T^3/4
            sage: G = 3^3*T^3*(alpha^4 - alpha)^2 - (4*alpha^3 - 1)^3
            sage: G = G/G.leading_coefficient()
            sage: pAdicValuation(K).mac_lane_approximants(G) # long time
            [[ Gauss valuation induced by 3-adic valuation, v((1 + O(3^20))*T + (2 + O(3^20))) = 1/9, v((1 + O(3^20))*T^9 + (2*3 + 2*3^2 + O(3^21))*T^8 + (3 + 3^5 + O(3^21))*T^7 + (2*3 + 2*3^2 + 3^3 + 2*3^4 + 2*3^5 + 3^6 + O(3^21))*T^6 + (2*3 + 2*3^2 + 3^4 + 3^6 + 2*3^7 + O(3^21))*T^5 + (3 + 3^2 + 3^3 + 2*3^6 + 2*3^7 + 3^8 + O(3^21))*T^4 + (2*3 + 2*3^2 + 3^3 + 2*3^5 + 2*3^6 + 2*3^7 + 2*3^8 + O(3^21))*T^3 + (2*3 + 2*3^2 + 3^3 + 2*3^4 + 3^5 + 2*3^6 + 2*3^7 + 2*3^8 + O(3^21))*T^2 + (3 + 2*3^2 + 2*3^3 + 2*3^4 + 2*3^7 + 3^8 + O(3^21))*T + (2 + 2*3 + 2*3^2 + 2*3^4 + 2*3^5 + 3^7 + O(3^20))) = 55/27 ]]

        A similar example::

          sage: R.<x> = QQ[]
          sage: v = pAdicValuation(QQ, 3)
          sage: G = (x^3 + 3)^3 - 81
          sage: v.mac_lane_approximants(G)
          [[ Gauss valuation induced by 3-adic valuation, v(x) = 1/3, v(x^3 + 3*x + 3) = 13/9 ]]

        Another problematic case::

            sage: R.<x>=QQ[] 
            sage: Delta=x^12 + 20*x^11 + 154*x^10 + 664*x^9 + 1873*x^8 + 3808*x^7 + 5980*x^6 + 7560*x^5 + 7799*x^4 + 6508*x^3 + 4290*x^2 + 2224*x + 887 
            sage: K.<theta>=NumberField(x^6+108) 
            sage: K.is_galois()
            True
            sage: vK=pAdicValuation(QQ,2).extension(K)
            sage: vK(2) 
            3 
            sage: vK(theta) 
            1 
            sage: G=Delta.change_ring(K) 
            sage: V=vK.mac_lane_approximants(G) 

        OUTPUT:

        A list of valuations over the parent of ``G``.

        .. NOTE::

            For an irreducible ``G`` over `K`, these extensions correspond to
            the extensions of this valuation to the field `K[x]/(G(x))`, with
            `K` the number field corresponding to the domain of this valuation,
            see Theorem 2.1 in [ML1936]

        REFERENCES:

        .. [ML1936] Mac Lane, S. (1936). A construction for prime ideals as absolute
        values of an algebraic field. Duke Mathematical Journal, 2(3), 492-510.

        """
        R = G.parent()
        if R.base_ring() is not self.domain():
            raise ValueError("G must be defined over the domain of this valuation")
        if not assume_squarefree and not G.is_squarefree():
            raise ValueError("G must be squarefree")
        if not G.is_monic():
            raise ValueError("G must be monic")

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

    def change_ring(self, ring):
        return pAdicValuation(ring, self.prime())

    def extensions(self, ring):
        if self.domain() is ring:
            return [self]
        if self.domain().fraction_field() is not self.domain():
            if self.domain().fraction_field().is_subring(ring):
                return pAdicValuation(self.domain().fraction_field(), self).extensions(ring)
        if self.domain().is_subring(ring):
            from sage.rings.number_field.number_field import is_NumberField
            if is_NumberField(ring):
                if ring.base() is self.domain():
                    from valuation_space import DiscreteValuationSpace
                    from limit_valuation import FiniteExtensionFromInfiniteValuation
                    parent = DiscreteValuationSpace(ring)
                    approximants = self.mac_lane_approximants(ring.relative_polynomial())
                    return [parent.__make_element_class__(pAdicFromLimitValuation)(parent, approximant, ring.relative_polynomial()) for approximant in approximants]
                if ring.base() is not ring and self.domain().is_subring(ring.base()):
                    return sum([w.extensions(ring) for w in self.extensions(ring.base())], [])
        raise NotImplementedError("extensions of %r from %r to %r not implemented"%(self, self.domain(), ring))

    def restriction(self, ring):
        if self.domain().has_coerce_map_from(ring):
            return pAdicValuation(ring, self.prime())

        K = self.domain()
        if L is K:
            return self
        if L is K.fraction_field():
            return pAdicValuation(L, self)
        if L.base_ring() is not K:
            raise ValueError("L must be a simple finite extension of %r but %r is not"%(K, L))

        if algorithm == "ideal":
            I = L.ideal(self.prime()).factor()
            if len(I) > 1:
                raise ValueError("extension to %s is not unique"%L)
            I = I[0]
            return pAdicValuation(L, I)
        elif algorithm == "mac_lane":
            W = self.mac_lane_approximants(L.defining_polynomial())
            if len(W) > 1:
                raise ValueError("extension to %s is not unique"%L)
            return pAdicValuation(L, W[0])
        else: raise ValueError

class pAdicValuation_padic(pAdicValuation_base):
    """
    The `p`-adic valuation of a `p`-adic ring.

    INPUT:

    - ``R`` -- a `p`-adic ring

    EXAMPLES::

        sage: pAdicValuation(Qp(2)) #indirect doctest
        2-adic valuation

    """
    def __init__(self, R):
        """
        Initialization.

        TESTS::

            sage: from sage.rings.padics.padic_valuation import padicValuation_padic # optional: integrated
            sage: isinstance(pAdicValuation(Qp(2)), pAdicValuation_padic)
            True

        """
        pAdicValuation_base.__init__(self, R, R.prime())

    def reduce(self, x):
        """
        Reduce ``x`` modulo the ideal of elements of positive valuation.

        INPUT:

        - ``x`` -- an element of the domain of this valuation

        OUTPUT:

        An element of the :meth:`residue_field`.

        EXAMPLES::

            sage: R = Zp(3)
            sage: pAdicValuation(Zp(3)).reduce(R(4))
            1

        """
        x = self.domain().coerce(x)
        return self.residue_field()(x.residue())

    def lift(self, x):
        """
        Lift ``x`` from the residue field to the domain of this valuation.

        INPUT:

        - ``x`` -- an element of the :meth:`residue_field`

        EXAMPLES::

            sage: R = Zp(3)
            sage: v = pAdicValuation(R)
            sage: xbar = v.reduce(R(4))
            sage: v.lift(xbar)
            1 + O(3^20)

        """
        if x.parent() is not self.residue_field():
            raise ValueError("x must be in the residue field of the valuation")

        return self.domain()(x).lift_to_precision()

    def uniformizer(self):
        """
        Return a uniformizer of this valuation, i.e., `p` as an element of the
        domain.

        EXAMPLES::

            sage: v = pAdicValuation(Zp(3))
            sage: v.uniformizer()
            3 + O(3^21)

        """
        return self.domain().uniformizer()

    def shift(self, c, v):
        """
        Multiply ``c`` by a ``v``th power of the uniformizer.

        INPUT:

        - ``c`` -- an element of the domain of this valuation

        - ``v`` -- an integer

        OUTPUT:

        If the resulting element has negative valation, then it will be brought
        to the fraction field, otherwise it is an element of the domain.

        EXAMPLES::

            sage: R = Zp(3)
            sage: v = pAdicValuation(Zp(3))
            sage: v.shift(R(2),3)
            2*3^3 + O(3^23)

        """
        from sage.rings.all import ZZ
        if c.parent() is not self.domain():
            raise ValueError("c must be in the domain of the valuation")
        if not v in ZZ:
            raise ValueError("v must be an integer")
        v = ZZ(v)

        return c<<v

class pAdicValuation_int(pAdicValuation_base):
    pass

class pAdicValuation_number_field(pAdicValuation_base):
    def __init__(self, R, valuation, approximants):
        p = valuation.residue_field().characteristic()
        assert(p>0)
        pAdicValuation_base.__init__(self, R, p)
        assert(valuation in approximants)
        self._valuation = valuation
        self._approximants = approximants

    def _mac_lane_step(self):
        E,F = self._valuation.E(),self._valuation.F()
        self._valuation = self._valuation.mac_lane_step(self.domain().relative_polynomial())
        assert len(self._valuation)== 1
        self._valuation = self._valuation[0]
        assert E == self._valuation.E()
        assert F == self._valuation.F()

    def _repr_(self):
        return "%r-adic valuation"%(self._valuation)

    def _call_(self, x):
        if x.parent() is not self.domain():
            raise ValueError("x must be in the domain of the valuation")

        if x.is_zero():
            from sage.rings.all import infinity
            return infinity
        
        f = self._valuation.domain()(x.list())
        while True:
            vals = self._valuation.valuations(f)
            ret = min(vals)
            argmin = [1 for t in vals if t == ret]
            if len(argmin)>1:
                self._mac_lane_step()
                continue
            return ret*self._valuation.E()

    def reduce(self, x):
        x = self.domain().coerce(x)

        if self(x)>0:
            return self.residue_field().zero()
        f = self._valuation.domain()(x.list())
        while True:
            ret = self._valuation.reduce(f)
            if ret.degree():
                self._mac_lane_step()
                continue
            return ret[0]

    def lift(self, x):
        if x.parent() is not self.residue_field():
            raise ValueError("x must be in the residue field of the valuation")
        if x.is_one():
            return self.domain().one()
        if x.is_zero():
            return self.domain().zero()

        return self._valuation.lift(x)(self.domain().fraction_field().gen())

    def _repr_(self):
        """
        Return a printable representation of this valuation.

        EXAMPLES::

            sage: pAdicValuation(ZZ, 3)._repr_()
            '3-adic valuation'

        """
        if isinstance(self._valuation, DiscretePseudoValuation):
            from augmented_valuation import AugmentedValuation
            # print the minimal information that singles out this valuation from all approximants
            approximants = [v._augmentations() for v in self._approximants]
            augmentations = self._valuation._augmentations()
            unique_approximant = self._valuation
            for l in range(len(augmentations)):
                if len([a for a in approximants if a[:l+1] == augmentations[:l+1]]) == 1:
                    unique_approximant = augmentations[:l+1]
                    break
            if unique_approximant[0].is_gauss_valuation():
                unique_approximant[0] = unique_approximant[0].constant_valuation()
            if len(unique_approximant) == 1:
                return repr(unique_approximant[0])
            return "[ %s ]-adic valuation"%(", ".join("v(%r) = %r" if (isinstance(v, AugmentedValuation) and v.domain() == self._valuation.domain()) else repr(v) for v in unique_approximant))
        return "%s-adic valuation"%(self._valuation)

class pAdicFromLimitValuation(FiniteExtensionFromLimitValuation):
    def _to_base_domain(self, f):
        return f.polynomial()(self._base_valuation.domain().gen())
    def _from_base_domain(self, f):
        return f(self.domain().gen())
