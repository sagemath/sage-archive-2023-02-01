# -*- coding: utf-8 -*-
"""
`p`-adic valuations on number fields and their subrings and completions.

AUTHORS:

- Julian Rüth (2013-03-16): initial version

"""
#*****************************************************************************
#       Copyright (C) 2013-2016 Julian Rüth <julian.rueth@fsfe.org>
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

from valuation import DiscreteValuation
from mapped_valuation import FiniteExtensionFromLimitValuation
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
        ValueError: The valuation Gauss valuation induced by 5-adic valuation does not approximate a unique extension of 5-adic valuation with respect to x^2 + 1

    When ``R`` is an absolute or relative number field, or a subring thereof,
    ``prime`` can also be specified by providing a valuation on the base ring
    that has a unique extension::

        sage: pAdicValuation(CyclotomicField(5), pAdicValuation(ZZ, 5))
        5-adic valuation

    When the extension is not unique, this does not work::

        sage: pAdicValuation(GaussianIntegers(), pAdicValuation(ZZ, 5))
        Traceback (most recent call last):
        ...
        ValueError: The valuation Gauss valuation induced by 5-adic valuation does not approximate a unique extension of 5-adic valuation with respect to x^2 + 1

    For a number field which is of the form `K[x]/(G)`, you can specify a
    valuation by providing a discrete pseudo-valuation on `K[x]` which sends
    `G` to `\infty`. This lets us specify which extension of the 5-adic
    valuation we care about in the above example::

        sage: R.<x> = QQ[]
        sage: v = pAdicValuation(GaussianIntegers(), GaussValuation(R, pAdicValuation(QQ, 5)).augmentation(x + 2, infinity))
        sage: w = pAdicValuation(GaussianIntegers(), GaussValuation(R, pAdicValuation(QQ, 5)).augmentation(x + 1/2, infinity))
        sage: v == w
        False

    Note that you get the same valuation, even if you write down the
    pseudo-valuation differently::

        sage: ww = pAdicValuation(GaussianIntegers(), GaussValuation(R, pAdicValuation(QQ, 5)).augmentation(x + 3, infinity))
        sage: w is ww
        True

    The valuation ``prime`` does not need to send the defining polynomial `G`
    to `\infty`. It is sufficient if it singles out one of the valuations on
    the number field.  This is important if the prime only factors over the
    completion, i.e., if it is not possible to write down one of the factors
    within the number field::

        sage: v = GaussValuation(R, pAdicValuation(QQ, 5)).augmentation(x + 3, 1)
        sage: pAdicValuation(GaussianIntegers().fraction_field(), v)
        [ 5-adic valuation, v(x + 3) = 1 ]-adic valuation

    Finally, ``prime`` can also be a fractional ideal of a number field if it
    singles out an extension of a `p`-adic valuation of the base field::

        sage: R = GaussianIntegers()
        sage: I = R.fraction_field().gen()
        sage: pAdicValuation(R, R.fractional_ideal(I + 1))
        2-adic valuation

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
        r"""
        Create a unique key identifying the valuation of ``R`` with respect to
        ``prime``.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: pAdicValuation(QQ, 2) # indirect doctest
            2-adic valuation

        """
        from sage.rings.all import ZZ
        if prime is None:
            raise ValueError("prime must be specified for this ring")
        from valuation import DiscretePseudoValuation
        if isinstance(prime, DiscretePseudoValuation):
            prime = prime.uniformizer()
        if prime not in ZZ or not ZZ(prime).is_prime():
            raise ValueError("prime must be a prime in the integers but %s is not"%(prime,))
        return R, prime

    def create_key_for_local_ring(self, R, prime):
        r"""
        Create a unique key identifying the valuation of ``R`` with respect to
        ``prime``.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: pAdicValuation(Qp(2)) # indirect doctest
            2-adic valuation

        """
        # We do not care much about the value of prime since there is only one
        # reasonable p-adic valuation here
        if prime is not None and R(prime).valuation() <= 0:
            raise ValueError("prime must be an element of positive valuation")

        return (R,)

    def create_key_and_extra_args_for_number_field(self, R, prime):
        r"""
        Create a unique key identifying the valuation of ``R`` with respect to
        ``prime``.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: pAdicValuation(GaussianIntegers(), 2) # indirect doctest
            2-adic valuation

        """
        from sage.rings.number_field.number_field_ideal import NumberFieldFractionalIdeal
        # To make our lives easier, we move prime to the fraction field of R
        # which we denote in the following as L = K[x]/(G), do all computations
        # there and then come back the original ring
        L = R.fraction_field()
        G = L.relative_polynomial()
        K = L.base_ring()

        from valuation import DiscretePseudoValuation
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
        Create a unique key identifying the valuation of ``R`` with respect to
        ``v``.

        .. NOTE::

            ``prime``, the original parameter that was passed to
            :meth:`create_key_and_extra_args``, is only used to provide more
            meaningful error messages

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: pAdicValuation(GaussianIntegers(), pAdicValuation(ZZ, 2)) # indirect doctest
            2-adic valuation

        """
        # To make our lives easier, we move prime to the fraction field of R
        # which we denote in the following as L = K[x]/(G), do all computations
        # there and then come back the original ring
        L = R.fraction_field()
        K = L.base_ring()
        G = L.relative_polynomial().change_ring(K)

        # Lift v to a pseudo-valuation on K[x]
        # First, we lift valuations defined on subrings of K to valuations on K[x]
        if v.domain().fraction_field() is not G.parent().fraction_field():
            if v.domain().is_subring(K):
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
        vK = v.restriction(v.domain().base_ring())
        approximants = vK.mac_lane_approximants(L.relative_polynomial())
        approximant = vK.mac_lane_approximant(L.relative_polynomial(), v, approximants=approximants)

        return (R, approximant, L.construction()), {'approximants': approximants}

    def create_key_and_extra_args_for_number_field_from_ideal(self, R, I, prime):
        r"""
        Create a unique key identifying the valuation of ``R`` with respect to
        ``I``.

        .. NOTE::

            ``prime``, the original parameter that was passed to
            :meth:`create_key_and_extra_args``, is only used to provide more
            meaningful error messages

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: pAdicValuation(GaussianIntegers(), GaussianIntegers().ideal(2)) # indirect doctest
            2-adic valuation

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
        p = I.relative_norm()
        F = p.factor()
        if len(F) != 1:
            raise ValueError("%r does not lie over a single prime of %r"%(I, K))
        vK = pAdicValuation(K, F[0][0])
        candidates = vK.mac_lane_approximants(G)

        candidates_for_I = [c for c in candidates if all(c(g.polynomial()) > 0 for g in I.gens())]
        assert(len(candidates_for_I) > 0) # This should not be possible, unless I contains a unit
        if len(candidates_for_I) > 1:
            raise ValueError("%s does not single out a unique extension of %s to %s"%(prime, vK, L))
        else:
            # equality of number fields has it quirks since it says that two
            # fields are == even if they are distinguishable (because they come
            # from different constructions.)
            # Including structure() into the key seems to be a way to distinguish such cases properly.
            # This used to be an issue but seems to be fixed, namely, the
            # absolute_field of a number field was deemed equivalent to the
            # directly created absolute field, even though the absolute_field
            # carried the information where it came from
            return (R, candidates_for_I[0], L.construction()), {'approximants': candidates}

    def create_object(self, version, key, **extra_args):
        r"""
        Create a `p`-adic valuation from ``key``.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: pAdicValuation(ZZ, 5) # indirect doctest
            5-adic valuation

        """
        from sage.rings.all import ZZ, QQ
        from sage.rings.padics.padic_generic import pAdicGeneric
        from valuation_space import DiscretePseudoValuationSpace
        R = key[0]
        K = R.fraction_field()
        parent = DiscretePseudoValuationSpace(R)
        if isinstance(R, pAdicGeneric):
            assert(len(key)==1)
            return parent.__make_element_class__(pAdicValuation_padic)(parent)
        elif R is ZZ or R is QQ:
            prime = key[1]
            assert(len(key)==2)
            return parent.__make_element_class__(pAdicValuation_int)(parent, prime)
        else: # Number field case
            v = key[1]
            _ = key[2] # ignored
            approximants = extra_args['approximants']
            parent = DiscretePseudoValuationSpace(R)
            return parent.__make_element_class__(pAdicFromLimitValuation)(parent, v, K.relative_polynomial().change_ring(R.base()), approximants)

pAdicValuation = PadicValuationFactory("pAdicValuation")

class pAdicValuation_base(DiscreteValuation):
    """
    Abstract base class for `p`-adic valuations.

    INPUT:

    - ``ring`` -- an integral domain

    - ``p`` -- a rational prime over which this valuation lies, not
      necessarily a uniformizer for the valuation

    EXAMPLES::

        sage: from mac_lane import * # optional: standalone
        sage: pAdicValuation(ZZ, 3)
        3-adic valuation

        sage: pAdicValuation(QQ, 5)
        5-adic valuation

     For `p`-adic rings, ``p`` has to match the `p` of the ring.

        sage: v = pAdicValuation(Zp(3), 2); v
        Traceback (most recent call last):
        ...
        ValueError: prime must be an element of positive valuation

    TESTS::

        sage: TestSuite(pAdicValuation(ZZ, 3)).run() # long time
        sage: TestSuite(pAdicValuation(QQ, 5)).run() # long time
        sage: TestSuite(pAdicValuation(Zp(5), 5)).run() # long time

    """
    def __init__(self, parent, p):
        """
        TESTS::

            sage: from mac_lane import * # optional: standalone
            sage: isinstance(pAdicValuation(ZZ, 3), pAdicValuation_base)
            True

        """
        DiscreteValuation.__init__(self, parent)

        from sage.rings.all import ZZ
        self._p = ZZ(p)

    def p(self):
        """
        Return the `p` of this `p`-adic valuation.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: pAdicValuation(GaussianIntegers(), 2).p()
            2

        """
        return self._p

    def reduce(self, x):
        """
        Reduce ``x`` modulo the ideal of elements of positive valuation.

        INPUT:

        - ``x`` -- an element in the domain of this valuation

        OUTPUT:

        An element of the :meth:`residue_field`.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: v = pAdicValuation(ZZ, 3)
            sage: v.reduce(4)
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
        x = self.residue_field().coerce(x)

        return self.domain()(x)

    def is_unramified(self, G, include_steps=False, assume_squarefree=False):
        """
        Return whether ``G`` defines a single unramified extension of the
        completion of the domain of this valuation.

        INPUT:

        - ``G`` -- a monic squarefree polynomial over the domain of this valuation

        - ``include_steps`` -- a boolean (default: ``False``); whether to
          include the approximate valuations that were used to determine the
          result in the return value.

        - ``assume_squarefree`` -- a boolean (default: ``False``); whether to
          assume that ``G`` is square-free over the completion of the domain of
          this valuation. Setting this to ``True`` can significantly improve
          the performance.

        EXAMPLES:

        We consider an extension as unramified if its ramification index is 1.
        Hence, a trivial extension is unramified::

            sage: from mac_lane import * # optional: standalone
            sage: R.<x> = QQ[]
            sage: v = pAdicValuation(QQ, 2)
            sage: v.is_unramified(x)
            True

        If ``G`` remains irreducible in reduction, then it defines an
        unramified extension::

            sage: v.is_unramified(x^2 + x + 1)
            True

        However, even if ``G`` factors, it might define an unramified
        extension::

            sage: v.is_unramified(x^2 + 2*x + 4)
            True

        """
        R = G.parent()

        from sage.rings.polynomial.polynomial_ring import is_PolynomialRing
        if not is_PolynomialRing(R) or R.base_ring() is not self.domain() or not len(R.gens()) == 1 or not G.is_monic():
            raise ValueError("G must be a monic univariate polynomial over the domain of this valuation")
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
        Return whether ``G`` defines a single totally ramified extension of the
        completion of the domain of this valuation.

        INPUT:

        - ``G`` -- a monic squarefree polynomial over the domain of this valuation

        - ``include_steps`` -- a boolean (default: ``False``); where to include
          the valuations produced during the process of checking whether ``G``
          is totally ramified in the return value

        - ``assume_squarefree`` -- a boolean (default: ``False``); whether to
          assume that ``G`` is square-free over the completion of the domain of
          this valuation. Setting this to ``True`` can significantly improve
          the performance.

        ALGORITHM:

        This is a simplified version of :meth:`mac_lane_approximants`.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
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

        from sage.rings.polynomial.polynomial_ring import is_PolynomialRing
        if not is_PolynomialRing(R) or R.base_ring() is not self.domain() or not len(R.gens()) == 1 or not G.is_monic():
            raise ValueError("G must be a monic univariate polynomial over the domain of this valuation")
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

    def change_domain(self, ring):
        r"""
        Change the domain of this valuation to ``ring`` if possible.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: v = pAdicValuation(ZZ, 2)
            sage: v.change_domain(QQ).domain()
            Rational Field

        """
        return pAdicValuation(ring, self.p())

    def extensions(self, ring):
        r"""
        Return the extensions of this valuation to ``ring``.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: v = pAdicValuation(ZZ, 2)
            sage: v.extensions(GaussianIntegers())
            [2-adic valuation]

        """
        if self.domain() is ring:
            return [self]
        if self.domain().fraction_field() is not self.domain():
            if self.domain().fraction_field().is_subring(ring):
                return pAdicValuation(self.domain().fraction_field(), self).extensions(ring)
        if self.domain().is_subring(ring):
            from sage.rings.number_field.number_field import is_NumberField
            if is_NumberField(ring.fraction_field()):
                if ring.base().fraction_field() is self.domain().fraction_field():
                    from valuation_space import DiscretePseudoValuationSpace
                    parent = DiscretePseudoValuationSpace(ring)
                    approximants = self.mac_lane_approximants(ring.fraction_field().relative_polynomial().change_ring(self.domain()))
                    return [pAdicValuation(ring, approximant) for approximant in approximants]
                if ring.base() is not ring and self.domain().is_subring(ring.base()):
                    return sum([w.extensions(ring) for w in self.extensions(ring.base())], [])
        raise NotImplementedError("can not compute extensions of %r from %r to %r"%(self, self.domain(), ring))

    def restriction(self, ring):
        r"""
        Return the restriction of this valuation to ``ring``.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: v = pAdicValuation(GaussianIntegers(), 2)
            sage: v.restriction(ZZ)
            2-adic valuation

        """
        if not ring.is_subring(self.domain()):
            raise ValueError("ring must be a subring of the domain of this valuation but %r is not a subring of %r"%(ring, self.domain()))

        return pAdicValuation(ring, self.p())

class pAdicValuation_padic(pAdicValuation_base):
    """
    The `p`-adic valuation of a complete `p`-adic ring.

    INPUT:

    - ``R`` -- a `p`-adic ring

    EXAMPLES::

        sage: from mac_lane import * # optional: standalone
        sage: v = pAdicValuation(Qp(2)); v #indirect doctest
        2-adic valuation

    TESTS::

        sage: TestSuite(v).run() # optional: integrated, long time

    """
    def __init__(self, parent):
        """
        TESTS::

            sage: from mac_lane import * # optional: standalone
            sage: from sage.rings.padics.padic_valuation import padicValuation_padic # optional: integrated
            sage: isinstance(pAdicValuation(Qp(2)), pAdicValuation_padic)
            True

        """
        pAdicValuation_base.__init__(self, parent, parent.domain().prime())

    def reduce(self, x):
        """
        Reduce ``x`` modulo the ideal of elements of positive valuation.

        INPUT:

        - ``x`` -- an element of the domain of this valuation

        OUTPUT:

        An element of the :meth:`residue_field`.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R = Zp(3)
            sage: pAdicValuation(Zp(3)).reduce(R(4))
            1

        """
        x = self.domain().coerce(x)
        return self.residue_field()(x.residue())

    def lift(self, x):
        """
        Lift ``x`` from the :meth:`residue_field` to the :meth:`domain` of this
        valuation.

        INPUT:

        - ``x`` -- an element of the :meth:`residue_field`

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R = Zp(3)
            sage: v = pAdicValuation(R)
            sage: xbar = v.reduce(R(4))
            sage: v.lift(xbar)
            1 + O(3^20)

        """
        x = self.residue_field().coerce(x)
        return self.domain()(x).lift_to_precision()

    def uniformizer(self):
        """
        Return a uniformizer of this valuation.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: v = pAdicValuation(Zp(3))
            sage: v.uniformizer()
            3 + O(3^21)

        """
        return self.domain().uniformizer()

    def shift(self, c, v):
        """
        Multiply ``c`` by a power of a uniformizer such that its valuation
        changes by ``v``.

        INPUT:

        - ``c`` -- an element of the domain of this valuation

        - ``v`` -- an integer

        OUTPUT:

        If the resulting element has negative valation, it will be brought to
        the fraction field, otherwise it is an element of the domain.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: R = Zp(3)
            sage: v = pAdicValuation(Zp(3))
            sage: v.shift(R(2),3)
            2*3^3 + O(3^23)

        TESTS:

        Make sure that we work around the following bug in p-adics::

            sage: R = Zp(3)
            sage: R(1) >> 1 # this should throw an exception
            O(3^19)

        However, our shift gets this right::

            sage: v = pAdicValuation(Zp(3))
            sage: v.shift(R(1), -1)
            Traceback (most recent call last):
            ...
            ValueError: can not shift down 1 + O(3^20) by 1 in 3-adic Ring with capped relative precision 20 with respect to 3-adic valuation

        Note that the same happens for ZpCA::

            sage: R = ZpCA(3)
            sage: R(1) >> 1 # this should throw an exception
            O(3^19)

        But we also detect that one::

            sage: v = pAdicValuation(ZpCA(3))
            sage: v.shift(R(1), -1)
            Traceback (most recent call last):
            ...
            ValueError: can not shift down 1 + O(3^20) by 1 in 3-adic Ring with capped absolute precision 20 with respect to 3-adic valuation

        """
        from sage.rings.all import QQ, ZZ
        c = self.domain().coerce(c)
        if c == 0:
            if v == 0:
                return c
            raise ValueError("can not shift a zero")

        v = QQ(v)
        if v not in self.value_group():
            raise ValueError("%r is not in the value group of %r"%(v, self))
        v = ZZ(v * self.domain().ramification_index())

        # Work around a bug in ZpCR/ZpCA p-adics
        # TODO: fix this upstream
        from sage.rings.padics.padic_capped_absolute_element import pAdicCappedAbsoluteElement
        from sage.rings.padics.padic_capped_relative_element import pAdicCappedRelativeElement
        from sage.categories.fields import Fields
        if (isinstance(c, pAdicCappedAbsoluteElement) or (isinstance(c, pAdicCappedRelativeElement) and c.parent() not in Fields())) and -v > self(c):
            raise ValueError("can not shift down %r by %r in %r with respect to %r"%(c, -v, self.domain(), self))

        return c<<v

    def _repr_(self):
        """
        Return a printable representation of this valuation.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: pAdicValuation(ZZ, 3)._repr_()
            '3-adic valuation'

        """
        return "%s-adic valuation"%(self.p())

    def _call_(self, x):
        r"""
        Evaluate this valuation at ``x``.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: K = Qp(3)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - 3)
            sage: pAdicValuation(L, 3)(3)
            1

        """
        return x.ordp()

    def residue_ring(self):
        r"""
        Return the residue field of this valuation.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: pAdicValuation(Qq(9, names='a'), 3).residue_ring()
            Finite Field in a0 of size 3^2

        """
        return self.domain().residue_field()

class pAdicValuation_int(pAdicValuation_base):
    r"""
    A `p`-adic valuation on the integers or the rationals.

    EXAMPLES::

        sage: from mac_lane import * # optional: standalone
        sage: v = pAdicValuation(ZZ, 3); v
        3-adic valuation

    TESTS::

        sage: TestSuite(v).run() # long time

    """
    def _repr_(self):
        """
        Return a printable representation of this valuation.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: pAdicValuation(ZZ, 3)._repr_()
            '3-adic valuation'

        """
        return "%s-adic valuation"%(self.p())

    def _call_(self, x):
        """
        Evaluate this valuation at ``x``.

        INPUT::

        - ``x`` --  an element in the domain of this valuation

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: pAdicValuation(ZZ, 3)(9)
            2

        """
        return x.valuation(self.p())

    def uniformizer(self):
        """
        Return a uniformizer of this `p`-adic valuation, i.e., `p` as an
        element of the domain.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: v = pAdicValuation(ZZ, 3)
            sage: v.uniformizer()
            3

        """
        return self.domain()(self.p())

    def residue_ring(self):
        """
        Return the residue field of this valuation.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: v = pAdicValuation(ZZ, 3)
            sage: v.residue_ring()
            Finite Field of size 3

        """
        from sage.rings.all import GF
        return GF(self.p())

    def _ge_(self, other):
        r"""
        Return whether this valuation is greater than or equal than ``other``
        everywhere.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: v = pAdicValuation(ZZ, 2)
            sage: w = TrivialValuation(ZZ)
            sage: v >= w
            True

        """
        if other.is_trivial():
            return other.is_discrete_valuation()
        if isinstance(other, pAdicValuation_int):
            return self.p() == other.p()
        return super(pAdicValuation_base, self)._ge_(other)


class pAdicFromLimitValuation(FiniteExtensionFromLimitValuation, pAdicValuation_base):
    r"""
    A `p`-adic valuation on a number field or a subring thereof, i.e., a
    valuation that extends the `p`-adic valuation on the integers.

    EXAMPLES::

        sage: from mac_lane import * # optional: standalone
        sage: v = pAdicValuation(GaussianIntegers(), 3); v
        3-adic valuation

    TESTS::

        sage: TestSuite(v).run() # long time

    """
    def __init__(self, parent, approximant, G, approximants):
        r"""
        TESTS::

            sage: from mac_lane import * # optional: standalone
            sage: v = pAdicValuation(GaussianIntegers(), 3)
            sage: isinstance(v, pAdicFromLimitValuation)
            True

        """
        FiniteExtensionFromLimitValuation.__init__(self, parent, approximant, G, approximants)
        pAdicValuation_base.__init__(self, parent, approximant.restriction(approximant.domain().base_ring()).p())

    def _to_base_domain(self, f):
        r"""
        Return ``f``, an element of the domain of this valuation, as an element
        of the domain of the underlying limit valuation.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: v = pAdicValuation(GaussianIntegers(), 3)
            sage: I = GaussianIntegers().fraction_field().gen()
            sage: v._to_base_domain(I)
            x

        """
        return f.polynomial()(self._base_valuation.domain().gen())

    def _from_base_domain(self, f):
        r"""
        Return ``f``, an element of the domain of this valuation, as an element
        of the domain of the underlying limit valuation.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: v = pAdicValuation(GaussianIntegers(), 3)
            sage: v._from_base_domain(v._base_valuation.domain().gen())
            I

        """
        return f(self.domain().fraction_field().gen())
