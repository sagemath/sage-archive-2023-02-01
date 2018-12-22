# -*- coding: utf-8 -*-
r"""
Partition Shifting and Raising Operator Algebras

This module contains families of operators that act on partitions or, more
generally, integer sequences. In particular, this includes Young's raising
operators `R_{ij}`, which act on integer sequences by adding `1` to the `i`th
entry and subtracting `1` to the `j`th entry. A special case is acting on
partitions.

AUTHORS:

- Matthew Lancellotti, George H. Seelinger (2018): Initial version

"""

#*****************************************************************************
#  Copyright (C) 2018 Matthew Lancellotti <mvlancellotti@gmail.com>
#                     George H. Seelinger <ghseeli@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.algebras.group_algebra import GroupAlgebra_class
from sage.categories.algebras import Algebras
from sage.categories.groups import Groups
from sage.combinat.composition import Composition
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.partition import Partitions, Partition
from sage.combinat.sf.sf import SymmetricFunctions
from sage.combinat.sf.sfa import SymmetricFunctionAlgebra_generic
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.misc_c import prod
from sage.rings.all import RationalField, NonNegativeIntegerSemiring
from sage.rings.integer_ring import IntegerRing

def free_group_elm_to_partition(elm):
    r"""
    Given an element of an abelian free group indexed by natural numbers,
    return an integer sequence of the powers of the generators.

    INPUT:

    -``elm`` -- an element of an abelian free group indexed by natural numbers

    EXAMPLES::

        sage: from sage.combinat.partition_shifting_algebras import free_group_elm_to_partition
        sage: F = Groups.Commutative().free(NN,'F')
        sage: elm = F.gens()[0]^2*F.gens()[1]^(-3)*F.gens()[5]
        sage: free_group_elm_to_partition(elm)
        (2, -3, 0, 0, 0, 1)
        sage: free_group_elm_to_partition(F.one())
        (0,)
    """
    power_dict = elm.dict()
    max_nonzero_entry = max(power_dict.keys()+[0])
    return tuple(power_dict.get(i,0) for i in range(max_nonzero_entry+1))

def shifting_operator_action_algebra_elm_to_partition_list(elm):
    r"""
    Given an element in the :class:`ShiftingOperatorActionAlgebra`, return
    a list of tuples of the form (powers, coefficient) representing
    inputted element.

    INPUT:

    -``elm`` -- an element in the :class:`ShiftingOperatorActionAlgebra`

    EXAMPLES::

        sage: from sage.combinat.partition_shifting_algebras import ShiftingOperatorActionAlgebra, shifting_operator_action_algebra_elm_to_partition_list
        sage: S = ShiftingOperatorActionAlgebra(ZZ)
        sage: elm = S([2,0,-1,4])+5*S([0,1,0])
        sage: Set(shifting_operator_action_algebra_elm_to_partition_list(elm)) == Set([((2, 0, -1, 4), 1), ((0, 1), 5)])
        True
        sage: shifting_operator_action_algebra_elm_to_partition_list(S.one())
        [((0,), 1)]
    """
    elm_list = list(elm)
    return [(free_group_elm_to_partition(supp),coeff) for (supp,coeff) in elm_list]

class ShiftingOperatorActionAlgebra(GroupAlgebra_class):
    r"""
    A shifting operator action algera.

    This object serves mainly as a superclass to :class:`ShiftingOperatorAlgebra`
    and :class:`RaisingOperatorAlgebra` to act on symmetric functions. In
    general, it is not meant to be used directly by the user. Instead, use
    :class:`ShiftingOperatorAlgebra` or :class:`RaisingOperatorAlgebra`.

    This is an implementation of the ring on which shifting
    operators formally act, and is meant to be isomorphic to
    `R[x_1^\pm, x_2^\pm, x_3^\pm, \ldots]` for some ring `R`, accomplished by
    creating the group algebra of the countably infinite commutative free group
    over ring `R`. In this class, the partition `\lambda=(\lambda_1,\lambda_2,\ldots,
    \lambda_\ell)` is encoded as `x_1^{\lambda_1} x_2^{\lambda_2} \cdots
    x_\ell^{\lambda_\ell}` and this notion generalizes for any sequence of
    integers with finite support. Then, for example, Young's raising operator
    `R_{ij}` is encoded in this space as simply `\frac{x_i}{x_j}`.

    To extend the action of the raising operators on a symmetric function basis,
    one defines an `R`-module homomorphism `\phi` from `R[x_1^\pm, x_2^\pm, \ldots]`
    to the symmetric functions. Then, for a symmetric function basis `b`, we
    compute `R_{ij} b_\lambda` by first computing
    `\frac{x_i}{x_j} x_1^{\lambda_1} \cdots x_\ell^{\lambda_\ell}` and then
    applying `\phi` to the result.

    This is consistent with how many references work formally with raising
    operators. For instance, see exposition surrounding [BMPS2018]_
    Equation (4.1)

    INPUT:

    -``base_ring`` -- the ring of coefficients.

    OPTIONAL ARGUMENTS:

    -``prefix`` -- (default ``"F"``) a label for the basis elements.

    EXAMPLES:

    We initialize the algebra::

        sage: from sage.combinat.partition_shifting_algebras import ShiftingOperatorActionAlgebra
        sage: A = ShiftingOperatorActionAlgebra(QQ); A
        Ring of polynomials in countably infinite variables with prefix F on which shifting operators act over Rational Field

    The algebra also comes equipped with homomorphisms to various
    symmetric function bases. Note, however, not all homomorphisms are
    equivalent. ::

        sage: from sage.combinat.partition_shifting_algebras import ShiftingOperatorActionAlgebra
        sage: A = ShiftingOperatorActionAlgebra(QQ, prefix='x')
        sage: elm = A([3,2,1]); elm
        x0^3*x1^2*x2^1
        sage: sym = SymmetricFunctions(QQ)
        sage: h = sym.h()
        sage: h_mor = h.convert_map_from(A)
        sage: h_mor(elm)
        h[3, 2, 1]
        sage: h(elm)
        h[3, 2, 1]
        sage: s = sym.s()
        sage: s_mor = s.convert_map_from(A)
        sage: s_mor(elm)
        s[3, 2, 1]
        sage: s(elm)
        s[3, 2, 1]
        sage: s_mor(elm) == s(h_mor(elm))
        False
    """
    def __init__(self, base_ring, prefix='F'):
        F = Groups.Commutative().free(NonNegativeIntegerSemiring(),prefix)
        self.prefix = prefix
        category = F.category().Algebras(base_ring)
        GroupAlgebra_class.__init__(self, base_ring, F, category=category)
        self._register_builtin_conversions()

    def _element_constructor_(self, seq):
        r"""
        Construct an element of ``self``.

        INPUT::

        -``seq`` -- Any finite iterable of integers representing the powers of
                   the generators.

        EXAMPLES::

            sage: from sage.combinat.partition_shifting_algebras import ShiftingOperatorActionAlgebra
            sage: A = ShiftingOperatorActionAlgebra(QQ, prefix='x')
            sage: A([1,0,-1])
            x0^1*x2^-1
            sage: 3*A([0,2,0])
            3*x1^2
        """
        F = self.group()
        return self.basis()[prod([F.gens()[i]**(seq[i]) for i in range(len(seq))],F.one())]

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.combinat.partition_shifting_algebras import ShiftingOperatorActionAlgebra
            sage: ShiftingOperatorActionAlgebra(QQ, prefix='x')
            Ring of polynomials in countably infinite variables with prefix x on which shifting operators act over Rational Field
        """
        return "Ring of polynomials in countably infinite variables with prefix " + self.prefix + " on which shifting operators act over " + repr(self.base_ring())

    def from_iterable_indexed_parent(self, iterable_indexed_elm):
        r"""
        Return an element of ``self`` from the underlying iterables of
        ``iterable_indexed_elm``.

        EXAMPLES::

            sage: from sage.combinat.partition_shifting_algebras import ShiftingOperatorActionAlgebra
            sage: S = ShiftingOperatorActionAlgebra(QQ)
            sage: sym = SymmetricFunctions(QQ)
            sage: h = sym.h()
            sage: S.from_iterable_indexed_parent(2*h[3,2,1] + 3*h[2,1]) == 2*S([3,2,1])+3*S([2,1])
            True
        """
        R = iterable_indexed_elm.parent()
        if self.base_ring().has_coerce_map_from(R.base_ring()):
            f = R.module_morphism(self._element_constructor_, codomain=self)
            return f(iterable_indexed_elm)
        else:
            return self._element_constructor_(iterable_indexed_elm)

    def _register_builtin_conversions(self):
        r"""
        Registers known conversions from :class:`ShiftingOperatorActionAlgebra`
        instance ``self`` to known symmetric function bases.

        Currently, the built-in bases for conversion are

        - The complete homogeneous symmetric function basis
        - The Schur symmetric function basis

        ..  WARNING::

            For internal use only during the initialization!

        """
        base = self.base_ring()
        sym = SymmetricFunctions(base)
        h = sym.h()
        h.register_conversion(self.module_morphism(lambda supp: self._supp_to_h(supp,h), codomain=h))
        s = sym.s()
        s.register_conversion(self.module_morphism(lambda supp: self._supp_to_s(supp,s), codomain=s))

    def build_and_register_conversion(self, support_map, codomain):
        r"""
        Builds a module homomorphism from a map sending integer sequences to
        ``codomain`` and registers the result into Sage's conversion model.

        Usually ``codomain`` is a basis of symmetric functions.

        The intended use is to define a morphism from
        ``self`` to a basis of symmetric functions that will be used by
        :class:`ShiftingOperatorAlgebra` or :class:`RaisingOperatorAlgebra` to
        correctly define the action of the operators on the symmetric functions.
        For examples of this, see the ``build_and_register_conversion``
        method of :class:`ShiftingOperatorAlgebra`.

        Note that the actions on the complete homogeneous symmetric functions
        and on the Schur functions by morphisms registered to the
        :class:`ShiftingOperatorActionAlgebra` by :meth:`_register_builtin_conversions`.

        INPUT::

        - ``support_map`` -- A map from integer sequences to ``codomain``.

        - ``codomain``` -- The codomain of ``support_map``, usually a basis of symmetric functions.

        EXAMPLES::

            sage: from sage.combinat.partition_shifting_algebras import ShiftingOperatorActionAlgebra
            sage: S = ShiftingOperatorActionAlgebra(QQ)
            sage: sym = SymmetricFunctions(QQ)
            sage: p = sym.p()
            sage: zero_map = lambda part: p.zero()
            sage: S.build_and_register_conversion(zero_map, p)
            sage: p(2*S([1,0,-1])+S([2,1,0])-3*S([0,1,3])) # indirect doctest
            0


        ..  SEEALSO::
            :class:`ShiftingOperatorAlgebra`, :class:`RaisingOperatorAlgebra`

        ..  WARNING::

            Because :class:`ShiftingOperatorActionAlgebra` ultimately inherits
            from :class:`UniqueRepresentation`, once you register a conversion,
            it will apply to all instances of
            :class:`ShiftingOperatorActionAlgebra` over the same base ring with
            the same prefix::

                sage: from sage.combinat.partition_shifting_algebras import ShiftingOperatorActionAlgebra
                sage: sym = SymmetricFunctions(QQ)
                sage: p = sym.p()
                sage: p._unset_coercions_used()
                sage: zero_map = lambda free_gp_elm: p.zero()
                sage: S = ShiftingOperatorActionAlgebra(QQ)
                sage: S.build_and_register_conversion(zero_map, p)
                sage: p(2*S([1,0,-2])+5*S([3,2])-4*S([2,3]))
                0
                sage: S = ShiftingOperatorActionAlgebra(QQ)
                sage: p(2*S([1,0,-2])+5*S([3,2])-4*S([2,3])) # conversion persists across instances
                0
                sage: St = ShiftingOperatorActionAlgebra(QQ['t'], prefix='x') # but not with different inputs
                sage: p(St([2,3]))
                Traceback (most recent call last):
                ...
                TypeError: do not know how to make x (= x0^2*x1^3) an element of self
        """
        precompose_map = lambda supp: support_map(free_group_elm_to_partition(supp))
        module_morphism = self.module_morphism(precompose_map, codomain=codomain)
        codomain.register_conversion(module_morphism)

    def _supp_to_h(self, supp, basis):
        r"""
        This is a helper funciton that is not meant to be called directly.

        Given the support of an element `x_1^{\gamma_1} x_2^{\gamma_2} \cdots
        x_\ell^{\gamma_\ell}` in the `` ShiftingOperatorActionAlgebra`` and a
        symmetric function algebra basis `b` generated by `\{b_1, b_2, b_3,
        \ldots\}`, return the element `b_{\gamma_1} b_{\gamma_2} \cdots
        b_{\gamma_\ell}` where `b_0 = 1` and `b_{-n} = 0` for all positive
        integers `n`. The canonical example for `b` in this case would be `h`.

        EXAMPLES::

            sage: from sage.combinat.partition_shifting_algebras import ShiftingOperatorActionAlgebra
            sage: S = ShiftingOperatorActionAlgebra(QQ)
            sage: sym = SymmetricFunctions(QQ)
            sage: h = sym.h()
            sage: S._supp_to_h(S([3,2,1]).support_of_term(), h)
            h[3, 2, 1]
            sage: S._supp_to_h(S([2,3,1]).support_of_term(), h)
            h[3, 2, 1]
            sage: S._supp_to_h(S([2,3,-1]).support_of_term(), h) == h.zero()
            True
            sage: S._supp_to_h(S([2,3,0]).support_of_term(), h)
            h[3, 2]
        """
        gamma = list(reversed(sorted(free_group_elm_to_partition(supp))))
        if gamma in Partitions():
            return basis(gamma)
        else:
            return basis.zero()

    def _supp_to_s(self, supp, basis):
        r"""
        This is a helper function that is not meant to be called directly.

        Given the support of an element `x_1^{\gamma_1} x_2^{\gamma_2} \cdots
        x_\ell^{\gamma_\ell}` in the `` ShiftingOperatorActionAlgebra``, return
        the appropriate `s_\gamma` in the Schur basis using
        "Schur function straightening" in [BMPS2018]_ Proposition 4.1.

        EXAMPLES::

            sage: from sage.combinat.partition_shifting_algebras import ShiftingOperatorActionAlgebra
            sage: S = ShiftingOperatorActionAlgebra(QQ)
            sage: sym = SymmetricFunctions(QQ)
            sage: s = sym.s()
            sage: S._supp_to_s(S([3,2,1]).support_of_term(), s)
            s[3, 2, 1]
            sage: S._supp_to_s(S([2,3,1]).support_of_term(), s)
            0
            sage: S._supp_to_s(S([2,4,-1,1]).support_of_term(), s)
            s[3, 3]
            sage: S._supp_to_s(S([3,2,0]).support_of_term(), s)
            s[3, 2]
        """
        def has_nonnegative_parts(lis):
            return all(e >= 0 for e in lis)

        def has_distinct_parts(lis):
            return len(set(lis)) == len(lis)

        def number_of_noninversions(lis):
            num = 0
            for i in range(len(lis)):
                for j in range(i + 1, len(lis)):
                    # i < j is already enforced
                    if lis[i] < lis[j]:
                        num += 1
            return num

        gamma = list(free_group_elm_to_partition(supp))
        rho = list(range(len(gamma) - 1, -1, -1))
        combined = [g + r for g, r in zip(gamma, rho)]
        if has_distinct_parts(combined) and has_nonnegative_parts(combined):
            sign = (-1)**number_of_noninversions(combined)
            sort_combined = reversed(sorted(combined))
            new_gamma = [sc - r for sc, r in zip(sort_combined, rho)]
            return sign * basis(Partition(new_gamma))
        else:
            return basis.zero()

    def _repr_term(self, term):
        r"""
        Return a string representation of ``term``.

        EXAMPLES::

            sage: from sage.combinat.partition_shifting_algebras import ShiftingOperatorActionAlgebra
            sage: S = ShiftingOperatorActionAlgebra(QQ, prefix='x')
            sage: S([1,-2,3]) # indirect doctest
            x0^1*x1^-2*x2^3
        """
        pow_dict = term.dict()
        if len(pow_dict) == 0:
            return '1'
        parts = [self.prefix+repr(k)+'^'+repr(v) for (k,v) in pow_dict.iteritems()]
        return reduce(lambda a,b: a+'*'+b, parts)

class ShiftingSequenceSpace():
    r"""
    A helper for :class:`ShiftingOperatorAlgebra` that contains all
    sequences (represented as tuples) in appropriate base with finite support.

    Helps :class:`ShiftingOperatorAlgebra` know which indices are valid and
    which indices are not for the basis.

    EXAMPLES::

        sage: from sage.combinat.partition_shifting_algebras import ShiftingSequenceSpace
        sage: S = ShiftingSequenceSpace()
        sage: (1, -1) in S
        True
        sage: (1, -1, 0, 9) in S
        True
        sage: [1, -1] in S
        False
        sage: (0.5, 1) in S
        False
    """
    def __init__(self, base=IntegerRing()):
        self.base = base

    def __contains__(self, seq):
        r"""
        Return ``True`` if and only if ``seq`` is a valid shifting sequence.

        EXAMPLES::

            sage: from sage.combinat.partition_shifting_algebras import ShiftingSequenceSpace
            sage: S = ShiftingSequenceSpace()
            sage: (1, -1) in S
            True
            sage: (1, -1, 0, 9) in S
            True
            sage: [1, -1] in S
            False
            sage: (0.5, 1) in S
            False
        """
        if not isinstance(seq, tuple):
            return False
        return not any(i not in self.base for i in seq)

    CHECK_ERROR_MESSAGE = 'Expected valid index (a tuple of {base}), but instead received {seq}.'

    def check(self, seq):
        r"""
        Verify that ``seq`` is a valid shifting sequence.

        If it is not, raise a ``ValueError``.

        EXAMPLES::

            sage: from sage.combinat.partition_shifting_algebras import ShiftingSequenceSpace
            sage: S = ShiftingSequenceSpace()
            sage: S.check((1, -1))
            sage: S.check((1, -1, 0, 9))
            sage: S.check([1, -1])
            Traceback (most recent call last):
            ...
            ValueError: Expected valid index (a tuple of Integer Ring), but instead received [1, -1].
            sage: S.check((0.5, 1))
            Traceback (most recent call last):
            ...
            ValueError: Expected valid index (a tuple of Integer Ring), but instead received (0.500000000000000, 1).
        """
        if not self.__contains__(seq):
            raise ValueError(self.CHECK_ERROR_MESSAGE.format(
                base=self.base, seq=seq))


class RaisingSequenceSpace(ShiftingSequenceSpace):
    r"""
    A helper for :class:`RaisingOperatorAlgebra` containing all integer
    sequences of finite support that sum to zero.

    Helps :class:`RaisingOperatorAlgebra` know which indices are valid and which
    indices are not for the basis.

    EXAMPLES::

        sage: from sage.combinat.partition_shifting_algebras import RaisingSequenceSpace
        sage: RS = RaisingSequenceSpace()
        sage: (1, -1) in RS
        True
        sage: (1, 0, -1) in RS
        True
        sage: (1, -1, 0, 9) in RS
        False
        sage: [1, -1] in RS
        False
    """

    CHECK_ERROR_MESSAGE = 'Expected valid index (a tuple of {base} elements, where every partial sum is nonnegative and every total sum is 0), but instead received {seq}.'

    def __contains__(self, seq):
        r"""
        Return ``True`` if and only if ``seq`` is a valid raising sequence.

        EXAMPLES::

            sage: from sage.combinat.partition_shifting_algebras import RaisingSequenceSpace
            sage: RS = RaisingSequenceSpace()
            sage: (1, -1) in RS
            True
            sage: (1, 0, -1) in RS
            True
            sage: (1, -1, 0, 9) in RS
            False
            sage: [1, -1] in RS
            False
        """
        # check that it is a shifting sequence
        if not ShiftingSequenceSpace.__contains__(self, seq):
            return False
        # check that every partial sum is nonnegative
        partial_sum = 0
        for term in seq:
            partial_sum += term
            if partial_sum < 0:
                return False
        # check that total sum is 0
        if partial_sum != 0:
            return False
        # finally, succeed
        return True


class ShiftingOperatorAlgebra(ShiftingOperatorActionAlgebra):
    r"""
    An algebra of shifting operators.

    We follow the following convention:

    ``S[(1, 0, -1, 2)]`` is the shifting operator that raises the first part by
    1, lowers the third part by 1, and raises the fourth part by 2.

    In addition to acting on partitions (or any integer sequence), the shifting
    operators can also act on
    symmetric functions in a basis `b` when a conversion from
    :class:`ShiftingOperatorActionAlgebra` to `b` has been registered,
    preferably using :meth:`build_and_register_conversion`.

    OPTIONAL ARGUMENTS:

    - ``base_ring`` -- (default ``RationalField()['t']``) the ring you will use on the raising operators.

    - ``prefix`` -- (default ``"S"``) the label for the shifting operators.

    EXAMPLES::

        sage: S = ShiftingOperatorAlgebra()
        sage: s = SymmetricFunctions(RationalField()['t']).s()
        sage: h = SymmetricFunctions(RationalField()['t']).h()

        sage: S[(1, -1, 2)]
        S(1, -1, 2)
        sage: S[(1, -1, 2)]([5, 4])
        [([6, 3, 2], 1)]
        sage: S[(1, -1, 2)](s[5, 4])
        s[6, 3, 2]
        sage: S[(1, -1, 2)](h[5, 4])
        h[6, 3, 2]

        sage: (1 - S[(1,-1)]) * (1 - S[(4,)]) == S([]) - S([1, -1]) - S([4]) + S([5, -1])
        True
        sage: ((1 - S[(1,-1)]) * (1 - S[(4,)]))(s[2, 2, 1])
        s[2, 2, 1] - s[3, 1, 1] - s[6, 2, 1] + s[7, 1, 1]

    ..  SEEALSO::
        :class:`RaisingOperatorAlgebra`
    """
    def __init__(self, base_ring=RationalField()['t'], prefix='S', basis_indices=ShiftingSequenceSpace()):
        self._prefix = prefix
        self._base_ring = base_ring
        # a single basis index looks like (1, 0, -1, 2), for example
        self._basis_indices = basis_indices
        ShiftingOperatorActionAlgebra.__init__(
            self,
            self._base_ring,
            prefix=self._prefix)

    def __getitem__(self, seq):
        r"""
        Return the shifting operator whose index is ``seq``.

        This method is only for basis indices.

        EXAMPLES::

            sage: S = ShiftingOperatorAlgebra()
            sage: S[1, 1, -9]
            S(1, 1, -9)
        """
        return ShiftingOperatorActionAlgebra._element_constructor_(self, seq)

    def _element_constructor_(self, seq):
        r"""
        Return the shifting operator whose index is ``seq``.

        This method is only for basis indices.

        EXAMPLES::

            sage: S = ShiftingOperatorAlgebra()
            sage: S._element_constructor_([1, 1, -9])
            S(1, 1, -9)
            sage: S[1, 1, -9]
            S(1, 1, -9)
        """
        seq = tuple(seq)
        self._basis_indices.check(seq)
        return ShiftingOperatorActionAlgebra._element_constructor_(self, seq)

    def _repr_(self):
        r"""
        Return a string describing ``self``.

        EXAMPLES::

            sage: S = ShiftingOperatorAlgebra()
            sage: S
            Shifting Operator Algebra over Univariate Polynomial Ring in t over Rational Field
        """
        return "Shifting Operator Algebra over {base_ring}".format(base_ring=self._base_ring)

    def _repr_term(self, term):
        r"""
        Return a string representation of the shifting operator.

        EXAMPLES::

            sage: S = ShiftingOperatorAlgebra(QQ)
            sage: S[1, 1, -9] # indirect doctest
            S(1, 1, -9)
            sage: S[1, 1, 0, -9, 0] # indirect doctest
            S(1, 1, 0, -9)
            sage: S.one() # indirect doctest
            S()
        """
        if term == self.group().one():
            return self.prefix + '()'
        else:
            return self.prefix + repr(free_group_elm_to_partition(term))

    def ambient(self):
        r"""
        Return the ambient :class:`ShiftingSequenceActionAlgebra` in which
        the :class:`ShiftingOperatorAlgebra` ``self`` is a subalgebra.

        EXAMPLES::

            sage: S = ShiftingOperatorAlgebra()
            sage: S.ambient()
            Ring of polynomials in countably infinite variables with prefix S on which shifting operators act over Univariate Polynomial Ring in t over Rational Field
        """
        return self.lift.codomain()

    @lazy_attribute
    def lift(self):
        r"""
        Return the lift map from the :class:`ShiftingOperatorAlgebra` ``self``
        to the ambient :class:`ShiftingSequenceActionAlgebra`.

        EXAMPLES::

            sage: S = ShiftingOperatorAlgebra()
            sage: elm = S[1, 1, -9]
            sage: lifted = S.lift(elm); lifted
            S0^1*S1^1*S2^-9
            sage: lifted.parent() is S.ambient()
            True
        """
        amb = ShiftingOperatorActionAlgebra(self._base_ring, prefix=self._prefix)
        phi = self.module_morphism(amb.monomial, codomain=amb)
        phi.register_as_coercion()
        return phi

    def build_and_register_conversion(self, support_map, codomain):
        r"""
        Builds a module homomorphism from a map sending integer sequences to
        ``codomain`` and registers the result into Sage's conversion model.

        EXAMPLES::

            sage: S = ShiftingOperatorAlgebra(QQ)
            sage: sym = SymmetricFunctions(QQ)
            sage: p = sym.p()
            sage: zero_map = lambda part: p.zero()
            sage: S.build_and_register_conversion(zero_map, p)
            sage: op = S((1, -1))
            sage: op(2*p[4,3]+5*p[2,2]+7*p[2]) == p.zero() # indirect doctest
            True

        For a more illustrative example, we can implement a simple
        (but not mathematically justified!) conversion on the monomial basis.::

            sage: S = RaisingOperatorAlgebra(QQ)
            sage: sym = SymmetricFunctions(QQ)
            sage: m = sym.m()
            sage: supp_map = lambda gamma: m(list(reversed(sorted(gamma)))) if list(reversed(sorted(gamma))) in Partitions() else m.zero()
            sage: S.build_and_register_conversion(supp_map, m)
            sage: ROA = RaisingOperatorAlgebra(QQ)
            sage: op = ROA.ij(0,1)
            sage: op(2*m[4,3]+5*m[2,2]+7*m[2]) == 2*m[5, 2] + 5*m[3, 1] # indirect doctest
            True

        ..  SEEALSO::
            :class:`RaisingOperatorAlgebra` and the
            ``build_and_register_conversion`` method of
            :class:`ShiftingOperatorActionAlgebra`
        """
        A = self.ambient()
        A.build_and_register_conversion(support_map, codomain)

    class Element(ShiftingOperatorActionAlgebra.Element):
        r"""
        An element of a :class`ShiftingOperatorAlgebra`.
        """

        def indices(self):
            r"""
            Return the indices of ``self`` since the support is a free group
            element but it is more natural to work with the powers.

            EXAMPLES::

                sage: S = ShiftingOperatorAlgebra()
                sage: (S[2, 1] + S[1, 1]).indices()
                [(1, 1), (2, 1)]
            """
            return [free_group_elm_to_partition(supp) for supp in self.support()]

        def index(self):
            r"""
            Return the index of ``self`` since the support is a free group
            element but it is more natural to work with the powers.

            This method is only for basis elements.

            EXAMPLES::

                sage: S = ShiftingOperatorAlgebra()
                sage: S[2, 1].index()
                (2, 1)
                sage: (S[2, 1] + S[1, 1]).index()
                Traceback (most recent call last):
                ...
                ValueError: This is only defined for basis elements.  For other elements, use indices() instead.
         """
            if len(self) != 1:
                raise ValueError(
                    "This is only defined for basis elements.  For other elements, use indices() instead.")
            return self.indices()[0]

        @staticmethod
        def _call_basis_on_index(seq, index):
            r"""
            For internal use only!

            Return the action of the basis element indexed by ``seq`` upon the composition ``index``.

            INPUTS:

            - ``seq`` -- The sequence of the basis element that acts.

            - ``index`` -- A sequence (typically a composition or a partition) that we act upon.

            EXAMPLES::

                sage: S = ShiftingOperatorAlgebra()
                sage: S[2, 1]._call_basis_on_index([1, 1], [1, 2, 3, 4, 5])
                [2, 3, 3, 4, 5]

            ..  SEEALSO::
                :meth:`__call__`
            """
            # pad sequence and index with 0's
            index = list(index) + [0] * (len(seq) - len(index))
            seq = tuple(seq) + (0,) * (len(index) - len(seq))
            # raise and drop
            return [v + s for v, s in zip(index, seq)]

        def __call__(self, operand):
            r"""
            Call method for shifting sequence operators to act on objects.

            EXAMPLES::

                sage: S = ShiftingOperatorAlgebra(QQ)
                sage: op = S([1,1]) + 2*S([0,1,0,1])
                sage: Set(op([1,1,1,1])) == Set([([2,2,1,1], 1), ([1,2,1,2], 2)])
                True
                sage: Set(op(Partition([1,1,1,1]))) == Set([([2,2,1,1], 1), ([1,2,1,2], 2)])
                True
                sage: sym = SymmetricFunctions(QQ)
                sage: h = sym.h()
                sage: op(h[1,1,1,1])
                3*h[2, 2, 1, 1]
                sage: s = sym.s()
                sage: op(s[1,1,1,1])
                s[2, 2, 1, 1]
                sage: m = sym.m()
                sage: Set(op(m[1,1,1,1])) == Set([((2,2,1,1), 1), ((1,2,1,2), 2)])
                True
            """
            A = self.parent().ambient()
            parent = self.parent()
            if isinstance(operand, (list, tuple, Composition, Partition)):
                self_terms = [(free_group_elm_to_partition(supp),coeff) for (supp,coeff) in self]
                return [(self._call_basis_on_index(index,operand),coeff) for (index,coeff) in self_terms]
            if hasattr(operand,"parent"):
                lift_operand = A.from_iterable_indexed_parent(operand)
                result = A(self)*lift_operand
                operand_parent = operand.parent()
                try:
                    return operand_parent(result)
                except TypeError:
                    return shifting_operator_action_algebra_elm_to_partition_list(result)
            else:
                return self(Composition(operand))

class RaisingOperatorAlgebra(ShiftingOperatorAlgebra):
    r"""
    An algebra of raising operators.

    This class subclasses :class:`ShiftingOperatorAlgebra` and inherits the
    large majority of its functionality from there.

    We follow the following convention!:

    ``R[(1, 0, -1)]`` is the raising operator that raises the first part by 1
    and lowers the third part by 1.

    For a definition of raising operators, see [cat]_ Definition 2.1, but be
    wary that the notation is different there.  See :meth:`ij` for a way to
    create operators using the notation in the paper.

    If you do NOT want any restrictions on the allowed sequences, use
    :class:`ShiftingOperatorAlgebra` instead of :class:`RaisingOperatorAlgebra`.

    OPTIONAL ARGUMENTS:

    - ``base_ring`` -- (default ``RationalField()['t']``) the ring you will use on the raising operators.

    - ``prefix`` -- (default ``"R"``) the label for the raising operators.

    EXAMPLES::

        sage: R = RaisingOperatorAlgebra()
        sage: s = SymmetricFunctions(RationalField()['t']).s()
        sage: h = SymmetricFunctions(RationalField()['t']).h()

        sage: R[(1, -1)]
        R(1, -1)
        sage: R[(1, -1)](s[5, 4])
        s[6, 3]
        sage: R[(1, -1)](h[5, 4])
        h[6, 3]

        sage: (1 - R[(1,-1)]) * (1 - R[(0,1,-1)]) == 1 - R[(0, 1, -1)] - R[(1, -1)] + R[(1, 0, -1)]
        True
        sage: ((1 - R[(1,-1)]) * (1 - R[(0,1,-1)]))(s[2, 2, 1])
        s[2, 2, 1] - s[3, 1, 1] + s[3, 2]

    We can also use raising operators to implement the Jacobi-Trudi identity::

        sage: op = (1-R[(1,-1)])*(1-R[(1,0,-1)])*(1-R[(0,1,-1)])
        sage: s(op(h[3,2,1]))
        s[3, 2, 1]

    ..  SEEALSO::
        :class:`ShiftingOperatorAlgebra`
    """

    def __init__(self, base_ring=RationalField()['t'], prefix='R'):
        ShiftingOperatorAlgebra.__init__(self,
                                         base_ring=base_ring,
                                         prefix=prefix,
                                         basis_indices=RaisingSequenceSpace())

    def ij(self, i, j):
        r"""
        Return the raising operator `R_{ij}` as notated in [BMPS2018]_ Definition 2.1.

        Shorthand element constructor that allows you to create raising
        operators using the familiar `R_{ij}` notation found in
        [BMPS2018]_ Definition 2.1, with the exception that indices here are 0-based,
        not 1-based.

        EXAMPLES:

        Create the raising operator which raises part 0 and lowers part 2
        (indices are 0-based)::

            sage: R = RaisingOperatorAlgebra()
            sage: R.ij(0, 2)
            R(1, 0, -1)

        ..  SEEALSO::
            :meth:`ShiftingOperatorAlgebra._element_constructor_`, :meth:`ShiftingOperatorAlgebra.__getitem__`
        """
        if not i in NonNegativeIntegerSemiring():
            raise ValueError(
                'i must be a natural number.  You input i = {i}.'.format(i=i))
        if not j in NonNegativeIntegerSemiring():
            raise ValueError(
                'j must be a natural number.  You input j = {j}.'.format(j=j))
        if not i < j:
            raise ValueError(
                'Index j must be greater than index i.  You input (i, j) = ({i}, {j}).'.format(i=i, j=j))
        seq = [0] * (max(i, j) + 1)
        seq[i] = 1
        seq[j] = -1
        seq = tuple(seq)
        return self._element_constructor_(seq)
