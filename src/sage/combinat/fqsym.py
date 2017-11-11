# -*- coding: utf-8 -*-
r"""
Free Quasi-symmetric functions

AUTHORS:

- Frédéric Chapoton, Darij Grinberg (2017)
"""

# ****************************************************************************
#       Copyright (C) 2010-2015 Frédéric Chapoton <chapoton@unistra.fr>,
#       Copyright (C) 2017      Darij Grinberg <dgrinber at umn.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.bindable_class import BindableClass
from sage.misc.lazy_attribute import lazy_attribute
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.hopf_algebras import HopfAlgebras
from sage.categories.realizations import Realizations, Category_realization_of_parent
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.permutation import Permutations, Permutation
from sage.groups.perm_gps.permgroup_element import PermutationGroupElement
from sage.combinat.words.word import Word
from sage.combinat.symmetric_group_algebra import SymmetricGroupAlgebra

class FQSymBasis_abstract(CombinatorialFreeModule, BindableClass):
    """
    Abstract base class for bases of FQSym.

    This must define two attributes:

    - ``_prefix`` -- the basis prefix
    - ``_basis_name`` -- the name of the basis and must match one
      of the names that the basis can be constructed from FQSym
    """
    def __init__(self, alg):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: TestSuite(algebras.FQSym(QQ).F()).run()
        """
        CombinatorialFreeModule.__init__(self, alg.base_ring(),
                                         Permutations(),
                                         category=FQSymBases(alg),
                                         bracket="", prefix=self._prefix)

    def _coerce_map_from_(self, R):
        r"""
        Return ``True`` if there is a coercion from ``R`` into ``self``
        and ``False`` otherwise.

        The things that coerce into ``self`` are

        - free quasi-symmetric functions over a base with
          a coercion map into ``self.base_ring()``

        EXAMPLES::

            sage: F = algebras.FQSym(GF(7)).F(); F
            Free Quasi-symmetric functions over Finite Field of size 7 in the F basis

        Elements of the free quasi-symmetric functions canonically coerce in::

            sage: x, y, z = F([1]), F([2,1]), F([1,3,2])
            sage: F.coerce(x+y) == x+y
            True

        The free quasi-symmetric functions over `\ZZ` coerces in,
        since `\ZZ` coerces to `\GF{7}`::

            sage: G = algebras.FQSym(ZZ).F()
            sage: Gx, Gy = G([1]), G([2,1])
            sage: z = F.coerce(Gx+Gy); z
            F[1] + F[2, 1]
            sage: z.parent() is F
            True

        However, `\GF{7}` does not coerce to `\ZZ`, so free
        quasi-symmetric functions over `\GF{7}` does not coerce
        to the same algebra over `\ZZ`::

            sage: G.coerce(y)
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Free Quasi-symmetric functions
             over Finite Field of size 7 in the F basis to
             Free Quasi-symmetric functions over Integer Ring in the F basis

        TESTS::

            sage: F = algebras.FQSym(ZZ).F()
            sage: G = algebras.FQSym(QQ).F()
            sage: F.has_coerce_map_from(G)
            False
            sage: G.has_coerce_map_from(F)
            True
            sage: F.has_coerce_map_from(QQ)
            False
            sage: G.has_coerce_map_from(QQ)
            True
            sage: F.has_coerce_map_from(PolynomialRing(ZZ, 3, 'x,y,z'))
            False
        """
        # free quasi-symmetric functions in the same variables
        # over any base that coerces in:
        if isinstance(R, FQSymBasis_abstract):
            if R.realization_of() == self.realization_of():
                return True
            if not self.base_ring().has_coerce_map_from(R.base_ring()):
                return False
            if self._basis_name == R._basis_name: # The same basis
                def coerce_base_ring(self, x):
                    return self._from_dict(x.monomial_coefficients())
                return coerce_base_ring
            # Otherwise lift that basis up and then coerce over
            target = getattr(self.realization_of(), R._basis_name)()
            return self._coerce_map_via([target], R)
        return super(FQSymBasis_abstract, self)._coerce_map_from_(R)

class FreeQuasisymmetricFunctions(UniqueRepresentation, Parent):
    r"""
    The free quasi-symmetric functions.

    The Hopf algebra `FQSym` of free quasi-symmetric functions
    over a commutative ring `R` is the free `R`-module with basis
    indexed by all permutations (i.e., the indexing set is
    the disjoint union of all symmetric groups).
    Its product is determined by the shifted shuffles of two
    permutations, whereas its coproduct is given by splitting
    a permutation (regarded as a word) into two (at every
    possible point) and standardizing the two pieces.
    This Hopf algebra was introduced in [MR]_.
    See [GriRei16]_ (Chapter 8) for a treatment using modern
    notations.

    In more detail:
    For each `n \geq 0`, consider the symmetric group `S_n`.
    Let `S` be the disjoint union of the `S_n` over all
    `n \geq 0`.
    Then, `FQSym` is the free `R`-module with basis
    `(F_w)_{w \in S}`.
    This `R`-module is graded, with the `n`-th graded
    component being spanned by all `F_w` for `w \in S_n`.
    A multiplication is defined on `FQSym` as follows:
    For any two permutations `u \in S_k` and `v \in S_l`,
    we set

    .. MATH::

        F_u F_v = \sum F_w ,

    where the sum is over all shuffles of `u` with `v[k]`.
    Here, the permutations `u` and `v` are regarded as words
    (by writing them in one-line notation), and `v[k]` means
    the word obtained from `v` by increasing each letter by
    `k` (for example, `(1,4,2,3)[5] = (6,9,7,8)`); and the
    shuffles `w` are translated back into permutations.
    This defines an associative multiplication on `FQSym`;
    its unity is `F_e`, where `e` is the identity
    permutation in `S_0`.

    As an associative algebra, `FQSym` has the richer structure
    of a dendriform algebra. This means that the associative
    product ``*`` is decomposed as a sum of two binary operations

    .. MATH::

        x y = x \succ y + x \prec y

    that satisfy the axioms:

    .. MATH::

        (x \succ y) \prec z = x \succ (y \prec z),

    .. MATH::

        (x \prec y) \prec z = x \prec (y z),

    .. MATH::

        (x y) \succ z = x \succ (y \succ z).

    These two binary operations are defined similarly to the
    (associative) product above: We set

    .. MATH::

        F_u \prec F_v = \sum F_w ,

    where the sum is now over all shuffles of `u` with `v[k]`
    whose first letter is taken from `u` (rather than from
    `v[k]`). Similarly,

    .. MATH::

        F_u \succ F_v = \sum F_w ,

    where the sum is over all remaining shuffles of `u` with
    `v[k]`.

    .. TODO::

        Explain what `1 \prec 1` and `1 \succ 1` are.

    .. TODO::

        Doctest all 6 possibilities involving `1` on one
        side of a `\prec` or `\succ`.

    .. NOTE::

        The usual binary operator ``*`` is used for the
        associative product.

    EXAMPLES::

        sage: F = algebras.FQSym(ZZ).F()
        sage: x,y,z = F([1]), F([1,2]), F([1,3,2])
        sage: (x * y) * z
        F[1, 2, 3, 4, 6, 5] + ...

    The product of `FQSym` is associative::

        sage: x * (y * z) == (x * y) * z
        True

    The associative product decomposes into two parts::

        sage: x * y == F.prec(x, y) + F.succ(x, y)
        True

    The axioms of dendriform algebra hold::

        sage: F.prec(F.succ(x, y), z) == F.succ(x, F.prec(y, z))
        True
        sage: F.prec(F.prec(x, y), z) == F.prec(x, y * z)
        True
        sage: F.succ(x * y, z) == F.succ(x, F.succ(y, z))
        True

    `FQSym` is also known as the Malvenuto-Reutenauer algebra::

        sage: algebras.MalvenutoReutenauer(ZZ)
        Free Quasi-symmetric functions over Integer Ring

    REFERENCES:

    - [MR]_
    - [LodayRonco]_
    - [GriRei16]_
    """

    def __init__(self, R):
        """
        Initialize ``self``.

        TESTS::

            sage: A = algebras.FQSym(QQ); A
            Free Quasi-symmetric functions over Rational Field
            sage: TestSuite(A).run()  # long time (3s)

            sage: F = algebras.FQSym(QQ)
            sage: TestSuite(F).run() # long time (3s)
        """
        self._category = HopfAlgebras(R).Graded().Connected()
        Parent.__init__(self, base=R, category=self._category.WithRealizations())

    def _repr_(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: algebras.FQSym(QQ)  # indirect doctest
            Free Quasi-symmetric functions over Rational Field
        """
        s = "Free Quasi-symmetric functions over {}"
        return s.format(self.base_ring())

    def a_realization(self):
        r"""
        Return a particular realization of ``self`` (the F-basis).

        EXAMPLES::

            sage: FQSym = algebras.FQSym(QQ)
            sage: FQSym.a_realization()
            Free Quasi-symmetric functions over Rational Field in the F basis
        """
        return self.F()

    class F(FQSymBasis_abstract):
        """
        The F-basis of `FQSym`.

        This is the basis `(F_w)`, with `w` ranging over all
        permutations. See the documentation of :class:`FQSym`
        for details.

        EXAMPLES::

            sage: FQSym = algebras.FQSym(QQ)
            sage: FQSym.F()
            Free Quasi-symmetric functions over Rational Field in the F basis
        """
        _prefix = "F"
        _basis_name = "F"

        def _element_constructor_(self, x):
            r"""
            Convert ``x`` into ``self``.

            EXAMPLES::

                sage: R = algebras.FQSym(QQ).F()
                sage: x, y, z = R([1]), R([2,1]), R([3,2,1])
                sage: R(x)
                F[1]
                sage: R(x+4*y)
                F[1] + 4*F[2, 1]
                sage: R(1)
                F[]

                sage: D = algebras.FQSym(ZZ).F()
                sage: X, Y, Z = D([1]), D([2,1]), D([3,2,1])
                sage: R(X-Y).parent()
                Free Quasi-symmetric functions over Rational Field in the F basis

                sage: R([1, 3, 2])
                F[1, 3, 2]
                sage: R(Permutation([1, 3, 2]))
                F[1, 3, 2]
                sage: R(SymmetricGroup(4)(Permutation([1,3,4,2])))
                F[1, 3, 4, 2]
            """
            if isinstance(x, (list, tuple, PermutationGroupElement)):
                x = Permutation(x)
            try:
                P = x.parent()
                if isinstance(P, FreeQuasisymmetricFunctions.F):
                    if P is self:
                        return x
                    return self.element_class(self, x.monomial_coefficients())
            except AttributeError:
                pass
            return CombinatorialFreeModule._element_constructor_(self, x)

        def degree_on_basis(self, t):
            """
            Return the degree of a permutation in
            the algebra of free quasi-symmetric functions.

            This is the size of the permutation (i.e., the `n`
            for which the permutation belongs to `S_n`).

            EXAMPLES::

                sage: A = algebras.FQSym(QQ).F()
                sage: u = Permutation([2,1])
                sage: A.degree_on_basis(u)
                2
            """
            return len(t)

        @cached_method
        def an_element(self):
            """
            Return an element of ``self``.

            EXAMPLES::

                sage: A = algebras.FQSym(QQ).F()
                sage: A.an_element()
                F[1] + 2*F[1, 2] + 2*F[2, 1]
            """
            o = self([1])
            return o + 2 * o * o

        def some_elements(self):
            """
            Return some elements of the free quasi-symmetric functions.

            EXAMPLES::

                sage: A = algebras.FQSym(QQ).F()
                sage: A.some_elements()
                [F[], F[1], F[1, 2] + F[2, 1],
                 F[] + F[1, 2] + F[2, 1]]
            """
            u = self.one()
            o = self([1])
            x = o * o
            y = u + x
            return [u, o, x, y]

        def one_basis(self):
            """
            Return the index of the unit.

            EXAMPLES::

                sage: A = algebras.FQSym(QQ).F()
                sage: A.one_basis()
                []
            """
            Perm = self.basis().keys()
            return Perm([])

        def product_on_basis(self, x, y):
            r"""
            Return the `*` associative product of two permutations.

            This is the shifted shuffle of `x` and `y`.

            .. SEEALSO::

                :meth:`succ_product_on_basis`, :meth:`prec_product_on_basis`

            EXAMPLES::

                sage: A = algebras.FQSym(QQ).F()
                sage: x = Permutation([1])
                sage: A.product_on_basis(x, x)
                F[1, 2] + F[2, 1]
            """
            n = len(x)
            basis = self.basis()
            return self.sum(basis[u] for u in x.shifted_shuffle(y))

        def succ_product_on_basis(self, x, y):
            r"""
            Return the `\succ` product of two permutations.

            This is the shifted shuffle of `x` and `y` with the additional
            condition that the first letter of the result comes from `y`.

            The usual symbol for this operation is `\succ`.

            .. SEEALSO::

                - :meth:`product_on_basis`, :meth:`prec_product_on_basis`

            EXAMPLES::

                sage: A = algebras.FQSym(QQ).F()
                sage: x = Permutation([1,2])
                sage: A.succ_product_on_basis(x, x)
                F[3, 1, 2, 4] + F[3, 1, 4, 2] + F[3, 4, 1, 2]
                sage: y = Permutation([])
                sage: A.succ_product_on_basis(x, y) == 0
                True
                sage: A.succ_product_on_basis(y, x) == A(x)
                True

            TESTS::

                sage: u = A.one().support()[0]
                sage: A.succ_product_on_basis(u, u)
                Traceback (most recent call last):
                ...
                ValueError: products | < | and | > | are not defined
            """
            if not y:
                if not x:
                    raise ValueError("products | < | and | > | are not defined")
                else:
                    return self.zero()
            basis = self.basis()
            if not x:
                return basis[y]
            K = basis.keys()
            n = len(x)
            shy = Word([a + n for a in y])
            shy0 = shy[0]
            return self.sum(basis[K([shy0] + list(u))]
                            for u in Word(x).shuffle(Word(shy[1:])))

        def prec_product_on_basis(self, x, y):
            r"""
            Return the `\prec` product of two permutations.

            This is the shifted shuffle of `x` and `y` with the additional
            condition that the first letter of the result comes from `x`.

            The usual symbol for this operation is `\prec`.

            .. SEEALSO::

                :meth:`product_on_basis`, :meth:`succ_product_on_basis`

            EXAMPLES::

                sage: A = algebras.FQSym(QQ).F()
                sage: x = Permutation([1,2])
                sage: A.prec_product_on_basis(x, x)
                F[1, 2, 3, 4] + F[1, 3, 2, 4] + F[1, 3, 4, 2]
                sage: y = Permutation([])
                sage: A.prec_product_on_basis(x, y) == A(x)
                True
                sage: A.prec_product_on_basis(y, x) == 0
                True

            TESTS::

                sage: u = A.one().support()[0]
                sage: A.prec_product_on_basis(u, u)
                Traceback (most recent call last):
                ...
                ValueError: products | < | and | > | are not defined
            """
            if not x and not y:
                raise ValueError("products | < | and | > | are not defined")
            if not x:
                return self.zero()
            basis = self.basis()
            if not y:
                return basis[x]
            K = basis.keys()
            n = len(x)
            shy = Word([a + n for a in y])
            x0 = x[0]
            return self.sum(basis[K([x0] + list(u))]
                            for u in Word(x[1:]).shuffle(shy))

        def coproduct_on_basis(self, x):
            r"""
            Return the coproduct of `F_{\sigma}` for `\sigma` a permutation.

            EXAMPLES::

                sage: A = algebras.FQSym(QQ).F()
                sage: x = A([1])
                sage: ascii_art(A.coproduct(A.one()))  # indirect doctest
                1 # 1

                sage: ascii_art(A.coproduct(x))  # indirect doctest
                1 # F    + F    # 1
                     [1]    [1]

                sage: A = algebras.FQSym(QQ).F()
                sage: x, y, z = A([1]), A([2,1]), A([3,2,1])
                sage: A.coproduct(z)
                F[] # F[3, 2, 1] + F[1] # F[2, 1] + F[2, 1] # F[1]
                + F[3, 2, 1] # F[]
            """
            if not len(x):
                return self.one().tensor(self.one())
            return sum(self(Word(x[:i]).standard_permutation()).tensor(
                                self(Word(x[i:]).standard_permutation()))
                        for i in range(len(x) + 1))

class FQSymBases(Category_realization_of_parent):
    r"""
    The category of bases of `FQSym`.
    """
    def __init__(self, base):
        r"""
        Initialize the bases of an `FQSym`

        INPUT:

        - ``base`` -- an instance of `FQSym`

        TESTS::

            sage: from sage.combinat.fqsym import FQSymBases
            sage: FQSym = algebras.FQSym(ZZ)
            sage: bases = FQSymBases(FQSym)
            sage: FQSym.F() in bases
            True
        """
        Category_realization_of_parent.__init__(self, base)

    def _repr_(self):
        r"""
        Return the representation of ``self``.

        EXAMPLES::

            sage: from sage.combinat.fqsym import FQSymBases
            sage: FQSym = algebras.FQSym(ZZ)
            sage: FQSymBases(FQSym)
            Category of bases of Free Quasi-symmetric functions over Integer Ring
        """
        return "Category of bases of {}".format(self.base())

    def super_categories(self):
        r"""
        The super categories of ``self``.

        EXAMPLES::

            sage: from sage.combinat.fqsym import FQSymBases
            sage: FQSym = algebras.FQSym(ZZ)
            sage: bases = FQSymBases(FQSym)
            sage: bases.super_categories()
            [Category of graded connected hopf algebras with basis over Integer Ring,
             Category of realizations of Free Quasi-symmetric functions over Integer Ring]
        """
        return [self.base()._category.WithBasis().Graded(),
                Realizations(self.base())]

    class ParentMethods:
        def _repr_(self):
            """
            Text representation of this basis of `FQSym`.

            EXAMPLES::

                sage: FQSym = algebras.FQSym(ZZ)
                sage: FQSym.F()
                Free Quasi-symmetric functions over Integer Ring in the F basis
            """
            return "{} in the {} basis".format(self.realization_of(), self._basis_name)

        def __getitem__(self, p):
            """
            Return the basis element indexed by ``p``.

            INPUT:

            - ``p`` -- a permutation

            EXAMPLES::

                sage: R = algebras.FQSym(QQ).F()
                sage: R[[1, 3, 2]]
                F[1, 3, 2]
                sage: R[Permutation([1, 3, 2])]
                F[1, 3, 2]
                sage: R[SymmetricGroup(4)(Permutation([1,3,4,2]))]
                F[1, 3, 4, 2]
            """
            return self.monomial(Permutation(p))

        def is_field(self, proof=True):
            """
            Return whether this `FQSym` is a field.

            EXAMPLES::

                sage: F = algebras.FQSym(QQ).F()
                sage: F.is_field()
                False
            """
            return False

        def is_commutative(self):
            """
            Return whether this `FQSym` is commutative.

            EXAMPLES::

                sage: F = algebras.FQSym(ZZ).F()
                sage: F.is_commutative()
                False
            """
            return self.base_ring().is_zero()

        @lazy_attribute
        def succ(self):
            r"""
            Return the `\succ` product.

            This is the shifted shuffle of `x` and `y` with the additional
            condition that the first letter of the result comes from `y`.

            The usual symbol for this operation is `\succ`.

            .. SEEALSO::

                :meth:`product`, :meth:`prec`, :meth:`over`, :meth:`under`

            EXAMPLES::

                sage: A = algebras.FQSym(QQ).F()
                sage: x = A([1])
                sage: A.succ(x, x)
                F[2, 1]
                sage: y = A([3,1,2])
                sage: A.succ(x, y)
                F[4, 1, 2, 3] + F[4, 2, 1, 3] + F[4, 2, 3, 1]
                sage: A.succ(y, x)
                F[4, 3, 1, 2]
            """
            suc = self.succ_product_on_basis
            return self._module_morphism(self._module_morphism(suc, position=0,
                                                               codomain=self),
                                         position=1)

        @lazy_attribute
        def prec(self):
            r"""
            Return the `\prec` product.

            This is the shifted shuffle of `x` and `y` with the additional
            condition that the first letter of the result comes from `x`.

            The usual symbol for this operation is `\prec`.

            .. SEEALSO::

                :meth:`product`, :meth:`succ`, :meth:`over`, :meth:`under`

            EXAMPLES::

                sage: A = algebras.FQSym(QQ).F()
                sage: x = A([2,1])
                sage: A.prec(x, x)
                F[2, 1, 4, 3] + F[2, 4, 1, 3] + F[2, 4, 3, 1]
                sage: y = A([2,1,3])
                sage: A.prec(x, y)
                F[2, 1, 4, 3, 5] + F[2, 4, 1, 3, 5] + F[2, 4, 3, 1, 5]
                 + F[2, 4, 3, 5, 1]
                sage: A.prec(y, x)
                F[2, 1, 3, 5, 4] + F[2, 1, 5, 3, 4] + F[2, 1, 5, 4, 3]
                 + F[2, 5, 1, 3, 4] + F[2, 5, 1, 4, 3] + F[2, 5, 4, 1, 3]
            """
            pre = self.prec_product_on_basis
            return self._module_morphism(self._module_morphism(pre, position=0,
                                                               codomain=self),
                                         position=1)

        def from_symmetric_group_algebra(self, x):
            """
            Return the element of `FQSym` corresponding to the element
            `x` of a symmetric group algebra.

            EXAMPLES::

                sage: A = algebras.FQSym(QQ).F()
                sage: SGA4 = SymmetricGroupAlgebra(QQ, 4)
                sage: x = SGA4([1,3,2,4]) + 5/2 * SGA4([1,2,4,3])
                sage: A.from_symmetric_group_algebra(x)
                5/2*F[1, 2, 4, 3] + F[1, 3, 2, 4]
                sage: A.from_symmetric_group_algebra(SGA4.zero())
                0
            """
            return self._from_dict({Permutation(key): c for (key, c) in x})

    class ElementMethods:
        def to_symmetric_group_algebra(self, n=None):
            """
            Return the element of a symmetric group algebra
            corresponding to the element ``self`` of `FQSym`.

            INPUT:

            - ``n`` -- integer (default: the maximal degree of ``self``);
              the rank of the target symmetric group algebra

            EXAMPLES::

                sage: A = algebras.FQSym(QQ).F()
                sage: x = A([1,3,2,4]) + 5/2 * A([1,2,4,3])
                sage: x.to_symmetric_group_algebra()
                5/2*[1, 2, 4, 3] + [1, 3, 2, 4]
                sage: x.to_symmetric_group_algebra(n=7)
                5/2*[1, 2, 4, 3, 5, 6, 7] + [1, 3, 2, 4, 5, 6, 7]
                sage: a = A.zero().to_symmetric_group_algebra(); a
                0
                sage: parent(a)
                Symmetric group algebra of order 0 over Rational Field

                sage: y = A([1,3,2,4]) + 5/2 * A([2,1])
                sage: y.to_symmetric_group_algebra()
                [1, 3, 2, 4] + 5/2*[2, 1, 3, 4]
                sage: y.to_symmetric_group_algebra(6)
                [1, 3, 2, 4, 5, 6] + 5/2*[2, 1, 3, 4, 5, 6]
            """
            if not self:
                if n is None:
                    n = 0
                return SymmetricGroupAlgebra(self.base_ring(), n).zero()
            m = self.maximal_degree()
            if n is None:
                n = m
            elif n < m:
                raise ValueError("n must be at least the maximal degree")

            SGA = SymmetricGroupAlgebra(self.base_ring(), n)
            return SGA._from_dict({Permutations(n)(key): c for (key, c) in self})

