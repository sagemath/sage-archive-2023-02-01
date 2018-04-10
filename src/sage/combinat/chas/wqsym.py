# -*- coding: utf-8 -*-
r"""
Word Quasi-symmetric functions

AUTHORS:

- Travis Scrimshaw (2018-04-09): initial implementation
"""

# ****************************************************************************
#       Copyright (C) 2018 Travis Scrimshaw <tcscrims at gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.bindable_class import BindableClass
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.hopf_algebras import HopfAlgebras
from sage.categories.realizations import Category_realization_of_parent
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.set_partition_ordered import OrderedSetPartitions
from sage.combinat.shuffle import ShuffleProduct_overlapping

class WQSymBasis_abstract(CombinatorialFreeModule, BindableClass):
    """
    Abstract base class for bases of `WQSym`.

    This must define two attributes:

    - ``_prefix`` -- the basis prefix
    - ``_basis_name`` -- the name of the basis (must match one
      of the names that the basis can be constructed from `WQSym`)
    """
    def __init__(self, alg, graded=True):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: M = algebras.WQSym(QQ).M()
            sage: TestSuite(M).run()  # long time
        """
        CombinatorialFreeModule.__init__(self, alg.base_ring(),
                                         OrderedSetPartitions(),
                                         category=WQSymBases(alg, graded),
                                         bracket="", prefix=self._prefix)

    def _coerce_map_from_(self, R):
        r"""
        Return ``True`` if there is a coercion from ``R`` into ``self``
        and ``False`` otherwise.

        The things that coerce into ``self`` are

        - word quasi-symmetric functions over a base with
          a coercion map into ``self.base_ring()``

        EXAMPLES::

            sage: M = algebras.WQSym(GF(7)).M(); M
            Word Quasi-symmetric functions over Finite Field of size 7 in the Monomial basis

        Elements of the word quasi-symmetric functions canonically coerce in::

            sage: x, y = M([[1]]), M([[2,1]])
            sage: M.coerce(x+y) == x+y
            True

        The word quasi-symmetric functions over `\ZZ` coerces in,
        since `\ZZ` coerces to `\GF{7}`::

            sage: N = algebras.WQSym(ZZ).M()
            sage: Nx, Ny = N([[1]]), N([[2,1]])
            sage: z = M.coerce(Nx+Ny); z
            M[{1}] + M[{1, 2}]
            sage: z.parent() is M
            True

        However, `\GF{7}` does not coerce to `\ZZ`, so word
        quasi-symmetric functions over `\GF{7}` does not coerce
        to the same algebra over `\ZZ`::

            sage: N.coerce(y)
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Word Quasi-symmetric functions
             over Finite Field of size 7 in the Monomial basis to
             Word Quasi-symmetric functions over Integer Ring in the Monomial basis

        TESTS::

            sage: M = algebras.WQSym(ZZ).M()
            sage: N = algebras.WQSym(QQ).M()
            sage: M.has_coerce_map_from(N)
            False
            sage: N.has_coerce_map_from(M)
            True
            sage: M.has_coerce_map_from(QQ)
            False
            sage: N.has_coerce_map_from(QQ)
            True
            sage: M.has_coerce_map_from(PolynomialRing(ZZ, 3, 'x,y,z'))
            False
        """
        # word quasi-symmetric functions in the same variables
        # over any base that coerces in:
        if isinstance(R, WQSymBasis_abstract):
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
        return super(WQSymBasis_abstract, self)._coerce_map_from_(R)

    @cached_method
    def an_element(self):
        """
        Return an element of ``self``.

        EXAMPLES::

            sage: M = algebras.WQSym(QQ).M()
            sage: M.an_element()
            M[{1}] + 2*M[{1}, {2}]
        """
        return self([[1]]) + 2*self([[1],[2]])

    def some_elements(self):
        """
        Return some elements of the word quasi-symmetric functions.

        EXAMPLES::

            sage: M = algebras.WQSym(QQ).M()
            sage: M.some_elements()
            [M[], M[{1}], M[{1, 2}],
             M[{1}] + M[{1}, {2}],
             M[] + 1/2*M[{1}]]
        """
        u = self.one()
        o = self([[1]])
        s = self.base_ring().an_element()
        return [u, o, self([[1,2]]), o + self([[1],[2]]), u + s*o]

class WordQuasisymmetricFunctions(UniqueRepresentation, Parent):
    r"""
    The word quasi-symmetric functions.

    The ring of word quasi-symmetric functions can be defined as a
    subring of the ring of all bounded-degree noncommutative power
    series in countably many indeterminates (i.e., elements in
    `R \langle \langle x_1, x_2, x_3, \ldots \rangle \rangle` of bounded
    degree). Namely, consider words over the alphabet `\{1, 2, 3, \ldots\}`;
    every noncommutative power series is an infinite `R`-linear
    combination of these words.
    For each such word `w`, we define the *packing* of `w` to be the
    word `\operatorname{pack}(w)` that is obtained from `w` by replacing
    the smallest letter that appears in `w` by `1`, the second-smallest
    letter that appears in `w` by `2`, etc. (for example,
    `\operatorname{pack}(4112774) = 3112443`).
    A word `w` is said to be *packed* if `\operatorname{pack}(w) = w`.
    For each packed word `u`, we define the noncommutative power series
    `\mathbf{M}_u = \sum w`, where the sum ranges over all words `w`
    satisfying `\operatorname{pack}(w) = u`.
    The span of these power series `\mathbf{M}_u` is a subring of the
    ring of all noncommutative power series; it is called the ring of
    word quasi-symmetric functions, and is denoted by `WQSym`.

    For each nonnegative integer `n`, there is a bijection between
    packed words of length `n` and ordered set partitions of
    `\{1, 2, \ldots, n\}`. Under this bijection, a packed word
    `u = (u_1, u_2, \ldots, u_n)` of length `n` corresponds to the
    ordered set partition `P = (P_1, P_2, \ldots, P_k)` of
    `\{1, 2, \ldots, n\}` whose `i`-th part `P_i` (for each `i`) is the
    set of all `j \in \{1, 2, \ldots, n\}` such that `u_j = i`.

    The basis element `\mathbf{M}_u` is also denoted as `\mathbf{M}_P`
    in this situation and is implemented using the latter indexing.
    The basis `(\mathbf{M}_P)_P` is called the *Monomial basis* and
    is implemented at
    :class:`~sage.combinat.chas.wqsym.WordQuasisymmetricFunctions.M`.

    `WQSym` is endowed with a connected graded Hopf algebra structure (see
    Section 2.2 of [NoThWi08]_, Section 1.1 of [FoiMal14]_ and
    Section 4.3.2 of [MeNoTh11]_) given by

    .. MATH::

        \Delta(\mathbf{M}_{(P_1,\ldots,P_{\ell})}) = \sum_{i=0}^{\ell}
            \mathbf{M}_{\operatorname{st}(P_1, \ldots, P_i)} \otimes
            \mathbf{M}_{\operatorname{st}(P_{i+1}, \ldots, P_{\ell})}.

    Here, for any ordered set partition `(Q_1, \ldots, Q_k)` of a
    finite set `Z` of integers, we let `\operatorname{st}(Q_1, \ldots, Q_k)`
    denote the set partition obtained from `Z` by replacing the smallest
    element appearing in it by `1`, the second-smallest element by `2`,
    and so on.

    A rule for multiplying elements of the monomial basis relies on the
    *quasi-shuffle product* of two ordered set partitions.
    The quasi-shuffle product `\Box` is given by
    :class:`~sage.combinat.shuffle.ShuffleProduct_overlapping` with ``+``
    being the union of the sets. The product `\mathbf{M}_P \mathbf{M}_Q`
    for two ordered set partitions `P` and `Q` of `[n]` and `[m]`
    is then given by

    .. MATH::

        \mathbf{M}_P \mathbf{M}_Q
        = \sum_{R \in P \Box Q^+} \mathbf{M}_R ,

    where `Q^+` means `Q` with all numbers shifted upwards by `n`.

    Sometimes, `WQSym` is also denoted as `NCQSym`.

    REFERENCES:

    - [FoiMal14]_
    - [MeNoTh11]_
    - [NoThWi08]_

    EXAMPLES::

        sage: WQSym = algebras.WQSym(ZZ)
        sage: WQSym
        Word Quasi-symmetric functions over Integer Ring
        sage: M = WQSym.M()
        sage: M
        Word Quasi-symmetric functions over Integer Ring in the Monomial basis
        sage: M[[]]
        M[]
        sage: M[[1,2,3]]
        M[{1, 2, 3}]
        sage: M[[1,2],[3]]
        M[{1, 2}, {3}]
        sage: M[[2, 3], [5], [6], [4], [1]].coproduct()
        M[] # M[{2, 3}, {5}, {6}, {4}, {1}] + M[{1, 2}] # M[{3}, {4}, {2}, {1}]
         + M[{1, 2}, {3}] # M[{3}, {2}, {1}] + M[{1, 2}, {3}, {4}] # M[{2}, {1}]
         + M[{1, 2}, {4}, {5}, {3}] # M[{1}] + M[{2, 3}, {5}, {6}, {4}, {1}] # M[]
        sage: M[[1,2,3]] * M[[1,2],[3]]
        M[{1, 2, 3}, {4, 5}, {6}] + M[{1, 2, 3, 4, 5}, {6}]
         + M[{4, 5}, {1, 2, 3}, {6}] + M[{4, 5}, {1, 2, 3, 6}]
         + M[{4, 5}, {6}, {1, 2, 3}]
        sage: M[[1,2,3]].antipode()
        -M[{1, 2, 3}]
        sage: M[[1], [2], [3]].antipode()
        -M[{1, 2, 3}] - M[{2, 3}, {1}] - M[{3}, {1, 2}] - M[{3}, {2}, {1}]
        sage: x = M[[1],[2],[3]] + 3*M[[2],[1]]
        sage: x.counit()
        0
        sage: x.antipode()
        3*M[{1}, {2}] + 3*M[{1, 2}] - M[{1, 2, 3}] - M[{2, 3}, {1}]
         - M[{3}, {1, 2}] - M[{3}, {2}, {1}]

    TESTS::

        sage: a = M[OrderedSetPartition([[1]])]
        sage: b = M[OrderedSetPartitions(1)([[1]])]
        sage: c = M[[1]]
        sage: a == b == c
        True

    .. TODO::

        Dendriform structure.
        Bergeron-Zabrocki/Menous-Novelli-Thibon basis.
    """
    def __init__(self, R):
        """
        Initialize ``self``.

        TESTS::

            sage: A = algebras.WQSym(QQ)
            sage: TestSuite(A).run()  # long time
        """
        category = HopfAlgebras(R).Graded().Connected()
        Parent.__init__(self, base=R, category=category.WithRealizations())

    def _repr_(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: algebras.WQSym(QQ)  # indirect doctest
            Word Quasi-symmetric functions over Rational Field
        """
        return "Word Quasi-symmetric functions over {}".format(self.base_ring())

    def a_realization(self):
        r"""
        Return a particular realization of ``self`` (the `M`-basis).

        EXAMPLES::

            sage: WQSym = algebras.WQSym(QQ)
            sage: WQSym.a_realization()
            Word Quasi-symmetric functions over Rational Field in the Monomial basis
        """
        return self.M()

    _shorthands = tuple(['M', 'X', 'C'])

    class M(WQSymBasis_abstract):
        r"""
        The Monomial basis of `WQSym`.

        The family `(\mathbf{M}_u)`, as defined in
        :class:`~sage.combinat.chas.wqsym.WordQuasiSymmetricFunctions`
        with `u` ranging over all packed words, is a basis for the
        free `R`-module `WQSym` and called the *Monomial basis*.
        Here it is labelled using ordered set partitions.

        EXAMPLES::

            sage: WQSym = algebras.WQSym(QQ)
            sage: WQSym.M()
            Word Quasi-symmetric functions over Rational Field in the Monomial basis
        """
        _prefix = "M"
        _basis_name = "Monomial"

        def product_on_basis(self, x, y):
            r"""
            Return the (associative) `*` product of the basis elements
            of ``self`` indexed by the ordered set partitions `x` and
            `y`.

            This is the shifted quasi-shuffle product of `x` and `y`.

            EXAMPLES::

                sage: A = algebras.WQSym(QQ).M()
                sage: x = OrderedSetPartition([[1],[2,3]])
                sage: y = OrderedSetPartition([[1,2]])
                sage: z = OrderedSetPartition([[1,2],[3]])
                sage: A.product_on_basis(x, y)
                M[{1}, {2, 3}, {4, 5}] + M[{1}, {2, 3, 4, 5}]
                 + M[{1}, {4, 5}, {2, 3}] + M[{1, 4, 5}, {2, 3}]
                 + M[{4, 5}, {1}, {2, 3}]
                sage: A.product_on_basis(x, z)
                M[{1}, {2, 3}, {4, 5}, {6}] + M[{1}, {2, 3, 4, 5}, {6}]
                 + M[{1}, {4, 5}, {2, 3}, {6}] + M[{1}, {4, 5}, {2, 3, 6}]
                 + M[{1}, {4, 5}, {6}, {2, 3}] + M[{1, 4, 5}, {2, 3}, {6}]
                 + M[{1, 4, 5}, {2, 3, 6}] + M[{1, 4, 5}, {6}, {2, 3}]
                 + M[{4, 5}, {1}, {2, 3}, {6}] + M[{4, 5}, {1}, {2, 3, 6}]
                 + M[{4, 5}, {1}, {6}, {2, 3}] + M[{4, 5}, {1, 6}, {2, 3}]
                 + M[{4, 5}, {6}, {1}, {2, 3}]
                sage: A.product_on_basis(y, y)
                M[{1, 2}, {3, 4}] + M[{1, 2, 3, 4}] + M[{3, 4}, {1, 2}]

            TESTS::

                sage: one = OrderedSetPartition([])
                sage: all(A.product_on_basis(one, z) == A(z) == A.basis()[z] for z in OrderedSetPartitions(3))
                True
                sage: all(A.product_on_basis(z, one) == A(z) == A.basis()[z] for z in OrderedSetPartitions(3))
                True
            """
            K = self.basis().keys()
            if not x:
                return self.monomial(y)
            m = max(max(part) for part in x) # The degree of x
            x = [set(part) for part in x]
            yshift = [[val + m for val in part] for part in y]
            def union(X,Y): return X.union(Y)
            return self.sum_of_monomials(ShuffleProduct_overlapping(x, yshift, K, union))

        def coproduct_on_basis(self, x):
            r"""
            Return the coproduct of ``self`` on the basis element
            indexed by the ordered set partition ``x``.

            EXAMPLES::

                sage: M = algebras.WQSym(QQ).M()

                sage: M.coproduct(M.one())  # indirect doctest
                M[] # M[]
                sage: M.coproduct( M([[1]]) )  # indirect doctest
                M[] # M[{1}] + M[{1}] # M[]
                sage: M.coproduct( M([[1,2]]) )
                M[] # M[{1, 2}] + M[{1, 2}] # M[]
                sage: M.coproduct( M([[1], [2]]) )
                M[] # M[{1}, {2}] + M[{1}] # M[{1}] + M[{1}, {2}] # M[]
            """
            if not len(x):
                return self.one().tensor(self.one())
            K = self.indices()
            def standardize(P): # standardize an ordered set partition
                base = sorted(sum((list(part) for part in P), []))
                # base is the ground set of P, as a sorted list.
                d = {val: i+1 for i,val in enumerate(base)}
                # d is the unique order isomorphism from base to
                # {1, 2, ..., |base|} (encoded as dict).
                return K([[d[x] for x in part] for part in P])
            T = self.tensor_square()
            return T.sum_of_monomials((standardize(x[:i]), standardize(x[i:]))
                                      for i in range(len(x) + 1))

    Monomial = M

    class X(WQSymBasis_abstract):
        r"""
        The Characteristic basis of `WQSym`.

        The *Characteristic basis* is a graded basis `(X_P)` of `WQSym`,
        indexed by ordered set partitions `P`. It is defined by

        .. MATH::

            X_P = (-1)^{\ell(P)} \mathbf{M}_P ,

        where `(\mathbf{M}_P)_P` denotes the Monomial basis,
        and where `\ell(P)` denotes the number of blocks in an ordered
        set partition `P`.

        EXAMPLES::

            sage: WQSym = algebras.WQSym(QQ)
            sage: X = WQSym.X(); X
            Word Quasi-symmetric functions over Rational Field in the Characteristic basis

            sage: X[[1,2,3]] * X[[1,2],[3]]
            X[{1, 2, 3}, {4, 5}, {6}] - X[{1, 2, 3, 4, 5}, {6}]
             + X[{4, 5}, {1, 2, 3}, {6}] - X[{4, 5}, {1, 2, 3, 6}]
             + X[{4, 5}, {6}, {1, 2, 3}]

            sage: X[[1, 4], [3], [2]].coproduct()
            X[] # X[{1, 4}, {3}, {2}] + X[{1, 2}] # X[{2}, {1}]
             + X[{1, 3}, {2}] # X[{1}] + X[{1, 4}, {3}, {2}] # X[]

            sage: M = WQSym.M()
            sage: M(X[[1, 2, 3]])
            -M[{1, 2, 3}]
            sage: M(X[[1, 3], [2]])
            M[{1, 3}, {2}]
            sage: X(M[[1, 2, 3]])
            -X[{1, 2, 3}]
            sage: X(M[[1, 3], [2]])
            X[{1, 3}, {2}]
        """
        _prefix = "X"
        _basis_name = "Characteristic"

        def __init__(self, alg):
            """
            Initialize ``self``.

            EXAMPLES::

                sage: X = algebras.WQSym(QQ).X()
                sage: TestSuite(X).run()  # long time
            """
            WQSymBasis_abstract.__init__(self, alg)

            M = self.realization_of().M()
            mone = -self.base_ring().one()
            def sgn(P): return mone**len(P)
            self.module_morphism(codomain=M, diagonal=sgn).register_as_coercion()
            M.module_morphism(codomain=self, diagonal=sgn).register_as_coercion()

    Characteristic = X

    class C(WQSymBasis_abstract):
        r"""
        The Cone basis of `WQSym`.

        Let `(X_P)_P` denote the Characteristic basis of `WQSym`.
        Denote the quasi-shuffle of two ordered set partitions `A` and
        `B` by `A \Box B`. For an ordered set partition
        `P = (P_1, \ldots, P_{\ell})`, we form a list of ordered set
        partitions `[P] := (P'_1, \ldots, P'_k)` as follows.
        Define a strictly decreasing sequence of integers
        `\ell + 1 = i_0 > i_1 > \cdots > i_k = 1` recursively by
        requiring that `\min P_{i_j} \leq \min P_a` for all `a < i_{j-1}`.
        Set `P'_j = (P_{i_j}, \ldots, P_{i_{j-1}-1})`.

        The *Cone basis* `(C_P)_P` is defined by

        .. MATH::

            C_P = \sum_Q X_Q,

        where the sum is over all elements `Q` of the quasi-shuffle
        product `P'_1 \Box P'_2 \Box \cdots \Box P'_k` with
        `[P] = (P'_1, \ldots, P'_k)`.

        EXAMPLES::

            sage: WQSym = algebras.WQSym(QQ)
            sage: C = WQSym.C()
            sage: C
            Word Quasi-symmetric functions over Rational Field in the Cone basis

            sage: X = WQSym.X()
            sage: X(C[[2,3],[1,4]])
            X[{1, 2, 3, 4}] + X[{1, 4}, {2, 3}] + X[{2, 3}, {1, 4}]
            sage: X(C[[1,4],[2,3]])
            X[{1, 4}, {2, 3}]
            sage: X(C[[2,3],[1],[4]])
            X[{1}, {2, 3}, {4}] + X[{1}, {2, 3, 4}] + X[{1}, {4}, {2, 3}]
             + X[{1, 2, 3}, {4}] + X[{2, 3}, {1}, {4}]
            sage: X(C[[3], [2, 5], [1, 4]])
            X[{1, 2, 3, 4, 5}] + X[{1, 2, 4, 5}, {3}] + X[{1, 3, 4}, {2, 5}]
             + X[{1, 4}, {2, 3, 5}] + X[{1, 4}, {2, 5}, {3}]
             + X[{1, 4}, {3}, {2, 5}] + X[{2, 3, 5}, {1, 4}]
             + X[{2, 5}, {1, 3, 4}] + X[{2, 5}, {1, 4}, {3}]
             + X[{2, 5}, {3}, {1, 4}] + X[{3}, {1, 2, 4, 5}]
             + X[{3}, {1, 4}, {2, 5}] + X[{3}, {2, 5}, {1, 4}]
            sage: C(X[[2,3],[1,4]])
            -C[{1, 2, 3, 4}] - C[{1, 4}, {2, 3}] + C[{2, 3}, {1, 4}]

        REFERENCES:

        - Section 4 of [Early2017]_
        """
        _prefix = "C"
        _basis_name = "Cone"

        def __init__(self, alg):
            """
            Initialize ``self``.

            EXAMPLES::

                sage: C = algebras.WQSym(QQ).C()
                sage: TestSuite(C).run()  # long time
            """
            WQSymBasis_abstract.__init__(self, alg)

            X = self.realization_of().X()
            phi = self.module_morphism(self._C_to_X, codomain=X, unitriangular="upper")
            phi.register_as_coercion()
            (~phi).register_as_coercion()

        def some_elements(self):
            """
            Return some elements of the word quasi-symmetric functions
            in the Cone basis.

            EXAMPLES::

                sage: C = algebras.WQSym(QQ).C()
                sage: C.some_elements()
                [C[], C[{1}], C[{1, 2}], C[] + 1/2*C[{1}]]
            """
            u = self.one()
            o = self([[1]])
            s = self.base_ring().an_element()
            return [u, o, self([[1,2]]), u + s*o]

        def _C_to_X(self, P):
            """
            Return the image of the basis element of ``self`` indexed
            by ``P`` in the Characteristic basis.

            EXAMPLES::

                sage: C = algebras.WQSym(QQ).C()
                sage: OSP = C.basis().keys()
                sage: C._C_to_X(OSP([[2,3],[1,4]]))
                X[{1, 2, 3, 4}] + X[{1, 4}, {2, 3}] + X[{2, 3}, {1, 4}]
            """
            X = self.realization_of().X()
            if not P:
                return X.one()

            OSP = self.basis().keys()

            # Convert to standard set of ordered set partitions
            temp = list(P)
            data = []
            while temp:
                i = min(min(X) for X in temp)
                for j,A in enumerate(temp):
                    if i in A:
                        data.append(OSP(temp[j:]))
                        temp = temp[:j]
                        break

            # Perform the quasi-shuffle product
            cur = {data[0]: 1}
            for B in data[1:]:
                ret = {}
                for A in cur:
                    for C in ShuffleProduct_overlapping(A, B, element_constructor=OSP):
                        if C in ret:
                            ret[C] += cur[A]
                        else:
                            ret[C] = cur[A]
                cur = ret

            # Return the result in the X basis
            return X._from_dict(cur, coerce=True)

    Cone = C

class WQSymBases(Category_realization_of_parent):
    r"""
    The category of bases of `WQSym`.
    """
    def __init__(self, base, graded):
        r"""
        Initialize ``self``.

        INPUT:

        - ``base`` -- an instance of `WQSym`
        - ``graded`` -- boolean; if the basis is graded or filtered

        TESTS::

            sage: from sage.combinat.chas.wqsym import WQSymBases
            sage: WQSym = algebras.WQSym(ZZ)
            sage: bases = WQSymBases(WQSym, True)
            sage: WQSym.M() in bases
            True
        """
        self._graded = graded
        Category_realization_of_parent.__init__(self, base)

    def _repr_(self):
        r"""
        Return the representation of ``self``.

        EXAMPLES::

            sage: from sage.combinat.chas.wqsym import WQSymBases
            sage: WQSym = algebras.WQSym(ZZ)
            sage: WQSymBases(WQSym, True)
            Category of graded bases of Word Quasi-symmetric functions over Integer Ring
            sage: WQSymBases(WQSym, False)
            Category of filtered bases of Word Quasi-symmetric functions over Integer Ring
        """
        if self._graded:
            type_str = "graded"
        else:
            type_str = "filtered"
        return "Category of {} bases of {}".format(type_str, self.base())

    def super_categories(self):
        r"""
        The super categories of ``self``.

        EXAMPLES::

            sage: from sage.combinat.chas.wqsym import WQSymBases
            sage: WQSym = algebras.WQSym(ZZ)
            sage: bases = WQSymBases(WQSym, True)
            sage: bases.super_categories()
            [Category of realizations of Word Quasi-symmetric functions over Integer Ring,
             Join of Category of realizations of hopf algebras over Integer Ring
                 and Category of graded algebras over Integer Ring,
             Category of graded connected hopf algebras with basis over Integer Ring]

            sage: bases = WQSymBases(WQSym, False)
            sage: bases.super_categories()
            [Category of realizations of Word Quasi-symmetric functions over Integer Ring,
             Join of Category of realizations of hopf algebras over Integer Ring
                 and Category of graded algebras over Integer Ring,
             Join of Category of filtered connected hopf algebras with basis over Integer Ring
                 and Category of graded algebras over Integer Ring]
        """
        R = self.base().base_ring()
        cat = HopfAlgebras(R).Graded().WithBasis()
        if self._graded:
            cat = cat.Graded()
        else:
            cat = cat.Filtered()
        return [self.base().Realizations(),
                HopfAlgebras(R).Graded().Realizations(),
                cat.Connected()]

    class ParentMethods:
        def _repr_(self):
            """
            Text representation of this basis of `WQSym`.

            EXAMPLES::

                sage: WQSym = algebras.WQSym(ZZ)
                sage: WQSym.M()
                Word Quasi-symmetric functions over Integer Ring in the Monomial basis
            """
            return "{} in the {} basis".format(self.realization_of(), self._basis_name)

        def __getitem__(self, p):
            """
            Return the basis element indexed by ``p``.

            INPUT:

            - ``p`` -- an ordered set partition

            EXAMPLES::

                sage: M = algebras.WQSym(QQ).M()
                sage: M[[1, 3, 2]]
                M[{1, 2, 3}]
                sage: M[[1,3],[2]]
                M[{1, 3}, {2}]
                sage: M[OrderedSetPartition([[2],[1,4],[3,5]])]
                M[{2}, {1, 4}, {3, 5}]
            """
            try:
                return self.monomial(self._indices(p))
            except TypeError:
                return self.monomial(self._indices([p]))

        def is_field(self, proof=True):
            """
            Return whether ``self`` is a field.

            EXAMPLES::

                sage: M = algebras.WQSym(QQ).M()
                sage: M.is_field()
                False
            """
            return False

        def is_commutative(self):
            """
            Return whether ``self`` is commutative.

            EXAMPLES::

                sage: M = algebras.WQSym(ZZ).M()
                sage: M.is_commutative()
                False
            """
            return self.base_ring().is_zero()

        def one_basis(self):
            """
            Return the index of the unit.

            EXAMPLES::

                sage: A = algebras.WQSym(QQ).M()
                sage: A.one_basis()
                []
            """
            OSP = self.basis().keys()
            return OSP([])

        def degree_on_basis(self, t):
            """
            Return the degree of an ordered set partition in
            the algebra of word quasi-symmetric functions.

            This is the sum of the sizes of the blocks of the
            ordered set partition.

            EXAMPLES::

                sage: A = algebras.WQSym(QQ).M()
                sage: u = OrderedSetPartition([[2,1]])
                sage: A.degree_on_basis(u)
                2
                sage: u = OrderedSetPartition([[2], [1]])
                sage: A.degree_on_basis(u)
                2
            """
            return sum(len(part) for part in t)

    class ElementMethods:
        pass

