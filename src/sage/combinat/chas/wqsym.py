# -*- coding: utf-8 -*-
r"""
Word Quasi-symmetric functions

AUTHORS:

- Travis Scrimshaw (2018-04-09)
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
    def __init__(self, alg):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: M = algebras.WQSym(QQ).M()
            sage: TestSuite(M).run()  # long time
        """
        CombinatorialFreeModule.__init__(self, alg.base_ring(),
                                         OrderedSetPartitions(),
                                         category=WQSymBases(alg),
                                         sorting_key=self._basis_key,
                                         bracket="", prefix=self._prefix)

    def _basis_key(self, X):
        """
        Return an object with a total order for sorting the basis
        for the output of an element.

        EXAMPLES::

            sage: M = algebras.WQSym(ZZ).M(); M
            Word Quasi-symmetric functions over Integer Ring in the M basis
            sage: M([[4,1], [3], [2]]) + M([[3,2], [4], [1]]) # indirect doctest
            M[{1, 4}, {3}, {2}] + M[{2, 3}, {4}, {1}]
        """
        return [sorted(part) for part in X]

    def _coerce_map_from_(self, R):
        r"""
        Return ``True`` if there is a coercion from ``R`` into ``self``
        and ``False`` otherwise.

        The things that coerce into ``self`` are

        - word quasi-symmetric functions over a base with
          a coercion map into ``self.base_ring()``

        EXAMPLES::

            sage: M = algebras.WQSym(GF(7)).M(); M
            Word Quasi-symmetric functions over Finite Field of size 7 in the M basis

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
             over Finite Field of size 7 in the M basis to
             Word Quasi-symmetric functions over Integer Ring in the M basis

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

class WordQuasisymmetricFunctions(UniqueRepresentation, Parent):
    r"""
    The word quasi-symmetric functions.

    EXAMPLES::

        sage: WQSym = algebras.WQSym(ZZ)
        sage: WQSym
        Word Quasi-symmetric functions over Integer Ring
        sage: M = WQSym.M()
        sage: M
        Word Quasi-symmetric functions over Integer Ring in the M basis
        sage: M[[]]
        M[]
        sage: M[[1,2,3]]
        M[{1, 2, 3}]
        sage: M[[1,2],[3]]
        M[{1, 2}, {3}]
        sage: M[[1,2,3]] * M[[1,2],[3]]
        M[{1, 2, 3}, {4, 5}, {6}] + M[{1, 2, 3, 4, 5}, {6}]
         + M[{4, 5}, {1, 2, 3}, {6}] + M[{4, 5}, {1, 2, 3, 6}]
         + M[{4, 5}, {6}, {1, 2, 3}]
        sage: x = M[[1],[2],[3]] + 3*M[[2],[1]]
        sage: x.counit()
        0
        sage: x.antipode()
        3*M[{1}, {2}] + 3*M[{1, 2}] - M[{1, 2, 3}] - M[{2, 3}, {1}]
         - M[{3}, {1, 2}] - M[{3}, {2}, {1}]
    """
    def __init__(self, R):
        """
        Initialize ``self``.

        TESTS::

            sage: A = algebras.WQSym(QQ)
            sage: TestSuite(A).run()  # long time
        """
        self._category = HopfAlgebras(R).Graded().Connected()
        Parent.__init__(self, base=R, category=self._category.WithRealizations())

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
            Word Quasi-symmetric functions over Rational Field in the M basis
        """
        return self.M()

    class M(WQSymBasis_abstract):
        r"""
        The `M`-basis of `WQSym`.

        This is the basis `(M_P)`, with `P` ranging over all
        ordered set partitions. See the documentation of
        :class:`WQSym` for details.

        EXAMPLES::

            sage: WQSym = algebras.WQSym(QQ)
            sage: WQSym.M()
            Word Quasi-symmetric functions over Rational Field in the M basis
        """
        _prefix = "M"
        _basis_name = "M"

        def degree_on_basis(self, t):
            """
            Return the degree of an ordered set partition in
            the algebra of word quasi-symmetric functions.

            This is the length of the ordered set partition.

            EXAMPLES::

                sage: A = algebras.WQSym(QQ).M()
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

        def product_on_basis(self, x, y):
            r"""
            Return the `*` associative product of two permutations.

            This is the shifted quasi-shuffle of `x` and `y`.

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
            """
            K = self.basis().keys()
            if not x:
                return self.monomial(y)
            m = max(max(part) for part in x)
            x = [set(part) for part in x]
            yshift = [[val + m for val in part] for part in y]
            def union(X,Y): return X.union(Y)
            return self.sum_of_monomials(ShuffleProduct_overlapping(x, yshift, K, union))

        def coproduct_on_basis(self, x):
            r"""
            Return the coproduct of the basis element indexed by ``x``.

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
            def standardize(P):
                base = sorted(sum((list(part) for part in P), []))
                d = {val: i+1 for i,val in enumerate(base)}
                K = self.indices()
                return K([[d[x] for x in part] for part in P])
            T = self.tensor_square()
            return T.sum_of_monomials((standardize(x[:i]), standardize(x[i:]))
                                      for i in range(len(x) + 1))

class WQSymBases(Category_realization_of_parent):
    r"""
    The category of bases of `WQSym`.
    """
    def __init__(self, base):
        r"""
        Initialize ``self``.

        INPUT:

        - ``base`` -- an instance of `WQSym`

        TESTS::

            sage: from sage.combinat.chas.wqsym import WQSymBases
            sage: WQSym = algebras.WQSym(ZZ)
            sage: bases = WQSymBases(WQSym)
            sage: WQSym.M() in bases
            True
        """
        Category_realization_of_parent.__init__(self, base)

    def _repr_(self):
        r"""
        Return the representation of ``self``.

        EXAMPLES::

            sage: from sage.combinat.chas.wqsym import WQSymBases
            sage: WQSym = algebras.WQSym(ZZ)
            sage: WQSymBases(WQSym)
            Category of bases of Word Quasi-symmetric functions over Integer Ring
        """
        return "Category of bases of {}".format(self.base())

    def super_categories(self):
        r"""
        The super categories of ``self``.

        EXAMPLES::

            sage: from sage.combinat.chas.wqsym import WQSymBases
            sage: WQSym = algebras.WQSym(ZZ)
            sage: bases = WQSymBases(WQSym)
            sage: bases.super_categories()
            [Category of graded connected hopf algebras with basis over Integer Ring,
             Category of realizations of Word Quasi-symmetric functions over Integer Ring]
        """
        return [self.base()._category.WithBasis().Graded(),
                self.base().Realizations()]

    class ParentMethods:
        def _repr_(self):
            """
            Text representation of this basis of `WQSym`.

            EXAMPLES::

                sage: WQSym = algebras.WQSym(ZZ)
                sage: WQSym.M()
                Word Quasi-symmetric functions over Integer Ring in the M basis
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

    class ElementMethods:
        pass

