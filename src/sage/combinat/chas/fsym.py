r"""
Poirier-Reutenauer Hopf algebra of standard tableau

AUTHORS:

- Franco Saliola (2012): initial implementation
- Travis Scrimshaw (2018-04-11): added missing doctests and reorganization
"""

# ****************************************************************************
#       Copyright (C) 2012 Franco Saliola   <saliola at gmail.com>
#                     2018 Travis Scrimshaw <tcscrims at gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.bindable_class import BindableClass
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.misc_c import prod
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.functions.other import factorial

from sage.categories.rings import Rings
from sage.categories.tensor import tensor
from sage.categories.realizations import Category_realization_of_parent
from sage.categories.graded_hopf_algebras_with_basis import GradedHopfAlgebrasWithBasis
from sage.categories.graded_hopf_algebras import GradedHopfAlgebras

from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.tableau import Tableau, StandardTableau, StandardTableaux
from sage.combinat.set_partition_ordered import OrderedSetPartitions
from sage.combinat.permutation import Permutations
from sage.combinat.sf.sf import SymmetricFunctions

class FSymBasis_abstract(CombinatorialFreeModule, BindableClass):
    """
    Abstract base class for bases of `FSym`.

    This must define two attributes:

    - ``_prefix`` -- the basis prefix
    - ``_basis_name`` -- the name of the basis (must match one
      of the names that the basis can be constructed from `FSym`)
    """
    def __init__(self, alg, graded=True):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: G = algebras.FSym(QQ).G()
            sage: TestSuite(G).run()  # long time

        Checks for the antipode::

            sage: FSym = algebras.FSym(QQ)
            sage: G = FSym.G()
            sage: for b in G.basis(degree=3):
            ....:     print("%s : %s" % (b, b.antipode()))
            G[123] : -G[1|2|3]
            G[13|2] : -G[13|2]
            G[12|3] : -G[12|3]
            G[1|2|3] : -G[123]

            sage: F = FSym.dual().F()
            sage: for b in F.basis(degree=3):
            ....:     print("%s : %s" % (b, b.antipode()))
            F[123] : -F[1|2|3]
            F[13|2] : -F[13|2]
            F[12|3] : -F[12|3]
            F[1|2|3] : -F[123]
        """
        CombinatorialFreeModule.__init__(self, alg.base_ring(),
                                         StandardTableaux(),
                                         category=FSymBases(alg),
                                         bracket="", prefix=self._prefix)

    def _coerce_map_from_(self, R):
        r"""
        Return ``True`` if there is a coercion from ``R`` into ``self``
        and ``False`` otherwise.

        The things that coerce into ``self`` are

        - word quasi-symmetric functions over a base with
          a coercion map into ``self.base_ring()``

        EXAMPLES::

            sage: G = algebras.FSym(GF(7)).G(); G
            Hopf algebra of standard tableaux over the Finite Field of size 7
             in the Fundamental basis

        Elements of the word quasi-symmetric functions canonically coerce in::

            sage: x, y = G([[1]]), G([[1,3],[2]])
            sage: G.coerce(x + y) == x + y
            True

        The word quasi-symmetric functions over `\ZZ` coerces in,
        since `\ZZ` coerces to `\GF{7}`::

            sage: H = algebras.FSym(ZZ).G()
            sage: Hx, Hy = H([[1]]), H([[1,3],[2]])
            sage: z = G.coerce(Hx+Hy); z
            G[1] + G[13|2]
            sage: z.parent() is G
            True

        However, `\GF{7}` does not coerce to `\ZZ`, so word
        quasi-symmetric functions over `\GF{7}` does not coerce
        to the same algebra over `\ZZ`::

            sage: H.coerce(y)
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Hopf algebra of standard tableaux
             over the Finite Field of size 7 in the Fundamental basis
             to Hopf algebra of standard tableaux over the Integer Ring in the Fundamental basis

        TESTS::

            sage: G = algebras.FSym(ZZ).G()
            sage: H = algebras.FSym(QQ).G()
            sage: G.has_coerce_map_from(H)
            False
            sage: H.has_coerce_map_from(G)
            True
            sage: G.has_coerce_map_from(QQ)
            False
            sage: H.has_coerce_map_from(QQ)
            True
            sage: G.has_coerce_map_from(PolynomialRing(ZZ, 3, 'x,y,z'))
            False
        """
        # word quasi-symmetric functions in the same variables
        # over any base that coerces in:
        if isinstance(R, FSymBasis_abstract):
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
        return super(FSymBasis_abstract, self)._coerce_map_from_(R)

    def some_elements(self):
        """
        Return some elements of the word quasi-symmetric functions.

        EXAMPLES::

            sage: G = algebras.FSym(QQ).G()
            sage: G.some_elements()
            [G[], G[1], G[12], G[1] + G[1|2], G[] + 1/2*G[1]]
        """
        u = self.one()
        o = self([[1]])
        s = self.base_ring().an_element()
        return [u, o, self([[1,2]]), o + self([[1],[2]]), u + s*o]

    def _repr_term(self, phi):
        r"""
        The string representation of a basis element.

        EXAMPLES:

        We use a compact notation for set partitions::

            sage: FSym = algebras.FSym(QQ)
            sage: G = FSym.G()
            sage: G.zero()
            0
            sage: G[[1],[2],[3],[4]]
            G[1|2|3|4]
            sage: G[[1,3,5],[2,4]]
            G[135|24]
        """
        return "{}[{}]".format(self._prefix,
                               "|".join("".join(map(str,block))
                                                for block in sorted(phi, key=min)))

class FSymBases(Category_realization_of_parent):
    r"""
    The category of bases of `FSym`.
    """
    def super_categories(self):
        """
        The super categories of ``self``.

        EXAMPLES::

            sage: from sage.combinat.chas.fsym import FSymBases
            sage: FSym = algebras.FSym(ZZ)
            sage: bases = FSymBases(FSym)
            sage: bases.super_categories()
            [Category of realizations of Hopf algebra of standard tableaux over the Integer Ring,
             Join of Category of realizations of hopf algebras over Integer Ring
                 and Category of graded algebras over Integer Ring,
             Category of graded connected hopf algebras with basis over Integer Ring]
        """
        R = self.base().base_ring()
        return [self.base().Realizations(),
                GradedHopfAlgebras(R).Realizations(),
                GradedHopfAlgebrasWithBasis(R).Connected()]

    class ParentMethods:
        def _repr_(self):
            """
            Text representation of this basis of `FSym`.

            EXAMPLES::

                sage: FSym = algebras.FSym(ZZ)
                sage: FSym.G()
                Hopf algebra of standard tableaux over the Integer Ring
                 in the Fundamental basis
            """
            return "{} in the {} basis".format(self.realization_of(), self._basis_name)

        def __getitem__(self, key):
            r"""
            Override the ``__getitem__`` method to allow passing of
            arguments to SetPartition.

            EXAMPLES:

            Construct the basis element indexed by a set partition by
            passing data that defines the set partition::

                sage: FSym = algebras.FSym(QQ)
                sage: G = FSym.G()
                sage: G[[1,3],[2]]
                G[13|2]
                sage: G[[1,3],[2]].leading_support() in StandardTableaux()
                True
            """
            try:
                return self.monomial(self._indices(list(key)))
            except (TypeError, ValueError):
                return self.monomial(self._indices([key]))

        def basis(self, degree=None):
            r"""
            The basis elements (optionally: of the specified degree).

            OUTPUT: Family

            EXAMPLES::

                sage: FSym = algebras.FSym(QQ)
                sage: TG = FSym.G()
                sage: TG.basis()
                Lazy family (Term map from Standard tableaux to Hopf algebra of standard tableaux over the Rational Field in the Fundamental basis(i))_{i in Standard tableaux}
                sage: TG.basis().keys()
                Standard tableaux
                sage: TG.basis(degree=3).keys()
                Standard tableaux of size 3
                sage: TG.basis(degree=3).list()
                [G[123], G[13|2], G[12|3], G[1|2|3]]
            """
            from sage.combinat.family import Family
            if degree is None:
                return Family(self._indices, self.monomial)
            else:
                return Family(StandardTableaux(degree), self.monomial)

        @cached_method
        def one_basis(self):
            r"""
            Return the basis index corresponding to `1`.

            EXAMPLES::

                sage: FSym = algebras.FSym(QQ)
                sage: TG = FSym.G()
                sage: TG.one_basis()
                []
            """
            return self._indices([])

        def duality_pairing(self, x, y):
            r"""
            Scalar product for which the Hopf algebra of tableaux is
            self-dual: the basis of tableaux is declared to be orthonormal.

            EXAMPLES::

                sage: FSym = algebras.FSym(QQ)
                sage: G = FSym.G()
                sage: F = G.dual_basis()
                sage: t1 = StandardTableau([[1,3,5],[2,4]])
                sage: t2 = StandardTableau([[1,3],[2,5],[4]])
                sage: G.duality_pairing(G[t1], F[t2])
                0
                sage: G.duality_pairing(G[t1], F[t1])
                1
                sage: G.duality_pairing(G[t2], F[t2])
                1

                sage: z = G[[1,3,5],[2,4]]
                sage: all(F.duality_pairing(F[p1] * F[p2], z)
                ....:     for ((p1, p2), c) in z.coproduct())
                True
            """
            x = self(x)
            y = self.dual_basis()(y)
            return sum(coeff * y[t] for (t, coeff) in x)

        def duality_pairing_matrix(self, basis, degree):
            r"""
            The matrix of scalar products between elements of `FSym` and
            elements of `FSym^*`.

            INPUT:

            - ``basis`` -- a basis of the dual Hopf algebra
            - ``degree`` -- a non-negative integer

            OUTPUT:

            - the matrix of scalar products between the basis ``self`` and the
              basis ``basis`` in the dual Hopf algebra of degree ``degree``

            EXAMPLES::

                sage: FSym = algebras.FSym(QQ)
                sage: G = FSym.G()
                sage: G.duality_pairing_matrix(G.dual_basis(), 3)
                [1 0 0 0]
                [0 1 0 0]
                [0 0 1 0]
                [0 0 0 1]
            """
            from sage.matrix.constructor import matrix
            keys = self.basis(degree=degree).keys()
            return matrix(self.base_ring(),
                          [[self.duality_pairing(self[t], basis[s])
                            for t in keys] for s in keys])

        def degree_on_basis(self, t):
            """
            Return the degree of a standard tableaux in the algebra
            of free symmetric functions.

            This is the size of the tableaux.

            EXAMPLES::

                sage: G = algebras.FSym(QQ).G()
                sage: t = StandardTableau([[1,3],[2]])
                sage: G.degree_on_basis(t)
                3
                sage: u = StandardTableau([[1,3,4,5],[2]])
                sage: G.degree_on_basis(u)
                5
            """
            return t.size()

    class ElementMethods:
        def duality_pairing(self, other):
            r"""
            Compute the pairing between ``self`` and an element ``other`` of the dual.

            EXAMPLES::

                sage: FSym = algebras.FSym(QQ)
                sage: G = FSym.G()
                sage: F = G.dual_basis()
                sage: elt = G[[1,3],[2]] - 3*G[[1,2],[3]]
                sage: elt.duality_pairing(F[[1,3],[2]])
                1
                sage: elt.duality_pairing(F[[1,2],[3]])
                -3
                sage: elt.duality_pairing(F[[1,2]])
                0
            """
            return self.parent().duality_pairing(self, other)

class FreeSymmetricFunctions(UniqueRepresentation, Parent):
    """
    The free symmetric functions.

    The *free symmetric functions* is the combinatorial Hopf algebra
    defined using tableaux and denoted `FSym`.
    """
    def __init__(self, base_ring):
        r"""
        TESTS:

        This Hopf algebra embeds as a Hopf-subalgebra of the Hopf algebra of
        permutations: the basis element indexed by a tableau `t` is sent to
        the sum of permutations whose `P`-tableau is `t`.

        ::

            sage: TG = algebras.FSym(QQ).G()
            sage: F = algebras.FQSym(QQ).F()
            sage: TG[[1,3],[2]]
            G[13|2]
            sage: TG[[1,3],[2]].to_fqsym()
            G[2, 1, 3] + G[3, 1, 2]
            sage: F(TG[[1,3],[2]].to_fqsym())
            F[2, 1, 3] + F[2, 3, 1]

        This embedding is a Hopf algebra morphism::

            sage: all((TG[t1] * TG[t2]).to_fqsym() == TG[t1].to_fqsym() * TG[t2].to_fqsym()
            ....:     for t1 in StandardTableaux(2)
            ....:     for t2 in StandardTableaux(3))
            True

            sage: T2 = TG.tensor_square()
            sage: phi = T2.module_morphism(
            ....:               lambda x: tensor([TG[x[0]].to_fqsym(), TG[x[1]].to_fqsym()]),
            ....:               codomain=F.tensor_square())
            sage: all(phi(TG[t].coproduct()) == F(TG[t].to_fqsym()).coproduct()
            ....:     for t in StandardTableaux(4))
            True

        There is a Hopf algebra map onto the Hopf algebra of symmetric
        functions, which maps a tableau `t` to the Schur function indexed
        by the shape of `t`::

            sage: TG = algebras.FSym(QQ).G()
            sage: t = StandardTableau([[1,3],[2,4],[5]])
            sage: TG[t]
            G[13|24|5]
            sage: TG[t].to_symmetric_function()
            s[2, 2, 1]
        """
        cat = GradedHopfAlgebras(base_ring).Connected()
        Parent.__init__(self, base=base_ring, category=cat.WithRealizations())

    _shorthands = ['G']

    def a_realization(self):
        r"""
        Return a particular realization of ``self`` (the Fundamental basis).

        EXAMPLES::

            sage: FSym = algebras.FSym(QQ)
            sage: FSym.a_realization()
            Hopf algebra of standard tableaux over the Rational Field
             in the Fundamental basis
        """
        return self.Fundamental()

    def dual(self):
        r"""
        Return the dual Hopf algebra of `FSym`.

        EXAMPLES::

            sage: algebras.FSym(QQ).dual()
            Dual Hopf algebra of standard tableaux over the Rational Field
        """
        return FreeSymmetricFunctions_Dual(self.base_ring())

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: algebras.FSym(QQ)
            Hopf algebra of standard tableaux over the Rational Field
        """
        return "Hopf algebra of standard tableaux over the %s" % self.base_ring()

    class Fundamental(FSymBasis_abstract):
        r"""
        The Hopf algebra of tableaux on the Fundamental basis.

        EXAMPLES::

            sage: FSym = algebras.FSym(QQ)
            sage: TG = FSym.G()
            sage: TG
            Hopf algebra of standard tableaux over the Rational Field
             in the Fundamental basis

        Elements of the algebra look like::

            sage: TG.an_element()
            2*G[] + 2*G[1] + 3*G[12]

        TESTS::

            sage: FSym = algebras.FSym(QQ)
            sage: TG = FSym.G()
            sage: TestSuite(TG).run()
        """
        _prefix = "G"
        _basis_name = "Fundamental"

        def dual_basis(self):
            r"""
            Return the dual basis to ``self``.

            EXAMPLES::

                sage: G = algebras.FSym(QQ).G()
                sage: G.dual_basis()
                Dual Hopf algebra of standard tableaux over the Rational Field
                 in the FundamentalDual basis
            """
            return self.realization_of().dual().F()

        @cached_method
        def product_on_basis(self, t1, t2):
            r"""
            EXAMPLES::

                sage: FSym = algebras.FSym(QQ)
                sage: G = FSym.G()
                sage: t1 = StandardTableau([[1,2], [3]])
                sage: t2 = StandardTableau([[1,2,3]])
                sage: G.product_on_basis(t1, t2)
                G[12456|3] + G[1256|3|4] + G[1256|34] + G[126|35|4]

                sage: t1 = StandardTableau([[1],[2]])
                sage: t2 = StandardTableau([[1,2]])
                sage: G.product_on_basis(t1, t2)
                G[134|2] + G[14|2|3]

                sage: t1 = StandardTableau([[1,2],[3]])
                sage: t2 = StandardTableau([[1],[2]])
                sage: G.product_on_basis(t1, t2)
                G[12|3|4|5] + G[12|34|5] + G[124|3|5] + G[124|35]
            """
            n = t1.size()
            m = t1.size() + t2.size()
            tableaux = []
            for t in StandardTableaux(m):
                if t.restrict(n) == t1 and standardize(t.anti_restrict(n).rectify()) == t2:
                    tableaux.append(t)
            return self.sum_of_monomials(tableaux)

        @cached_method
        def coproduct_on_basis(self, t):
            r"""
            EXAMPLES::

                sage: FSym = algebras.FSym(QQ)
                sage: G = FSym.G()
                sage: t = StandardTableau([[1,2,5], [3,4]])
                sage: G.coproduct_on_basis(t)
                G[] # G[125|34] + G[1|2] # G[13|2] + G[12] # G[12|3]
                 + G[12] # G[123] + G[1] # G[12|34] + G[1] # G[124|3]
                 + G[123] # G[12] + G[13|2] # G[1|2] + G[13|2] # G[12]
                 + G[12|34] # G[1] + G[134|2] # G[1] + G[125|34] # G[]
            """
            # we use the duality to compute this
            n = t.size()
            L = []
            dual_basis = self.dual_basis()
            for i in range(n + 1):
                for t1 in StandardTableaux(i):
                    for t2 in StandardTableaux(n - i):
                        coeff = (dual_basis[t1] * dual_basis[t2])[t]
                        if coeff:
                            L.append(((t1, t2), coeff))
            TT = self.tensor_square()
            return TT.sum_of_terms(L)

        class Element(FSymBasis_abstract.Element):
            def to_fqsym(self):
                """
                Return the image of ``self`` in the natural map to `FQSym`.

                EXAMPLES::

                    sage: FSym = algebras.FSym(QQ)
                    sage: G = FSym.G()
                    sage: t = StandardTableau([[1,3],[2,4],[5]])
                    sage: G[t].to_fqsym()
                    G[2, 1, 5, 4, 3] + G[3, 1, 5, 4, 2] + G[3, 2, 5, 4, 1]
                     + G[4, 1, 5, 3, 2] + G[4, 2, 5, 3, 1]
                """
                from sage.combinat.fqsym import FreeQuasisymmetricFunctions
                R = self.parent().base_ring()
                G = FreeQuasisymmetricFunctions(R).G()
                P = G._indices
                def on_basis(t):
                    return G.sum_of_monomials(P(sigma) for sigma in Permutations(t.size())
                                              if sigma.right_tableau() == t)
                return G.linear_combination((on_basis(t), coeff) for t, coeff in self)

            def to_symmetric_function(self):
                """
                Return the image of ``self`` in the natural map to `Sym`.

                EXAMPLES::

                    sage: FSym = algebras.FSym(QQ)
                    sage: G = FSym.G()
                    sage: t = StandardTableau([[1,3],[2,4],[5]])
                    sage: G[t].to_symmetric_function()
                    s[2, 2, 1]
                """
                s = SymmetricFunctions(self.parent().base_ring()).s()
                return s.sum_of_terms((t.shape(), coeff) for t, coeff in self)

    G = Fundamental

class FreeSymmetricFunctions_Dual(UniqueRepresentation, Parent):
    """
    The Hopf dual `FSym^*` of the free symmetric functions `FSym`.
    """
    def __init__(self, base_ring):
        r"""
        Initialize ``self``.

        TESTS::

            sage: FSymD = algebras.FSym(QQ).dual()
            sage: TestSuite(FSymD).run()
        """
        # TODO: the following line won't be needed when CategoryObject won't override base_ring
        self._base = base_ring
        Parent.__init__(self, category=GradedHopfAlgebras(base_ring).WithRealizations())

    _shorthands = ['F']

    def a_realization(self):
        r"""
        Return a particular realization of ``self`` (the Fundamental
        dual basis).

        EXAMPLES::

            sage: FSym = algebras.FSym(QQ).dual()
            sage: FSym.a_realization()
            Dual Hopf algebra of standard tableaux over the Rational Field
             in the FundamentalDual basis
        """
        return self.FundamentalDual()

    def dual(self):
        r"""
        Return the dual Hopf algebra of ``self``, which is `FSym`.

        EXAMPLES::

            sage: D = algebras.FSym(QQ).dual()
            sage: D.dual()
            Hopf algebra of standard tableaux over the Rational Field
        """
        return FreeSymmetricFunctions(self.base_ring())

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: algebras.FSym(QQ).dual()
            Dual Hopf algebra of standard tableaux over the Rational Field
        """
        return "Dual Hopf algebra of standard tableaux over the %s" % self.base_ring()

    class FundamentalDual(FSymBasis_abstract):
        r"""
        The Hopf algebra of tableaux on the fundamental dual basis.

        EXAMPLES::

            sage: FSym = algebras.FSym(QQ)
            sage: TF = FSym.dual().F()
            sage: TF
            Dual Hopf algebra of standard tableaux over the Rational Field
             in the FundamentalDual basis

        Elements of the algebra look like::

            sage: TF.an_element()
            3*F[12] + 2*F[1] + 2*F[]

        TESTS::

            sage: FSym = algebras.FSym(QQ)
            sage: TF = FSym.dual().F()
            sage: TestSuite(TF).run()
        """
        _prefix = "F"
        _basis_name = "FundamentalDual"

        def _coerce_map_from_(self, R):
            r"""
            Return ``True`` if there is a coercion from ``R`` into ``self``
            and ``False`` otherwise.

            The things that coerce into ``self`` are

            - word quasi-symmetric functions over a base with
              a coercion map into ``self.base_ring()``
            - symmetric functions over a base with a coercion
              map into ``self.base_ring()``
            - free quasi-symmetric functions over a base with
              a coercion map into ``self.base_ring()``

            EXAMPLES:

            `FSym` is a Hopf-quotient algbera of `FQSym`: the basis
            element indexed by a permutation `\sigma` is mapped to
            the tableau `Q(\sigma)`::

                sage: TF = algebras.FSym(QQ).dual().F()
                sage: SF = algebras.FQSym(QQ).F()
                sage: TF(SF([3,1,4,5,2]))
                F[134|25]

            This mapping is a Hopf algebra morphism::

                sage: all(TF(SF[p1] * SF[p2]) == TF(SF[p1]) * TF(SF[p2])
                ....:     for p1 in Permutations(2)
                ....:     for p2 in Permutations(3))
                True

                sage: SSym2 = SF.tensor_square()
                sage: phi = SSym2.module_morphism(
                ....:               lambda x: tensor([TF(SF[x[0]]), TF(SF[x[1]])]),
                ....:               codomain=TF.tensor_square())
                sage: all(phi(SF[p].coproduct()) == TF(SF[p]).coproduct()
                ....:     for p in Permutations(4))
                True

            There is a Hopf algebra morphism from `Sym`::

                sage: Sym = SymmetricFunctions(QQ)
                sage: s = Sym.schur()
                sage: TF = algebras.FSym(QQ).dual().F()
                sage: TF(s[2,1])
                F[12|3] + F[13|2]
                sage: TF(s[2,2,1])
                F[12|34|5] + F[12|35|4] + F[13|24|5] + F[13|25|4] + F[14|25|3]
            """
            if hasattr(R, "realization_of"):
                if not self.base_ring().has_coerce_map_from(R.base_ring()):
                    return False
                A = R.realization_of()
                # FQSym to FSym^*
                from sage.combinat.fqsym import FreeQuasisymmetricFunctions
                if isinstance(A, FreeQuasisymmetricFunctions):
                    F = A.F()
                    def F_to_SF_on_basis(sigma):
                        return self.monomial(sigma.right_tableau())
                    phi = F.module_morphism(F_to_SF_on_basis, codomain=self)
                    if R is F:
                        return phi
                    return phi * F.coerce_map_from(R)

                # Sym to FSym^*
                if isinstance(A, SymmetricFunctions):
                    s = A.s()
                    def s_to_F_on_basis(mu):
                        return self.sum_of_monomials(StandardTableaux(mu))
                    phi = s.module_morphism(s_to_F_on_basis, codomain=self)
                    if R is s:
                        return phi
                    return phi * s.coerce_map_from(R)
            return super(FreeSymmetricFunctions_Dual.FundamentalDual, self)._coerce_map_from_(R)

        def dual_basis(self):
            r"""
            Return the dual basis to ``self``.

            EXAMPLES::

                sage: F = algebras.FSym(QQ).dual().F()
                sage: F.dual_basis()
                Hopf algebra of standard tableaux over the Rational Field
                 in the Fundamental basis
            """
            return self.realization_of().dual().G()

        @cached_method
        def product_on_basis(self, t1, t2):
            r"""
            EXAMPLES::

                sage: FSym = algebras.FSym(QQ)
                sage: TF = FSym.dual().F()
                sage: t1 = StandardTableau([[1,2]])
                sage: TF.product_on_basis(t1, t1)
                F[12|34] + F[123|4] + F[1234] + F[124|3] + F[13|24] + F[134|2]
            """
            z = []
            n = t1.size()
            m = t2.size()
            for (I, J) in OrderedSetPartitions(range(1, n + m + 1), [n, m]):
                tt1 = [[I[x-1] for x in row] for row in t1]
                tt2 = [tuple([J[x-1] for x in row]) for row in t2]
                z.append(StandardTableau(Tableau(tt1).slide_multiply(tt2)))
            return self.sum_of_monomials(z)

        @cached_method
        def coproduct_on_basis(self, t):
            r"""
            EXAMPLES::

                sage: FSym = algebras.FSym(QQ)
                sage: TF = FSym.dual().F()
                sage: t = StandardTableau([[1,2,5], [3,4]])
                sage: TF.coproduct_on_basis(t)
                F[] # F[125|34] + F[1] # F[134|2] + F[12] # F[123] + F[12|3] # F[12] + F[12|34] # F[1] + F[125|34] # F[]
            """
            terms = [(t.restrict(i), standardize(t.anti_restrict(i).rectify()))
                        for i in range(t.size()+1)]
            return self.tensor_square().sum_of_monomials(terms)

        class Element(FSymBasis_abstract.Element):
            def to_quasisymmetric_function(self):
                """
                Return the image of ``self`` as a quasi-symmetric function
                in the Fundamental basis.

                EXAMPLES::

                    sage: F = algebras.FSym(QQ).dual().F()
                    sage: F[[1,3,5],[2,4]].to_quasisymmetric_function()
                    F[1, 2, 2]
                """
                from sage.combinat.ncsf_qsym.qsym import QuasiSymmetricFunctions
                QF = QuasiSymmetricFunctions(self.base_ring()).Fundamental()
                return QF.sum_of_terms((descent_composition(t), coeff)
                                       for t, coeff in self)

    F = FundamentalDual

##### utility functions for tableaux

# TODO: add these methods to tableaux!
from sage.combinat.composition import Composition

def standardize(t):
    A = sorted(t.to_word())
    std = {j:(i+1) for (i, j) in enumerate(A)}
    return StandardTableau([[std[i] for i in row] for row in t])


def ascent_set(t):
    locations = {}
    for (i,row) in enumerate(t):
        for (j,entry) in enumerate(row):
            locations[entry] = (i,j)
    ascents = [t.size()]
    for i in range(1,t.size()):
        # ascent means i+1 appears to the right or above
        x, y = locations[i]
        u, v = locations[i+1]
        if u <= x:
            ascents.append(i)
    return sorted(ascents)

def descent_set(t):
    ascents = set(ascent_set(t))
    return [i for i in range(1,1+t.size()) if i not in ascents]

def descent_composition(t):
    return Composition(from_subset=(descent_set(t), t.size()))


r"""
There are various morphisms between the Hopf algebras below.
We test that the diagram of morphisms is commutative::

    sage: NSym = NonCommutativeSymmetricFunctions(QQ)
    sage: QSym = QuasiSymmetricFunctions(QQ)
    sage: FQSym = algebras.FQSym(QQ)
    sage: FSym = algebras.FSym(QQ)
    sage: FSymDual = FSym.dual()
    sage: Sym = SymmetricFunctions(QQ)

    sage: def go(composition):
    ....:     x = NSym.a_realization()[composition]
    ....:     if QSym(FQSym(x)) != QSym(Sym(x)): return False
    ....:     if Sym(FSym(x)) != Sym(x): return False
    ....:     if FQSym(FSym(x)) != FQSym(x): return False
    ....:     if FSymDual(Sym(x)) != FSymDual(FQSym(x)): return False
    ....:     return True

    sage: go([2,1,2])
    True
    sage: all(all(go(comp) for comp in Compositions(n)) for n in range(5))
    True

    sage: def go2(n):
    ....:     for sigma in Permutations(n):
    ....:          x = FQSym.F()[sigma]
    ....:          if QSym(FSymDual(x)) != QSym(x): return False
    ....:     for mu in Partitions(n):
    ....:          x = Sym.s()[mu]
    ....:          if QSym(FSymDual(x)) != QSym(x): return False
    ....:     return True

    sage: all(go2(n) for n in range(6))
    True
"""
