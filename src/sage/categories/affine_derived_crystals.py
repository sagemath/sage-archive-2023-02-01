r"""
Affine Derived Subalgebra Crystals
"""
#*****************************************************************************
#  Copyright (C) 2015   Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.categories.category_singleton import Category_singleton
from sage.categories.crystals import Crystals
from sage.categories.regular_crystals import RegularCrystals
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.tensor import TensorProductsCategory
from sage.graphs.dot2tex_utils import have_dot2tex
from sage.functions.other import ceil
from sage.rings.all import ZZ


class AffineDerivedSubalgebraCrystals(Category_singleton):
    r"""
    The category of `U_q'(\mathfrak{g})`-crystals, where `\mathfrak{g}`
    is of affine type.

    EXAMPLES::

        sage: from sage.categories.affine_derived_crystals import AffineDerivedSubalgebraCrystals
        sage: C = AffineDerivedSubalgebraCrystals()
        sage: C
        Category of affine derived subalgebra crystals
        sage: C.super_categories()
        [Category of crystals]
        sage: C.example()
        Kirillov-Reshetikhin crystal of type ['A', 3, 1] with (r,s)=(1,1)

    TESTS::

        sage: TestSuite(C).run()
        sage: B = FiniteCrystals().example()
        sage: TestSuite(B).run(verbose = True)
        running ._test_an_element() . . . pass
        running ._test_cardinality() . . . pass
        running ._test_category() . . . pass
        running ._test_elements() . . .
          Running the test suite of self.an_element()
          running ._test_category() . . . pass
          running ._test_eq() . . . pass
          running ._test_not_implemented_methods() . . . pass
          running ._test_pickling() . . . pass
          running ._test_stembridge_local_axioms() . . . pass
          pass
        running ._test_elements_eq_reflexive() . . . pass
        running ._test_elements_eq_symmetric() . . . pass
        running ._test_elements_eq_transitive() . . . pass
        running ._test_elements_neq() . . . pass
        running ._test_enumerated_set_contains() . . . pass
        running ._test_enumerated_set_iter_cardinality() . . . pass
        running ._test_enumerated_set_iter_list() . . . pass
        running ._test_eq() . . . pass
        running ._test_fast_iter() . . . pass
        running ._test_not_implemented_methods() . . . pass
        running ._test_pickling() . . . pass
        running ._test_some_elements() . . . pass
        running ._test_stembridge_local_axioms() . . . pass
    """
    @cached_method
    def super_categories(self):
        r"""
        EXAMPLES::

            sage: from sage.categories.affine_derived_crystals import AffineDerivedSubalgebraCrystals
            sage: AffineDerivedSubalgebraCrystals().super_categories()
            [Category of crystals]
        """
        return [Crystals()]

    def example(self, n = 3):
        """
        Returns an example of Kirillov-Reshetikhin crystals, as per
        :meth:`Category.example`.

        EXAMPLES::

            sage: from sage.categories.affine_derived_crystals import AffineDerivedSubalgebraCrystals
            sage: B = AffineDerivedSubalgebraCrystals().example(); B
            Kirillov-Reshetikhin crystal of type ['A', 3, 1] with (r,s)=(1,1)
        """
        from sage.combinat.crystals.kirillov_reshetikhin import KirillovReshetikhinCrystal
        return KirillovReshetikhinCrystal(['A', n, 1], 1, 1)

    class ParentMethods:
        def weight_lattice_realization(self):
            """
            Return the weight lattice realization used to express weights
            of elements in ``self``.

            The default is to use the non-extended affine weight lattice.

            EXAMPLES::

                sage: C = crystals.Letters(['A', 5])
                sage: C.weight_lattice_realization()
                Ambient space of the Root system of type ['A', 5]
                sage: K = crystals.KirillovReshetikhin(['A',2,1], 1, 1)
                sage: K.weight_lattice_realization()
                Weight lattice of the Root system of type ['A', 2, 1]
            """
            F = self.cartan_type().root_system()
            return F.weight_lattice(extended=False)

        def digraph(self, subset=None, index_set=None):
            r"""
            Return the :class:`DiGraph` associated to ``self``.

            INPUT:

            - ``subset`` -- (optional) a subset of vertices for
              which the digraph should be constructed

            - ``index_set`` -- (optional) the index set to draw arrows

            .. SEEALSO::

                :meth:`sage.categories.crystals.Crystals.ParentMethods.digraph`

            EXAMPLES::

                sage: C = crystals.KirillovReshetikhin(['D',4,1], 2, 1)
                sage: G = C.digraph()
                sage: G.latex_options()  # optional - dot2tex
                LaTeX options for Digraph on 29 vertices:
                {...'edge_options': <function <lambda> at 0x...>,...}
                sage: view(G, tightpage=True)  # optional - dot2tex graphviz, not tested (opens external window)
            """
            G = Crystals().parent_class.digraph(self, subset, index_set)
            if have_dot2tex():
                f = lambda u_v_label: ({"backward": u_v_label[2] == 0})
                G.set_latex_options(edge_options=f)
            return G

# TODO: Should we make "regular" an axiom?
class RegularAffineDerivedSubalgebraCrystals(Category_singleton):
    r"""
    The category of regular `U_q'(\mathfrak{g})`-crystals, where
    `\mathfrak{g}` is of affine type.
    """
    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: from sage.categories.affine_derived_crystals import RegularAffineDerivedSubalgebraCrystals
            sage: RegularAffineDerivedSubalgebraCrystals().super_categories()
            [Category of regular crystals,
             Category of affine derived subalgebra crystals]
        """
        return [RegularCrystals(), AffineDerivedSubalgebraCrystals()]

    class ElementMethods:
        def classical_weight(self):
            """
            Return the classical weight of ``self``.

            EXAMPLES::

                sage: R = RootSystem(['A',2,1])
                sage: La = R.weight_space().basis()
                sage: LS = crystals.ProjectedLevelZeroLSPaths(2*La[1])
                sage: hw = LS.classically_highest_weight_vectors()
                sage: [(v.weight(), v.classical_weight()) for v in hw]
                [(-2*Lambda[0] + 2*Lambda[1], (2, 0, 0)),
                 (-Lambda[0] + Lambda[2], (1, 1, 0))]
            """
            CT = self.cartan_type().classical()
            I0 = CT.index_set()
            La = CT.root_system().ambient_space().fundamental_weights()
            return sum(La[i] * (self.phi(i) - self.epsilon(i)) for i in I0)

class KirillovReshetikhinCrystals(Category_singleton):
    """
    Category of Kirillov-Reshetikhin crystals.
    """
    @cached_method
    def super_categories(self):
        r"""
        EXAMPLES::

            sage: from sage.categories.affine_derived_crystals import KirillovReshetikhinCrystals
            sage: KirillovReshetikhinCrystals().super_categories()
            [Category of finite regular affine derived subalgebra crystals]
        """
        return [RegularAffineDerivedSubalgebraCrystals().Finite()]

    class ParentMethods:
        @abstract_method
        def r(self):
            r"""
            Return the value `r` in ``self`` written as `B^{r,s}`.

            EXAMPLES::

                sage: K = crystals.KirillovReshetikhin(['A',3,1], 2,4)
                sage: K.r()
                2
            """

        @abstract_method
        def s(self):
            r"""
            Return the value `s` in ``self`` written as `B^{r,s}`.

            EXAMPLES::

                sage: K = crystals.KirillovReshetikhin(['A',3,1], 2,4)
                sage: K.s()
                4
            """

        @abstract_method(optional=True)
        def classical_decomposition(self):
            """
            Return the classical decomposition of ``self``.

            EXAMPLES::

                sage: K = crystals.KirillovReshetikhin(['A',3,1], 2,2)
                sage: K.classical_decomposition()
                The crystal of tableaux of type ['A', 3] and shape(s) [[2, 2]]
            """

        @cached_method
        def classically_highest_weight_vectors(self):
            """
            Return the classically highest weight elements of ``self``.

            EXAMPLES::

                sage: K = crystals.KirillovReshetikhin(['E',6,1],1,1)
                sage: K.classically_highest_weight_vectors()
                ([(1,)],)
            """
            I0 = self.cartan_type().classical().index_set()
            return tuple([x for x in self if x.is_highest_weight(I0)])

        # TODO: This is duplicated in tensor product category
        def cardinality(self):
            """
            Return the cardinality of ``self``.

            EXAMPLES::

                sage: K = crystals.KirillovReshetikhin(['E',6,1], 1,1)
                sage: K.cardinality()
                27
                sage: K = crystals.KirillovReshetikhin(['C',6,1], 4,3)
                sage: K.cardinality()
                4736732
            """
            CWLR = self.cartan_type().classical().root_system().ambient_space()
            return sum(CWLR.weyl_dimension(mg.classical_weight())
                       for mg in self.classically_highest_weight_vectors())

        @cached_method
        def maximal_vector(self):
            r"""
            Return the unique element of classical weight `s \Lambda_r`
            in ``self``.

            EXAMPLES::

                sage: K = crystals.KirillovReshetikhin(['C',2,1],1,2)
                sage: K.maximal_vector()
                [[1, 1]]
                sage: K = crystals.KirillovReshetikhin(['E',6,1],1,1)
                sage: K.maximal_vector()
                [(1,)]

                sage: K = crystals.KirillovReshetikhin(['D',4,1],2,1)
                sage: K.maximal_vector()
                [[1], [2]]
            """
            R = self.weight_lattice_realization()
            Lambda = R.fundamental_weights()
            r = self.r()
            s = self.s()
            weight = s*Lambda[r] - s*Lambda[0] * Lambda[r].level() / Lambda[0].level()

            # First check the module generators as it is likely to be in here
            for b in self.module_generators:
                if b.weight() == weight:
                    return b

            # Otherwise check all of the elements
            for b in self:
                if b not in self.module_generators and b.weight() == weight:
                    return b

            assert False, "BUG: invalid Kirillov-Reshetikhin crystal"

        # TODO: Should this be in one of the super categories?
        def affinization(self):
            """
            Return the corresponding affinization crystal of ``self``.

            EXAMPLES::

                sage: K = crystals.KirillovReshetikhin(['A',2,1], 1, 1)
                sage: K.affinization()
                Affinization of Kirillov-Reshetikhin crystal of type ['A', 2, 1] with (r,s)=(1,1)

                sage: K = crystals.KirillovReshetikhin(['A',2,1], 1, 1, model='KR')
                sage: K.affinization()
                Affinization of Kirillov-Reshetikhin tableaux of type ['A', 2, 1] and shape (1, 1)
            """
            from sage.combinat.crystals.affinization import AffinizationOfCrystal
            return AffinizationOfCrystal(self)

        @cached_method
        def R_matrix(self, K):
            r"""
            Return the combinatorial `R`-matrix of ``self`` to ``K``.

            The *combinatorial* `R`-*matrix* is the affine crystal
            isomorphism `R : L \otimes K \to K \otimes L` which maps
            `u_{L} \otimes u_K` to `u_K \otimes u_{L}`, where `u_K`
            is the unique element in `K = B^{r,s}` of weight
            `s\Lambda_r - s c \Lambda_0` (see :meth:`module_generator`).

            INPUT:

            - ``self`` -- a crystal `L`
            - ``K`` -- a Kirillov-Reshetikhin crystal of the same type as `L`

            EXAMPLES::

                sage: K = crystals.KirillovReshetikhin(['A',2,1],1,1)
                sage: L = crystals.KirillovReshetikhin(['A',2,1],1,2)
                sage: f = K.R_matrix(L)
                sage: [[b,f(b)] for b in crystals.TensorProduct(K,L)]
                [[[[[1]], [[1, 1]]], [[[1, 1]], [[1]]]],
                 [[[[1]], [[1, 2]]], [[[1, 1]], [[2]]]],
                 [[[[1]], [[2, 2]]], [[[1, 2]], [[2]]]],
                 [[[[1]], [[1, 3]]], [[[1, 1]], [[3]]]],
                 [[[[1]], [[2, 3]]], [[[1, 2]], [[3]]]],
                 [[[[1]], [[3, 3]]], [[[1, 3]], [[3]]]],
                 [[[[2]], [[1, 1]]], [[[1, 2]], [[1]]]],
                 [[[[2]], [[1, 2]]], [[[2, 2]], [[1]]]],
                 [[[[2]], [[2, 2]]], [[[2, 2]], [[2]]]],
                 [[[[2]], [[1, 3]]], [[[2, 3]], [[1]]]],
                 [[[[2]], [[2, 3]]], [[[2, 2]], [[3]]]],
                 [[[[2]], [[3, 3]]], [[[2, 3]], [[3]]]],
                 [[[[3]], [[1, 1]]], [[[1, 3]], [[1]]]],
                 [[[[3]], [[1, 2]]], [[[1, 3]], [[2]]]],
                 [[[[3]], [[2, 2]]], [[[2, 3]], [[2]]]],
                 [[[[3]], [[1, 3]]], [[[3, 3]], [[1]]]],
                 [[[[3]], [[2, 3]]], [[[3, 3]], [[2]]]],
                 [[[[3]], [[3, 3]]], [[[3, 3]], [[3]]]]]

                sage: K = crystals.KirillovReshetikhin(['D',4,1],1,1)
                sage: L = crystals.KirillovReshetikhin(['D',4,1],2,1)
                sage: f = K.R_matrix(L)
                sage: T = crystals.TensorProduct(K,L)
                sage: b = T( K(rows=[[1]]), L(rows=[]) )
                sage: f(b)
                [[[2], [-2]], [[1]]]

            Alternatively, one can compute the combinatorial `R`-matrix
            using the isomorphism method of digraphs::

                sage: K1 = crystals.KirillovReshetikhin(['A',2,1],1,1)
                sage: K2 = crystals.KirillovReshetikhin(['A',2,1],2,1)
                sage: T1 = crystals.TensorProduct(K1,K2)
                sage: T2 = crystals.TensorProduct(K2,K1)
                sage: T1.digraph().is_isomorphic(T2.digraph(), edge_labels = True, certify = True) #todo: not implemented (see #10904 and #10549)
                (True, {[[[1]], [[2], [3]]]: [[[1], [3]], [[2]]], [[[3]], [[2], [3]]]: [[[2], [3]], [[3]]],
                [[[3]], [[1], [3]]]: [[[1], [3]], [[3]]], [[[1]], [[1], [3]]]: [[[1], [3]], [[1]]], [[[1]],
                [[1], [2]]]: [[[1], [2]], [[1]]], [[[2]], [[1], [2]]]: [[[1], [2]], [[2]]], [[[3]],
                [[1], [2]]]: [[[2], [3]], [[1]]], [[[2]], [[1], [3]]]: [[[1], [2]], [[3]]], [[[2]], [[2], [3]]]: [[[2], [3]], [[2]]]})
            """
            from sage.combinat.crystals.tensor_product import TensorProductOfCrystals
            T1 = TensorProductOfCrystals(self, K)
            T2 = TensorProductOfCrystals(K, self)
            gen1 = T1( self.maximal_vector(), K.maximal_vector() )
            gen2 = T2( K.maximal_vector(), self.maximal_vector() )
            g = { gen1 : gen2 }
            return T1.crystal_morphism(g, check=False)

        def is_perfect(self, ell=None):
            r"""
            Check if ``self`` is a perfect crystal of level ``ell``.

            A crystal `\mathcal{B}` is perfect of level `\ell` if:

            #. `\mathcal{B}` is isomorphic to the crystal graph of a
               finite-dimensional `U_q^{'}(\mathfrak{g})`-module.
            #. `\mathcal{B}\otimes \mathcal{B}` is connected.
            #. There exists a `\lambda\in X`, such that
               `\mathrm{wt}(\mathcal{B}) \subset \lambda + \sum_{i\in I}
                \ZZ_{\le 0} \alpha_i` and there is a unique element in
               `\mathcal{B}` of classical weight `\lambda`.
            #. For all `b \in \mathcal{B}`,
               `\mathrm{level}(\varepsilon (b)) \geq \ell`.
            #. For all `\Lambda` dominant weights of level `\ell`, there
               exist unique elements `b_{\Lambda}, b^{\Lambda} \in
               \mathcal{B}`, such that `\varepsilon ( b_{\Lambda}) =
               \Lambda = \varphi( b^{\Lambda})`.

            Points (1)-(3) are known to hold. This method checks
            points (4) and (5).

            If ``self`` is the Kirillov-Reshetikhin crystal `B^{r,s}`,
            then it was proven for non-exceptional types in [FOS2010]_
            that it is perfect if and only if `s/c_r` is an integer
            (where `c_r` is a constant related to the type of the crystal).

            It is conjectured this is true for all affine types.

            INPUT:

            - ``ell`` -- (default: `s / c_r`) integer; the level

            REFERENCES:

            .. [FOS2010] G. Fourier, M. Okado, A. Schilling. *Perfectness of
               Kirillov-Reshetikhin crystals for nonexceptional types*.
               Contemp. Math. 506 (2010) 127-143 ( :arxiv:`0811.1604` )

            EXAMPLES::

                sage: K = crystals.KirillovReshetikhin(['A',2,1], 1, 1)
                sage: K.is_perfect()
                True

                sage: K = crystals.KirillovReshetikhin(['C',2,1], 1, 1)
                sage: K.is_perfect()
                False

                sage: K = crystals.KirillovReshetikhin(['C',2,1], 1, 2)
                sage: K.is_perfect()
                True

            .. TODO::

                Implement a version for tensor products of KR crystals.
            """
            if ell is None:
                ell = self.s() / self.cartan_type().c()[self.r()]
                if ell not in ZZ:
                    return False

            if ell not in ZZ:
                raise ValueError("perfectness not defined for non-integral levels")

            # [FOS2010]_ check
            if self.cartan_type().classical().type() not in ['E','F','G']:
                return ell == self.s() / self.cartan_type().c()[self.r()]

            # Check by definition
            # TODO: This is duplicated from ProjectedLevelZeroLSPaths, combine the two methods.
            # TODO: Similarly, don't duplicate in the tensor product category, maybe
            #   move this to the derived affine category?
            MPhi = []
            for b in self:
                p = b.Phi().level()
                assert p == b.Epsilon().level()
                if p < level:
                    return False
                if p == level:
                    MPhi += [b]
            weights = []
            I = self.index_set()
            rank = len(I)
            La = self.weight_lattice_realization().basis()
            from sage.combinat.integer_vector import IntegerVectors
            for n in range(1,level+1):
                for c in IntegerVectors(n, rank):
                    w = sum(c[i]*La[i] for i in I)
                    if w.level() == level:
                        weights.append(w)
            return sorted(b.Phi() for b in MPhi) == sorted(weights)

        def level(self):
            r"""
            Return the level of ``self`` when ``self`` is a perfect crystal.

            .. SEEALSO::

                :meth:`~sage.categories.affine_derived_crystals.KirillovReshetikhinCrystals.ParentMethods.is_perfect`

            EXAMPLES::

                sage: K = crystals.KirillovReshetikhin(['A',2,1], 1, 1)
                sage: K.level()
                1
                sage: K = crystals.KirillovReshetikhin(['C',2,1], 1, 2)
                sage: K.level()
                1
                sage: K = crystals.KirillovReshetikhin(['D',4,1], 1, 3)
                sage: K.level()
                3

                sage: K = crystals.KirillovReshetikhin(['C',2,1], 1, 1)
                sage: K.level()
                Traceback (most recent call last):
                ...
                ValueError: this crystal is not perfect
            """
            if not self.is_perfect():
                raise ValueError("this crystal is not perfect")
            return self.s() / self.cartan_type().c()[self.r()]

        def q_dimension(self, q=None, prec=None, use_product=False):
            """
            Return the `q`-dimension of ``self``.

            The `q`-dimension of a KR crystal is defined as the `q`-dimension of
            the underlying classical crystal.

            EXAMPLES::

                sage: KRC = crystals.KirillovReshetikhin(['A',2,1], 2,2)
                sage: KRC.q_dimension()
                q^4 + q^3 + 2*q^2 + q + 1
                sage: KRC = crystals.KirillovReshetikhin(['D',4,1], 2,1)
                sage: KRC.q_dimension()
                q^10 + q^9 + 3*q^8 + 3*q^7 + 4*q^6 + 4*q^5 + 4*q^4 + 3*q^3 + 3*q^2 + q + 2
            """
            return self.classical_decomposition().q_dimension(q, prec, use_product)

    class ElementMethods:
        def lusztig_involution(self):
            r"""
            Return the result of the classical Lusztig involution on ``self``.

            EXAMPLES::

                sage: KRT = crystals.KirillovReshetikhin(['D',4,1], 2, 3, model='KR')
                sage: mg = KRT.module_generators[1]
                sage: mg.lusztig_involution()
                [[-2, -2, 1], [-1, -1, 2]]
                sage: elt = mg.f_string([2,1,3,2]); elt
                [[3, -2, 1], [4, -1, 2]]
                sage: elt.lusztig_involution()
                [[-4, -2, 1], [-3, -1, 2]]
            """
            Cl = self.parent().cartan_type().classical()
            I = Cl.index_set()
            aut = Cl.opposition_automorphism()
            hw = self.to_highest_weight(I)[1]
            hw.reverse()
            return self.to_lowest_weight(I)[0].e_string(aut[i] for i in hw)

    class TensorProducts(TensorProductsCategory):
        """
        The category of tensor products of Kirillov-Reshetikhin crystals.
        """
        @cached_method
        def extra_super_categories(self):
            """
            EXAMPLES::

                sage: from sage.categories.affine_derived_crystals import KirillovReshetikhinCrystals
                sage: KirillovReshetikhinCrystals().TensorProducts().extra_super_categories()
                [Category of finite regular affine derived subalgebra crystals]
            """
            return [RegularAffineDerivedSubalgebraCrystals().Finite()]

        class ParentMethods:
            @cached_method
            def maximal_vector(self):
                """
                Return the maximal vector of ``self``.

                EXAMPLES::

                    sage: K = crystals.KirillovReshetikhin(['A',2,1],1,1)
                    sage: T = crystals.TensorProduct(K,K,K)
                    sage: T.maximal_vector()
                    [[[1]], [[1]], [[1]]]
                """
                return self(*[K.maximal_vector() for K in self.crystals])

            @cached_method
            def classically_highest_weight_vectors(self):
                """
                Return the classically highest weight elements of ``self``.

                This works by using a backtracking algorithm since if
                `b_2 \otimes b_1` is classically highest weight then `b_1`
                is classically highest weight.

                EXAMPLES::

                    sage: K = crystals.KirillovReshetikhin(['A',2,1],1,1)
                    sage: T = crystals.TensorProduct(K,K,K)
                    sage: T.classically_highest_weight_vectors()
                    ([[[1]], [[1]], [[1]]],
                     [[[2]], [[1]], [[1]]],
                     [[[1]], [[2]], [[1]]],
                     [[[3]], [[2]], [[1]]])
                """
                n = len(self.crystals)
                I0 = self.cartan_type().classical().index_set()
                it = [ iter(self.crystals[-1].classically_highest_weight_vectors()) ]
                path = []
                ret = []
                while it:
                    try:
                        x = next(it[-1])
                    except StopIteration:
                        it.pop()
                        if path:
                            path.pop(0)
                        continue

                    b = self.element_class(self, [x] + path)
                    if not b.is_highest_weight(index_set=I0):
                        continue
                    path.insert(0, x)
                    if len(path) == n:
                        ret.append(b)
                        path.pop(0)
                    else:
                        it.append( iter(self.crystals[-len(path)-1]) )
                return tuple(ret)

            # TODO: This is duplicated in KR crystals category
            def cardinality(self):
                """
                Return the cardinality of ``self``.

                EXAMPLES::

                    sage: RC = RiggedConfigurations(['A', 3, 1], [[3, 2], [1, 2]])
                    sage: RC.cardinality()
                    100
                    sage: len(RC.list())
                    100

                    sage: RC = RiggedConfigurations(['E', 7, 1], [[1,1]])
                    sage: RC.cardinality()
                    134
                    sage: len(RC.list())
                    134

                    sage: RC = RiggedConfigurations(['B', 3, 1], [[2,2],[1,2]])
                    sage: RC.cardinality()
                    5130
                """
                CWLR = self.cartan_type().classical().root_system().ambient_space()
                return sum(CWLR.weyl_dimension(mg.classical_weight())
                           for mg in self.classically_highest_weight_vectors())

            def one_dimensional_configuration_sum(self, q=None, group_components=True):
                r"""
                Compute the one-dimensional configuration sum of ``self``.

                INPUT:

                - ``q`` -- (default: ``None``) a variable or ``None``;
                  if ``None``, a variable `q` is set in the code
                - ``group_components`` -- (default: ``True``) boolean; if
                  ``True``, then the terms are grouped by classical component

                The one-dimensional configuration sum is the sum of the
                weights of all elements in the crystal weighted by the
                energy function.

                EXAMPLES::

                    sage: K = crystals.KirillovReshetikhin(['A',2,1],1,1)
                    sage: T = crystals.TensorProduct(K,K)
                    sage: T.one_dimensional_configuration_sum()
                    B[-2*Lambda[1] + 2*Lambda[2]] + (q+1)*B[-Lambda[1]]
                     + (q+1)*B[Lambda[1] - Lambda[2]] + B[2*Lambda[1]]
                     + B[-2*Lambda[2]] + (q+1)*B[Lambda[2]]
                    sage: R.<t> = ZZ[]
                    sage: T.one_dimensional_configuration_sum(t, False)
                    B[-2*Lambda[1] + 2*Lambda[2]] + (t+1)*B[-Lambda[1]]
                     + (t+1)*B[Lambda[1] - Lambda[2]] + B[2*Lambda[1]]
                     + B[-2*Lambda[2]] + (t+1)*B[Lambda[2]]

                    sage: R = RootSystem(['A',2,1])
                    sage: La = R.weight_space().basis()
                    sage: LS = crystals.ProjectedLevelZeroLSPaths(2*La[1])
                    sage: LS.one_dimensional_configuration_sum() == T.one_dimensional_configuration_sum() # long time
                    True

                TESTS::

                    sage: K1 = crystals.KirillovReshetikhin(['A',2,1],1,1)
                    sage: K2 = crystals.KirillovReshetikhin(['A',2,1],2,1)
                    sage: T = crystals.TensorProduct(K1,K2)
                    sage: T.one_dimensional_configuration_sum() == T.one_dimensional_configuration_sum(group_components=False)
                    True
                """
                if q is None:
                    from sage.rings.all import QQ
                    q = QQ['q'].gens()[0]
                P0 = self.weight_lattice_realization().classical()
                B = P0.algebra(q.parent())
                if group_components:
                    G = self.digraph(index_set=self.cartan_type().classical().index_set())
                    C = G.connected_components()
                    return B.sum(q**(c[0].energy_function())*B.sum(B(P0(b.weight())) for b in c)
                                 for c in C)
                return B.sum(q**(b.energy_function())*B(P0(b.weight())) for b in self)

        class ElementMethods:
            def energy_function(self):
                r"""
                Return the energy function of ``self``.

                The energy is only defined when ``self`` is an element of a
                tensor product of affine Kirillov-Reshetikhin crystals. In this
                implementation, it is assumed that ``self`` is an element of a
                tensor product of perfect crystals of the same level, see
                Theorem 7.5 in [SchillingTingley2011]_.

                INPUT:

                - ``self`` -- an element of a tensor product of perfect
                  Kirillov-Reshetkhin crystals of the same level

                OUTPUT: an integer

                REFERENCES:

                .. [SchillingTingley2011] A. Schilling, P. Tingley.
                   *Demazure crystals, Kirillov-Reshetikhin crystals, and
                   the energy function*.
                   Electronic Journal of Combinatorics. **19(2)**. 2012.
                   :arXiv:`1104.2359`

                EXAMPLES::

                    sage: K = crystals.KirillovReshetikhin(['A',2,1],1,1)
                    sage: T = crystals.TensorProduct(K,K,K)
                    sage: hw = [b for b in T if all(b.epsilon(i)==0 for i in [1,2])]
                    sage: for b in hw:
                    ....:    print b, b.energy_function()
                    [[[1]], [[1]], [[1]]] 0
                    [[[1]], [[2]], [[1]]] 2
                    [[[2]], [[1]], [[1]]] 1
                    [[[3]], [[2]], [[1]]] 3

                    sage: K = crystals.KirillovReshetikhin(['C',2,1],1,2)
                    sage: T = crystals.TensorProduct(K,K)
                    sage: hw = [b for b in T if all(b.epsilon(i)==0 for i in [1,2])]
                    sage: for b in hw:  # long time (5s on sage.math, 2011)
                    ....:     print b, b.energy_function()
                    [[], []] 4
                    [[], [[1, 1]]] 1
                    [[[1, 1]], []] 3
                    [[[1, 1]], [[1, 1]]] 0
                    [[[1, 2]], [[1, 1]]] 1
                    [[[2, 2]], [[1, 1]]] 2
                    [[[-1, -1]], [[1, 1]]] 2
                    [[[1, -1]], [[1, 1]]] 2
                    [[[2, -1]], [[1, 1]]] 2

                    sage: K = crystals.KirillovReshetikhin(['C',2,1],1,1)
                    sage: T = crystals.TensorProduct(K)
                    sage: t = T.module_generators[0]
                    sage: t.energy_function()
                    Traceback (most recent call last):
                    ...
                    NotImplementedError: all crystals in the tensor product need to be perfect of the same level
                """
                C = self.parent().crystals[0]
                ell = ceil(C.s()/C.cartan_type().c()[C.r()])
                if any(ell != K.s()/K.cartan_type().c()[K.r()] for K in self.parent().crystals):
                    raise NotImplementedError("all crystals in the tensor product need to be perfect of the same level")
                t = self.parent()(*[K.module_generator() for K in self.parent().crystals])
                d = t.affine_grading()
                return d - self.affine_grading()

            def affine_grading(self):
                r"""
                Return the affine grading of ``self``.

                The affine grading is only defined when ``self`` is an
                element of a tensor product of Kirillov-Reshetikhin
                crystals. It is calculated by finding a path from ``self``
                to a ground state path using the helper method
                :meth:`e_string_to_ground_state` and counting the number
                of affine Kashiwara operators `e_0` applied on the way.

                OUTPUT: an integer

                EXAMPLES::

                    sage: K = crystals.KirillovReshetikhin(['A',2,1],1,1)
                    sage: T = crystals.TensorProduct(K,K)
                    sage: t = T.module_generators[0]
                    sage: t.affine_grading()
                    1

                    sage: K = crystals.KirillovReshetikhin(['A',2,1],1,1)
                    sage: T = crystals.TensorProduct(K,K,K)
                    sage: hw = [b for b in T if all(b.epsilon(i)==0 for i in [1,2])]
                    sage: for b in hw:
                    ....:     print b, b.affine_grading()
                    [[[1]], [[1]], [[1]]] 3
                    [[[1]], [[2]], [[1]]] 1
                    [[[2]], [[1]], [[1]]] 2
                    [[[3]], [[2]], [[1]]] 0

                    sage: K = crystals.KirillovReshetikhin(['C',2,1],1,1)
                    sage: T = crystals.TensorProduct(K,K,K)
                    sage: hw = [b for b in T if all(b.epsilon(i)==0 for i in [1,2])]
                    sage: for b in hw:
                    ....:     print b, b.affine_grading()
                    [[[1]], [[1]], [[1]]] 2
                    [[[1]], [[2]], [[1]]] 1
                    [[[1]], [[-1]], [[1]]] 0
                    [[[2]], [[1]], [[1]]] 1
                    [[[-2]], [[2]], [[1]]] 0
                    [[[-1]], [[1]], [[1]]] 1
                """
                return self.e_string_to_ground_state().count(0)

            @cached_method
            def e_string_to_ground_state(self):
                r"""
                Return a string of integers in the index set
                `(i_1, \ldots, i_k)` such that `e_{i_k} \cdots e_{i_1}`
                of ``self`` is the ground state.

                This method calculates a path from ``self`` to a ground
                state path using Demazure arrows as defined in Lemma 7.3
                in [SchillingTingley2011]_.

                OUTPUT: a tuple of integers `(i_1, \ldots, i_k)`

                EXAMPLES::

                    sage: K = crystals.KirillovReshetikhin(['A',2,1],1,1)
                    sage: T = crystals.TensorProduct(K,K)
                    sage: t = T.module_generators[0]
                    sage: t.e_string_to_ground_state()
                    (0, 2)

                    sage: K = crystals.KirillovReshetikhin(['C',2,1],1,1)
                    sage: T = crystals.TensorProduct(K,K)
                    sage: t = T.module_generators[0]; t
                    [[[1]], [[1]]]
                    sage: t.e_string_to_ground_state()
                    (0,)
                    sage: x=t.e(0)
                    sage: x.e_string_to_ground_state()
                    ()
                    sage: y=t.f_string([1,2,1,1,0]); y
                    [[[2]], [[1]]]
                    sage: y.e_string_to_ground_state()
                    ()
                """
                I = self.cartan_type().classical().index_set()
                ell = max(ceil(K.s()/K.cartan_type().c()[K.r()]) for K in self.parent().crystals)
                for i in I:
                    if self.epsilon(i) > 0:
                        return (i,) + (self.e(i)).e_string_to_ground_state()
                if self.epsilon(0) > ell:
                    return (0,) + (self.e(0)).e_string_to_ground_state()
                return ()

