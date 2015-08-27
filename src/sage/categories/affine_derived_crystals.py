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
#from sage.categories.category_with_axiom import CategoryWithAxiom
from sage.categories.category_singleton import Category_singleton
from sage.categories.crystals import Crystals
from sage.categories.regular_crystals import RegularCrystals
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.tensor import TensorProductsCategory
from sage.graphs.dot2tex_utils import have_dot2tex

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
        [Category of finite regular crystals]
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
            [Category of finite regular crystals]
        """
        return [RegularCrystals().Finite()]

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
            [Category of affine derived subalgebra crystals]
        """
        return [AffineDerivedSubalgebraCrystals().Finite()]

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

        @abstract_method
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
            """
            I0 = self.cartan_type().classical().index_set()
            return tuple([x for x in self if x.is_highest_weight(I0)])

        @cached_method
        def maximal_vector(self):
            r"""
            Return the unique element of classical weight `s \Lambda_r`
            in ``self``.

            EXAMPLES::

                sage: K = crystals.KirillovReshetikhin(['C',2,1],1,2)
                sage: K.module_generator()
                [[1, 1]]
                sage: K = crystals.KirillovReshetikhin(['E',6,1],1,1)
                sage: K.module_generator()
                [(1,)]

                sage: K = crystals.KirillovReshetikhin(['D',4,1],2,1)
                sage: K.module_generator()
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
        class ParentMethods:
            @cached_method
            def classically_highest_weight_vectors(self):
                """
                Return the classically highest weight elements of ``self``.

                This works by using a backtracking algorithm since if
                `b_2 \otimes b_1` is classically highest weight then `b_1`
                is classically highest weight.

                EXAMPLES::

                    sage: C = crystals.Tableaux(['D',4], shape=[2,2])
                    sage: D = crystals.Tableaux(['D',4], shape=[1])
                    sage: T = crystals.TensorProduct(D, C)
                    sage: T.highest_weight_vectors()
                    ([[[1]], [[1, 1], [2, 2]]],
                     [[[3]], [[1, 1], [2, 2]]],
                     [[[-2]], [[1, 1], [2, 2]]])
                    sage: L = filter(lambda x: x.is_highest_weight(), T)
                    sage: tuple(L) == T.highest_weight_vectors()
                    True
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
                    if not b.is_highest_weight():
                        continue
                    path.insert(0, x)
                    if len(path) == n:
                        ret.append(b)
                        path.pop(0)
                    else:
                        it.append( iter(self.crystals[-len(path)-1]) )
                return tuple(ret)

            def cardinality(self):
                """
                Return the cardinality of ``self``.

                EXAMPLES::

                    sage: RC = RiggedConfigurations(['A', 3, 1], [[3, 2], [1, 2]])
                    sage: RC.cardinality()
                    100
                    sage: RC = RiggedConfigurations(['B', 3, 1], [[2,2],[1,2]])
                    sage: RC.cardinality()
                    5130
                    sage: RC = RiggedConfigurations(['E', 7, 1], [[1,1]])
                    sage: RC.cardinality()
                    134
                """
                CWLR = self.cartan_type().classical().root_system().ambient_space()
                return sum(CWLR.weyl_dimension(mg.classical_weight())
                           for mg in self.classically_highest_weight_vectors())

