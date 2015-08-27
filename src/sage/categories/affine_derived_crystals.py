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

        sage: C = FiniteCrystals()
        sage: C
        Category of finite crystals
        sage: C.super_categories()
        [Category of crystals, Category of finite enumerated sets]
        sage: C.example()
        Highest weight crystal of type A_3 of highest weight omega_1

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

            sage: FiniteCrystals().extra_super_categories()
            [Category of finite enumerated sets]
        """
        return [RegularCrystals().Finite()]

    def example(self, n = 3):
        """
        Returns an example of Kirillov-Reshetikhin crystals, as per
        :meth:`Category.example`.

        EXAMPLES::

            sage: B = KirillovReshetikhinCrystals().example(); B
            Highest weight crystal of type A_3 of highest weight omega_1
        """
        from sage.combinat.crystals.kirillov_reshetikhin import KirillovReshetikhinCrystal
        return KirillovReshetkhinCrystal(['A', n, 1], 1, 1)

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

            sage: RegularCrystals().super_categories()
            [Category of crystals]
        """
        return [AffineDerivedSubalgebraCrystals()]

    class ParentMethods:
        @abstract_method
        def r(self):
            """
            Return the value `r` in ``self`` written as `B^{r,s}`.
            """

        @abstract_method
        def s(self):
            """
            Return the value `s` in ``self`` written as `B^{r,s}`.
            """

    class TensorProducts(TensorProductsCategory):
        """
        The category of tensor products of Kirillov-Reshetikhin crystals.
        """

