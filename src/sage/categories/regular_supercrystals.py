r"""
Regular Supercrystals
"""

#*****************************************************************************
#       Copyright (C) 2017 Franco Saliola <saliola@gmail.com>
#                     2017 Anne Schilling <anne at math.ucdavis.edu>
#                     2017 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import print_function

from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.categories.category_singleton import Category_singleton
from sage.categories.crystals import Crystals
from sage.categories.tensor import TensorProductsCategory
from sage.combinat.subset import Subsets
from sage.graphs.dot2tex_utils import have_dot2tex

class RegularSuperCrystals(Category_singleton):
    r"""
    The category of crystals for super Lie algebras.

    EXAMPLES::

        sage: from sage.categories.regular_supercrystals import RegularSuperCrystals
        sage: C = RegularSuperCrystals()
        sage: C
        Category of regular super crystals
        sage: C.super_categories()
        [Category of finite crystals]

    Parents in this category should implement the following methods:

    - either an attribute ``_cartan_type`` or a method ``cartan_type``

    - ``module_generators``: a list (or container) of distinct elements
      which generate the crystal using `f_i`

    Furthermore, their elements ``x`` should implement the following
    methods:

    - ``x.e(i)`` (returning `e_i(x)`)

    - ``x.f(i)`` (returning `f_i(x)`)

    - ``x.weight()`` (returning `weight(x)`)

    EXAMPLES::

        sage: from sage.misc.abstract_method import abstract_methods_of_class
        sage: from sage.categories.regular_supercrystals import RegularSuperCrystals
        sage: abstract_methods_of_class(RegularSuperCrystals().element_class)
        {'optional': [], 'required': ['e', 'f', 'weight']}

    TESTS::

        sage: from sage.categories.regular_supercrystals import RegularSuperCrystals
        sage: C = RegularSuperCrystals()
        sage: TestSuite(C).run()
        sage: B = crystals.Letters(['A',[1,1]]); B
        The crystal of letters for type ['A', [1, 1]]
        sage: TestSuite(B).run(verbose = True)
        running ._test_an_element() . . . pass
        running ._test_cardinality() . . . pass
        running ._test_category() . . . pass
        running ._test_elements() . . .
          Running the test suite of self.an_element()
          running ._test_category() . . . pass
          running ._test_eq() . . . pass
          running ._test_new() . . . pass
          running ._test_not_implemented_methods() . . . pass
          running ._test_pickling() . . . pass
          pass
        running ._test_elements_eq_reflexive() . . . pass
        running ._test_elements_eq_symmetric() . . . pass
        running ._test_elements_eq_transitive() . . . pass
        running ._test_elements_neq() . . . pass
        running ._test_enumerated_set_contains() . . . pass
        running ._test_enumerated_set_iter_cardinality() . . . pass
        running ._test_enumerated_set_iter_list() . . . pass
        running ._test_eq() . . . pass
        running ._test_new() . . . pass
        running ._test_not_implemented_methods() . . . pass
        running ._test_pickling() . . . pass
        running ._test_some_elements() . . . pass
    """
    def super_categories(self):
        r"""
        EXAMPLES::

            sage: from sage.categories.regular_supercrystals import RegularSuperCrystalCategory
            sage: RegularSuperCrystalCategory().super_categories()
            [Category of finite crystals]
        """
        return [Crystals().Finite()]

    class ParentMethods:
        @cached_method
        def digraph(self):
            r"""
            Return the :class:`DiGraph` associated to ``self``.

            INPUT:

            - ``subset`` -- (optional) a subset of vertices for
              which the digraph should be constructed

            - ``index_set`` -- (optional) the index set to draw arrows

            EXAMPLES::

                sage: B = crystals.Letters(['A', [1,3]])
                sage: G = B.digraph(); G
                Digraph on 6 vertices

            The edges of the crystal graph are by default colored using
            blue for edge 1, red for edge 2, green for edge 3, and dashed with
            the corresponding color for barred edges. Edge 0 is dotted black::

                sage: view(G)  # optional - dot2tex graphviz, not tested (opens external window)
            """
            from sage.graphs.digraph import DiGraph
            from sage.misc.latex import LatexExpr
            from sage.combinat.root_system.cartan_type import CartanType

            d = {x: {} for x in self}
            for i in self.index_set():
                for x in d:
                    y = x.f(i)
                    if y is not None:
                        d[x][y] = i
            G = DiGraph(d, format='dict_of_dicts')

            def edge_options((u, v, l)):
                edge_opts = { 'edge_string': '->', 'color': 'black' }
                if l > 0:
                    edge_opts['color'] = CartanType._colors.get(l, 'black')
                    edge_opts['label'] = LatexExpr(str(l))
                elif l < 0:
                    edge_opts['color'] = "dashed," + CartanType._colors.get(-l, 'black')
                    edge_opts['label'] = LatexExpr("\\overline{%s}" % str(-l))
                else:
                    edge_opts['color'] = "dotted," + CartanType._colors.get(l, 'black')
                    edge_opts['label'] = LatexExpr(str(l))
                return edge_opts

            G.set_latex_options(format="dot2tex", edge_labels=True, edge_options=edge_options)

            return G

        def connected_components_generators(self):
            r"""
            Return tuple of generators for each of the connected components of ``self``.

            EXAMPLES::

                sage: B = crystals.Letters(['A',[1,2]])
                sage: B.connected_components_generators()
                [(-2,)]

                sage: T = B.tensor(B)
                sage: T.connected_components_generators()
                [([-2, -1],), ([-2, -2],)]
                sage: s1, s2 = T.connected_components()
                sage: s = s1 + s2
                sage: s.connected_components_generators()
                [([-2, -1],), ([-2, -2],)]
            """
            # NOTE: we compute the connected components of the digraph,
            # then highest weight elements in each connected components
            X = []
            for connected_component_vertices in self.digraph().connected_components():
                gens = [g for g in connected_component_vertices if g.is_highest_weight()]
                X.append(tuple(gens))
            return X

        def connected_components(self):
            r"""
            Return the connected components of ``self`` as subcrystals.

            EXAMPLES::

                sage: B = crystals.Letters(['A',[1,2]])
                sage: B.connected_components()
                [Subcrystal of The crystal of letters for type ['A', [1, 2]]]

                sage: T = B.tensor(B)
                sage: T.connected_components()
                [Subcrystal of Full tensor product of the crystals 
                [The crystal of letters for type ['A', [1, 2]], 
                 The crystal of letters for type ['A', [1, 2]]],
                 Subcrystal of Full tensor product of the crystals 
                 [The crystal of letters for type ['A', [1, 2]], 
                 The crystal of letters for type ['A', [1, 2]]]]
            """
            category = RegularSuperCrystals()
            index_set = self.index_set()
            cartan_type = self.cartan_type()
            CCs = []

            for mg in self.connected_components_generators():
                if not isinstance(mg, tuple):
                    mg = (mg,)
                subcrystal = self.subcrystal(generators=mg,
                                             index_set=index_set,
                                             cartan_type=cartan_type,
                                             category=category)
                CCs.append(subcrystal)

            return CCs

        def tensor(self, *crystals, **options):
            """
            Return the tensor product of ``self`` with the crystals ``B``.

            EXAMPLES::

                sage: B = crystals.Letters(['A',[1,2]])
                sage: C = crystals.Tableaux(['A',[1,2]], shape = [2,1])
                sage: T = C.tensor(B); T
                Full tensor product of the crystals [Crystal of BKK tableaux of shape [2, 1] of gl(2|3), 
                The crystal of letters for type ['A', [1, 2]]]
                sage: S = B.tensor(C); S
                Full tensor product of the crystals [The crystal of letters for type ['A', [1, 2]], 
                Crystal of BKK tableaux of shape [2, 1] of gl(2|3)]
                sage: G = T.digraph()
                sage: H = S.digraph()
                sage: G.is_isomorphic(H, edge_labels= True)
                True
            """
            cartan_type = self.cartan_type()
            from sage.combinat.crystals.tensor_product import FullTensorProductOfSuperCrystals
            if any(c.cartan_type() != cartan_type for c in crystals):
                raise ValueError("all crystals must be of the same Cartan type")
            return FullTensorProductOfSuperCrystals((self,) + tuple(crystals), **options)

        def character(self):
            """
            Returns the character of ``self``.

            TODO: once the `WeylCharacterRing` is implemented, make this consistent
                  with the implementation in `sage.categories.classical_crystals.character`.

            EXAMPLES::

                sage: B = crystals.Letters(['A',[1,2]])
                sage: B.character()
                B[(1, 0, 0, 0, 0)] + B[(0, 1, 0, 0, 0)] + B[(0, 0, 1, 0, 0)] + B[(0, 0, 0, 1, 0)] 
                + B[(0, 0, 0, 0, 1)]
            """
            from sage.rings.all import ZZ
            A = self.weight_lattice_realization().algebra(ZZ)
            return A.sum(A(x.weight()) for x in self)

    class ElementMethods:

        def epsilon(self, i):
            """
            Return `\varepsilon_i` of ``self``.

            """
            string_length = 0
            x = self
            while True:
                x = x.e(i)
                if x is None:
                    return string_length
                else:
                    string_length += 1

        def phi(self, i):
            string_length = 0
            x = self
            while True:
                x = x.f(i)
                if x is None:
                    return string_length
                else:
                    string_length += 1

    class TensorProducts(TensorProductsCategory):
        """
        The category of regular crystals constructed by tensor
        product of regular crystals.
        """
        @cached_method
        def extra_super_categories(self):
            """
            EXAMPLES::

                sage: RegularCrystals().TensorProducts().extra_super_categories()
                [Category of regular crystals]
            """
            return [self.base_category()]

