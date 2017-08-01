r"""
Regular SuperCrystals
"""
#*****************************************************************************
#  Copyright (C) 2017    Travis Scrimshaw <tcscrims at gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#****************************************************************************
from __future__ import print_function

from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.categories.category_singleton import Category_singleton
from sage.categories.crystals import Crystals
from sage.categories.tensor import TensorProductsCategory
from sage.combinat.subset import Subsets
from sage.graphs.dot2tex_utils import have_dot2tex

class RegularSuperCrystals(Category_singleton):
    def super_categories(self):
        r"""
        EXAMPLES::

            sage: from sage.categories.regular_supercrystals import RegularSuperCrystalCategory
            sage: RegularSuperCrystalCategory().super_categories()
            [Category of finite crystals]
        """
        return [Crystals().Finite()]

    class ParentMethods:
        def m(self):
            r"""
            EXAMPLES::

                sage: from bkk_crystals import BKKOneBoxCrystal
                sage: c = BKKOneBoxCrystal(2, 3)
                sage: c.m()
                2

            Test for tensor products::

                sage: t = c.tensor(c)
                sage: t.m()
                2

            Test for direct sums::

                sage: s1, s2 = t.connected_components()
                sage: s = s1 + s2
                sage: s.m()
                2
            """
            return 1 - min(self.index_set())

        def n(self):
            r"""
            EXAMPLES::

                sage: from bkk_crystals import BKKOneBoxCrystal
                sage: c = BKKOneBoxCrystal(2, 3)
                sage: c.n()
                3

            Test for tensor products::

                sage: t = c.tensor(c)
                sage: t.n()
                3

            Test for direct sums::

                sage: s1, s2 = t.connected_components()
                sage: s = s1 + s2
                sage: s.n()
                3
            """
            return 1 + max(self.index_set())

        @cached_method
        def digraph(self):
            r"""
            EXAMPLES::

                sage: from bkk_crystals import BKKOneBoxCrystal
                sage: c = BKKOneBoxCrystal(2, 3)
                sage: c.digraph()
                Digraph on 5 vertices
            """
            from sage.graphs.digraph import DiGraph
            from sage.misc.latex import LatexExpr

            G = DiGraph()
            for i in self.index_set():
                for x in self:
                    y = x.f(i)
                    if y is not None:
                        G.add_edge(x, y, i)

            def edge_options((u, v, l)):
                edge_color_by_label = { 0: 'black', 1: 'blue', 2: 'red', 3: 'green' }
                edge_opts = { 'edge_string': '->', 'color': 'black' }
                if l > 0:
                    edge_opts['color'] = edge_color_by_label[l]
                    edge_opts['label'] = LatexExpr(str(l))
                elif l < 0:
                    edge_opts['color'] = "dashed," + edge_color_by_label[-l]
                    edge_opts['label'] = LatexExpr("\\overline{%s}" % str(-l))
                else:
                    edge_opts['color'] = "dotted," + edge_color_by_label[l]
                    edge_opts['label'] = LatexExpr(str(l))
                return edge_opts

            G.set_latex_options(format="dot2tex", edge_labels=True, edge_options=edge_options)

            return G

        def connected_components_generators(self):
            r"""
            EXAMPLES::

                sage: from bkk_crystals import BKKOneBoxCrystal
                sage: c = BKKOneBoxCrystal(2, 3)
                sage: c.connected_components_generators()
                [(-2,)]

                sage: t = c.tensor(c)
                sage: t.connected_components_generators()
                [([-2, -1],), ([-2, -2],)]

                sage: t = c.tensor(c)
                sage: s1, s2 = t.connected_components()
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
            EXAMPLES::

                sage: from bkk_crystals import BKKOneBoxCrystal
                sage: c = BKKOneBoxCrystal(2, 3)
                sage: c.connected_components()
                [Subcrystal of BKK crystal on semistandard tableaux of shape [1] with entries in (-2, -1, 1, 2, 3)]
                sage: t = c.tensor(c)
                sage: t.connected_components()
                [Subcrystal of <class 'bkk_crystals.TensorProductOfSuperCrystals_with_category'>,
                 Subcrystal of <class 'bkk_crystals.TensorProductOfSuperCrystals_with_category'>]
            """
            category = RegularSuperCrystalCategory()
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
            r"""
            EXAMPLES::

                sage: from bkk_crystals import BKKOneBoxCrystal
                sage: c = BKKOneBoxCrystal(2, 3)
                sage: c.tensor(c)
                <class 'bkk_crystals.TensorProductOfSuperCrystals_with_category'>
            """
            from sage.combinat.crystals.bkk_crystals import TensorProductOfSuperCrystals
            return TensorProductOfSuperCrystals((self,) + crystals, **options)

        def direct_sum(self, X):
            r"""
            EXAMPLES::

                sage: from bkk_crystals import BKKOneBoxCrystal
                sage: c = BKKOneBoxCrystal(2, 3)
                sage: t = c.tensor(c)
                sage: s1, s2 = t.connected_components()
                sage: s1 + s2
                Direct sum of the crystals Family (Subcrystal of <class 'bkk_crystals.TensorProductOfSuperCrystals_with_category'>, Subcrystal of <class 'bkk_crystals.TensorProductOfSuperCrystals_with_category'>)
            """
            from sage.combinat.crystals.direct_sum import DirectSumOfCrystals
            return DirectSumOfCrystals([self, X])

    class ElementMethods:
        def epsilon(self, i):
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

