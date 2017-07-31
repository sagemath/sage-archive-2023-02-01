#######################################################################
#                   tensor product of BKK crystals                    #
#######################################################################

from sage.categories.cartesian_product import cartesian_product
from sage.categories.crystals import Crystals
from sage.categories.enumerated_sets import EnumeratedSets
from sage.combinat.crystals.direct_sum import DirectSumOfCrystals
from sage.combinat.root_system.cartan_type import CartanType
from sage.graphs.digraph import DiGraph
from sage.misc.cachefunc import cached_method
from sage.misc.latex import LatexExpr
from sage.misc.misc import uniq
from sage.structure.element_wrapper import ElementWrapper
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.lazy_attribute import lazy_attribute


from sage.combinat.crystals.tensor_product import TensorProductOfCrystals, TensorProductOfCrystalsElement
from sage.categories.category_singleton import Category_singleton

class BKKCrystalCategory(Category_singleton):
    def super_categories(self):
        r"""
        EXAMPLES::

            sage: from bkk_crystals import BKKCrystalCategory
            sage: BKKCrystalCategory().super_categories()
            [Category of finite crystals]

        """
        return [Crystals().Finite()]

    class ParentMethods:

        # NOTE: we need to define _index_set because some classes override the
        # index_set method (such as subcrystals, which returns self._index_set)
        @lazy_attribute
        def _index_set(self):
            r"""
            EXAMPLES::

                sage: from bkk_crystals import BKKOneBoxCrystal
                sage: c = BKKOneBoxCrystal(2, 3)
                sage: c._index_set
                (-1, 0, 1, 2)

            Test that the index set is computed correctly for tensor products::

                sage: t = c.tensor(c)
                sage: t._index_set
                (-1, 0, 1, 2)

            Test that the index set is computed correctly for direct sums::

                sage: s1, s2 = t.connected_components()
                sage: s = s1 + s2
                sage: s._index_set
                (-1, 0, 1, 2)
            """
            if hasattr(self, 'crystals'):
                # FIXME: where should this code go? it works for tensor
                # products and direct sum (and anything with a crystals method)
                ems = uniq([c.m() for c in self.crystals])
                assert len(ems) == 1
                m = ems[0]

                ens = uniq([c.n() for c in self.crystals])
                assert len(ens) == 1
                n = ens[0]

                return tuple(range(-m + 1, n))
            else:
                return NotImplemented

        def index_set(self):
            r"""
            EXAMPLES::

                sage: from bkk_crystals import BKKOneBoxCrystal
                sage: c = BKKOneBoxCrystal(2, 3)
                sage: c.index_set()
                (-1, 0, 1, 2)

            Test that the index set is computed correctly for tensor products::

                sage: t = c.tensor(c)
                sage: t.index_set()
                (-1, 0, 1, 2)

            Test that the index set is computed correctly for direct sums::

                sage: s1, s2 = t.connected_components()
                sage: s = s1 + s2
                sage: s.index_set()
                (-1, 0, 1, 2)
            """
            return tuple(self._index_set)

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
                [Subcrystal of <class 'bkk_crystals.TensorProductOfBKKCrystals_with_category'>,
                 Subcrystal of <class 'bkk_crystals.TensorProductOfBKKCrystals_with_category'>]
            """
            category = BKKCrystalCategory()
            index_set = self.index_set()
            CCs = []

            # FIXME: setting the cartan type and then deleting it is a hack!
            # the subcrystal init method insists on the cartan type, but
            # I don't want to define it; hopefully, this doesn't break anything
            cartan_type = CartanType("A", max(self.n(), self.m()) - 1)
            for mg in self.connected_components_generators():
                if not isinstance(mg, tuple):
                    mg = (mg,)
                subcrystal = self.subcrystal(generators=mg,
                                             index_set=index_set,
                                             cartan_type=cartan_type,
                                             category=category)
                subcrystal._cartan_type = None
                CCs.append(subcrystal)

            return CCs

        def tensor(self, *crystals, **options):
            r"""
            EXAMPLES::

                sage: from bkk_crystals import BKKOneBoxCrystal
                sage: c = BKKOneBoxCrystal(2, 3)
                sage: c.tensor(c)
                <class 'bkk_crystals.TensorProductOfBKKCrystals_with_category'>
            """
            return TensorProductOfBKKCrystals((self,) + crystals, **options)

        def direct_sum(self, X):
            r"""
            EXAMPLES::

                sage: from bkk_crystals import BKKOneBoxCrystal
                sage: c = BKKOneBoxCrystal(2, 3)
                sage: t = c.tensor(c)
                sage: s1, s2 = t.connected_components()
                sage: s1 + s2
                Direct sum of the crystals Family (Subcrystal of <class 'bkk_crystals.TensorProductOfBKKCrystals_with_category'>, Subcrystal of <class 'bkk_crystals.TensorProductOfBKKCrystals_with_category'>)
            """
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

class BKKCrystalForVectorRepresentation(UniqueRepresentation, Parent):

    def __init__(self, m, n):
        Parent.__init__(self, category = BKKCrystalCategory())
        self._index_set = tuple(range(-m + 1, n))
        self.module_generators = [ t for t in self ]

    def __iter__(self):
        for t in range(-self.m(), self.n() + 1):
            if t != 0:
                yield self(t)

    def _repr_(self):
        return "BKK crystal on semistandard tableaux of shape [1] with entries in {}".format(tuple(self))

    # temporary workaround while an_element is overriden by Parent
    _an_element_ = EnumeratedSets.ParentMethods._an_element_

    class Element(ElementWrapper):

        def e(self, i):
            assert i in self.parent().index_set()
            b = self.value
            if i < 0:
                if b == i:
                    return self.parent()(b - 1)
            elif 0 < i:
                if b == i + 1:
                    return self.parent()(b - 1)
            elif i == 0 and b == 1:
                return self.parent()(-1)

        def f(self, i):
            assert i in self.parent().index_set()
            b = self.value
            if 0 < i:
                if b == i:
                    return self.parent()(b + 1)
            elif i < 0:
                if b == i - 1:
                    return self.parent()(b + 1)
            elif i == 0 and b == -1:
                return self.parent()(1)

        def weight(self):
            from sage.modules.free_module_element import vector
            elements = list(self.parent())
            v = vector([0]*len(elements))
            i = elements.index(self)
            v[i] = 1
            return v

BKKOneBoxCrystal = BKKCrystalForVectorRepresentation

class TensorProductOfBKKCrystals(TensorProductOfCrystals):
    r"""
    TESTS:

    In the BKK paper, they point out that there are elements that are highest
    weight vectors, but not "genuine". Here we verify their examples.

    First some import statements::

        sage: from bkk_crystals import BKKOneBoxCrystal
        sage: from bkk_tableaux import japanese_reading_word

    Construct the 6-fold tensor product of the BKK crystal for shape [1]::

        sage: c = BKKOneBoxCrystal(2, 2)
        sage: # TODO: extend tensor product defintion so that tensor([c]*6) works
        sage: c6 = reduce(lambda x, y: tensor((x,y)), [c]*6)

    Here is a function that finds the element of ``c6`` whose reading word is
    the Japanese reading word of the given tableau::

        sage: def find_element_by_japanese_reading_word(t):
        ....:     w = japanese_reading_word(t)
        ....:     s = str(reduce(lambda x, y : [x, y], w))
        ....:     for x in c6:
        ....:         if str(x) == s:
        ....:             return x

        sage: t = Tableau([ [-2,-2,-2], [-1,-1], [1] ])
        sage: ascii_art(t)
         -2 -2 -2
         -1 -1
          1
        sage: japanese_reading_word(t)
        word: -2,-2,-1,-2,-1,1
        sage: find_element_by_japanese_reading_word(t)
        [[[[[-2, -2], -1], -2], -1], 1]

    The following three elements are all highest weight vectors,
    and they all belong to the same connected component of ``c6``::

        sage: H0 = Tableau([ [-2,-2,-2], [-1,-1], [1] ])
        sage: x0 = find_element_by_japanese_reading_word(H0)
        sage: x0.is_highest_weight()
        True
        sage: x0.weight()
        (3, 2, 1, 0)

        sage: H1 = Tableau([ [-2,-2,-2], [-1,2], [1] ])
        sage: x1 = find_element_by_japanese_reading_word(H1)
        sage: x1.is_highest_weight()
        True
        sage: x1.weight()
        (3, 1, 1, 1)

        sage: H2 = Tableau([ [-2,-2,2], [-1,-1], [1] ])
        sage: x2 = find_element_by_japanese_reading_word(H2)
        sage: x2.is_highest_weight()
        True
        sage: x2.weight()
        (2, 2, 1, 1)

        sage: for cc in c6.connected_components():
        ....:    if x0 in cc:
        ....:        break
        sage: (x0 in cc) and (x1 in cc) and (x2 in cc)
        True

    Similarly, the following three elements are all lowest weight vectors,
    and they all belong to the same connected component of ``c6``::

        sage: L0 = Tableau([ [-1,1,2], [1,2], [2] ])
        sage: y0 = find_element_by_japanese_reading_word(L0)
        sage: y0.is_lowest_weight()
        True
        sage: y0.weight()
        (0, 1, 2, 3)

        sage: L1 = Tableau([ [-2,1,2], [-1,2], [2] ])
        sage: y1 = find_element_by_japanese_reading_word(L1)
        sage: y1.is_lowest_weight()
        True
        sage: y1.weight()
        (1, 1, 1, 3)

        sage: L2 = Tableau([ [-2,1,2], [-1,2], [1] ])
        sage: y2 = find_element_by_japanese_reading_word(L2)
        sage: y2.is_lowest_weight()
        True
        sage: y2.weight()
        (1, 1, 2, 2)

        sage: for cc in c6.connected_components():
        ....:    if y0 in cc:
        ....:        break
        sage: (y0 in cc) and (y1 in cc) and (y2 in cc)
        True

    """
    def __init__(self, crystals):
        assert isinstance(crystals, tuple)

        Parent.__init__(self, category=BKKCrystalCategory())
        self.crystals = crystals

        self.module_generators = self.list()

    def __iter__(self):
        for x in cartesian_product([c.list() for c in self.crystals]):
            yield self(*x)

class TensorProductOfBKKCrystalsElement(TensorProductOfCrystalsElement):

    def f(self, i):
        # FIXME: generalize this to deal with more than two tensor factors
        (b1, b2) = self

        # Proposition 2.8.ii.a
        if i < 0:
            if b1.phi(i) > b2.epsilon(i):
                x = b1.f(i)
                if x is not None:
                    return self.parent()(x, b2)
            else:
                y = b2.f(i)
                if y is not None:
                    return self.parent()(b1, y)

        # Proposition 2.8.ii.b
        elif 0 < i:
            if b2.phi(i) > b1.epsilon(i):
                y = b2.f(i)
                if y is not None:
                    return self.parent()(b1, y)
            else:
                x = b1.f(i)
                if x is not None:
                    return self.parent()(x, b2)

        # Proposition 2.8.ii.c
        elif i == 0:
            if b1.phi(i) + b1.epsilon(i):
                x = b1.f(i)
                if x is not None:
                    return self.parent()(x, b2)
            else:
                y = b2.f(i)
                if y is not None:
                    return self.parent()(b1, y)

    def e(self, i):
        (b1, b2) = self

        # Proposition 2.8.ii.a
        if i < 0:
            if b1.phi(i) >= b2.epsilon(i):
                x = b1.e(i)
                if x is not None:
                    return self.parent()(x, b2)
            else:
                y = b2.e(i)
                if y is not None:
                    return self.parent()(b1, y)

        # Proposition 2.8.ii.b
        elif 0 < i:
            if b2.phi(i) >= b1.epsilon(i):
                y = b2.e(i)
                if y is not None:
                    return self.parent()(b1, y)
            else:
                x = b1.e(i)
                if x is not None:
                    return self.parent()(x, b2)

        # Proposition 2.8.ii.c
        elif i == 0:
            if b1.phi(i) + b1.epsilon(i):
                x = b1.e(i)
                if x is not None:
                    return self.parent()(x, b2)
            else:
                y = b2.e(i)
                if y is not None:
                    return self.parent()(b1, y)

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

TensorProductOfBKKCrystals.Element = TensorProductOfBKKCrystalsElement
