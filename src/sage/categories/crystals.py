r"""
Crystals
"""
#*****************************************************************************
#  Copyright (C) 2010    Anne Schilling <anne at math.ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.abstract_method import abstract_method
from sage.misc.lazy_import import LazyImport
from sage.categories.category_singleton import Category_singleton
from sage.categories.enumerated_sets import EnumeratedSets
from sage.categories.tensor import TensorProductsCategory
from sage.categories.category import HomCategory
from sage.categories.morphism import Morphism
from sage.categories.homset import Hom, Homset
from sage.misc.latex import latex
from sage.combinat import ranker
from sage.graphs.dot2tex_utils import have_dot2tex
from sage.rings.integer import Integer
from sage.sets.recursively_enumerated_set import RecursivelyEnumeratedSet

class Crystals(Category_singleton):
    r"""
    The category of crystals.

    See :mod:`sage.combinat.crystals.crystals` for an introduction to crystals.

    EXAMPLES::

        sage: C = Crystals()
        sage: C
        Category of crystals
        sage: C.super_categories()
        [Category of... enumerated sets]
        sage: C.example()
        Highest weight crystal of type A_3 of highest weight omega_1

    Parents in this category should implement the following methods:

    - either an attribute ``_cartan_type`` or a method ``cartan_type``

    - ``module_generators``: a list (or container) of distinct elements
      which generate the crystal using `f_i`

    Furthermore, their elements should implement the following methods:

    - ``x.e(i)`` (returning `e_i(x)`)

    - ``x.f(i)`` (returning `f_i(x)`)

    - ``x.epsilon(i)`` (returning `\varepsilon_i(x)`)

    - ``x.phi(i)`` (returning `\varphi_i(x)`)

    EXAMPLES::

        sage: from sage.misc.abstract_method import abstract_methods_of_class
        sage: abstract_methods_of_class(Crystals().element_class)
        {'optional': [], 'required': ['e', 'epsilon', 'f', 'phi', 'weight']}

    TESTS::

        sage: TestSuite(C).run()
        sage: B = Crystals().example()
        sage: TestSuite(B).run(verbose = True)
        running ._test_an_element() . . . pass
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

    def super_categories(self):
        r"""
        EXAMPLES::

            sage: Crystals().super_categories()
            [Category of enumerated sets]
        """
        return [EnumeratedSets()]

    def example(self, choice="highwt", **kwds):
        r"""
        Returns an example of a crystal, as per
        :meth:`Category.example()
        <sage.categories.category.Category.example>`.

        INPUT::

         - ``choice`` -- str [default: 'highwt']. Can be either 'highwt'
           for the highest weight crystal of type A, or 'naive' for an
           example of a broken crystal.

         - ``**kwds`` -- keyword arguments passed onto the constructor for the
           chosen crystal.

        EXAMPLES::

            sage: Crystals().example(choice='highwt', n=5)
            Highest weight crystal of type A_5 of highest weight omega_1
            sage: Crystals().example(choice='naive')
            A broken crystal, defined by digraph, of dimension five.
        """
        import sage.categories.examples.crystals as examples
        if choice == "naive":
            return examples.NaiveCrystal(**kwds)
        else:
            from sage.rings.integer import Integer
            if isinstance(choice, Integer):
                return examples.HighestWeightCrystalOfTypeA(n=choice, **kwds)
            else:
                return examples.HighestWeightCrystalOfTypeA(**kwds)

    class HomCategory(HomCategory):
        """
        The category of homomorphisms sets `Hom(X,Y)` for crystals `X, Y`.
        """
        class ElementMethods:
            @abstract_method(optional=True)
            def is_isomorphism(self):
                """
                Check if ``self`` is a crystal isomorphism.

                EXAMPLES::

                    sage: B = crystals.Tableaux(['C',2], shape=[1,1])
                    sage: C = crystals.Tableaux(['C',2], ([2,1], [1,1]))
                    sage: psi = B.crystal_morphism(C.module_generators[1:], codomain=C)
                    sage: psi.is_isomorphism()
                    False
                """

            @abstract_method(optional=True)
            def is_embedding(self):
                """
                Check if ``self`` is an injective crystal morphism.

                EXAMPLES::

                    sage: B = crystals.Tableaux(['C',2], shape=[1,1])
                    sage: C = crystals.Tableaux(['C',2], ([2,1], [1,1]))
                    sage: psi = B.crystal_morphism(C.module_generators[1:], codomain=C)
                    sage: psi.is_embedding()
                    True
                """

            @abstract_method(optional=True)
            def is_strict(self):
                """
                Check if ``self`` is a strict crystal morphism.

                EXAMPLES::

                    sage: B = crystals.Tableaux(['C',2], shape=[1,1])
                    sage: C = crystals.Tableaux(['C',2], ([2,1], [1,1]))
                    sage: psi = B.crystal_morphism(C.module_generators[1:], codomain=C)
                    sage: psi.is_strict()
                    True
                """

    class ParentMethods:

        def an_element(self):
            """
            Returns an element of ``self``

                sage: C = crystals.Letters(['A', 5])
                sage: C.an_element()
                1
            """
            return self.first()

        def weight_lattice_realization(self):
            """
            Returns the weight lattice realization used to express weights.

            This default implementation uses the ambient space of the
            root system for (non relabelled) finite types and the
            weight lattice otherwise. This is a legacy from when
            ambient spaces were partially implemented, and may be
            changed in the future.

            EXAMPLES::

                sage: C = crystals.Letters(['A', 5])
                sage: C.weight_lattice_realization()
                Ambient space of the Root system of type ['A', 5]
                sage: K = crystals.KirillovReshetikhin(['A',2,1], 1, 1)
                sage: K.weight_lattice_realization()
                Weight lattice of the Root system of type ['A', 2, 1]
            """
            F = self.cartan_type().root_system()
            if F.is_finite() and F.ambient_space() is not None:
                return F.ambient_space()
            else:
                return F.weight_lattice()

        def cartan_type(self):
            """
            Returns the Cartan type of the crystal

            EXAMPLES::

                sage: C = crystals.Letters(['A',2])
                sage: C.cartan_type()
                ['A', 2]
            """
            return self._cartan_type

        @cached_method
        def index_set(self):
            """
            Returns the index set of the Dynkin diagram underlying the crystal

            EXAMPLES::

                sage: C = crystals.Letters(['A', 5])
                sage: C.index_set()
                (1, 2, 3, 4, 5)
            """
            return self.cartan_type().index_set()

        def Lambda(self):
            """
            Returns the fundamental weights in the weight lattice
            realization for the root system associated with the crystal

            EXAMPLES::

                sage: C = crystals.Letters(['A', 5])
                sage: C.Lambda()
                Finite family {1: (1, 0, 0, 0, 0, 0), 2: (1, 1, 0, 0, 0, 0), 3: (1, 1, 1, 0, 0, 0), 4: (1, 1, 1, 1, 0, 0), 5: (1, 1, 1, 1, 1, 0)}
            """
            return self.weight_lattice_realization().fundamental_weights()

        def __iter__(self, index_set=None, max_depth=float('inf')):
            """
            Return an iterator over the elements of ``self``.

            INPUT:

            - ``index_set`` -- (Default: ``None``) The index set; if ``None``
              then use the index set of the crystal

            - ``max_depth`` -- (Default: infinity) The maximum depth to build

            The iteration order is not specified except that, if
            ``max_depth`` is finite, then the iteration goes depth by
            depth.

            EXAMPLES::

                sage: C = crystals.LSPaths(['A',2,1],[-1,0,1])
                sage: C.__iter__.__module__
                'sage.categories.crystals'
                sage: g = C.__iter__()
                sage: next(g)
                (-Lambda[0] + Lambda[2],)
                sage: next(g)
                (Lambda[0] - Lambda[1] + delta,)
                sage: next(g)
                (Lambda[1] - Lambda[2] + delta,)
                sage: next(g)
                (-Lambda[0] + Lambda[2] + delta,)
                sage: next(g)
                (Lambda[1] - Lambda[2],)

                sage: sorted(C.__iter__(index_set=[1,2]), key=str)
                [(-Lambda[0] + Lambda[2],),
                 (Lambda[0] - Lambda[1],),
                 (Lambda[1] - Lambda[2],)]

                sage: sorted(C.__iter__(max_depth=1), key=str)
                [(-Lambda[0] + Lambda[2],),
                 (Lambda[0] - Lambda[1] + delta,),
                 (Lambda[1] - Lambda[2],)]

            """
            if index_set is None:
                index_set = self.index_set()
            if max_depth < float('inf'):
                return RecursivelyEnumeratedSet(self.module_generators,
                                                lambda x: [x.f(i) for i in index_set]
                                                        + [x.e(i) for i in index_set],
                                                structure=None, enumeration='breadth',
                                                max_depth=max_depth).__iter__()
            return RecursivelyEnumeratedSet(self.module_generators,
                                            lambda x: [x.f(i) for i in index_set]
                                                    + [x.e(i) for i in index_set],
                                            structure=None,
                                            enumeration='naive').__iter__()

        def subcrystal(self, index_set=None, generators=None, max_depth=float("inf"),
                       direction="both", contained=None,
                       virtualization=None, scaling_factors=None,
                       cartan_type=None, category=None):
            r"""
            Construct the subcrystal from ``generators`` using `e_i` and/or
            `f_i` for all `i` in ``index_set``.

            INPUT:

            - ``index_set`` -- (default: ``None``) the index set; if ``None``
              then use the index set of the crystal

            - ``generators`` -- (default: ``None``) the list of generators; if
              ``None`` then use the module generators of the crystal

            - ``max_depth`` -- (default: infinity) the maximum depth to build

            - ``direction`` -- (default: ``'both'``) the direction to build
              the subcrystal; it can be one of the following:

              - ``'both'`` - using both `e_i` and `f_i`
              - ``'upper'`` - using `e_i`
              - ``'lower'`` - using `f_i`

            - ``contained`` -- (optional) a set (or function) defining the
              containment in the subcrystal

            - ``virtualization``, ``scaling_factors`` -- (optional)
              dictionaries whose key `i` corresponds to the sets `\sigma_i`
              and `\gamma_i` respectively used to define virtual crystals; see
              :class:`~sage.combinat.crystals.virtual_crystal.VirtualCrystal`

            - ``cartan_type`` -- (optional) specify the Cartan type of the
              subcrystal

            - ``category`` -- (optional) specify the category of the subcrystal

            EXAMPLES::

                sage: C = crystals.KirillovReshetikhin(['A',3,1], 1, 2)
                sage: S = list(C.subcrystal(index_set=[1,2])); S
                [[[1, 1]], [[1, 2]], [[1, 3]], [[2, 2]], [[2, 3]], [[3, 3]]]
                sage: C.cardinality()
                10
                sage: len(S)
                6
                sage: list(C.subcrystal(index_set=[1,3], generators=[C(1,4)]))
                [[[1, 4]], [[2, 4]], [[2, 3]], [[1, 3]]]
                sage: list(C.subcrystal(index_set=[1,3], generators=[C(1,4)], max_depth=1))
                [[[1, 4]], [[2, 4]], [[1, 3]]]
                sage: list(C.subcrystal(index_set=[1,3], generators=[C(1,4)], direction='upper'))
                [[[1, 4]], [[1, 3]]]
                sage: list(C.subcrystal(index_set=[1,3], generators=[C(1,4)], direction='lower'))
                [[[1, 4]], [[2, 4]]]
            """
            from sage.combinat.crystals.subcrystal import Subcrystal
            from sage.categories.finite_crystals import FiniteCrystals

            if cartan_type is None:
                cartan_type = self.cartan_type()
            else:
                from sage.combinat.root_system.cartan_type import CartanType
                cartan_type = CartanType(cartan_type)
            if index_set is None:
                index_set = self.index_set()
            if generators is None:
                generators = self.module_generators

            if max_depth == float('inf'):
                if self not in FiniteCrystals():
                    if (contained is None and index_set == self.index_set()
                            and generators == self.module_generators
                            and scaling_factors is None and virtualization is None):
                        return self
                    return Subcrystal(self, contained, generators,
                                      virtualization, scaling_factors,
                                      cartan_type, index_set, category)
                if contained is None and direction == 'both':
                    if category is None:
                        category = FiniteCrystals()
                    return Subcrystal(self, contained, generators,
                                      virtualization, scaling_factors,
                                      cartan_type, index_set, category)

            # TODO: Make this work for virtual crystals as well
            if direction == 'both':
                subset = RecursivelyEnumeratedSet(generators, 
                                                  lambda x: [x.f(i) for i in index_set]
                                                          + [x.e(i) for i in index_set],
                                                  structure=None, enumeration='breadth',
                                                  max_depth=max_depth)
            elif direction == 'upper':
                subset = RecursivelyEnumeratedSet(generators, 
                                                  lambda x: [x.e(i) for i in index_set],
                                                  structure=None, enumeration='breadth',
                                                  max_depth=max_depth)
            elif direction == 'lower':
                subset = RecursivelyEnumeratedSet(generators, 
                                                  lambda x: [x.f(i) for i in index_set],
                                                  structure=None, enumeration='breadth',
                                                  max_depth=max_depth)
            else:
                raise ValueError("direction must be either 'both', 'upper', or 'lower'")

            # We perform the filtering here since checking containment
            #   in a frozenset should be fast
            if contained is not None:
                subset = frozenset(x for x in subset if contained(x))
            else:
                subset = frozenset(subset)

            if category is None:
                category = FiniteCrystals()
            else:
               category = FiniteCrystals().join(category)

            if self in FiniteCrystals() and len(subset) == self.cardinality():
                if index_set == self.index_set():
                    return self
                return Subcrystal(self, subset, generators,
                                  virtualization, scaling_factors,
                                  cartan_type, index_set, category)

            return Subcrystal(self, subset, generators,
                              virtualization, scaling_factors,
                              cartan_type, index_set, category)

        def _Hom_(self, Y, category=None, **options):
            r"""
            Return the homset from ``self`` to ``Y`` in the
            category ``category``.

            INPUT::

            - ``Y`` -- a crystal
            - ``category`` -- a subcategory of :class:`Crysals`() or ``None``

            The sole purpose of this method is to construct the homset
            as a :class:`~sage.categories.crystals.CrystalHomset`. If
            ``category`` is specified and is not a subcategory of
            :class:`Crystals`, a ``TypeError`` is raised instead

            This method is not meant to be called directly. Please use
            :func:`sage.categories.homset.Hom` instead.

            EXAMPLES::

                sage: B = crystals.elementary.B(['A',2], 1)
                sage: H = B._Hom_(B); H
                Set of Crystal Morphisms from The 1-elementary crystal of type ['A', 2]
                 to The 1-elementary crystal of type ['A', 2]
            """
            if category is None:
                category = self.category()
            elif not category.is_subcategory(Crystals()):
                raise TypeError("{} is not a subcategory of Crystals()".format(category))
            if Y not in Crystals():
                raise TypeError("{} is not a crystal".format(Y))
            return CrystalHomset(self, Y, category=category, **options)

        def crystal_morphism(self, on_gens, codomain=None,
                             cartan_type=None, index_set=None, generators=None,
                             automorphism=None,
                             virtualization=None, scaling_factors=None,
                             category=None, check=True):
            r"""
            Construct a crystal morphism from ``self`` to another crystal
            ``codomain``.

            INPUT:

            - ``on_gens`` -- a function that determines the image of the
              generators of ``self`` under the crystal morphism
            - ``codomain`` -- (default: ``self``) the codomain of the morphism
            - ``cartan_type`` -- (optional) the Cartan type of the morphism;
              the default is the Cartan type of ``self``
            - ``index_set`` -- (optional) the index set of the morphism;
              the default is the index set of the Cartan type
            - ``generators`` -- (optional) the generators to define the
              morphism; the default is the generators of ``self``
            - ``automorphism`` -- (optional) the automorphism to perform the
              twisting
            - ``virtualization`` -- a dictionary whose keys are in the index
              set of the domain and whose values are lists of entries in the
              index set of the codomain
            - ``scaling_factors`` -- a dictionary whose keys are in the index
              set of the domain and whose values are scaling factors for the
              weight, `\varepsilon` and `\varphi`.
            - ``category`` -- (optional) the category for the crystal morphism;
              the default is the category of :class:`Crystals`.
            - ``check`` -- (default: ``True``) check if the crystal morphism
              is valid

            .. SEEALSO::

                For more examples, see
                :class:`sage.categories.crystals.CrystalHomset`.

            EXAMPLES:

            We construct the natural embedding of a crystal using tableaux
            into the tensor product of single boxes via the reading word::

                sage: B = crystals.Tableaux(['A',2], shape=[2,1])
                sage: F = crystals.Tableaux(['A',2], shape=[1])
                sage: T = crystals.TensorProduct(F, F, F)
                sage: mg = T.highest_weight_vectors()[2]; mg
                [[[2]], [[1]], [[1]]]
                sage: psi = B.crystal_morphism([mg], codomain=T); psi
                ['A', 2] Crystal morphism:
                  From: The crystal of tableaux of type ['A', 2] and shape(s) [[2, 1]]
                  To:   Full tensor product of the crystals
                         [The crystal of tableaux of type ['A', 2] and shape(s) [[1]],
                          The crystal of tableaux of type ['A', 2] and shape(s) [[1]],
                          The crystal of tableaux of type ['A', 2] and shape(s) [[1]]]
                  Defn: [2, 1, 1] |--> [[[2]], [[1]], [[1]]]
                sage: b = B.module_generators[0]
                sage: b.pp()
                  1  1
                  2
                sage: psi(b)
                [[[2]], [[1]], [[1]]]
                sage: psi(b.f(2))
                [[[3]], [[1]], [[1]]]
                sage: psi(b.f_string([2,1,1]))
                [[[3]], [[2]], [[2]]]
                sage: lw = b.to_lowest_weight()[0]
                sage: lw.pp()
                  2  3
                  3
                sage: psi(lw)
                [[[3]], [[2]], [[3]]]
                sage: psi(lw) == mg.to_lowest_weight()[0]
                True

            We now take the other isomorphic factor in the tensor product::

                sage: mg = T.highest_weight_vectors()[1]; mg
                [[[1]], [[2]], [[1]]]
                sage: psi = B.crystal_morphism([mg], codomain=T)
                sage: psi(lw)
                [[[3]], [[3]], [[2]]]

            We construct a crystal morphism as classical crystals using a
            Kirillov-Reshetikhin crystal::

                sage: B = crystals.Tableaux(['D', 4], shape=[1,1])
                sage: K = crystals.KirillovReshetikhin(['D',4,1], 2,2)
                sage: K.module_generators
                [[], [[1], [2]], [[1, 1], [2, 2]]]
                sage: v = K.module_generators[1]
                sage: psi = B.crystal_morphism([v], codomain=K, cartan_type=['D',4], category=FiniteCrystals())
                sage: psi
                ['D', 4] Crystal morphism:
                  From: The crystal of tableaux of type ['D', 4] and shape(s) [[1, 1]]
                  To:   Kirillov-Reshetikhin crystal of type ['D', 4, 1] with (r,s)=(2,2)
                  Defn: [2, 1] |--> [[1], [2]]
                sage: b = B.module_generators[0]
                sage: psi(b)
                [[1], [2]]
                sage: psi(b.to_lowest_weight()[0])
                [[-2], [-1]]
            """
            # Determine the codomain
            if codomain is None:
                if hasattr(on_gens, 'codomain'):
                    codomain = on_gens.codomain()
                elif isinstance(on_gens, (list, tuple)):
                    if on_gens:
                        codomain = on_gens[0].parent()
                elif isinstance(on_gens, dict):
                    if on_gens:
                        codomain = on_gens.values()[0].parent()
                else:
                    for x in self.module_generators:
                        y = on_gens(x)
                        if y is not None:
                            codomain = y.parent()
                            break
            if codomain is None:
                codomain = self
            elif codomain not in Crystals():
                raise ValueError("the codomain must be a crystal")

            if category is None:
                category = self.category()
            homset = Hom(self, codomain, category=category)
            return homset(on_gens, cartan_type, index_set, generators,
                          automorphism, virtualization, scaling_factors, check)

        def digraph(self, subset=None, index_set=None):
            """
            Returns the DiGraph associated to ``self``.

            INPUT:

            - ``subset`` -- (Optional) A subset of vertices for
              which the digraph should be constructed

            - ``index_set`` -- (Optional) The index set to draw arrows

            EXAMPLES::

                sage: C = Crystals().example(5)
                sage: C.digraph()
                Digraph on 6 vertices

            The edges of the crystal graph are by default colored using blue for edge 1, red for edge 2,
            and green for edge 3::

                sage: C = Crystals().example(3)
                sage: G = C.digraph()
                sage: view(G, pdflatex=True, tightpage=True)  #optional - dot2tex graphviz

            One may also overwrite the colors::

                sage: C = Crystals().example(3)
                sage: G = C.digraph()
                sage: G.set_latex_options(color_by_label = {1:"red", 2:"purple", 3:"blue"})
                sage: view(G, pdflatex=True, tightpage=True)  #optional - dot2tex graphviz

            Or one may add colors to yet unspecified edges::

                sage: C = Crystals().example(4)
                sage: G = C.digraph()
                sage: C.cartan_type()._index_set_coloring[4]="purple"
                sage: view(G, pdflatex=True, tightpage=True)  #optional - dot2tex graphviz

            Here is an example of how to take the top part up to a given depth of an infinite dimensional
            crystal::

                sage: C = CartanType(['C',2,1])
                sage: La = C.root_system().weight_lattice().fundamental_weights()
                sage: T = crystals.HighestWeight(La[0])
                sage: S = T.subcrystal(max_depth=3)
                sage: G = T.digraph(subset=S); G
                Digraph on 5 vertices
                sage: sorted(G.vertices(), key=str)
                [(-Lambda[0] + 2*Lambda[1] - delta,),
                 (1/2*Lambda[0] + Lambda[1] - Lambda[2] - 1/2*delta, -1/2*Lambda[0] + Lambda[1] - 1/2*delta),
                 (1/2*Lambda[0] - Lambda[1] + Lambda[2] - 1/2*delta, -1/2*Lambda[0] + Lambda[1] - 1/2*delta),
                 (Lambda[0] - 2*Lambda[1] + 2*Lambda[2] - delta,),
                 (Lambda[0],)]

            Here is a way to construct a picture of a Demazure crystal using
            the ``subset`` option::

                sage: B = crystals.Tableaux(['A',2], shape=[2,1])
                sage: C = CombinatorialFreeModule(QQ,B)
                sage: t = B.highest_weight_vector()
                sage: b = C(t)
                sage: D = B.demazure_operator(b,[2,1]); D
                B[[[1, 1], [2]]] + B[[[1, 2], [2]]] + B[[[1, 3], [2]]] + B[[[1, 1], [3]]] + B[[[1, 3], [3]]]
                sage: G = B.digraph(subset=D.support())
                sage: G.vertices()
                [[[1, 1], [2]], [[1, 2], [2]], [[1, 3], [2]], [[1, 1], [3]], [[1, 3], [3]]]
                sage: view(G, pdflatex=True, tightpage=True)  #optional - dot2tex graphviz

            We can also choose to display particular arrows using the
            ``index_set`` option::

                sage: C = crystals.KirillovReshetikhin(['D',4,1], 2, 1)
                sage: G = C.digraph(index_set=[1,3])
                sage: len(G.edges())
                20
                sage: view(G, pdflatex=True, tightpage=True)  #optional - dot2tex graphviz

            TODO: add more tests
            """
            from sage.graphs.all import DiGraph
            from sage.categories.highest_weight_crystals import HighestWeightCrystals
            d = {}
            if self in HighestWeightCrystals:
                f = lambda u_v_label: ({})
            else:
                f = lambda u_v_label: ({"backward": u_v_label[2] == 0})

            # Parse optional arguments
            if subset is None:
                subset = self
            if index_set is None:
                index_set = self.index_set()

            for x in subset:
                d[x] = {}
                for i in index_set:
                    child = x.f(i)
                    if child is None or child not in subset:
                        continue
                    d[x][child]=i
            G = DiGraph(d)
            if have_dot2tex():
                G.set_latex_options(format="dot2tex",
                                    edge_labels = True,
                                    color_by_label = self.cartan_type()._index_set_coloring,
                                    edge_options = f)
            return G

        def latex_file(self, filename):
            r"""
            Exports a file, suitable for pdflatex, to 'filename'. This requires
            a proper installation of ``dot2tex`` in sage-python. For more
            information see the documentation for ``self.latex()``.

            EXAMPLES::

                sage: C = crystals.Letters(['A', 5])
                sage: C.latex_file('/tmp/test.tex') #optional - dot2tex
            """
            header = r"""\documentclass{article}
            \usepackage[x11names, rgb]{xcolor}
            \usepackage[utf8]{inputenc}
            \usepackage{tikz}
            \usetikzlibrary{snakes,arrows,shapes}
            \usepackage{amsmath}
            \usepackage[active,tightpage]{preview}
            \newenvironment{bla}{}{}
            \PreviewEnvironment{bla}

            \begin{document}
            \begin{bla}"""

            footer = r"""\end{bla}
            \end{document}"""

            f = open(filename, 'w')
            f.write(header + self.latex() + footer)
            f.close()

        def _latex_(self, **options):
            r"""
            Returns the crystal graph as a latex string. This can be exported
            to a file with self.latex_file('filename').

            EXAMPLES::

                sage: T = crystals.Tableaux(['A',2],shape=[1])
                sage: T._latex_()   #optional - dot2tex
                '...tikzpicture...'
                sage: view(T, pdflatex = True, tightpage = True) #optional - dot2tex graphviz

            One can for example also color the edges using the following options::

                sage: T = crystals.Tableaux(['A',2],shape=[1])
                sage: T._latex_(color_by_label = {0:"black", 1:"red", 2:"blue"})   #optional - dot2tex graphviz
                '...tikzpicture...'
            """
            if not have_dot2tex():
                print "dot2tex not available.  Install after running \'sage -sh\'"
                return
            G=self.digraph()
            G.set_latex_options(**options)
            return G._latex_()

        latex = _latex_

        def metapost(self, filename, thicklines=False, labels=True, scaling_factor=1.0, tallness=1.0):
            r"""
            Use C.metapost("filename.mp",[options]), where options can be:

            thicklines = True (for thicker edges) labels = False (to suppress
            labeling of the vertices) scaling_factor=value, where value is a
            floating point number, 1.0 by default. Increasing or decreasing the
            scaling factor changes the size of the image. tallness=1.0.
            Increasing makes the image taller without increasing the width.

            Root operators e(1) or f(1) move along red lines, e(2) or f(2)
            along green. The highest weight is in the lower left. Vertices with
            the same weight are kept close together. The concise labels on the
            nodes are strings introduced by Berenstein and Zelevinsky and
            Littelmann; see Littelmann's paper Cones, Crystals, Patterns,
            sections 5 and 6.

            For Cartan types B2 or C2, the pattern has the form

            a2 a3 a4 a1

            where c\*a2 = a3 = 2\*a4 =0 and a1=0, with c=2 for B2, c=1 for C2.
            Applying e(2) a1 times, e(1) a2 times, e(2) a3 times, e(1) a4 times
            returns to the highest weight. (Observe that Littelmann writes the
            roots in opposite of the usual order, so our e(1) is his e(2) for
            these Cartan types.) For type A2, the pattern has the form

            a3 a2 a1

            where applying e(1) a1 times, e(2) a2 times then e(3) a1 times
            returns to the highest weight. These data determine the vertex and
            may be translated into a Gelfand-Tsetlin pattern or tableau.

            EXAMPLES::

                sage: C = crystals.Letters(['A', 2])
                sage: C.metapost('/tmp/test.mp') #optional

            ::

                sage: C = crystals.Letters(['A', 5])
                sage: C.metapost('/tmp/test.mp')
                Traceback (most recent call last):
                ...
                NotImplementedError
            """
            # FIXME: those tests are not robust
            # Should use instead self.cartan_type() == CartanType(['B',2])
            if self.cartan_type()[0] == 'B' and self.cartan_type()[1] == 2:
                word = [2,1,2,1]
            elif self.cartan_type()[0] == 'C' and self.cartan_type()[1] == 2:
                word = [2,1,2,1]
            elif self.cartan_type()[0] == 'A' and self.cartan_type()[1] == 2:
                word = [1,2,1]
            else:
                raise NotImplementedError
            size = self.cardinality()
            string_data = []
            for i in range(size):
                turtle = self.list()[i]
                string_datum = []
                for j in word:
                    turtlewalk = 0
                    while turtle.e(j) is not None:
                        turtle = turtle.e(j)
                        turtlewalk += 1
                    string_datum.append(turtlewalk)
                string_data.append(string_datum)

            if self.cartan_type()[0] == 'A':
                if labels:
                    c0 = int(55*scaling_factor)
                    c1 = int(-25*scaling_factor)
                    c2 = int(45*tallness*scaling_factor)
                    c3 = int(-12*scaling_factor)
                    c4 = int(-12*scaling_factor)
                else:
                    c0 = int(45*scaling_factor)
                    c1 = int(-20*scaling_factor)
                    c2 = int(35*tallness*scaling_factor)
                    c3 = int(12*scaling_factor)
                    c4 = int(-12*scaling_factor)
                outstring = "verbatimtex\n\\magnification=600\netex\n\nbeginfig(-1);\nsx:=35; sy:=30;\n\nz1000=(%d,0);\nz1001=(%d,%d);\nz1002=(%d,%d);\nz2001=(-3,3);\nz2002=(3,3);\nz2003=(0,-3);\nz2004=(7,0);\nz2005=(0,7);\nz2006=(-7,0);\nz2007=(0,7);\n\n"%(c0,c1,c2,c3,c4)
            else:
                if labels:
                    outstring = "verbatimtex\n\\magnification=600\netex\n\nbeginfig(-1);\n\nsx := %d;\nsy=%d;\n\nz1000=(2*sx,0);\nz1001=(-sx,sy);\nz1002=(-16,-10);\n\nz2001=(0,-3);\nz2002=(-5,3);\nz2003=(0,3);\nz2004=(5,3);\nz2005=(10,1);\nz2006=(0,10);\nz2007=(-10,1);\nz2008=(0,-8);\n\n"%(int(scaling_factor*40),int(tallness*scaling_factor*40))
                else:
                    outstring = "beginfig(-1);\n\nsx := %d;\nsy := %d;\n\nz1000=(2*sx,0);\nz1001=(-sx,sy);\nz1002=(-5,-5);\n\nz1003=(10,10);\n\n"%(int(scaling_factor*35),int(tallness*scaling_factor*35))
            for i in range(size):
                if self.cartan_type()[0] == 'A':
                    [a1,a2,a3] = string_data[i]
                else:
                    [a1,a2,a3,a4] = string_data[i]
                shift = 0
                for j in range(i):
                    if self.cartan_type()[0] == 'A':
                        [b1,b2,b3] = string_data[j]
                        if b1+b3 == a1+a3 and b2 == a2:
                            shift += 1
                    else:
                        [b1,b2,b3,b4] = string_data[j]
                        if b1+b3 == a1+a3 and b2+b4 == a2+a4:
                            shift += 1
                if self.cartan_type()[0] == 'A':
                    outstring = outstring +"z%d=%d*z1000+%d*z1001+%d*z1002;\n"%(i,a1+a3,a2,shift)
                else:
                    outstring = outstring +"z%d=%d*z1000+%d*z1001+%d*z1002;\n"%(i,a1+a3,a2+a4,shift)
            outstring = outstring + "\n"
            if thicklines:
                outstring = outstring +"pickup pencircle scaled 2\n\n"
            for i in range(size):
                for j in range(1,3):
                    dest = self.list()[i].f(j)
                    if dest is not None:
                        dest = self.list().index(dest)
                        if j == 1:
                            col = "red;"
                        else:
                            col = "green;  "
                        if self.cartan_type()[0] == 'A':
                            [a1,a2,a3] = string_data[i] # included to facilitate hand editing of the .mp file
                            outstring = outstring+"draw z%d--z%d withcolor %s   %% %d %d %d\n"%(i,dest,col,a1,a2,a3)
                        else:
                            [a1,a2,a3,a4] = string_data[i]
                            outstring = outstring+"draw z%d--z%d withcolor %s   %% %d %d %d %d\n"%(i,dest,col,a1,a2,a3,a4)
            outstring += "\npickup pencircle scaled 3;\n\n"
            for i in range(self.cardinality()):
                if labels:
                    if self.cartan_type()[0] == 'A':
                        outstring = outstring+"pickup pencircle scaled 15;\nfill z%d+z2004..z%d+z2006..z%d+z2006..z%d+z2007..cycle withcolor white;\nlabel(btex %d etex, z%d+z2001);\nlabel(btex %d etex, z%d+z2002);\nlabel(btex %d etex, z%d+z2003);\npickup pencircle scaled .5;\ndraw z%d+z2004..z%d+z2006..z%d+z2006..z%d+z2007..cycle;\n"%(i,i,i,i,string_data[i][2],i,string_data[i][1],i,string_data[i][0],i,i,i,i,i)
                    else:
                        outstring = outstring+"%%%d %d %d %d\npickup pencircle scaled 1;\nfill z%d+z2005..z%d+z2006..z%d+z2007..z%d+z2008..cycle withcolor white;\nlabel(btex %d etex, z%d+z2001);\nlabel(btex %d etex, z%d+z2002);\nlabel(btex %d etex, z%d+z2003);\nlabel(btex %d etex, z%d+z2004);\npickup pencircle scaled .5;\ndraw z%d+z2005..z%d+z2006..z%d+z2007..z%d+z2008..cycle;\n\n"%(string_data[i][0],string_data[i][1],string_data[i][2],string_data[i][3],i,i,i,i,string_data[i][0],i,string_data[i][1],i,string_data[i][2],i,string_data[i][3],i,i,i,i,i)
                else:
                    outstring += "drawdot z%d;\n"%i
            outstring += "\nendfig;\n\nend;\n\n"

            f = open(filename, 'w')
            f.write(outstring)
            f.close()

        def dot_tex(self):
            r"""
            Return a dot_tex string representation of ``self``.

            EXAMPLES::

                sage: C = crystals.Letters(['A',2])
                sage: C.dot_tex()
                'digraph G { \n  node [ shape=plaintext ];\n  N_0 [ label = " ", texlbl = "$1$" ];\n  N_1 [ label = " ", texlbl = "$2$" ];\n  N_2 [ label = " ", texlbl = "$3$" ];\n  N_0 -> N_1 [ label = " ", texlbl = "1" ];\n  N_1 -> N_2 [ label = " ", texlbl = "2" ];\n}'
            """
            import re
            from sage.combinat import ranker
            rank = ranker.from_list(self.list())[0]
            vertex_key = lambda x: "N_"+str(rank(x))

            # To do: check the regular expression
            # Removing %-style comments, newlines, quotes
            # This should probably be moved to sage.misc.latex
            quoted_latex = lambda x: re.sub("\"|\r|(%[^\n]*)?\n","", latex(x))

            result = "digraph G { \n  node [ shape=plaintext ];\n"

            for x in self:
                result += "  " + vertex_key(x) + " [ label = \" \", texlbl = \"$"+quoted_latex(x)+"$\" ];\n"
            for x in self:
                for i in self.index_set():
                    child = x.f(i)
                    if child is None:
                        continue
    #                result += "  " + vertex_key(x) + " -> "+vertex_key(child)+ " [ label = \" \", texlbl = \""+quoted_latex(i)+"\" ];\n"
                    if i == 0:
                        option = "dir = back, "
                        (source, target) = (child, x)
                    else:
                        option = ""
                        (source, target) = (x, child)
                    result += "  " + vertex_key(source) + " -> "+vertex_key(target)+ " [ "+option+"label = \" \", texlbl = \""+quoted_latex(i)+"\" ];\n"
            result+="}"
            return result

        def plot(self, **options):
            """
            Return the plot of ``self`` as a directed graph.

            EXAMPLES::

                sage: C = crystals.Letters(['A', 5])
                sage: print(C.plot())
                Graphics object consisting of 17 graphics primitives
            """
            return self.digraph().plot(edge_labels=True,vertex_size=0,**options)

        def plot3d(self, **options):
            """
            Return the 3-dimensional plot of ``self`` as a directed graph.

            EXAMPLES::

                sage: C = crystals.KirillovReshetikhin(['A',3,1],2,1)
                sage: print(C.plot3d())
                Graphics3d Object
            """
            G = self.digraph(**options)
            return G.plot3d()

        def tensor(self, *crystals, **options):
            """
            Return the tensor product of ``self`` with the crystals ``B``.

            EXAMPLES::

                sage: C = crystals.Letters(['A', 3])
                sage: B = crystals.infinity.Tableaux(['A', 3])
                sage: T = C.tensor(C, B); T
                Full tensor product of the crystals
                 [The crystal of letters for type ['A', 3],
                  The crystal of letters for type ['A', 3],
                  The infinity crystal of tableaux of type ['A', 3]]
                sage: tensor([C, C, B]) is T
                True

                sage: C = crystals.Letters(['A',2])
                sage: T = C.tensor(C, C, generators=[[C(2),C(1),C(1)],[C(1),C(2),C(1)]]); T
                The tensor product of the crystals
                 [The crystal of letters for type ['A', 2],
                  The crystal of letters for type ['A', 2],
                  The crystal of letters for type ['A', 2]]
                sage: T.module_generators
                [[2, 1, 1], [1, 2, 1]]
            """
            from sage.combinat.crystals.tensor_product import TensorProductOfCrystals
            return TensorProductOfCrystals(self, *crystals, **options)

        def direct_sum(self, X):
            """
            Return the direct sum of ``self`` with ``X``.

            EXAMPLES::

                sage: B = crystals.Tableaux(['A',2], shape=[2,1])
                sage: C = crystals.Letters(['A',2])
                sage: B.direct_sum(C)
                Direct sum of the crystals Family
                (The crystal of tableaux of type ['A', 2] and shape(s) [[2, 1]],
                 The crystal of letters for type ['A', 2])

            As a shorthand, we can use ``+``::

                sage: B + C
                Direct sum of the crystals Family
                (The crystal of tableaux of type ['A', 2] and shape(s) [[2, 1]],
                 The crystal of letters for type ['A', 2])
            """
            if X not in Crystals():
                raise ValueError("{} is not a crystal".format(X))
            from sage.combinat.crystals.direct_sum import DirectSumOfCrystals
            return DirectSumOfCrystals([self, X])

        __add__ = direct_sum

        @abstract_method(optional=True)
        def connected_components_generators(self):
            """
            Return a tuple of generators for each of the connected components
            of ``self``.

            EXAMPLES::

                sage: B = crystals.Tableaux(['A',2], shape=[2,1])
                sage: C = crystals.Letters(['A',2])
                sage: T = crystals.TensorProduct(B,C)
                sage: T.connected_components_generators()
                ([[[1, 1], [2]], 1], [[[1, 2], [2]], 1], [[[1, 2], [3]], 1])
            """

        def connected_components(self):
            """
            Return the connected components of ``self`` as subcrystals.

            EXAMPLES::

                sage: B = crystals.Tableaux(['A',2], shape=[2,1])
                sage: C = crystals.Letters(['A',2])
                sage: T = crystals.TensorProduct(B,C)
                sage: T.connected_components()
                [Subcrystal of Full tensor product of the crystals
                 [The crystal of tableaux of type ['A', 2] and shape(s) [[2, 1]],
                  The crystal of letters for type ['A', 2]],
                 Subcrystal of Full tensor product of the crystals
                 [The crystal of tableaux of type ['A', 2] and shape(s) [[2, 1]],
                  The crystal of letters for type ['A', 2]],
                 Subcrystal of Full tensor product of the crystals
                 [The crystal of tableaux of type ['A', 2] and shape(s) [[2, 1]],
                  The crystal of letters for type ['A', 2]]]
            """
            return map(lambda mg: self.subcrystal(generators=[mg]),
                       self.connected_components_generators())

        def number_of_connected_components(self):
            """
            Return the number of connected components of ``self``.

            EXAMPLES::

                sage: B = crystals.Tableaux(['A',2], shape=[2,1])
                sage: C = crystals.Letters(['A',2])
                sage: T = crystals.TensorProduct(B,C)
                sage: T.number_of_connected_components()
                3
            """
            return len(self.connected_components_generators())

        def is_connected(self):
            """
            Return ``True`` if ``self`` is a connected crystal.

            EXAMPLES::

                sage: B = crystals.Tableaux(['A',2], shape=[2,1])
                sage: C = crystals.Letters(['A',2])
                sage: T = crystals.TensorProduct(B,C)
                sage: B.is_connected()
                True
                sage: T.is_connected()
                False
            """
            return self.number_of_connected_components() == 1

    class ElementMethods:

        @cached_method
        def index_set(self):
            """
            EXAMPLES::

                sage: C = crystals.Letters(['A',5])
                sage: C(1).index_set()
                (1, 2, 3, 4, 5)
            """
            return self.parent().index_set()

        def cartan_type(self):
            """
            Returns the Cartan type associated to ``self``

            EXAMPLES::

                sage: C = crystals.Letters(['A', 5])
                sage: C(1).cartan_type()
                ['A', 5]
            """
            return self.parent().cartan_type()

        @abstract_method
        def e(self, i):
            r"""
            Returns `e_i` of ``self`` if it exists or ``None`` otherwise.

            This method should be implemented by the element class of
            the crystal.

            EXAMPLES::

                sage: C = Crystals().example(5)
                sage: x = C[2]; x
                3
                sage: x.e(1), x.e(2), x.e(3)
                (None, 2, None)
            """

        @abstract_method
        def f(self, i):
            r"""
            Returns `f_i` of ``self`` if it exists or ``None`` otherwise.

            This method should be implemented by the element class of
            the crystal.

            EXAMPLES::

                sage: C = Crystals().example(5)
                sage: x = C[1]; x
                2
                sage: x.f(1), x.f(2), x.f(3)
                (None, 3, None)
            """

        @abstract_method
        def epsilon(self, i):
            r"""
            EXAMPLES::

                sage: C = crystals.Letters(['A',5])
                sage: C(1).epsilon(1)
                0
                sage: C(2).epsilon(1)
                1
            """

        @abstract_method
        def phi(self, i):
            r"""
            EXAMPLES::

                sage: C = crystals.Letters(['A',5])
                sage: C(1).phi(1)
                1
                sage: C(2).phi(1)
                0
            """

        @abstract_method
        def weight(self):
            r"""
            Return the weight of this crystal element.

            This method should be implemented by the element class of
            the crystal.

            EXAMPLES::

                sage: C = crystals.Letters(['A',5])
                sage: C(1).weight()
                (1, 0, 0, 0, 0, 0)
            """

        def phi_minus_epsilon(self, i):
            """
            Return `\phi_i - \epsilon_i` of ``self``.

            There are sometimes better implementations using the weight
            for this. It is used for reflections along a string.

            EXAMPLES::

                sage: C = crystals.Letters(['A',5])
                sage: C(1).phi_minus_epsilon(1)
                1
            """
            return self.phi(i) - self.epsilon(i)

        def Epsilon(self):
            """
            EXAMPLES::

                sage: C = crystals.Letters(['A',5])
                sage: C(0).Epsilon()
                (0, 0, 0, 0, 0, 0)
                sage: C(1).Epsilon()
                (0, 0, 0, 0, 0, 0)
                sage: C(2).Epsilon()
                (1, 0, 0, 0, 0, 0)
            """
            Lambda = self.parent().Lambda()
            return sum(self.epsilon(i) * Lambda[i] for i in self.index_set())

        def Phi(self):
            """
            EXAMPLES::

                sage: C = crystals.Letters(['A',5])
                sage: C(0).Phi()
                (0, 0, 0, 0, 0, 0)
                sage: C(1).Phi()
                (1, 0, 0, 0, 0, 0)
                sage: C(2).Phi()
                (1, 1, 0, 0, 0, 0)
            """
            Lambda = self.parent().Lambda()
            return sum(self.phi(i) * Lambda[i] for i in self.index_set())

        def f_string(self, list):
            r"""
            Applies `f_{i_r} \cdots f_{i_1}` to self for ``list`` as
            `[i_1, ..., i_r]`

            EXAMPLES::

                sage: C = crystals.Letters(['A',3])
                sage: b = C(1)
                sage: b.f_string([1,2])
                3
                sage: b.f_string([2,1])
            """
            b = self
            for i in list:
                b = b.f(i)
                if b is None:
                    return None
            return b

        def e_string(self, list):
            r"""
            Applies `e_{i_r} \cdots e_{i_1}` to self for ``list`` as
            `[i_1, ..., i_r]`

            EXAMPLES::

                sage: C = crystals.Letters(['A',3])
                sage: b = C(3)
                sage: b.e_string([2,1])
                1
                sage: b.e_string([1,2])
            """
            b = self
            for i in list:
                b = b.e(i)
                if b is None:
                    return None
            return b

        def s(self, i):
            r"""
            Return the reflection of ``self`` along its `i`-string.

            EXAMPLES::

                sage: C = crystals.Tableaux(['A',2], shape=[2,1])
                sage: b = C(rows=[[1,1],[3]])
                sage: b.s(1)
                [[2, 2], [3]]
                sage: b = C(rows=[[1,2],[3]])
                sage: b.s(2)
                [[1, 2], [3]]
                sage: T = crystals.Tableaux(['A',2],shape=[4])
                sage: t = T(rows=[[1,2,2,2]])
                sage: t.s(1)
                [[1, 1, 1, 2]]
            """
            d = self.phi_minus_epsilon(i)
            b = self
            if d > 0:
                for j in range(d):
                    b = b.f(i)
            else:
                for j in range(-d):
                    b = b.e(i)
            return b

        def is_highest_weight(self, index_set = None):
            r"""
            Returns ``True`` if ``self`` is a highest weight.
            Specifying the option ``index_set`` to be a subset `I` of the
            index set of the underlying crystal, finds all highest
            weight vectors for arrows in `I`.

            EXAMPLES::

                sage: C = crystals.Letters(['A',5])
                sage: C(1).is_highest_weight()
                True
                sage: C(2).is_highest_weight()
                False
                sage: C(2).is_highest_weight(index_set = [2,3,4,5])
                True
            """
            if index_set is None:
                index_set = self.index_set()
            return all(self.e(i) is None for i in index_set)

        def is_lowest_weight(self, index_set = None):
            r"""
            Returns ``True`` if ``self`` is a lowest weight.
            Specifying the option ``index_set`` to be a subset `I` of the
            index set of the underlying crystal, finds all lowest
            weight vectors for arrows in `I`.

            EXAMPLES::

                sage: C = crystals.Letters(['A',5])
                sage: C(1).is_lowest_weight()
                False
                sage: C(6).is_lowest_weight()
                True
                sage: C(4).is_lowest_weight(index_set = [1,3])
                True
            """
            if index_set is None:
                index_set = self.index_set()
            return all(self.f(i) is None for i in index_set)

        def to_highest_weight(self, index_set = None):
            r"""
            Yields the highest weight element `u` and a list `[i_1,...,i_k]`
            such that `self = f_{i_1} ... f_{i_k} u`, where `i_1,...,i_k` are
            elements in `index_set`. By default the index set is assumed to be
            the full index set of self.

            EXAMPLES::

                sage: T = crystals.Tableaux(['A',3], shape = [1])
                sage: t = T(rows = [[3]])
                sage: t.to_highest_weight()
                [[[1]], [2, 1]]
                sage: T = crystals.Tableaux(['A',3], shape = [2,1])
                sage: t = T(rows = [[1,2],[4]])
                sage: t.to_highest_weight()
                [[[1, 1], [2]], [1, 3, 2]]
                sage: t.to_highest_weight(index_set = [3])
                [[[1, 2], [3]], [3]]
                sage: K = crystals.KirillovReshetikhin(['A',3,1],2,1)
                sage: t = K(rows=[[2],[3]]); t.to_highest_weight(index_set=[1])
                [[[1], [3]], [1]]
                sage: t.to_highest_weight()
                Traceback (most recent call last):
                ...
                ValueError: This is not a highest weight crystals!
            """
            from sage.categories.highest_weight_crystals import HighestWeightCrystals
            if index_set is None:
                if HighestWeightCrystals() not in self.parent().categories():
                    raise ValueError("This is not a highest weight crystals!")
                index_set = self.index_set()
            for i in index_set:
                next = self.e(i)
                if next is not None:
                    hw = next.to_highest_weight(index_set = index_set)
                    return [hw[0], [i] + hw[1]]
            return [self, []]

        def to_lowest_weight(self, index_set = None):
            r"""
            Yields the lowest weight element `u` and a list `[i_1,...,i_k]`
            such that `self = e_{i_1} ... e_{i_k} u`, where `i_1,...,i_k` are
            elements in `index_set`. By default the index set is assumed to be
            the full index set of self.

            EXAMPLES::

                sage: T = crystals.Tableaux(['A',3], shape = [1])
                sage: t = T(rows = [[3]])
                sage: t.to_lowest_weight()
                [[[4]], [3]]
                sage: T = crystals.Tableaux(['A',3], shape = [2,1])
                sage: t = T(rows = [[1,2],[4]])
                sage: t.to_lowest_weight()
                [[[3, 4], [4]], [1, 2, 2, 3]]
                sage: t.to_lowest_weight(index_set = [3])
                [[[1, 2], [4]], []]
                sage: K = crystals.KirillovReshetikhin(['A',3,1],2,1)
                sage: t = K.module_generator(); t
                [[1], [2]]
                sage: t.to_lowest_weight(index_set=[1,2,3])
                [[[3], [4]], [2, 1, 3, 2]]
                sage: t.to_lowest_weight()
                Traceback (most recent call last):
                ...
                ValueError: This is not a highest weight crystals!
            """
            from sage.categories.highest_weight_crystals import HighestWeightCrystals
            if index_set is None:
                if HighestWeightCrystals() not in self.parent().categories():
                    raise ValueError, "This is not a highest weight crystals!"
                index_set = self.index_set()
            for i in index_set:
                next = self.f(i)
                if next is not None:
                    lw = next.to_lowest_weight(index_set = index_set)
                    return [lw[0], [i] + lw[1]]
            return [self, []]

        def subcrystal(self, index_set=None, max_depth=float("inf"), direction="both",
                       contained=None, cartan_type=None, category=None):
            r"""
            Construct the subcrystal generated by ``self`` using `e_i` and/or
            `f_i` for all `i` in ``index_set``.

            INPUT:

            - ``index_set`` -- (default: ``None``) the index set; if ``None``
              then use the index set of the crystal

            - ``max_depth`` -- (default: infinity) the maximum depth to build

            - ``direction`` -- (default: ``'both'``) the direction to build
              the subcrystal; it can be one of the following:

              - ``'both'`` - using both `e_i` and `f_i`
              - ``'upper'`` - using `e_i`
              - ``'lower'`` - using `f_i`

            - ``contained`` -- (optional) a set (or function) defining the
              containment in the subcrystal

            - ``cartan_type`` -- (optional) specify the Cartan type of the
              subcrystal

            - ``category`` -- (optional) specify the category of the subcrystal

            .. SEEALSO::

                - :meth:`Crystals.ParentMethods.subcrystal()`

            EXAMPLES::

                sage: C = crystals.KirillovReshetikhin(['A',3,1], 1, 2)
                sage: elt = C(1,4)
                sage: list(elt.subcrystal(index_set=[1,3]))
                [[[1, 4]], [[2, 4]], [[2, 3]], [[1, 3]]]
                sage: list(elt.subcrystal(index_set=[1,3], max_depth=1))
                [[[1, 4]], [[2, 4]], [[1, 3]]]
                sage: list(elt.subcrystal(index_set=[1,3], direction='upper'))
                [[[1, 4]], [[1, 3]]]
                sage: list(elt.subcrystal(index_set=[1,3], direction='lower'))
                [[[1, 4]], [[2, 4]]]
            """
            return self.parent().subcrystal(generators=[self], index_set=index_set,
                                            max_depth=max_depth, direction=direction)

    class SubcategoryMethods:
        """
        Methods for all subcategories.
        """
        def TensorProducts(self):
            r"""
            Return the full subcategory of objects of ``self`` constructed
            as tensor products.

            .. SEEALSO::

                - :class:`.tensor.TensorProductsCategory`
                - :class:`~.covariant_functorial_construction.RegressiveCovariantFunctorialConstruction`.

            EXAMPLES::

                sage: HighestWeightCrystals().TensorProducts()
                Category of tensor products of highest weight crystals
            """
            return TensorProductsCategory.category_of(self)

    class TensorProducts(TensorProductsCategory):
        """
        The category of crystals constructed by tensor product of crystals.
        """
        @cached_method
        def extra_super_categories(self):
            """
            EXAMPLES::

                sage: Crystals().TensorProducts().extra_super_categories()
                [Category of crystals]
            """
            return [self.base_category()]

    Finite = LazyImport('sage.categories.finite_crystals', 'FiniteCrystals')

###############################################################################
## Morphisms

class CrystalMorphism(Morphism):
    r"""
    A crystal morphism.

    INPUT:

    - ``parent`` -- a homset
    - ``cartan_type`` -- (optional) the Cartan type of the morphism; the
      default is the Cartan type of the domain
    - ``index_set`` -- (optional) the index set of the morphism; the
      default is the index set of the Cartan type
    """
    def __init__(self, parent, cartan_type=None, index_set=None):
        """
        Initialize ``self``.

        TESTS::

            sage: B = crystals.Tableaux(['A',2], shape=[2,1])
            sage: H = Hom(B, B)
            sage: psi = H.an_element()
        """
        if cartan_type is None:
            cartan_type = parent.domain().cartan_type()
        from sage.combinat.root_system.cartan_type import CartanType
        self._cartan_type = CartanType(cartan_type)

        if index_set is None:
            index_set = self._cartan_type.index_set()
        self._index_set = tuple(index_set)
        Morphism.__init__(self, parent)

    def _repr_type(self):
        """
        Used internally in printing this morphism.

        TESTS::

            sage: B = crystals.Tableaux(['A',2], shape=[2,1])
            sage: H = Hom(B, B)
            sage: psi = H.an_element()
            sage: psi._repr_type()
            "['A', 2] Crystal"

            sage: psi = H(lambda x: None, index_set=[1])
            sage: psi._repr_type()
            "['A', 2] [1]-Crystal"
        """
        if self._index_set != self._cartan_type.index_set():
            return "{} {}-Crystal".format(self._cartan_type, list(self._index_set))
        return "{} Crystal".format(self._cartan_type)

    def index_set(self):
        """
        Return the index set of ``self``.

        EXAMPLES::

            sage: B = crystals.Tableaux(['A',2], shape=[2,1])
            sage: H = Hom(B, B)
            sage: psi = H.an_element()
            sage: psi.index_set()
            (1, 2)
            sage: psi = H(lambda x: None, index_set=[1])
            sage: psi.index_set()
            (1,)
        """
        return self._index_set

    def cartan_type(self):
        """
        Return the Cartan type of ``self``.

        EXAMPLES::

            sage: B = crystals.Tableaux(['A',2], shape=[2,1])
            sage: psi = Hom(B, B).an_element()
            sage: psi.cartan_type()
            ['A', 2]
        """
        return self._cartan_type

    def is_injective(self):
        """
        Return if ``self`` is an injective crystal morphism.

        EXAMPLES::

            sage: B = crystals.Tableaux(['A',2], shape=[2,1])
            sage: psi = Hom(B, B).an_element()
            sage: psi.is_injective()
            False
        """
        # This is needed because is_injective is defined in a superclass, so
        #   we can't overwrite it with the category
        return self.is_embedding()

    def __call__(self, x, *args, **kwds):
        """
        Apply this map to ``x``. We need to do special processing
        for ``None``.

        EXAMPLES::

            sage: B = crystals.Tableaux(['A',2], shape=[2,1])
            sage: F = crystals.Tableaux(['A',2], shape=[1])
            sage: T = crystals.TensorProduct(F, F, F)
            sage: H = Hom(T, B)
            sage: b = B.module_generators[0]
            sage: psi = H((None, b, b, None), generators=T.highest_weight_vectors())
            sage: psi(None)
            sage: map(psi, T.highest_weight_vectors())
            [None, [[1, 1], [2]], [[1, 1], [2]], None]
        """
        if x is None:
            return None
        return super(CrystalMorphism, self).__call__(x, *args, **kwds)

class CrystalMorphismByGenerators(CrystalMorphism):
    """
    A crystal morphism defined by a set of generators.

    INPUT:

    - ``parent`` -- a homset
    - ``on_gens`` -- a function for the images of the generators
    - ``cartan_type`` -- (optional) the Cartan type of the morphism; the
      default is the Cartan type of the domain
    - ``index_set`` -- (optional) the index set of the morphism; the
      default is the index set of the Cartan type
    - ``gens`` -- a list of generators to define the morphism; the
      default is to use the module generators of the crystal
    - ``check`` -- (default: ``True``) check if the crystal morphism is valid

    If there is a fast way to find the module generators, overwrite
    :meth:`to_module_generator()` which currently looks for a module generator
    by a breath first search.
    """
    def __init__(self, parent, on_gens, cartan_type=None, index_set=None,
                 gens=None, check=True):
        """
        Construct a crystal morphism.

        TESTS::

            sage: B = crystals.Tableaux(['A',2], shape=[2,1])
            sage: H = Hom(B, B)
            sage: psi = H.an_element()
        """
        if gens is None:
            if isinstance(on_gens, dict):
                gens = on_gens.keys()
            else:
                gens = parent.domain().module_generators
        self._gens = tuple(gens)

        # Make sure on_gens is a function
        if isinstance(on_gens, dict):
            f = lambda x: on_gens[x]
        elif isinstance(on_gens, (list, tuple)):
            if len(self._gens) != len(on_gens):
                raise ValueError("not all module generators have an image")
            d = {x: y for x,y in zip(self._gens, on_gens)}
            f = lambda x: d[x]
        else:
            f = on_gens
        self._on_gens = f
        self._path_mg_cache = {x: (x, [], []) for x in self._gens}

        CrystalMorphism.__init__(self, parent, cartan_type, index_set)

        if check:
            self._check()

    def _check(self):
        """
        Check that ``self`` is a valid crystal morphism.

        EXAMPLES::

            sage: B = crystals.Tableaux(['A',2], shape=[2,1])
            sage: C = crystals.Tableaux(['A',2], shape=[4,2])
            sage: B.crystal_morphism(C.module_generators, codomain=C) # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: invalid crystal morphism: weights do not match

            sage: La = RootSystem(['A',2]).ambient_space().fundamental_weights()
            sage: T = crystals.TensorProduct(crystals.elementary.T(['A',2], La[1]+La[2]), B)
            sage: psi = T.crystal_morphism(C.module_generators, codomain=C)
        """
        index_set = self._index_set
        acx = self.domain().weight_lattice_realization().simple_coroots()
        acy = self.codomain().weight_lattice_realization().simple_coroots()
        for x in self._gens:
            y = self._on_gens(x)
            if y is None:
                continue
            xwt = x.weight()
            ywt = y.weight()
            if any(xwt.scalar(acx[i]) != ywt.scalar(acy[i]) for i in index_set):
                raise ValueError("invalid crystal morphism: weights do not match")
            if any(x.epsilon(i) != y.epsilon(i) for i in index_set):
                raise ValueError("invalid crystal morphism: epsilons do not match")
            if any(x.phi(i) != y.phi(i) for i in index_set):
                raise ValueError("invalid crystal morphism: phis do not match")

    def _call_(self, x):
        """
        Return the image of ``x`` under ``self``.

        EXAMPLES::

            sage: B = crystals.Tableaux(['A',2], shape=[2,1])
            sage: H = Hom(B, B)
            sage: psi = H(B.module_generators)
            sage: psi(B.highest_weight_vector())
            [[1, 1], [2]]
        """
        mg, ef, indices = self.to_module_generator(x)
        cur = self._on_gens(mg)
        for op,i in reversed(zip(ef, indices)):
            if cur is None:
                return None

            if op == 'e':
                cur = cur.f(i)
            elif op == 'f':
                cur = cur.e(i)
        return cur

    # TODO: Does this belong in the element_class of the Crystals() category?
    def to_module_generator(self, x):
        """
        Return a generator ``mg`` and a path of `e_i` and `f_i` operations
        to ``mg``.

        OUTPUT:

        A tuple consisting of:

        - a module generator,
        - a list of ``'e'`` and ``'f'`` to denote which operation, and
        - a list of matching indices.

        EXAMPLES::

            sage: B = crystals.elementary.Elementary(['A',2], 2)
            sage: psi = B.crystal_morphism(B.module_generators)
            sage: psi.to_module_generator(B(4))
            (0, ['f', 'f', 'f', 'f'], [2, 2, 2, 2])
            sage: psi.to_module_generator(B(-2))
            (0, ['e', 'e'], [2, 2])
        """
        if x in self._path_mg_cache:
            return self._path_mg_cache[x]

        mg = set(self._path_mg_cache.keys())
        visited = set([None, x])
        index_set = self._index_set
        todo = [x]
        ef = [[]]
        indices = [[]]

        while len(todo) > 0:
            cur = todo.pop(0)
            cur_ef = ef.pop(0)
            cur_indices = indices.pop(0)
            for i in index_set:
                next = cur.e(i)
                if next in mg:
                    gen,ef,indices = self._path_mg_cache[next]
                    ef = cur_ef + ['e'] + ef
                    indices = cur_indices + [i] + indices
                    self._path_mg_cache[x] = (gen, ef, indices)
                    return (gen, ef, indices)
                if next not in visited:
                    todo.append(next)
                    ef.append(cur_ef + ['e'])
                    indices.append(cur_indices + [i])
                    visited.add(next)

                # Now for f's
                next = cur.f(i)
                if next in mg:
                    gen,ef,indices = self._path_mg_cache[next]
                    ef = cur_ef + ['f'] + ef
                    indices = cur_indices + [i] + indices
                    self._path_mg_cache[x] = (gen, ef, indices)
                    return (gen, ef, indices)
                if next not in visited:
                    todo.append(next)
                    ef.append(cur_ef + ['f'])
                    indices.append(cur_indices + [i])
                    visited.add(next)
        raise ValueError("no module generator in the component of {}".format(x))

    def _repr_defn(self):
        """
        Used in constructing string representation of ``self``.

        EXAMPLES::

            sage: B = crystals.Tableaux(['A',2], shape=[2,1])
            sage: F = crystals.Tableaux(['A',2], shape=[1])
            sage: T = crystals.TensorProduct(F, F, F)
            sage: H = Hom(T, B)
            sage: b = B.highest_weight_vector()
            sage: psi = H((None, b, b, None), generators=T.highest_weight_vectors())
            sage: print psi._repr_defn()
            [[[1]], [[1]], [[1]]] |--> None
            [[[1]], [[2]], [[1]]] |--> [2, 1, 1]
            [[[2]], [[1]], [[1]]] |--> [2, 1, 1]
            [[[3]], [[2]], [[1]]] |--> None
        """
        return '\n'.join(['{} |--> {}'.format(mg, im)
                          for mg,im in zip(self._gens, self.im_gens())])

    @cached_method
    def im_gens(self):
        """
        Return the image of the generators of ``self`` as a tuple.

        EXAMPLES::

            sage: B = crystals.Tableaux(['A',2], shape=[2,1])
            sage: F = crystals.Tableaux(['A',2], shape=[1])
            sage: T = crystals.TensorProduct(F, F, F)
            sage: H = Hom(T, B)
            sage: b = B.highest_weight_vector()
            sage: psi = H((None, b, b, None), generators=T.highest_weight_vectors())
            sage: psi.im_gens()
            (None, [[1, 1], [2]], [[1, 1], [2]], None)
        """
        return tuple(map(self._on_gens, self._gens))

class TwistedCrystalMorphismByGenerators(CrystalMorphismByGenerators):
    """
    A crystal morphism defined by a set of generators twisted by
    an automorphism.

    INPUT:

    - ``parent`` -- a homset
    - ``on_gens`` -- a function for the images of the generators
    - ``automorphism`` -- the twisting automorphism
    - ``cartan_type`` -- (optional) the Cartan type of the morphism; the
      default is the Cartan type of the domain
    - ``index_set`` -- (optional) the index set of the morphism; the
      default is the index set of the Cartan type
    - ``gens`` -- a list of generators to define the morphism; the
      default is to use the module generators of the crystal
    - ``check`` -- (default: ``True``) check if the crystal morphism is valid
    """
    def __init__(self, parent, on_gens, automorphism, cartan_type=None,
                 index_set=None, gens=None, check=True):
        """
        Construct a twisted crystal morphism.

        TESTS::

            sage: B = crystals.Tableaux(['D',4], shape=[1])
            sage: H = Hom(B, B)
            sage: d = {1:1, 2:2, 3:4, 4:3}
            sage: psi = H(B.module_generators, automorphism=lambda x: d[x])
        """
        self._twist = automorphism
        CrystalMorphismByGenerators.__init__(self, parent, on_gens, cartan_type,
                                             index_set, gens, check)

    def _check(self):
        """
        Check if ``self`` is a valid twisted crystal morphism.

        TESTS::

            sage: B = crystals.Tableaux(['D',4], shape=[1])
            sage: H = Hom(B, B)
            sage: d = {1:1, 2:2, 3:4, 4:3}
            sage: psi = H(B.module_generators, automorphism=lambda x: d[x]) # indirect doctest
        """
        index_set = self._index_set
        acx = self.domain().weight_lattice_realization().simple_coroots()
        acy = self.codomain().weight_lattice_realization().simple_coroots()
        for x in self._gens:
            y = self._on_gens(x)
            if y is None:
                continue
            xwt = x.weight()
            ywt = y.weight()
            if any(xwt.scalar(acx[self._twist(i)]) != ywt.scalar(acy[i]) for i in index_set):
                raise ValueError("invalid crystal morphism: weights do not match")
            if any(x.epsilon(self._twist(i)) != y.epsilon(i) for i in index_set):
                raise ValueError("invalid crystal morphism: epsilons do not match")
            if any(x.phi(self._twist(i)) != y.phi(i) for i in index_set):
                raise ValueError("invalid crystal morphism: phis do not match")

    def _call_(self, x):
        """
        Return the image of ``x`` under ``self``.

        EXAMPLES::

            sage: B = crystals.Tableaux(['D',4], shape=[1])
            sage: H = Hom(B, B)
            sage: d = {1:1, 2:2, 3:4, 4:3}
            sage: psi = H(B.module_generators, automorphism=lambda x: d[x])
            sage: b = B.highest_weight_vector()
            sage: psi(b.f_string([1,2,3]))
            [[-4]]
            sage: psi(b.f_string([1,2,4]))
            [[4]]
        """
        mg, ef, indices = self.to_module_generator(x)
        cur = self._on_gens(mg)
        for op,i in reversed(zip(ef, indices)):
            if cur is None:
                return None
            if op == 'e':
                cur = cur.f(self._twist(i))
            elif op == 'f':
                cur = cur.e(self._twist(i))
        return cur

class VirtualCrystalMorphismByGenerators(CrystalMorphismByGenerators):
    r"""
    A crystal morphism defined by a set of generators which create a virtual
    crystal inside the codomain.

    INPUT:

    - ``parent`` -- a homset
    - ``on_gens`` -- a function for the images of the generators
    - ``virtualization`` -- a dictionary whose keys are in the index set of
      the domain and whose values are lists of entries in the index set of the
      codomain
    - ``scaling_factors`` -- a dictionary whose keys are in the index set of
      the domain and whose values are scaling factors for the weight,
      `\varepsilon` and `\varphi`.
    - ``cartan_type`` -- (optional) the Cartan type of the morphism; the
      default is the Cartan type of the domain
    - ``index_set`` -- (optional) the index set of the morphism; the
      default is the index set of the Cartan type
    - ``gens`` -- a list of generators to define the morphism; the
      default is to use the module generators of the crystal
    - ``check`` -- (default: ``True``) check if the crystal morphism is valid
    """
    def __init__(self, parent, on_gens, virtualization=None, scaling_factors=None,
                 cartan_type=None, index_set=None, gens=None, check=True):
        """
        Construct a virtual crystal morphism.

        TESTS::

            sage: B = crystals.Tableaux(['B',3], shape=[1])
            sage: C = crystals.Tableaux(['D',4], shape=[2])
            sage: H = Hom(B, C)
            sage: psi = H(C.module_generators)
        """
        # We skip the check on the initialization because we may not have initialized everything
        CrystalMorphismByGenerators.__init__(self, parent, on_gens, cartan_type,
                                             index_set, gens, False)
        if scaling_factors is None:
            scaling_factors = {i: 1 for i in index_set}
        if virtualization is None:
            virtualization = {i: (i,) for i in index_set}
        self._virtualization = virtualization
        self._scaling_factors = scaling_factors

        # Now that everything is initialized, run the check (if it is wanted)
        if check:
            self._check()

    def _check(self):
        """
        Check if ``self`` is a valid virtual crystal morphism.

        TESTS::

            sage: B = crystals.Tableaux(['B',3], shape=[1])
            sage: C = crystals.Tableaux(['D',4], shape=[2])
            sage: H = Hom(B, C)
            sage: psi = H(C.module_generators) # indirect doctest
        """
        index_set = self._index_set
        acx = self.domain().weight_lattice_realization().simple_coroots()
        La = self.codomain().weight_lattice_realization().fundamental_weights()
        v = self._virtualization
        sf = self._scaling_factors
        for x in self._gens:
            y = self._on_gens(x)
            if y is None:
                continue
            xwt = x.weight()
            xwt = {i: xwt.scalar(acx[i]) for i in index_set}
            if y.weight() != sum(xwt[i] * sf[i] * La[j] for i in index_set for j in v[i]):
                raise ValueError("invalid crystal morphism: weights do not match")
            if any(sf[i] * x.epsilon(i) != y.epsilon(j) for i in index_set for j in v[i]):
                raise ValueError("invalid crystal morphism: epsilons are not aligned")
            if any(sf[i] * x.phi(i) != y.phi(j) for i in index_set for j in v[i]):
                raise ValueError("invalid crystal morphism: phis are not aligned")

    def _call_(self, x):
        """
        Return the image of ``x`` under ``self``.

        EXAMPLES::

            sage: B = crystals.Tableaux(['B',3], shape=[1])
            sage: C = crystals.Tableaux(['D',4], shape=[2])
            sage: H = Hom(B, C)
            sage: psi = H(C.module_generators)
            sage: psi(B.highest_weight_vector())
            [[1, 1]]
        """
        mg, ef, indices = self.to_module_generator(x)
        cur = self._on_gens(mg)
        for op,i in reversed(zip(ef, indices)):
            if cur is None:
                return None

            s = []
            sf = self._scaling_factors[i]
            for j in self._virtualization[i]:
                s += [j]*sf
            if op == 'e':
                cur = cur.f_string(s)
            elif op == 'f':
                cur = cur.e_string(s)
        return cur

    def image(self):
        """
        Return the image of ``self`` in the codomain as a
        :class:`~sage.combinat.crystals.virtual_crystal.VirtualCrystal`.

        .. WARNING::

            This assumes that ``self`` is a strict crystal morphism.

        .. TODO::

            A more general image as a (virtual) subcrystal and provide a
            better mechanic to construct virtual crystals.

        EXAMPLES::

            sage: B = crystals.Tableaux(['B',3], shape=[1])
            sage: C = crystals.Tableaux(['D',4], shape=[2])
            sage: H = Hom(B, C)
            sage: psi = H(C.module_generators)
            sage: psi.image()
            Virtual crystal of The crystal of tableaux of type ['D', 4] and shape(s) [[2]] of type ['B', 3]
        """
        #if not self.is_strict():
        #    raise NotImplementedError
        from sage.combinat.crystals.virtual_crystal import VirtualCrystal
        return VirtualCrystal(self.codomain(), self._virtualization, self._scaling_factors,
                              None, self.im_gens(), self.cartan_type(), self.index_set(),
                              self.domain().category())

###############################################################################
## Homset

class CrystalHomset(Homset):
    r"""
    The set of crystal morphisms from one crystal to another.

    An `U_q(\mathfrak{g})` `I`-crystal morphism `\Psi : B \to C` is a map
    `\Psi : B \cup \{ 0 \} \to C \cup \{ 0 \}` such that:

    - `\Psi(0) = 0`.
    - If `b \in B` and `\Psi(b) \in C`, then
      `\mathrm{wt}(\Psi(b)) = \mathrm{wt}(b)`,
      `\varepsilon_i(\Psi(b)) = \varepsilon_i(b)`, and
      `\varphi_i(\Psi(b)) = \varphi_i(b)` for all `i \in I`.
    - If `b, b^{\prime} \in B`, `\Psi(b), \Psi(b^{\prime}) \in C` and
      `f_i b = b^{\prime}`, then `f_i \Psi(b) = \Psi(b^{\prime})` and
      `\Psi(b) = e_i \Psi(b^{\prime})` for all `i \in I`.

    If the Cartan type is unambiguous, it is often surpressed from the
    notation. Similarly if the full index set is used, it is also surpressed.

    We can also generalize the definition of a crystal morphism by considering
    a map of `\sigma` of the (now possibly different) Dynkin diagrams `B` and
    `C` a along with scaling factors `\gamma_i \in \ZZ` for `i \in I`.
    Let `\sigma_i` denote the orbit of `i` under `\sigma`. Then a *virtual*
    crystal morphism `\Psi` is a map `\Psi` such that the following holds:

    - `\Psi(0) = 0`.
    - If `b \in B` and `\Psi(b) \in C`, then for all `j \in \sigma_i`:

    .. MATH::

        \varepsilon_i(b) = \frac{\varepsilon_j}{\gamma_j} \Psi(b), \quad
        \widehat{\varphi}_i(b) = \frac{\varphi_j}{\gamma_j} \Psi(b), \quad
        \mathrm{wt}(\Phi(b)) = c_i \sum_{j \in \sigma_i} \gamma_j \Lambda_j.

    where `\mathrm{wt}(b) = \sum_i c_i \Lambda_i`.

    - If `b, b^{\prime} \in B`, `\Psi(b), \Psi(b^{\prime}) \in C` and
      `f_i b = b^{\prime}`, then independent of the ordering of `\sigma_i`
      we have:

      .. MATH::

          \Psi(b^{\prime}) = \widehat{e}_i \Psi(b) =
              \prod_{j \in \sigma_i} e_j^{\gamma_i} \Psi(b), \quad
          \Psi(b^{\prime}) = \widehat{f}_i \Psi(b) =
              \prod_{j \in \sigma_i} f_j^{\gamma_i} \Psi(b).

    If `\gamma_i = 1` for all `I` and the Dynkin diagrams are the same, then
    we call `\Psi` a *twisted* crystal morphism.

    INPUT:

    - ``X`` -- the domain
    - ``Y`` -- the codomain
    - ``category`` -- (optional) the category of the crystal morphisms

    EXAMPLES:

    We begin with the natural embedding of `B(2\Lambda_1)` into
    `B(\Lambda_1) \otimes B(\Lambda_1)` in type `A_1`::

        sage: B = crystals.Tableaux(['A',1], shape=[2])
        sage: F = crystals.Tableaux(['A',1], shape=[1])
        sage: T = crystals.TensorProduct(F, F)
        sage: v = T.highest_weight_vectors()[0]; v
        [[[1]], [[1]]]
        sage: H = Hom(B, T)
        sage: psi = H([v])
        sage: b = B.highest_weight_vector(); b
        [[1, 1]]
        sage: psi(b)
        [[[1]], [[1]]]
        sage: b.f(1)
        [[1, 2]]
        sage: psi(b.f(1))
        [[[1]], [[2]]]

    We now look at the decomposition of `B(\Lambda_1) \otimes B(\Lambda_1)`
    into `B(2\Lambda_1) \oplus B(0)`::

        sage: B0 = crystals.Tableaux(['A',1], shape=[])
        sage: D = crystals.DirectSum([B, B0])
        sage: H = Hom(T, D)
        sage: psi = H(D.module_generators)
        sage: psi
        ['A', 1] Crystal morphism:
          From: Full tensor product of the crystals
           [The crystal of tableaux of type ['A', 1] and shape(s) [[1]],
            The crystal of tableaux of type ['A', 1] and shape(s) [[1]]]
          To:   Direct sum of the crystals Family
           (The crystal of tableaux of type ['A', 1] and shape(s) [[2]],
            The crystal of tableaux of type ['A', 1] and shape(s) [[]])
          Defn: [[[1]], [[1]]] |--> [1, 1]
                [[[2]], [[1]]] |--> []
        sage: psi.is_isomorphism()
        True

    We can always construct the trivial morphism which sends
    everything to `0`::

        sage: Binf = crystals.infinity.Tableaux(['B', 2])
        sage: B = crystals.Tableaux(['B',2], shape=[1])
        sage: H = Hom(Binf, B)
        sage: psi = H(lambda x: None)
        sage: psi(Binf.highest_weight_vector())

    For Kirillov-Reshetikhin crystals, we consider the map to the
    corresponding classical crystal::

        sage: K = crystals.KirillovReshetikhin(['D',4,1], 2,1)
        sage: B = K.classical_decomposition()
        sage: H = Hom(K, B)
        sage: psi = H(lambda x: x.lift(), cartan_type=['D',4])
        sage: L = map(psi, K.module_generators); L
        [[], [[1], [2]]]
        sage: all(x.parent() == B for x in L)
        True

    Next we consider a type `D_4` crystal morphism where we twist by
    `3 \leftrightarrow 4`::

        sage: B = crystals.Tableaux(['D',4], shape=[1])
        sage: H = Hom(B, B)
        sage: d = {1:1, 2:2, 3:4, 4:3}
        sage: psi = H(B.module_generators, automorphism=lambda x: d[x])
        sage: b = B.highest_weight_vector()
        sage: b.f_string([1,2,3])
        [[4]]
        sage: b.f_string([1,2,4])
        [[-4]]
        sage: psi(b.f_string([1,2,3]))
        [[-4]]
        sage: psi(b.f_string([1,2,4]))
        [[4]]

    We construct the natural virtual embedding of a type `B_3` into a type
    `D_4` crystal::

        sage: B = crystals.Tableaux(['B',3], shape=[1])
        sage: C = crystals.Tableaux(['D',4], shape=[2])
        sage: H = Hom(B, C)
        sage: psi = H(C.module_generators)
        sage: psi
        ['B', 3] Crystal morphism:
          From: The crystal of tableaux of type ['B', 3] and shape(s) [[1]]
          To:   The crystal of tableaux of type ['D', 4] and shape(s) [[2]]
          Defn: [1] |--> [1, 1]
        sage: for b in B: print "{} |--> {}".format(b, psi(b))
        [1] |--> [1, 1]
        [2] |--> [2, 2]
        [3] |--> [3, 3]
        [0] |--> [3, -3]
        [-3] |--> [-3, -3]
        [-2] |--> [-2, -2]
        [-1] |--> [-1, -1]
    """
    def __init__(self, X, Y, category=None):
        """
        Initialize ``self``.

        TESTS::

            sage: B = crystals.Tableaux(['A', 2], shape=[2,1])
            sage: H = Hom(B, B)
            sage: Binf = crystals.infinity.Tableaux(['B',2])
            sage: H = Hom(Binf, B)
        """
        if category is None:
            category = Crystals()
        # TODO: Should we make one of the types of morphisms into the self.Element?
        Homset.__init__(self, X, Y, category)

    _generic = CrystalMorphismByGenerators
    _twisted = TwistedCrystalMorphismByGenerators
    _virtual = VirtualCrystalMorphismByGenerators

    def _repr_(self):
        """
        TESTS::

            sage: B = crystals.Tableaux(['A', 2], shape=[2,1])
            sage: Hom(B, B)
            Set of Crystal Morphisms from The crystal of tableaux of type ['A', 2] and shape(s) [[2, 1]]
             to The crystal of tableaux of type ['A', 2] and shape(s) [[2, 1]]
        """
        return "Set of Crystal Morphisms from {} to {}".format(self.domain(), self.codomain())

    def _coerce_impl(self, x):
        """
        Check to see if we can coerce ``x`` into a morphism with the
        correct parameters.

        EXAMPLES::

            sage: B = crystals.Tableaux(['B',3], shape=[2,1])
            sage: H = Hom(B, B)
            sage: H(H.an_element()) # indirect doctest
            ['B', 3] Crystal endomorphism of The crystal of tableaux of type ['B', 3] and shape(s) [[2, 1]]
              Defn: [2, 1, 1] |--> None
        """
        if not isinstance(x, CrystalMorphism):
            raise TypeError

        if x.parent() is self:
            return x

        # Case 1: the parent fits
        if x.parent() == self:
            if isinstance(x, CrystalMorphismByGenerators):
                return self._generic(self, x._on_gens, x._cartan_type, x._index_set, x._gens)
            if isinstance(x, TwistedCrystalMorphismByGenerators):
                return self._twisted(self, x._on_gens, x._automorphism,
                                     x._cartan_type, x._index_set, x._gens)
            if isinstance(x, VirtualCrystalMorphismByGenerators):
                return self._virtual(self, x._on_gens, x._virtualization, x._scaling_factors,
                                     x._cartan_type, x._index_set, x._gens)

        # TODO: Should we try extraordinary measures (like twisting)?
        raise ValueError

    def __call__(self, on_gens, cartan_type=None, index_set=None, generators=None,
                 automorphism=None, virtualization=None, scaling_factors=None, check=True):
        """
        Construct a crystal morphism.

        EXAMPLES::

            sage: B = crystals.Tableaux(['A', 2], shape=[2,1])
            sage: H = Hom(B, B)
            sage: psi = H(B.module_generators)

            sage: F = crystals.Tableaux(['A',3], shape=[1])
            sage: T = crystals.TensorProduct(F, F, F)
            sage: H = Hom(B, T)
            sage: v = T.highest_weight_vectors()[2]
            sage: psi = H([v], cartan_type=['A',2])
        """
        if isinstance(on_gens, CrystalMorphism):
            return self._coerce_impl(on_gens)

        if cartan_type is None:
            cartan_type = self.domain().cartan_type()
        else:
            from sage.combinat.root_system.cartan_type import CartanType
            cartan_type = CartanType(cartan_type)
        if index_set is None:
            index_set = cartan_type.index_set()

        # Check for twisted morphisms
        if automorphism is not None:
            return self._twisted(self, on_gens, automorphism,
                                 cartan_type, index_set, generators, check)

        # Check for virtual morphisms
        if virtualization is not None or scaling_factors is not None:
            return self._virtual(self, on_gens, virtualization, scaling_factors,
                                 cartan_type, index_set, generators, check)

        # Try as a natural folding
        if cartan_type != self.codomain().cartan_type():
            fct = self.domain().cartan_type().as_folding()
            if fct.folding_of() == self.codomain().cartan_type():
                return self._virtual(self, on_gens, fct.folding_orbit(), fct.scaling_factors(),
                                     cartan_type, index_set, generators, check)

        # It is a regular morphism
        return self._generic(self, on_gens, cartan_type, index_set, generators, check)

    def _an_element_(self):
        """
        Return an element of ``self``. Every homset has the crystal morphism
        which sends all elements to ``None``.

        EXAMPLES::

            sage: B = crystals.Tableaux(['A', 2], shape=[2,1])
            sage: C = crystals.infinity.Tableaux(['A', 2])
            sage: H = Hom(B, C)
            sage: H.an_element()
            ['A', 2] Crystal morphism:
              From: The crystal of tableaux of type ['A', 2] and shape(s) [[2, 1]]
              To:   The infinity crystal of tableaux of type ['A', 2]
              Defn: [2, 1, 1] |--> None
        """
        return self._generic(self, lambda x: None)

