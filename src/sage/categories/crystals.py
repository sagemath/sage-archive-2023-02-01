r"""
Crystals
"""
#*****************************************************************************
#  Copyright (C) 2010    Anne Schilling <anne at math.ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import CachedFunction, cached_method
from sage.misc.abstract_method import abstract_method
from sage.categories.category_singleton import Category_singleton
from sage.categories.enumerated_sets import EnumeratedSets
from sage.categories.category import HomCategory
from sage.categories.morphism import Morphism
from sage.categories.homset import Hom, Homset
from sage.misc.latex import latex
from sage.combinat import ranker
from sage.graphs.dot2tex_utils import have_dot2tex
from sage.rings.integer import Integer

class Crystals(Category_singleton):
    r"""
    The category of crystals

    See :mod:`sage.combinat.crystals` for an introduction to crystals.

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
        {'required': ['e', 'epsilon', 'f', 'phi', 'weight'], 'optional': []}

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
                """

            @abstract_method(optional=True)
            def is_injective(self):
                """
                Check if ``self`` is an ijective crystal morphism.
                """

            @abstract_method(optional=True)
            def is_strict(self):
                """
                Check if ``self`` is a strict crystal morphism.
                """

    class ParentMethods:

        def an_element(self):
            """
            Returns an element of ``self``

                sage: C = CrystalOfLetters(['A', 5])
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

                sage: C = CrystalOfLetters(['A', 5])
                sage: C.weight_lattice_realization()
                Ambient space of the Root system of type ['A', 5]
                sage: K = KirillovReshetikhinCrystal(['A',2,1], 1, 1)
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

                sage: C = CrystalOfLetters(['A',2])
                sage: C.cartan_type()
                ['A', 2]
            """
            return self._cartan_type

        @cached_method
        def index_set(self):
            """
            Returns the index set of the Dynkin diagram underlying the crystal

            EXAMPLES::

                sage: C = CrystalOfLetters(['A', 5])
                sage: C.index_set()
                (1, 2, 3, 4, 5)
            """
            return self.cartan_type().index_set()

        def Lambda(self):
            """
            Returns the fundamental weights in the weight lattice
            realization for the root system associated with the crystal

            EXAMPLES::

                sage: C = CrystalOfLetters(['A', 5])
                sage: C.Lambda()
                Finite family {1: (1, 0, 0, 0, 0, 0), 2: (1, 1, 0, 0, 0, 0), 3: (1, 1, 1, 0, 0, 0), 4: (1, 1, 1, 1, 0, 0), 5: (1, 1, 1, 1, 1, 0)}
            """
            return self.weight_lattice_realization().fundamental_weights()

        def __iter__(self, index_set=None, max_depth=float('inf')):
            """
            Returns the iterator of ``self``.

            INPUT:

            - ``index_set`` -- (Default: ``None``) The index set; if ``None``
              then use the index set of the crystal

            - ``max_depth`` -- (Default: infinity) The maximum depth to build

            EXAMPLES::

                sage: C = CrystalOfLSPaths(['A',2,1],[-1,0,1])
                sage: C.__iter__.__module__
                'sage.categories.crystals'
                sage: g = C.__iter__()
                sage: g.next()
                (-Lambda[0] + Lambda[2],)
                sage: g.next()
                (Lambda[0] - Lambda[1] + delta,)
                sage: g.next()
                (Lambda[1] - Lambda[2] + delta,)
                sage: g.next()
                (Lambda[1] - Lambda[2],)
                sage: g.next()
                (Lambda[0] - Lambda[1],)
                sage: h = C.__iter__(index_set=[1,2])
                sage: h.next()
                (-Lambda[0] + Lambda[2],)
                sage: h.next()
                (Lambda[1] - Lambda[2],)
                sage: h.next()
                (Lambda[0] - Lambda[1],)
                sage: h.next()
                Traceback (most recent call last):
                ...
                StopIteration
                sage: g = C.__iter__(max_depth=1)
                sage: g.next()
                (-Lambda[0] + Lambda[2],)
                sage: g.next()
                (Lambda[1] - Lambda[2],)
                sage: g.next()
                (Lambda[0] - Lambda[1] + delta,)
                sage: h.next()
                Traceback (most recent call last):
                ...
                StopIteration

            """
            if index_set is None:
                index_set = self.index_set()
            if max_depth < float('inf'):
                from sage.combinat.backtrack import TransitiveIdealGraded
                return TransitiveIdealGraded(lambda x: [x.f(i) for i in index_set]
                                                     + [x.e(i) for i in index_set],
                                             self.module_generators, max_depth).__iter__()
            from sage.combinat.backtrack import TransitiveIdeal
            return TransitiveIdeal(lambda x: [x.f(i) for i in index_set]
                                           + [x.e(i) for i in index_set],
                                   self.module_generators).__iter__()

        def subcrystal(self, index_set=None, generators=None, max_depth=float("inf"),
                       direction="both", contained=None, category=None):
            r"""
            Construct the subcrystal from ``generators`` using `e_i` and `f_i`
            for all `i` in ``index_set``.

            INPUT:

            - ``index_set`` -- (Default: ``None``) The index set; if ``None``
              then use the index set of the crystal

            - ``generators`` -- (Default: ``None``) The list of generators; if
              ``None`` then use the module generators of the crystal

            - ``max_depth`` -- (Default: infinity) The maximum depth to build

            - ``direction`` -- (Default: ``'both'``) The direction to build
              the subcrystal. It can be one of the following:

              - ``'both'`` - Using both `e_i` and `f_i`
              - ``'upper'`` - Using `e_i`
              - ``'lower'`` - Using `f_i`

            - ``contained`` -- (Optional) A function defining the containment
              in the subcrystal

            - ``category`` -- (Optional) Specify the category of the subcrystal

            EXAMPLES::

                sage: C = KirillovReshetikhinCrystal(['A',3,1], 1, 2)
                sage: S = list(C.subcrystal(index_set=[1,2])); S
                [[[1, 1]], [[1, 2]], [[1, 3]], [[2, 2]], [[2, 3]], [[3, 3]]]
                sage: C.cardinality()
                10
                sage: len(S)
                6
                sage: list(C.subcrystal(index_set=[1,3], generators=[C(1,4)]))
                [[[1, 4]], [[2, 4]], [[1, 3]], [[2, 3]]]
                sage: list(C.subcrystal(index_set=[1,3], generators=[C(1,4)], max_depth=1))
                [[[1, 4]], [[2, 4]], [[1, 3]]]
                sage: list(C.subcrystal(index_set=[1,3], generators=[C(1,4)], direction='upper'))
                [[[1, 4]], [[1, 3]]]
                sage: list(C.subcrystal(index_set=[1,3], generators=[C(1,4)], direction='lower'))
                [[[1, 4]], [[2, 4]]]
            """
            from sage.combinat.crystals.subcrystal import Subcrystal
            from sage.categories.finite_crystals import FiniteCrystals

            if index_set is None:
                index_set = self.index_set()
            if generators is None:
                generators = self.module_generators

            if max_depth == float('inf'):
                if self not in FiniteCrystals():
                    if contained is None and index_set == self.index_set() \
                            and generators == self.module_generators:
                        return self
                    return Subcrystal(self, contained, generators,
                                      index_set=index_set, category=category)
                if contained is None and direction == 'both':
                    if category is None:
                        category = FiniteCrystals()
                    return Subcrystal(self, contained, generators,
                                      index_set=index_set, category=category)

            from sage.combinat.backtrack import TransitiveIdealGraded
            if direction == 'both':
                subset = TransitiveIdealGraded(lambda x: [x.f(i) for i in index_set]
                                                     + [x.e(i) for i in index_set],
                                             generators, max_depth)
            elif direction == 'upper':
                subset = TransitiveIdealGraded(lambda x: [x.e(i) for i in index_set],
                                             generators, max_depth)
            elif direction == 'lower':
                subset = TransitiveIdealGraded(lambda x: [x.f(i) for i in index_set],
                                             generators, max_depth)
            else:
                raise ValueError("direction must be either 'both', 'upper', or 'lower'")

            subset = frozenset(subset)
            if category is None:
                category = FiniteCrystals()

            if self in FiniteCrystals() and len(subset) == self.cardinality():
                if contained is None and index_set == self.index_set():
                    return self
                return Subcrystal(self, contained, generators,
                                  index_set=index_set, category=category)

            if contained is not None:
                contained = lambda x: x in subset and contained(x)
            else:
                contained = lambda x: x in subset
            return Subcrystal(self, contained, generators, index_set=index_set, category=category)

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

                sage: B = CrystalOfTableaux(['A',2], shape=[2,1])
                sage: H = B._Hom_(B, category=Crystals()); H
                Set of Homomorphisms from Rational Field to Rational Field
            """
            if category is None:
                category = self.category()
            elif not category.is_subcategory(Crystals()):
                raise TypeError("{} is not a subcategory of Crystals()".format(category))
            if Y not in Crystals():
                raise TypeError("{} is not a crystal".format(Y))
            return CrystalHomset(self, Y, category=category, **options)

        def crystal_morphism(self, g, index_set = None, automorphism = lambda i : i, direction = 'down', direction_image = 'down',
                             similarity_factor = None, similarity_factor_domain = None, cached = False, acyclic = True):
            r"""
            Constructs a morphism from the crystal ``self`` to another crystal.
            The input `g` can either be a function of a (sub)set of elements of self to
            element in another crystal or a dictionary between certain elements.
            Usually one would map highest weight elements or crystal generators to each
            other using g.
            Specifying index_set gives the opportunity to define the morphism as `I`-crystals
            where `I =` index_set. If index_set is not specified, the index set of self is used.
            It is also possible to define twisted-morphisms by specifying an automorphism on the
            nodes in te Dynkin diagram (or the index_set).
            The option direction and direction_image indicate whether to use `f_i` or `e_i` in
            self or the image crystal to construct the morphism, depending on whether the direction
            is set to 'down' or 'up'.
            It is also possible to set a similarity_factor. This should be a dictionary between
            the elements in the index set and positive integers. The crystal operator `f_i` then gets
            mapped to `f_i^{m_i}` where `m_i =` similarity_factor[i].
            Setting similarity_factor_domain to a dictionary between the index set and positive integers
            has the effect that `f_i^{m_i}` gets mapped to `f_i` where `m_i =` similarity_factor_domain[i].
            Finally, it is possible to set the option `acyclic = False`. This calculates an isomorphism
            for cyclic crystals (for example finite affine crystals). In this case the input function `g`
            is supposed to be given as a dictionary.

            EXAMPLES::

                sage: C2 = CrystalOfLetters(['A',2])
                sage: C3 = CrystalOfLetters(['A',3])
                sage: g = {C2.module_generators[0] : C3.module_generators[0]}
                sage: g_full = C2.crystal_morphism(g)
                sage: g_full(C2(1))
                1
                sage: g_full(C2(2))
                2
                sage: g = {C2(1) : C2(3)}
                sage: g_full = C2.crystal_morphism(g, automorphism = lambda i : 3-i, direction_image = 'up')
                sage: [g_full(b) for b in C2]
                [3, 2, 1]
                sage: T = CrystalOfTableaux(['A',2], shape = [2])
                sage: g = {C2(1) : T(rows=[[1,1]])}
                sage: g_full = C2.crystal_morphism(g, similarity_factor = {1:2, 2:2})
                sage: [g_full(b) for b in C2]
                [[[1, 1]], [[2, 2]], [[3, 3]]]
                sage: g = {T(rows=[[1,1]]) : C2(1)}
                sage: g_full = T.crystal_morphism(g, similarity_factor_domain = {1:2, 2:2})
                sage: g_full(T(rows=[[2,2]]))
                2

                sage: B1=KirillovReshetikhinCrystal(['A',2,1],1,1)
                sage: B2=KirillovReshetikhinCrystal(['A',2,1],1,2)
                sage: T=TensorProductOfCrystals(B1,B2)
                sage: T1=TensorProductOfCrystals(B2,B1)
                sage: La = T.weight_lattice_realization().fundamental_weights()
                sage: t = [b for b in T if b.weight() == -3*La[0] + 3*La[1]][0]
                sage: t1 = [b for b in T1 if b.weight() == -3*La[0] + 3*La[1]][0]
                sage: g={t:t1}
                sage: f=T.crystal_morphism(g,acyclic = False)
                sage: [[b,f(b)] for b in T]
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
            """
            if index_set is None:
                index_set = self.index_set()
            if similarity_factor is None:
                similarity_factor = dict( (i,1) for i in index_set )
            if similarity_factor_domain is None:
                similarity_factor_domain = dict( (i,1) for i in index_set )
            if direction == 'down':
                e_string = 'e_string'
            else:
                e_string = 'f_string'
            if direction_image == 'down':
                f_string = 'f_string'
            else:
                f_string = 'e_string'

            if acyclic:
                if type(g) == dict:
                    g = g.__getitem__

                def morphism(b):
                    for i in index_set:
                        c = getattr(b, e_string)([i for k in range(similarity_factor_domain[i])])
                        if c is not None:
                            d = getattr(morphism(c), f_string)([automorphism(i) for k in range(similarity_factor[i])])
                            if d is not None:
                                return d
                            else:
                                raise ValueError, "This is not a morphism!"
                            #now we know that b is hw
                    return g(b)

            else:
                import copy
                morphism = copy.copy(g)
                known = set( g.keys() )
                todo = copy.copy(known)
                images = set( [g[x] for x in known] )
                # Invariants:
                # - images contains all known morphism(x)
                # - known contains all elements x for which we know morphism(x)
                # - todo  contains all elements x for which we haven't propagated to each child
                while todo <> set( [] ):
                    x = todo.pop()
                    for i in index_set:
                        eix  = getattr(x, f_string)([i for k in range(similarity_factor_domain[i])])
                        eigx = getattr(morphism[x], f_string)([automorphism(i) for k in range(similarity_factor[i])])
                        if bool(eix is None) <> bool(eigx is None):
                            # This is not a crystal morphism!
                            raise ValueError, "This is not a morphism!" #, print("x="x,"g(x)="g(x),"i="i)
                        if (eix is not None) and (eix not in known):
                            todo.add(eix)
                            known.add(eix)
                            morphism[eix] = eigx
                            images.add(eigx)
                # Check that the module generators are indeed module generators
                assert(len(known) == self.cardinality())
                # We may not want to systematically run those tests,
                # to allow for non bijective crystal morphism
                # Add an option CheckBijective?
                if not ( len(known) == len(images) and len(images) == images.pop().parent().cardinality() ):
                    return(None)
                return morphism.__getitem__

            if cached:
                return morphism
            else:
                return CachedFunction(morphism)

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
                sage: T = HighestWeightCrystal(La[0])
                sage: S = T.subcrystal(max_depth=3)
                sage: G = T.digraph(subset=S); G
                Digraph on 5 vertices
                sage: G.vertices()
                [(1/2*Lambda[0] + Lambda[1] - Lambda[2] - 1/2*delta, -1/2*Lambda[0] + Lambda[1] - 1/2*delta),
                (-Lambda[0] + 2*Lambda[1] - delta,), (Lambda[0] - 2*Lambda[1] + 2*Lambda[2] - delta,),
                (1/2*Lambda[0] - Lambda[1] + Lambda[2] - 1/2*delta, -1/2*Lambda[0] + Lambda[1] - 1/2*delta), (Lambda[0],)]

            Here is a way to construct a picture of a Demazure crystal using
            the ``subset`` option::

                sage: B = CrystalOfTableaux(['A',2], shape=[2,1])
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

                sage: C = KirillovReshetikhinCrystal(['D',4,1], 2, 1)
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
                f = lambda (u,v,label): ({})
            else:
                f = lambda (u,v,label): ({"backward":label ==0})

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

                sage: C = CrystalOfLetters(['A', 5])
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

                sage: T = CrystalOfTableaux(['A',2],shape=[1])
                sage: T._latex_()   #optional - dot2tex
                '...tikzpicture...'
                sage: view(T, pdflatex = True, tightpage = True) #optional - dot2tex graphviz

            One can for example also color the edges using the following options::

                sage: T = CrystalOfTableaux(['A',2],shape=[1])
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

                sage: C = CrystalOfLetters(['A', 2])
                sage: C.metapost('/tmp/test.mp') #optional

            ::

                sage: C = CrystalOfLetters(['A', 5])
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
                    while not turtle.e(j) == None:
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
                    if not dest == None:
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
            Returns a dot_tex string representation of ``self``.

            EXAMPLES::

                sage: C = CrystalOfLetters(['A',2])
                sage: C.dot_tex()
                'digraph G { \n  node [ shape=plaintext ];\n  N_0 [ label = " ", texlbl = "$1$" ];\n  N_1 [ label = " ", texlbl = "$2$" ];\n  N_2 [ label = " ", texlbl = "$3$" ];\n  N_0 -> N_1 [ label = " ", texlbl = "1" ];\n  N_1 -> N_2 [ label = " ", texlbl = "2" ];\n}'
            """
            import re
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
            Returns the plot of self as a directed graph.

            EXAMPLES::

                sage: C = CrystalOfLetters(['A', 5])
                sage: print(C.plot())
                Graphics object consisting of 17 graphics primitives
            """
            return self.digraph().plot(edge_labels=True,vertex_size=0,**options)

        def plot3d(self, **options):
            """
            Returns the 3-dimensional plot of self as a directed graph.

            EXAMPLES::

                sage: C = KirillovReshetikhinCrystal(['A',3,1],2,1)
                sage: print(C.plot3d())
                Graphics3d Object
            """
            G = self.digraph(**options)
            return G.plot3d()

        def direct_sum(self, X):
            """
            Return the direct sum of ``self`` with ``X``.

            EXAMPLES::

                sage: B = CrystalOfTableaux(['A',2], shape=[2,1])
                sage: C = CrystalOfLetters(['A',2])
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

                sage: B = CrystalOfTableaux(['A',2], shape=[2,1])
                sage: C = CrystalOfLetters(['A',2])
                sage: T = TensorProductOfCrystals(B,C)
                sage: T.connected_components_generators()
                ([[[1, 1], [2]], 1], [[[1, 2], [2]], 1], [[[1, 2], [3]], 1])
            """

        def connected_components(self):
            """
            Return the connected components of ``self`` as subcrystals.

            EXAMPLES::

                sage: B = CrystalOfTableaux(['A',2], shape=[2,1])
                sage: C = CrystalOfLetters(['A',2])
                sage: T = TensorProductOfCrystals(B,C)
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

        def connected_components_number(self):
            """
            Return the number of connected components of ``self``.

            EXAMPLES::

                sage: B = CrystalOfTableaux(['A',2], shape=[2,1])
                sage: C = CrystalOfLetters(['A',2])
                sage: T = TensorProductOfCrystals(B,C)
                sage: T.connected_components_number()
                3
            """
            return len(self.connected_components_generators())

        def is_connected(self):
            """
            Return ``True`` if ``self`` is a connected crystal.

            EXAMPLES::

                sage: B = CrystalOfTableaux(['A',2], shape=[2,1])
                sage: C = CrystalOfLetters(['A',2])
                sage: T = TensorProductOfCrystals(B,C)
                sage: B.is_connected()
                True
                sage: T.is_connected()
                False
            """
            return self.connected_components_number() == 1

    class ElementMethods:

        @cached_method
        def index_set(self):
            """
            EXAMPLES::

                sage: C = CrystalOfLetters(['A',5])
                sage: C(1).index_set()
                (1, 2, 3, 4, 5)
            """
            return self.parent().index_set()

        def cartan_type(self):
            """
            Returns the cartan type associated to ``self``

            EXAMPLES::

                sage: C = CrystalOfLetters(['A', 5])
                sage: C(1).cartan_type()
                ['A', 5]
            """
            return self.parent().cartan_type()

        @abstract_method
        def e(self, i):
            r"""
            Returns `e_i(x)` if it exists or ``None`` otherwise.

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
            Returns `f_i(x)` if it exists or ``None`` otherwise.

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

                sage: C = CrystalOfLetters(['A',5])
                sage: C(1).epsilon(1)
                0
                sage: C(2).epsilon(1)
                1
            """

        @abstract_method
        def phi(self, i):
            r"""
            EXAMPLES::

                sage: C = CrystalOfLetters(['A',5])
                sage: C(1).phi(1)
                1
                sage: C(2).phi(1)
                0
            """

        @abstract_method
        def weight(self):
            r"""
            Returns the weight of this crystal element

            This method should be implemented by the element class of
            the crystal.

            EXAMPLES::

                sage: C = CrystalOfLetters(['A',5])
                sage: C(1).weight()
                (1, 0, 0, 0, 0, 0)
            """

        def phi_minus_epsilon(self, i):
            """
            Returns `\phi_i - \epsilon_i` of self. There are sometimes
            better implementations using the weight for this. It is used
            for reflections along a string.

            EXAMPLES::

                sage: C = CrystalOfLetters(['A',5])
                sage: C(1).phi_minus_epsilon(1)
                1
            """
            return self.phi(i) - self.epsilon(i)

        def Epsilon(self):
            """
            EXAMPLES::

                sage: C = CrystalOfLetters(['A',5])
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

                sage: C = CrystalOfLetters(['A',5])
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
            Applies `f_{i_r} ... f_{i_1}` to self for `list = [i_1, ..., i_r]`

            EXAMPLES::

                sage: C = CrystalOfLetters(['A',3])
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
            Applies `e_{i_r} ... e_{i_1}` to self for `list = [i_1, ..., i_r]`

            EXAMPLES::

                sage: C = CrystalOfLetters(['A',3])
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
            Returns the reflection of ``self`` along its `i`-string

            EXAMPLES::

                sage: C = CrystalOfTableaux(['A',2], shape=[2,1])
                sage: b=C(rows=[[1,1],[3]])
                sage: b.s(1)
                [[2, 2], [3]]
                sage: b=C(rows=[[1,2],[3]])
                sage: b.s(2)
                [[1, 2], [3]]
                sage: T=CrystalOfTableaux(['A',2],shape=[4])
                sage: t=T(rows=[[1,2,2,2]])
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

                sage: C = CrystalOfLetters(['A',5])
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

                sage: C = CrystalOfLetters(['A',5])
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

                sage: T = CrystalOfTableaux(['A',3], shape = [1])
                sage: t = T(rows = [[3]])
                sage: t.to_highest_weight()
                [[[1]], [2, 1]]
                sage: T = CrystalOfTableaux(['A',3], shape = [2,1])
                sage: t = T(rows = [[1,2],[4]])
                sage: t.to_highest_weight()
                [[[1, 1], [2]], [1, 3, 2]]
                sage: t.to_highest_weight(index_set = [3])
                [[[1, 2], [3]], [3]]
                sage: K = KirillovReshetikhinCrystal(['A',3,1],2,1)
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

                sage: T = CrystalOfTableaux(['A',3], shape = [1])
                sage: t = T(rows = [[3]])
                sage: t.to_lowest_weight()
                [[[4]], [3]]
                sage: T = CrystalOfTableaux(['A',3], shape = [2,1])
                sage: t = T(rows = [[1,2],[4]])
                sage: t.to_lowest_weight()
                [[[3, 4], [4]], [1, 2, 2, 3]]
                sage: t.to_lowest_weight(index_set = [3])
                [[[1, 2], [4]], []]
                sage: K = KirillovReshetikhinCrystal(['A',3,1],2,1)
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

        def subcrystal(self, index_set=None, max_depth=float("inf"),
                       direction="both", contained=None, category=None):
            r"""
            Construct the subcrystal generated by ``self`` using `e_i` and/or
            `f_i` for all `i` in ``index_set``.

            INPUT:

            - ``index_set`` -- (Default: ``None``) The index set; if ``None``
              then use the index set of the crystal

            - ``max_depth`` -- (Default: infinity) The maximum depth to build

            - ``direction`` -- (Default: ``'both'``) The direction to build
              the subcrystal. It can be one of the following:

              - ``'both'`` - Using both `e_i` and `f_i`
              - ``'upper'`` - Using `e_i`
              - ``'lower'`` - Using `f_i`

            - ``contained`` -- (Optional) A function defining the containment
              in the subcrystal

            - ``category`` -- (Optional) Specify the category of the subcrystal

            .. SEEALSO::

                - :meth:`Crystals.ParentMethods.subcrystal()`.

            EXAMPLES::

                sage: C = KirillovReshetikhinCrystal(['A',3,1], 1, 2)
                sage: elt = C(1,4)
                sage: list(elt.subcrystal(index_set=[1,3]))
                [[[1, 4]], [[2, 4]], [[1, 3]], [[2, 3]]]
                sage: list(elt.subcrystal(index_set=[1,3], max_depth=1))
                [[[1, 4]], [[2, 4]], [[1, 3]]]
                sage: list(elt.subcrystal(index_set=[1,3], direction='upper'))
                [[[1, 4]], [[1, 3]]]
                sage: list(elt.subcrystal(index_set=[1,3], direction='lower'))
                [[[1, 4]], [[2, 4]]]
            """
            return self.parent().subcrystal(generators=[self], index_set=index_set,
                                            max_depth=max_depth, direction=direction,
                                            contained=contained, category=category)


#########################################################################
## Morphisms

class CrystalHomset(Homset):
    """
    The set of morphisms from one crystal to another.
    """
    def __init__(self, X, Y, category=None, cartan_type=None, index_set=None):
        """
        Initialize ``self``.
        """
        if category is None:
            category = Crystals()
        if cartan_type is None:
            cartan_type = X.cartan_type()
        if index_set is None:
            index_set = cartan_type.index_set()

        from sage.combinat.root_system.cartan_type import CartanType
        self._cartan_type = CartanType(cartan_type)
        self._index_set = tuple(index_set)
        Homset.__init__(self, X, Y, category)

    def _repr_(self):
        """
        TESTS::

            sage: B = CrystalOfTableaux(['A', 2], shape=[2,1])
            sage: Hom(B, B)
            sage: Hom(B, B, index_set=[1])
        """
        i_set = ''
        if self._index_set is not None:
            i_set = ' ' + repr(self._index_set)
        return "Set of{} {}-Crystal Morphisms from {} to {}".format(i_set,
            self.domain().cartan_type(), self.domain(), self.codomain())

    def _coerce_impl(self, x):
        """
        Check to see if we can coerce ``x`` into a morphism with the
        correct parameters.

        EXAMPLES::
        """
        if not isinstance(x, CrystalMorphism):
            raise TypeError

        if x.parent() is self:
            return x

        # Case 1: the parent fits
        if x.parent() == self:
            if isinstance(x, CrystalMorphismByGenerators):
                return CrystalMorphismByGenerators(self, x._on_gens, x._to_mg)

        if self.cartan_type() != x.cartan_type():
            raise ValueError

        # Case 2: it has a smaller index set
        if set(self._index_set).is_subset(set(x.parent()._index_set)):
            if isinstance(x, HighestWeightCrystalMorphism):
                return HighestWeightCrystalMorphism(self, x._on_gens, self._index_set)
            elif isinstance(x, CrystalMorphismByGenerators):
                return CrystalMorphismByGenerators(self, x._on_gens, x._to_mg)

        raise ValueError

    def __call__(self, on_gens=None, to_module_generator=None, check=True):
        """
        Construct a crystal morphism.
        """
        if isinstance(on_gens, CrystalMorphism):
            return self._coerce_impl(on_gens)

        if on_gens is None:
            # Match generators by weight
            cod_wt = dict((x.weight(), x) for x in self.codomain().module_generators)
            if len(cod_wt) != self.codomain().module_generators:
                raise ValueError("if no map is specified,"
                                 " the highest weights of the codomain's"
                                 " module generators must be unique")
            on_gens = lambda x: cod_wt.get(x.weight(), None)

        from sage.categories.highest_weight_crystals import HighestWeightCrystals
        if self.domain() in HighestWeightCrystals():
            from sage.categories.highest_weight_crystals import HighestWeightCrystalMorphism
            return HighestWeightCrystalMophism(self, on_gens, self._index_set)
        return CrystalMorphismByGenerators(self, on_gens, to_module_generator)

    def cartan_type(self):
        """
        Return the Cartan type of this Hom space.
        """
        return self._cartan_type

    def index_set(self):
        """
        Return the index set of this Hom space.
        """
        return self._index_set

    # TODO: Combine with the crystal_morphism in ParentMethods
    def crystal_morphism(self, domain, on_gens, codomain=None, index_set=None):
        if codomain is None:
            if hasattr(on_gens, 'codomain'):
                codomain = on_gens.codomain()
            else:
                for x in self.module_generators:
                    y = on_gens(x)
                    if y is not None:
                        codomain = y.parent()
                        break
        if codomain is None:
            raise ValueError("no codomain specified")

        if category is None:
            if codomain not in Crystals():
                raise ValueError("the codomain must be a crystal")
            category = Crystals()
        parent = Hom(domain, codomain, category=category, index_set=index_set)
        return CrystalMorphismByGenerators(parent, on_gens)

class CrystalMorphism(Morphism):
    """
    A crystal morphism.
    """
    def _repr_type(self):
        """
        Used internally in printing this morphism.

        TESTS::

            sage: B = CrystalOfTableaux(['A',2], shape=[2,1])
            sage: phi = B.hom(B)
            sage: phi._repr_type()
            "['A', 2]-Crystal"
        """
        return "{}-Crystal".format(self.parent()._cartan_type)

class CrystalMorphismByGenerators(CrystalMorphism):
    """
    A crystal morphism defined by a set of generators.

    INPUT:

    - ``parent`` -- a Hom category
    - ``on_gens`` -- a function for the images of the generators
    - ``gens`` -- a list of generators to define the morphism; the
      default is to use the module generators of the crystal

    If there is a fast way to find the module generators, overwrite
    :meth:`to_module_generator()` which currently looks for a module generator
    by a breath first search.
    """
    def __init__(self, parent, on_gens, gens=None):
        """
        Construct a crystal morphism.

        EXAMPLES::
        """
        if gens is None:
            gens = dom.module_generators
        self._gens = gens

        # Make sure on_gens is a function
        if isinstance(on_gens, dict):
            on_gens = lambda x: on_gens[x]
        elif isinstance(on_gens, (list, tuple)):
            dom = parent.domain()
            if len(self._gens) != len(on_gens):
                raise ValueError("not all module generators have an image")
            on_gens = lambda x: on_gens[self._gens.index(x)]
        self._on_gens = on_gens
        self._path_mg_cache = {x: (x, [], []) for x in self._gens}

        CrystalMorphism.__init__(self, parent, mor_type)

    def _call_(self, x):
        """
        Return the image of ``x`` under ``self``.

        EXAMPLES::
        """
        if x is None:
            return None
        mg, ef, indices = self.to_module_generator(x)
        cur = self._on_gens(mg)
        for op,i in zip(ef, indices):
            if cur is None:
                return None
            cur = getattr(cur, op)(i)
        return cur

    def to_module_generator(self, x):
        """
        Return a generator ``mg`` and a path of `e_i` and `f_i` operations
        to ``mg``.

        OUTPUT:

        A tuple:

        - a module generator,
        - a list of 'e' and 'f' to denote which operation,
        - a list of matching indices.
        """
        mg = set(self._path_mg_cache.keys())
        if x in mg:
            return mg[x]

        visited = set([None, x])
        index_set = self.parent().index_set()
        ef = [[]]
        indices = [[]]
        todo = [x]

        while len(todo) > 0:
            for i in index_set:
                cur = todo.pop(0)
                cur_ef = ef.pop(0)
                cur_indices = indices.pop(0)

                next = cur.e(i)
                if next in mg:
                    gen,ef,indices = self._path_mg_cache[next]
                    ef = cur_ef + ['e'] + ef
                    indices = cur_indices + [i] + indices
                    self._path_mg_cache[x] = (gen, ef, indices)
                    return (gen, ef, indices)
                if next not in visited:
                    ef.append(cur_ef + ['e'])
                    indices.append(cur_indices + [i])
                    visited.add(next)
                    todo.append(next)

                # Now for f's
                next = cur.f(i)
                if next in mg:
                    gen,ef,indices = self._path_mg_cache[next]
                    ef = cur_ef + ['f'] + ef
                    indices = cur_indices + [i] + indices
                    self._path_mg_cache[x] = (gen, ef, indices)
                    return (gen, ef, indices)
                if next not in visited:
                    ef.append(cur_ef + ['f'])
                    indices.append(cur_indices + [i])
                    visited.add(next)
                    todo.append(next)
        raise ValueError("no module generator in the component of {}".format(x))

    def _repr_defn(self):
        """
        Used in constructing string representation of ``self``.

        EXAMPLES::

            sage: B = CrystalOfTableaux(['A',2], shape=[2,1])
            sage: C = CrystalOfTableaux(['A',2], shape=[4,2])
            sage: f = B.hom(B, [C(2,1,3,1,1,2)])
            sage: print f._repr_defn()
            x |--> x^2
            y |--> x + y
        """
        D = self.domain()
        return '\n'.join(['{} |--> {}'.format(mg, im)
                          for mg,im in zip(self._gens, self.im_gens())])

    @cached_method
    def im_gens(self):
        """
        Return the image of the generators of ``self`` as a tuple.
        """
        return tuple(map(self._on_gens, self._gens))

class TwistedCrystalMorphismByGenerators(CrystalMorphismByGenerators):
    """
    A crystal morphism defined by a set of generators twisted by
    an automorphism.

    INPUT:

    - ``parent`` -- a Hom category
    - ``on_gens`` -- a function for the images of the generators
    - ``automorphism`` -- the twisting automorphism
    - ``gens`` -- a list of generators to define the morphism; the
      default is to use the module generators of the crystal
    """
    def __init__(self, parent, on_gens, automorphism, gens=None):
        """
        Construct a twisted crystal morphism.

        EXAMPLES::
        """
        self._twist = automorphism
        CrystalMorphismByGenerators.__init__(self, parent, on_gens, to_module_generator, gens)

    def _call_(self, x):
        """
        Return the image of ``x`` under ``self``.

        EXAMPLES::
        """
        if x is None:
            return None
        mg, ef, indices = self._to_mg(x)
        cur = self._on_gens(mg)
        for op,i in zip(ef, indices):
            if cur is None:
                return None
            cur = getattr(cur, op)(self._twist(i))
        return cur

class VirtualCrystalMorphismByGenerators(CrystalMorphismByGenerators):
    """
    A crystal morphism defined by a set of generators which create a virtual
    crystal inside the codomain.

    INPUT:

    - ``parent`` -- a Hom category
    - ``on_gens`` -- a function for the images of the generators
    - ``virtualization`` -- a dictionary whose keys are in the index set of
      the domain and whose values are lists of entries in the index set of the
      codomain
    - ``gens`` -- a list of generators to define the morphism; the
      default is to use the module generators of the crystal
    """
    def __init__(self, parent, on_gens, virtualization, gens=None):
        """
        Construct a virtual crystal morphism.

        EXAMPLES::
        """
        self._virtualization = virtualization
        CrystalMorphismByGenerators.__init__(self, parent, on_gens, to_module_generator, gens)

    def _call_(self, x):
        """
        Return the image of ``x`` under ``self``.

        EXAMPLES::
        """
        if x is None:
            return None
        mg, ef, indices = self._to_mg(x)
        cur = self._on_gens(mg)
        for op,i in zip(ef, indices):
            if cur is None:
                return None
            cur = getattr(cur, op + "_string")(self._virtualization[i])
        return cur

