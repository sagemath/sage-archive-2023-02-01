"""
Highest Weight Crystals
"""
#*****************************************************************************
#  Copyright (C) 2010    Anne Schilling <anne at math.ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from sage.categories.category_singleton import Category_singleton
from sage.categories.crystals import Crystals, CrystalHomset, \
    CrystalMorphismByGenerators, TwistedCrystalMorphismByGenerators, \
    VirtualCrystalMorphismByGenerators

class HighestWeightCrystals(Category_singleton):
    """
    The category of highest weight crystals.

    A crystal is highest weight if it is acyclic; in particular, every
    connected component has a unique highest weight element, and that
    element generate the component.

    EXAMPLES::

        sage: C = HighestWeightCrystals()
        sage: C
        Category of highest weight crystals
        sage: C.super_categories()
        [Category of crystals]
        sage: C.example()
        Highest weight crystal of type A_3 of highest weight omega_1

    TESTS::

        sage: TestSuite(C).run()
        sage: B = HighestWeightCrystals().example()
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

            sage: HighestWeightCrystals().super_categories()
            [Category of crystals]
        """
        return [Crystals()]

    def example(self):
        """
        Returns an example of highest weight crystals, as per
        :meth:`Category.example`.

        EXAMPLES::

            sage: B = HighestWeightCrystals().example(); B
            Highest weight crystal of type A_3 of highest weight omega_1
        """
        from sage.categories.crystals import Crystals
        return Crystals().example()

    class ParentMethods:

        @cached_method
        def highest_weight_vectors(self):
            r"""
            Returns the highest weight vectors of ``self``

            This default implementation selects among the module
            generators those that are highest weight, and caches the result.
            A crystal element `b` is highest weight if `e_i(b)=0` for all `i` in the
            index set.


            EXAMPLES::

                sage: C = CrystalOfLetters(['A',5])
                sage: C.highest_weight_vectors()
                [1]

            ::

                sage: C = CrystalOfLetters(['A',2])
                sage: T = TensorProductOfCrystals(C,C,C,generators=[[C(2),C(1),C(1)],[C(1),C(2),C(1)]])
                sage: T.highest_weight_vectors()
                [[2, 1, 1], [1, 2, 1]]
            """
            return [g for g in self.module_generators if g.is_highest_weight()]

        def highest_weight_vector(self):
            r"""
            Returns the highest weight vector if there is a single one;
            otherwise, raises an error.

            Caveat: this assumes that :meth:`.highest_weight_vectors`
            returns a list or tuple.

            EXAMPLES::

                sage: C = CrystalOfLetters(['A',5])
                sage: C.highest_weight_vector()
                1
            """
            hw = self.highest_weight_vectors();
            if len(hw) == 1:
                return hw[0]
            else:
                raise RuntimeError("The crystal does not have exactly one highest weight vector")

        # TODO: Not every highest weight crystal is a lowest weight crystal
        @cached_method
        def lowest_weight_vectors(self):
            r"""
            Return the lowest weight vectors of ``self``.

            This default implementation selects among all elements of the crystal
            those that are lowest weight, and cache the result.
            A crystal element `b` is lowest weight if `f_i(b)=0` for all `i` in the
            index set.

            EXAMPLES::

                sage: C = CrystalOfLetters(['A',5])
                sage: C.lowest_weight_vectors()
                [6]

            ::

                sage: C = CrystalOfLetters(['A',2])
                sage: T = TensorProductOfCrystals(C,C,C,generators=[[C(2),C(1),C(1)],[C(1),C(2),C(1)]])
                sage: T.lowest_weight_vectors()
                [[3, 2, 3], [3, 3, 2]]
            """
            return [g for g in self if g.is_lowest_weight()]

        def __iter__(self, index_set=None, max_depth = float("inf")):
            """
            Returns the iterator of ``self``.

            INPUT:

            - ``index_set`` -- (Default: ``None``) The index set; if ``None``
              then use the index set of the crystal

            - ``max_depth`` -- (Default: infinity) The maximum depth to build

            EXAMPLES::

                sage: C = CrystalOfLSPaths(['A',2,1],[0,1,0])
                sage: [p for p in C.__iter__(max_depth=3)]
                [(Lambda[1],), (Lambda[0] - Lambda[1] + Lambda[2],), (2*Lambda[0] - Lambda[2],),
                (-Lambda[0] + 2*Lambda[2] - delta,),
                (1/2*Lambda[0] + Lambda[1] - Lambda[2] - 1/2*delta, -1/2*Lambda[0] + Lambda[2] - 1/2*delta),
                (-Lambda[0] + Lambda[1] + 1/2*Lambda[2] - delta, Lambda[0] - 1/2*Lambda[2])]
                sage: [p for p in C.__iter__(index_set=[0, 1], max_depth=3)]
                [(Lambda[1],), (Lambda[0] - Lambda[1] + Lambda[2],), (-Lambda[0] + 2*Lambda[2] - delta,)]
            """
            if index_set is None:
                index_set = self.index_set()
            from sage.combinat.backtrack import TransitiveIdealGraded
            return TransitiveIdealGraded(lambda x: [x.f(i) for i in index_set],
                                         self.module_generators, max_depth).__iter__()

        # TODO: This is not correct if a factor has multiple heads (i.e., we
        #   should have a category for uniqueness of highest/lowest weights)
        connected_components_generators = highest_weight_vectors

        def _Hom_(self, Y, category=None, **options):
            r"""
            Return the homset from ``self`` to ``Y`` in the
            category ``category``.

            INPUT::

            - ``Y`` -- a crystal
            - ``category`` -- a subcategory of :class:`HighestWeightCrysals`()
              or ``None``

            The sole purpose of this method is to construct the homset as a
            :class:`~sage.categories.highest_weight_crystals.HighestWeightCrystalHomset`.
            If ``category`` is specified and is not a subcategory of
            :class:`HighestWeightCrystals`, a ``TypeError`` is raised instead

            This method is not meant to be called directly. Please use
            :func:`sage.categories.homset.Hom` instead.

            EXAMPLES::

                sage: B = CrystalOfTableaux(['A',2], shape=[2,1])
                sage: H = B._Hom_(B)
                sage: H
                Set of Crystal Morphisms from The crystal of tableaux of type ['A', 2] and shape(s) [[2, 1]]
                 to The crystal of tableaux of type ['A', 2] and shape(s) [[2, 1]]
                sage: type(H)
                <class 'sage.categories.highest_weight_crystals.HighestWeightCrystalHomset_with_category'>
            """
            if category is None:
                category = self.category()
            elif not category.is_subcategory(HighestWeightCrystals()):
                raise TypeError("{} is not a subcategory of HighestWeightCrystals()".format(category))
            if Y not in Crystals():
                raise TypeError("{} is not a crystal".format(Y))
            return HighestWeightCrystalHomset(self, Y, category=category, **options)

    class ElementMethods:
        pass

###############################################################################
## Morphisms

class HighestWeightCrystalMorphism(CrystalMorphismByGenerators):
    """
    A crystal morphism whose domain is a highest weight crystal.

    INPUT:

    - ``parent`` -- a homset
    - ``on_gens`` -- a function for the images of the module generators
    - ``cartan_type`` -- (optional) the Cartan type of the morphism; the
      default is the Cartan type of the domain
    - ``index_set`` -- (optional) the index set of the morphism; the
      default is the index set of the Cartan type
    - ``gens`` -- a list of generators to define the morphism; the
      default is to use the highest weight vectors of the crystal
    - ``check`` -- (default: ``True``) check if the crystal morphism is valid
    """
    def __init__(self, parent, on_gens, cartan_type=None, index_set=None, gens=None, check=True):
        """
        Initialize ``self``.

        TESTS::

            sage: B = InfinityCrystalOfTableaux(['B',2])
            sage: C = InfinityCrystalOfNakajimaMonomials(['B',2])
            sage: psi = B.crystal_morphism(C.module_generators)
        """
        if gens is None:
            gens = parent.domain().highest_weight_vectors()
        return CrystalMorphismByGenerators.__init__(self, parent, on_gens, cartan_type,
                                                    index_set, gens, check)

    def _call_(self, x):
        """
        Return the image of ``x`` under ``self``.

        EXAMPLES::

            sage: B = InfinityCrystalOfTableaux(['B',2])
            sage: C = InfinityCrystalOfNakajimaMonomials(['B',2])
            sage: psi = B.crystal_morphism(C.module_generators)
            sage: b = B.highest_weight_vector()
            sage: psi(b)
            1
            sage: c = psi(b.f_string([1,1,1,2,2,1,2,2])); c
            Y(1,0)^-4 Y(2,0)^4 Y(2,1)^-4 
            sage: c == C.highest_weight_vector().f_string([1,1,1,2,2,1,2,2])
            True
        """
        mg, path = x.to_highest_weight(self._index_set)
        y = self._on_gens(mg)
        if y is None:
            return None
        return y.f_string(reversed(path))

class HighestWeightTwistedCrystalMorphism(TwistedCrystalMorphismByGenerators):
    """
    A twisted crystal morphism whose domain is a highest weight crystal.

    INPUT:

    - ``parent`` -- a homset
    - ``on_gens`` -- a function for the images of the generators
    - ``automorphism`` -- the twisting automorphism
    - ``cartan_type`` -- (optional) the Cartan type of the morphism; the
      default is the Cartan type of the domain
    - ``index_set`` -- (optional) the index set of the morphism; the
      default is the index set of the Cartan type
    - ``gens`` -- a list of generators to define the morphism; the
      default is to use the highest weight vectors of the crystal
    - ``check`` -- (default: ``True``) check if the crystal morphism is valid
    """
    def __init__(self, parent, on_gens, automorphism, cartan_type=None,
                 index_set=None, gens=None, check=True):
        """
        Initialize ``self``.

        TESTS::

            sage: B = CrystalOfTableaux(['D',4], shape=[1])
            sage: H = Hom(B, B)
            sage: d = {1:1, 2:2, 3:4, 4:3}
            sage: psi = H(B.module_generators, automorphism=lambda x: d[x])
        """
        if gens is None:
            gens = parent.domain().highest_weight_vectors()
        return TwistedCrystalMorphismByGenerators.__init__(self, parent, on_gens, automorphism,
                                                           cartan_type, index_set, gens, check)

    def _call_(self, x):
        """
        Return the image of ``x`` under ``self``.

        TESTS::

            sage: B = CrystalOfTableaux(['D',4], shape=[1])
            sage: H = Hom(B, B)
            sage: d = {1:1, 2:2, 3:4, 4:3}
            sage: psi = H(B.module_generators, automorphism=lambda x: d[x])
            sage: b = B.highest_weight_vector()
            sage: psi(b.f_string([1,2,3]))
            [[-4]]
        """
        mg, path = x.to_highest_weight(self._index_set)
        y = self._on_gens(mg)
        if y is None:
            return None
        return y.f_string(map(self._twist, reversed(path)))

class HighestWeightVirtualCrystalMorphism(VirtualCrystalMorphismByGenerators):
    r"""
    A virtual crystal morphism whose domain is a highest weight crystal.

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
      default is to use the highest weight vectors of the crystal
    - ``check`` -- (default: ``True``) check if the crystal morphism is valid
    """
    def __init__(self, parent, on_gens, virtualization, scaling_factors,
                 cartan_type=None, index_set=None, gens=None, check=True):
        """
        Construct a virtual crystal morphism.

        TESTS::

            sage: B = CrystalOfTableaux(['B',3], shape=[1])
            sage: C = CrystalOfTableaux(['D',4], shape=[2])
            sage: H = Hom(B, C)
            sage: psi = H(C.module_generators)
        """
        if gens is None:
            gens = parent.domain().highest_weight_vectors()
        VirtualCrystalMorphismByGenerators.__init__(self, parent, on_gens,
                                                    virtualization, scaling_factors,
                                                    cartan_type, index_set, gens, check)

    def _call_(self, x):
        """
        Return the image of ``x`` under ``self``.

        TESTS::

            sage: B = CrystalOfTableaux(['B',3], shape=[1])
            sage: C = CrystalOfTableaux(['D',4], shape=[2])
            sage: H = Hom(B, C)
            sage: psi = H(C.module_generators)
            sage: psi(B.module_generators[0])
            [[1, 1]]
        """
        mg, path = x.to_highest_weight(self._index_set)
        cur = self._on_gens(mg)
        for i in reversed(path):
            if cur is None:
                return None
            s = []
            sf = self._scaling_factors[i]
            for j in self._virtualization[i]:
                s += [j]*sf
            cur = cur.f_string(s)
        return cur

class HighestWeightCrystalHomset(CrystalHomset):
    """
    The set of crystal morphisms from a highest weight crystal to
    another crystal.

    .. SEEALSO::

        See :class:`sage.categories.crystals.CrystalHomset` for more
        information.
    """
    def __init__(self, X, Y, category=None):
        """
        Initialize ``self``.

        TESTS::

            sage: B = CrystalOfTableaux(['A', 2], shape=[2,1])
            sage: H = Hom(B, B)
            sage: B = InfinityCrystalOfTableaux(['B',2])
            sage: H = Hom(B, B)
        """
        if category is None:
            category = HighestWeightCrystals()
        CrystalHomset.__init__(self, X, Y, category)

    _generic = HighestWeightCrystalMorphism
    _twisted = HighestWeightTwistedCrystalMorphism
    _virtual = HighestWeightVirtualCrystalMorphism
