r"""
Finite Crystals
"""
#*****************************************************************************
#  Copyright (C) 2010    Anne Schilling <anne at math.ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from sage.categories.category import HomCategory
from sage.categories.category import Category
from sage.categories.crystals import Crystals
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets

class FiniteCrystals(Category):
    """
    The category of finite crystals.

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

            sage: FiniteCrystals().super_categories()
            [Category of crystals, Category of finite enumerated sets]
        """
        return [Crystals(), FiniteEnumeratedSets()]

    def example(self, n = 3):
        """
        Returns an example of highest weight crystals, as per
        :meth:`Category.example`.

        EXAMPLES::

            sage: B = FiniteCrystals().example(); B
            Highest weight crystal of type A_3 of highest weight omega_1
        """
        from sage.categories.crystals import Crystals
        return Crystals().example(n)

    class HomCategory(HomCategory):
        """
        The category of homomorphisms sets `Hom(X,Y)` for a finite crystal `X`
        and an arbitrary crystal `Y`.
        """
        class ElementMethods:
            # Should be is_injective, but that is defined in the class itself
            @cached_method
            def is_embedding(self):
                """
                Return if ``self`` is an injective crystal isomorphism by a
                brute force check.

                EXAMPLES::

                    sage: B = crystals.Tableaux(['A',2], shape=[2,1])
                    sage: C = crystals.Tableaux(['A',2], ([2,1], [1,1]))
                    sage: psi = B.crystal_morphism(C.module_generators[:1], codomain=C)
                    sage: psi.is_embedding()
                    True
                    sage: psi = C.crystal_morphism([B.module_generators[0], None], codomain=B)
                    sage: psi.is_embedding()
                    False
                """
                image = set([None])
                for x in self.domain():
                    y = self(x)
                    if y in image:
                        return False
                    image.add(y)
                return True

            @cached_method
            def is_strict(self):
                """
                Return if ``self`` is a strict crystal morphism by a
                brute force check.

                EXAMPLES::

                    sage: B = crystals.Tableaux(['A',2], shape=[2,1])
                    sage: C = crystals.Tableaux(['A',2], ([2,1], [1,1]))
                    sage: psi = B.crystal_morphism(C.module_generators[:1], codomain=C)
                    sage: psi.is_strict()
                    True
                    sage: psi = C.crystal_morphism([B.module_generators[0], None], codomain=B)
                    sage: psi.is_strict()
                    False
                """
                for x in self.domain():
                    for i in self._index_set:
                        y = self(x)
                        if y is None:
                            return False
                        if self(x.e(i)) != self(x).e(i) or self(x.f(i)) != self(x).f(i):
                            return False
                return True

            def is_isomorphism(self):
                """
                Check if ``self`` is a crystal isomorphism, which is true
                if and only if this is a strict embedding.

                EXAMPLES::

                    sage: B = crystals.Tableaux(['A',2], shape=[1,1])
                    sage: C = crystals.Tableaux(['A',2], ([2,1], [1,1]))
                    sage: psi = B.crystal_morphism(C.module_generators[1:], codomain=C)
                    sage: psi.is_isomorphism()
                    False
                    sage: K = crystals.KirillovReshetikhin(['A',2,1], 2,1)
                    sage: psi = K.crystal_morphism(B.module_generators, codomain=B, cartan_type=['A',2])
                    sage: psi.is_isomorphism()
                    True
                """
                return self.domain().is_isomorphic(self.codomain(), self.index_set())

            def image(self):
                """
                Return the image of ``self`` in the codomain as a subcrystal.

                EXAMPLES::

                    sage: B = crystals.Tableaux(['D',4], shape=[])
                    sage: K = crystals.KirillovReshetikhin(['D',4,1], 2,1, model='KR')
                    sage: psi = B.crystal_morphism(K.module_generators[:1], codomain=K, category=FiniteCrystals())
                    sage: S = psi.image()
                    sage: S
                    Subcrystal of Kirillov-Reshetikhin tableaux of type ['D', 4, 1] and shape (2, 1)
                    sage: list(S)
                    [[[1], [-1]]]
                """
                image = set([])
                for x in self.domain():
                    y = self(x)
                    if y is not None:
                        image.add(y)
                gens = filter(lambda x: x is not None, self.im_gens())
                from sage.combinat.crystals.subcrystal import Subcrystal
                return Subcrystal(self.codomain(), image, gens, cartan_type=self.cartan_type(),
                                  index_set=self.index_set(), category=self.domain().category())

    class ParentMethods:
        @cached_method
        def connected_components_generators(self):
            """
            Return a tuple of generators for each of the connected components
            of ``self``.

            EXAMPLES::

                sage: K = crystals.KirillovReshetikhin(['D',4,1], 2,1)
                sage: K.connected_components_generators()
                ([[1], [2]],)
            """
            return tuple(map(lambda x: x[0], self.digraph().connected_components()))

        def is_isomorphic(self, X, index_set=None):
            """
            Check if ``self`` is isomorphic as `I`-crystals to another
            crystal ``X`` by a brute force check of their digraphs.

            INPUT:

            - ``X`` -- another crystal
            - ``index_set`` -- (optional) the set `I`; the default is to use
              the index set of ``self``.

            EXAMPLES::

                sage: K = crystals.KirillovReshetikhin(['D',4,1], 2,1)
                sage: B = K.classical_decomposition()
                sage: B.is_isomorphic(K)
                False
                sage: K.is_isomorphic(B)
                False
                sage: K.is_isomorphic(B, index_set=[1,2,3,4])
                True
                sage: B.is_isomorphic(K, index_set=[1,2,3,4])
                True
            """
            from sage.rings.infinity import infinity
            if X.cardinality() == infinity:
                return False
            G = self.digraph(index_set=index_set)
            GX = X.digraph(index_set=index_set)
            return G.is_isomorphic(GX, edge_labels=True)

