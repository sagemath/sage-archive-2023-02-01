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
from sage.categories.crystals import Crystals, CrystalMorphismByGenerators
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
            @cached_method
            def is_injective(self):
                """
                Return if ``self`` is an injective crystal isomorphism by a
                brute force check.
                """
                image = set([])
                for x in self.domain:
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
                """
                index_set = self.domain
                for x in self.domain:
                    for i in index_set:
                        if self(x.e(i)) != self(x).e(i) or self(x.f(i)) != self(x).f(i):
                            return False
                return True

            def is_isomorphism(self):
                """
                Check if ``self`` is a crystal isomorphism, which is true
                if and only if this is a strict embedding.
                """
                return self.domain.is_isomorphic(self.codomain)

    class ParentMethods:
        @cached_method
        def connected_components_generators(self):
            """
            Return a tuple of generators for each of the connected components
            of ``self``.
            """
            return tuple(map(lambda x: x[0], self.digraph().connected_components()))

        def is_isomorphic(self, X):
            """
            Check if ``self`` is isomorphic to another crystal ``X`` by
            brute force check of their digraphs.
            """
            from sage.rings.infinity import infinity
            if X.cardinality() == infinity:
                return False
            return self.digraph().is_ismorphic(X.digraph(), edge_labels=True)

