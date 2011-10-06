r"""
Highest Weight Crystals
"""
#*****************************************************************************
#  Copyright (C) 2010    Anne Schilling <anne at math.ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from sage.categories.category import Category
from sage.categories.category_singleton import Category_singleton
from sage.categories.crystals import Crystals

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
        running ._test_elements_eq() . . . pass
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

        @cached_method
        def lowest_weight_vectors(self):
            r"""
            Returns the lowest weight vectors of ``self``

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

    class ElementMethods:

        pass
