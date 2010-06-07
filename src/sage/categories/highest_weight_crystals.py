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
from sage.categories.crystals import Crystals

class HighestWeightCrystals(Category):
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
            generators those that are highest weight, and cache the
            result.

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

            Caveat: this assumes that :meth:`.highest_weight_vector`
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


    class ElementMethods:

        def to_highest_weight(self, list = [], index_set = None):
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
                sage: t.to_highest_weight()
                [[[1]], [2, 1]]
                sage: T = CrystalOfTableaux(['A',3], shape = [2,1])
                sage: t = T(rows = [[1,2],[4]])
                sage: t.to_highest_weight()
                [[[1, 1], [2]], [1, 3, 2]]
                sage: t.to_highest_weight(index_set = [3])
                [[[1, 2], [3]], [3]]
            """
            if index_set is None:
                index_set = self.index_set()
            for i in index_set:
                if self.epsilon(i) <> 0:
                    self = self.e(i)
                    return self.to_highest_weight(list = list + [i], index_set = index_set)
            return [self, list]
