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

    class ParentMethods:

        def list(self):
            """
            Returns a list of the elements of ``self`` obtained by
            repeatedly applying the `f_i` operators to the module
            generators of ``self``.

            EXAMPLES::

                sage: C = FiniteCrystals().example(5)
                sage: l = C._list_brute_force()
                sage: l.sort(); l
                [1, 2, 3, 4, 5, 6]
            """
            # Should use transitiveIdeal
            # should be transformed to __iter__ instead of list
            # To be moved in a super category CombinatorialModule
            result = set(self.module_generators)
            todo = result.copy()
            while len(todo) > 0:
                x = todo.pop()
                for i in self.index_set():
                    y = x.f(i)
                    if y == None or y in result:
                        continue
                    todo.add(y)
                    result.add(y)
            return list(result)

        _list_brute_force = list
