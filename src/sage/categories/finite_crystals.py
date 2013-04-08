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

        def _test_stembridge_local_axioms(self, index_set=None, verbose=False, complete=False, **options):
            r"""
            This implements tests for the Stembridge local characterization on the finite crystal ``self``.
            The current implementation only uses the rules for simply-laced types.  Crystals
            of other types should still pass the test, but expansion of this test to
            non-simply laced type would be desirable.

            One can specify an index set smaller than the full index set of the crystal,
            using the option ``index_set``.

            Running with ``verbose=True`` will print each node for which a local axiom
            test applies.

            Running with ``complete=True`` will continue to run the test past the first
            failure of the local axioms.  This is probably only useful in conjunction
            with the verbose option, to see all places where the local axioms fail.

            EXAMPLES::

                sage: T = CrystalOfTableaux(['A',3], shape=[2,1])
                sage: T._test_stembridge_local_axioms()
                True
                sage: T._test_stembridge_local_axioms(verbose=True)
                True
                sage: T._test_stembridge_local_axioms(index_set=[1,3])
                True

                sage: B=Crystals().example(choice='naive')
                sage: B._test_stembridge_local_axioms()
                Traceback (most recent call last):
                ...
                AssertionError: None
            """
            tester = self._tester(**options)
            goodness=True

            for x in self:
                goodness=x._test_stembridge_local_axioms(index_set, verbose)
                if goodness==False and not complete:
                    tester.fail()
            tester.assertTrue(goodness)
            return goodness
