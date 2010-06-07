r"""
Example of a crystal
"""
#*****************************************************************************
#  Copyright (C) 2010 Anne Schilling <anne at math.ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.structure.parent import Parent
from sage.structure.element_wrapper import ElementWrapper
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.classical_crystals import ClassicalCrystals
from sage.categories.enumerated_sets import EnumeratedSets
from sage.combinat.root_system.cartan_type import CartanType

class HighestWeightCrystalOfTypeA(UniqueRepresentation, Parent):
    r"""
    An example of a crystal: the highest weight crystal of type `A_n`
    of highest weight `\omega_1`.

    The purpose of this class is to provide a minimal template for
    implementing crystals. See
    :class:`~sage.combinat.crystals.letters.CrystalOfLetters` for a
    full featured and optimized implementation.

    EXAMPLES::

        sage: C = Crystals().example()
        sage: C
        Highest weight crystal of type A_3 of highest weight omega_1
        sage: C.category()
        Category of classical crystals

    The elements of this crystal are in the set `\{1,\ldots,n+1\}`::

        sage: C.list()
        [1, 2, 3,  4]
        sage: C.module_generators[0]
        1

    The crystal operators themselves correspond to the elementary
    transpositions::

        sage: b = C.module_generators[0]
        sage: b.f(1)
        2
        sage: b.f(1).e(1) == b
        True

    TESTS::

        sage: TestSuite(C).run(verbose = True)
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

    Only the following basic operations are implemented:
     - :meth:`.cartan_type` or an attribute _cartan_type
     - an attribute `module_generators`
     - :meth:`.Element.e`
     - :meth:`.Element.f`
     - :meth:`.Element.weight`

    All the other usual crystal operations are inherited from the
    categories; for example::

        sage: C.cardinality()
        4
    """

    def __init__(self, n = 3):
        """
        EXAMPLES::

            sage: C = sage.categories.examples.crystals.HighestWeightCrystalOfTypeA(4)
            sage: C == Crystals().example(4)
            True
        """
        Parent.__init__(self, category = ClassicalCrystals())
        self.n = n
        self._cartan_type = CartanType(['A',n])
        self.module_generators = [ self(1) ]

    def _repr_(self):
        """
        EXAMPLES::

            sage: Crystals().example()
            Highest weight crystal of type A_3 of highest weight omega_1
        """
        return "Highest weight crystal of type A_%s of highest weight omega_1"%(self.n)

    # temporary woraround while an_element is overriden by Parent
    _an_element_ = EnumeratedSets.ParentMethods._an_element_

    class Element(ElementWrapper):

        def e(self, i):
            r"""
            Returns the action of `e_i` on ``self``.

            EXAMPLES::

                sage: C = Crystals().example(4)
                sage: [[c,i,c.e(i)] for i in C.index_set() for c in C if c.e(i) is not None]
                [[2, 1, 1], [3, 2, 2], [4, 3, 3], [5, 4, 4]]
            """
            assert i in self.index_set()
            if self.value == i+1:
                return self.parent()(self.value-1)
            else:
                return None

        def f(self, i):
            r"""
            Returns the action of `f_i` on ``self``.

            EXAMPLES::

                sage: C = Crystals().example(4)
                sage: [[c,i,c.f(i)] for i in C.index_set() for c in C if c.f(i) is not None]
                [[1, 1, 2], [2, 2, 3], [3, 3, 4], [4, 4, 5]]
            """
            assert i in self.index_set()
            if self.value == i:
                return self.parent()(self.value+1)
            else:
                return None

Example = HighestWeightCrystalOfTypeA
