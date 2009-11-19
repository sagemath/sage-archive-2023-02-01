"""
Finite Enumerated Sets
"""
#*****************************************************************************
#  Copyright (C) 2009 Florent Hivert <Florent.Hivert@univ-rouen.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.rings.integer import Integer

#################################################################
class FiniteEnumeratedSet(UniqueRepresentation, Parent):
    """
    A class for finite enumerated set.

    Returns the finite enumerated set with elements in ``elements``
    where ``element`` is any (finite) iterable object.

    The main purpose is to provide a variant of ``list`` or ``tuple``,
    which is a parent with an interface consistent with
    ``EnumeratedSets`` and has unique representation.
    The list of the elements is expanded in memory.


    EXAMPLES::

        sage: S = FiniteEnumeratedSet([1, 2, 3])
        sage: S
        {1, 2, 3}
        sage: S.list()
        [1, 2, 3]
        sage: S.cardinality()
        3
        sage: S.random_element()
        1
        sage: S.first()
        1
        sage: S.category()
        Category of finite enumerated sets
        sage: TestSuite(S).run(verbose = True)
        running ._test_an_element() ... done
        running ._test_element_pickling() ... done
        running ._test_enumerated_set_contains() ... done
        running ._test_enumerated_set_iter_cardinality() ... done
        running ._test_enumerated_set_iter_list() ... done
        running ._test_not_implemented_methods() ... done
        running ._test_pickling() ... done
        running ._test_some_elements() ... done

    Note that being and enumerated set, the result depends on the order::

        sage: S1 = FiniteEnumeratedSet((1, 2, 3))
        sage: S1
        {1, 2, 3}
        sage: S1.list()
        [1, 2, 3]
        sage: S1 == S
        True
        sage: S2 = FiniteEnumeratedSet((2, 1, 3))
        sage: S2 == S
        False

    As an abuse, repeated entries in ``elements`` are allowed to model
    multisets::

        sage: S1 = FiniteEnumeratedSet((1, 2, 1, 2, 2, 3))
        sage: S1
        {1, 2, 1, 2, 2, 3}

    Finaly the elements are not aware of their parent::

        sage: S.first().parent()
        Integer Ring
    """

    @staticmethod
    def __classcall__(cls, iterable):
        """
        Standard trick to expand the iterable upon input, and
        guarantees unique representation, independently of the type of
        the iterable. See ``UniqueRepresentation``.

        TESTS::

            sage: S1 = FiniteEnumeratedSet([1, 2, 3])
            sage: S2 = FiniteEnumeratedSet((1, 2, 3))
            sage: S3 = FiniteEnumeratedSet((x for x in range(1,4)))
            sage: S1 is S2
            True
            sage: S2 is S3
            True
        """
        return super(FiniteEnumeratedSet, cls).__classcall__(cls, tuple(iterable))

    def __init__(self, elements):
        """
        TESTS::

            sage: S = FiniteEnumeratedSet([1,2,3])
            sage: S == loads(dumps(S))
            True
        """
        self._elements = elements
        Parent.__init__(self, category = FiniteEnumeratedSets())

    def __repr__(self):
        """
        TESTS::

            sage: S = FiniteEnumeratedSet([1,2,3])
            sage: repr(S)
            '{1, 2, 3}'
            sage: S = FiniteEnumeratedSet([1])
            sage: repr(S)
            '{1}'
        """
        if len(self._elements) == 1: # avoid printing '{1,}'
            return "{" + str(self._elements[0]) + '}'
        return "{" + str(self._elements)[1:-1] + '}'

    def __contains__(self, x):
        """
        TESTS::

            sage: S = FiniteEnumeratedSet([1,2,3])
            sage: 1 in S
            True
            sage: 2 in S
            True
            sage: 4 in S
            False
            sage: ZZ in S
            False
        """
        return x in self._elements

    def list(self):
        """
        TESTS::

            sage: S = FiniteEnumeratedSet([1,2,3])
            sage: S.list()
            [1, 2, 3]
        """
        return list(self._elements)

    def an_element(self):
        """
        TESTS::

            sage: S = FiniteEnumeratedSet([1,2,3])
            sage: S.an_element()
            1
        """
        return self._elements[0]

    def cardinality(self):
        """
        TESTS::

            sage: S = FiniteEnumeratedSet([1,2,3])
            sage: S.cardinality()
            3
        """
        return Integer(len(self._elements))

