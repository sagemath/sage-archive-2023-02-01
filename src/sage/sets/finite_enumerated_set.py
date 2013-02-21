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
from sage.categories.sets_cat import EmptySetError
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
        Category of facade finite enumerated sets
        sage: TestSuite(S).run()

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
        return super(FiniteEnumeratedSet, cls).__classcall__(
                cls,
                tuple(iterable))

    def __init__(self, elements):
        """
        TESTS::

            sage: TestSuite(FiniteEnumeratedSet([1,2,3])).run()
            sage: TestSuite(FiniteEnumeratedSet([])).run()
        """
        self._elements = elements
        Parent.__init__(self, facade = True, category = FiniteEnumeratedSets())

    def __nonzero__(self):
        r"""
        Conversion to boolean.

        EXAMPLES::

            sage: bool(FiniteEnumeratedSet('abc'))
            True
            sage: bool(FiniteEnumeratedSet([]))
            False
        """
        return bool(self._elements)

    def _repr_(self):
        """
        TESTS::

            sage: S = FiniteEnumeratedSet([1,2,3])
            sage: repr(S)
            '{1, 2, 3}'
            sage: S = FiniteEnumeratedSet(['1','2','3'])
            sage: repr(S)
            "{'1', '2', '3'}"
            sage: S = FiniteEnumeratedSet([1])
            sage: repr(S)
            '{1}'
        """
        return '{' + ', '.join(repr(e) for e in self._elements) + '}'

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

            sage: S.is_parent_of(2)
            True
            sage: S.is_parent_of(4)
            False
        """
        return x in self._elements

    is_parent_of = __contains__

    def __iter__(self):
        r"""
        Iterator over the element of self.

        EXAMPLES::

            sage: for i in FiniteEnumeratedSet([1,2,3]): print i
            1
            2
            3
        """
        return iter(self._elements)

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
        if not self._elements:
            raise EmptySetError
        return self._elements[0]

    def first(self):
        r"""
        Return the first element of the enumeration or raise an EmptySetError if
        the set is empty.

        EXAMPLES::

            sage: S = FiniteEnumeratedSet('abc')
            sage: S.first()
            'a'
        """
        if not self._elements:
            raise EmptySetError
        return self._elements[0]

    def last(self):
        r"""
        Returns the last element of the iteration or raise an EmptySetError if
        the set is empty.

        EXAMPLES::

            sage: S = FiniteEnumeratedSet([0,'a',1.23, 'd'])
            sage: S.last()
            'd'
        """
        if not self._elements:
            raise EmptySetError
        return self._elements[-1]

    def random_element(self):
        r"""
        Return a random element.

        EXAMPLES::

            sage: S = FiniteEnumeratedSet('abc')
            sage: S.random_element()   # random
            'b'
        """
        if not self._elements:
            raise EmptySetError
        from sage.misc.prandom import choice
        return choice(self._elements)

    def cardinality(self):
        """
        TESTS::

            sage: S = FiniteEnumeratedSet([1,2,3])
            sage: S.cardinality()
            3
        """
        return Integer(len(self._elements))

    def rank(self, x):
        """
        Returns the index of ``x`` in this finite enumerated set.

        EXAMPLES::

            sage: S = FiniteEnumeratedSet(['a','b','c'])
            sage: S.index('b')
            1
        """
        return self._elements.index(x)

    index = rank

    def unrank(self,i):
        r"""
        Return the element at position i.

        EXAMPLES::

            sage: S = FiniteEnumeratedSet([1,'a',-51])
            sage: S[0], S[1], S[2]
            (1, 'a', -51)
            sage: S[3]
            Traceback (most recent call last):
            ...
            IndexError: list index out of range
            sage: S[-1], S[-2], S[-3]
            (-51, 'a', 1)
            sage: S[-4]
            Traceback (most recent call last):
            ...
            IndexError: list index out of range
        """
        return self._elements[i]

    def _element_constructor_(self, el):
        """
        TESTS::

            sage: S = FiniteEnumeratedSet([1,2,3])
            sage: one,two,three = sorted(S)
            sage: S(1) is one
            True
            sage: x = 3
            sage: x is three
            False
            sage: S(x) is three
            True
            sage: S(0)
            Traceback (most recent call last):
            ...
            ValueError: 0 not in {1, 2, 3}
        """
        try:
            return self._elements[self.rank(el)]
        except (ValueError,KeyError):
            raise ValueError("%s not in %s"%(el, self))
