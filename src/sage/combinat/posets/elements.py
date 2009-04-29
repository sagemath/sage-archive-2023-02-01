r"""
Elements of posets, lattices, semilattices, etc.
"""
#*****************************************************************************
#       Copyright (C) 2008 Peter Jipsen <jipsen@chapman.edu>,
#                          Franco Saliola <saliola@gmail.com>
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
#*****************************************************************************
from sage.structure.element import Element

class PosetElement(Element):
    def __init__(self, poset, element, vertex):
        r"""
        Establishes the parent-child relationship between ``poset``
        and ``element``, where ``element`` is associated to the
        vertex ``vertex`` of the Hasse diagram of the poset.

        INPUT:

        - ``poset`` - a poset object

        - ``element`` - any object

        - ``vertex`` - a vertex of the Hasse diagram of the poset

        TESTS::

            sage: from sage.combinat.posets.elements import PosetElement
            sage: P = Poset([[1,2],[4],[3],[4],[]])
            sage: e = PosetElement(P, "elm", 0)
            sage: e.parent() is P
            True
            sage: e == loads(dumps(e))
            True
        """
        Element.__init__(self, poset)
        if isinstance(element, self.parent()._element_type):
            self.element = element.element
        else:
            self.element = element
        self.vertex = vertex

    def __repr__(self):
        """
        TESTS::

            sage: repr(Poset([[1,2],[4],[3],[4],[]])(0))
            '0'
        """
        return "%s" %str(self.element)

    def __eq__(self,other):
        """
        TESTS::

            sage: P = Poset([[1,2],[4],[3],[4],[]])
            sage: P(0).__eq__(P(4))
            False
            sage: from sage.combinat.posets.elements import PosetElement
            sage: PosetElement(P,0,3) == PosetElement(P,0,3)
            True
            sage: PosetElement(P,1,3) == PosetElement(P,0,3)
            False
            sage: PosetElement(P,0,2) == PosetElement(P,0,3)
            False
        """
        return self.parent() == other.parent() \
                and self.element == other.element \
                and self.vertex == other.vertex

    def _cmp(self,other):
        """
        TESTS::

            sage: P = Poset([[1,2],[4],[3],[4],[]])
            sage: P(0)._cmp(P(4))
            -1
            sage: P(4)._cmp(P(0))
            1
            sage: P(0)._cmp(P(0))
            0
        """
        if self.parent() != other.parent():
            raise ValueError, "The elements are not in the same poset."
        if self == other: return 0
        if self.parent().is_less_than(self,other): return -1
        if self.parent().is_less_than(other,self): return  1
        return None

    def __lt__(self,other):
        """
        TESTS

        ::

            sage: dag = DiGraph({0:[2,3], 1:[3,4], 2:[5], 3:[5], 4:[5]})
            sage: P = Poset(dag)
            sage: P(0) < P(1)
            False
            sage: P(4) < P(1)
            False
            sage: P(0) < P(0)
            False
        """
        return self._cmp(other) == -1 or False

    def __le__(self,other):
        """
        TESTS

        ::

            sage: dag = DiGraph({0:[2,3], 1:[3,4], 2:[5], 3:[5], 4:[5]})
            sage: P = Poset(dag)
            sage: P(1) <= P(0)
            False
            sage: P(0) <= P(1)
            False
            sage: P(0) <= P(3)
            True
            sage: P(0) <= P(0)
            True
        """
        return self == other or self._cmp(other) == -1 or False

    def __gt__(self,other):
        """
        TESTS

        ::

            sage: dag = DiGraph({0:[2,3], 1:[3,4], 2:[5], 3:[5], 4:[5]})
            sage: P = Poset(dag)
            sage: P(0).__gt__(P(5))
            False
            sage: P(5).__gt__(P(0))
            True
            sage: P(0).__gt__(P(0))
            False
        """
        return self._cmp(other) == 1 or False

    def __ge__(self,other):
        """
        TESTS

        ::

            sage: dag = DiGraph({0:[2,3], 1:[3,4], 2:[5], 3:[5], 4:[5]})
            sage: P = Poset(dag)
            sage: P(0).__ge__(P(5))
            False
            sage: P(5).__ge__(P(0))
            True
            sage: P(0).__ge__(P(0))
            True
        """
        return self == other or self._cmp(other) == 1 or False

class MeetSemilatticeElement(PosetElement):
    def __mul__(self,other):
        r"""
        Return the meet of ``self`` and ``other`` in the lattice.

        EXAMPLES::

            sage: D = Posets.DiamondPoset(5)
            sage: D(1) * D(2)
            0
            sage: D(1) * D(1)
            1
            sage: D(1) * D(0)
            0
            sage: D(1) * D(4)
            1
        """
        return self.parent().meet(self,other)

class JoinSemilatticeElement(PosetElement):
    def __add__(self,other):
        r"""
        Return the join of ``self`` and ``other`` in the lattice.

        EXAMPLES::

            sage: D = Posets.DiamondPoset(5)
            sage: D(1) + D(2)
            4
            sage: D(1) + D(1)
            1
            sage: D(1) + D(4)
            4
            sage: D(1) + D(0)
            1
        """
        return self.parent().join(self,other)

class LatticePosetElement(MeetSemilatticeElement,JoinSemilatticeElement):
    pass
