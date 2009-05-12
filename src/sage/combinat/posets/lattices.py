r"""
SemiLattices and Lattices
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
from sage.combinat.posets.posets import Poset, FinitePoset
from sage.combinat.posets.elements import (LatticePosetElement,
                                           MeetSemilatticeElement,
                                           JoinSemilatticeElement)
import copy

####################################################################################

def MeetSemilattice(data):
    r"""
    Construct a meet semi-lattice from various forms of input data.

    INPUT:

    - ``data`` - any data that defines a poset that is also a meet
      semilattice. See the documentation for ``Poset``.

    EXAMPLES:

    Using data that defines a poset::

          sage: MeetSemilattice([[1,2],[3],[3]])
          Finite meet-semilattice containing 4 elements

    Using a previously constructed poset::

          sage: P = Poset([[1,2],[3],[3]])
          sage: L = MeetSemilattice(P); L
          Finite meet-semilattice containing 4 elements
          sage: type(L)
          <class 'sage.combinat.posets.lattices.FiniteMeetSemilattice'>

    If the data is not a lattice, then an error is raised::

          sage: elms = [1,2,3,4,5,6,7]
          sage: rels = [[1,2],[3,4],[4,5],[2,5]]
          sage: MeetSemilattice((elms, rels))
          Traceback (most recent call last):
          ...
          ValueError: Not a meet semilattice.
    """
    if isinstance(data,FiniteMeetSemilattice):
        return data
    else:
        P = Poset(data)
        if P.is_meet_semilattice():
            M = copy.copy(P)
            M.__class__ = FiniteMeetSemilattice
            return M
        else:
            raise ValueError, "Not a meet semilattice."

class FiniteMeetSemilattice(FinitePoset):
    """
    ..note::
        We assume that the argument passed to MeetSemilattice is the poset
        of a meet-semilattice (i.e. a poset with greatest lower bound for
        each pair of elements).

    TESTS::

        sage: M = MeetSemilattice([[1,2],[3],[3]])
        sage: M == loads(dumps(M))
        True

    ::

        sage: P = Poset([[1,2],[3],[3]])
        sage: M = MeetSemilattice(P)
        sage: M == loads(dumps(M))
        True

    """
    _element_type = MeetSemilatticeElement

    def __repr__(self):
        r"""
        TESTS::

            sage: M = MeetSemilattice([[1,2],[3],[3]])
            sage: M.__repr__()
            'Finite meet-semilattice containing 4 elements'

        ::

            sage: P = Poset([[1,2],[3],[3]])
            sage: M = MeetSemilattice(P)
            sage: M.__repr__()
            'Finite meet-semilattice containing 4 elements'
        """
        return "Finite meet-semilattice containing %s elements"\
                %self._hasse_diagram.order()

    def meet(self,x,y):
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
        return self._vertex_to_element(self._hasse_diagram._meet[x.vertex,y.vertex])

####################################################################################

def JoinSemilattice(data):
    r"""
    Construct a join semi-lattice from various forms of input data.

    INPUT:

    - ``data`` - any data that defines a poset that is also a join
      semilattice. See the documentation for ``Poset``.

    EXAMPLES:

    Using data that defines a poset::

          sage: JoinSemilattice([[1,2],[3],[3]])
          Finite join-semilattice containing 4 elements

    Using a previously constructed poset::

          sage: P = Poset([[1,2],[3],[3]])
          sage: J = JoinSemilattice(P); J
          Finite join-semilattice containing 4 elements
          sage: type(J)
          <class 'sage.combinat.posets.lattices.FiniteJoinSemilattice'>

    If the data is not a lattice, then an error is raised::

          sage: elms = [1,2,3,4,5,6,7]
          sage: rels = [[1,2],[3,4],[4,5],[2,5]]
          sage: JoinSemilattice((elms, rels))
          Traceback (most recent call last):
          ...
          ValueError: Not a join semilattice.
    """
    if isinstance(data,FiniteJoinSemilattice):
        return data
    else:
        P = Poset(data)
        if P.is_join_semilattice():
            J = copy.copy(P)
            J.__class__ = FiniteJoinSemilattice
            return J
        else:
            raise ValueError, "Not a join semilattice."

class FiniteJoinSemilattice(FinitePoset):
    """
    We assume that the argument passed to FiniteJoinSemilattice is the
    poset of a join-semilattice (i.e. a poset with least upper bound
    for each pair of elements).

    TESTS::

        sage: J = JoinSemilattice([[1,2],[3],[3]])
        sage: J == loads(dumps(J))
        True

    ::

        sage: P = Poset([[1,2],[3],[3]])
        sage: J = JoinSemilattice(P)
        sage: J == loads(dumps(J))
        True

    """
    _element_type = JoinSemilatticeElement

    def __repr__(self):
        r"""
        TESTS::

            sage: J = JoinSemilattice([[1,2],[3],[3]])
            sage: J.__repr__()
            'Finite join-semilattice containing 4 elements'

        ::

            sage: P = Poset([[1,2],[3],[3]])
            sage: J = JoinSemilattice(P)
            sage: J.__repr__()
            'Finite join-semilattice containing 4 elements'
        """
        return "Finite join-semilattice containing %s elements"\
                %self._hasse_diagram.order()

    def join(self,x,y):
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
        return self._vertex_to_element(self._hasse_diagram._join[x.vertex,y.vertex])

####################################################################################

def LatticePoset(data):
    r"""
    Construct a lattice from various forms of input data.

    INPUT:

    - ``data`` - any data that defines a poset. See the documentation
      for ``Poset``.

    OUTPUT:

        FiniteLatticePoset -- an instance of FiniteLatticePoset

    EXAMPLES:

    Using data that defines a poset::

          sage: LatticePoset([[1,2],[3],[3]])
          Finite lattice containing 4 elements

    Using a previously constructed poset::

          sage: P = Poset([[1,2],[3],[3]])
          sage: L = LatticePoset(P); L
          Finite lattice containing 4 elements
          sage: type(L)
          <class 'sage.combinat.posets.lattices.FiniteLatticePoset'>

    If the data is not a lattice, then an error is raised::

          sage: elms = [1,2,3,4,5,6,7]
          sage: rels = [[1,2],[3,4],[4,5],[2,5]]
          sage: LatticePoset((elms, rels))
          Traceback (most recent call last):
          ...
          ValueError: Not a lattice.
    """
    if isinstance(data,FiniteLatticePoset):
        return data
    else:
        P = Poset(data)
        if P.is_meet_semilattice() and P.is_join_semilattice():
            L = copy.copy(P)
            L.__class__ = FiniteLatticePoset
            return L
        else:
            raise ValueError, "Not a lattice."

class FiniteLatticePoset(FiniteMeetSemilattice, FiniteJoinSemilattice):
    """
    We assume that the argument passed to FiniteLatticePoset is the
    poset of a lattice (i.e. a poset with greatest lower bound and
    least upper bound for each pair of elements).

    TESTS::

        sage: L = LatticePoset([[1,2],[3],[3]])
        sage: L == loads(dumps(L))
        True

    ::

        sage: P = Poset([[1,2],[3],[3]])
        sage: L = LatticePoset(P)
        sage: L == loads(dumps(L))
        True

    """
    _element_type = LatticePosetElement

    def __repr__(self):
        r"""
        TESTS::

            sage: L = LatticePoset([[1,2],[3],[3]])
            sage: L.__repr__()
            'Finite lattice containing 4 elements'

        ::

            sage: P = Poset([[1,2],[3],[3]])
            sage: L = LatticePoset(P)
            sage: L.__repr__()
            'Finite lattice containing 4 elements'
        """
        return "Finite lattice containing %s elements"%self._hasse_diagram.order()

    def is_distributive(self):
        r"""
        Returns ``True`` if the lattice is distributive, and ``False``
        otherwise.

        EXAMPLES::

            sage: L = LatticePoset({0:[1,2],1:[3],2:[3]})
            sage: L.is_distributive()
            True
            sage: L = LatticePoset({0:[1,2,3],1:[4],2:[4],3:[4]})
            sage: L.is_distributive()
            False
        """
        return self._hasse_diagram.is_distributive_lattice()

    def is_complemented(self):
        r"""
        Returns ``True`` if ``self`` is a complemented lattice, and
        ``False`` otherwise.

        EXAMPLES::

            sage: L = LatticePoset({0:[1,2,3],1:[4],2:[4],3:[4]})
            sage: L.is_complemented()
            True

            sage: L = LatticePoset({0:[1,2],1:[3],2:[3],3:[4]})
            sage: L.is_complemented()
            False
        """
        return self._hasse_diagram.is_complemented_lattice()

    def complements(self):
        r"""
        Returns a list of the elements of the lattice.

        A complement of ``x`` is an element ``y`` such that the meet
        of ``x`` and ``y`` is the bottom element of ``self`` and the
        join of ``x`` and ``y`` is the top element of ``self``.

        EXAMPLES::

            sage: L = LatticePoset({0:[1,2,3],1:[4],2:[4],3:[4]})
            sage: L.complements()
            [4, 3, 3, 2, 0]

            sage: L = LatticePoset({0:[1,2],1:[3],2:[3],3:[4]})
            sage: L.complements()
            [4, None, None, None, 0]
        """
        return self._hasse_diagram.complements()

####################################################################################
