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

def LatticePoset(data):
    # Data type 0: if we have a lattice, then we're happy.
    if isinstance(data,FiniteLatticePoset):
        return data
    else:
        P = Poset(data)
        if P.is_meet_semilattice() and P.is_join_semilattice():
            return FiniteLatticePoset(P)
        else:
            raise ValueError, "Not a lattice."

class FiniteLatticePoset(FinitePoset):
    """
    We assume that the argument passed to FiniteLatticePoset is the
    poset of a lattice (i.e. a poset with greatest lower bound and
    least upper bound for each pair of elements).
    """
    def __repr__(self):
        return "Finite lattice containing %s elements"%self.hasse_diagram().order()

    def meet(self,x,y):
        return self._vertex_to_element(self.hasse_diagram()._meet[x.vertex,y.vertex])

    def join(self,x,y):
        return self._vertex_to_element(self.hasse_diagram()._join[x.vertex,y.vertex])

    def is_distributive(self):
        return self.hasse_diagram().is_distributive_lattice()

    def is_complemented(self):
        return self.hasse_diagram().is_complemented_lattice()

    def complements(self):
        return self.hasse_diagram().complements()

####################################################################################

def MeetSemilattice(data):
    if isinstance(data,FiniteMeetSemilattice):
        return data
    else:
        P = Poset(data)
        if P.is_meet_semilattice():
            return FiniteMeetSemilattice(P)
        else:
            raise ValueError, "Not a meet semilattice."

class FiniteMeetSemilattice(FinitePoset):
    """
    We assume that the argument passed to MeetSemilattice is the poset
    of a meet-semilattice (i.e. a poset with greatest lower bound for
    each pair of elements).
    """
    def __repr__(self):
        return "Finite meet-semilattice containing %s elements"\
                %self.hasse_diagram().order()

    def meet(self,x,y):
        return self._vertex_to_element(self.hasse_diagram()._meet[x.vertex,y.vertex])

####################################################################################

def JoinSemilattice(data):
    if isinstance(data,FiniteJoinSemilattice):
        return data
    else:
        P = Poset(data)
        if P.is_join_semilattice():
            return FiniteJoinSemilattice(P)
        else:
            raise ValueError, "Not a join semilattice."

class FiniteJoinSemilattice(FinitePoset):
    """
    We assume that the argument passed to FiniteJoinSemilattice is the
    poset of a join-semilattice (i.e. a poset with least upper bound
    for each pair of elements).
    """
    def __repr__(self):
        return "Finite join-semilattice containing %s elements"\
                %self.hasse_diagram().order()

    def join(self,x,y):
        return self._vertex_to_element(self.hasse_diagram()._join[x.vertex,y.vertex])

####################################################################################
