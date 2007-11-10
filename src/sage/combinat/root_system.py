#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
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

from cartan_type import CartanType
from cartan_matrix import cartan_matrix
from sage.combinat.combinatorial_algebra import CombinatorialAlgebra

def RootSystem(t):
    """

    """
    ct = CartanType(t)
    return RootSystem_generic(ct)

class RootSystem_generic:
    def __init__(self, ct):
        """
        TESTS:
            #sage: R = RootSystem(['A',3])
            #sage: R == loads(dumps(R))
            #True
        """
        self.ct = ct


    def cartan_type(self):
        """
        Returns the Cartan type of the root system.

        EXAMPLES:
            sage: R = RootSystem(['A',3])
            sage: R.cartan_type()
            ['A', 3]
        """
        return self.ct

    def cartan_matrix(self):
        """

        """
        return cartan_matrix(self.cartan_type())

    def index_set(self):
        """

        """
        return self.cartan_type().index_set()

    def coroot_lattice(self):
        raise NotImplementedError

    def root_lattice(self):
        raise NotImplementedError

    def weight_lattice(self):
        raise NotImplementedError

    def dual_weight_lattice(self):
        return self.coroot_lattice()

    def ambient_lattice(self):
        raise NotImplementedError

class RootSystem_affine(RootSystem_generic):
    def weight_lattice(self):
        return affine_weight_lattice(self)

    def dual_weight_lattice(self):
        raise NotImplementedError



class CorootLattice_generic(CombinatorialAlgebra):
    def __init__(self, ct):
        self.ct = ct
        self._prefix = "alphacheck"
        self._combinatorial_class = None


    def cartan_type(self):
        return self.ct




class RootSystemRealization_generic:
    """

    """
    ## Realization of a root system (that is of the
    ## (positive/negative/simple) roots inside some space, not necessarily
    ## with any particular structure

