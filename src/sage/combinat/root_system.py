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
from sage.modules.free_module import FreeModule
from sage.rings.all import ZZ

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



class AmbientLattice_generic:
    def __init__(self, ct):
        self.ct = ct
        self.n = ct.rank()+1

        self._free_module = FreeModule(ZZ, self.n)

    def __call__(self, i):
        raise NotImplementedError

class AmbientLattice_a(AmbientLattice_generic):

    def root(self, i, j):
        return self._free_module.gen(i) - self._free_module.gen(j)

    def simple_roots(self):
        return [ self.root(i, i+1) for i in range(self.n-1) ]

    def roots(self):
        return self.positive_roots() + self.negative_roots()

    def negative_roots(self):
        res = []
        for j in range(self.n):
            for i in range(j,self.n+1):
                res.append(  self.root(i,j) )
        return res

    def positive_roots(self):
        res = []
        for j in range(self.n):
            for i in range(j):
                res.append(  self.root(i,j) )
        return res

    def fundamental_weights(self):
        return [ sum([self(j) for j in range(i)]) for i in range(self.n)]
