#*****************************************************************************
#       Copyright (C) 2007 Nicolas M. Thiery <nthiery at users.sf.net>
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

from sage.combinat.free_module import CombinatorialFreeModule, CombinatorialFreeModuleElement
from root_lattice_realization import RootLatticeRealization, RootLatticeRealizationElement
from sage.rings.all import ZZ


class RootSpace(CombinatorialFreeModule, RootLatticeRealization):

    def __init__(self, root_system, base_ring):
        self.root_system = root_system
        basis_name = "alphacheck" if root_system.dualSide else "alpha"
        CombinatorialFreeModule.__init__(self, base_ring,\
                                         root_system.index_set(),\
                                         element_class = RootSpaceElement,\
                                         prefix=basis_name)
        # This would be more efficient, but this is not used intensively,
        # and we would have to be careful with pickling:
        # self.simple_root = self.basis().__getitem__

    def simple_root(self, i):
        return self.basis()[i]

    def __repr__(self):
        """
        TEST:
            sage: RootSystem(['A',4]).root_lattice()
            The root lattice of the Root system of type ['A', 4]
            sage: RootSystem(['B',4]).root_space()
            The root space over the Rational Field of the Root system of type ['B', 4]
            sage: RootSystem(['A',4]).coroot_lattice()
            The coroot lattice of the Root system of type ['A', 4]
            sage: RootSystem(['B',4]).coroot_space()
            The coroot space over the Rational Field of the Root system of type ['B', 4]

        """
        return "The %s"%self.root_system.dualString+\
               ("root lattice " if self.base_ring() == ZZ else "root space over the %s "%self.base_ring())+\
               "of the %s"%(self.root_system.dual if self.root_system.dualSide else self.root_system)

    simple_roots = CombinatorialFreeModule.basis


class RootSpaceElement(CombinatorialFreeModuleElement, RootLatticeRealizationElement):

    def scalar(self, lambdacheck):
        """
        The scalar product between the root lattice and
        the coroot lattice.
        """
        zero = self.parent().base_ring().zero_element()
        cartan_matrix = self.parent().dynkin_diagram()
        return sum( (sum( (lambdacheck[i]*s for (i,s) in cartan_matrix.column(j)), zero) * c
                         for (j,c) in self), zero)

