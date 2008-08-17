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

class RootSpace(CombinatorialFreeModule, RootLatticeRealization):
    def __init__(self, root_system, base_ring):
        """
        EXAMPLES:
            sage: r = RootSystem(['A',4]).root_lattice()
            sage: r == loads(dumps(r))
            True
        """
        self.root_system = root_system
        basis_name = "alphacheck" if root_system.dual_side else "alpha"
        CombinatorialFreeModule.__init__(self, base_ring,
                                         root_system.index_set(),
                                         element_class = RootSpaceElement,
                                         prefix=basis_name)
    def simple_root(self, i):
        """
        EXAMPLES:
            sage: r = RootSystem(['A',4]).root_lattice()
            sage: r.simple_root(1)
            alpha[1]
        """
        # This would be more efficient, but this is not used intensively,
        # and we would have to be careful with pickling:
        # self.simple_root = self.basis().__getitem__
        return self.simple_roots()[i]

    def __repr__(self):
        """
        TESTS:
            sage: RootSystem(['A',4]).root_lattice()
            Root lattice of the Root system of type ['A', 4]
            sage: RootSystem(['B',4]).root_space()
            Root space over the Rational Field of the Root system of type ['B', 4]
            sage: RootSystem(['A',4]).coroot_lattice()
            Coroot lattice of the Root system of type ['A', 4]
            sage: RootSystem(['B',4]).coroot_space()
            Coroot space over the Rational Field of the Root system of type ['B', 4]

        """
        return self._name_string()

    def _name_string(self, capitalize=True, base_ring=True, type=True):
        """
        EXAMPLES:
            sage: RootSystem(['A',4]).root_space()._name_string()
            "Root space over the Rational Field of the Root system of type ['A', 4]"
        """
        return self._name_string_helper("root", capitalize=capitalize, base_ring=base_ring, type=type)


    simple_roots = CombinatorialFreeModule.basis


class RootSpaceElement(CombinatorialFreeModuleElement, RootLatticeRealizationElement):
    def scalar(self, lambdacheck):
        """
        The scalar product between the root lattice and
        the coroot lattice.

        EXAMPLES:
            sage: r = RootSystem(['A',4]).root_lattice()
            sage: cr = RootSystem(['A',4]).coroot_lattice()
            sage: a1 = r.simple_root(1)
            sage: ac1 = cr.simple_root(1)
            sage: a1.scalar(ac1)
            2

        """
        zero = self.parent().base_ring().zero_element()
        cartan_matrix = self.parent().dynkin_diagram()
        return sum( (sum( (lambdacheck[i]*s for i,s in cartan_matrix.column(j)), zero) * c for j,c in self), zero)

