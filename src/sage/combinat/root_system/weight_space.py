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
from sage.combinat.family import Family
from root_lattice_realization import RootLatticeRealizationElement
from weight_lattice_realization import WeightLatticeRealization
from sage.rings.all import ZZ


class WeightSpace(CombinatorialFreeModule, WeightLatticeRealization):

    def __init__(self, root_system, base_ring):
        self.root_system = root_system
        basis_name = "Lambdacheck" if root_system.dualSide else "Lambda"
        CombinatorialFreeModule.__init__(self, base_ring,\
                                         root_system.index_set(),\
                                         element_class = WeightSpaceElement,\
                                         prefix=basis_name)

    def __repr__(self):
        """
        TEST:
            sage: RootSystem(['A',4]).weight_lattice()
            The weight lattice of the Root system of type ['A', 4]
            sage: RootSystem(['B',4]).weight_space()
            The weight space over the Rational Field of the Root system of type ['B', 4]
            sage: RootSystem(['A',4]).coweight_lattice()
            The coweight lattice of the Root system of type ['A', 4]
            sage: RootSystem(['B',4]).coweight_space()
            The coweight space over the Rational Field of the Root system of type ['B', 4]

        """
        return "The %s"%self.root_system.dualString+\
               ("weight lattice " if self.base_ring() == ZZ else "weight space over the %s "%self.base_ring())+\
               "of the %s"%(self.root_system.dual if self.root_system.dualSide else self.root_system)

    fundamental_weights = CombinatorialFreeModule.basis

    def simple_root(self, j):
        """
        Returns the $j$-th simple root

        TESTS:
            sage: R = RootSystem(["C",4])
            sage: R.weight_lattice().simple_root(3)
            -Lambda[2] + 2*Lambda[3] - Lambda[4]

            sage: R.weight_lattice().simple_roots()
            Finite family {1: 2*Lambda[1] - Lambda[2], 2: -Lambda[1] + 2*Lambda[2] - Lambda[3], 3: -Lambda[2] + 2*Lambda[3] - Lambda[4], 4: -2*Lambda[3] + 2*Lambda[4]}

        """
        assert(j in self.index_set())
        return self._from_dict(dict([(i,c) for (i,c) in self.root_system.dynkin_diagram().column(j)]))


class WeightSpaceElement(CombinatorialFreeModuleElement, RootLatticeRealizationElement):

    def scalar(self, lambdacheck):
        """
        The canonical scalar product between the weight lattice and
        the coroot lattice.

        TODO: merge with_apply_multi_module_morphism
        """
        zero = self.parent().base_ring().zero_element()
        if len(self) < len(lambdacheck):
            return sum( (lambdacheck[i]*c for (i,c) in self), zero)
        else:
            return sum( (self[i]*c for (i,c) in lambdacheck), zero)
