r"""
Crystal of Rigged Configurations

AUTHORS:

- Travis Scrimshaw (2010-09-26): Initial version

We only consider the highest weight crystal structure, not the
Kirillov-Reshetikhin structure, and we extend this to all types.

INPUT:

- ``cartan_type`` -- A Cartan type

- ``wt`` -- the highest weight in the weight lattice
        
EXAMPLES::
"""

#*****************************************************************************
#       Copyright (C) 2013 Travis Scrimshaw <tscrim at ucdavis.edu>
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

from sage.misc.cachefunc import cached_method
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.categories.highest_weight_crystals import HighestWeightCrystals
from sage.categories.regular_crystals import RegularCrystals
from sage.categories.classical_crystals import ClassicalCrystals
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.rigged_configurations.rigged_configurations import RiggedConfigurationOptions
from sage.combinat.rigged_configurations.rigged_configuration_element import RCHighestWeightElement
from sage.combinat.rigged_configurations.rigged_partition import RiggedPartition, RiggedPartitionTypeB

# Note on implementation, this class is used for simply-laced types only
class CrystalOfRiggedConfigurations(Parent, UniqueRepresentation):
    r"""
    A highest weight crystal of rigged configurations.

    The crystal structure is given in [CrysStructSchilling06]_.

    INPUT:

    - ``cartan_type`` -- a (simply-laced) Cartan type

    - ``wt`` -- the highest weight vector in the weight lattice
    """
    @staticmethod
    def __classcall_private__(cls, cartan_type, wt):
        r"""
        Normalize the input arguments to ensure unique representation.

        EXAMPLES::
        """
        cartan_type = CartanType(cartan_type)
        wt_lattice = cartan_type.root_system().weight_lattice()
        wt = wt_lattice(wt)
        return super(CrystalOfRiggedConfigurations, cls).__classcall__(cls, cartan_type, wt)

    def __init__(self, cartan_type, wt):
        r"""
        Initialize ``self``.

        EXAMPLES::
        """
        self._cartan_type = cartan_type
        self._wt = wt
        # We store the cartan matrix for the vacancy number calculations for speed
        self._cartan_matrix = self._cartan_type.cartan_matrix()
        if self._cartan_type.is_finite():
            category = ClassicalCrystals()
        else:
            category = (RegularCrystals(), HighestWeightCrystals())
        Parent.__init__(self, category=category)
        self.module_generators = (self.element_class( self,
                                  partition_list=[[] for i in cartan_type.index_set()] ),)

    global_options = RiggedConfigurationOptions

    def _repr_(self):
        """
        Return a string representation of ``self``.
        """
        return "Crystal of rigged configurations of type {0} and weight {1}".format(
                self._cartan_type, self._wt)

    def _element_constructor_(self, *lst, **options):
        """
        Construct a ``RiggedConfigurationElement``.

        Typically the user should not call this method since it does not check
        if it is an actual configuration in the crystal. Instead the user
        should use the iterator.

        EXAMPLES::
        """
        return self.element_class(self, list(lst), **options)

    def _calc_vacancy_number(self, partitions, a, i, **options):
        r"""
        Calculate the vacancy number of the `i`-th row of the `a`-th rigged
        partition.

        This assumes that `\gamma_a = 1` for all `a` and `(\alpha_a \mid
        \alpha_b ) = A_{ab}`.

        INPUT:

        - ``partitions`` -- the list of rigged partitions we are using

        - ``a`` -- the rigged partition index

        - ``i`` -- the row index of the `a`-th rigged partition

        TESTS::
        """
        vac_num = self._wt[self.index_set()[a]]

        if i is None:
            row_len = float('inf')
        else:
            row_len = partitions[a][i]

        for b, value in enumerate(self._cartan_matrix.row(a)):
            vac_num -= value * partitions[b].get_num_cells_to_column(row_len)

        return vac_num

    Element = RCHighestWeightElement

