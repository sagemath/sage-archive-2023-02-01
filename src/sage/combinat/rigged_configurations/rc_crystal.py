r"""
Crystal of Rigged Configurations

AUTHORS:

- Travis Scrimshaw (2010-09-26): Initial version

We only consider the highest weight crystal structure, not the
Kirillov-Reshetikhin structure, and we extend this to most
symmetrizable types.
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
from sage.misc.lazy_attribute import lazy_attribute
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.categories.highest_weight_crystals import HighestWeightCrystals
from sage.categories.regular_crystals import RegularCrystals
from sage.categories.classical_crystals import ClassicalCrystals
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.rigged_configurations.rigged_configurations import RiggedConfigurationOptions
from sage.combinat.rigged_configurations.rigged_configuration_element import (
     RCHighestWeightElement, RCHWNonSimplyLacedElement)

# Note on implementation, this class is used for simply-laced types only
class CrystalOfRiggedConfigurations(Parent, UniqueRepresentation):
    r"""
    A highest weight crystal of rigged configurations.

    The crystal structure for finite simply-laced types
    is given in [CrysStructSchilling06]_.

    INPUT:

    - ``cartan_type`` -- a (simply-laced) Cartan type

    - ``wt`` -- the highest weight vector in the weight lattice
    """
    @staticmethod
    def __classcall_private__(cls, cartan_type, wt=None):
        r"""
        Normalize the input arguments to ensure unique representation.

        EXAMPLES::
        """
        if wt is None:
            wt = cartan_type
            cartan_type = wt.parent().cartan_type()
        else:
            cartan_type = CartanType(cartan_type)
            wt_lattice = cartan_type.root_system().weight_lattice()
            wt = wt_lattice(wt)

        if not cartan_type.is_simply_laced():
            vct = cartan_type.as_folding()
            return CrystalOfNonSimplyLacedRC(vct, wt)

        return super(CrystalOfRiggedConfigurations, cls).__classcall__(cls, wt)

    def __init__(self, wt):
        r"""
        Initialize ``self``.

        EXAMPLES::
        """
        self._cartan_type = wt.parent().cartan_type()
        self._wt = wt
        self._rc_index = self._cartan_type.index_set()
        # We store the cartan matrix for the vacancy number calculations for speed
        self._cartan_matrix = self._cartan_type.cartan_matrix()
        if self._cartan_type.is_finite():
            category = ClassicalCrystals()
        else:
            category = (RegularCrystals(), HighestWeightCrystals(), InfiniteEnumeratedSets())
        Parent.__init__(self, category=category)
        n = self._cartan_type.rank() #== len(self._cartan_type.index_set())
        self.module_generators = (self.element_class( self, partition_list=[[] for i in range(n)] ),)

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

class CrystalOfNonSimplyLacedRC(CrystalOfRiggedConfigurations):
    """
    Highest weight crystal of rigged configurations in non-simply-laced type.
    """
    def __init__(self, vct, wt):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: La = RootSystem(['C',2,1]).weight_lattice().fundamental_weights()
            sage: RC = CrystalOfRiggedConfigurations(['C',2,1], La[0]); RC
            Rigged configurations of type ['C', 2, 1]
         """
        self._folded_ct = vct
        CrystalOfRiggedConfigurations.__init__(self, wt)

    @lazy_attribute
    def virtual(self):
        """
        Return the corresponding virtual crystal.

        EXAMPLES::

            sage: La = RootSystem(['C',2,1]).weight_lattice().fundamental_weights()
            sage: RC = crystals.RiggedConfigurations(['C',2,1], La[0])
            sage: RC
            B infinity rigged configurations of type ['C', 3]
            sage: RC.virtual
            B infinity rigged configurations of type ['A', 3]
        """
        P = self._folded_ct._folding.root_system().weight_lattice()
        gamma = self._folded_ct.scaling_factors()
        sigma = self._folded_ct.folding_orbit()
        vwt = P.sum_of_terms((b, gamma[a]*c) for a,c in self._wt for b in sigma[a])
        return CrystalOfRiggedConfigurations(vwt)

    def _calc_vacancy_number(self, partitions, a, i, **options):
        r"""
        Calculate the vacancy number of the `i`-th row of the `a`-th rigged
        partition.

        INPUT:

        - ``partitions`` -- the list of rigged partitions we are using

        - ``a`` -- the rigged partition index

        - ``i`` -- the row index of the `a`-th rigged partition

        TESTS::

            sage: La = RootSystem(['C',2]).weight_lattice().fundamental_weights()
            sage: RC = crystals.RiggedConfigurations(La[1])
            sage: elt = RC(partition_list=[[1], [1]])
            sage: RC._calc_vacancy_number(elt.nu(), 0, 0)
            0
            sage: RC._calc_vacancy_number(elt.nu(), 1, 0)
            -1
        """
        I = self.index_set()
        ia = I[a]
        vac_num = self._wt[ia]

        if i is None:
            row_len = float("inf")
        else:
            row_len = partitions[a][i]

        gamma = self._folded_ct.scaling_factors()
        for b, value in enumerate(self._cartan_matrix.row(a)):
            ib = I[b]
            q = partitions[b].get_num_cells_to_column(gamma[ia]*row_len, gamma[ib])
            vac_num -= value * q / gamma[ib]

        return vac_num

    def to_virtual(self, rc):
        """
        Convert ``rc`` into a rigged configuration in the virtual crystal.

        INPUT:

        - ``rc`` -- a rigged configuration element

        EXAMPLES::

            sage: La = RootSystem(['C',2,1]).weight_lattice().fundamental_weights()
            sage: RC = crystals.RiggedConfigurations(['C',2,1], La[0])
        """
        gamma = map(int, self._folded_ct.scaling_factors())
        sigma = self._folded_ct._orbit
        n = self._folded_ct._folding.rank()
        vindex = self._folded_ct._folding.index_set()
        partitions = [None] * n
        riggings = [None] * n
        vac_nums = [None] * n
        # -1 for indexing
        for a, rp in enumerate(rc):
            for i in sigma[a]:
                k = vindex.index(i)
                partitions[k] = [row_len*gamma[a] for row_len in rp._list]
                riggings[k] = [rig_val*gamma[a] for rig_val in rp.rigging]
                vac_nums[k] = [vac_num*gamma[a] for vac_num in rp.vacancy_numbers]
        return self.virtual.element_class(self.virtual, partition_list=partitions,
                            rigging_list=riggings,
                            vacancy_numbers_list=vac_nums)

    def from_virtual(self, vrc):
        """
        Convert ``vrc`` in the virtual crystal into a rigged configution of
        the original Cartan type.

        INPUT:

        - ``vrc`` -- a virtual rigged configuration

        EXAMPLES::

            sage: La = RootSystem(['C',2,1]).weight_lattice().fundamental_weights()
            sage: RC = crystals.RiggedConfigurations(['C',2,1], La[0])
        """
        gamma = list(self._folded_ct.scaling_factors()) #map(int, self._folded_ct.scaling_factors())
        sigma = self._folded_ct._orbit
        n = self._cartan_type.rank()
        partitions = [None] * n
        riggings = [None] * n
        vac_nums = [None] * n
        vindex = self._folded_ct._folding.index_set()
        # TODO: Handle special cases for A^{(2)} even and its dual?
        for a in range(n):
            index = vindex.index(sigma[a][0])
            partitions[a] = [row_len // gamma[a] for row_len in vrc[index]._list]
            riggings[a] = [rig_val / gamma[a] for rig_val in vrc[index].rigging]
            vac_nums[a] = [vac_val / gamma[a] for vac_val in vrc[index].vacancy_numbers]
        return self.element_class(self, partition_list=partitions,
                                  rigging_list=riggings, vacancy_numbers_list=vac_nums)

    Element = RCHWNonSimplyLacedElement

