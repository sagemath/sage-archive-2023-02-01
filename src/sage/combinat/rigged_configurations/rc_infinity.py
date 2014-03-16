r"""
Rigged Configurations of `\mathcal{B}(\infty)`

AUTHORS:

- Travis Scrimshaw (2013-04-16): Initial version
"""

#*****************************************************************************
#       Copyright (C) 2013 Travis Scrimshaw <tscrim@ucdavis.edu>
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
from sage.rings.all import QQ
from sage.categories.highest_weight_crystals import HighestWeightCrystals
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.rigged_configurations.rigged_configuration_element import RiggedConfigurationElement
from sage.combinat.rigged_configurations.rigged_configurations import RiggedConfigurationOptions
from sage.combinat.rigged_configurations.rigged_partition import RiggedPartition

# Note on implementation, this class is used for simply-laced types only
class InfinityCrystalOfRiggedConfigurations(Parent, UniqueRepresentation):
    r"""
    Class of rigged configurations modeling `\mathcal{B}(\infty)`.

    INPUT:

    - ``cartan_type`` -- a Cartan type
    """
    @staticmethod
    def __classcall_private__(cls, cartan_type):
        r"""
        Normalize the input arguments to ensure unique representation.

        EXAMPLES::

            sage: RC1 = InfinityCrystalOfRiggedConfigurations(CartanType(['A',3]))
            sage: RC2 = InfinityCrystalOfRiggedConfigurations(['A',3])
            sage: RC2 is RC1
            True
        """
        cartan_type = CartanType(cartan_type)
        if not cartan_type.is_simply_laced():
            vct = cartan_type.as_folding()
            return InfinityCrystalOfVirtualRC(vct)

        return super(InfinityCrystalOfRiggedConfigurations, cls).__classcall__(cls, cartan_type)

    def __init__(self, cartan_type):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: RC = InfinityCrystalOfRiggedConfigurations(['A',3])
            sage: TestSuite(RC).run()
        """
        self._cartan_type = cartan_type
        Parent.__init__(self, category=HighestWeightCrystals())
        # We store the cartan matrix for the vacancy number calculations for speed
        self._cartan_matrix = self._cartan_type.cartan_matrix()
        self.module_generators = (self.element_class(self, rigging_list=[[]]*cartan_type.rank()),)

    global_options = RiggedConfigurationOptions

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: InfinityCrystalOfRiggedConfigurations(['A',3])
            The infinity crystal of rigged configurations of type ['A', 3]
        """
        return "The infinity crystal of rigged configurations of type {}".format(self._cartan_type)

    def __iter__(self):
        """
        Returns the iterator of ``self``.

        EXAMPLES::

            sage: RC = InfinityCrystalOfRiggedConfigurations(['A', 3])
            sage: g = RC.__iter__()
            sage: g.next()
            <BLANKLINE>
            (/)
            <BLANKLINE>
            (/)
            <BLANKLINE>
            (/)
            <BLANKLINE>
            sage: g.next()
            <BLANKLINE>
            (/)
            <BLANKLINE>
            -1[ ]-1
            <BLANKLINE>
            (/)
            <BLANKLINE>
        """
        index_set = self._cartan_type.index_set()
        from sage.combinat.backtrack import TransitiveIdeal
        return TransitiveIdeal(lambda x: [x.f(i) for i in index_set],
                               self.module_generators).__iter__()

    def _element_constructor_(self, lst=None, **options):
        """
        Construct an element of ``self`` from ``lst``.
        """
        return self.element_class(self, lst, **options)

    def _calc_vacancy_number(self, partitions, a, i, **options):
        r"""
        Calculate the vacancy number of the `i`-th row of the `a`-th rigged
        partition.

        This assumes that `\gamma_a = 1` for all `a` and `(\alpha_a \mid
        \alpha_b ) = A_{ab}`.

        INPUT:

        - ``partitions`` -- The list of rigged partitions we are using

        - ``a``          -- The rigged partition index

        - ``i``          -- The row index of the `a`-th rigged partition

        TESTS::

            sage: RC = InfinityCrystalOfRiggedConfigurations(['A', 4, 1])
            sage: elt = RC(partition_list=[[1], [1], [], []])
            sage: RC._calc_vacancy_number(elt.nu(), 1, 0)
            0
        """
        row_len = partitions[a][i]
        vac_num = 0
        for b, value in enumerate(self._cartan_matrix.row(a)):
            vac_num -= value * partitions[b].get_num_cells_to_column(row_len)

        return vac_num

    class Element(RiggedConfigurationElement):
        """
        A rigged configuration for simply-laced types.

        Typically to create a specific rigged configuration, the user will pass
        in the optional argument **partition_list** and if the user wants to
        specify the rigging values, give the optional argument **rigging_list**
        as well. If **rigging_list** is not passed, the rigging values are set
        to the corresponding vacancy numbers.

        INPUT:

        - ``parent``            -- The parent of this element

        - ``rigged_partitions`` -- A list of rigged partitions

        There are two optional arguments to explicitly construct a rigged
        configuration. The first is **partition_list** which gives a list of
        partitions, and the second is **rigging_list** which is a list of
        corresponding lists of riggings. If only partition_list is specified,
        then it sets the rigging equal to the calculated vacancy numbers.

        EXAMPLES:

        Type `A_n^{(1)}` examples::

            sage: RC = InfinityCrystalOfRiggedConfigurations(['A', 4, 1])
            sage: RC(partition_list=[[2], [2, 2], [2], [2]])
            <BLANKLINE>
            0[ ][ ]0
            <BLANKLINE>
            -2[ ][ ]-2
            -2[ ][ ]-2
            <BLANKLINE>
            2[ ][ ]2
            <BLANKLINE>
            -2[ ][ ]-2
            <BLANKLINE>

            sage: RC = InfinityCrystalOfRiggedConfigurations(['A', 4, 1])
            sage: RC(partition_list=[[], [], [], []])
            <BLANKLINE>
            (/)
            <BLANKLINE>
            (/)
            <BLANKLINE>
            (/)
            <BLANKLINE>
            (/)
            <BLANKLINE>

        Type `D_n^{(1)}` examples::

            sage: RC = InfinityCrystalOfRiggedConfigurations(['D', 4, 1])
            sage: RC(partition_list=[[3], [3,2], [4], [3]])
            <BLANKLINE>
            -1[ ][ ][ ]-1
            <BLANKLINE>
            1[ ][ ][ ]1
            0[ ][ ]0
            <BLANKLINE>
            -3[ ][ ][ ][ ]-3
            <BLANKLINE>
            -1[ ][ ][ ]-1
            <BLANKLINE>

            sage: RC = InfinityCrystalOfRiggedConfigurations(['D', 4, 1])
            sage: RC(partition_list=[[1], [1,1], [1], [1]])
            <BLANKLINE>
            1[ ]1
            <BLANKLINE>
            0[ ]0
            0[ ]0
            <BLANKLINE>
            0[ ]0
            <BLANKLINE>
            0[ ]0
            <BLANKLINE>
            sage: RC(partition_list=[[1], [1,1], [1], [1]], rigging_list=[[0], [0,0], [0], [0]])
            <BLANKLINE>
            1[ ]0
            <BLANKLINE>
            0[ ]0
            0[ ]0
            <BLANKLINE>
            0[ ]0
            <BLANKLINE>
            0[ ]0
            <BLANKLINE>

        TESTS::

            sage: RC = InfinityCrystalOfRiggedConfigurations(['A', 4, 1])
            sage: elt = RC(partition_list=[[1], [1], [], []], rigging_list=[[-1], [0], [], []]); elt
            <BLANKLINE>
            -1[ ]-1
            <BLANKLINE>
            0[ ]0
            <BLANKLINE>
            (/)
            <BLANKLINE>
            (/)
            <BLANKLINE>
            sage: TestSuite(elt).run()
        """

        def epsilon(self, a):
            r"""
            Return `\varepsilon_a` of ``self``.
            """
            ep = 0
            cur = self.e(a)
            while cur is not None:
                cur = cur.e(a)
                ep += 1
            return ep

        def phi(self, a):
            r"""
            Return `\varphi_a` of ``self``.

            Let `x \in \mathcal{RC}(\infty)` Define `\varphi_a(x) :=
            \varepsilon_a(x) + \langle h_a, \mathrm{wt}(x) \rangle`,
            where `h_a` is the `a`-th simple coroot and `\mathrm{wt}(x)` is
            the :meth:`weight` of `x`.

            INPUT:

            - ``i`` -- An element of the index set

            EXAMPLES::
            """
            P = self.parent().weight_lattice_realization()
            h = P.simple_coroots()
            return self.epsilon(a) + P(self.weight()).scalar(h[a])

        def weight(self):
            """
            Return the weight of ``self``.
            """
            P = self.parent().weight_lattice_realization()
            alpha = list(P.simple_roots())
            return sum(sum(x) * alpha[i] for i,x in enumerate(self))

class InfinityCrystalOfVirtualRC(InfinityCrystalOfRiggedConfigurations):
    r"""
    Rrigged configurations for `\mathcal{B}(\infty)` in non-simply-laced types.
    """
    def __init__(self, vct):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: RC = InfinityCrystalOfRiggedConfigurations(['C',2,1]); RC
            Rigged configurations of type ['C', 2, 1]
         """
        self._folded_ct = vct
        InfinityCrystalOfRiggedConfigurations.__init__(self, vct._cartan_type)

    @lazy_attribute
    def virtual(self):
        """
        Return the corresponding virtual crystal.

        EXAMPLES::

            sage: RC = InfinityCrystalOfRiggedConfigurations(['C',2])
            sage: RC
            B infinity rigged configurations of type ['C', 3]
            sage: RC.virtual
            B infinity rigged configurations of type ['A', 3]
        """
        return InfinityCrystalOfRiggedConfigurations(self._folded_ct._folding)

    def to_virtual(self, rc):
        """
        Convert ``rc`` into a rigged configuration in the virtual crystal.

        INPUT:

        - ``rc`` -- A rigged configuration element

        EXAMPLES::

            sage: RC = InfinityCrystalOfRiggedConfigurations(['C',2])
            sage: elt = RC(partition_list=[[3],[2]]); elt
            <BLANKLINE>
            0[ ][ ][ ]0
            <BLANKLINE>
            0[ ][ ]0
            sage: velt = RC.to_virtual(elt); velt
            <BLANKLINE>
            0[ ][ ][ ]0
            <BLANKLINE>
            0[ ][ ][ ][ ]0
            <BLANKLINE>
            0[ ][ ][ ]0
            sage: velt.parent()
            B infinity rigged configurations of type ['A', 3]
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

        - ``vrc`` -- A virtual rigged configuration

        EXAMPLES::

            sage: RC = InfinityCrystalOfRiggedConfigurations(['C',2])
            sage: elt = RC(partition_list=[[3],[2]])
            sage: vrc_elt = RC.to_virtual(elt)
            sage: ret = RC.from_virtual(vrc_elt); ret
            <BLANKLINE>
            0[ ][ ][ ]0
            <BLANKLINE>
            0[ ][ ]0
            sage: ret == elt
            True
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

    class Element(InfinityCrystalOfRiggedConfigurations.Element):
        """
        Virtual rigged configuration elements.

        TESTS::
        """
        def to_virtual_configuration(self):
            """
            Return the corresponding rigged configuration in the virtual crystal.

            EXAMPLES::
            """
            return self.parent().to_virtual(self)

        def e(self, a):
            """
            Return the action of `e_a` on ``self``.

            This works by lifting into the virtual configuration, then applying

            .. MATH::

                \hat{e}_a = \prod_{j \in \iota(a)} e_j^{\gamma_j}

            and pulling back.

            EXAMPLES::
            """
            vct = self.parent()._folded_ct
            gamma = vct.scaling_factors()
            L = []
            for i in vct.folding_orbit()[a]:
                L.extend([i]*int(gamma[a]))
            virtual_rc = self.parent().to_virtual(self).e_string(L)
            if virtual_rc is None:
                return None
            return self.parent().from_virtual(virtual_rc)

        def f(self, a):
            """
            Return the action of `f_a` on ``self``.

            This works by lifting into the virtual configuration, then applying

            .. MATH::

                \hat{f}_a = \prod_{j \in \iota(a)} f_j^{\gamma_j}

            and pulling back.

            EXAMPLES::
            """
            vct = self.parent()._folded_ct
            gamma = vct.scaling_factors()
            L = []
            for i in vct.folding_orbit()[a]:
                L.extend([i]*int(gamma[a]))
            virtual_rc = self.parent().to_virtual(self).f_string(L)
            if virtual_rc is None:
                return None
            return self.parent().from_virtual(virtual_rc)

class InfinityCrystalOfRCA2Even(InfinityCrystalOfVirtualRC):
    """
    Infinity crystal of rigged configurations for type `A_{2n}^{(2)}`.
    """
    def to_virtual(self, rc):
        """
        Convert ``rc`` into a rigged configuration in the virtual crystal.

        INPUT:

        - ``rc`` -- a rigged configuration element

        EXAMPLES::

            sage: from sage.combinat.rigged_configurations.rc_infinity import InfinityCrystalOfRCA2Even
            sage: RC = InfinityCrystalOfRCA2Even(CartanType(['A',4,2]).as_folding())
            sage: elt = RC(partition_list=[[1],[1]]); elt
            <BLANKLINE>
            -1[ ]-1
            <BLANKLINE>
            1[ ]1
            <BLANKLINE>
            sage: velt = RC.to_virtual(elt); velt
            <BLANKLINE>
            -1[ ]-1
            <BLANKLINE>
            2[ ]2
            <BLANKLINE>
            -1[ ]-1
            <BLANKLINE>
            sage: velt.parent()
            Rigged configurations of type ['A', 3, 1] and factor(s) ((2, 2), (2, 2))
        """
        gamma = self._folded_ct.scaling_factors()
        sigma = self._folded_ct.folding_orbit()
        n = self._folded_ct._folding.rank()
        partitions = [None] * n
        riggings = [None] * n
        vac_nums = [None] * n
        # +/- 1 for indexing
        for a in range(len(rc)):
            for i in sigma[a]:
                partitions[i] = [row_len for row_len in rc[a]._list]
                riggings[i] = [rig_val*gamma[a] for rig_val in rc[a].rigging]
                vac_nums[i] = [vac_num*gamma[a] for vac_num in rc[a].vacancy_numbers]
        return self.virtual.element_class(self.virtual, partition_list=partitions,
                            rigging_list=riggings,
                            vacancy_numbers_list=vac_nums)

    def from_virtual(self, vrc):
        """
        Convert ``vrc`` in the virtual crystal into a rigged configution of
        the original Cartan type.

        INPUT:

        - ``vrc`` -- a virtual rigged configuration element

        EXAMPLES::

            sage: from sage.combinat.rigged_configurations.rc_infinity import InfinityCrystalOfRCA2Even
            sage: RC = InfinityCrystalOfRCA2Even(CartanType(['A',4,2]).as_folding())
            sage: elt = RC(partition_list=[[1],[1]])
            sage: velt = RC.to_virtual(elt)
            sage: ret = RC.from_virtual(velt); ret
            <BLANKLINE>
            -1[ ]-1
            <BLANKLINE>
            1[ ]1
            <BLANKLINE>
            sage: ret == elt
            True
        """
        gamma = self._folded_ct.scaling_factors()
        sigma = self._folded_ct.folding_orbit()
        n = self._cartan_type.rank()
        partitions = [None] * n
        riggings = [None] * n
        vac_nums = [None] * n
        # +/- 1 for indexing
        for a in range(n):
            index = sigma[a][0]
            partitions[index] = [row_len for row_len in vrc[index]._list]
            riggings[index] = [rig_val//gamma[a] for rig_val in vrc[index].rigging]
            vac_nums[a] = [vac_val//gamma[a] for vac_val in vrc[index].vacancy_numbers]
        return self.element_class(self, partition_list=partitions,
                                  rigging_list=riggings, vacancy_numbers_list=vac_nums)

