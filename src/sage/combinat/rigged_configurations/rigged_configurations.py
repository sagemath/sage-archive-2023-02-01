r"""
Rigged Configurations

AUTHORS:

- Travis Scrimshaw (2010-09-26): Initial version

Rigged configurations form combinatorial objects first introduced by Kerov, Kirillov and Reshetikhin
that arose from studies of statistical mechanical models using the Bethe Ansatz.
They are sequences of rigged partitions. A rigged partition is a partition together
with a label associated to each part that satisfy certain constraints. The labels
are also called riggings.

Rigged configurations exist for all affine Kac-Moody Lie algebras. See for example
[HKOTT2002]_. In Sage they are specified by providing a Cartan type and a list of
rectangular shapes `R`. The list of all (highest weight) rigged configurations
for given `R` is computed via the (virtual) Kleber algorithm (see also
:class:`~sage.combinat.rigged_configurations.kleber_tree.KleberTree` and
:class:`~sage.combinat.rigged_configurations.kleber_tree.VirtualKleberTree`).

There exists a crystal structure on the set of rigged configurations, see
[CrysStructSchilling06]_ and [OSS03]_. The highest weight rigged
configurations are those where all riggings are nonnegative. The list of
all rigged configurations is computed from the highest weight ones using the
crystal operators.

Rigged configurations are in bijection with tensor products of Kirillov-Reshetikhin tableaux,
see [RigConBijection]_, [BijectionDn]_, and [BijectionLRT]_. The list of rectangles `R` corresponds
to the shape of the Kirillov-Reshetikhin crystals appearing in the tensor product.
Kirillov-Reshetikhin crystals are implemented in Sage, see :class:`KirillovReshetikhinCrystal`, however,
in the bijection with rigged configurations a different realization of the elements in the
crystal are obtained, which are coined Kirillov-Reshetkihin tableaux, see :class:`KirillovReshetikhinTableaux`.
For more details see [AffineRigConDn]_.

INPUT:

- ``cartan_type`` -- A Cartan type

- ``B`` -- A list of positive integer pairs `(r,s)` specifying the width `s`
  and height `r` of the sequence of rectangles

EXAMPLES::

    sage: RC = RiggedConfigurations(['A', 3, 1], [[3, 2], [1, 2], [1, 1]])
    sage: RC
    Rigged configurations of type ['A', 3, 1] and factor(s) ((3, 2), (1, 2), (1, 1))

    sage: RC = RiggedConfigurations(['A', 3, 1], [[2,1]]); RC
    Rigged configurations of type ['A', 3, 1] and factor(s) ((2, 1),)
    sage: RC.cardinality()
    6
    sage: len(RC.list()) == RC.cardinality()
    True
    sage: RC.list()    # random
    [
    (/)
    <BLANKLINE>
    (/)
    <BLANKLINE>
    (/)
    ,
    (/)
    <BLANKLINE>
    -1[ ]-1
    <BLANKLINE>
    (/)
    ,
    (/)
    <BLANKLINE>
    0[ ]0
    <BLANKLINE>
    -1[ ]-1
    ,
    -1[ ]-1
    <BLANKLINE>
    0[ ]0
    <BLANKLINE>
    (/)
    ,
    -1[ ]-1
    <BLANKLINE>
    1[ ]1
    <BLANKLINE>
    -1[ ]-1
    ,
    0[ ]0
    <BLANKLINE>
    -1[ ]-1
    -1[ ]-1
    <BLANKLINE>
    0[ ]0
    ]

A rigged configuration element with all riggings equal to the vacancy numbers can be created as follows::

    sage: RC = RiggedConfigurations(['A', 3, 1], [[3,2], [2,1], [1,1], [1,1]]); RC
    Rigged configurations of type ['A', 3, 1] and factor(s) ((3, 2), (2, 1), (1, 1), (1, 1))
    sage: elt = RC(partition_list=[[1],[],[]]); elt
    <BLANKLINE>
    0[ ]0
    <BLANKLINE>
    (/)
    <BLANKLINE>
    (/)
    <BLANKLINE>

If on the other hand we also want to specify the riggings, this can be achieved as follows::

    sage: RC = RiggedConfigurations(['D', 7, 1], [[3,3],[5,2],[4,3],[2,3],[4,4],[3,1],[1,4],[2,2]])
    sage: elt = RC(partition_list=[[2],[3,2,1],[2,2,1,1],[2,2,1,1,1,1],[3,2,1,1,1,1],[2,1,1],[2,2]],
    ....:          rigging_list=[[2],[1,0,0],[4,1,2,1],[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,0],[0,0]])
    sage: elt
    <BLANKLINE>
    3[ ][ ]2
    <BLANKLINE>
    1[ ][ ][ ]1
    2[ ][ ]0
    1[ ]0
    <BLANKLINE>
    4[ ][ ]4
    4[ ][ ]1
    3[ ]2
    3[ ]1
    <BLANKLINE>
    2[ ][ ]1
    2[ ][ ]0
    0[ ]0
    0[ ]0
    0[ ]0
    0[ ]0
    <BLANKLINE>
    0[ ][ ][ ]0
    2[ ][ ]1
    0[ ]0
    0[ ]0
    0[ ]0
    0[ ]0
    <BLANKLINE>
    0[ ][ ]0
    0[ ]0
    0[ ]0
    <BLANKLINE>
    0[ ][ ]0
    0[ ][ ]0
    <BLANKLINE>

To obtain the Kirillov-Reshetikhin (KR) tableaux under the bijection between rigged configurations and KR
tableaux, we can type the following. This example was checked against Reiho Sakamoto's Mathematica program
on rigged configurations::

    sage: output = elt.to_tensor_product_of_kirillov_reshetikhin_tableaux(); output
    [[1, 1, 1], [2, 3, 3], [3, 4, -5]] (X) [[1, 1], [2, 2], [3, 3], [5, -6], [6, -5]] (X) [[1, 1, 2], [2, 2, 3], [3, 3, 7], [4, 4, -7]] (X) [[1, 1, 1], [2, 2, 2]] (X) [[1, 1, 1, 3], [2, 2, 3, 4], [3, 3, 4, 5], [4, 4, 5, 6]] (X) [[1], [2], [3]] (X) [[1, 1, 1, 1]] (X) [[1, 1], [2, 2]]
    sage: elt.to_tensor_product_of_kirillov_reshetikhin_tableaux().to_rigged_configuration() == elt
    True
    sage: output.to_rigged_configuration().to_tensor_product_of_kirillov_reshetikhin_tableaux() == output
    True

To get the highest weight rigged configurations, use the attribute
``module_generators``::

    sage: RC = RiggedConfigurations(['D',4,1], [[2,1]])
    sage: for x in RC.module_generators: x
    <BLANKLINE>
    (/)
    <BLANKLINE>
    (/)
    <BLANKLINE>
    (/)
    <BLANKLINE>
    (/)
    <BLANKLINE>
    0[ ]0
    <BLANKLINE>
    0[ ]0
    0[ ]0
    <BLANKLINE>
    0[ ]0
    <BLANKLINE>
    0[ ]0

TESTS::

    sage: RC = RiggedConfigurations(['A', 3, 1], [[3,2], [2,1], [1,1], [1,1]])
    sage: len(RC.module_generators)
    17
    sage: RC = RiggedConfigurations(['D', 4, 1], [[1, 1]])
    sage: RC.cardinality()
    8
    sage: RC = RiggedConfigurations(['D', 4, 1], [[2, 1]])
    sage: c = RC.cardinality(); c
    29
    sage: K = KirillovReshetikhinCrystal(['D',4,1],2,1)
    sage: K.cardinality() == c
    True

REFERENCES:

.. [HKOTT2002] G. Hatayama, A. Kuniba, M. Okado, T. Takagi, Z. Tsuboi.
   Paths, Crystals and Fermionic Formulae
   Prog.Math.Phys. 23 (2002) 205-272

.. [CrysStructSchilling06] Anne Schilling.
   Crystal structure on rigged configurations.
   International Mathematics Research Notices.
   Volume 2006. 2006. Article ID 97376. Pages 1-27.

.. [RigConBijection] Masato Okado, Anne Schilling, Mark Shimozono.
   A crystal to rigged configuration bijection for non-exceptional affine
   algebras.
   Algebraic Combinatorics and Quantum Groups.
   Edited by N. Jing. World Scientific. 2003. Pages 85-124.

.. [BijectionDn] Anne Schilling.
   A bijection between type `D_n^{(1)}` crystals and rigged configurations.
   J. Algebra. 285. 2005. 292-334

.. [BijectionLRT] Anatol N. Kirillov, Anne Schilling, Mark Shimozono.
   A bijection between Littlewood-Richardson tableaux and rigged
   configurations.
   Selecta Mathematica (N.S.). 8. 2002. 67-135 (:mathscinet:`MR1890195`).

.. [AffineRigConDn] Masato Okado, Reiho Sakamoto, Anne Schilling.
   Affine crystal structure on rigged configurations of type `D_n^{(1)}`.
   J. Algebraic Combinatorics, to appear, :doi:`10.1007/s10801-012-0383-z` (:arXiv:`1109.3523`)
"""

#*****************************************************************************
#       Copyright (C) 2010-2012 Travis Scrimshaw <tscrim@ucdavis.edu>
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
from sage.combinat.misc import IterableFunctionCall
from sage.rings.all import ZZ, QQ
from sage.categories.finite_crystals import FiniteCrystals
from sage.categories.regular_crystals import RegularCrystals
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.cartesian_product import CartesianProduct
from sage.combinat.rigged_configurations.kleber_tree import KleberTree, VirtualKleberTree
from sage.combinat.rigged_configurations.rigged_configuration_element import RiggedConfigurationElement, \
  RCVirtualElement

# Note on implementation, this class is used for simply-laced types only
class RiggedConfigurations(Parent, UniqueRepresentation):
    r"""
    Class of rigged configurations.

    Let `\overline{I}` denote the classical index set associated to the Cartan type
    of the rigged configurations. A rigged configuration of multiplicity
    array `L_i^{(a)}` and dominant weight `\Lambda` is a sequence of partitions
    `\{ \nu^{(a)} \mid a \in \overline{I} \}` such that

    .. MATH::

        \sum_{\overline{I} \times \mathbb{Z}_{>0}} i m_i^{(a)} \alpha_a
        = \sum_{\overline{I} \times \mathbb{Z}_{>0}} i L_i^{(a)} \Lambda_a
        - \Lambda

    where `\alpha_a` is a simple root, `\Lambda_a` is a fundamental weight,
    and `m_i^{(a)}` is the number of rows of length `i` in the partition
    `\nu^{(a)}`.

    Each sequence of partitions also comes with a sequence of values called
    vacancy numbers and riggings. Vacancy numbers are computed based upon the
    partitions and `L_i^{(a)}`, and the riggings must satisfy certain
    conditions. For more, see [RigConBijection]_ [CrysStructSchilling06]_
    [BijectionLRT]_.

    They are in bijection with
    :class:`TensorProductOfKirillovReshetikhinTableaux` of non-exceptional
    affine types (see [RigConBijection]_, [BijectionLRT]_, [BijectionDn]_) and all
    admit a classical crystal structure [CrysStructSchilling06]_.

    .. NOTE::

        All non-simply-laced rigged configurations have not been proven to
        give rise to aligned virtual crystals (i.e. have the correct crystal
        structure or ismorphic as affine crystals to the tensor product of
        Kirillov-Reshetikhin tableaux).

    INPUT:

    - ``cartan_type`` -- a Cartan type

    - ``B`` -- a list of positive integer tuples `(r,s)` corresponding to the
      tensor factors in the bijection with tensor product of
      Kirillov-Reshetikhin tableaux

    A rigged configuration element with all riggings equal to the vacancy numbers can be created as follows::

        sage: RC = RiggedConfigurations(['A', 3, 1], [[3,2], [2,1], [1,1], [1,1]]); RC
        Rigged configurations of type ['A', 3, 1] and factor(s) ((3, 2), (2, 1), (1, 1), (1, 1))
        sage: elt = RC(partition_list=[[1],[],[]]); elt
        <BLANKLINE>
        0[ ]0
        <BLANKLINE>
        (/)
        <BLANKLINE>
        (/)
        <BLANKLINE>

    If on the other hand we also want to specify the riggings, this can be achieved as follows::

        sage: RC = RiggedConfigurations(['D', 7, 1], [[3,3],[5,2],[4,3],[2,3],[4,4],[3,1],[1,4],[2,2]])
        sage: elt = RC(partition_list=[[2],[3,2,1],[2,2,1,1],[2,2,1,1,1,1],[3,2,1,1,1,1],[2,1,1],[2,2]],
        ....:          rigging_list=[[2],[1,0,0],[4,1,2,1],[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,0],[0,0]])
        sage: elt
        <BLANKLINE>
        3[ ][ ]2
        <BLANKLINE>
        1[ ][ ][ ]1
        2[ ][ ]0
        1[ ]0
        <BLANKLINE>
        4[ ][ ]4
        4[ ][ ]1
        3[ ]2
        3[ ]1
        <BLANKLINE>
        2[ ][ ]1
        2[ ][ ]0
        0[ ]0
        0[ ]0
        0[ ]0
        0[ ]0
        <BLANKLINE>
        0[ ][ ][ ]0
        2[ ][ ]1
        0[ ]0
        0[ ]0
        0[ ]0
        0[ ]0
        <BLANKLINE>
        0[ ][ ]0
        0[ ]0
        0[ ]0
        <BLANKLINE>
        0[ ][ ]0
        0[ ][ ]0
        <BLANKLINE>

    To obtain the Kirillov-Reshetikhin (KR) tableau under the bijection between rigged configurations and KR
    tableaux, we can type the following. This example was checked against Reiho Sakamoto's Mathematica program
    on rigged configurations::

        sage: output = elt.to_tensor_product_of_kirillov_reshetikhin_tableaux(); output
        [[1, 1, 1], [2, 3, 3], [3, 4, -5]] (X) [[1, 1], [2, 2], [3, 3], [5, -6], [6, -5]] (X)
        [[1, 1, 2], [2, 2, 3], [3, 3, 7], [4, 4, -7]] (X) [[1, 1, 1], [2, 2, 2]] (X)
        [[1, 1, 1, 3], [2, 2, 3, 4], [3, 3, 4, 5], [4, 4, 5, 6]] (X) [[1], [2], [3]] (X) [[1, 1, 1, 1]] (X) [[1, 1], [2, 2]]
        sage: elt.to_tensor_product_of_kirillov_reshetikhin_tableaux().to_rigged_configuration() == elt
        True
        sage: output.to_rigged_configuration().to_tensor_product_of_kirillov_reshetikhin_tableaux() == output
        True

    We can also convert between rigged configurations and tensor products of
    Kirillov-Reshetikhin crystals::

        sage: RC = RiggedConfigurations(['D', 4, 1], [[2, 1]])
        sage: elt = RC(partition_list=[[1],[1,1],[1],[1]])
        sage: tp_krc = elt.to_tensor_product_of_kirillov_reshetikhin_crystals(); tp_krc
        [[]]
        sage: ret = RC(tp_krc)
        sage: ret == elt
        True

    ::

        sage: RC = RiggedConfigurations(['D', 4, 1], [[4,1], [3,3]])
        sage: KR1 = KirillovReshetikhinCrystal(['D', 4, 1], 4, 1)
        sage: KR2 = KirillovReshetikhinCrystal(['D', 4, 1], 3, 3)
        sage: T = TensorProductOfCrystals(KR1, KR2)
        sage: t = T[1]; t
        [[++++, []], [+++-, [[1], [2], [4], [-4]]]]
        sage: ret = RC(t)
        sage: ret.to_tensor_product_of_kirillov_reshetikhin_crystals()
        [[++++, []], [+++-, [[1], [2], [4], [-4]]]]
    """
    @staticmethod
    def __classcall_private__(cls, cartan_type, B):
        r"""
        Normalize the input arguments to ensure unique representation.

        EXAMPLES::

            sage: RC1 = RiggedConfigurations(CartanType(['A',3,1]), [[2,2]])
            sage: RC2 = RiggedConfigurations(['A',3,1], [(2,2)])
            sage: RC3 = RiggedConfigurations(['A',3,1], ((2,2),))
            sage: RC2 is RC1, RC3 is RC1
            (True, True)
        """
        cartan_type = CartanType(cartan_type)
        if not cartan_type.is_affine():
            raise ValueError("The Cartan type must be affine")

        # Standardize B input into a tuple of tuples
        B = tuple(tuple(factor) for factor in B)

        if cartan_type.type() == 'BC': # Type `A_{2n}^{(2)}`
            return RCTypeA2Even(cartan_type, B)
        if cartan_type.dual().type() == 'BC': # Type 'A_{2n}^{(2)\dagger`
            return RCTypeA2Dual(cartan_type, B)
        # We check the classical type to account for A^{(1)}_1 which is not
        #    a virtual rigged configuration.
        if not cartan_type.classical().is_simply_laced():
            return RCVirtual(cartan_type, B)

        return super(RiggedConfigurations, cls).__classcall__(cls, cartan_type, B)

    def __init__(self, cartan_type, B):
        r"""
        Initialize the RiggedConfigurations class.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A', 3, 1], [[3, 2], [1, 2], [1, 1]])
            sage: RC(partition_list=[[2],[2],[2]])
            <BLANKLINE>
            1[ ][ ]1
            <BLANKLINE>
            0[ ][ ]0
            <BLANKLINE>
            0[ ][ ]0
            sage: RC(partition_list=[[2],[2],[2]], rigging_list=[[0],[0],[0]])
            <BLANKLINE>
            1[ ][ ]0
            <BLANKLINE>
            0[ ][ ]0
            <BLANKLINE>
            0[ ][ ]0
            sage: RC
            Rigged configurations of type ['A', 3, 1] and factor(s) ((3, 2), (1, 2), (1, 1))
            sage: TestSuite(RC).run()  # long time (4s on sage.math, 2012)
            sage: RC = RiggedConfigurations(['A',1,1], [[1,1], [1,1]])
            sage: TestSuite(RC).run()
            sage: RC = RiggedConfigurations(['A',2,1], [[1,1], [2,1]])
            sage: TestSuite(RC).run()
            sage: RC = RiggedConfigurations(['D', 4, 1], [[2,2]])
            sage: TestSuite(RC).run() # long time
            sage: RC = RiggedConfigurations(['D', 4, 1], [[3,1]])
            sage: TestSuite(RC).run() # long time
            sage: RC = RiggedConfigurations(['D', 4, 1], [[4,3]])
            sage: TestSuite(RC).run() # long time
        """

        self._cartan_type = cartan_type
        self.dims = B
        # We store the cartan matrix for the vacancy number calculations for speed
        self._cartan_matrix = self._cartan_type.classical().cartan_matrix()
        Parent.__init__(self, category=(RegularCrystals(), FiniteCrystals()))
        self.rename("Rigged configurations of type %s and factor(s) %s" % (cartan_type, B))

    def __iter__(self):
        """
        Returns the iterator of ``self``.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A', 3, 1], [[2,1], [1,1]])
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
        index_set = self._cartan_type.classical().index_set()
        from sage.combinat.backtrack import TransitiveIdeal
        return TransitiveIdeal(lambda x: [x.f(i) for i in index_set],
                               self.module_generators).__iter__()

    use_half_width_boxes_type_B = True

    @lazy_attribute
    def module_generators(self):
        r"""
        Module generators for this set of rigged configurations.

        Iterate over the highest weight rigged configurations by moving
        through the :class:`KleberTree` and then setting appropriate
        values of the partitions.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['D', 4, 1], [[2,1]])
            sage: for x in RC.module_generators: x
            <BLANKLINE>
            (/)
            <BLANKLINE>
            (/)
            <BLANKLINE>
            (/)
            <BLANKLINE>
            (/)
            <BLANKLINE>
            <BLANKLINE>
            0[ ]0
            <BLANKLINE>
            0[ ]0
            0[ ]0
            <BLANKLINE>
            0[ ]0
            <BLANKLINE>
            0[ ]0
            <BLANKLINE>
        """
        module_gens = []
        n = self._cartan_type.classical().rank()

        for tree_node in self.kleber_tree():
            shapes = []
            cur = tree_node
            path_lambda = [cur.up_root.to_vector()] # Build the lambda values
            # Note that these are not same lambda as in the paper,
            #   but a less computational version.
            while cur.parent_node is not None:
                path_lambda.insert(0, (cur.parent_node.up_root - cur.up_root).to_vector())
                cur = cur.parent_node

            for a in range(n):
                shapes.append([])
                for i, cur_lambda in enumerate(path_lambda):
                    for j in range(cur_lambda[a]):
                        shapes[-1].insert(0, i)

            # Start with a base to calculate the vacancy numbers
            # Make a copy just to be safe
            base = self.element_class(self, partition_list=shapes[:])

            # Build out the blocks for the partition values
            vac_nums = []
            blocks = []
            L = []

            for partition in base:
                vac_nums.append(partition.vacancy_numbers)
                blocks.append([[]])

                # If the partition is empty, there's nothing to do
                if len(partition) <= 0:
                    L.append([[]])
                    continue

                # Setup the first block
                block_len = partition[0]
                for i, rowLen in enumerate(partition):
                    # If we've gone to a different sized block, then update the
                    #   values which change when moving to a new block size
                    if block_len != rowLen:
                        blocks[-1].append([])
                        block_len = rowLen

                    blocks[-1][-1].append(partition.vacancy_numbers[i])

                L2 = []
                for block in blocks[-1]:
                    L2.append(IterableFunctionCall(self._block_iterator, block))
                L.append(CartesianProduct(*L2))

            # TODO Find a more efficient method without appealing to the CartesianProduct
            C = CartesianProduct(*L)
            for curBlocks in C:
                module_gens.append( self.element_class(self, KT_constructor=[shapes[:],
                                         self._blocks_to_values(curBlocks[:]),
                                         vac_nums[:]]) )

        return tuple(module_gens)

    def _block_iterator(self, container):
        r"""
        Iterate over all possible riggings for a particular block.

        Helper iterator which iterates over all possible partitions contained
        within the container.

        INPUT:

        - ``container`` -- A list of widths of the rows of the container

        TESTS::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2, 2]])
            sage: for x in RC._block_iterator([]): x
            []
            sage: for x in RC._block_iterator([2,3]): x
            [0, 0]
            [1, 0]
            [1, 1]
            [2, 0]
            [2, 1]
            [2, 2]
        """
        if len(container) == 0:
            yield []
            return

        pos = 0
        length = len(container)
        ret_part = [-1] * length
        while pos >= 0:
            ret_part[pos] += 1

            if ret_part[pos] > container[pos] or (pos != 0 and ret_part[pos] > ret_part[pos - 1]):
                ret_part[pos] = -1
                pos -= 1
            else:
                pos += 1

            if pos == length:
                yield ret_part[:]
                pos -= 1

    def _blocks_to_values(self, blocks):
        r"""
        Convert an array of blocks into a list of partition values.

        INPUT:

        - ``blocks`` -- The (2-dim) array blocks of the partition values.

        TESTS::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2, 2]])
            sage: RC._blocks_to_values([[[2, 1]]])
            [[2, 1]]
        """
        values = []
        for part_block in blocks:
            if len(part_block) == 0:
                values.append([])
            else:
                values.append(part_block[0][:]) # Need to make a copy
                for block in part_block[1:]:
                    values[-1].extend(block)
        return values

    def _element_constructor_(self, *lst, **options):
        """
        Construct a ``RiggedConfigurationElement``.

        Typically the user should not call this method since it does not check
        if it is a valid configuration. Instead the user should use the
        iterator methods.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2, 1]])
            sage: RC(partition_list=[[1], [1], [], []], rigging_list=[[-1], [0], [], []]) # indirect doctest
            <BLANKLINE>
            -1[ ]-1
            <BLANKLINE>
            0[ ]0
            <BLANKLINE>
            (/)
            <BLANKLINE>
            (/)
            <BLANKLINE>
        """
        from sage.combinat.rigged_configurations.tensor_product_kr_tableaux_element import TensorProductOfKirillovReshetikhinTableauxElement
        if isinstance(lst[0], TensorProductOfKirillovReshetikhinTableauxElement):
            if self != lst[0].parent().rigged_configurations():
                raise ValueError("incorrect bijection image")
            return lst[0].to_rigged_configuration()

        from sage.combinat.crystals.tensor_product import TensorProductOfRegularCrystalsElement
        if isinstance(lst[0], TensorProductOfRegularCrystalsElement):
            lst = lst[0]
        from sage.combinat.crystals.kirillov_reshetikhin import KirillovReshetikhinGenericCrystalElement
        if isinstance(lst[0], KirillovReshetikhinGenericCrystalElement):
            KRT = self.tensor_product_of_kirillov_reshetikhin_tableaux()
            krt_elt = KRT(*[x.to_kirillov_reshetikhin_tableau() for x in lst])
            return krt_elt.to_rigged_configuration()

        return self.element_class(self, list(lst), **options)

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

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2, 1]])
            sage: elt = RC(partition_list=[[1], [1], [], []])
            sage: RC._calc_vacancy_number(elt.nu(), 1, 0)
            0
        """
        row_len = partitions[a][i]

        vac_num = 0
        if "B" in options:
            for tableau in options["B"]:
                # The + 1 is to convert between indexing methods
                if len(tableau) == a + 1:
                    vac_num += min(row_len, len(tableau[0]))
        elif "L" in options:
            L = options["L"]
            if L.has_key(a):
                for kvp in L[a].items():
                    vac_num += min(kvp[0], row_len) * kvp[1]
        elif "dims" in options:
            for dim in options["dims"]:
                if dim[0] == a + 1:
                    vac_num += min(dim[1], row_len)
        else:
            for dim in self.dims:
                if dim[0] == a + 1:
                    vac_num += min(dim[1], row_len)

        for b, value in enumerate(self._cartan_matrix.row(a)):
            vac_num -= value * partitions[b].get_num_cells_to_column(row_len)

        return vac_num

    def kleber_tree(self):
        r"""
        Return the underlying Kleber tree used to generate all highest
        weight rigged configurations.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A',3,1], [[1,1], [2,1]])
            sage: RC.kleber_tree()
            Kleber tree of Cartan type ['A', 3, 1] and B = ((1, 1), (2, 1))
        """
        return KleberTree(self._cartan_type, self.dims)

    @cached_method
    def tensor_product_of_kirillov_reshetikhin_tableaux(self):
        """
        Return the corresponding tensor product of Kirillov-Reshetikhin
        tableaux.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A', 3, 1], [[3, 2], [1, 2]])
            sage: RC.tensor_product_of_kirillov_reshetikhin_tableaux()
            Tensor product of Kirillov-Reshetikhin tableaux of type ['A', 3, 1] and factor(s) ((3, 2), (1, 2))
        """
        from sage.combinat.rigged_configurations.tensor_product_kr_tableaux import TensorProductOfKirillovReshetikhinTableaux
        return TensorProductOfKirillovReshetikhinTableaux(self._cartan_type, self.dims)

    def cardinality(self):
        """
        Return the cardinality of ``self``.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A', 3, 1], [[3, 2], [1, 2]])
            sage: RC.cardinality()
            100
            sage: RC = RiggedConfigurations(['B', 3, 1], [[2,2],[1,2]])
            sage: RC.cardinality()
            5130
            sage: RC = RiggedConfigurations(['E', 7, 1], [[1,1]])
            sage: RC.cardinality()
            134
        """
        try:
            return self.tensor_product_of_kirillov_reshetikhin_tableaux().cardinality()
        except NotImplementedError: # Backup by direct counting
            return ZZ(sum(1 for x in self))

    @cached_method
    def tensor_product_of_kirillov_reshetikhin_crystals(self):
        """
        Return the corresponding tensor product of Kirillov-Reshetikhin
        crystals.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A', 3, 1], [[3,1],[2,2]])
            sage: RC.tensor_product_of_kirillov_reshetikhin_crystals()
            Full tensor product of the crystals
            [Kirillov-Reshetikhin crystal of type ['A', 3, 1] with (r,s)=(3,1),
             Kirillov-Reshetikhin crystal of type ['A', 3, 1] with (r,s)=(2,2)]
        """
        return self.tensor_product_of_kirillov_reshetikhin_tableaux().tensor_product_of_kirillov_reshetikhin_crystals()

    def fermionic_formula(self, q=None):
        r"""
        Return the fermoinic formula associated to ``self``.

        Given a set of rigged configurations `RC(\Lambda, L)`, the fermonic
        formula is defined as:

        .. MATH::

            M(\lambda, L; q) = \sum_{(\nu,J)} q^{cc(\nu, J)}

        where we sum over all classically highest weight rigged
        configurations where `cc` is the :meth:`cocharge statistic
        <sage.combinat.rigged_configurations.rigged_configuration_element.RiggedConfigurationElement.cc>`.
        This is known to reduce to

        .. MATH::

            M(\lambda, L; q) = \sum_{\nu} q^{cc(\nu)} \prod_{(a,i) \in
            I \times \ZZ} \begin{bmatrix} p_i^{(a)} + m_i^{(a)} \\ m_i^{(a)}
            \end{bmatrix}_q.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A', 3, 1], [[3,1],[2,2]])
            sage: RC.fermionic_formula()
            q + 1
            sage: RC = RiggedConfigurations(['D', 4, 1], [[3,2],[4,1],[2,2]])
            sage: RC.fermionic_formula()
            q^6 + 6*q^5 + 11*q^4 + 14*q^3 + 11*q^2 + 4*q + 1
            sage: RC = RiggedConfigurations(['E', 6, 1], [[2,2]])
            sage: RC.fermionic_formula()
            q^2 + q + 1
            sage: RC = RiggedConfigurations(['B', 3, 1], [[3,1],[2,2]])
            sage: RC.fermionic_formula()
            q^3 + 3*q^2 + 2*q + 1
            sage: RC = RiggedConfigurations(['C', 3, 1], [[3,1],[2,2]])
            sage: RC.fermionic_formula()
            q^4 + 2*q^3 + 3*q^2 + 2*q + 1
            sage: RC = RiggedConfigurations(['D', 4, 2], [[3,1],[2,2]])
            sage: RC.fermionic_formula()
            q^6 + 2*q^5 + 4*q^4 + 3*q^3 + 3*q^2 + q + 1
        """
        if q is None:
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            q = PolynomialRing(ZZ, 'q').gen(0)
        return sum([q**x.cc() for x in self.module_generators], ZZ.zero())

    def _test_bijection(self, **options):
        r"""
        Test function to make sure that the bijection between rigged
        configurations and Kirillov-Reshetikhin tableaux is correct.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A', 3, 1], [[2,1],[1,1]])
            sage: RC._test_bijection()
        """
        tester = self._tester(**options)
        rejects = []
        for x in self:
            y = x.to_tensor_product_of_kirillov_reshetikhin_tableaux()
            z = y.to_rigged_configuration()
            if z != x:
                rejects.append((x, z))

        tester.assertTrue(len(rejects) == 0, "Bijection is not correct: %s"%rejects)
        if len(rejects) != 0:
            return rejects

RiggedConfigurations.Element = RiggedConfigurationElement

class RCVirtual(RiggedConfigurations):
    r"""
    Class of virtual rigged configurations.

    These are rigged configurations which lift to a virtual configuration.

    For more on rigged configurations, see :class:`RiggedConfigurations`.
    """
    @staticmethod
    def __classcall_private__(cls, cartan_type, B):
        r"""
        Normalize the input arguments to ensure unique representation.

        EXAMPLES::

            sage: RC1 = RiggedConfigurations(CartanType(['A',4,2]), [[2,2]])
            sage: RC2 = RiggedConfigurations(['A',4,2], [(2,2)])
            sage: RC3 = RiggedConfigurations(['BC',2,2], ((2,2),))
            sage: RC2 is RC1, RC3 is RC1
            (True, True)
        """
        cartan_type = CartanType(cartan_type)

        # Standardize B input into a tuple of tuples
        B = tuple(map(tuple, B))
        return super(RCVirtual, cls).__classcall__(cls, cartan_type, B)

    def __init__(self, cartan_type, dims):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['C',2,1], [[1,1]])
            sage: TestSuite(RC).run()
            sage: RC = RiggedConfigurations(['C',2,1], [[1,2],[1,1],[2,1]]); RC
            Rigged configurations of type ['C', 2, 1] and factor(s) ((1, 2), (1, 1), (2, 1))
            sage: TestSuite(RC).run() # long time
            sage: RC = RiggedConfigurations(['B',3,1], [[3,1],[1,1]])
            sage: TestSuite(RC).run() # long time
            sage: RC = RiggedConfigurations(['D',4,2], [[2,1]])
            sage: TestSuite(RC).run() # long time
            sage: RC = RiggedConfigurations(['A',5,2], [[2,1]])
            sage: TestSuite(RC).run() # long time
         """
        self._folded_ct = cartan_type.as_folding()
        RiggedConfigurations.__init__(self, cartan_type, dims)

    def _calc_vacancy_number(self, partitions, a, i, **options):
        r"""
        Calculate the vacancy number of the `i`-th row of the `a`-th rigged
        partition.

        INPUT:

        - ``partitions`` -- The list of rigged partitions we are using

        - ``a``          -- The rigged partition index

        - ``i``          -- The row index of the `a`-th rigged partition

        TESTS::

            sage: RC = RiggedConfigurations(['C', 4, 1], [[2, 1]])
            sage: elt = RC(partition_list=[[1], [2], [2], [1]])
            sage: RC._calc_vacancy_number(elt.nu(), 1, 0)
            0
        """
        row_len = partitions[a][i]

        vac_num = 0
        if "B" in options:
            for tableau in options["B"]:
                # The + 1 is to convert between indexing methods
                if len(tableau) == a + 1:
                    vac_num += min(row_len, len(tableau[0]))
        elif "L" in options:
            L = options["L"]
            if L.has_key(a):
                for kvp in L[a].items():
                    vac_num += min(kvp[0], row_len) * kvp[1]
        elif "dims" in options:
            for dim in options["dims"]:
                if dim[0] == a + 1:
                    vac_num += min(dim[1], row_len)
        else:
            for dim in self.dims:
                if dim[0] == a + 1:
                    vac_num += min(dim[1], row_len)

        gamma = self._folded_ct.scaling_factors()
        for b, value in enumerate(self._cartan_matrix.row(a)):
            vac_num -= value * partitions[b].get_num_cells_to_column(gamma[a+1]*row_len, gamma[b+1]) // gamma[b+1]

        return vac_num

    @lazy_attribute
    def module_generators(self):
        r"""
        Module generators for this set of rigged configurations.

        Iterate over the highest weight rigged configurations by moving
        through the :class:`KleberTree` and then setting appropriate
        values of the partitions.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['C', 3, 1], [[1,2]])
            sage: for x in RC.module_generators: x
            <BLANKLINE>
            (/)
            <BLANKLINE>
            (/)
            <BLANKLINE>
            (/)
            <BLANKLINE>
            <BLANKLINE>
            0[ ][ ]0
            <BLANKLINE>
            0[ ][ ]0
            <BLANKLINE>
            0[ ]0
            <BLANKLINE>
        """
        module_gens = []
        vec_len = len(self.kleber_tree().root.up_root.to_vector())

        for tree_node in self.kleber_tree():
            shapes = []
            cur = tree_node
            path_lambda = [cur.up_root.to_vector()] # Build the lambda values
            # Note that these are not same lambda as in the paper,
            #   but a less computational version.
            while cur.parent_node is not None:
                path_lambda.insert(0, (cur.parent_node.up_root - cur.up_root).to_vector())
                cur = cur.parent_node

            for a in range(vec_len):
                shapes.append([])
                for i, cur_lambda in enumerate(path_lambda):
                    for j in range(cur_lambda[a]):
                        shapes[-1].insert(0, i)

            # Convert from the virtual rigged configuration
            # As a special case, we do not need to do anything for type `A_{2n}^{(2)}`
            sigma = self._folded_ct.folding_orbit()
            shapes = shapes[:len(sigma)-1]
            if self._cartan_type.type() != 'BC':
                gammatilde = self._folded_ct.scaling_factors()
                for a in range(1, len(sigma)):
                    index = sigma[a][0] - 1 # for indexing
                    for i in range(len(shapes[index])):
                        shapes[index][i] = shapes[index][i] // gammatilde[a]

            # Start with a base to calculate the vacancy numbers
            # Make a copy just to be safe
            base = self.element_class(self, partition_list=shapes[:])

            # Build out the blocks for the partition values
            vac_nums = []
            blocks = []
            L = []

            for partition in base:
                vac_nums.append(partition.vacancy_numbers)
                blocks.append([[]])

                # If the partition is empty, there's nothing to do
                if len(partition) <= 0:
                    L.append([[]])
                    continue

                # Setup the first block
                block_len = partition[0]
                for i, rowLen in enumerate(partition):
                    # If we've gone to a different sized block, then update the
                    #   values which change when moving to a new block size
                    if block_len != rowLen:
                        blocks[-1].append([])
                        block_len = rowLen

                    blocks[-1][-1].append(partition.vacancy_numbers[i])

                L2 = []
                for block in blocks[-1]:
                    L2.append(IterableFunctionCall(self._block_iterator, block))
                L.append(CartesianProduct(*L2))

            # TODO Find a more efficient method without appealing to the CartesianProduct
            C = CartesianProduct(*L)
            for curBlocks in C:
                module_gens.append( self.element_class(self, KT_constructor=[shapes[:],
                                         self._blocks_to_values(curBlocks[:]), vac_nums[:]]) )

        return tuple(module_gens)

    def kleber_tree(self):
        r"""
        Return the underlying (virtual) Kleber tree used to generate all
        highest weight rigged configurations.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['C',3,1], [[1,1], [2,1]])
            sage: RC.kleber_tree()
            Virtual Kleber tree of Cartan type ['C', 3, 1] and B = ((1, 1), (2, 1))
        """
        return VirtualKleberTree(self._cartan_type, self.dims)

    @lazy_attribute
    def virtual(self):
        """
        Return the corresponding virtual crystal.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['C',2,1], [[1,2],[1,1],[2,1]])
            sage: RC
            Rigged configurations of type ['C', 2, 1] and factor(s) ((1, 2), (1, 1), (2, 1))
            sage: RC.virtual
            Rigged configurations of type ['A', 3, 1] and factor(s) ((1, 2), (3, 2), (1, 1), (3, 1), (2, 2))
        """
        gamma = self._folded_ct.scaling_factors()
        sigma = self._folded_ct.folding_orbit()
        virtual_dims = []
        for r,s in self.dims:
            for a in sigma[r]:
                virtual_dims.append([a, s*gamma[r]])
        return RiggedConfigurations(self._folded_ct._folding, virtual_dims)

    def to_virtual(self, rc):
        """
        Convert ``rc`` into a rigged configuration in the virtual crystal.

        INPUT:

        - ``rc`` -- A rigged configuration element

        EXAMPLES::

            sage: RC = RiggedConfigurations(['C',2,1], [[1,2],[1,1],[2,1]])
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
            Rigged configurations of type ['A', 3, 1] and factor(s) ((1, 2), (3, 2), (1, 1), (3, 1), (2, 2))
        """
        gamma = self._folded_ct.scaling_factors()
        sigma = self._folded_ct.folding_orbit()
        n = self._folded_ct._folding.classical().rank()
        partitions = [None] * n
        riggings = [None] * n
        vac_nums = [None] * n
        # +/- 1 for indexing
        for a in range(len(rc)):
            for i in sigma[a+1]:
                partitions[i-1] = [row_len*gamma[a+1] for row_len in rc[a]._list]
                riggings[i-1] = [rig_val*gamma[a+1] for rig_val in rc[a].rigging]
                vac_nums[i-1] = [vac_num*gamma[a+1] for vac_num in rc[a].vacancy_numbers]
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

            sage: RC = RiggedConfigurations(['C',2,1], [[1,2],[1,1],[2,1]])
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
        gamma = self._folded_ct.scaling_factors()
        sigma = self._folded_ct.folding_orbit()
        n = self._cartan_type.classical().rank()
        partitions = [None] * n
        riggings = [None] * n
        vac_nums = [None] * n
        # +/- 1 for indexing
        for a in range(n):
            index = sigma[a+1][0]-1
            partitions[a] = [row_len//gamma[a+1] for row_len in vrc[index]._list]
            riggings[a] = [rig_val//gamma[a+1] for rig_val in vrc[index].rigging]
            vac_nums[a] = [vac_val//gamma[a+1] for vac_val in vrc[index].vacancy_numbers]
        return self.element_class(self, partition_list=partitions,
                                  rigging_list=riggings, vacancy_numbers_list=vac_nums)

    def _test_virtual_vacancy_numbers(self, **options):
        """
        Test to make sure that the vacancy numbers obtained from the virtual
        rigged configuration agree with the explicit computation of the
        vacancy numbers done here.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['B', 3, 1], [[2,1]])
            sage: RC._test_virtual_vacancy_numbers()
        """
        tester = self._tester(**options)
        for x in self:
            parts_list = [p._list[:] for p in x]
            elt = self.element_class(self, partition_list=parts_list)
            for i, p in enumerate(elt):
                for j, vac_num in enumerate(p.vacancy_numbers):
                    tester.assertTrue(vac_num == x[i].vacancy_numbers[j],
                      "Incorrect vacancy number: %s\nComputed: %s\nFor: %s"\
                      %(x[i].vacancy_numbers[j],vac_num, x))

RCVirtual.Element = RCVirtualElement

class RCTypeA2Even(RCVirtual):
    """
    Rigged configurations for type `A_{2n}^{(2)}`.

    For more on rigged configurations, see :class:`RiggedConfigurations`.

    EXAMPLES::

        sage: RC = RiggedConfigurations(['A',4,2], [[2,1], [1,2]])
        sage: RC.cardinality()
        150
        sage: RC = RiggedConfigurations(['A',2,2], [[1,1]])
        sage: RC.cardinality()
        3
        sage: RC = RiggedConfigurations(['A',2,2], [[1,2],[1,1]])
        sage: TestSuite(RC).run()
        sage: RC = RiggedConfigurations(['A',4,2], [[2,1]])
        sage: TestSuite(RC).run() # long time
    """
    @lazy_attribute
    def virtual(self):
        """
        Return the corresponding virtual crystal.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A',4,2], [[1,2],[1,1],[2,1]])
            sage: RC
            Rigged configurations of type ['BC', 2, 2] and factor(s) ((1, 2), (1, 1), (2, 1))
            sage: RC.virtual
            Rigged configurations of type ['A', 3, 1] and factor(s) ((1, 2), (3, 2), (1, 1), (3, 1), (2, 1), (2, 1))
        """
        sigma = self._folded_ct.folding_orbit()
        n = len(sigma) - 1
        virtual_dims = []
        for r,s in self.dims:
            if r == n:
                virtual_dims.extend([[n, s], [n, s]])
            else:
                for a in sigma[r]:
                    virtual_dims.append([a, s])
        return RiggedConfigurations(self._folded_ct._folding, virtual_dims)

    def _calc_vacancy_number(self, partitions, a, i, **options):
        r"""
        Calculate the vacancy number of the `i`-th row of the `a`-th rigged
        partition.

        This is a special implementation for type `A_{2n}^{(2)}`.

        INPUT:

        - ``partitions`` -- The list of rigged partitions we are using

        - ``a``          -- The rigged partition index

        - ``i``          -- The row index of the `a`-th rigged partition

        TESTS::

            sage: RC = RiggedConfigurations(['A', 4, 2], [[2, 1]])
            sage: elt = RC(partition_list=[[1], [2]])
            sage: RC._calc_vacancy_number(elt.nu(), 1, 0)
            0
        """
        row_len = partitions[a][i]

        vac_num = 0
        if "B" in options:
            for tableau in options["B"]:
                # The + 1 is to convert between indexing methods
                if len(tableau) == a + 1:
                    vac_num += min(row_len, len(tableau[0]))
        elif "L" in options:
            L = options["L"]
            if L.has_key(a):
                for kvp in L[a].items():
                    vac_num += min(kvp[0], row_len) * kvp[1]
        elif "dims" in options:
            for dim in options["dims"]:
                if dim[0] == a + 1:
                    vac_num += min(dim[1], row_len)
        else:
            for dim in self.dims:
                if dim[0] == a + 1:
                    vac_num += min(dim[1], row_len)

        gamma = self._folded_ct.scaling_factors()
        for b, value in enumerate(self._cartan_matrix.row(a)):
            vac_num -= value * partitions[b].get_num_cells_to_column(row_len) // gamma[b+1]

        return vac_num

    def to_virtual(self, rc):
        """
        Convert ``rc`` into a rigged configuration in the virtual crystal.

        INPUT:

        - ``rc`` -- A rigged configuration element

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A',4,2], [[2,2]])
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
        n = self._folded_ct._folding.classical().rank()
        partitions = [None] * n
        riggings = [None] * n
        vac_nums = [None] * n
        # +/- 1 for indexing
        for a in range(len(rc)):
            for i in sigma[a+1]:
                partitions[i-1] = [row_len for row_len in rc[a]._list]
                riggings[i-1] = [rig_val*gamma[a+1] for rig_val in rc[a].rigging]
                vac_nums[i-1] = [vac_num*gamma[a+1] for vac_num in rc[a].vacancy_numbers]
        return self.virtual.element_class(self.virtual, partition_list=partitions,
                            rigging_list=riggings,
                            vacancy_numbers_list=vac_nums)

    def from_virtual(self, vrc):
        """
        Convert ``vrc`` in the virtual crystal into a rigged configution of
        the original Cartan type.

        INPUT:

        - ``vrc`` -- A virtual rigged configuration element

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A',4,2], [[2,2]])
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
        n = self._cartan_type.classical().rank()
        partitions = [None] * n
        riggings = [None] * n
        vac_nums = [None] * n
        # +/- 1 for indexing
        for a in range(n):
            index = sigma[a+1][0]-1
            partitions[index] = [row_len for row_len in vrc[index]._list]
            riggings[index] = [rig_val//gamma[a+1] for rig_val in vrc[index].rigging]
            vac_nums[a] = [vac_val//gamma[a+1] for vac_val in vrc[index].vacancy_numbers]
        return self.element_class(self, partition_list=partitions,
                                  rigging_list=riggings, vacancy_numbers_list=vac_nums)

class RCTypeA2Dual(RCTypeA2Even):
    r"""
    Rigged configurations of type `A_{2n}^{(2)\dagger}`.

    For more on rigged configurations, see :class:`RiggedConfigurations`.

    EXAMPLES::

        sage: RC = RiggedConfigurations(CartanType(['A',4,2]).dual(), [[1,2],[1,1],[2,1]])
        sage: RC
        Rigged configurations of type ['BC', 2, 2]^* and factor(s) ((1, 2), (1, 1), (2, 1))
        sage: RC.cardinality()
        750
        sage: RC.virtual
        Rigged configurations of type ['A', 3, 1] and factor(s) ((1, 2), (3, 2), (1, 1), (3, 1), (2, 1), (2, 1))
        sage: RC = RiggedConfigurations(CartanType(['A',2,2]).dual(), [[1,1]])
        sage: RC.cardinality()
        3
        sage: RC = RiggedConfigurations(CartanType(['A',2,2]).dual(), [[1,2],[1,1]])
        sage: TestSuite(RC).run()
        sage: RC = RiggedConfigurations(CartanType(['A',4,2]).dual(), [[2,1]])
        sage: TestSuite(RC).run() # long time
    """
    def _calc_vacancy_number(self, partitions, a, i, **options):
        r"""
        Calculate the vacancy number. A special case is needed for the `n`-th
        partition for type `A_{2n}^{(2)\dagger}`.

        INPUT:

        - ``partitions`` -- The list of rigged partitions we are using

        - ``a``          -- The rigged partition index

        - ``i``          -- The row index of the `a`-th rigged partition

        TESTS::

            sage: RC = RiggedConfigurations(CartanType(['A', 6, 2]).dual(), [[2,1]])
            sage: elt = RC(partition_list=[[1], [2], [2]])
            sage: RC._calc_vacancy_number(elt.nu(), 0, 0)
            -1
        """
        if a != self._cartan_type.classical().rank()-1:
            return RCTypeA2Even._calc_vacancy_number(self, partitions, a, i, **options)
        row_len = partitions[a][i]

        vac_num = 0
        if "B" in options:
            for tableau in options["B"]:
                # The + 1 is to convert between indexing methods
                if len(tableau) == a + 1:
                    vac_num += min(row_len, len(tableau[0]))
        elif "L" in options:
            L = options["L"]
            if L.has_key(a):
                for kvp in L[a].items():
                    vac_num += min(kvp[0], row_len) * kvp[1]
        elif "dims" in options:
            for dim in options["dims"]:
                if dim[0] == a + 1:
                    vac_num += min(dim[1], row_len)
        else:
            for dim in self.dims:
                if dim[0] == a + 1:
                    vac_num += min(dim[1], row_len)

        for b, value in enumerate(self._cartan_matrix.row(a)):
            vac_num -= value * partitions[b].get_num_cells_to_column(row_len) / 2

        return vac_num

    @lazy_attribute
    def module_generators(self):
        r"""
        Module generators for rigged configurations of type
        `A_{2n}^{(2)\dagger}`.

        Iterate over the highest weight rigged configurations by moving
        through the :class:`KleberTree` and then setting appropriate
        values of the partitions. This also skips rigged configurations where
        `P_i^{(n)} < 1` when `i` is odd.

        EXAMPLES::

            sage: RC = RiggedConfigurations(CartanType(['A', 4, 2]).dual(), [[1,1]])
            sage: for x in RC.module_generators: x
            <BLANKLINE>
            (/)
            <BLANKLINE>
            (/)
            <BLANKLINE>
        """
        module_gens = []
        # This is for the non-simply-laced types
        vec_len = len(self.kleber_tree().root.up_root.to_vector())

        for tree_node in self.kleber_tree():
            shapes = []
            cur = tree_node
            path_lambda = [cur.up_root.to_vector()] # Build the lambda values
            # Note that these are not same lambda as in the paper,
            #   but a less computational version.
            while cur.parent_node is not None:
                path_lambda.insert(0, (cur.parent_node.up_root - cur.up_root).to_vector())
                cur = cur.parent_node

            for a in range(vec_len):
                shapes.append([])
                for i, cur_lambda in enumerate(path_lambda):
                    for j in range(cur_lambda[a]):
                        shapes[-1].insert(0, i)

            # We are not simply-laced, so convert from the virtual rigged configuration
            sigma = self._folded_ct.folding_orbit()
            shapes = shapes[:len(sigma)-1]
            # Nothing more to do since gamma[i] == 1 for all i >= 1

            # Start with a base to calculate the vacancy numbers
            # Make a copy just to be safe
            base = self.element_class(self, partition_list=shapes[:])

            # Check the special condition of odd rows in the n-th partition
            invalid_RC = False
            for i in range(len(base[-1]._list)):
                if base[-1]._list[i] % 2 == 1 and base[-1].vacancy_numbers[i] < 1:
                    invalid_RC = True
                    break
            # If it is invalid, skip it
            if invalid_RC:
                continue

            # Build out the blocks for the partition values
            vac_nums = []
            blocks = []
            L = []

            for partition in base[:-1]:
                vac_nums.append(partition.vacancy_numbers)
                blocks.append([[]])

                # If the partition is empty, there's nothing to do
                if len(partition) <= 0:
                    L.append([[]])
                    continue

                # Setup the first block
                block_len = partition[0]
                for i, row_len in enumerate(partition):
                    # If we've gone to a different sized block, then update the
                    #   values which change when moving to a new block size
                    if block_len != row_len:
                        blocks[-1].append([])
                        block_len = row_len

                    blocks[-1][-1].append(partition.vacancy_numbers[i])

                L2 = []
                for block in blocks[-1]:
                    L2.append(IterableFunctionCall(self._block_iterator, block))
                L.append(CartesianProduct(*L2))

            # Special case for the final tableau
            partition = base[-1]
            vac_nums.append(partition.vacancy_numbers)

            # If the partition is empty, there's nothing to do
            if len(partition) <= 0:
                L.append([[]])
            else:
                # Setup the first block
                block_len = partition[0]
                blocks = [[]]
                odd_block = []
                for i, row_len in enumerate(partition):
                    # If we've gone to a different sized block, then update the
                    #   values which change when moving to a new block size
                    if block_len != row_len:
                        blocks.append([])
                        odd_block.append(block_len % 2 == 1)
                        block_len = row_len

                    blocks[-1].append(partition.vacancy_numbers[i])
                odd_block.append(block_len % 2 == 1)

                L2 = []
                for i, block in enumerate(blocks):
                    if odd_block[i]:
                        L2.append(IterableFunctionCall(self._block_iterator_n_odd, block))
                    else:
                        L2.append(IterableFunctionCall(self._block_iterator, block))
                L.append(CartesianProduct(*L2))

            # TODO Find a more efficient method without appealing to the CartesianProduct
            C = CartesianProduct(*L)
            for curBlocks in C:
                module_gens.append( self.element_class(self, KT_constructor=[shapes[:],
                                         self._blocks_to_values(curBlocks[:]), vac_nums[:]]) )

        return tuple(module_gens)

    def _block_iterator_n_odd(self, container):
        r"""
        Iterate over all possible riggings for a block of odd length in the
        `n`-th rigged partition for type `A_{2n}^{(2)\dagger}`.

        Helper iterator which iterates over all possible partitions of
        `\frac{2k+1}{2}` sizes contained within the container.

        INPUT:

        - ``container`` -- A list the widths of the rows of the container

        TESTS::

            sage: RC = RiggedConfigurations(CartanType(['A', 4, 2]).dual(), [[2, 2]])
            sage: for x in RC._block_iterator_n_odd([]): x
            []
            sage: for x in RC._block_iterator_n_odd([2,2]): x
            [1/2, 1/2]
            [3/2, 1/2]
            [3/2, 3/2]
        """
        if len(container) == 0:
            yield []
            return

        half = lambda x: QQ(x) / QQ(2)
        pos = 0
        length = len(container)
        ret_part = [-1] * length
        while pos >= 0:
            ret_part[pos] += 2

            if ret_part[pos] > container[pos]*2 or (pos != 0 and ret_part[pos] > ret_part[pos - 1]):
                ret_part[pos] = -1
                pos -= 1
            else:
                pos += 1

            if pos == length:
                yield map(half, ret_part[:])
                pos -= 1

    def to_virtual(self, rc):
        """
        Convert ``rc`` into a rigged configuration in the virtual crystal.

        INPUT:

        - ``rc`` -- A rigged configuration element

        EXAMPLES::

            sage: RC = RiggedConfigurations(CartanType(['A',4,2]).dual(), [[2,2]])
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
        gammatilde = list(self._folded_ct.scaling_factors())
        gammatilde[-1] = 2
        sigma = self._folded_ct.folding_orbit()
        n = self._folded_ct._folding.classical().rank()
        partitions = [None] * n
        riggings = [None] * n
        vac_nums = [None] * n
        # +/- 1 for indexing
        for a in range(len(rc)):
            for i in sigma[a+1]:
                partitions[i-1] = [row_len for row_len in rc[a]._list]
                riggings[i-1] = [rig_val*gammatilde[a+1] for rig_val in rc[a].rigging]
                vac_nums[i-1] = [vac_num for vac_num in rc[a].vacancy_numbers]
        return self.virtual.element_class(self.virtual, partition_list=partitions,
                            rigging_list=riggings,
                            vacancy_numbers_list=vac_nums)

    def from_virtual(self, vrc):
        """
        Convert ``vrc`` in the virtual crystal into a rigged configution of
        the original Cartan type.

        INPUT:

        - ``vrc`` -- A virtual rigged configuration element

        EXAMPLES::

            sage: RC = RiggedConfigurations(CartanType(['A',4,2]).dual(), [[2,2]])
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
        gammatilde = list(self._folded_ct.scaling_factors())
        gammatilde[-1] = QQ(2)
        sigma = self._folded_ct.folding_orbit()
        n = self._cartan_type.classical().rank()
        partitions = [None] * n
        riggings = [None] * n
        vac_nums = [None] * n
        # +/- 1 for indexing
        for a in range(n):
            index = sigma[a+1][0]-1
            partitions[index] = [row_len for row_len in vrc[index]._list]
            riggings[index] = [rig_val/gammatilde[a+1] for rig_val in vrc[index].rigging]
            vac_nums[a] = [vac_val for vac_val in vrc[index].vacancy_numbers]
        return self.element_class(self, partition_list=partitions,
                                  rigging_list=riggings, vacancy_numbers_list=vac_nums)

# For experimentation purposes only.
# I'm keeping this for when we implement the R-matrix in general
def R_matrix(ct, RS1, RS2, only_highest_weight=False):
    r"""
    Output pairs of Kirillov-Reshetikhin tableaux under the action of the
    (combinatorial) `R`-matrix.

    INPUT:

    - ``ct`` -- A Cartan type
    - ``RS1``, ``RS2`` -- Pairs of `(r_i, s_i)` for `i = 1,2` for the two
      tensor factors `B^{r_1, s_1} \otimes B^{r_2, s_2}`
    - ``only_highest_weight`` -- Output only the highest weight elements

    EXAMPLES::

        sage: from sage.combinat.rigged_configurations.rigged_configurations import R_matrix
        sage: L = R_matrix(['D',4,1], [2,1], [1,1])
        sage: len(L)
        232

    Check that the `R`-matrix is the identity on `B \otimes B`::

        sage: L = R_matrix(['A',2,1], [2,3], [2,3])
        sage: len(L)
        100
        sage: all(x == y for x,y in L)
        True
    """
    ct = CartanType(ct)
    RC = RiggedConfigurations(ct, [RS1, RS2])
    RC2 = RiggedConfigurations(ct, [RS2, RS1])
    ret_list = []
    if only_highest_weight:
        L = RC.module_generators
    else:
        L = RC
    for x in L:
        x2 = RC2(*x)
        ret_list.append([x.to_tensor_product_of_kirillov_reshetikhin_tableaux(),
                         x2.to_tensor_product_of_kirillov_reshetikhin_tableaux()])
    return ret_list

