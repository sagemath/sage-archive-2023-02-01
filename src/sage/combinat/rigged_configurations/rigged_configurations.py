r"""
Rigged configurations

AUTHORS:

- Travis Scrimshaw (2010-09-26): Initial version

Rigged configurations form combinatorial objects first introduced by Kirillov and Reshetikhin
that arose from studies of statistical mechanical models using the Bethe Ansatz.
They are sequences of rigged partitions. A rigged partition is a partition together
with a label associated to each part that satisfy certain constraints. The labels
are also called riggings.

Rigged configurations exist for all affine Kac-Moody Lie algebras. See for example
[HKOTT2002]_. In Sage they are specified by providing a Cartan type and a list of
rectangular shapes `R`. Currently, they are only implemented for types `A_n^{(1)}` and
`D_n^{(1)}`. The list of all (highest weight) rigged configurations for given `R` is computed
via the Kleber algorithm (see also :class:`KleberTree`).

There exists a crystal structure on the set of rigged configurations, see [CrysStructSchilling06]_.
The highest weight rigged configurations are those where all riggings are nonnegative.
The list of all rigged configurations is computed from the highest weight ones using the
crystal operators.

Rigged configurations are in bijection with tensor products of Kirillov-Reshetikhin tableaux,
see [RigConBijection]_, [BijectionDn]_, and [BijectionLRT]_. The list of rectangles `R` corresponds
to the shape of the Kirillov-Reshetikhin crystals appearing in the tensor product.
Kirillov-Reshetikhin crystals are implemented in Sage, see :class:`KirillovReshetikhinCrystal`, however,
in the bijection with rigged configurations a different realization of the elements in the
crystal are obtained, which are coined Kirillov-Reshetkihin tableaux, see :class:`KirillovReshetikhinTableaux`.
For more details see [AffineRigConDn]_.

INPUT:

- ``cartan_type`` -- a Cartan type (currently only types `A_n^{(1)}` and `D_n^{(1)}` are supported)
- ``B`` -- a list of positive integer tuples `[r,s]` specifying the width `s` and height `r` of the sequence of rectangles

EXAMPLES::

    sage: RC = RiggedConfigurations(['A', 3, 1], [[3, 2], [1, 2], [1, 1]])
    sage: RC
    Rigged configurations of type ['A', 3, 1] and factors ((3, 2), (1, 2), (1, 1))

    sage: RC = RiggedConfigurations(['A', 3, 1], [[2,1]]); RC
    Rigged configurations of type ['A', 3, 1] and factors ((2, 1),)
    sage: RC.cardinality()
    6
    sage: len(RC.list()) == RC.cardinality()
    True
    sage: sorted(RC.list())    # random
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
    Rigged configurations of type ['A', 3, 1] and factors ((3, 2), (2, 1), (1, 1), (1, 1))
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
    ...            rigging_list=[[2],[1,0,0],[4,1,2,1],[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,0],[0,0]])
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

    sage: output = elt.to_tensor_product_of_Kirillov_Reshetikhin_tableaux(); output
    [[1, 1, 1], [2, 3, 3], [3, 4, -5]] (X) [[1, 1], [2, 2], [3, 3], [5, -6], [6, -5]] (X) [[1, 1, 2], [2, 2, 3], [3, 3, 7], [4, 4, -7]] (X) [[1, 1, 1], [2, 2, 2]] (X) [[1, 1, 1, 3], [2, 2, 3, 4], [3, 3, 4, 5], [4, 4, 5, 6]] (X) [[1], [2], [3]] (X) [[1, 1, 1, 1]] (X) [[1, 1], [2, 2]]
    sage: elt.to_tensor_product_of_Kirillov_Reshetikhin_tableaux().to_rigged_configuration() == elt
    True
    sage: output.to_rigged_configuration().to_tensor_product_of_Kirillov_Reshetikhin_tableaux() == output
    True

TESTS::

    sage: RC = HighestWeightRiggedConfigurations(['A', 3, 1], [[3,2], [2,1], [1,1], [1,1]])
    sage: RC.cardinality()
    17
    sage: RC = RiggedConfigurations(['D', 4, 1], [[1, 1]])
    sage: RC.cardinality()
    8
    sage: RC = RiggedConfigurations(['D', 4, 1], [[2, 1]])
    sage: c= RC.cardinality(); c
    29
    sage: K = KirillovReshetikhinCrystal(['D',4,1],2,1)
    sage: K.cardinality() == c
    True

.. TODO:: Implement affine crystal structure and remove conversion to classical
    Cartan type.

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
from sage.structure.element import Element, parent
from sage.combinat.misc import IterableFunctionCall
from sage.categories.finite_crystals import FiniteCrystals
from sage.categories.regular_crystals import RegularCrystals
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.partition import Partition
from sage.combinat.cartesian_product import CartesianProduct
from sage.combinat.rigged_configurations.kleber_tree import KleberTree
from sage.combinat.rigged_configurations.rigged_configuration_element import RiggedConfigurationElement

class AbstractRiggedConfigurations(UniqueRepresentation, Parent):
    r"""
    Abstract root class for all rigged configurations.

    This class should never be created directly. Instead see
    :class:`RiggedConfigurations`.
    """
    # TODO: Have the option to have L_{rc} as a dictionary to reduce the memory
    #   footprint in case it is sparse (as a matrix)

    def __init__(self, cartan_type, B, biject_class):
        r"""
        Initialize all common elements for rigged configurations.

        INPUT:

        - ``cartan_type``  -- The Cartan type
        - ``B``            -- A list of dimensions `[r,s]` for rectangles of height `r` and width `s`
        - ``biject_class`` -- The class the bijection creates

        TESTS::

            sage: RC = HighestWeightRiggedConfigurations(['A', 3, 1], [[3, 2], [1, 2], [1, 1]]) # indirect doctest
            sage: RC
            Highest weight rigged configurations of type ['A', 3, 1] and factors ((3, 2), (1, 2), (1, 1))
            sage: RC = RiggedConfigurations(['A', 3, 1], [[3, 2], [1, 2], [1, 1]]) # indirect doctest
            sage: RC
            Rigged configurations of type ['A', 3, 1] and factors ((3, 2), (1, 2), (1, 1))
        """
        assert cartan_type.is_affine()
        self._affine_ct = cartan_type

        # Force to classical since we do not have e_0 and f_0 yet
        self._cartan_type = cartan_type.classical()
        self.dims = B
        self._bijection_class = biject_class
        self.rename("Rigged configurations of type %s and factors %s" % (cartan_type, B))
        Parent.__init__(self, category=(RegularCrystals(), FiniteCrystals()))

    def _highest_weight_iter(self):
        r"""
        Iterate over the highest weight rigged configurations by moving through
        the :class:`KleberTree` and then setting appropriate
        values of the partitions.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A', 2, 1], [[1, 1]]) # indirect doctest
            sage: len([x for x in RC])
            3
        """
        n = self._cartan_type.n

        for tree_node in self.kleber_tree:
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
            base = self(partition_list=shapes[:])

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
                blockLen = partition[0]
                for i, rowLen in enumerate(partition):
                    # If we've gone to a different sized block, then update the
                    #   values which change when moving to a new block size
                    if blockLen != rowLen:
                        blocks[-1].append([])
                        blockLen = rowLen

                    blocks[-1][-1].append(partition.vacancy_numbers[i])

                L2 = []
                for block in blocks[-1]:
                    L2.append(IterableFunctionCall(self._block_iterator, block))
                L.append(CartesianProduct(*L2))

            # TODO Find a more efficient method without appealing to the CartesianProduct
            C = CartesianProduct(*L)
            for curBlocks in C:
                yield self(KT_constructor=[shapes[:], self._blocks_to_values(
                                curBlocks[:]), vac_nums[:]])

    def _block_iterator(self, container):
        r"""
        Iterate over all possible partitions.

        Helper iterator which iterates over all possible partitions contained
        within the container.

        INPUT:

        - ``container`` -- A list of widths of the rows of the container

        TESTS::

            sage: HWRC = HighestWeightRiggedConfigurations(['A', 4, 1], [[2, 2]])
            sage: for x in HWRC._block_iterator([]): x
            ...
            []
            sage: for x in HWRC._block_iterator([2,3]): x
            ...
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
        retPart = [-1] * length
        while pos >= 0:
            retPart[pos] += 1

            if retPart[pos] > container[pos] or (pos != 0 and retPart[pos] > retPart[pos - 1]):
                retPart[pos] = -1
                pos -= 1
            else:
                pos += 1

            if pos == length:
                yield retPart[:]
                pos -= 1

    def _blocks_to_values(self, blocks):
        r"""
        Convert an array of blocks into a list of partition values.

        INPUT:

        - ``blocks`` -- The (2-dim) array blocks of the partition values.

        TESTS::

            sage: HWRC = HighestWeightRiggedConfigurations(['A', 4, 1], [[2, 2]])
            sage: HWRC._blocks_to_values([[[2, 1]]])
            [[2, 1]]
        """
        values = []
        for partBlock in blocks:
            if len(partBlock) == 0:
                values.append([])
            else:
                values.append(partBlock[0][:]) # Need to make a copy
                for block in partBlock[1:]:
                    values[-1].extend(block)
        return values

    def _element_constructor_(self, *lst, **options):
        """
        Construct a RiggedConfigurationElement.

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
                raise ValueError("Incorrect bijection image.")
            return lst[0].to_rigged_configuration()

        from sage.combinat.crystals.tensor_product import TensorProductOfRegularCrystalsElement
        if isinstance(lst[0], TensorProductOfRegularCrystalsElement):
            lst = lst[0]
        from sage.combinat.crystals.kirillov_reshetikhin import KirillovReshetikhinGenericCrystalElement
        if isinstance(lst[0], KirillovReshetikhinGenericCrystalElement):
            KRT = self.tensor_product_of_Kirillov_Reshetikhin_tableaux()
            krt_elt = KRT(*[x.to_Kirillov_Reshetikhin_tableau() for x in lst])
            return krt_elt.to_rigged_configuration()

        return self.element_class(self, *lst, **options)

    def _calc_vacancy_number(self, partitions, a, i, **options):
        r"""
        Calculate the vacancy number of the `i`-th row of the `a`-th rigged
        partition.

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

        #for j, value in enumerate(self._cartan_type.classical().cartan_matrix().row(a)):
        for j, value in enumerate(self._cartan_type.cartan_matrix().row(a)):
            vac_num -= value * partitions[j].get_num_cells_to_column(row_len)

        return vac_num

    def cartan_type(self):
        r"""
        Return the classical Cartan type of this set of rigged configurations.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2, 2]])
            sage: RC.cartan_type()
            ['A', 4]
        """
        return(self._cartan_type)

    @lazy_attribute
    def kleber_tree(self):
        r"""
        Return the underlying Kleber tree used to generate all admissible configurations.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2, 2]])
            sage: RC.kleber_tree
            Kleber tree of Cartan type ['A', 4, 1] and B = ((2, 2),)
        """
        return KleberTree(self._affine_ct, self.dims)

    def tensor_product_of_Kirillov_Reshetikhin_tableaux(self):
        """
        Return the corresponding tensor product of Kirillov-Reshetikhin
        tableaux.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A', 3, 1], [[3, 2], [1, 2]])
            sage: RC.tensor_product_of_Kirillov_Reshetikhin_tableaux()
            Tensor product of Kirillov-Reshetikhin tableaux of type ['A', 3, 1] and tableau shape(s) [[2, 2, 2], [2]]
        """
        return self._bijection_class(self._affine_ct, self.dims)

class HighestWeightRiggedConfigurations(AbstractRiggedConfigurations):
    r"""
    Class for all highest weight rigged configurations.

    A rigged configuration is highest weight if all of its vacancy numbers and
    rigging values are non-negative.

    For more on rigged configurations see :class:`RiggedConfigurations`.
    """

    @staticmethod
    def __classcall_private__(cls, cartan_type, B):
        r"""
        Normalize the input arguments to ensure unique representation.

        EXAMPLES::

            sage: RC1 = HighestWeightRiggedConfigurations(CartanType(['A',3,1]), [[2,2]])
            sage: RC2 = HighestWeightRiggedConfigurations(['A',3,1], [(2,2)])
            sage: RC3 = HighestWeightRiggedConfigurations(['A',3,1], ((2,2),))
            sage: RC2 is RC1, RC3 is RC1
            (True, True)
        """
        cartan_type = CartanType(cartan_type)

        # Standardize B input into a tuple of tuples
        assert B is not None
        B = tuple(tuple(factor) for factor in B)
        return super(HighestWeightRiggedConfigurations, cls).__classcall__(cls, cartan_type, B)

    def __init__(self, cartan_type, B):
        r"""
        Initialize the HighestWeightRiggedConfigurations class.

        INPUT:

        - ``cartan_type`` -- The (affine) Cartan type
        - ``B``           -- A list of dimensions `[r,s]` for rectangles of height `r` and width `s`

        EXAMPLES::

            sage: RC = HighestWeightRiggedConfigurations(['A', 3, 1], [[3, 2], [1, 2], [1, 1]])
            sage: RC
            Highest weight rigged configurations of type ['A', 3, 1] and factors ((3, 2), (1, 2), (1, 1))
            sage: TestSuite(RC).run()
        """

        from tensor_product_kr_tableaux import HighestWeightTensorProductOfKirillovReshetikhinTableaux # For circular imports
        AbstractRiggedConfigurations.__init__(self, cartan_type, B, HighestWeightTensorProductOfKirillovReshetikhinTableaux)
        self.rename("Highest weight rigged configurations of type %s and factors %s" % (cartan_type, B))

    __iter__ = AbstractRiggedConfigurations._highest_weight_iter

    def list(self):
        r"""
        Create a list of the elements by using the iterator.

        EXAMPLES::

            sage: RC = HighestWeightRiggedConfigurations(['D', 4, 1], [[2, 2]])
            sage: len(RC.list()) == RC.cardinality()
            True
            sage: sorted(RC.list())    # random
            [
            (/)
            <BLANKLINE>
            (/)
            <BLANKLINE>
            (/)
            <BLANKLINE>
            (/)
            ,
            0[ ]0
            <BLANKLINE>
            0[ ]0
            0[ ]0
            <BLANKLINE>
            0[ ]0
            <BLANKLINE>
            0[ ]0
            ,
            0[ ][ ]0
            <BLANKLINE>
            0[ ][ ]0
            0[ ][ ]0
            <BLANKLINE>
            0[ ][ ]0
            <BLANKLINE>
            0[ ][ ]0
            ]
        """
        # This is needed to overwrite the list method from the FiniteCrystals
        #   category which generates the list via f_a applications.
        return [x for x in self]

    @lazy_attribute
    def module_generators(self):
        r"""
        Module generators for this set of rigged configurations.

        The module generators for rigged configurations are the highest weight
        rigged configurations. For more info, see :class:`HighestWeightRiggedConfigurations`.

        EXAMPLES::

            sage: RC = HighestWeightRiggedConfigurations(['D', 4, 1], [[2,1]])
            sage: for x in RC.module_generators: x
            ...
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
        # This was separated out for speed; only needed for __iter__
        # lazy_attribute does not inherit
        return [x for x in self._highest_weight_iter()]


HighestWeightRiggedConfigurations.Element = RiggedConfigurationElement

class RiggedConfigurations(AbstractRiggedConfigurations):
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

    INPUT:

    - ``cartan_type`` -- a Cartan type (currently only types `A_n^{(1)}` and `D_n^{(1)}` are supported)
    - ``B`` -- a list of positive integer tuples `[r,s]` specifying the width `s` and height `r` of the sequence of rectangles

    A rigged configuration element with all riggings equal to the vacancy numbers can be created as follows::

        sage: RC = RiggedConfigurations(['A', 3, 1], [[3,2], [2,1], [1,1], [1,1]]); RC
        Rigged configurations of type ['A', 3, 1] and factors ((3, 2), (2, 1), (1, 1), (1, 1))
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
        ...            rigging_list=[[2],[1,0,0],[4,1,2,1],[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,0],[0,0]])
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

        sage: output = elt.to_tensor_product_of_Kirillov_Reshetikhin_tableaux(); output
        [[1, 1, 1], [2, 3, 3], [3, 4, -5]] (X) [[1, 1], [2, 2], [3, 3], [5, -6], [6, -5]] (X)
        [[1, 1, 2], [2, 2, 3], [3, 3, 7], [4, 4, -7]] (X) [[1, 1, 1], [2, 2, 2]] (X)
        [[1, 1, 1, 3], [2, 2, 3, 4], [3, 3, 4, 5], [4, 4, 5, 6]] (X) [[1], [2], [3]] (X) [[1, 1, 1, 1]] (X) [[1, 1], [2, 2]]
        sage: elt.to_tensor_product_of_Kirillov_Reshetikhin_tableaux().to_rigged_configuration() == elt
        True
        sage: output.to_rigged_configuration().to_tensor_product_of_Kirillov_Reshetikhin_tableaux() == output
        True

    We can also convert between rigged configurations and tensor products of
    Kirillov-Reshetikhin crystals::

        sage: RC = RiggedConfigurations(['D', 4, 1], [[2, 1]])
        sage: elt = RC(partition_list=[[1],[1,1],[1],[1]])
        sage: tp_krc = elt.to_tensor_product_of_Kirillov_Reshetikhin_crystals(); tp_krc
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
        sage: ret.to_tensor_product_of_Kirillov_Reshetikhin_crystals()
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

        # Standardize B input into a tuple of tuples
        assert B is not None
        B = tuple(tuple(factor) for factor in B)
        return super(RiggedConfigurations, cls).__classcall__(cls, cartan_type, B)

    def __init__(self, cartan_type, B):
        r"""
        Initialize the RiggedConfigurations class.

        INPUT:

        - ``cartan_type``    -- The Cartan type (the crystal type and n value)
        - ``B``              -- A list of dimensions.

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
            Rigged configurations of type ['A', 3, 1] and factors ((3, 2), (1, 2), (1, 1))
            sage: TestSuite(RC).run()  # long time (4s on sage.math, 2012)
            sage: RC = RiggedConfigurations(['D', 4, 1], [[2,2]])
            sage: TestSuite(RC).run() # long time
            sage: RC = RiggedConfigurations(['D', 4, 1], [[3,1]])
            sage: TestSuite(RC).run() # long time
            sage: RC = RiggedConfigurations(['D', 4, 1], [[4,3]])
            sage: TestSuite(RC).run() # long time
        """
        from tensor_product_kr_tableaux import TensorProductOfKirillovReshetikhinTableaux # For circular imports
        AbstractRiggedConfigurations.__init__(self, cartan_type, B, TensorProductOfKirillovReshetikhinTableaux)
        self.rename("Rigged configurations of type %s and factors %s" % (cartan_type, B))

    @lazy_attribute
    def module_generators(self):
        r"""
        Module generators for this set of rigged configurations.

        The module generators for rigged configurations are the highest weight
        rigged configurations. For more info, see :class:`HighestWeightRiggedConfigurations`.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['D', 4, 1], [[2,1]])
            sage: for x in RC.module_generators: x
            ...
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
        # This was separated out for speed; only needed for __iter__
        # lazy_attribute does not inherit
        return [x for x in self._highest_weight_iter()]

    def _test_bijection(self, **options):
        r"""
        Test function to make sure that the bijection between rigged
        configurations and Kirillov-Reshetikhin tableaux is correct.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[3,2],[4,1]])
            sage: RC._test_bijection()
        """
        tester = self._tester(**options)
        rejects = []
        for x in self:
            y = x.to_tensor_product_of_Kirillov_Reshetikhin_tableaux()
            z = y.to_rigged_configuration()
            if z != x:
                rejects.append((x, z))

        tester.assertTrue(len(rejects) == 0, "Bijection is not correct: %s"%rejects)
        if len(rejects) != 0:
            return rejects

    def tensor_product_of_Kirillov_Reshetikhin_crystals(self):
        """
        Return the corresponding tensor product of Kirillov-Reshetikhin
        crystals.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A', 3, 1], [[3,1],[2,2]])
            sage: RC.tensor_product_of_Kirillov_Reshetikhin_crystals()
            Full tensor product of the crystals [Kirillov-Reshetikhin crystal of type ['A', 3, 1] with (r,s)=(3,1),
            Kirillov-Reshetikhin crystal of type ['A', 3, 1] with (r,s)=(2,2)]
        """
        return self._bijection_class(self._affine_ct, self.dims).tensor_product_of_Kirillov_Reshetikhin_crystals()

RiggedConfigurations.Element = RiggedConfigurationElement

# For experimentation purposes only.
# I'm keeping this for when we implement the R-matrix in general
#def R_matrix_test(n, L1, L2):
#    RC = HighestWeightRiggedConfigurations(['D', n, 1], [L1, L2])
#    RC2 = HighestWeightRiggedConfigurations(['D', n, 1], [L2, L1])
#    ret_list = []
#    for x in RC:
#        x2 = RC2(*x)
#        ret_list.append([x.to_tensor_product_of_Kirillov_Reshetikhin_tableaux(),
#                         x2.to_tensor_product_of_Kirillov_Reshetikhin_tableaux()])
#    return ret_list
