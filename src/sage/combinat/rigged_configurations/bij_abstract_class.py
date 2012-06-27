r"""
Abstract classes for the rigged configuration bijections

This file contains two sets of classes, one for the bijection from KR tableaux to
rigged configurations and the other for the reverse bijection. We do this for
two reasons, one is because we can store a state in the bijection locally, so
we do not have to constantly pass it around between functions. The other is because
it makes the code easier to read in the \*_element.py files.

These classes are not meant to be used by the user and are only supposed to be
used internally to perform the bijections between
:class:`TensorProductOfKirillovReshetikhinTableaux` and
:class:`RiggedConfigurations`.

AUTHORS:

- Travis Scrimshaw (2011-04-15): Initial version
"""

#*****************************************************************************
#       Copyright (C) 2011, 2012 Travis Scrimshaw <tscrim@ucdavis.edu>
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

from copy import deepcopy
from sage.misc.abstract_method import abstract_method

class KRTToRCBijectionAbstract:
    """
    Root abstract class for the bijection from KR tableaux to rigged configurations.

    This class holds the state of the bijection and generates the next state.
    This class should never be created directly.
    """

    def __init__(self, krt):
        """
        Initialize the bijection by obtaining the important information from
        the KR tableaux.

        INPUT:

        - ``krt`` -- The tensor product of KR tableaux

        EXAMPLES::

            sage: KRT = TensorProductOfKirillovReshetikhinTableaux(['A', 4, 1], [[2,1]])
            sage: from sage.combinat.rigged_configurations.bij_abstract_class import KRTToRCBijectionAbstract
            sage: bijection = KRTToRCBijectionAbstract(KRT(pathlist=[[4,3]]))
            sage: bijection.cur_path
            []
            sage: bijection.ret_rig_con
            <BLANKLINE>
            (/)
            <BLANKLINE>
            (/)
            <BLANKLINE>
            (/)
            <BLANKLINE>
            (/)
            <BLANKLINE>
            sage: TestSuite(bijection).run()
        """
        # L = [[0]] * crystalPath.parent().cartan_type().n
        # for dim in crystalPath.parent().dims:
        #     L[dim[0]-1][0] += 1

        self.ret_rig_con = krt.parent()._bijection_class(
          krt.parent().affine_ct,
          krt.parent().dims)(partition_list=[[]] *
          krt.parent().cartan_type().n)
        # We allow this to be mutable to make the bijection easier to program.
        # Upon completing the bijection, this will be set to immutable.
        # Do not call this, the object could be in a mutable state and ultimately
        #   be placed in an unstable state.
        # The user will (and should) never know about this temporary mutable state.
        self.ret_rig_con._set_mutable()
        self.cur_path = []

    def __eq__(self, rhs):
        r"""
        Check equality.

        This is only here for pickling check. This is a temporary placeholder
        class, and as such, should never be compared.

        TESTS::

            sage: KRT = TensorProductOfKirillovReshetikhinTableaux(['A', 4, 1], [[2,1]])
            sage: from sage.combinat.rigged_configurations.bij_type_A import KRTToRCBijectionTypeA
            sage: bijection = KRTToRCBijectionTypeA(KRT(pathlist=[[4,3]]))
            sage: bijection2 = KRTToRCBijectionTypeA(KRT(pathlist=[[4,3]]))
            sage: bijection == bijection2
            True
        """
        return isinstance(rhs, KRTToRCBijectionAbstract)

    @abstract_method
    def next_state(self, val, tableau_height):
        r"""
        Build the next state in the bijection.

        INPUT:

        - ``val``            -- The value we are adding
        - ``tableau_height`` -- The height of the tableau

        TESTS::

            sage: KRT = TensorProductOfKirillovReshetikhinTableaux(['A', 4, 1], [[2,1]])
            sage: from sage.combinat.rigged_configurations.bij_type_A import KRTToRCBijectionTypeA
            sage: bijection = KRTToRCBijectionTypeA(KRT(pathlist=[[4,3]]))
            sage: bijection.next_state(3, 0)
            sage: bijection.ret_rig_con
            <BLANKLINE>
            -1[ ]-1
            <BLANKLINE>
            -1[ ]-1
            <BLANKLINE>
            (/)
            <BLANKLINE>
            (/)
            <BLANKLINE>
        """

    def _update_vacancy_nums(self, a):
        r"""
        Update the vacancy numbers of a rigged partition.

        Helper function to (batch) update the vacancy numbers of the rigged
        partition at position `a` in the rigged configuration stored by this
        bijection.

        INPUT:

        - ``a`` -- The index of the partition to update

        TESTS::

            sage: KRT = TensorProductOfKirillovReshetikhinTableaux(['A', 4, 1], [[2,1]])
            sage: from sage.combinat.rigged_configurations.bij_abstract_class import KRTToRCBijectionAbstract
            sage: bijection = KRTToRCBijectionAbstract(KRT(pathlist=[[4,3]]))
            sage: bijection._update_vacancy_nums(2)
        """
        # Check to make sure we have a valid index (currently removed)
        # If the current tableau is empty, there is nothing to do
        if len(self.ret_rig_con[a]) == 0: # Check to see if we have vacancy numbers
            return

        # Setup the first block
        blockLen = self.ret_rig_con[a][0]
        vac_num = self.ret_rig_con.parent()._calc_vacancy_number(self.ret_rig_con.nu(),
                                                                 a, 0, B=self.cur_path)

        for i, row_len in enumerate(self.ret_rig_con[a]):
            # If we've gone to a different sized block, then update the
            #   values which change when moving to a new block size
            if blockLen != row_len:
                vac_num = self.ret_rig_con.parent()._calc_vacancy_number(self.ret_rig_con.nu(),
                                                                         a, i, B=self.cur_path)
                blockLen = row_len
            self.ret_rig_con[a].vacancy_numbers[i] = vac_num

    def _update_partition_values(self, a):
        r"""
        Update the partition values of a rigged partition.

        Helper function to update the partition values of a given rigged
        partition row. This will go through all of our partition values and set
        them to our vacancy number if the corresponding row has been changed
        (indicated by being set to None).

        INPUT:

        - ``a`` -- The index of the partition to update

        TESTS::

            sage: KRT = TensorProductOfKirillovReshetikhinTableaux(['A', 4, 1], [[2,1]])
            sage: from sage.combinat.rigged_configurations.bij_abstract_class import KRTToRCBijectionAbstract
            sage: bijection = KRTToRCBijectionAbstract(KRT(pathlist=[[4,3]]))
            sage: bijection._update_partition_values(2)
        """
        rigged_partition = self.ret_rig_con[a]
        for index, value in enumerate(rigged_partition.rigging):
            if value == None:
                rigged_partition.rigging[index] = rigged_partition.vacancy_numbers[index]
                if index > 0  and rigged_partition[index - 1] == rigged_partition[index] \
                  and rigged_partition.rigging[index - 1] < rigged_partition.rigging[index]:
                    # If we need to reorder
                    pos = 0
                    width = rigged_partition[index]
                    #width = rigged_partition._get_item(index) # ...New way... - Change all
                    val = rigged_partition.rigging[index]
                    for i in reversed(range(index - 1)):
                        if rigged_partition[i] > width or rigged_partition.rigging[i] >= val:
                            pos = i
                            break

                    rigged_partition.rigging.pop(index)
                    rigged_partition.rigging.insert(pos, val)

class RCToKRTBijectionAbstract:
    """
    Root abstract class for the bijection from rigged configurations to crystal paths.

    This class holds the state of the bijection and generates the next state.
    This class should never be created directly.
    """

    def __init__(self, RC_element):
        """
        Initialize the bijection helper.

        INPUT:

        - ``RC_element`` -- The rigged configuration

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2, 1]])
            sage: from sage.combinat.rigged_configurations.bij_abstract_class import RCToKRTBijectionAbstract
            sage: bijection = RCToKRTBijectionAbstract(RC(partition_list=[[1],[1],[1],[1]]))
            sage: bijection.rigged_con
            <BLANKLINE>
            -1[ ]-1
            <BLANKLINE>
            1[ ]1
            <BLANKLINE>
            0[ ]0
            <BLANKLINE>
            -1[ ]-1
            <BLANKLINE>
            sage: TestSuite(bijection).run()
        """
        # L = [[0]] * crystalPath.parent().cartan_type().n
        # for dim in crystalPath.parent().dims:
        #     L[dim[0]-1][0] += 1

        # Make a mutable clone of the rigged configuration for the bijection
        # This will be deleted when the bijection is completed
        #self.rigged_con = deepcopy(RC_element)
        self.rigged_con = RC_element.__copy__()

        # Build an empty path for the vacancy numbers
        self.rem_path = []
        for dim in self.rigged_con.parent().dims:
            self.rem_path.append([])
            for i in range(dim[0]):
                self.rem_path[-1].append([None] * dim[1])

        # Note that this implementation of the bijection is destructive to cur_partitions,
        #   therefore we will make a (deep) copy of the partitions.
        # TODO: Convert from cur_partitions to rigged_con
        self.cur_partitions = deepcopy(list(self.rigged_con)[:])

    def __eq__(self, rhs):
        r"""
        Check equality.

        This is only here for pickling check. This is a temporary placeholder
        class, and as such, should never be compared.

        TESTS::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2, 1]])
            sage: from sage.combinat.rigged_configurations.bij_type_A import RCToKRTBijectionTypeA
            sage: bijection = RCToKRTBijectionTypeA(RC(partition_list=[[1],[1],[1],[1]]))
            sage: bijection2 = RCToKRTBijectionTypeA(RC(partition_list=[[1],[1],[1],[1]]))
            sage: bijection == bijection2
            True
        """
        return isinstance(rhs, RCToKRTBijectionAbstract)

    @abstract_method
    def next_state(self):
        """
        Build the next state in the bijection.

        TESTS::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2, 1]])
            sage: from sage.combinat.rigged_configurations.bij_type_A import RCToKRTBijectionTypeA
            sage: bijection = RCToKRTBijectionTypeA(RC(partition_list=[[1],[1],[1],[1]]))
            sage: bijection.tj(1)
            sage: bijection.next_state()
            5
            sage: bijection.cur_partitions
            [(/)
            , (/)
            , (/)
            , (/)
            ]
        """

    def _update_vacancy_numbers(self, a):
        r"""
        Update the vacancy numbers during the bijection.

        INPUT:

        - ``a`` -- The index of the partition to update

        TESTS::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2, 1]])
            sage: from sage.combinat.rigged_configurations.bij_abstract_class import RCToKRTBijectionAbstract
            sage: bijection = RCToKRTBijectionAbstract(RC(partition_list=[[1],[1],[1],[1]]))
            sage: bijection._update_vacancy_numbers(2)
        """

        # Nothing to do if there the rigged partition is empty
        if len(self.cur_partitions[a]) == 0:
            return

        partition = self.cur_partitions[a]

        # Setup the first block
        blockLen = partition[0]
        vacNum = self.rigged_con.parent()._calc_vacancy_number(self.cur_partitions,
                                                               a, 0, B=self.rem_path)

        for i, rowLen in enumerate(self.cur_partitions[a]):
            # If we've gone to a different sized block, then update the
            #   values which change when moving to a new block size
            if blockLen != rowLen:
                vacNum = self.rigged_con.parent()._calc_vacancy_number(self.cur_partitions,
                                                                       a, i, B=self.rem_path)
                blockLen = rowLen

            partition.vacancy_numbers[i] = vacNum

    @abstract_method
    def tj(self, k):
        r"""
        Perform the map `tj` to "move" our path to one with a leading box of
        height 1.

        TESTS::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2, 1]])
            sage: from sage.combinat.rigged_configurations.bij_type_A import RCToKRTBijectionTypeA
            sage: bijection = RCToKRTBijectionTypeA(RC(partition_list=[[1],[1],[1],[1]]))
            sage: bijection.tj(1)
            sage: bijection.cur_partitions
            [-1[ ]-1
            , 1[ ]1
            , 0[ ]0
            , -1[ ]-1
            ]
        """

    def _find_singular_string(self, partition, last_size):
        r"""
        Return the index of the singular string or None if not found.

        Helper method to find a singular string at least as long as last_size.

        INPUT:

        - ``partition`` -- The partition to look in
        - ``last_size`` -- The last size found

        TESTS::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2, 1]])
            sage: from sage.combinat.rigged_configurations.bij_abstract_class import RCToKRTBijectionAbstract
            sage: bijection = RCToKRTBijectionAbstract(RC(partition_list=[[1],[1],[1],[1]]))
            sage: bijection._find_singular_string(bijection.cur_partitions[2], 2)
            sage: bijection._find_singular_string(bijection.cur_partitions[2], 0)
            0
        """
        for i in reversed(range(0, len(partition))):
            if partition[i] >= last_size and \
              partition.vacancy_numbers[i] == partition.rigging[i]:
                return i

