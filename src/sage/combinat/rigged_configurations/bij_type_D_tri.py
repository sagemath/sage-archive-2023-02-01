r"""
Bijection classes for type `D_4^{(3)}`.

Part of the (internal) classes which runs the bijection between rigged
configurations and KR tableaux of type `D_4^{(3)}`.

AUTHORS:

- Travis Scrimshaw (2014-09-10): Initial version

TESTS::

    sage: KRT = crystals.TensorProductOfKirillovReshetikhinTableaux(['D', 4, 3], [[2,1]])
    sage: from sage.combinat.rigged_configurations.bij_type_D_tri import KRTToRCBijectionTypeDTri
    sage: bijection = KRTToRCBijectionTypeDTri(KRT(pathlist=[[-1,2]]))
    sage: TestSuite(bijection).run()
    sage: RC = RiggedConfigurations(['D', 4, 3], [[2, 1]])
    sage: from sage.combinat.rigged_configurations.bij_type_D_tri import RCToKRTBijectionTypeDTri
    sage: bijection = RCToKRTBijectionTypeDTri(RC(partition_list=[[],[],[]]))
    sage: TestSuite(bijection).run()
"""

#*****************************************************************************
#       Copyright (C) 2014 Travis Scrimshaw <tscrim@ucdavis.edu>
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

from sage.combinat.rigged_configurations.bij_type_A import KRTToRCBijectionTypeA
from sage.combinat.rigged_configurations.bij_type_A import RCToKRTBijectionTypeA

class KRTToRCBijectionTypeDTri(KRTToRCBijectionTypeA):
    r"""
    Specific implementation of the bijection from KR tableaux to rigged
    configurations for type `D_4^{(3)}`.

    This inherits from type `A_n^{(1)}` because we use the same methods in
    some places.
    """

    def next_state(self, val):
        r"""
        Build the next state for type `D_4^{(3)}`.

        TESTS::

            sage: KRT = crystals.TensorProductOfKirillovReshetikhinTableaux(['D', 4, 3], [[2,1]])
            sage: from sage.combinat.rigged_configurations.bij_type_D_tri import KRTToRCBijectionTypeDTri
            sage: bijection = KRTToRCBijectionTypeDTri(KRT(pathlist=[[-1,2]]))
            sage: bijection.cur_path.insert(0, [])
            sage: bijection.cur_dims.insert(0, [0, 1])
            sage: bijection.cur_path[0].insert(0, [2])
            sage: bijection.next_state(2)
        """

        n = self.n
        tableau_height = len(self.cur_path[0]) - 1

        if val == 'E':
            self.ret_rig_con[0].insert_cell(max_width)
            self.ret_rig_con[1].insert_cell(max_width)
            if tableau_height == 0:
                self.ret_rig_con[0].insert_cell(max_width)
            self._update_vacancy_nums(0)
            self._update_vacancy_nums(1)
            self._update_partition_values(0)
            self._update_partition_values(1)
            return

        if val > 0:
            # If it is a regular value, we follow the A_n rules
            KRTToRCBijectionTypeA.next_state(self, val)
            return

        pos_val = -val

        if pos_val == 0:
            if len(self.ret_rig_con[pos_val - 1]) > 0:
                max_width = self.ret_rig_con[0][0]
            else:
                max_width = 1
            max_width = self.ret_rig_con[0].insert_cell(max_width)
            width_n = max_width + 1

            # Follow regular A_n rules
            for a in reversed(range(tableau_height, 1)):
                max_width = self.ret_rig_con[a].insert_cell(max_width)
                self._update_vacancy_nums(a + 1)
                self._update_partition_values(a + 1)
            self._update_vacancy_nums(tableau_height)
            self._update_partition_values(tableau_height)
            if tableau_height > 0:
                self._update_vacancy_nums(tableau_height-1)
                self._update_partition_values(tableau_height-1)

            # Make the largest string at \nu^{(1)} quasi-singular
            p = self.ret_rig_con[0]
            num_rows = len(p)
            for i in range(num_rows):
                if p._list[i] == width_n:
                    j = i+1
                    while j < num_rows and p._list[j] == width_n \
                      and p.vacancy_numbers[j] == p.rigging[j]:
                        j += 1
                    p.rigging[j-1] -= 1
                    break
            return

        case_S = [None] * n
        pos_val = -val

        if pos_val > 3:
            # Always add a cell to the first singular value in the first
            #   tableau we are updating.
            if len(self.ret_rig_con[pos_val - 1]) > 0:
                max_width = self.ret_rig_con[pos_val - 1][0]
            else:
                max_width = 1

            # Add cells similar to type A_n but we move to the right
            for a in range(pos_val - 1, 2):
                max_width = self.ret_rig_con[a].insert_cell(max_width)
                case_S[a] = max_width
        else:
            if len(self.ret_rig_con[0]) > 0:
                max_width = self.ret_rig_con[0][0]
            else:
                max_width = 1

        # Special case for going through 0
        # If we find a quasi-singular string first, then we are in case (Q, S)
        #   otherwise we will find a singular string and insert 2 cells
        P = self.ret_rig_con[0]
        num_rows = len(P)
        case_QS = False
        for i in range(num_rows + 1):
            if i == num_rows:
                max_width = 0
                if case_QS:
                    P._list.append(1)
                    P.vacancy_numbers.append(None)
                    # Go through our partition until we find a length of greater than 1
                    j = len(P._list) - 1
                    while j >= 0 and P._list[j] == 1:
                        j -= 1
                    P.rigging.insert(j + 1, None)
                    width_n = 1
                else:
                    # Go through our partition until we find a length of greater than 2
                    j = len(P._list) - 1
                    while j >= 0 and P._list[j] <= 2:
                        j -= 1
                    P._list.insert(j+1, 2)
                    P.vacancy_numbers.insert(j+1, None)
                    P.rigging.insert(j+1, None)
                break
            elif P._list[i] <= max_width:
                if P.vacancy_numbers[i] == P.rigging[i]:
                    max_width = P._list[i]
                    if case_QS:
                        P._list[i] += 1
                        width_n = P._list[i]
                        P.rigging[i] = None
                    else:
                        j = i - 1
                        while j >= 0 and P._list[j] <= max_width + 2:
                            P.rigging[j+1] = P.rigging[j] # Shuffle it along
                            j -= 1
                        P._list.pop(i)
                        P._list.insert(j+1, max_width + 2)
                        P.rigging[j+1] = None
                    break
                elif P.vacancy_numbers[i] - 1 == P.rigging[i] and not case_QS:
                    case_QS = True
                    P._list[i] += 1
                    P.rigging[i] = None
                    # No need to set max_width here since we will find a singular string

        # Now go back following the regular C_n (ish) rules
        for a in reversed(range(tableau_height, 2)):
            if case_S[a] == max_width:
                self._insert_cell_case_S(self.ret_rig_con[a])
            else:
                max_width = self.ret_rig_con[a].insert_cell(max_width)
            self._update_vacancy_nums(a + 1)
            self._update_partition_values(a + 1)

        # Update the final rigged partitions
        self._update_vacancy_nums(tableau_height)
        if tableau_height >= n - 1:
            self._correct_vacancy_nums()
        self._update_partition_values(tableau_height)

        if pos_val <= tableau_height:
            for a in range(pos_val-1, tableau_height):
                self._update_vacancy_nums(a)
                self._update_partition_values(a)
            if pos_val > 1:
                self._update_vacancy_nums(pos_val - 2)
                self._update_partition_values(pos_val - 2)
        elif tableau_height > 0:
            self._update_vacancy_nums(tableau_height - 1)
            self._update_partition_values(tableau_height - 1)

        if case_QS:
            # Make the new string quasi-singular
            num_rows = len(partition)
            for i in range(num_rows):
                if partition._list[i] == width_n:
                    j = i+1
                    while j < num_rows and partition._list[j] == width_n \
                      and partition.vacancy_numbers[j] == partition.rigging[j]:
                        j += 1
                    partition.rigging[j-1] -= 1
                    break

    def _insert_cell_case_S(self, partition):
        """
        Insert a cell when case `(S)` holds.

        TESTS::

            sage: RC = RiggedConfigurations(['C', 2, 1], [[2, 2]])
            sage: RP = RC(partition_list=[[2],[2,2]])[1]
            sage: RP
            -4[ ][ ]-4
            -4[ ][ ]-4
            <BLANKLINE>
            sage: RP.rigging[0] = None
            sage: from sage.combinat.rigged_configurations.bij_type_D_tri import KRTToRCBijectionTypeDTri
            sage: KRT = crystals.TensorProductOfKirillovReshetikhinTableaux(['D', 4, 3], [[2,1]])
            sage: bijection = KRTToRCBijectionTypeDTri(KRT(pathlist=[[-1,2]]))
            sage: bijection._insert_cell_case_S(RP)
            sage: RP
            -4[ ][ ][ ]None
            -4[ ][ ]-4
            <BLANKLINE>
        """
        # Special case when adding twice to the first row
        if partition.rigging[0] is None:
            partition._list[0] += 1
            return

        num_rows = len(partition)
        for i in reversed(range(1, num_rows)):
            if partition.rigging[i] is None:
                j = i - 1
                while j >= 0 and partition._list[j] == partition._list[i]:
                    partition.rigging[j+1] = partition.rigging[j] # Shuffle it along
                    j -= 1
                partition._list[j+1] += 1
                partition.rigging[j+1] = None
                return

class RCToKRTBijectionTypeDTri(RCToKRTBijectionTypeA):
    r"""
    Specific implementation of the bijection from rigged configurations to
    tensor products of KR tableaux for type `D_4^{(3)}`.
    """

    def next_state(self, height):
        r"""
        Build the next state for type `D_4^{(3)}`.

        TESTS::

            sage: RC = RiggedConfigurations(['D', 4, 3], [[2, 1]])
            sage: from sage.combinat.rigged_configurations.bij_type_D_tri import RCToKRTBijectionTypeDTri
            sage: bijection = RCToKRTBijectionTypeDTri(RC(partition_list=[[2],[2],[1]]))
            sage: bijection.next_state(0)
            -1
        """
        ell = [None] * 6
        case_S = [False] * 3
        case_Q = False
        b = None

        # Calculate the rank and ell values

        last_size = 0
        for a in range(height, 2):
            ell[a] = self._find_singular_string(self.cur_partitions[a], last_size)

            if ell[a] is None:
                b = a + 1
                break
            else:
                last_size = self.cur_partitions[a][ell[a]]

        if b is None:
            partition = self.cur_partitions[0]
            # Modified version of _find_singular_string()
            for i in reversed(range(len(partition))):
                if partition[i] >= last_size:
                    if partition.vacancy_numbers[i] == partition.rigging[i] and i != ell[0]:
                        if partition[i] == 1:
                            b = 'E'
                        else:
                            last_size = partition[i]
                        case_S[2] = True
                        ell[3] = i
                        break
                    elif partition.vacancy_numbers[i] - 1 == partition.rigging[i] and not case_Q:
                        case_Q = True
                        # Check if the block is singular
                        block_size = partition[i]
                        for j in reversed(range(i)):
                            if partition[j] != block_size:
                                break
                            elif partition.vacancy_numbers[j] == partition.rigging[j] and j != ell[0]:
                                case_Q = False
                                break
                        if case_Q:
                            last_size = partition[i] + 1
                            ell[2] = i

            if ell[3] is None:
                if not case_Q:
                    b = 3
                else:
                    b = 0

        if b is None: # Going back
            if self.cur_partitions[1][ell[1]] == last_size:
                ell[4] = ell[1]
                case_S[1] = True
            else:
                ell[4] = self._find_singular_string(self.cur_partitions[1], last_size)

                if ell[4] is None:
                    b = -3
                else:
                    last_size = self.cur_partitions[1][ell[4]]

        if b is None: # Final partition
            P = self.cur_partitions[0]
            if ell[0] is not None and P[ell[0]] == last_size:
                ell[5] = ell[0]
                case_S[0] = True
            else:
                # Modified form of _find_singular_string
                end = ell[3]
                for i in reversed(range(end)):
                    if P[i] >= last_size and P.vacancy_numbers[i] == P.rigging[i]:
                        ell[5] = i
                        break

                if ell[5] is None:
                    b = -2

        if b is None:
            b = -1

        # Determine the new rigged configuration by removing boxes from the
        #   selected string and then making the new string singular
        if case_S[1]:
            row1 = [self.cur_partitions[1].remove_cell(ell[4], 2)]
        else:
            row1 = [self.cur_partitions[1].remove_cell(ell[1]), 
                    self.cur_partitions[1].remove_cell(ell[4])]

        if case_S[0]:
            row0 = [self.cur_partitions[0].remove_cell(ell[5], 2)]
            row0.append( self.cur_partitions[0].remove_cell(ell[3], 2) )
        else:
            if case_Q:
                if ell[0] < ell[2]:
                    row0 = [self.cur_partitions[0].remove_cell(ell[2]),
                            self.cur_partitions[0].remove_cell(ell[0])]
                else:
                    row0 = [self.cur_partitions[0].remove_cell(ell[0]),
                            self.cur_partitions[0].remove_cell(ell[2])]
                if case_S[2]:
                    quasi = self.cur_partitions[0].remove_cell(ell[3])
            else:
                row0 = [self.cur_partitions[0].remove_cell(ell[0])]
                if case_S[2]:
                    row0.append( self.cur_partitions[0].remove_cell(ell[3], 2) )

            row0.append( self.cur_partitions[0].remove_cell(ell[5]) )

        self._update_vacancy_numbers(0)
        self._update_vacancy_numbers(1)

        for l in row1:
            if l is not None:
                self.cur_partitions[1].rigging[l] = self.cur_partitions[1].vacancy_numbers[l]
        for l in row0:
            if l is not None:
                self.cur_partitions[0].rigging[l] = self.cur_partitions[0].vacancy_numbers[l]

        # If case (Q,S) holds, then we must make the larger string quasisingular    
        if case_Q and case_S[2]:
            P = self.cur_partitions[0]
            vac_num = P.vacancy_numbers[quasi]
            P.rigging[quasi] = vac_num
            block_len = P[quasi]
            j = quasi + 1
            length = len(P)
            # Find the place for the quasisingular rigging
            while j < length and P[j] == block_len and P.rigging[j] == vac_num:
                j += 1
            P.rigging[j-1] = vac_num - 1

        return(b)

