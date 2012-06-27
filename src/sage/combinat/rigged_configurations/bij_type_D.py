r"""
Specific implementations of the bijection classes for type `D_n^{(1)}`.

Part of the (internal) classes which runs the bijection between rigged
configurations and KR tableaux of type `D_n^{(1)}`.

AUTHORS:

- Travis Scrimshaw (2011-04-15): Initial version

TESTS::

    sage: KRT = TensorProductOfKirillovReshetikhinTableaux(['D', 4, 1], [[2,1]])
    sage: from sage.combinat.rigged_configurations.bij_type_D import KRTToRCBijectionTypeD
    sage: bijection = KRTToRCBijectionTypeD(KRT(pathlist=[[-2, 3]]))
    sage: TestSuite(bijection).run()
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

from sage.combinat.rigged_configurations.bij_type_A import KRTToRCBijectionTypeA
from sage.combinat.rigged_configurations.bij_type_A import RCToKRTBijectionTypeA

class KRTToRCBijectionTypeD(KRTToRCBijectionTypeA):
    r"""
    Specific implementation of the bijection from KR tableaux to rigged configurations for type `D_n^{(1)}`.

    This inherits from type `A_n^{(1)}` because we use the same methods in
    some places.
    """

    def next_state(self, val, tableau_height):
        r"""
        Build the next state for type `D_n^{(1)}`.

        EXAMPLES::

            sage: KRT = TensorProductOfKirillovReshetikhinTableaux(['D', 4, 1], [[2,1]])
            sage: from sage.combinat.rigged_configurations.bij_type_D import KRTToRCBijectionTypeD
            sage: bijection = KRTToRCBijectionTypeD(KRT(pathlist=[[-2, 3]]))
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

        # Note that type D_n only contains the (absolute) values between 1 and n
        #   (unlike type A_n which is 1 to n+1). Thus we only need to subtract 1
        #   to match the indices.
        # Also note that we must subtract 1 from n to match the indices as well.

        # If it is a regular value, we follow the A_n rules
        if val > 0:
            KRTToRCBijectionTypeA.next_state(self, val, tableau_height)

            # If we are inserting n-1 (which will only update nu[a] for a <=
            #   n-1), we need to update the vacancy number of nu[n] as well
            if val == self.ret_rig_con.parent()._cartan_type.n - 1:
                self._update_vacancy_nums(val)
            return

        n = self.ret_rig_con.parent()._cartan_type.n
        posVal = -val

        if posVal == n:
            if len(self.ret_rig_con[n - 1]) > 0:
                maxWidth = self.ret_rig_con[n - 1][0] + 1
            else:
                maxWidth = 1

            # Update the last one and skip a rigged partition
            maxWidth = self.ret_rig_con[n - 1].insert_cell(maxWidth)

            for a in reversed(range(tableau_height, n - 2)):
                maxWidth = self.ret_rig_con[a].insert_cell(maxWidth)
                self._update_vacancy_nums(a + 1)
                self._update_partition_values(a + 1)

            # Update all other effected remaining values

            self._update_vacancy_nums(n - 1)
            self._update_partition_values(n - 1)

            self._update_vacancy_nums(tableau_height)
            self._update_partition_values(tableau_height)

            self._update_vacancy_nums(tableau_height - 1)
            self._update_partition_values(tableau_height - 1)

            self._update_vacancy_nums(n - 2)
            return

        # Always add a cell to the first singular value in the first
        #   tableau we are updating.
        if len(self.ret_rig_con[posVal - 1]) > 0:
            maxWidth = self.ret_rig_con[posVal - 1][0] + 1
        else:
            maxWidth = 1
        # Special case for `\overline{n-1}` to take the larger of the last two
        if posVal == n - 1 and len(self.ret_rig_con[n - 1]) > 0 and \
          self.ret_rig_con[n - 1][0] + 1 > maxWidth:
            maxWidth = self.ret_rig_con[n - 1][0] + 1


        # Add cells similiar to type A_n but we move to the right until we reach
        #   the value of n-2
        for a in range(posVal - 1, n - 2):
            maxWidth = self.ret_rig_con[a].insert_cell(maxWidth)

        # Handle the special behavior near values of n
        maxWidth2 = self.ret_rig_con[n - 2].insert_cell(maxWidth)
        maxWidth = self.ret_rig_con[n - 1].insert_cell(maxWidth)

        if maxWidth2 < maxWidth:
            maxWidth = maxWidth2

        # Go back following the regular A_n rules
        maxWidth = self.ret_rig_con[n - 3].insert_cell(maxWidth)

        self._update_vacancy_nums(n - 2)
        self._update_partition_values(n - 2)
        self._update_vacancy_nums(n - 1)
        self._update_partition_values(n - 1)

        for a in reversed(range(tableau_height, n - 3)):
            maxWidth = self.ret_rig_con[a].insert_cell(maxWidth)
            self._update_vacancy_nums(a + 1)
            self._update_partition_values(a + 1)

        # Update the final rigged tableau
        self._update_vacancy_nums(tableau_height)
        self._update_partition_values(tableau_height)

        if posVal < tableau_height:
            for a in range(posVal - 1, tableau_height):
                self._update_vacancy_nums(a)
                self._update_partition_values(a)
        elif tableau_height > 0:
            self._update_vacancy_nums(tableau_height - 1)
            self._update_partition_values(tableau_height - 1)

class RCToKRTBijectionTypeD(RCToKRTBijectionTypeA):
    r"""
    Specific implementation of the bijection from rigged configurations to tensor products of KR tableaux for type `D_n^{(1)}`.
    """

    def next_state(self):
        r"""
        Build the next state for type `D_n^{(1)}`.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['D', 4, 1], [[2, 1]])
            sage: from sage.combinat.rigged_configurations.bij_type_D import RCToKRTBijectionTypeD
            sage: bijection = RCToKRTBijectionTypeD(RC(partition_list=[[],[1,1],[1],[1]]))
            sage: bijection.tj(1)
            sage: bijection.next_state()
            1
            sage: bijection.cur_partitions
            [(/)
            , -1[ ]-1
            -1[ ]-1
            , 0[ ]0
            , 0[ ]0
            ]

        """
        n = self.rigged_con.parent()._cartan_type.n
        ell = [None] * (2 * n - 2) # No `\bar{\ell}^{n-1}` and `\bar{\ell}^n`
        b = None

        # Calculate the rank and ell values

        last_size = 0
        for a in range(0, n - 2):
            ell[a] = self._find_singular_string(self.cur_partitions[a], last_size)

            if ell[a] is None:
                b = a + 1
                break
            else:
                last_size = self.cur_partitions[a][ell[a]]

        if b is None:
            # Do the special cases when we've reached n - 2
            ell[n - 2] = self._find_singular_string(self.cur_partitions[n - 2], last_size)
            ell[n - 1] = self._find_singular_string(self.cur_partitions[n - 1], last_size)

            if ell[n - 2] is not None:
                temp_size = self.cur_partitions[n - 2][ell[n - 2]]
                if ell[n - 1] is not None:
                    last_size = self.cur_partitions[n - 1][ell[n - 1]]
                    if temp_size > last_size:
                        last_size = temp_size
                else:
                    b = n
            else:
                if ell[n - 1] is not None:
                    b = -n
                else:
                    b = n - 1

        if b is None:
            # Now go back
            for a in reversed(range(0, n - 2)):
                # Modified form of _find_singular_string
                for i in reversed(range(0, ell[a])):
                    if self.cur_partitions[a][i] >= last_size and \
                      self.cur_partitions[a].vacancy_numbers[i] == self.cur_partitions[a].rigging[i]:
                            ell[n + a] = i
                            break

                if ell[n + a] is None:
                    b = -(a + 2)
                    break
                else:
                    last_size = self.cur_partitions[a][ell[n + a]]

        if b is None:
            b = -1

        # Determine the new rigged configuration by removing a box from the selected
        #   string and then making the new string singular
        ret_row = self.cur_partitions[0].remove_cell(ell[0])
        ret_row_bar = self.cur_partitions[0].remove_cell(ell[n])
        for a in range(1, n - 2):
            ret_row_next = self.cur_partitions[a].remove_cell(ell[a])
            ret_row_bar_next = self.cur_partitions[a].remove_cell(ell[n + a])

            self._update_vacancy_numbers(a - 1)
            if ret_row is not None:
                self.cur_partitions[a - 1].rigging[ret_row] = \
                  self.cur_partitions[a - 1].vacancy_numbers[ret_row]
            if ret_row_bar is not None:
                self.cur_partitions[a - 1].rigging[ret_row_bar] = \
                  self.cur_partitions[a - 1].vacancy_numbers[ret_row_bar]

            ret_row = ret_row_next
            ret_row_bar = ret_row_bar_next

        # Special behavior for all a > n-3
        ret_row_next = self.cur_partitions[n - 2].remove_cell(ell[n - 2])
        ret_row_bar_next = self.cur_partitions[n - 1].remove_cell(ell[n - 1])

        self._update_vacancy_numbers(n - 3)
        if ret_row is not None:
            self.cur_partitions[n - 3].rigging[ret_row] = \
              self.cur_partitions[n - 3].vacancy_numbers[ret_row]
        if ret_row_bar is not None:
            self.cur_partitions[n - 3].rigging[ret_row_bar] = \
              self.cur_partitions[n - 3].vacancy_numbers[ret_row_bar]

        self._update_vacancy_numbers(n - 2)
        if ret_row_next is not None:
            self.cur_partitions[n - 2].rigging[ret_row_next] = \
              self.cur_partitions[n - 2].vacancy_numbers[ret_row_next]

        self._update_vacancy_numbers(n - 1)
        if ret_row_bar_next is not None:
            self.cur_partitions[n - 1].rigging[ret_row_bar_next] = \
              self.cur_partitions[n - 1].vacancy_numbers[ret_row_bar_next]

        return(b)

    def tj(self, k):
        r"""
        The map `tj` for crystals of type `D_n^{(1)}`.

        INPUT:

        - ``k`` -- The height of the crystal `B^{k,1}`.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['D', 4, 1], [[2, 1]])
            sage: from sage.combinat.rigged_configurations.bij_type_D import RCToKRTBijectionTypeD
            sage: bijection = RCToKRTBijectionTypeD(RC(partition_list=[[1],[2,1],[2],[1]]))
            sage: bijection.tj(2)
            sage: bijection.cur_partitions
            [0[ ]0
            0[ ]0
            , -1[ ][ ]-1
            0[ ]0
            , -1[ ][ ]-1
            , 0[ ]0
            ]
        """
        n = len(self.cur_partitions)

        if k == n - 1:
            # TODO
            pass
        elif k == n:
            # TODO
            pass
        else:
            RCToKRTBijectionTypeA.tj(self, k)

