r"""
Bijection classes for type `E_{6,7}^{(1)}`

Part of the (internal) classes which runs the bijection between rigged
configurations and KR tableaux of type `E_{6,7}^{(1)}`.

AUTHORS:

- Travis Scrimshaw (2011-04-15): Initial version

TESTS::

    sage: KRT = crystals.TensorProductOfKirillovReshetikhinTableaux(['D', 4, 1], [[2,1]])
    sage: from sage.combinat.rigged_configurations.bij_type_D import KRTToRCBijectionTypeD
    sage: bijection = KRTToRCBijectionTypeD(KRT(pathlist=[[3, 2]]))
    sage: TestSuite(bijection).run()
    sage: RC = RiggedConfigurations(['D', 4, 1], [[2, 1]])
    sage: from sage.combinat.rigged_configurations.bij_type_D import RCToKRTBijectionTypeD
    sage: bijection = RCToKRTBijectionTypeD(RC(partition_list=[[],[],[],[]]))
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

from sage.combinat.rigged_configurations.bij_abstract_class import KRTToRCBijectionAbstract
from sage.combinat.rigged_configurations.bij_abstract_class import RCToKRTBijectionAbstract

class KRTToRCBijectionTypeE67(KRTToRCBijectionAbstract):
    r"""
    Specific implementation of the bijection from KR tableaux to rigged
    configurations for type `E_{6,7}^{(1)}`.
    """
    def next_state(self, val):
        r"""
        Build the next state for type `E_{6,7}^{(1)}`.

        TESTS::

            sage: KRT = crystals.TensorProductOfKirillovReshetikhinTableaux(['D', 4, 1], [[2,1]])
            sage: from sage.combinat.rigged_configurations.bij_type_D import KRTToRCBijectionTypeD
            sage: bijection = KRTToRCBijectionTypeD(KRT(pathlist=[[5,3]]))
            sage: bijection.cur_path.insert(0, [])
            sage: bijection.cur_dims.insert(0, [0, 1])
            sage: bijection.cur_path[0].insert(0, [3])
            sage: bijection.next_state(3)
        """
        # Note that type D_n only contains the (absolute) values between 1 and n
        #   (unlike type A_n which is 1 to n+1). Thus we only need to subtract 1
        #   to match the indices.
        # Also note that we must subtract 1 from n to match the indices as well.
        n = self.n
        tableau_height = len(self.cur_path[0]) - 1

        # If it is a regular value, we follow the A_n rules
        if val > 0:
            KRTToRCBijectionTypeA.next_state(self, val)

            # If we are inserting n-1 (which will only update nu[a] for a <=
            #   n-1), we need to update the vacancy number of nu[n] as well
            if val == n - 1:
                self._update_vacancy_nums(val)
            if tableau_height >= n - 2:
                self._correct_vacancy_nums()
            return

        pos_val = -val

        if pos_val == n:
            # Special case for `\overline{n}` and adding to make height `n`
            #    where we only update the vacancy numbers
            # This only occurs with `r = n - 1`
            if self.cur_dims[0][0] == n - 1 and tableau_height == n - 1:
                self._update_vacancy_nums(n-2)
                self._update_vacancy_nums(n-1)
                self._correct_vacancy_nums()
                return

            if len(self.ret_rig_con[n - 1]) > 0:
                max_width = self.ret_rig_con[n - 1][0] + 1
            else:
                max_width = 1

            # Update the last one and skip a rigged partition
            max_width = self.ret_rig_con[n - 1].insert_cell(max_width)

            for a in reversed(range(tableau_height, n - 2)):
                max_width = self.ret_rig_con[a].insert_cell(max_width)
                self._update_vacancy_nums(a + 1)
                self._update_partition_values(a + 1)

            # Update all other effected remaining values

            self._update_vacancy_nums(n - 1)
            if tableau_height >= n - 2:
                self._correct_vacancy_nums()
            self._update_partition_values(n - 1)

            self._update_vacancy_nums(tableau_height)
            self._update_partition_values(tableau_height)

            if tableau_height > 0:
                self._update_vacancy_nums(tableau_height - 1)
                self._update_partition_values(tableau_height - 1)

            self._update_vacancy_nums(n - 2)
            return

        # Always add a cell to the first singular value in the first
        #   tableau we are updating.
        if len(self.ret_rig_con[pos_val - 1]) > 0:
            max_width = self.ret_rig_con[pos_val - 1][0] + 1
        else:
            max_width = 1
        # Special case for `\overline{n-1}` to take the larger of the last two
        if pos_val == n - 1 and len(self.ret_rig_con[n - 1]) > 0 and \
          self.ret_rig_con[n - 1][0] + 1 > max_width:
            max_width = self.ret_rig_con[n - 1][0] + 1

        # Add cells similar to type A_n but we move to the right until we reach
        #   the value of n-2
        for a in range(pos_val - 1, n - 2):
            max_width = self.ret_rig_con[a].insert_cell(max_width)

        # Handle the special behavior near values of n
        if tableau_height <= n - 2:
            max_width2 = self.ret_rig_con[n - 2].insert_cell(max_width)
            max_width = self.ret_rig_con[n - 1].insert_cell(max_width)
            if max_width2 < max_width:
                max_width = max_width2
        elif pos_val <= self.cur_dims[0][0]:
            # Special case when the height will become n
            max_width = self.ret_rig_con[self.cur_dims[0][0] - 1].insert_cell(max_width)

        # Go back following the regular A_n rules
        if tableau_height <= n - 3:
            max_width = self.ret_rig_con[n - 3].insert_cell(max_width)

        self._update_vacancy_nums(n - 2)
        self._update_vacancy_nums(n - 1)
        if tableau_height >= n - 2:
            self._correct_vacancy_nums()
        self._update_partition_values(n - 2)
        self._update_partition_values(n - 1)

        for a in reversed(range(tableau_height, n - 3)):
            max_width = self.ret_rig_con[a].insert_cell(max_width)
            self._update_vacancy_nums(a + 1)
            self._update_partition_values(a + 1)

        # Update the final rigged partitions
        if tableau_height < n - 2:
            self._update_vacancy_nums(tableau_height)
            self._update_partition_values(tableau_height)

            if pos_val <= tableau_height:
                for a in range(pos_val - 1, tableau_height):
                    self._update_vacancy_nums(a)
                    self._update_partition_values(a)
                if pos_val > 1:
                    self._update_vacancy_nums(pos_val-2)
                    self._update_partition_values(pos_val-2)
            elif 0 < tableau_height:
                self._update_vacancy_nums(tableau_height - 1)
                self._update_partition_values(tableau_height - 1)
        elif pos_val <= n - 1:
            for a in range(pos_val - 1, n - 2):
                self._update_vacancy_nums(a)
                self._update_partition_values(a)
            if pos_val > 1:
                self._update_vacancy_nums(pos_val-2)
                self._update_partition_values(pos_val-2)

class RCToKRTBijectionTypeE67(RCToKRTBijectionAbstract):
    r"""
    Specific implementation of the bijection from rigged configurations
    to tensor products of KR tableaux for type `E_{6,7}^{(1)}`.
    """
    def next_state(self, r):
        r"""
        Build the next state for type `E_{6,7}^{(1)}`.

        TESTS::

            sage: RC = RiggedConfigurations(['D', 4, 1], [[2, 1]])
            sage: from sage.combinat.rigged_configurations.bij_type_D import RCToKRTBijectionTypeD
            sage: bijection = RCToKRTBijectionTypeD(RC(partition_list=[[],[1,1],[1],[1]]))
            sage: bijection.next_state(0)
            1
        """
        r += 1 # Hack for now
        C = self.KRT.letters
        if r == 1:
            b = C.module_generators[0]
        elif r == 2: # remaining factor is r = 1
            b = C((-1, 2))
        elif r == 3: # remaining factor is r = 2
            b = C((-1, 3))
        elif r == 4: # remaining factor is r = 3
            b = C((-3, 4))
        elif r == 5: # remaining factor is r = 2??
            b = C((-2, 5))
        elif r == 6: # remaining factor is r = 1
            b = C((-1,6))

        last_size = 0
        pos = [x for x in b.value if x > 0]
        found = True
        while found:
            found = False
            assert len(pos) <= 2
            data = [self._find_singular_string(self.cur_partitions[a-1], last_size)
                    for a in pos]
            data = [(val, pos[i], self.cur_partitions[pos[i]-1][val])
                    for i,val in enumerate(data) if val is not None]
            if not data:
                break

            min_val = min(l for i,a,l in data)
            for i,a,l in data:
                if l == min_val:
                    found = True
                    last_size = l
                    self.cur_partitions[a-1].remove_cell(i)
                    b = b.f(a)
                    pos = [x for x in b.value if x > 0]
                    break

        for a in range(6):
            self._update_vacancy_numbers(a)
            for i in range(len(self.cur_partitions[a])):
                if self.cur_partitions[a].rigging[i] is None:
                    self.cur_partitions[a].rigging[i] = self.cur_partitions[a].vacancy_numbers[i]

        return(b)

