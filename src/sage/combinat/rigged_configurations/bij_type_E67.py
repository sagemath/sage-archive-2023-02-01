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

        def find_singular_string(p, max_width):
            max_pos = -1
            if max_width > 0:
                for i, vac_num in enumerate(p.vacancy_numbers):
                    if p[i] <= max_width and vac_num == p.rigging[i]:
                        max_pos = i
                        break
            if max_pos == -1:
                return 0
            return p[max_pos]

        b = self.tp_krt.parent().letters(val)
        max_width = float("inf")
        found = True
        while found:
            found = False
            data = [(-a, find_singular_string(self.ret_rig_con[-a-1], max_width))
                    for a in b.value if a < 0]
            if not data:
                break

            max_val = max(l for a,l in data)
            for a,l in data:
                if l == max_val:
                    found = True
                    self.ret_rig_con[a-1].insert_cell(max_width)
                    max_width = l
                    b = b.e(a)
                    break

        for a in range(6):
            self._update_vacancy_nums(a)
            self._update_partition_values(a)

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
        found = True
        while found:
            found = False
            data = [(a, self._find_singular_string(self.cur_partitions[a-1], last_size))
                    for a in b.value if a > 0]
            data = [(val, a, self.cur_partitions[a-1][val])
                    for a,val in data if val is not None]
            if not data:
                break

            min_val = min(l for i,a,l in data)
            for i,a,l in data:
                if l == min_val:
                    found = True
                    last_size = l
                    self.cur_partitions[a-1].remove_cell(i)
                    b = b.f(a)
                    break

        for a in range(6):
            self._update_vacancy_numbers(a)
            for i in range(len(self.cur_partitions[a])):
                if self.cur_partitions[a].rigging[i] is None:
                    self.cur_partitions[a].rigging[i] = self.cur_partitions[a].vacancy_numbers[i]

        return(b)

