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
from sage.combinat.crystals.letters import CrystalOfLetters
from sage.misc.cachefunc import cached_method

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
        end = self._endpoint(self.cur_dims[0][0])
        max_width = float("inf")
        found = (b != end)
        while found:
            found = False
            data = [(-a, find_singular_string(self.ret_rig_con[-a-1], max_width))
                    for a in b.value if a < 0]
            if not data:
                break

            max_val = max(l for a,l in data)
            for a,l in data:
                if l == max_val:
                    self.ret_rig_con[a-1].insert_cell(max_width)
                    max_width = l
                    b = b.e(a)
                    found = (b != end)
                    break

        for a in range(6):
            self._update_vacancy_nums(a)
            self._update_partition_values(a)

    def _next_index(self, r, target):
        """
        Return the next index after ``r`` when performing a step
        in the bijection going towards ``target``.

        TESTS::

            sage: KRT = crystals.TensorProductOfKirillovReshetikhinTableaux(['E', 6, 1], [[5,1]])
            sage: from sage.combinat.rigged_configurations.bij_type_E67 import KRTToRCBijectionTypeE67
            sage: bijection = KRTToRCBijectionTypeE67(KRT.module_generators[0])
            sage: bijection._next_index(3, 5)
            2
            sage: bijection._next_index(2, 5)
            5
            sage: bijection._next_index(3, 4)
            4
            sage: bijection._next_index(1, 5)
            3
            sage: bijection._next_index(1, 4)
            3
            sage: bijection._next_index(1, 6)
            6
        """
        #       3-2-5
        #      / \
        # 0 - 1   4
        #      \
        #       6
        if self.tp_krt.cartan_type().classical().rank() == 6:
            if r == 0:
                return 1
            if r == 1:
                if target == 6:
                    return 6
                return 3
            if r == 3:
                if target == 4:
                    return 4
                return 2
            if r == 2:
                return 5
        else: # rank == 7
            if r == 0:
                return 7

    @cached_method
    def _endpoint(self, r):
        r"""
        Return the endpoint for the bijection in type `E_6^{(1)}`.
        """
        if self.tp_krt.cartan_type().classical().rank() == 6:
            return _endpoint6(r)
        else:
            return _endpoint7(r)

class RCToKRTBijectionTypeE67(RCToKRTBijectionAbstract):
    r"""
    Specific implementation of the bijection from rigged configurations
    to tensor products of KR tableaux for type `E_{6,7}^{(1)}`.
    """
    def next_state(self, r):
        r"""
        Build the next state for type `E_{6,7}^{(1)}`.

        TESTS::

            sage: RC = RiggedConfigurations(['E', 6, 1], [[2, 1]])
            sage: from sage.combinat.rigged_configurations.bij_type_E67 import RCToKRTBijectionTypeE67
            sage: bijection = RCToKRTBijectionTypeE67(RC(partition_list=[[1],[1,1],[1,1],[1,1,1], [1,1], [1]]))
            sage: bijection.next_state(1)
            (-2, 1)
        """
        last_size = 0
        found = True
        b = self._endpoint(r)
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

    def _next_index(self, r):
        """
        Return the next index after ``r`` when performing a step
        in the bijection.

        TESTS::

            sage: RC = RiggedConfigurations(['E', 6, 1], [[2, 1]])
            sage: from sage.combinat.rigged_configurations.bij_type_E67 import RCToKRTBijectionTypeE67
            sage: bijection = RCToKRTBijectionTypeE67(RC(partition_list=[[1],[1,1],[1,1],[1,1,1], [1,1], [1]]))
            sage: bijection._next_index(2)
            3
        """
        if self.KRT.cartan_type().classical().rank() == 6:
            if r == 1:
                return 0
            if r == 2:
                return 3
            if r == 3:
                return 1
            if r == 4:
                return 3
            if r == 5:
                return 2
            if r == 6:
                return 1
        else: # rank == 7
            if r == 7:
                return 0

    @cached_method
    def _endpoint(self, r):
        r"""
        Return the endpoint for the bijection in type `E_6^{(1)}`.
        """
        if self.KRT.cartan_type().classical().rank() == 6:
            return _endpoint6(r)
        else:
            return _endpoint7(r)

def _endpoint6(r):
    r"""
    Return the endpoint for the bijection in type `E_6^{(1)}`.
    """
    C = CrystalOfLetters(['E',6])
    if r == 1:
        return C.module_generators[0]  # C((1,))
    elif r == 2:
        return C((-3, 2))
    elif r == 3:
        return C((-1, 3))
    elif r == 4:
        return C((-3, 4))
    elif r == 5:
        return C((-2, 5))
    elif r == 6:
        return C((-1,6))

