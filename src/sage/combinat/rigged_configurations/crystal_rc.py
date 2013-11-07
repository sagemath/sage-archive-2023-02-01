r"""
Crystal of Rigged Configurations

AUTHORS:

- Travis Scrimshaw (2010-09-26): Initial version

We only consider the highest weight crystal structure, not the
Kirillov-Reshetikhin structure, and we extend this to all types.

INPUT:

- ``cartan_type`` -- A Cartan type

- ``wt`` -- the highest weight in the weight lattice
        
EXAMPLES::
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
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.structure.list_clone import ClonableArray
from sage.categories.highest_weight_crystals import HighestWeightCrystals
from sage.categories.regular_crystals import RegularCrystals
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.rigged_configurations.rigged_configuration_element import RiggedConfigurationElement
from sage.combinat.rigged_configurations.rigged_partition import RiggedPartition, RiggedPartitionTypeB

# Note on implementation, this class is used for simply-laced types only
class CrystalOfRiggedConfigurations(Parent, UniqueRepresentation):
    r"""
    A highest weight crystal of rigged configurations.

    The crystal structure is given in [CrysStructSchilling06]_.

    INPUT:

    - ``cartan_type`` -- a (simply-laced) Cartan type

    - ``wt`` -- the highest weight vector in the weight lattice
    """
    @staticmethod
    def __classcall_private__(cls, cartan_type, wt):
        r"""
        Normalize the input arguments to ensure unique representation.

        EXAMPLES::
        """
        cartan_type = CartanType(cartan_type)
        wt_lattice = cartan_type.root_system().weight_lattice()
        wt = wt_lattice(wt)
        return super(CrystalOfRiggedConfigurations, cls).__classcall__(cls, cartan_type, wt)

    def __init__(self, cartan_type, wt):
        r"""
        Initialize ``self``.

        EXAMPLES::
        """
        self._cartan_type = cartan_type
        self._wt = wt
        # We store the cartan matrix for the vacancy number calculations for speed
        self._cartan_matrix = self._cartan_type.cartan_matrix()
        Parent.__init__(self, category=(RegularCrystals(), HighestWeightCrystals()))
        self.module_generators = (self.element_class( self,
                                  partition_list=[[] for i in cartan_type.index_set()] ),)

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
        """
        vac_num = 0
        vac_num += self._wt[a]

        a = list(self.index_set()).index(a)
        row_len = partitions[a][i]

        for b, value in enumerate(self._cartan_matrix.row(a)):
            vac_num -= value * partitions[b].get_num_cells_to_column(row_len)

        return vac_num

    class Element(RiggedConfigurationElement):
        def __init__(self, parent, rigged_partitions=[], **options):
            r"""
            Construct a rigged configuration element.

            EXAMPLES::

                sage: RC = RiggedConfigurations(['A', 4, 1], [[2, 1]])
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
                sage: RC(partition_list=[[1], [1], [], []])
                <BLANKLINE>
                -1[ ]-1
                <BLANKLINE>
                0[ ]0
                <BLANKLINE>
                (/)
                <BLANKLINE>
                (/)
                <BLANKLINE>
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
            if "partition_list" in options:
                data = options["partition_list"]
                n = parent._cartan_type.rank()
                if len(data) == 0:
                    # Create a size n array of empty rigged tableau since no tableau
                    #   were given
                    nu = []
                    for i in range(n):
                        nu.append(RiggedPartition())
                else:
                    if len(data) != n: # otherwise n should be equal to the number of tableaux
                        raise ValueError("Incorrect number of partitions")

                    nu = []
                    if "rigging_list" in options:
                        rigging_data = options["rigging_list"]

                        if len(rigging_data) != n:
                            raise ValueError("Incorrect number of riggings")

                        for i in range(n):
                           nu.append(RiggedPartition(tuple(data[i]), \
                              list(rigging_data[i])))
                    else:
                        for partition_data in data:
                            nu.append(RiggedPartition(tuple(partition_data)))
            elif parent._cartan_type.rank() == len(rigged_partitions) and \
                isinstance(rigged_partitions[0], RiggedPartition):
                # The isinstance check is to make sure we are not in the n == 1 special case because
                #   Parent's __call__ always passes at least 1 argument to the element constructor

                # Special display case
                if parent.cartan_type().type() == 'B':
                    rigged_partitions[-1] = RiggedPartitionTypeB(rigged_partitions[-1])
                ClonableArray.__init__(self, parent, rigged_partitions)
                return
            else:
                # Otherwise we did not receive any info, create a size n array of
                #   empty rigged partitions
                nu = []
                for i in range(parent._cartan_type.rank()):
                    nu.append(RiggedPartition())
                #raise ValueError("Invalid input")
                #raise ValueError("Incorrect number of rigged partitions")

            # Set the vacancy numbers
            for a, partition in enumerate(nu):
                # If the partition is empty, there's nothing to do
                if len(partition) <= 0:
                    continue

                # Setup the first block
                block_len = partition[0]
                index_set = self.parent().index_set()
                vac_num = parent._calc_vacancy_number(nu, index_set[a], 0)

                for i, row_len in enumerate(partition):
                    # If we've gone to a different sized block, then update the
                    #   values which change when moving to a new block size
                    if block_len != row_len:
                        vac_num = parent._calc_vacancy_number(nu, index_set[a], i)
                        block_len = row_len

                    partition.vacancy_numbers[i] = vac_num
                    if partition.rigging[i] is None:
                        partition.rigging[i] = partition.vacancy_numbers[i]

            # Special display case
            if parent.cartan_type().type() == 'B':
                nu[-1] = RiggedPartitionTypeB(nu[-1])

            ClonableArray.__init__(self, parent, nu)

        def e(self, a):
            r"""
            Action of the crystal operator `e_a` on ``self``.
            """
            index_set = list(self.index_set())
            a = index_set.index(a)
            new_list = self[a][:]
            new_vac_nums = self[a].vacancy_numbers[:]
            new_rigging = self[a].rigging[:]

            # Find k and perform e_a
            k = None
            num_rows = len(new_list)
            cur_rigging = -1
            rigging_index = None
            for i in range(num_rows):
                if new_rigging[i] <= cur_rigging:
                    cur_rigging = new_rigging[i]
                    rigging_index = i

            # If we've not found a valid k
            if rigging_index is None:
                return None

            # Note that because the riggings are weakly decreasing, we will always
            #   remove the last box on of a block
            k = new_list[rigging_index]
            set_vac_num = False
            if k == 1:
                new_list.pop()
                new_vac_nums.pop()
                new_rigging.pop()
            else:
                new_list[rigging_index] -= 1
                cur_rigging += 1
                # Properly sort the riggings
                j = rigging_index + 1
                # Update the vacancy number if the row lengths are the same
                if j < num_rows and new_list[j] == new_list[rigging_index]:
                    new_vac_nums[rigging_index] = new_vac_nums[j]
                    set_vac_num = True
                while j < num_rows and new_list[j] == new_list[rigging_index] \
                  and new_rigging[j] > cur_rigging:
                    new_rigging[j-1] = new_rigging[j] # Shuffle it along
                    j += 1
                new_rigging[j-1] = cur_rigging

            new_partitions = []
            for b in range(len(self)):
                if b != a:
                    new_partitions.append(self._generate_partition_e(a, b, k))
                else:
                    # Update the vacancy numbers and the rigging
                    for i in range(len(new_vac_nums)):
                        if new_list[i] < k:
                            break

                        new_vac_nums[i] += 2
                        new_rigging[i] += 2

                    
                    if k != 1 and not set_vac_num: # If we did not remove a row nor found another row of length k-1
                        new_vac_nums[rigging_index] += 2

                    new_partitions.append(RiggedPartition(new_list, new_rigging, new_vac_nums))

            ret_RC = self.__class__(self.parent(), new_partitions)
            if k != 1 and not set_vac_num: # If we did not remove a row nor found another row of length k-1
                # Update that row's vacancy number
                ret_RC[a].vacancy_numbers[rigging_index] = \
                  self.parent()._calc_vacancy_number(ret_RC.nu(), index_set[a], rigging_index)
            return(ret_RC)

        def f(self, a):
            r"""
            Action of crystal operator `f_a` on ``self``.
            """
            index_set = list(self.index_set())
            a = index_set.index(a)

            new_list = self[a][:]
            new_vac_nums = self[a].vacancy_numbers[:]
            new_rigging = self[a].rigging[:]

            # Find k and perform f_a
            k = None
            add_index = -1 # Index where we will add our row too
            rigging_index = None # Index which we will pull the rigging from 
            cur_rigging = 0
            num_rows = len(new_list)
            for i in reversed(range(num_rows)):
                # If we need to increment a row, look for when we change rows for
                #   the correct index.
                if add_index is None and new_list[i] != new_list[rigging_index]:
                    add_index = i+1

                if new_rigging[i] <= cur_rigging:
                    cur_rigging = new_rigging[i]
                    k = new_list[i]
                    rigging_index = i
                    add_index = None

            # If we've not found a valid k
            if k is None:
                new_list.append(1)
                new_rigging.append(-1)
                new_vac_nums.append(None)
                k = 0
                add_index = num_rows
                num_rows += 1 # We've added a row
            else:
                if add_index is None: # We are adding to the first row in the list
                    add_index = 0
                new_list[add_index] += 1
                new_rigging.insert(add_index, new_rigging[rigging_index] - 1)
                new_vac_nums.insert(add_index, None)
                new_rigging.pop(rigging_index + 1) # add 1 for the insertion
                new_vac_nums.pop(rigging_index + 1)

            new_partitions = []
            for b in range(len(self)):
                if b != a:
                    new_partitions.append(self._generate_partition_f(a, b, k))
                else:
                    # Update the vacancy numbers and the rigging
                    for i in range(num_rows):
                        if new_list[i] <= k:
                            break

                        if i != add_index:
                            new_vac_nums[i] -= 2
                            new_rigging[i] -= 2

                    new_partitions.append(RiggedPartition(new_list, new_rigging, new_vac_nums))

            new_partitions[a].vacancy_numbers[add_index] = \
              self.parent()._calc_vacancy_number(new_partitions, index_set[a], add_index)
            if new_partitions[a].rigging[add_index] > new_partitions[a].vacancy_numbers[add_index]:
                return None

            # Note that we do not need to sort the rigging since if there was a
            #   smaller rigging in a larger row, then `k` would be larger.
            return self.__class__(self.parent(), new_partitions)

