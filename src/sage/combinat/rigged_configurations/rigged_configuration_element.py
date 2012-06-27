r"""
A specific rigged configuration

A rigged configuration element is a sequence of :class:`RiggedPartition`
objects.

AUTHORS:

- Travis Scrimshaw (2010-09-26): Initial version

.. TODO:: Implement crystal operators `e_0` and `f_0`.

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

from sage.misc.abstract_method import abstract_method
from sage.structure.list_clone import ClonableArray
from sage.combinat.crystals.letters import CrystalOfLetters
from sage.combinat.rigged_configurations.rigged_partition import RiggedPartition
from sage.combinat.rigged_configurations.bijection import RCToKRTBijection

class RiggedConfigurationElement(ClonableArray):
    """
    The RiggedConfigurationElement class is a list of :class:`RiggedPartition`
    objects.

    For more information on rigged configurations, see
    :class:`RiggedConfigurations`.

    Typically to create a specific rigged configuration, the user will pass in
    the optional argument **partition_list** and if the user wants to specify
    the rigging values, give the optional argument **rigging_list** as well.
    If **rigging_list** is not passed, the rigging values are set to the
    corresponding vacancy numbers.

    EXAMPLES:

    Type `A_n^{(1)}` examples::

        sage: RC = RiggedConfigurations(['A', 4, 1], [[2, 2]])
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

        sage: RC = HighestWeightRiggedConfigurations(['A', 4, 1], [[1, 1], [1, 1]])
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

        sage: RC = RiggedConfigurations(['D', 4, 1], [[2, 2]])
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

        sage: RC = HighestWeightRiggedConfigurations(['D', 4, 1], [[1,1], [2, 1]])
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
    """

    def __init__(self, parent, *rigged_partitions, **options):
        r"""
        Construct a rigged configuration element.

        INPUT:

        - ``parent``            -- The parent of this element
        - ``rigged_partitions`` -- A list of rigged partitions

        There are two optional arguments to explicitly construct a rigged
        configuration. The first is **partition_list** which gives a list of
        partitions, and the second is **rigging_list** which is a list of
        corresponding lists of riggings. If only partition_list is specified,
        then it sets the rigging equal to the calculated vacancy numbers.

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
            n = parent._cartan_type.n
            if len(data) == 0:
                # Create a size n array of empty rigged tableau since no tableau
                #   were given
                nu = []
                for i in range(parent._cartan_type.n):
                    nu.append(RiggedPartition())
            else:
                if len(data) != n: # otherwise n should be equal to the number of tableaux
                    raise ValueError

                nu = []
                if "rigging_list" in options:
                    rigging_data = options["rigging_list"]

                    if len(rigging_data) != n:
                        raise ValueError

                    for i in range(parent._cartan_type.n):
                       nu.append(RiggedPartition(tuple(data[i]), \
                          list(rigging_data[i])))
                else:
                    for partition_data in data:
                        nu.append(RiggedPartition(tuple(partition_data)))
        elif len(list(rigged_partitions)) == 0:
            # Create a size n array of empty rigged tableau since no tableau
            #   were given
            nu = []
            for i in range(parent._cartan_type.n):
                nu.append(RiggedPartition())
        elif "KT_constructor" in options:
            # Used only by the Kleber tree
            # Not recommended to be called by the user since it avoids safety
            #   checks for speed
            data = options["KT_constructor"]
            shape_data = data[0]
            rigging_data = data[1]
            vac_data = data[2]
            nu = []
            for i in range(parent._cartan_type.n):
                nu.append(RiggedPartition(shape_data[i], rigging_data[i], vac_data[i]))
            ClonableArray.__init__(self, parent, nu)
            return
        elif parent._cartan_type.n == len(list(rigged_partitions)):
            ClonableArray.__init__(self, parent, list(rigged_partitions))
            return
        else:
            # Otherwise we did not receive any info, create a size n array of
            #   empty rigged partitions
            nu = []
            for i in range(parent._cartan_type.n):
                nu.append(RiggedPartition())

        # Set the vacancy numbers
        for a, partition in enumerate(nu):
            # If the partition is empty, there's nothing to do
            if len(partition) <= 0:
                continue

            # Setup the first block
            block_len = partition[0]
            vac_num = parent._calc_vacancy_number(nu, a, 0)

            for i, row_len in enumerate(partition):
                # If we've gone to a different sized block, then update the
                #   values which change when moving to a new block size
                if block_len != row_len:
                    vac_num = parent._calc_vacancy_number(nu, a, i)
                    block_len = row_len

                partition.vacancy_numbers[i] = vac_num
                if partition.rigging[i] is None:
                    partition.rigging[i] = partition.vacancy_numbers[i]

        ClonableArray.__init__(self, parent, nu)

    def _repr_(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: RC = HighestWeightRiggedConfigurations(['D', 4, 1], [[2, 2]])
            sage: RC(partition_list=[[2], [3,1], [3], [3]]) # indirect doctest
            <BLANKLINE>
            -1[ ][ ]-1
            <BLANKLINE>
            2[ ][ ][ ]2
            0[ ]0
            <BLANKLINE>
            -2[ ][ ][ ]-2
            <BLANKLINE>
            -2[ ][ ][ ]-2
            <BLANKLINE>
            sage: RC(partition_list=[[],[],[],[]]) # indirect doctest
            <BLANKLINE>
            (/)
            <BLANKLINE>
            (/)
            <BLANKLINE>
            (/)
            <BLANKLINE>
            (/)
            <BLANKLINE>
        """
        retStr = ""
        for tableau in self:
            retStr += "\n" + repr(tableau)
        return(retStr)

    def _latex_(self):
        """
        Return the LaTeX representation of ``self``.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['D', 4, 1], [[2, 2]])
            sage: latex(RC(partition_list=[[2], [3,1], [3], [3]])) # indirect doctest
            {
            \begin{array}[t]{r|c|c|l}
            \cline{2-3} -1 &\phantom{x}&\phantom{x}& -1 \\
            \cline{2-3}
            \end{array}
            }
            \quad
            {
            \begin{array}[t]{r|c|c|c|l}
            \cline{2-4} 2 &\phantom{x}&\phantom{x}&\phantom{x}& 2 \\
            \cline{2-4} 0 &\phantom{x}& \multicolumn{3}{l}{0} \\
            \cline{2-2}
            \end{array}
            }
            \quad
            {
            \begin{array}[t]{r|c|c|c|l}
            \cline{2-4} -2 &\phantom{x}&\phantom{x}&\phantom{x}& -2 \\
            \cline{2-4}
            \end{array}
            }
            \quad
            {
            \begin{array}[t]{r|c|c|c|l}
            \cline{2-4} -2 &\phantom{x}&\phantom{x}&\phantom{x}& -2 \\
            \cline{2-4}
            \end{array}
            }
            sage: latex(RC(partition_list=[[],[],[],[]])) # indirect doctest
            {\emptyset}
            \quad
            {\emptyset}
            \quad
            {\emptyset}
            \quad
            {\emptyset}
        """
        ret_string = self[0]._latex_()

        for partition in self[1:]:
            ret_string += "\n\quad\n" + partition._latex_()

        return ret_string

    def check(self):
        """
        Make sure all of the riggings are less than or equal to the vacancy number.

        TESTS::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2, 1]])
            sage: elt = RC(partition_list=[[1], [1], [], []])
            sage: elt.check()
        """
        for partition in self:
            for i, vac_num in enumerate(partition.vacancy_numbers):
                assert vac_num >= partition.rigging[i], "rigging can be at most the vacancy number"

    def to_Kirillov_Reshetikhin_tableaux(self, display_steps=False, display_all=False, **options):
        r"""
        Perform the bijection from this rigged configuration to a tensor
        product of Kirillov-Reshetikhin tableaux given in [RigConBijection]_
        for single boxes and with [BijectionLRT]_ and [BijectionDn]_ for
        multiple columns and rows.

        INPUT:

        - ``display_steps`` -- (default: False) Boolean which indicates if we want to output each step in the algorithm.
        - ``display_all``   -- (default: False) Boolean which indicates if we to output from each step, `tj`, and the column splittings.

        OUTPUT:

        - The KR tableaux corresponding to this rigged configuration.

        EXAMPLES::

            sage: RC = HighestWeightRiggedConfigurations(['A', 4, 1], [[2, 2]])
            sage: RC(partition_list=[[2], [2,2], [2], [2]]).to_Kirillov_Reshetikhin_tableaux()
            [[3, 3], [5, 5]]
            sage: RC = HighestWeightRiggedConfigurations(['D', 4, 1], [[2, 2]])
            sage: RC(partition_list=[[2], [2,2], [1], [1]]).to_Kirillov_Reshetikhin_tableaux()
            [[2, 3], [3, -2]]
        """
        #Letters = CrystalOfLetters(self.parent()._cartan_type.classical())
        Letters = CrystalOfLetters(self.parent()._cartan_type)
        n = self.parent()._cartan_type.n

        # Pass in a copy of our partitions since the bijection is destructive to it.
        bijection = RCToKRTBijection(self)

        # This is technically bad, but because the first thing we do is append
        #   an empty list to ret_crystal_path, we correct this. We do it this
        #   way so that we do not have to remove an empty list after the
        #   bijection has been performed.
        ret_crystal_path = []

        for tableau in self.parent().dims:
            ret_crystal_path.append([])

            split_column = False
            while len(bijection.rem_path[0]) > 0 or split_column:
                # If we have finished a splited column
                if len(bijection.rem_path[0]) == 0:
                    bijection.rem_path.pop(0)
                    split_column = False

                if display_steps or display_all:
                    print "===================="
                    print repr(self.parent()(*bijection.cur_partitions))
                    print "--------------------"
                    print ret_crystal_path
                    print "--------------------\n"

                # Build the next state
                # Check to see if we need to split off a column
                if len(bijection.rem_path[0][0]) > 1:
                    for row in bijection.rem_path[0]:
                        row.pop()
                    bijection.rem_path.insert(0, [[None]] * len(bijection.rem_path[0]))

                    # Perform the corresponding splitting map on rigged configurations
                    # All it does is update the vacancy numbers
                    for a in range(n):
                        bijection._update_vacancy_numbers(a)
                    split_column = True

                    if display_all:
                        print "Split column:"
                        print repr(self.parent()(*bijection.cur_partitions))
                        print "--------------------\n"

                bijection.rem_path[0].pop()

                # Check to see if we needed to pull off a box
                if len(bijection.rem_path[0]) > 0:
                    bijection.tj(len(bijection.rem_path[0]) + 1)

                if display_all:
                    print "tj:"
                    print repr(self.parent()(*bijection.cur_partitions))
                    print "--------------------\n"

                b = bijection.next_state()

                # Make sure we have a crystal letter
                ret_crystal_path[-1].append(Letters(b)) # Append the rank

            bijection.rem_path.pop(0)

        # If you're curious about this, see the note in AbstractTensorProductOfKRTableaux._highest_weight_iter().
        # You should never call this option.
        if "KRT_init_hack" in options:
            return options["KRT_init_hack"](pathlist=ret_crystal_path)

        #return self.parent()._bijection_class(self.parent()._cartan_type,
        return self.parent()._bijection_class(self.parent()._affine_ct,
          self.parent().dims)(pathlist=ret_crystal_path)

    def nu(self):
        r"""
        Return the list `\nu` of rigged partitions of this rigged
        configuration element.

        OUTPUT:

        - The `\nu` array as a list

        EXAMPLES::

            sage: RC = HighestWeightRiggedConfigurations(['A', 4, 1], [[2, 2]])
            sage: RC(partition_list=[[2], [2,2], [2], [2]]).nu()
            [0[ ][ ]0
            , -2[ ][ ]-2
            -2[ ][ ]-2
            , 2[ ][ ]2
            , -2[ ][ ]-2
            ]
        """
        return list(self)

    def e(self, a):
        r"""
        Action of the crystal operator `e_a` on this rigged configuration element.

        This implements the method defined in [CrysStructSchilling06]_ which
        finds the value `k` which is  the length of the string with the
        smallest negative rigging of smallest length. Then it removes a box
        from a string of length `k` in the `a`-th rigged partition, keeping all
        colabels fixed and increasing the new label by one. If no such string
        exists, then `e_a` is undefined.

        INPUT:

        - ``a`` -- The index of the partition to remove a box.

        OUTPUT:

        - The resulting rigged configuration element.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2,1]])
            sage: elt = RC(partition_list=[[1], [1], [1], [1]])
            sage: elt.e(3)
            sage: elt.e(1)
            <BLANKLINE>
            (/)
            <BLANKLINE>
            0[ ]0
            <BLANKLINE>
            0[ ]0
            <BLANKLINE>
            -1[ ]-1
            <BLANKLINE>
        """
        assert a in self.parent()._cartan_type.index_set()
        if a == 0:
            raise NotImplementedError("Only classical crystal operators implemented")

        a -= 1 # For indexing

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
        if k == 1:
            new_list.pop()
            new_vac_nums.pop()
            new_rigging.pop()
        else:
            new_list[rigging_index] -= 1
            cur_rigging += 1
            # Properly sort the riggings
            j = rigging_index + 1
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

                if k != 1: # If we did not remove a row
                    new_vac_nums[rigging_index] += 2

                new_partitions.append(RiggedPartition(new_list, new_rigging, new_vac_nums))

        from sage.combinat.rigged_configurations.rigged_configurations import RiggedConfigurations
        #ret_RC = self.__class__(RiggedConfigurations(self.parent()._cartan_type, self.parent().dims),
        ret_RC = self.__class__(RiggedConfigurations(self.parent()._affine_ct, self.parent().dims),
                              *new_partitions)
        if k != 1: # If we did not remove a row
            # Update that row's vacancy number
            ret_RC[a].vacancy_numbers[rigging_index] = \
              self.parent()._calc_vacancy_number(ret_RC.nu(), a, rigging_index)
        return(ret_RC)

    def _generate_partition_e(self, a, b, k):
        r"""
        Generate a new partition for a given value of `a` by updating the
        vacancy numbers and preserving co-labels for the map `e_a`.

        INPUT:

        - ``a`` -- The index of the partition we operated on.
        - ``b`` -- The index of the partition to generate.
        - ``k`` -- The length of the string with the smallest negative rigging of smallest length.

        OUTPUT:

        - The constructed rigged partition.

        TESTS::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2,1]])
            sage: RC(partition_list=[[1], [1], [1], [1]])._generate_partition_e(1, 2, 1)
            -1[ ]-1
            <BLANKLINE>
        """
        #cartan_matrix = self.parent()._cartan_type.classical().cartan_matrix()
        cartan_matrix = self.parent()._cartan_type.cartan_matrix()
        # Check to make sure we will do something
        if cartan_matrix[a][b] == 0:
            return self[b]

        new_list = self[b][:]
        new_vac_nums = self[b].vacancy_numbers[:]
        new_rigging = self[b].rigging[:]

        # Update the vacancy numbers and the rigging
        value = cartan_matrix[a][b]
        for i in range(len(new_vac_nums)):
            if new_list[i] < k:
                break

            new_vac_nums[i] += value
            new_rigging[i] += value

        return(RiggedPartition(new_list, new_rigging, new_vac_nums))

    def f(self, a):
        r"""
        Action of crystal operator `f_a` on this rigged configuration element.

        This implements the method defined in [CrysStructSchilling06]_ which
        finds the value `k` which is  the length of the string with the
        smallest nonpositive rigging of largest length. Then it adds a box from
        a string of length `k` in the `a`-th rigged partition, keeping all
        colabels fixed and decreasing the new label by one. If no such string
        exists, then it adds a new string of length 1 with label `-1`. If any
        of the resulting vacancy numbers are larger than the labels (i.e. it
        is an invalid rigged configuration), then `f_a` is undefined.

        INPUT:

        - ``a`` -- The index of the partition to add a box.

        OUTPUT:

        - The resulting rigged configuration element.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2,1]])
            sage: elt = RC(partition_list=[[1], [1], [1], [1]])
            sage: elt.f(1)
            sage: elt.f(2)
            <BLANKLINE>
            0[ ]0
            <BLANKLINE>
            -1[ ]-1
            -1[ ]-1
            <BLANKLINE>
            1[ ]1
            <BLANKLINE>
            -1[ ]-1
            <BLANKLINE>
        """
        assert a in self.parent()._cartan_type.index_set()
        if a == 0:
            raise NotImplementedError("Only classical crystal operators implemented")

        a -= 1 # For indexing

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
          self.parent()._calc_vacancy_number(new_partitions, a, add_index)
        if new_partitions[a].rigging[add_index] > new_partitions[a].vacancy_numbers[add_index]:
            return None

        # Note that we do not need to sort the rigging since if there was a
        #   smaller rigging in a larger row, then `k` would be larger.
        from sage.combinat.rigged_configurations.rigged_configurations import RiggedConfigurations
        #return self.__class__(RiggedConfigurations(self.parent()._cartan_type, self.parent().dims),
        return self.__class__(RiggedConfigurations(self.parent()._affine_ct, self.parent().dims),
                              *new_partitions)

    def _generate_partition_f(self, a, b, k):
        r"""
        Generate a new partition for a given value of `a` by updating the
        vacancy numbers and preserving co-labels for the map `f_a`.

        INPUT:

        - ``a`` -- The index of the partition we operated on.
        - ``b`` -- The index of the partition to generate.
        - ``k`` -- The length of the string with smallest nonpositive rigging of largest length.

        OUTPUT:

        - The constructed rigged partition.

        TESTS::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2,1]])
            sage: RC(partition_list=[[1], [1], [1], [1]])._generate_partition_f(1, 2, 1)
            0[ ]0
            <BLANKLINE>
        """
        #cartan_matrix = self.parent()._cartan_type.classical().cartan_matrix()
        cartan_matrix = self.parent()._cartan_type.cartan_matrix()
        # Check to make sure we will do something
        if cartan_matrix[a][b] == 0:
            return self[b]

        new_list = self[b][:]
        new_vac_nums = self[b].vacancy_numbers[:]
        new_rigging = self[b].rigging[:]

        # Update the vacancy numbers and the rigging
        value = cartan_matrix[a][b]
        for i in range(len(new_vac_nums)):
            if new_list[i] <= k:
                break

            new_vac_nums[i] -= value
            new_rigging[i] -= value

        return(RiggedPartition(new_list, new_rigging, new_vac_nums))

    def cc(self):
        r"""
        Compute the cocharge statistic.

        Computes the cocharge statistic [CrysStructSchilling06]_ on this
        rigged configuration `(\nu, J)`. The cocharge statistic is defined as:

        .. MATH::

            cc(\nu, J) = \frac{1}{2} \sum_{a, b \in J}
            \sum_{j,k \geq 1} \left( \alpha_a \mid \alpha_b \right)
            \min(j, k) m_j^{(a)} m_k^{(b)}
            + \sum_{a, i} \left\lvert J^{(a, i)} \right\rvert.

        EXAMPLES::

            sage: HWRC = HighestWeightRiggedConfigurations(['A', 3, 1], [[3, 2], [2,1], [1,1]])
            sage: HWRC(partition_list=[[1], [1], []]).cc()
            1
        """
        cc = 0
        rigging_sum = 0
        num_partitions = len(self)
        #cartan_matrix = self.parent().cartan_type().classical().cartan_matrix()
        cartan_matrix = self.parent().cartan_type().cartan_matrix()

        for a in range(num_partitions):
            num_rows = len(self[a])
            for i in range(num_rows):
                rigging_sum += self[a].rigging[i]

                for b in range(a, num_partitions):
                    temp_cc = 0
                    num_rows_b = len(self[b])
                    # TODO - Optimize this by using the fact that a partition is weakly decreasing
                    for j in range(num_rows_b):
                        temp_cc += min(self[a][i], self[b][j])

                    if a != b:
                        cc += (cartan_matrix[a][b] + cartan_matrix[b][a]) * temp_cc
                    else: # Assuming A_{ii} = 2 for all types
                        cc += 2 * temp_cc

        return cc / 2 + rigging_sum

    def get_vacancy_numbers(self, a):
        r"""
        Return the list of all vacancy numbers of the rigged partition
        `\nu^{(a)}` (with duplicates).

        INPUT:

        - ``a`` -- The index of the rigged partition.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2, 2]])
            sage: RC(partition_list=[[1], [2,1], [1], []]).get_vacancy_numbers(2)
            [-2, -1]
        """
        return self[a-1].vacancy_numbers # -1 for indexing

    def get_vacancy_number(self, a, i):
        r"""
        Return the vacancy number `p_i^{(a)}`.

        INPUT:

        - ``a`` -- The index of the rigged partition.
        - ``i`` -- The row of the rigged partition.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2, 2]])
            sage: elt = RC(partition_list=[[1], [2,1], [1], []])
            sage: elt.get_vacancy_number(2, 3)
            sage: elt.get_vacancy_number(2, 2)
            -2
            sage: elt.get_vacancy_number(2, 1)
            -1
        """
        partition = self[a-1] # -1 for indexing
        for k, val in enumerate(partition):
            if val == i:
                return partition.vacancy_numbers[k]
            elif val < i:
                return None

        return None

