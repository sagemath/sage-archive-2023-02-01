r"""
Rigged Configuration Elements

A rigged configuration element is a sequence of :class:`RiggedPartition`
objects.

AUTHORS:

- Travis Scrimshaw (2010-09-26): Initial version
- Travis Scrimshaw (2012-10-25): Added virtual rigged confingurations
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
from sage.structure.list_clone import ClonableArray
from sage.combinat.rigged_configurations.rigged_partition import RiggedPartition, \
  RiggedPartitionTypeB

class RiggedConfigurationElement(ClonableArray):
    """
    A rigged configuration for simply-laced types.

    For more information on rigged configurations, see
    :class:`RiggedConfigurations`. For rigged configurations for
    non-simply-laced types, use :class:`RCNonSimplyLacedElement`.

    Typically to create a specific rigged configuration, the user will pass in
    the optional argument **partition_list** and if the user wants to specify
    the rigging values, give the optional argument **rigging_list** as well.
    If **rigging_list** is not passed, the rigging values are set to the
    corresponding vacancy numbers.

    INPUT:

    - ``parent`` -- the parent of this element

    - ``rigged_partitions`` -- a list of rigged partitions

    There are two optional arguments to explicitly construct a rigged
    configuration. The first is **partition_list** which gives a list of
    partitions, and the second is **rigging_list** which is a list of
    corresponding lists of riggings. If only partition_list is specified,
    then it sets the rigging equal to the calculated vacancy numbers.

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

        sage: RC = RiggedConfigurations(['A', 4, 1], [[1, 1], [1, 1]])
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

        sage: RC = RiggedConfigurations(['D', 4, 1], [[1,1], [2, 1]])
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

    We can go between
    :class:`tensor products of KR tableaux<TensorProductOfKirillovReshetikhinTableaux>`
    and tensor products of
    :mod:`KR crystals <sage.combinat.crystals.kirillov_reshetikhin>`::

        sage: RC = RiggedConfigurations(['D', 4, 1], [[1,1], [2,1]])
        sage: rc_elt = RC(partition_list=[[1], [1,1], [1], [1]])
        sage: tp_krtab = rc_elt.to_tensor_product_of_kirillov_reshetikhin_tableaux(); tp_krtab
        [[-2]] (X) [[1], [2]]
        sage: tp_krcrys = rc_elt.to_tensor_product_of_kirillov_reshetikhin_crystals(); tp_krcrys
        [[[-2]], [[1], [2]]]
        sage: tp_krcrys == tp_krtab.to_tensor_product_of_kirillov_reshetikhin_crystals()
        True
        sage: RC(tp_krcrys) == rc_elt
        True
        sage: RC(tp_krtab) == rc_elt
        True
        sage: tp_krtab.to_rigged_configuration() == rc_elt
        True
    """
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
        if "KT_constructor" in options:
            # Used only by the Kleber tree
            # Not recommended to be called by the user since it avoids safety
            #   checks for speed
            data = options["KT_constructor"]
            shape_data = data[0]
            rigging_data = data[1]
            vac_data = data[2]
            nu = []
            for i in range(parent._cartan_type.classical().rank()):
                nu.append(RiggedPartition(shape_data[i], rigging_data[i], vac_data[i]))
            # Special display case
            if parent.cartan_type().type() == 'B':
                nu[-1] = RiggedPartitionTypeB(nu[-1])
            ClonableArray.__init__(self, parent, nu)
            return
        elif "partition_list" in options:
            data = options["partition_list"]
            n = parent._cartan_type.classical().rank()
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
        elif parent._cartan_type.classical().rank() == len(rigged_partitions) and \
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
            for i in range(parent._cartan_type.classical().rank()):
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

        # Special display case
        if parent.cartan_type().type() == 'B':
            nu[-1] = RiggedPartitionTypeB(nu[-1])

        ClonableArray.__init__(self, parent, nu)

    def _repr_(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['D', 4, 1], [[2, 2]])
            sage: RC(partition_list=[[2], [3,1], [3], [3]])
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
            sage: RC(partition_list=[[],[],[],[]])
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
        ret_str = ""
        for tableau in self:
            ret_str += "\n" + repr(tableau)
        return(ret_str)

    def _latex_(self):
        r"""
        Return the LaTeX representation of ``self``.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['D', 4, 1], [[2, 2]])
            sage: latex(RC(partition_list=[[2], [3,1], [3], [3]]))
            {
            \begin{array}[t]{r|c|c|l}
            \cline{2-3} -1 &\phantom{|}&\phantom{|}& -1 \\
             \cline{2-3} 
            \end{array}
            } 
            \quad
             {
            \begin{array}[t]{r|c|c|c|l}
            \cline{2-4} 2 &\phantom{|}&\phantom{|}&\phantom{|}& 2 \\
             \cline{2-4} 0 &\phantom{|}& \multicolumn{3 }{l}{ 0 } \\
             \cline{2-2} 
            \end{array}
            } 
            \quad
             {
            \begin{array}[t]{r|c|c|c|l}
            \cline{2-4} -2 &\phantom{|}&\phantom{|}&\phantom{|}& -2 \\
             \cline{2-4} 
            \end{array}
            } 
            \quad
             {
            \begin{array}[t]{r|c|c|c|l}
            \cline{2-4} -2 &\phantom{|}&\phantom{|}&\phantom{|}& -2 \\
             \cline{2-4} 
            \end{array}
            }
            sage: latex(RC(partition_list=[[],[],[],[]]))
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

    def _ascii_art_(self):
        """
        Return an ASCII art representation of ``self``.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['D', 4, 1], [[2, 2]])
            sage: ascii_art(RC(partition_list=[[2], [3,1], [3], [3]]))
            -1[ ][ ]-1  2[ ][ ][ ]2  -2[ ][ ][ ]-2  -2[ ][ ][ ]-2
                        0[ ]0
            sage: ascii_art(RC(partition_list=[[],[],[],[]]))
            (/)  (/)  (/)  (/)
            sage: RC = RiggedConfigurations(['D', 7, 1], [[3,3],[5,2],[4,3],[2,3],[4,4],[3,1],[1,4],[2,2]])
            sage: elt = RC(partition_list=[[2],[3,2,1],[2,2,1,1],[2,2,1,1,1,1],[3,2,1,1,1,1],[2,1,1],[2,2]],
            ....:          rigging_list=[[2],[1,0,0],[4,1,2,1],[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,0],[0,0]])
            sage: ascii_art(elt)
            3[ ][ ]2  1[ ][ ][ ]1  4[ ][ ]4  2[ ][ ]1  0[ ][ ][ ]0  0[ ][ ]0  0[ ][ ]0
                      2[ ][ ]0     4[ ][ ]1  2[ ][ ]0  2[ ][ ]1     0[ ]0     0[ ][ ]0
                      1[ ]0        3[ ]2     0[ ]0     0[ ]0        0[ ]0
                                   3[ ]1     0[ ]0     0[ ]0
                                             0[ ]0     0[ ]0
                                             0[ ]0     0[ ]0
            sage: Partitions.global_options(convention='French')
            sage: ascii_art(elt)
                                             0[ ]0     0[ ]0
                                             0[ ]0     0[ ]0
                                   3[ ]1     0[ ]0     0[ ]0
                      1[ ]0        3[ ]2     0[ ]0     0[ ]0        0[ ]0
                      2[ ][ ]0     4[ ][ ]1  2[ ][ ]0  2[ ][ ]1     0[ ]0     0[ ][ ]0
            3[ ][ ]2  1[ ][ ][ ]1  4[ ][ ]4  2[ ][ ]1  0[ ][ ][ ]0  0[ ][ ]0  0[ ][ ]0
            sage: Partitions.global_options.reset()
        """
        from sage.combinat.partition import PartitionOptions
        if PartitionOptions['convention'] == "French":
            baseline = lambda s: 0
        else:
            baseline = lambda s: len(s)
        from sage.misc.ascii_art import AsciiArt
        s = repr(self[0]).splitlines()
        ret = AsciiArt(s, baseline=baseline(s))
        for tableau in self[1:]:
            s = repr(tableau).splitlines()
            ret += AsciiArt(["  "], baseline=baseline(s)) + AsciiArt(s, baseline=baseline(s))
        return ret

    def check(self):
        """
        Make sure all of the riggings are less than or equal to the
        vacancy number.

        TESTS::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2, 1]])
            sage: elt = RC(partition_list=[[1], [1], [], []])
            sage: elt.check()
        """
        for partition in self:
            for i, vac_num in enumerate(partition.vacancy_numbers):
                if vac_num < partition.rigging[i]:
                    raise ValueError("rigging can be at most the vacancy number")

    def to_tensor_product_of_kirillov_reshetikhin_tableaux(self, display_steps=False):
        r"""
        Perform the bijection from this rigged configuration to a tensor
        product of Kirillov-Reshetikhin tableaux given in [RigConBijection]_
        for single boxes and with [BijectionLRT]_ and [BijectionDn]_ for
        multiple columns and rows.

        .. NOTE::

            This is only proven to be a bijection in types `A_n^{(1)}`
            and `D_n^{(1)}`, as well as `\bigotimes_i B^{r_i,1}` and
            `\bigotimes_i B^{1,s_i}` for general affine types.

        INPUT:

        - ``display_steps`` -- (default: ``False``) boolean which indicates
          if we want to output each step in the algorithm

        OUTPUT:

        - The tensor product of KR tableaux element corresponding to this
          rigged configuration.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2, 2]])
            sage: RC(partition_list=[[2], [2,2], [2], [2]]).to_tensor_product_of_kirillov_reshetikhin_tableaux()
            [[3, 3], [5, 5]]
            sage: RC = RiggedConfigurations(['D', 4, 1], [[2, 2]])
            sage: elt = RC(partition_list=[[2], [2,2], [1], [1]])
            sage: tp_krt = elt.to_tensor_product_of_kirillov_reshetikhin_tableaux(); tp_krt
            [[2, 3], [3, -2]]

        This is invertible by calling
        :meth:`~sage.combinat.rigged_configurations.tensor_product_kr_tableaux_element.TensorProductOfKirillovReshetikhinTableauxElement.to_rigged_configuration()`::

            sage: ret = tp_krt.to_rigged_configuration(); ret
            <BLANKLINE>
            0[ ][ ]0
            <BLANKLINE>
            -2[ ][ ]-2
            -2[ ][ ]-2
            <BLANKLINE>
            0[ ]0
            <BLANKLINE>
            0[ ]0
            <BLANKLINE>
            sage: elt == ret
            True
        """
        from sage.combinat.rigged_configurations.bijection import RCToKRTBijection
        return RCToKRTBijection(self).run(display_steps)

    def to_tensor_product_of_kirillov_reshetikhin_crystals(self, display_steps=False):
        r"""
        Return the corresponding tensor product of Kirillov-Reshetikhin
        crystals.

        This is a composition of the map to a tensor product of KR tableaux,
        and then to a tensor product of KR crystals.

        INPUT:

        - ``display_steps`` -- (default: ``False``) boolean which indicates
          if we want to output each step in the algorithm

        EXAMPLES::

            sage: RC = RiggedConfigurations(['D', 4, 1], [[2, 2]])
            sage: elt = RC(partition_list=[[2], [2,2], [1], [1]])
            sage: krc = elt.to_tensor_product_of_kirillov_reshetikhin_crystals(); krc
            [[[2, 3], [3, -2]]]

        We can recover the rigged configuration::

            sage: ret = RC(krc); ret
            <BLANKLINE>
            0[ ][ ]0
            <BLANKLINE>
            -2[ ][ ]-2
            -2[ ][ ]-2
            <BLANKLINE>
            0[ ]0
            <BLANKLINE>
            0[ ]0
            <BLANKLINE>
            sage: elt == ret
            True
        """
        kr_tab = self.to_tensor_product_of_kirillov_reshetikhin_tableaux(display_steps)
        return kr_tab.to_tensor_product_of_kirillov_reshetikhin_crystals()

    def nu(self):
        r"""
        Return the list `\nu` of rigged partitions of this rigged
        configuration element.

        OUTPUT:

        The `\nu` array as a list.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2, 2]])
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
        Action of the crystal operator `e_a` on ``self``.

        This implements the method defined in [CrysStructSchilling06]_ which
        finds the value `k` which is  the length of the string with the
        smallest negative rigging of smallest length. Then it removes a box
        from a string of length `k` in the `a`-th rigged partition, keeping all
        colabels fixed and increasing the new label by one. If no such string
        exists, then `e_a` is undefined.

        .. TODO::

            Implement `f_0` without appealing to tensor product of
            KR tableaux.

        INPUT:

        - ``a`` -- the index of the partition to remove a box

        OUTPUT:

        The resulting rigged configuration element.

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
        if a not in self.parent()._cartan_type.index_set():
            raise ValueError("{} is not in the index set".format(a))
        if a == 0:
            try:
                ret = self.to_tensor_product_of_kirillov_reshetikhin_tableaux().e(0)
                if ret is None:
                    return None
                return ret.to_rigged_configuration()
            except NotImplementedError:
                # We haven't implemented the bijection yet, so return None
                # This is to make sure we can at least view it as a classical
                #   crystal if there is no bijection.
                return None

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
              self.parent()._calc_vacancy_number(ret_RC.nu(), a, rigging_index)
        return(ret_RC)

    def _generate_partition_e(self, a, b, k):
        r"""
        Generate a new partition for a given value of `a` by updating the
        vacancy numbers and preserving co-labels for the map `e_a`.

        INPUT:

        - ``a`` -- the index of the partition we operated on
        - ``b`` -- the index of the partition to generate
        - ``k`` -- the length of the string with the smallest negative
          rigging of smallest length

        OUTPUT:

        The constructed rigged partition.

        TESTS::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2,1]])
            sage: RC(partition_list=[[1], [1], [1], [1]])._generate_partition_e(1, 2, 1)
            -1[ ]-1
            <BLANKLINE>
        """
        # Check to make sure we will do something
        if self.parent()._cartan_matrix[a][b] == 0:
            return self[b]

        new_list = self[b][:]
        new_vac_nums = self[b].vacancy_numbers[:]
        new_rigging = self[b].rigging[:]

        # Update the vacancy numbers and the rigging
        value = self.parent()._cartan_matrix[a][b]
        for i in range(len(new_vac_nums)):
            if new_list[i] < k:
                break

            new_vac_nums[i] += value
            new_rigging[i] += value

        return(RiggedPartition(new_list, new_rigging, new_vac_nums))

    def f(self, a):
        r"""
        Action of the crystal operator `f_a` on ``self``.

        This implements the method defined in [CrysStructSchilling06]_ which
        finds the value `k` which is  the length of the string with the
        smallest nonpositive rigging of largest length. Then it adds a box from
        a string of length `k` in the `a`-th rigged partition, keeping all
        colabels fixed and decreasing the new label by one. If no such string
        exists, then it adds a new string of length 1 with label `-1`. If any
        of the resulting vacancy numbers are larger than the labels (i.e. it
        is an invalid rigged configuration), then `f_a` is undefined.

        .. TODO::

            Implement `f_0` without appealing to tensor product of
            KR tableaux.

        INPUT:

        - ``a`` -- the index of the partition to add a box

        OUTPUT:

        The resulting rigged configuration element.

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
        if a not in self.parent()._cartan_type.index_set():
            raise ValueError("{} is not in the index set".format(a))
        if a == 0:
            try:
                ret = self.to_tensor_product_of_kirillov_reshetikhin_tableaux().f(0)
                if ret is None:
                    return None
                return ret.to_rigged_configuration()
            except NotImplementedError:
                # We haven't implemented the bijection yet, so return None
                # This is to make sure we can at least view it as a classical
                #   crystal if there is no bijection.
                return None

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
        return self.__class__(self.parent(), new_partitions)

    def _generate_partition_f(self, a, b, k):
        r"""
        Generate a new partition for a given value of `a` by updating the
        vacancy numbers and preserving co-labels for the map `f_a`.

        INPUT:

        - ``a`` -- the index of the partition we operated on
        - ``b`` -- the index of the partition to generate
        - ``k`` -- the length of the string with smallest nonpositive rigging
          of largest length

        OUTPUT:

        The constructed rigged partition.

        TESTS::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2,1]])
            sage: RC(partition_list=[[1], [1], [1], [1]])._generate_partition_f(1, 2, 1)
            0[ ]0
            <BLANKLINE>
        """
        # Check to make sure we will do something
        if self.parent()._cartan_matrix[a][b] == 0:
            return self[b]

        new_list = self[b][:]
        new_vac_nums = self[b].vacancy_numbers[:]
        new_rigging = self[b].rigging[:]

        # Update the vacancy numbers and the rigging
        value = self.parent()._cartan_matrix[a][b]
        for i in range(len(new_vac_nums)):
            if new_list[i] <= k:
                break

            new_vac_nums[i] -= value
            new_rigging[i] -= value

        return(RiggedPartition(new_list, new_rigging, new_vac_nums))

#    def epsilon(self, a):
#        r"""
#        Return `\varepsilon_a` of ``self``.
#
#        Let `x_{\ell}` be the smallest string of `\nu^{(a)}` or `0` if
#        `\nu^{(a)} = \emptyset`, then we have
#        `\varepsilon_a = -\min(0, x_{\ell})`.
#
#        EXAMPLES::
#
#            sage: RC = RiggedConfigurations(['A', 4, 1], [[2,1]])
#        """
#        if len(self[a-1]) == 0:
#            return 0
#        return -min(0, min(self[a-1].rigging))
#
#    def phi(self, a):
#        r"""
#        Return `\varphi_a` of ``self``.
#
#        Let `x_{\ell}` be the smallest string of `\nu^{(a)}` or `0` if
#        `\nu^{(a)} = \emptyset`, then we have
#        `\varepsilon_a = p_{\infty}^{(a)} - min(0, x_{\ell})`.
#
#        EXAMPLES::
#
#            sage: RC = RiggedConfigurations(['A', 4, 1], [[2,1]])
#        """
#        p_inf = self.parent()._calc_vacancy_number(self, a, None)
#        if len(self[a-1]) == 0:
#            return p_inf
#        return p_inf - min(0, min(self[a-1].rigging))

    @cached_method
    def cocharge(self):
        r"""
        Compute the cocharge statistic of ``self``.

        Computes the cocharge statistic [CrysStructSchilling06]_ on this
        rigged configuration `(\nu, J)`. The cocharge statistic is defined as:

        .. MATH::

            cc(\nu, J) = \frac{1}{2} \sum_{a, b \in J}
            \sum_{j,k \geq 1} \left( \alpha_a \mid \alpha_b \right)
            \min(j, k) m_j^{(a)} m_k^{(b)}
            + \sum_{a, i} \left\lvert J^{(a, i)} \right\rvert.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A', 3, 1], [[3, 2], [2,1], [1,1]])
            sage: RC(partition_list=[[1], [1], []]).cocharge()
            1
        """
        cc = 0
        rigging_sum = 0
        for a, p in enumerate(self):
            for pos, i in enumerate(p._list):
                # Add the rigging
                rigging_sum += p.rigging[pos]
                # Add the L matrix contribution
                for dim in self.parent().dims:
                    if dim[0] == a + 1:
                        cc += min(dim[1], i)
                # Subtract the vacancy number
                cc -= p.vacancy_numbers[pos]
        return cc / 2 + rigging_sum

    cc = cocharge

    @cached_method
    def charge(self):
        r"""
        Compute the charge statistic of ``self``.

        Let `B` denote a set of rigged configurations. The *charge* `c` of
        a rigged configuration `b` is computed as

        .. MATH::

            c(b) = \max(cc(b) \mid b \in B) - cc(b).

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A', 3, 1], [[3, 2], [2,1], [1,1]])
            sage: RC(partition_list=[[],[],[]]).charge()
            2
            sage: RC(partition_list=[[1], [1], []]).charge()
            1
        """
        B = self.parent()
        if not hasattr(B, "_max_charge"):
            B._max_charge = max(b.cocharge() for b in B)
        return B._max_charge - self.cocharge()

    @cached_method
    def classical_weight(self):
        r"""
        Return the classical weight of ``self``.

        The classical weight `\Lambda` of a rigged configuration is

        .. MATH::

            \Lambda = \sum_{a \in \overline{I}} \sum_{i > 0}
            i L_i^{(a)} \Lambda_a - \sum_{a \in \overline{I}} \sum_{i > 0}
            i m_i^{(a)} \alpha_a

        EXAMPLES::

            sage: RC = RiggedConfigurations(['D',4,1], [[2,2]])
            sage: elt = RC(partition_list=[[2],[2,1],[1],[1]])
            sage: elt.classical_weight()
            (0, 1, 1, 0)

        This agrees with the corresponding classical weight as KR tableaux::

            sage: krt = elt.to_tensor_product_of_kirillov_reshetikhin_tableaux(); krt
            [[2, 1], [3, -1]]
            sage: krt.classical_weight() == elt.classical_weight()
            True

        TESTS:

        We check the classical weights agree in an entire crystal::
        
            sage: RC = RiggedConfigurations(['A',2,1], [[2,1], [1,1]])
            sage: passed_test = True
            sage: for x in RC:
            ....:    y = x.to_tensor_product_of_kirillov_reshetikhin_tableaux()
            ....:    if x.classical_weight() != y.classical_weight():
            ....:        passed_test = False
            ....:        break
            sage: passed_test
            True
        """
        F = self.cartan_type().classical().root_system()
        if F.ambient_space() is None:
            WLR = F.weight_lattice()
        else:
            WLR = F.ambient_space()
        la = WLR.fundamental_weights()
        sum = WLR.zero()
        for dim in self.parent().dims:
            sum += dim[1] * la[dim[0]]

        alpha = WLR.simple_roots()
        for a, partition in enumerate(self):
            for row_len in partition:
                sum -= row_len * alpha[a+1] # +1 for indexing
        return sum

    def get_vacancy_numbers(self, a):
        r"""
        Return the list of all vacancy numbers of the rigged partition
        `\nu^{(a)}` (with duplicates).

        INPUT:

        - ``a`` -- the index of the rigged partition

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

        - ``a`` -- the index of the rigged partition

        - ``i`` -- the row of the rigged partition

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

    def partition_rigging_lists(self):
        """
        Return the list of partitions and the associated list of riggings
        of ``self``.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A',3,1], [[1,2],[2,2]])
            sage: rc = RC(partition_list=[[2],[1],[1]], rigging_list=[[-1],[0],[-1]]); rc
            <BLANKLINE>
            -1[ ][ ]-1
            <BLANKLINE>
            1[ ]0
            <BLANKLINE>
            -1[ ]-1
            <BLANKLINE>
            sage: rc.partition_rigging_lists()
            [[[2], [1], [1]], [[-1], [0], [-1]]]
        """
        partitions = []
        riggings = []
        for p in self:
            partitions.append(list(p))
            riggings.append(list(p.rigging))
        return [partitions, riggings]

class RCNonSimplyLacedElement(RiggedConfigurationElement):
    """
    Rigged configuration elements for non-simply-laced types.

    TESTS::

        sage: RC = RiggedConfigurations(['C',2,1], [[1,2],[1,1],[2,1]])
        sage: elt = RC(partition_list=[[3],[2]]); elt
        <BLANKLINE>
        0[ ][ ][ ]0
        <BLANKLINE>
        0[ ][ ]0
        sage: TestSuite(elt).run()
    """
    def to_virtual_configuration(self):
        """
        Return the corresponding rigged configuration in the virtual crystal.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['C',2,1], [[1,2],[1,1],[2,1]])
            sage: elt = RC(partition_list=[[3],[2]]); elt
            <BLANKLINE>
            0[ ][ ][ ]0
            <BLANKLINE>
            0[ ][ ]0
            sage: elt.to_virtual_configuration()
            <BLANKLINE>
            0[ ][ ][ ]0
            <BLANKLINE>
            0[ ][ ][ ][ ]0
            <BLANKLINE>
            0[ ][ ][ ]0
        """
        return self.parent().to_virtual(self)

    def e(self, a):
        """
        Return the action of `e_a` on ``self``.

        This works by lifting into the virtual configuration, then applying

        .. MATH::

            \hat{e}_a = \prod_{j \in \iota(a)} e_j^{\gamma_j}

        and pulling back.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A',6,2], [[1,1]]*7)
            sage: elt = RC(partition_list=[[1]*5,[2,1,1],[3,2]])
            sage: elt.e(3)
            <BLANKLINE>
            0[ ]0
            0[ ]0
            0[ ]0
            0[ ]0
            0[ ]0
            <BLANKLINE>
            0[ ][ ]0
            1[ ]1
            1[ ]1
            <BLANKLINE>
            1[ ][ ]1
            1[ ]0
            <BLANKLINE>
        """
        vct = self.parent()._folded_ct
        L = []
        gamma = vct.scaling_factors()
        for i in vct.folding_orbit()[a]:
            L.extend([i]*gamma[a])
        virtual_rc = self.parent().to_virtual(self).e_string(L)
        if virtual_rc is None:
            return None
        return self.parent().from_virtual(virtual_rc)

    def f(self, a):
        """
        Return the action of `f_a` on ``self``.

        This works by lifting into the virtual configuration, then applying

        .. MATH::

            \hat{f}_a = \prod_{j \in \iota(a)} f_j^{\gamma_j}

        and pulling back.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A',6,2], [[1,1]]*7)
            sage: elt = RC(partition_list=[[1]*5,[2,1,1],[2,1]], rigging_list=[[0]*5,[0,1,1],[1,0]])
            sage: elt.f(3)
            <BLANKLINE>
            0[ ]0
            0[ ]0
            0[ ]0
            0[ ]0
            0[ ]0
            <BLANKLINE>
            1[ ][ ]1
            1[ ]1
            1[ ]1
            <BLANKLINE>
            -1[ ][ ][ ]-1
            0[ ][ ]0
            <BLANKLINE>
        """
        vct = self.parent()._folded_ct
        L = []
        gamma = vct.scaling_factors()
        for i in vct.folding_orbit()[a]:
            L.extend([i]*gamma[a])
        virtual_rc = self.parent().to_virtual(self).f_string(L)
        if virtual_rc is None:
            return None
        return self.parent().from_virtual(virtual_rc)

    @cached_method
    def cocharge(self):
        r"""
        Compute the cocharge statistic.

        Computes the cocharge statistic [OSS03]_ on this
        rigged configuration `(\nu, J)` by computing the cocharge as a virtual
        rigged configuration `(\hat{\nu}, \hat{J})` and then using the
        identity `cc(\hat{\nu}, \hat{J}) = \gamma_0 cc(\nu, J)`.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['C', 3, 1], [[2,1], [1,1]])
            sage: RC(partition_list=[[1,1],[2,1],[1,1]]).cocharge()
            1
        """
        #return self.to_virtual_configuration().cocharge() / self.parent()._folded_ct.gamma[0]
        vct = self.parent()._folded_ct
        cc = 0
        rigging_sum = 0
        sigma = vct.folding_orbit()
        gamma = vct.scaling_factors()
        for a, p in enumerate(self):
            t_check = len(sigma[a+1]) * gamma[a+1] / gamma[0]
            for pos, i in enumerate(p._list):
                # Add the rigging
                rigging_sum += t_check * p.rigging[pos]
                # Add the L matrix contribution
                for dim in self.parent().dims:
                    if dim[0] == a + 1:
                        cc += t_check * min(dim[1], i)
                # Subtract the vacancy number
                cc -= t_check * p.vacancy_numbers[pos]
        return cc / 2 + rigging_sum

    cc = cocharge

