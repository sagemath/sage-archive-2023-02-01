r"""
Rigged Configurations of `\mathcal{B}(\infty)`

AUTHORS:

- Travis Scrimshaw (2013-04-16): Initial version
"""

#*****************************************************************************
#       Copyright (C) 2013 Travis Scrimshaw <tscrim@ucdavis.edu>
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
from sage.structure.list_clone import ClonableArray
from sage.rings.all import QQ
from sage.categories.highest_weight_crystals import HighestWeightCrystals
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.rigged_configurations.rigged_partition import RiggedPartition

# Note on implementation, this class is used for simply-laced types only
class InfinityCrystalOfRiggedConfigurations(Parent, UniqueRepresentation):
    r"""
    Class of rigged configurations modeling `\mathcal{B}(\infty)`.

    INPUT:

    - ``cartan_type`` -- a Cartan type
    """
    @staticmethod
    def __classcall_private__(cls, cartan_type):
        r"""
        Normalize the input arguments to ensure unique representation.

        EXAMPLES::

            sage: RC1 = InfinityCrystalOfRiggedConfigurations(CartanType(['A',3]))
            sage: RC2 = InfinityCrystalOfRiggedConfigurations(['A',3])
            sage: RC2 is RC1
            True
        """
        cartan_type = CartanType(cartan_type)
        if not cartan_type.is_simply_laced():
            vct = cartan_type.as_folding()
            return InfinityCrystalOfVirtualRC(vct)

        return super(InfinityCrystalOfRiggedConfigurations, cls).__classcall__(cls, cartan_type)

    def __init__(self, cartan_type):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: RC = InfinityCrystalOfRiggedConfigurations(['A',3])
            sage: TestSuite(RC).run()
        """
        self._cartan_type = cartan_type
        # We store the cartan matrix for the vacancy number calculations for speed
        self._cartan_matrix = self._cartan_type.cartan_matrix()
        self.module_generators = (self.element_class(self, rigging_list=[[]]*cartan_type.rank()),)
        Parent.__init__(self, category=HighestWeightCrystals())

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: InfinityCrystalOfRiggedConfigurations(['A',3])
            The infinity crystal of rigged configurations of type ['A', 3]
        """
        return "The infinity crystal of rigged configurations of type %s"%self._cartan_type

    def __iter__(self):
        """
        Returns the iterator of ``self``.

        EXAMPLES::

            sage: RC = InfinityCrystalOfRiggedConfigurations(['A', 3])
            sage: g = RC.__iter__()
            sage: g.next()
            <BLANKLINE>
            (/)
            <BLANKLINE>
            (/)
            <BLANKLINE>
            (/)
            <BLANKLINE>
            sage: g.next()
            <BLANKLINE>
            (/)
            <BLANKLINE>
            -1[ ]-1
            <BLANKLINE>
            (/)
            <BLANKLINE>
        """
        index_set = self._cartan_type.index_set()
        from sage.combinat.backtrack import TransitiveIdeal
        return TransitiveIdeal(lambda x: [x.f(i) for i in index_set],
                               self.module_generators).__iter__()

    def _element_constructor_(self, lst=None, **options):
        """
        Construct an element of ``self`` from ``lst``.
        """
        return self.element_class(self, lst, **options)

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

            sage: RC = InfinityCrystalOfRiggedConfigurations(['A', 4, 1])
            sage: elt = RC(partition_list=[[1], [1], [], []])
            sage: RC._calc_vacancy_number(elt.nu(), 1, 0)
            0
        """
        row_len = partitions[a][i]
        vac_num = 0
        for b, value in enumerate(self._cartan_matrix.row(a)):
            vac_num -= value * partitions[b].get_num_cells_to_column(row_len)

        return vac_num

    class Element(ClonableArray):
        """
        A rigged configuration for simply-laced types.

        Typically to create a specific rigged configuration, the user will pass
        in the optional argument **partition_list** and if the user wants to
        specify the rigging values, give the optional argument **rigging_list**
        as well. If **rigging_list** is not passed, the rigging values are set
        to the corresponding vacancy numbers.

        INPUT:

        - ``parent``            -- The parent of this element

        - ``rigged_partitions`` -- A list of rigged partitions

        There are two optional arguments to explicitly construct a rigged
        configuration. The first is **partition_list** which gives a list of
        partitions, and the second is **rigging_list** which is a list of
        corresponding lists of riggings. If only partition_list is specified,
        then it sets the rigging equal to the calculated vacancy numbers.

        EXAMPLES:

        Type `A_n^{(1)}` examples::

            sage: RC = InfinityCrystalOfRiggedConfigurations(['A', 4, 1])
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

            sage: RC = InfinityCrystalOfRiggedConfigurations(['A', 4, 1])
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

            sage: RC = InfinityCrystalOfRiggedConfigurations(['D', 4, 1])
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

            sage: RC = InfinityCrystalOfRiggedConfigurations(['D', 4, 1])
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
        def __init__(self, parent, rigged_partitions=[], **options):
            r"""
            Construct a rigged configuration element.

            EXAMPLES::

                sage: RC = InfinityCrystalOfRiggedConfigurations(['A', 4, 1])
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

                sage: RC = InfinityCrystalOfRiggedConfigurations(['D', 4, 1])
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

                sage: RC = InfinityCrystalOfRiggedConfigurations(['D', 4, 1])
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

        def check(self):
            """
            Nothing to check.

            TESTS::

                sage: RC = RiggedConfigurations(['A', 4, 1], [[2, 1]])
                sage: elt = RC(partition_list=[[1], [1], [], []])
                sage: elt.check()
            """
            pass # No checks

        def nu(self):
            r"""
            Return the list `\nu` of rigged partitions of this rigged
            configuration element.

            OUTPUT:

            - The `\nu` array as a list

            EXAMPLES::

                sage: RC = InfinityCrystalOfRiggedConfigurations(['A', 4, 1])
                sage: RC(partition_list=[[2], [2,2], [2], [2]]).nu()
                [0[ ][ ]0
                , -2[ ][ ]-2
                -2[ ][ ]-2
                , 2[ ][ ]2
                , -2[ ][ ]-2
                ]
            """
            return list(self)

        def epsilon(self, a):
            r"""
            Return `\varepsilon_a` of ``self``.
            """
            ep = 0
            cur = self.e(a)
            while cur is not None:
                cur = cur.e(a)
                ep += 1
            return ep

        def phi(self, a):
            r"""
            Return `\varphi_a` of ``self``.

            Let `x \in \mathcal{RC}(\infty)` Define `\varphi_a(x) :=
            \varepsilon_a(x) + \langle h_a, \mathrm{wt}(x) \rangle`,
            where `h_a` is the `a`-th simple coroot and `\mathrm{wt}(x)` is
            the :meth:`weight` of `x`.

            INPUT:

            - ``i`` -- An element of the index set

            EXAMPLES::
            """
            P = self.parent().weight_lattice_realization()
            h = P.simple_coroots()
            return self.epsilon(a) + P(self.weight()).scalar(h[a])

        def weight(self):
            """
            Return the weight of ``self``.
            """
            P = self.parent().weight_lattice_realization()
            alpha = list(P.simple_roots())
            return sum(sum(x) * alpha[i] for i,x in enumerate(self))

        def e(self, a):
            r"""
            Action of the crystal operator `e_a` on this rigged configuration
            element.

            This implements the method defined in [CrysStructSchilling06]_ which
            finds the value `k` which is  the length of the string with the
            smallest negative rigging of smallest length. Then it removes a box
            from a string of length `k` in the `a`-th rigged partition, keeping all
            colabels fixed and increasing the new label by one. If no such string
            exists, then `e_a` is undefined.

            INPUT:

            - ``a`` -- The index of the partition to remove a box.

            OUTPUT:

            The resulting rigged configuration element.

            EXAMPLES::

                sage: RC = InfinityCrystalOfRiggedConfigurations(['A', 4, 1])
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
            a = self.index_set().index(a)

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

            - ``a`` -- The index of the partition we operated on.
            - ``b`` -- The index of the partition to generate.
            - ``k`` -- The length of the string with the smallest negative rigging of smallest length.

            OUTPUT:

            - The constructed rigged partition.

            TESTS::

                sage: RC = InfinityCrystalOfRiggedConfigurations(['A', 4, 1])
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

            The resulting rigged configuration element.

            EXAMPLES::

                sage: RC = InfinityCrystalOfRiggedConfigurations(['A', 4, 1])
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
            a = self.index_set().index(a)

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

            # Note that we do not need to sort the rigging since if there was a
            #   smaller rigging in a larger row, then `k` would be larger.
            return self.__class__(self.parent(), new_partitions)

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

                sage: RC = InfinityCrystalOfRiggedConfigurations(['A', 4, 1])
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

        @cached_method
        def cocharge(self):
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
            """
            cc = 0
            rigging_sum = 0
            for a, p in enumerate(self):
                for pos, i in enumerate(p._list):
                    # Add the rigging
                    rigging_sum += p.rigging[pos]
                    # Subtract the vacancy number
                    cc -= p.vacancy_numbers[pos]
            return cc / 2 + rigging_sum

        cc = cocharge

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
            """
            F = self.parent()._cartan_type.root_system()
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

class InfinityCrystalOfVirtualRC(InfinityCrystalOfRiggedConfigurations):
    r"""
    Rrigged configurations for `\mathcal{B}(\infty)` in non-simply-laced types.
    """
    def __init__(self, vct):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: RC = InfinityCrystalOfRiggedConfigurations(['C',2,1]); RC
            Rigged configurations of type ['C', 2, 1]
         """
        self._folded_ct = vct
        InfinityCrystalOfRiggedConfigurations.__init__(self, vct._cartan_type)

    @lazy_attribute
    def virtual(self):
        """
        Return the corresponding virtual crystal.

        EXAMPLES::

            sage: RC = InfinityCrystalOfRiggedConfigurations(['C',2])
            sage: RC
            B infinity rigged configurations of type ['C', 3]
            sage: RC.virtual
            B infinity rigged configurations of type ['A', 3]
        """
        return InfinityCrystalOfRiggedConfigurations(self._folded_ct._folding)

    def to_virtual(self, rc):
        """
        Convert ``rc`` into a rigged configuration in the virtual crystal.

        INPUT:

        - ``rc`` -- A rigged configuration element

        EXAMPLES::

            sage: RC = InfinityCrystalOfRiggedConfigurations(['C',2])
            sage: elt = RC(partition_list=[[3],[2]]); elt
            <BLANKLINE>
            0[ ][ ][ ]0
            <BLANKLINE>
            0[ ][ ]0
            sage: velt = RC.to_virtual(elt); velt
            <BLANKLINE>
            0[ ][ ][ ]0
            <BLANKLINE>
            0[ ][ ][ ][ ]0
            <BLANKLINE>
            0[ ][ ][ ]0
            sage: velt.parent()
            B infinity rigged configurations of type ['A', 3]
        """
        gamma = map(int, self._folded_ct.scaling_factors())
        sigma = self._folded_ct._orbit
        n = self._folded_ct._folding.rank()
        vindex = self._folded_ct._folding.index_set()
        partitions = [None] * n
        riggings = [None] * n
        vac_nums = [None] * n
        # -1 for indexing
        for a, rp in enumerate(rc):
            for i in sigma[a]:
                k = vindex.index(i)
                partitions[k] = [row_len*gamma[a] for row_len in rp._list]
                riggings[k] = [rig_val*gamma[a] for rig_val in rp.rigging]
                vac_nums[k] = [vac_num*gamma[a] for vac_num in rp.vacancy_numbers]
        return self.virtual.element_class(self.virtual, partition_list=partitions,
                            rigging_list=riggings,
                            vacancy_numbers_list=vac_nums)

    def from_virtual(self, vrc):
        """
        Convert ``vrc`` in the virtual crystal into a rigged configution of
        the original Cartan type.

        INPUT:

        - ``vrc`` -- A virtual rigged configuration

        EXAMPLES::

            sage: RC = InfinityCrystalOfRiggedConfigurations(['C',2])
            sage: elt = RC(partition_list=[[3],[2]])
            sage: vrc_elt = RC.to_virtual(elt)
            sage: ret = RC.from_virtual(vrc_elt); ret
            <BLANKLINE>
            0[ ][ ][ ]0
            <BLANKLINE>
            0[ ][ ]0
            sage: ret == elt
            True
        """
        gamma = list(self._folded_ct.scaling_factors()) #map(int, self._folded_ct.scaling_factors())
        sigma = self._folded_ct._orbit
        n = self._cartan_type.rank()
        partitions = [None] * n
        riggings = [None] * n
        vac_nums = [None] * n
        vindex = self._folded_ct._folding.index_set()
        # TODO: Handle special cases for A^{(2)} even and its dual?
        for a in range(n):
            index = vindex.index(sigma[a][0])
            partitions[a] = [row_len // gamma[a] for row_len in vrc[index]._list]
            riggings[a] = [rig_val / gamma[a] for rig_val in vrc[index].rigging]
            vac_nums[a] = [vac_val / gamma[a] for vac_val in vrc[index].vacancy_numbers]
        return self.element_class(self, partition_list=partitions,
                                  rigging_list=riggings, vacancy_numbers_list=vac_nums)

    class Element(InfinityCrystalOfRiggedConfigurations.Element):
        """
        Virtual rigged configuration elements.

        TESTS::
        """
        def to_virtual_configuration(self):
            """
            Return the corresponding rigged configuration in the virtual crystal.

            EXAMPLES::
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
            """
            vct = self.parent()._folded_ct
            gamma = vct.scaling_factors()
            L = []
            for i in vct.folding_orbit()[a]:
                L.extend([i]*int(gamma[a]))
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
            """
            vct = self.parent()._folded_ct
            gamma = vct.scaling_factors()
            L = []
            for i in vct.folding_orbit()[a]:
                L.extend([i]*int(gamma[a]))
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
            """
            #return self.to_virtual_configuration().cocharge() / self.parent()._folded_ct.gamma[0]
            vct = self.parent()._folded_ct
            gamma = vct.scaling_factors()
            sigma = vct.folding_orbit()
            cc = 0
            i_set = vct.cartan_type().index_set()
            rigging_sum = 0
            for a, p in enumerate(self):
                a = i_set[a]
                t_check = len(sigma[a]) * gamma[a] # / gamma[0]
                for pos, i in enumerate(p._list):
                    # Add the rigging
                    rigging_sum += t_check * p.rigging[pos]
                    # Subtract the vacancy number
                    cc -= t_check * p.vacancy_numbers[pos]
            return cc / 2 + rigging_sum

        cc = cocharge

class InfinityCrystalOfRCA2Even(InfinityCrystalOfVirtualRC):
    """
    Infinity crystal of rigged configurations for type `A_{2n}^{(2)}`.
    """
    def to_virtual(self, rc):
        """
        Convert ``rc`` into a rigged configuration in the virtual crystal.

        INPUT:

        - ``rc`` -- a rigged configuration element

        EXAMPLES::

            sage: from sage.combinat.rigged_configurations.rc_infinity import InfinityCrystalOfRCA2Even
            sage: RC = InfinityCrystalOfRCA2Even(CartanType(['A',4,2]).as_folding())
            sage: elt = RC(partition_list=[[1],[1]]); elt
            <BLANKLINE>
            -1[ ]-1
            <BLANKLINE>
            1[ ]1
            <BLANKLINE>
            sage: velt = RC.to_virtual(elt); velt
            <BLANKLINE>
            -1[ ]-1
            <BLANKLINE>
            2[ ]2
            <BLANKLINE>
            -1[ ]-1
            <BLANKLINE>
            sage: velt.parent()
            Rigged configurations of type ['A', 3, 1] and factor(s) ((2, 2), (2, 2))
        """
        gamma = self._folded_ct.scaling_factors()
        sigma = self._folded_ct.folding_orbit()
        n = self._folded_ct._folding.rank()
        partitions = [None] * n
        riggings = [None] * n
        vac_nums = [None] * n
        # +/- 1 for indexing
        for a in range(len(rc)):
            for i in sigma[a]:
                partitions[i] = [row_len for row_len in rc[a]._list]
                riggings[i] = [rig_val*gamma[a] for rig_val in rc[a].rigging]
                vac_nums[i] = [vac_num*gamma[a] for vac_num in rc[a].vacancy_numbers]
        return self.virtual.element_class(self.virtual, partition_list=partitions,
                            rigging_list=riggings,
                            vacancy_numbers_list=vac_nums)

    def from_virtual(self, vrc):
        """
        Convert ``vrc`` in the virtual crystal into a rigged configution of
        the original Cartan type.

        INPUT:

        - ``vrc`` -- a virtual rigged configuration element

        EXAMPLES::

            sage: from sage.combinat.rigged_configurations.rc_infinity import InfinityCrystalOfRCA2Even
            sage: RC = InfinityCrystalOfRCA2Even(CartanType(['A',4,2]).as_folding())
            sage: elt = RC(partition_list=[[1],[1]])
            sage: velt = RC.to_virtual(elt)
            sage: ret = RC.from_virtual(velt); ret
            <BLANKLINE>
            -1[ ]-1
            <BLANKLINE>
            1[ ]1
            <BLANKLINE>
            sage: ret == elt
            True
        """
        gamma = self._folded_ct.scaling_factors()
        sigma = self._folded_ct.folding_orbit()
        n = self._cartan_type.rank()
        partitions = [None] * n
        riggings = [None] * n
        vac_nums = [None] * n
        # +/- 1 for indexing
        for a in range(n):
            index = sigma[a][0]
            partitions[index] = [row_len for row_len in vrc[index]._list]
            riggings[index] = [rig_val//gamma[a] for rig_val in vrc[index].rigging]
            vac_nums[a] = [vac_val//gamma[a] for vac_val in vrc[index].vacancy_numbers]
        return self.element_class(self, partition_list=partitions,
                                  rigging_list=riggings, vacancy_numbers_list=vac_nums)

