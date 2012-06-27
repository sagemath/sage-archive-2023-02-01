r"""
Rigged partition

Class and methods of the rigged partition which are used by the rigged
configuration class. This is an internal class used by the rigged
configurations and KR tableaux during the bijection, and is not to be used by
the end-user.

We hold the partitions as an 1-dim array of positive integers where each
value corresponds to the length of the row. This is the shape of the
partition which can be accessed by the regular index.

The data for the vacancy number is also stored in a 1-dim array which each
entry corresponds to the row of the tableau, and similarly for the
partition values.

AUTHORS:

- Travis Scrimshaw (2010-09-26): Initial version

.. TODO:: Convert this to using `m_i` (multiplicities) with a dictionary?
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

from sage.combinat.partition import Partition_class

class RiggedPartition(Partition_class):
    r"""
    The RiggedPartition class which is the data structure of a rigged (i.e.
    marked or decorated) Young diagram of a partition.

    Note that this class as a stand-alone object does not make sense since the
    vacancy numbers are calculated using the entire rigged configuration. For
    more, see :class:`RiggedConfigurations`.

    EXAMPLES::

        sage: RC = RiggedConfigurations(['A', 4, 1], [[2, 2]])
        sage: RP = RC(partition_list=[[2],[2,2],[2,1],[2]])[2]
        sage: RP
        0[ ][ ]0
        -1[ ]-1
        <BLANKLINE>
        sage: type(RP)
        <class 'sage.combinat.rigged_configurations.rigged_partition.RiggedPartition'>
    """

    def __init__(self, shape=None, rigging_list=None, vacancy_nums=None):
        r"""
        Initialize by the rigged partition.

        Note that this only performs checks to see that the sizes match up.

        INPUT:

        - ``shape`` -- (Default: None) The shape
        - ``rigging_list`` -- (Default: None) The riggings
        - ``vacancy_nums`` -- (Default: None) The vacancy numbers

        TESTS::

            sage: sage.combinat.rigged_configurations.rigged_partition.RiggedPartition()
            (/)
            <BLANKLINE>
            sage: RP = sage.combinat.rigged_configurations.rigged_partition.RiggedPartition([2,1], [0,0], [1, 0])
            sage: RP
            1[ ][ ]0
            0[ ]0
            <BLANKLINE>
            sage: TestSuite(RP).run()
        """
        if shape is None:
            Partition_class.__init__(self, [])
            self.vacancy_numbers = []
            self.rigging = []
        else:
            Partition_class.__init__(self, shape)

            if vacancy_nums is not None:
                if len(shape) != len(vacancy_nums):
                    raise ValueError

                self.vacancy_numbers = list(vacancy_nums)
            else:
                self.vacancy_numbers = [None] * len(shape)

            if rigging_list is not None:

                if len(shape) != len(rigging_list):
                    raise ValueError

                self.rigging = list(rigging_list)
            else:
                self.rigging = [None] * len(shape)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2, 2]])
            sage: RC(partition_list=[[2],[2,2],[2,1],[2]])[2] # indirect doctest
            0[ ][ ]0
            -1[ ]-1
            <BLANKLINE>
        """
        # If it is empty, return saying so
        if len(self._list) == 0:
            return("(/)\n")

        retStr = ""
        for i, val in enumerate(self._list):
            retStr += str(self.vacancy_numbers[i])
            retStr += "[ ]"*val
            retStr += str(self.rigging[i])
            retStr += "\n"
        return(retStr)

    def _latex_(self):
        """
        Returns LaTeX representation of ``self``.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2, 2]])
            sage: latex(RC(partition_list=[[2],[2,2],[2,1],[2]])[2]) # indirect doctest
            {
            \begin{array}[t]{r|c|c|l}
            \cline{2-3} 0 &\phantom{x}&\phantom{x}& 0 \\
            \cline{2-3} -1 &\phantom{x}& \multicolumn{2}{l}{-1} \\
            \cline{2-2}
            \end{array}
            }

        """
        num_rows = len(self._list)
        if num_rows == 0:
            return "{\\emptyset}"

        num_cols = self._list[0]
        ret_string = ("{\n\\begin{array}[t]{r|" + "c|"*num_cols + "l}\n"
                       + "\\cline{2-" + repr(1 + num_cols) + "} "
                       + repr(self.vacancy_numbers[0]))
        for i, row_len in enumerate(self._list):

            ret_string += " &" + "\\phantom{x}&"*row_len

            if num_cols == row_len:
                ret_string += " " + repr(self.rigging[i])
            else:
                ret_string += " \\multicolumn{" + repr(num_cols - row_len + 1)
                ret_string += "}{l}{" + repr(self.rigging[i]) + "}"

            ret_string += " \\\\\n"

            ret_string += "\\cline{2-" + repr(1 + row_len) + "} "
            if i != num_rows - 1 and row_len != self._list[i + 1]:
                ret_string += repr(self.vacancy_numbers[i + 1])
        ret_string += "\n\\end{array}\n}"

        return ret_string

    def _clone(self):
        r"""
        Makes a (deep) copy of this rigged partition.

        TESTS::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2, 2]])
            sage: RP = RC(partition_list=[[2],[2,2],[2,1],[2]])[2]; RP
            0[ ][ ]0
            -1[ ]-1
            <BLANKLINE>
            sage: RP2 = RP._clone(); RP2
            0[ ][ ]0
            -1[ ]-1
            <BLANKLINE>
            sage: RP == RP2
            True
            sage: RP is RP2
            False
        """
        return(RiggedPartition(self._list[:], self.rigging[:], self.vacancy_numbers[:]))

    def __eq__(self, rhs):
        r"""
        Return true if ``self`` equals ``rhs``.

        TESTS::

            sage: RC = RiggedConfigurations(['A',2,1], [[1,1],[1,1],[1,1]])
            sage: x = RC(partition_list=[[1], []], rigging_list=[[0], []])
            sage: y = RC(partition_list=[[1], []], rigging_list=[[1], []])
            sage: x == y
            False
        """
        if isinstance(rhs, RiggedPartition):
            return(self._list.__eq__(rhs._list) and self.rigging.__eq__(rhs.rigging))

        return False

    # Should we move these functions to the CP -> RC bijections?

    def get_num_cells_to_column(self, end_column):
        r"""
        Get the number of cells in all columns before the ``end_column``.

        INPUT:

        - ``end_column`` -- The index of the column to end at

        OUTPUT:

        - The number of cells

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2, 2]])
            sage: RP = RC(partition_list=[[2],[2,2],[2,1],[2]])[2]
            sage: RP.get_num_cells_to_column(1)
            2
            sage: RP.get_num_cells_to_column(2)
            3
            sage: RP.get_num_cells_to_column(3)
            3
        """
        sum_cells = 0
        # Sum up from the reverse (the smallest row sizes)
        i = len(self._list) - 1
        while i >= 0 and self._list[i] < end_column:
            sum_cells += self._list[i]
            i -= 1

        # Add the remaining cells
        sum_cells += end_column * (i + 1)

        return sum_cells

    def insert_cell(self, max_width):
        r"""
        Insert a cell given at a singular value as long as its less than the
        specified width.

        Note that :meth:`insert_cell` does not update riggings or vacancy numbers, but
        it does prepare the space for them. Returns the width of the row we
        inserted at.

        INPUT:

        - ``max_width`` -- The maximum width (i.e. row length) that we can insert the cell at

        OUTPUT:

        - The width of the row we inserted at.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2, 2]])
            sage: RP = RC(partition_list=[[2],[2,2],[2,1],[2]])[2]
            sage: RP.insert_cell(2)
            2
            sage: RP
            0[ ][ ][ ]None
            -1[ ]-1
            <BLANKLINE>
        """
        maxPos = -1
        for i, vacNum in enumerate(self.vacancy_numbers):
            if self._list[i] <= max_width and vacNum == self.rigging[i]:
              # and vacNum > maxSingular:
                # maxSingular = vacNum
                maxPos = i
                break

        if maxPos == -1: # No singular values, then add a new row
            self._list.append(1)
            self.vacancy_numbers.append(None)
            # Go through our partition until we find a length of greater than 1
            i = len(self._list) - 1
            while i >= 0 and self._list[i] == 1:
                i -= 1
            self.rigging.insert(i + 1, None)
            return

        self._list[maxPos] += 1
        self.rigging[maxPos] = None # State that we've changed this row
        return self._list[maxPos] - 1

    def remove_cell(self, row):
        r"""
        Removes a cell at the specified ``row``.

        Note that :meth:`remove_cell` does not set/update the vacancy numbers or the
        riggings, but guarantees that the location has been allocated in the
        returned index.

        INPUT:

        - ``row`` -- The row to remove the cell from.

        OUTPUT:

        - The location of the newly constructed row or None if unable to remove row or if deleted a row.

        EXAMPLES::

            sage: RC = RiggedConfigurations(['A', 4, 1], [[2, 2]])
            sage: RP = RC(partition_list=[[2],[2,2],[2,1],[2]])[2]
            sage: RP.remove_cell(0)
            0
            sage: RP
            0[ ]0
            -1[ ]-1
            <BLANKLINE>
        """
        if row is None:
            return None

        if self._list[row] == 1:
            self._list.pop(row)
            self.vacancy_numbers.pop(row)
            self.rigging.pop(row)
            return None
        else:
            # Find the beginning of the next block
            block_len = self._list[row]
            if row + 1 == len(self._list):
                # If we are at the end, just do a simple remove
                self._list[row] -= 1
                return row
            else:
                for i in range(row + 1, len(self._list)):
                    if self._list[i] != block_len:
                        if i == row + 1:
                            # If the next row is a block change, just reduce by 1
                            self._list[row] -= 1
                            return row

                        # Otherwise we need to "move" the row
                        self._list.insert(i, block_len - 1)
                        # These should be updated (so there should be no need to carry them over)
                        self.vacancy_numbers.insert(i, None)
                        self.rigging.insert(i, None)

                        self._list.pop(row)
                        self.vacancy_numbers.pop(row)
                        self.rigging.pop(row)
                        return i - 1

                # We need to "move" the row to the end of the partition
                self._list.pop(row)
                self.vacancy_numbers.pop(row)
                self.rigging.pop(row)

                self._list.append(block_len - 1)
                # Placeholders as above
                self.vacancy_numbers.append(None)
                self.rigging.append(None)
                return len(self._list) - 1

