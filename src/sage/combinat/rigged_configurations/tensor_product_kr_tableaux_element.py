r"""
A specific tensor product of Kirillov-Reshetikhin tableaux

A tensor product of :class:`KirillovReshetikhinTableauxElement`.

AUTHORS:

- Travis Scrimshaw (2010-09-26): Initial version

.. TODO:: A proper is_highest_weight() function without needing to change the index set.

EXAMPLES:

Type `A_n` examples::

    sage: KRT = TensorProductOfKirillovReshetikhinTableaux(['A', 3, 1], [[1,1], [2,1], [1,1], [2,1], [2,1], [2,1]])
    sage: T = KRT(pathlist=[[2], [4,1], [3], [4,2], [3,1], [2,1]])
    sage: T
    [[2]] (X) [[1], [4]] (X) [[3]] (X) [[2], [4]] (X) [[1], [3]] (X) [[1], [2]]
    sage: T.to_rigged_configuration()
    <BLANKLINE>
    0[ ][ ]0
    1[ ]1
    <BLANKLINE>
    1[ ][ ]0
    1[ ]0
    1[ ]0
    <BLANKLINE>
    0[ ][ ]0
    <BLANKLINE>
    sage: T = KRT(pathlist=[[1], [2,1], [1], [4,1], [3,1], [2,1]])
    sage: T
    [[1]] (X) [[1], [2]] (X) [[1]] (X) [[1], [4]] (X) [[1], [3]] (X) [[1], [2]]
    sage: T.to_rigged_configuration()
    <BLANKLINE>
    (/)
    <BLANKLINE>
    1[ ]0
    1[ ]0
    <BLANKLINE>
    0[ ]0
    <BLANKLINE>

Type `D_n` examples::

    sage: KRT = TensorProductOfKirillovReshetikhinTableaux(['D', 4, 1], [[1,1], [1,1], [1,1], [1,1]])
    sage: T = KRT(pathlist=[[-1], [-1], [1], [1]])
    sage: T
    [[-1]] (X) [[-1]] (X) [[1]] (X) [[1]]
    sage: T.to_rigged_configuration()
    <BLANKLINE>
    0[ ][ ]0
    0[ ][ ]0
    <BLANKLINE>
    0[ ][ ]0
    0[ ][ ]0
    <BLANKLINE>
    0[ ][ ]0
    <BLANKLINE>
    0[ ][ ]0
    <BLANKLINE>
    sage: KRT = TensorProductOfKirillovReshetikhinTableaux(['D', 4, 1], [[2,1], [1,1], [1,1], [1,1]])
    sage: T = KRT(pathlist=[[3,2], [1], [-1], [1]])
    sage: T
    [[2], [3]] (X) [[1]] (X) [[-1]] (X) [[1]]
    sage: T.to_rigged_configuration()
    <BLANKLINE>
    0[ ]0
    0[ ]0
    0[ ]0
    <BLANKLINE>
    0[ ]0
    0[ ]0
    0[ ]0
    <BLANKLINE>
    1[ ]0
    <BLANKLINE>
    1[ ]0
    <BLANKLINE>
    sage: T.to_rigged_configuration().to_Kirillov_Reshetikhin_tableaux()
    [[2], [3]] (X) [[1]] (X) [[-1]] (X) [[1]]
"""

#*****************************************************************************
#       Copyright (C) 2010, 2011, 2012 Travis Scrimshaw <tscrim@ucdavis.edu>
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

from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.crystals.tensor_product import TensorProductOfCrystalsElement
from sage.combinat.rigged_configurations.bijection import KRTToRCBijection

class TensorProductOfKirillovReshetikhinTableauxElement(TensorProductOfCrystalsElement):
    """
    An element in a tensor product of Kirillov-Reshetikhin tableaux.

    For more on tensor product of Kirillov-Reshetikhin tableaux, see
    :class:`TensorProductOfKirillovReshetikhinTableaux`.
    """

    # Functions

    def __init__(self, parent, *path, **options):
        r"""
        Construct a TensorProductOfKirillovReshetikhinTableauxElement.

        INPUT:

        - ``parent`` -- Parent for this element
        - ``path``   -- The list of KR tableaux elements

        EXAMPLES::

            sage: KRT = TensorProductOfKirillovReshetikhinTableaux(['A', 3, 1], [[1, 1], [2, 1], [1, 1], [2, 1], [2, 1], [2, 1]])
            sage: T = KRT(pathlist=[[2], [4, 1], [3], [4, 2], [3, 1], [2, 1]])
            sage: T
            [[2]] (X) [[1], [4]] (X) [[3]] (X) [[2], [4]] (X) [[1], [3]] (X) [[1], [2]]
            sage: KRT = TensorProductOfKirillovReshetikhinTableaux(['A', 3, 1], [[3,3], [2,1]])
            sage: T = KRT(pathlist=[[3, 2, 1, 4, 2, 1, 4, 3, 1], [2, 1]])
            sage: T
            [[1, 1, 1], [2, 2, 3], [3, 4, 4]] (X) [[1], [2]]
            sage: KRT = TensorProductOfKirillovReshetikhinTableaux(['D', 4, 1], [[2, 1], [1, 1], [1, 1], [1, 1]])
            sage: T = KRT(pathlist=[[3,2], [1], [-1], [1]])
            sage: T
            [[2], [3]] (X) [[1]] (X) [[-1]] (X) [[1]]
            sage: TestSuite(T).run()
        """

        if "pathlist" in options:
            pathlist = options["pathlist"]
            TensorProductOfCrystalsElement.__init__(self, parent,
              [parent.crystals[i](list(tab)) for i, tab in enumerate(pathlist)])
        elif "list" in options:
            TensorProductOfCrystalsElement.__init__(self, parent, **options)
        else:
            TensorProductOfCrystalsElement.__init__(self, parent, list(path))

    def _repr_(self):
        """
        Return the string representation for self.

        EXAMPLES::

            sage: KRT = TensorProductOfKirillovReshetikhinTableaux(['D', 4, 1], [[2, 1], [1, 1], [1, 1], [1, 1]])
            sage: T = KRT(pathlist=[[3,2], [1], [-1], [1]])
            sage: T # indirect doctest
            [[2], [3]] (X) [[1]] (X) [[-1]] (X) [[1]]
        """
        retStr = repr(self[0])
        for i in range(1, len(self)):
            retStr += " (X) " + repr(self[i])
        return(retStr)

    def e(self, i):
        r"""
        Return the action of `e_i` on self.

        EXAMPLES::

            sage: KRT = TensorProductOfKirillovReshetikhinTableaux(['D', 4, 1], [[2,1]])
            sage: T = KRT(pathlist=[[4,3]])
            sage: T.e(1)
            sage: T.e(2)
            [[2], [4]]
        """

        if i != 0:
            return TensorProductOfCrystalsElement.e(self, i)

        return None

    def f(self, i):
        r"""
        Return the action of `f_i` on self.

        EXAMPLES::

            sage: KRT = TensorProductOfKirillovReshetikhinTableaux(['D', 4, 1], [[2,1]])
            sage: T = KRT(pathlist=[[4,3]])
            sage: T.f(1)
            sage: T.f(4)
            [[-4], [4]]
        """
        if i != 0:
            return TensorProductOfCrystalsElement.f(self, i)

        return None

    def to_rigged_configuration(self, display_steps=False):
        r"""
        Perform the bijection from this to a
        :class:`RiggedConfiguration` which is described
        in [RigConBijection]_, [BijectionLRT]_, and [BijectionDn]_.

        INPUT:

        - ``display_steps`` -- (default: False) Boolean which indicates if we want to output each step in the algorithm.

        OUTPUT:

        The rigged configuration corresponding to self.

        EXAMPLES::

            sage: KRT = HighestWeightTensorProductOfKirillovReshetikhinTableaux(['A', 3, 1], [[2,1], [2,1], [2,1]])
            sage: T = KRT(pathlist=[[4, 2], [3, 1], [2, 1]])
            sage: T
            [[2], [4]] (X) [[1], [3]] (X) [[1], [2]]
            sage: T.to_rigged_configuration()
            <BLANKLINE>
            0[ ]0
            <BLANKLINE>
            1[ ]1
            1[ ]0
            <BLANKLINE>
            0[ ]0
            <BLANKLINE>
            sage: KRT = TensorProductOfKirillovReshetikhinTableaux(['D', 4, 1], [[2,2]])
            sage: T = KRT(pathlist=[[2,1,4,3]])
            sage: T
            [[1, 3], [2, 4]]
            sage: T.to_rigged_configuration()
            <BLANKLINE>
            0[ ]0
            <BLANKLINE>
            -1[ ]-1
            -1[ ]-1
            <BLANKLINE>
            0[ ]0
            <BLANKLINE>
            (/)
        """

        bijection = KRTToRCBijection(self)

        for cur_crystal in reversed(self):
            # Iterate through the columns
            for col_number, cur_column in enumerate(reversed(cur_crystal.to_array(False))):
                bijection.cur_path.insert(0, []) # Prepend an empty list

                # Note that we do not need to worry about iterating over columns
                #   (see previous note about the data structure).
                # height is the height of the current tableau
                for height, letter in enumerate(reversed(cur_column)):
                    val = letter.value # Convert from a CrystalOfLetter to an Integer

                    if display_steps:
                        print "===================="
                        print repr(TensorProductOfKirillovReshetikhinTableauxElement(self.parent(), *bijection.cur_path))
                        print "--------------------"
                        print repr(bijection.ret_rig_con)
                        print "--------------------\n"

                    # Build the next state
                    # FIXME: Remove the single object list around letter
                    bijection.cur_path[0].insert(0, [letter]) # Prepend the value
                    bijection.next_state(val, height)

                # If we've split off a column, we need to merge the current column
                #   to the current crystal tableau
                if col_number > 0:
                    for i, letter_singleton in enumerate(bijection.cur_path[0]):
                        bijection.cur_path[1][i].insert(0, letter_singleton[0])
                    bijection.cur_path.pop(0)

                    # And perform the inverse column splitting map on the RC
                    for a in range(self.parent().cartan_type().n):
                        bijection._update_vacancy_nums(a)

        bijection.ret_rig_con.set_immutable() # Return it to immutable
        return(bijection.ret_rig_con)

#    FIXME
#    def is_highest_weight(self, **options):
#        r"""
#        Checks to see if the crystal path is a highest weight crystal path.
#        """
#        return super(CrystalsElement, self).is_highest_weight( \
#          index_set=self.parent().cartan_type().classical().index_set(), **options )
#    # is_highest_weight()
