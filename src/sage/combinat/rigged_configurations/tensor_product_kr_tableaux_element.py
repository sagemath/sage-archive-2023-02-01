r"""
An element of a tensor product of Kirillov-Reshetikhin tableaux

A tensor product of
:class:`~sage.combinat.rigged_configurations.kr_tableaux.KirillovReshetikhinTableauxElement`.

AUTHORS:

- Travis Scrimshaw (2010-09-26): Initial version

EXAMPLES:

Type `A_n^{(1)}` examples::

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

Type `D_n^{(1)}` examples::

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
    sage: T.to_rigged_configuration().to_tensor_product_of_Kirillov_Reshetikhin_tableaux()
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
        Return the string representation for ``self``.

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
        Return the action of `e_i` on ``self``.

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
        Return the action of `f_i` on ``self``.

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
        Perform the bijection from ``self`` to a
        :class:`rigged configuration<sage.combinat.rigged_configurations.rigged_configuration_element.RiggedConfigurationElement>`
        which is described in [RigConBijection]_, [BijectionLRT]_, and
        [BijectionDn]_.

        INPUT:

        - ``display_steps`` -- (default: ``False``) Boolean which indicates
          if we want to output each step in the algorithm.

        OUTPUT:

        The rigged configuration corresponding to ``self``.

        EXAMPLES:

        Type `A_n^{(1)}` example::

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

        Type `D_n^{(1)}` example::

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

        Type `D_n^{(1)}` spinor example::

            sage: CP = TensorProductOfKirillovReshetikhinTableaux(['D', 5, 1], [[5,1],[2,1],[1,1],[1,1],[1,1]])
            sage: elt = CP(pathlist=[[-2,-5,4,3,1],[-1,2],[1],[1],[1]])
            sage: elt
            [[1], [3], [4], [-5], [-2]] (X) [[2], [-1]] (X) [[1]] (X) [[1]] (X) [[1]]
            sage: elt.to_rigged_configuration()
            <BLANKLINE>
            2[ ][ ]1
            <BLANKLINE>
            0[ ][ ]0
            0[ ]0
            <BLANKLINE>
            0[ ][ ]0
            0[ ]0
            <BLANKLINE>
            0[ ]0
            <BLANKLINE>
            0[ ][ ]0
            <BLANKLINE>

        This is invertible by calling
        :meth:`~sage.combinat.rigged_configurations.rigged_configuration_element.RiggedConfigurationElement.to_tensor_product_of_Kirillov_Reshetikhin_tableaux()`::

            sage: KRT = TensorProductOfKirillovReshetikhinTableaux(['D', 4, 1], [[2,2]])
            sage: T = KRT(pathlist=[[2,1,4,3]])
            sage: rc = T.to_rigged_configuration()
            sage: ret = rc.to_tensor_product_of_Kirillov_Reshetikhin_tableaux(); ret
            [[1, 3], [2, 4]]
            sage: ret == T
            True
        """

        bijection = KRTToRCBijection(self)
        type = self.parent().cartan_type().letter
        n = self.parent().cartan_type().n

        for cur_crystal in reversed(self):
            r = cur_crystal.parent().r()
            # Iterate through the columns
            for col_number, cur_column in enumerate(reversed(cur_crystal.to_array(False))):
                bijection.cur_path.insert(0, []) # Prepend an empty list

                # Check to see if we are a spinor column
                if type == 'D' and r >= n-1:
                    if display_steps:
                        print "===================="
                        print repr(TensorProductOfKirillovReshetikhinTableauxElement(self.parent(), *bijection.cur_path))
                        print "--------------------"
                        print repr(bijection.ret_rig_con)
                        print "--------------------\n"
                        print "Applied doubling map"
                    bijection.doubling_map()

                bijection.cur_dims.insert(0, [0, 1])

                # Note that we do not need to worry about iterating over columns
                #   (see previous note about the data structure).
                for letter in reversed(cur_column):
                    if bijection.cur_dims[0][0] < r:
                        bijection.cur_dims[0][0] += 1
                    val = letter.value # Convert from a CrystalOfLetter to an Integer

                    if display_steps:
                        print "===================="
                        print repr(TensorProductOfKirillovReshetikhinTableauxElement(self.parent(), *bijection.cur_path))
                        print "--------------------"
                        print repr(bijection.ret_rig_con)
                        print "--------------------\n"

                    # Build the next state
                    bijection.cur_path[0].insert(0, [letter]) # Prepend the value
                    bijection.next_state(val)

                # Check to see if we are a spinor column
                if type == 'D' and r >= n-1:
                    if display_steps:
                        print "===================="
                        print repr(TensorProductOfKirillovReshetikhinTableauxElement(self.parent(), *bijection.cur_path))
                        print "--------------------"
                        print repr(bijection.ret_rig_con)
                        print "--------------------\n"
                        print "Applied halving map"
                    bijection.halving_map()

                # If we've split off a column, we need to merge the current column
                #   to the current crystal tableau
                if col_number > 0:
                    for i, letter_singleton in enumerate(bijection.cur_path[0]):
                        bijection.cur_path[1][i].insert(0, letter_singleton[0])
                    bijection.cur_dims.pop(0)
                    bijection.cur_dims[0][1] += 1

                    # And perform the inverse column splitting map on the RC
                    for a in range(n):
                        bijection._update_vacancy_nums(a)

        bijection.ret_rig_con.set_immutable() # Return it to immutable
        return(bijection.ret_rig_con)

    def to_tensor_product_of_Kirillov_Reshetikhin_crystals(self):
        """
        Return a tensor product of Kirillov-Reshetikhin crystals corresponding
        to ``self``.

        This works by performing the filling map on each individual factor.
        For more on the filling map, see :class:`KirillovReshetikhinTableaux`.

        EXAMPLES::

            sage: KRT = TensorProductOfKirillovReshetikhinTableaux(['D',4,1], [[1,1],[2,2]])
            sage: elt = KRT(pathlist=[[-1],[-1,2,-1,1]]); elt
            [[-1]] (X) [[2, 1], [-1, -1]]
            sage: tp_krc = elt.to_tensor_product_of_Kirillov_Reshetikhin_crystals(); tp_krc
            [[[-1]], [[2], [-1]]]

        We can recover the original tensor product of KR tableaux::

            sage: KRT(tp_krc)
            [[-1]] (X) [[2, 1], [-1, -1]]
            sage: ret = KRT(*tp_krc); ret
            [[-1]] (X) [[2, 1], [-1, -1]]
            sage: ret == elt
            True
        """
        TP = self.parent().tensor_product_of_Kirillov_Reshetikhin_crystals()
        return TP(*[x.to_Kirillov_Reshetikhin_crystal() for x in self])

#    FIXME
#    def is_highest_weight(self, **options):
#        r"""
#        Checks to see if the crystal path is a highest weight crystal path.
#        """
#        return super(CrystalsElement, self).is_highest_weight( \
#          index_set=self.parent().cartan_type().classical().index_set(), **options )
#    # is_highest_weight()
