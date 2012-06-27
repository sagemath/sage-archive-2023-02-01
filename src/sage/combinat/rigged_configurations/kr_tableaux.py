r"""
Kirillov-Reshetikhin tableaux

Kirillov-Reshetikhin tableaux are rectangular tableaux with `r` rows and `s` columns
that naturally arise under the bijection between rigged configurations and tableaux
[RigConBijection]_. They are in bijection with the elements of the Kirillov-Reshetikhin crystal
`B^{r,s}` under the filling map defined in [AffineRigConDn]_.

AUTHORS:

- Travis Scrimshaw (2012-01-03): Initial version

.. TODO:: Implement bijection to Kirillov-Reshetikhin crystals.

EXAMPLES::

    sage: KRT = KirillovReshetikhinTableaux(['A', 4, 1], 2, 1)
    sage: elt = KRT([4, 3]); elt
    [[3], [4]]

    sage: KRT = KirillovReshetikhinTableaux(['D', 4, 1], 2, 1)
    sage: elt = KRT([-1, 1]); elt
    [[1], [-1]]
"""

#*****************************************************************************
#       Copyright (C) 2012 Travis Scrimshaw <tscrim@ucdavis.edu>
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

# This contains both the parent and element classes. These should be split if
#   the classes grow larger.

from sage.misc.cachefunc import cached_method
from sage.misc.flatten import flatten

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent

from sage.categories.finite_crystals import FiniteCrystals
from sage.categories.classical_crystals import ClassicalCrystals

from sage.rings.integer import Integer

from sage.combinat.crystals.letters import CrystalOfLetters
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.crystals.tensor_product import CrystalOfWords
from sage.combinat.crystals.tensor_product import TensorProductOfCrystalsElement
from sage.combinat.crystals.kirillov_reshetikhin import horizontal_dominoes_removed

class KirillovReshetikhinTableaux(CrystalOfWords):
    r"""
    Class of Kirillov-Reshetikhin tableaux.

    Kirillov-Reshetikhin tableaux are rectangular tableaux with `r` rows and `s` columns
    that naturally arise under the bijection between rigged configurations and tableaux
    [RigConBijection]_. They are in bijection with the elements of the Kirillov-Reshetikhin crystal
    `B^{r,s}` under the filling map defined in [AffineRigConDn]_.

    For more information, see
    :class:`TensorProductOfKirillovReshetikhinTableaux`.
    """

    @staticmethod
    def __classcall_private__(cls, cartan_type, r, s):
        """
        Normalize the input arguments to ensure unique representation.

        EXAMPLES::

            sage: KRT1 = KirillovReshetikhinTableaux(CartanType(['A',3,1]), 2, 3)
            sage: KRT2 = KirillovReshetikhinTableaux(['A',3,1], 2, 3)
            sage: KRT1 is KRT2
            True
        """
        cartan_type = CartanType(cartan_type)
        return super(KirillovReshetikhinTableaux, cls).__classcall__(cls, cartan_type, r, s)

    def __init__(self, cartan_type, r, s):
        r"""
        Initialize the KirillovReshetikhinTableaux class.

        INPUT:

        - ``cartan_type`` -- The Cartan type
        - ``r``           -- The number of rows
        - ``s``           -- The number of columns

        EXAMPLES::

            sage: KRT = KirillovReshetikhinTableaux(['A', 4, 1], 2, 3); KRT
            Kirillov-Reshetikhin tableaux of type ['A', 4, 1] and shape (2, 3)
            sage: TestSuite(KRT).run()
        """
        assert cartan_type.is_affine()

        self._r = r
        self._s = s
        Parent.__init__(self, category=FiniteCrystals())
        self.rename("Kirillov-Reshetikhin tableaux of type %s and shape (%d, %d)" % (cartan_type, r, s))

        self._cartan_type = cartan_type.classical()
        self.letters = CrystalOfLetters(cartan_type)

        if self._cartan_type.letter == 'A':
            self.module_generators = (self._fill([r] * s),)
        elif self._cartan_type.letter == 'D':
            self.module_generators = tuple(self._fill(shape) for
                                           shape in horizontal_dominoes_removed(s, r))

    def _element_constructor_(self, list):
        """
        Construct a KirillovReshetikhinTableauxElement.

        EXAMPLES::

            sage: KRT = KirillovReshetikhinTableaux(['A', 4, 1], 2, 1)
            sage: KRT([3, 4]) # indirect doctest
            [[4], [3]]
            sage: KRT([4, 3]) # indirect doctest
            [[3], [4]]
        """
        return self.element_class(self, list)

    def _fill(self, shape):
        r"""
        Return the highest weight KR tableau of weight ``shape``.

        INPUT:

        - ``shape`` -- The weight of the highest weight KR tableau.

        OUTPUT:

        - A `r \times s` tableau.

        EXAMPLES::

            sage: CPF = KirillovReshetikhinTableaux(['A', 3, 1], 2, 3)
            sage: CPF._fill([2, 2, 2])
            [[1, 1, 1], [2, 2, 2]]
            sage: CPF = KirillovReshetikhinTableaux(['D', 4, 1], 2, 1)
            sage: CPF._fill([])
            [[1], [-1]]
            sage: CPF = KirillovReshetikhinTableaux(['D', 14, 1], 12, 7)
            sage: CPF._fill([10,10,8,2,2,2])
            [[1, 1, 1, 1, 1, 7, 1], [2, 2, 2, 2, 2, 8, 2], [3, 3, 7, 9, 7, 9, 3], [4, 4, 8, 10, 8, 10, 4], [5, 5, 9, 11, 9, 11, 5], [6, 6, 10, 12, 10, 12, 6], [7, 7, 11, -12, 11, -12, 7], [8, 8, 12, -11, 12, -11, 8], [9, 9, -12, -10, -12, -10, 9], [10, 10, -11, -9, -11, -9, -9], [-12, 11, -10, -8, -10, -8, -8], [-11, 12, -9, -7, -9, -7, -7]]
            sage: CPF._fill([10,10,6,2,2,2])
            [[1, 1, 1, 1, 1, 5, 1], [2, 2, 2, 2, 2, 6, 2], [3, 3, 9, 7, 9, 7, 3], [4, 4, 10, 8, 10, 8, 4], [5, 5, 11, 9, 11, 9, 5], [6, 6, 12, 10, 12, 10, 6], [7, 7, -12, 11, -12, 11, 7], [8, 8, -11, 12, -11, 12, 8], [9, 9, -10, -12, -10, -12, -8], [10, 10, -9, -11, -9, -11, -7], [-12, 11, -8, -10, -8, -10, -6], [-11, 12, -7, -9, -7, -9, -5]]
        """
        # Add zeros until the shape has length s
        shape_list = list(shape) # Make sure we have a list
        while len(shape_list) != self._s:
            shape_list.append(0)

        tableau = []
        i = 0
        # Step 0 - Fill first columns of height r
        while i < self._s and shape_list[i] == self._r:
            tableau.append( [self._r - j for j in range(self._r)] )
            i += 1

        # Step 1 - Add the alternating columns until we hit an odd number of columns
        c = -1
        while i < self._s:
            # If it is an odd number of columns
            if i == self._s - 1 or shape_list[i] != shape_list[i+1]:
                c = shape_list[i]
                i += 1
                break
            temp_list = [-(shape_list[i] + j + 1) for j in range(self._r - shape_list[i])]
            for j in range(shape_list[i]):
                temp_list.append(shape_list[i] - j)
            tableau.append(temp_list)
            tableau.append( [self._r - j for j in range(self._r)] )
            i += 2

        # Step 2 - Add the x dependent columns
        x = c + 1
        while i < self._s:
            temp_list = [-x - j for j in range(self._r - x + 1)] # +1 for indexing
            for j in range(x - shape_list[i] - 1): # +1 for indexing
                temp_list.append(self._r - j)
            x = temp_list[-1] # This is the h+1 entry of the column
            for j in range(shape_list[i]):
                temp_list.append(shape_list[i] - j)

            tableau.append(temp_list)
            i += 1

        # Step 3 - Add the final column
        if c > -1:
            val = (self._r + x - 1) / 2
            temp_list = [-x - j for j in range(self._r - val)]
            for j in range(val):
                temp_list.append(val - j)
            tableau.append(temp_list)

        return self([self.letters(x) for x in flatten(tableau)])

    def r(self):
        """
        Return the number of rows in this tableaux class' rectangle.

        EXAMPLES::

            sage: KRT = KirillovReshetikhinTableaux(['A', 4, 1], 2, 1)
            sage: KRT.r()
            2
        """
        return self._r

    def s(self):
        """
        Return the number of columns in this tableaux class' rectangle.

        EXAMPLES::

            sage: KRT = KirillovReshetikhinTableaux(['A', 4, 1], 2, 1)
            sage: KRT.s()
            1
        """
        return self._s

class KirillovReshetikhinTableauxElement(TensorProductOfCrystalsElement):
    r"""
    A Kirillov-Reshetikhin tableau.

    For more information, see :class:`KirillovReshetikhinTableaux` and
    :class:`TensorProductOfKirillovReshetikhinTableaux`.
    """

    def __init__(self, parent, list):
        r"""
        Construct a KirillovReshetikhinTableauxElement.

        EXAMPLES::

            sage: KRT = KirillovReshetikhinTableaux(['A', 4, 1], 2, 1)
            sage: elt = KRT([4, 3]); elt
            [[3], [4]]
            sage: TestSuite(elt).run()
        """
        # Make sure we are a list of letters
        if list != [] and type(list[0]) is Integer:
            list = [parent.letters(x) for x in list]
        TensorProductOfCrystalsElement.__init__(self, parent, list)

    def _repr_(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: KRT = KirillovReshetikhinTableaux(['A', 4, 1], 2, 1)
            sage: KRT([3,2]) # indirect doctest
            [[2], [3]]
        """
        return repr(self.to_array())

    @cached_method
    def to_array(self, rows=True):
        """
        Return a 2-dimensional array representation of this Kirillov-Reshetikhin element.

        If the output is in rows, then it outputs the top row first (in the
        English convention) from left to right.

        For example: if the reading word is [2, 1, 4, 3], so as a 2x2 tableau\n
        1 3\n
        2 4\n
        we output [[1, 3], [2, 4]].

        If the output is in columns, then it outputs the leftmost column first
        with the largest element first. In other words this parses the reading
        word into its columns.

        Continuing with the previous example, the output would be
        [[2, 1], [4, 3]].

        INPUT:

        - ``rows`` -- (Default: True) Set to true if the resulting array is by row, otherwise it is by column.

        EXAMPLES::

            sage: KRT = KirillovReshetikhinTableaux(['A', 4, 1], 2, 2)
            sage: elt = KRT([4, 3, 2, 1])
            sage: elt.to_array()
            [[3, 1], [4, 2]]
            sage: elt.to_array(False)
            [[4, 3], [2, 1]]
        """

        retList = []
        r = self.parent().r()
        s = self.parent().s()
        if rows:
            for i in reversed(range(r)):
                row = []
                for j in range(s):
                    row.append(self[j * r + i])
                retList.append(row)
        else:
            for j in range(s):
                col = []
                for i in range(r):
                    col.append(self[j * r + i])
                retList.append(col)

        return retList

KirillovReshetikhinTableaux.Element = KirillovReshetikhinTableauxElement
