r"""
Kirillov-Reshetikhin tableaux

Kirillov-Reshetikhin tableaux are rectangular tableaux with `r` rows and
`s` columns that naturally arise under the bijection between rigged
configurations and tableaux [RigConBijection]_. They are in bijection with
the elements of the Kirillov-Reshetikhin crystal `B^{r,s}` under the (inverse)
filling map. For more information, see :class:`KirillovReshetikhinTableaux`
and :class:`TensorProductOfKirillovReshetikhinTableaux`.

AUTHORS:

- Travis Scrimshaw (2012-01-03): Initial version

- Travis Scrimshaw (2012-11-14): Added bijection to KR crystals
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
from sage.misc.abstract_method import abstract_method
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
from sage.combinat.crystals.kirillov_reshetikhin import horizontal_dominoes_removed, \
  KirillovReshetikhinCrystal, KirillovReshetikhinGenericCrystalElement
from sage.combinat.partition import Partition

class KirillovReshetikhinTableaux(CrystalOfWords):
    r"""
    Kirillov-Reshetikhin tableaux.

    Kirillov-Reshetikhin tableaux are rectangular tableaux with `r` rows and
    `s` columns that naturally arise under the bijection between rigged
    configurations and tableaux [RigConBijection]_. They are in bijection with
    the elements of the Kirillov-Reshetikhin crystal `B^{r,s}` under the
    (inverse) filling map.

    When the Kirillov-Reshetkihin crystal is a full `r \times s` rectangle
    (such as in type `A^{(1)}_n` or `C^{(1)}_n` for `B^{n,s}`), the filling
    map is trivial. Thus the highest weight module is just filled with columns
    `[1, 2, \ldots, n]`.

    For type `D^{(1)}_n` for `B^{r,s}` when `r \leq n-2`, the filling map is
    defined in [AffineRigConDn]_.

    For the spinor cases, the crystal `B^{k,s}` is isomorphic to the classical
    crystal `B(\Lambda_k)`, and here we consider the Kirillov-Reshetikhin
    tableaux as living in `B(2\Lambda_k)` under the natural doubling map.
    In this case, `e_i` and `f_i` act as `e_i^2` and `f_i^2` respectively.
    See [BijectionDn]_.

    .. WARNING::

        The module generators for all types except `A^{(1)}_n`, `D^{(1)}_n`,
        and full rectangles have not been tested, much less proven, to be the
        correct output from the bijection.

    For more information about the bijection between rigged configurations
    and tensor products of Kirillov-Reshetikhin tableaux, see
    :class:`TensorProductOfKirillovReshetikhinTableaux`.

    EXAMPLES::

        sage: KRT = KirillovReshetikhinTableaux(['A', 4, 1], 2, 1)
        sage: elt = KRT([4, 3]); elt
        [[3], [4]]

        sage: KRT = KirillovReshetikhinTableaux(['D', 4, 1], 2, 1)
        sage: elt = KRT([-1, 1]); elt
        [[1], [-1]]

    We can create highest weight crystals from a given shape::

        sage: KRT = KirillovReshetikhinTableaux(['D', 4, 1], 2, 1)
        sage: KRT(shape=[])
        [[1], [-1]]
        sage: KRT = KirillovReshetikhinTableaux(['D', 4, 1], 2, 2)
        sage: KRT(shape=[1,1])
        [[1, 1], [2, -1]]

    We can go between :func:`KirillovReshetikhinCrystal` and
    :class:`KirillovReshetikhinTableaux` elements::

        sage: KRCrys = KirillovReshetikhinCrystal(['D', 4, 1], 2, 2)
        sage: KRTab = KirillovReshetikhinTableaux(['D', 4, 1], 2, 2)
        sage: elt = KRCrys(3, 2); elt
        [[2], [3]]
        sage: k = KRTab(elt); k
        [[2, 1], [3, -1]]
        sage: k.to_Kirillov_Reshetikhin_crystal()
        [[2], [3]]
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
        ct = CartanType(cartan_type)
        assert ct.is_affine()

        if ct.is_untwisted_affine():
            if ct.letter == 'D':
                if r == ct.n or r == ct.n - 1:
                    return KRTableauxSpin(ct, r, s)
                return KRTableauxTypeVertical(ct, r, s)

            if ct.letter == 'B':
                if r == ct.n:
                    return KRTableauxBn(ct, r, s)
                return KRTypeVertical(ct, r, s)

            if ct.letter == 'A' or (ct.letter == 'C' and r == ct.n):
                return KRTableauxRectangle(ct, r, s)
        else:
            if ct.dual().letter == 'B':
                return KRTableauxTypeVertical(ct, r, s)

        raise NotImplementedError
        #return super(KirillovReshetikhinTableaux, cls).__classcall__(cls, ct, r, s)

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
            sage: TestSuite(KRT).run()  # long time (4s on sage.math, 2013)
            sage: KRT = KirillovReshetikhinTableaux(['D', 4, 1], 2, 3); KRT
            Kirillov-Reshetikhin tableaux of type ['D', 4, 1] and shape (2, 3)
            sage: TestSuite(KRT).run()  # long time (53s on sage.math, 2013)
            sage: KRT = KirillovReshetikhinTableaux(['D', 4, 1], 4, 1); KRT
            Kirillov-Reshetikhin tableaux of type ['D', 4, 1] and shape (4, 1)
            sage: TestSuite(KRT).run()
        """
        self._r = r
        self._s = s
        Parent.__init__(self, category=FiniteCrystals())
        self.rename("Kirillov-Reshetikhin tableaux of type %s and shape (%d, %d)" % (cartan_type, r, s))

        self._cartan_type = cartan_type.classical()
        self.letters = CrystalOfLetters(self._cartan_type)

        self.module_generators = self._build_module_generators()

    @abstract_method
    def _build_module_generators(self):
        """
        Build the module generators.

        EXAMPLES::

            sage: KRT = KirillovReshetikhinTableaux(['A', 4, 1], 2, 3)
            sage: KRT._build_module_generators()
            ([[1, 1, 1], [2, 2, 2]],)
        """
#        shapes = self.Kirillov_Reshetikhin_crystal().classical_decomposition().shapes
#        self.module_generators = tuple(self._fill(Partition(shape).conjugate()) for shape in shapes)

    def _element_constructor_(self, list, **options):
        """
        Construct a :class:`KirillovReshetikhinTableauxElement`.

        EXAMPLES::

            sage: KRT = KirillovReshetikhinTableaux(['A', 4, 1], 2, 1)
            sage: KRT([3, 4]) # indirect doctest
            [[4], [3]]
            sage: KRT([4, 3])
            [[3], [4]]
        """
        return self.element_class(self, list, **options)

    def r(self):
        """
        Return the value `r` for this tableaux class.

        EXAMPLES::

            sage: KRT = KirillovReshetikhinTableaux(['A', 4, 1], 2, 1)
            sage: KRT.r()
            2
        """
        return self._r

    def s(self):
        """
        Return the value `s` for this tableaux class.

        EXAMPLES::

            sage: KRT = KirillovReshetikhinTableaux(['A', 4, 1], 2, 1)
            sage: KRT.s()
            1
        """
        return self._s

    def Kirillov_Reshetikhin_crystal(self):
        """
        Return the corresponding
        :func:`Kirillov-Reshetikhin crystal<KirillovReshetikhinCrystal>`.

        EXAMPLES::

            sage: KirillovReshetikhinTableaux(['A', 4, 1], 2, 1).Kirillov_Reshetikhin_crystal()
            Kirillov-Reshetikhin crystal of type ['A', 4, 1] with (r,s)=(2,1)
        """
        return KirillovReshetikhinCrystal(self._cartan_type.affine(),
                                          self._r, self._s)

class KRTableauxRectangle(KirillovReshetikhinTableaux):
    r"""
    Kirillov-Reshetkhin tableaux `B^{r,s}` whose module generator is a single
    `r \times s` rectangle.

    These are Kirillov-Reshetkhin tableaux `B^{r,s}` of type:

    - `A_n^{(1)}` for all `1 \leq r \leq n`,
    - `C_n^{(1)}` when `r = n`.
    """
    def _build_module_generators(self):
        r"""
        Build the module generators.

        There is only one module generator which corresponds to a single
        `r \times s` rectangle.

        EXAMPLES::

            sage: KRT = KirillovReshetikhinTableaux(['A', 4, 1], 2, 3)
            sage: KRT._build_module_generators()
            ([[1, 1, 1], [2, 2, 2]],)
        """
        tableau = []
        for i in range(self._s):
            tableau.append( [self._r - j for j in range(self._r)] )

        return (self([self.letters(x) for x in flatten(tableau)]),)

    def _element_constructor_(self, list, **options):
        """
        Construct a :class:`KirillovReshetikhinTableauxElement`.

        EXAMPLES::

            sage: KRT = KirillovReshetikhinTableaux(['A', 4, 1], 2, 1)
            sage: KRT([3, 4]) # indirect doctest
            [[4], [3]]
            sage: KRT([4, 3])
            [[3], [4]]
        """
        if isinstance(list, KirillovReshetikhinGenericCrystalElement):
            # Check to make sure it can be converted
            if list.cartan_type() != self.cartan_type().affine() \
              or list.parent().r() != self._r or list.parent().s() != self._s:
                raise ValueError("The Kirillov-Reshetikhin crystal must have the same Cartan type and shape")

            # To build a KR tableau from a KR crystal:
            # 1 - start with the highest weight KR tableau
            # 2 - determine a path from the KR crystal to its highest weight
            # 3 - apply the inverse path to the highest weight KR tableau
            f_str = reversed(list.lift().to_highest_weight()[1])
            return self.module_generators[0].f_string(f_str)

        return KirillovReshetikhinTableaux._element_constructor_(self, list, **options)

class KRTableauxTypeVertical(KirillovReshetikhinTableaux):
    r"""
    Kirillov-Reshetkihn tableaux `B^{r,s}` of type:

    - `D_n^{(1)}` for all `1 \leq r < n-1`,
    - `B_n^{(1)}` for all `1 \leq r < n`,
    - `A_{2n-1}^{(2)}` for all `1 \leq r \leq n`.
    """

    def _fill(self, shape):
        r"""
        Return the highest weight KR tableau of weight ``shape``.

        INPUT:

        - ``shape`` -- The weight of the highest weight KR tableau (the
          conjugate of the shape of the KR crystal's tableau)

        OUTPUT:

        - A `r \times s` tableau

        EXAMPLES::

            sage: KRT = KirillovReshetikhinTableaux(['D', 4, 1], 2, 1)
            sage: KRT._fill([])
            [[1], [-1]]
            sage: KRT = KirillovReshetikhinTableaux(['D', 14, 1], 12, 7)
            sage: KRT._fill([10,10,8,2,2,2])
            [[1, 1, 1, 1, 1, 7, 1], [2, 2, 2, 2, 2, 8, 2], [3, 3, 7, 9, 7, 9, 3], [4, 4, 8, 10, 8, 10, 4], [5, 5, 9, 11, 9, 11, 5], [6, 6, 10, 12, 10, 12, 6], [7, 7, 11, -12, 11, -12, 7], [8, 8, 12, -11, 12, -11, 8], [9, 9, -12, -10, -12, -10, 9], [10, 10, -11, -9, -11, -9, -9], [-12, 11, -10, -8, -10, -8, -8], [-11, 12, -9, -7, -9, -7, -7]]
            sage: KRT._fill([10,10,6,2,2,2])
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

    def _build_module_generators(self):
        """
        Build the module generators.

        EXAMPLES::

            sage: KRT = KirillovReshetikhinTableaux(['D',4,1], 2, 3)
            sage: KRT._build_module_generators()
            ([[-2, 1, 1], [-1, 2, -1]], [[1, -2, 1], [2, -1, 2]], [[1, 1, 1], [2, 2, -1]], [[1, 1, 1], [2, 2, 2]])
        """
        return tuple(self._fill(shape) for shape in
                     horizontal_dominoes_removed(self._s, self._r))

    def _element_constructor_(self, list, **options):
        """
        Construct a :class:`KirillovReshetikhinTableauxElement`.

        EXAMPLES::

            sage: KRT = KirillovReshetikhinTableaux(['A', 4, 1], 2, 1)
            sage: KRT([3, 4]) # indirect doctest
            [[4], [3]]
            sage: KRT([4, 3])
            [[3], [4]]
        """
        if "shape" in options:
            return self._fill(Partition(options["shape"]).conjugate())

        if isinstance(list, KirillovReshetikhinGenericCrystalElement):
            # Check to make sure it can be converted
            if list.cartan_type() != self.cartan_type().affine() \
              or list.parent().r() != self._r or list.parent().s() != self._s:
                raise ValueError("The Kirillov-Reshetikhin crystal must have the same Cartan type and shape")

            # To build a KR tableau from a KR crystal:
            # 1 - start with a highest weight KR tableau generated from the
            #  shape of the KR crystal
            # 2 - determine a path from the KR crystal to its highest weight
            # 3 - apply the inverse path to the highest weight KR tableau
            lifted = list.lift()
            shape = lifted.to_tableau().shape().conjugate()
            f_str = reversed(lifted.to_highest_weight()[1])
            return self._fill(shape).f_string(f_str)

        return KirillovReshetikhinTableaux._element_constructor_(self, list, **options)

class KRTableauxSpin(KRTableauxRectangle):
    r"""
    Kirillov-Reshetikhin tableaux `B^{r,s}` of type `D_n^{(1)}` with
    `r = n, n-1`.
    """
    def _build_module_generators(self):
        r"""
        Build the module generators.

        There is only one module generator which corresponds to a single
        `n \times s` rectangle.

        EXAMPLES::

            sage: KRT = KirillovReshetikhinTableaux(['D', 4, 1], 3, 3)
            sage: KRT._build_module_generators()
            ([[1, 1, 1], [2, 2, 2], [3, 3, 3], [-4, -4, -4]],)
        """
        if self._r == self.cartan_type().n:
            return KRTableauxRectangle._build_module_generators(self)

        tableau = []
        for i in range(self._s):
            tableau.append( [-4] + [self._r - j for j in range(self._r)] )

        return (self([self.letters(x) for x in flatten(tableau)]),)

# Placeholder for type B_n spinors
#class KRTableauxBn(KRTableauxSpin):
#    """
#    Kirillov-Reshetkhin tableaux `B^{n,s}` of type `B_n^{(1)}`.
#    """
#    def _build_module_generators(self):
#        """
#        Build the module generators.
#        """
#        shapes = KirillovReshetikhinCrystal(cartan_type, r, s).classical_decomposition().shapes
#        shapes = [Partition(map(lambda x: Integer(x*2), shape)).conjugate() for shape in shapes]
#        self.module_generators = tuple(self._fill(Partition(shape).conjugate()) for shape in shapes)

class KirillovReshetikhinTableauxElement(TensorProductOfCrystalsElement):
    r"""
    A Kirillov-Reshetikhin tableau.

    For more information, see :class:`KirillovReshetikhinTableaux` and
    :class:`TensorProductOfKirillovReshetikhinTableaux`.
    """

    def __init__(self, parent, list, **options):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: KRT = KirillovReshetikhinTableaux(['A', 4, 1], 2, 1)
            sage: elt = KRT([4, 3]); elt
            [[3], [4]]
            sage: TestSuite(elt).run()
        """
        # Make sure we are a list of letters
        if list != [] and type(list[0]) is not CrystalOfLetters:
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

    def _latex_(self):
        """
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: KRT = KirillovReshetikhinTableaux(['A', 4, 1], 2, 3)
            sage: latex(KRT([3,2,4,2,4,3])) # indirect doctest
            {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
            \raisebox{-.6ex}{$\begin{array}[b]{*{3}c}\cline{1-3}
            \lr{2}&\lr{2}&\lr{3}\\\cline{1-3}
            \lr{3}&\lr{4}&\lr{4}\\\cline{1-3}
            \end{array}$}
            }
        """
        from sage.combinat.output import tex_from_array
        return tex_from_array(self.to_array())

    def to_Kirillov_Reshetikhin_crystal(self):
        r"""
        Construct a :func:`KirillovReshetikhinCrystal` element from ``self``.

        We construct the Kirillov-Reshetikhin crystal element as follows:

        1. Determine the shape `\lambda` of the KR crystal from the weight.
        2. Determine a path `e_{i_1} e_{i_2} \cdots e_{i_k}` to the highest
           weight.
        3. Apply `f_{i_k} \cdots f_{i_2} f_{i_1}` to a highest weight KR
           crystal of shape `\lambda`.

        EXAMPLES::

            sage: KRT = KirillovReshetikhinTableaux(['D',4,1], 2,2)
            sage: elt = KRT([3,2,-1,1]); elt
            [[2, 1], [3, -1]]
            sage: elt.to_Kirillov_Reshetikhin_crystal()
            [[2], [3]]

        TESTS:

        Spinor tests::

            sage: KRT = KirillovReshetikhinTableaux(['D',4,1], 4, 3)
            sage: KRC = KirillovReshetikhinCrystal(['D',4,1], 4, 3)
            sage: elt = KRT([-3,-4,2,1,-3,-4,2,1,-2,-4,3,1]); elt
            [[1, 1, 1], [2, 2, 3], [-4, -4, -4], [-3, -3, -2]]
            sage: ret = elt.to_Kirillov_Reshetikhin_crystal(); ret
            [++--, [[1], [3], [-4], [-3]]]
            sage: test = KRT(ret); test
            [[1, 1, 1], [2, 2, 3], [-4, -4, -4], [-3, -3, -2]]
            sage: test == elt
            True
        """
        return self.parent().Kirillov_Reshetikhin_crystal()(self)

    @cached_method
    def to_array(self, rows=True):
        r"""
        Return a 2-dimensional array representation of this
        Kirillov-Reshetikhin element.

        If the output is in rows, then it outputs the top row first (in the
        English convention) from left to right.

        For example: if the reading word is `[2, 1, 4, 3]`, so as a
        `2 \times 2` tableau::

            1 3
            2 4

        we output ``[[1, 3], [2, 4]]``.

        If the output is in columns, then it outputs the leftmost column first
        with the bottom element first. In other words this parses the reading
        word into its columns.

        Continuing with the previous example, the output would be
        ``[[2, 1], [4, 3]]``.

        INPUT:

        - ``rows`` -- (Default: ``True``) Set to ``True`` if the resulting
          array is by row, otherwise it is by column.

        EXAMPLES::

            sage: KRT = KirillovReshetikhinTableaux(['A', 4, 1], 2, 2)
            sage: elt = KRT([4, 3, 2, 1])
            sage: elt.to_array()
            [[3, 1], [4, 2]]
            sage: elt.to_array(False)
            [[4, 3], [2, 1]]
        """

        ret_list = []
        h = self.parent()._r
        s = self.parent()._s
        if rows:
            for i in reversed(range(h)):
                row = []
                for j in range(s):
                    row.append(self[j * h + i])
                ret_list.append(row)
        else:
            for j in range(s):
                col = []
                for i in range(h):
                    col.append(self[j * h + i])
                ret_list.append(col)

        return ret_list

    def to_classical_highest_weight(self, index_set=None):
        r"""
        Return the classical highest weight element corresponding to ``self``.

        INPUT:

        - ``index_set`` -- (Default: ``None``) Return the highest weight
          with respect to the index set. If ``None`` is passed in, then this
          uses the classical index set.

        OUTPUT:

        A pair ``[H, f_str]`` where ``H`` is the highest weight element and
        ``f_str`` is a list of `a_i` of `f_{a_i}` needed to reach ``H``.

        EXAMPLES::

            sage: KRTab = KirillovReshetikhinTableaux(['D',4,1], 2,2)
            sage: elt = KRTab([3,2,-1,1]); elt
            [[2, 1], [3, -1]]
            sage: elt.to_classical_highest_weight()
            [[[1, 1], [2, -1]], [1, 2]]
        """
        if index_set is None:
            index_set = self.index_set()
        for i in index_set:
            if self.epsilon(i) != 0:
                self = self.e(i)
                hw = self.to_classical_highest_weight(index_set=index_set)
                return [hw[0], [i] + hw[1]]
        return [self, []]

    def classical_weight(self):
        r"""
        Return the classical weight of ``self``.

        EXAMPLES::

            sage: KRT = KirillovReshetikhinTableaux(['D',4,1], 2,2)
            sage: elt = KRT([3,2,-1,1]); elt
            [[2, 1], [3, -1]]
            sage: elt.classical_weight()
            (0, 1, 1, 0)
        """
        return self.Phi() - self.Epsilon()

    def e(self, i):
        """
        Perform the action of `e_i` on ``self``.

        .. TODO::

            Implement a direct action of `e_0` without moving to KR crystals.

        EXAMPLES::

            sage: KRT = KirillovReshetikhinTableaux(['D',4,1], 2,2)
            sage: KRT.module_generators[0].e(0)
            [[-2, 1], [-1, -1]]
        """
        if i == 0:
            ret = self.to_Kirillov_Reshetikhin_crystal().e0()
            if ret is None:
                return None
            return ret.to_Kirillov_Reshetikhin_tableau()
        return TensorProductOfCrystalsElement.e(self, i)

    def f(self, i):
        """
        Perform the action of `f_i` on ``self``.

        .. TODO::

            Implement a direct action of `f_0` without moving to KR crystals.

        EXAMPLES::

            sage: KRT = KirillovReshetikhinTableaux(['D',4,1], 2,2)
            sage: KRT.module_generators[0].f(0)
            [[1, 1], [2, -1]]
        """
        if i == 0:
            ret = self.to_Kirillov_Reshetikhin_crystal().f0()
            if ret is None:
                return None
            return ret.to_Kirillov_Reshetikhin_tableau()
        return TensorProductOfCrystalsElement.f(self, i)

KirillovReshetikhinTableaux.Element = KirillovReshetikhinTableauxElement

class KRTableauxSpinElement(KirillovReshetikhinTableauxElement):
    r"""
    Kirillov-Reshetikhin tableau for spinors.

    Here we are in the embedding `B(\Lambda_n) \to B(2 \Lambda_n)`, so `e_i`
    and `f_i` act by `e_i^2` and `f_i^2` respectively.
    """
    def e(self, i):
        r"""
        Calculate the action of `e_i` on ``self``.

        EXAMPLES::

            sage: KRT = KirillovReshetikhinTableaux(['D',4,1], 4, 1)
            sage: KRT([-1, -4, 3, 2]).e(1)
            [[1], [3], [-4], [-2]]
            sage: KRT([-1, -4, 3, 2]).e(3)
        """
        half = TensorProductOfCrystalsElement.e(self, i)
        if half is None:
            return None
        return TensorProductOfCrystalsElement.e(half, i)

    def f(self, i):
        r"""
        Calculate the action of `f_i` on ``self``.

        EXAMPLES::

            sage: KRT = KirillovReshetikhinTableaux(['D',4,1], 4, 1)
            sage: KRT([-1, -4, 3, 2]).f(1)
            sage: KRT([-1, -4, 3, 2]).f(3)
            [[2], [4], [-3], [-1]]
        """
        half = TensorProductOfCrystalsElement.f(self, i)
        if half is None:
            return None

        return TensorProductOfCrystalsElement.f(half, i)

    def epsilon(self, i):
        r"""
        Compute `\epsilon_i` of ``self``.

        EXAMPLES::

            sage: KRT = KirillovReshetikhinTableaux(['D',4,1], 4, 1)
            sage: KRT([-1, -4, 3, 2]).epsilon(1)
            1
            sage: KRT([-1, -4, 3, 2]).epsilon(3)
            0
        """
        return TensorProductOfCrystalsElement.epsilon(self, i) / 2

    def phi(self, i):
        r"""
        Compute `\phi_i` of ``self``.

        EXAMPLES::

            sage: KRT = KirillovReshetikhinTableaux(['D',4,1], 4, 1)
            sage: KRT([-1, -4, 3, 2]).phi(1)
            0
            sage: KRT([-1, -4, 3, 2]).phi(3)
            1
        """
        return TensorProductOfCrystalsElement.phi(self, i) / 2

    @cached_method
    def to_array(self, rows=True):
        r"""
        Return a 2-dimensional array representation of this
        Kirillov-Reshetikhin element.

        If the output is in rows, then it outputs the top row first (in the
        English convention) from left to right.

        For example: if the reading word is `[2, 1, 4, 3]`, so as a
        `2 \times 2` tableau::

            1 3
            2 4

        we output ``[[1, 3], [2, 4]]``.

        If the output is in columns, then it outputs the leftmost column first
        with the bottom element first. In other words this parses the reading
        word into its columns.

        Continuing with the previous example, the output would be
        ``[[2, 1], [4, 3]]``.

        INPUT:

        - ``rows`` -- (Default: ``True``) Set to ``True`` if the resulting
          array is by row, otherwise it is by column.

        EXAMPLES::

            sage: KRT = KirillovReshetikhinTableaux(['D', 4, 1], 4, 3)
            sage: elt = KRT([-3,-4,2,1,-3,-4,2,1,-2,-4,3,1])
            sage: elt.to_array()
            [[1, 1, 1], [2, 2, 3], [-4, -4, -4], [-3, -3, -2]]
            sage: elt.to_array(False)
            [[-3, -4, 2, 1], [-3, -4, 2, 1], [-2, -4, 3, 1]]
        """

        ret_list = []
        h = self.parent()._cartan_type.n
        s = self.parent()._s
        if rows:
            for i in reversed(range(h)):
                row = []
                for j in range(s):
                    row.append(self[j * h + i])
                ret_list.append(row)
        else:
            for j in range(s):
                col = []
                for i in range(h):
                    col.append(self[j * h + i])
                ret_list.append(col)

        return ret_list

KRTableauxSpin.Element = KRTableauxSpinElement
