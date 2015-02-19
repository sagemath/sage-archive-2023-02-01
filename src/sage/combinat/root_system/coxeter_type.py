"""
Coxeter Types
"""
#*****************************************************************************
#       Copyright (C) 2015 Travis Scrimshaw <tscrim at ucdavis.edu>,
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
from cartan_type import CartanType
from sage.matrix.all import MatrixSpace
from sage.rings.all import ZZ

class CoxeterType_abstract(object):
    """
    Abstract class for Coxeter types.
    """
    @abstract_method
    def rank(self):
        """
        Return the rank of ``self``.

        This is the number of nodes of the associated Coxeter diagram.

        EXAMPLES::

            sage: CartanType(['A', 4]).rank()
            4
            sage: CartanType(['A', 7, 2]).rank()
            5
            sage: CartanType(['I', 8]).rank()
            2
        """

    @abstract_method
    def index_set(self):
        """
        Return the index set for ``self``.

        This is the list of the nodes of the associated Coxeter diagram.

        EXAMPLES::

            sage: CartanType(['A', 3, 1]).index_set()
            (0, 1, 2, 3)
            sage: CartanType(['D', 4]).index_set()
            (1, 2, 3, 4)
            sage: CartanType(['A', 7, 2]).index_set()
            (0, 1, 2, 3, 4)
            sage: CartanType(['A', 7, 2]).index_set()
            (0, 1, 2, 3, 4)
            sage: CartanType(['A', 6, 2]).index_set()
            (0, 1, 2, 3)
            sage: CartanType(['D', 6, 2]).index_set()
            (0, 1, 2, 3, 4, 5)
            sage: CartanType(['E', 6, 1]).index_set()
            (0, 1, 2, 3, 4, 5, 6)
            sage: CartanType(['E', 6, 2]).index_set()
            (0, 1, 2, 3, 4)
            sage: CartanType(['A', 2, 2]).index_set()
            (0, 1)
            sage: CartanType(['G', 2, 1]).index_set()
            (0, 1, 2)
            sage: CartanType(['F', 4, 1]).index_set()
            (0, 1, 2, 3, 4)
        """

    @abstract_method
    def coxeter_matrix(self):
        """
        Return the Coxeter matrix associated to ``self``.

        EXAMPLES::

            sage: CartanType(['A', 3]).is_affine()
            False
            sage: CartanType(['A', 3, 1]).is_affine()
            True
        """

    @abstract_method
    def coxeter_diagram(self):
        """
        Return the Coxeter diagram associated to ``self``.

        EXAMPLES::

            sage: CartanType(['A', 3]).is_affine()
            False
            sage: CartanType(['A', 3, 1]).is_affine()
            True
        """

    @abstract_method
    def is_finite(self):
        """
        Return whether ``self`` is finite.

        EXAMPLES::

            sage: from sage.combinat.root_system.cartan_type import CartanType_abstract
            sage: C = CartanType_abstract()
            sage: C.is_finite()
            Traceback (most recent call last):
            ...
            NotImplementedError: <abstract method is_finite at ...>

        ::

            sage: CartanType(['A',4]).is_finite()
            True
            sage: CartanType(['A',4, 1]).is_finite()
            False
        """

    @abstract_method
    def is_affine(self):
        """
        Return whether ``self`` is affine.

        EXAMPLES::

            sage: CartanType(['A', 3]).is_affine()
            False
            sage: CartanType(['A', 3, 1]).is_affine()
            True
        """

    def is_crystallographic(self):
        """
        Return whether ``self`` is crystallographic.

        This returns ``False`` by default. Derived class should override this
        appropriately.

        EXAMPLES::

            sage: [ [t, t.is_crystallographic() ] for t in CartanType.samples(finite=True) ]
            [[['A', 1], True], [['A', 5], True],
             [['B', 1], True], [['B', 5], True],
             [['C', 1], True], [['C', 5], True],
             [['D', 2], True], [['D', 3], True], [['D', 5], True],
             [['E', 6], True], [['E', 7], True], [['E', 8], True],
             [['F', 4], True], [['G', 2], True],
             [['I', 5], False], [['H', 3], False], [['H', 4], False]]
        """
        return False

    def is_simply_laced(self):
        """
        Return whether ``self`` is simply laced.

        This returns ``False`` by default. Derived class should override this
        appropriately.

        EXAMPLES::

            sage: [ [t, t.is_simply_laced() ] for t in CartanType.samples() ]
            [[['A', 1], True], [['A', 5], True],
             [['B', 1], True], [['B', 5], False],
             [['C', 1], True], [['C', 5], False],
             [['D', 2], True], [['D', 3], True], [['D', 5], True],
             [['E', 6], True], [['E', 7], True], [['E', 8], True],
             [['F', 4], False], [['G', 2], False], [['I', 5], False], [['H', 3], False], [['H', 4], False],
             [['A', 1, 1], False], [['A', 5, 1], True],
             [['B', 1, 1], False], [['B', 5, 1], False],
             [['C', 1, 1], False], [['C', 5, 1], False],
             [['D', 3, 1], True], [['D', 5, 1], True],
             [['E', 6, 1], True], [['E', 7, 1], True], [['E', 8, 1], True],
             [['F', 4, 1], False], [['G', 2, 1], False],
             [['BC', 1, 2], False], [['BC', 5, 2], False],
             [['B', 5, 1]^*, False], [['C', 4, 1]^*, False], [['F', 4, 1]^*, False], [['G', 2, 1]^*, False],
             [['BC', 1, 2]^*, False], [['BC', 5, 2]^*, False]]
        """
        return False

class CoxeterTypeFromCartanType(CoxeterType_abstract):
    """
    A Coxeter type associated to a Cartan type.
    """
    def __init__(self, cartan_type):
        """
        Initialize ``self``.
        """
        self._cartan_type = cartan_type

    def _repr_(self):
        """
        Return a string representation of ``self``.
        """
        return "Coxeter type of {}".format(self._cartan_type)

    def coxeter_matrix(self):
        """
        Return the Coxeter matrix associated to ``self``.

        EXAMPLES::
        """
        return self._cartan_type.coxeter_matrix()

    def coxeter_diagram(self):
        """
        Return the Coxeter digramh of ``self``.
        """
        return self._cartan_type.coxeter_diagram()

    def cartan_type(self):
        """
        Return the Cartan type used to construct ``self``.
        """
        return self._cartan_type

    def rank(self):
        """
        Return the rank of ``self``.
        """
        return self._cartan_type.rank()

    def index_set(self):
        """
        Return the index set of ``self``.
        """
        return self._cartan_type.index_set()

    def is_finite(self):
        """
        Return if ``self`` is a finite type.
        """
        return self._cartan_type.is_finite()

    def is_affine(self):
        """
        Return if ``self`` is an affine type.
        """
        return self._cartan_type.is_affine()

    def is_crystallographic(self):
        """
        Return if ``self`` is crystallographic.
        """
        return self._cartan_type.is_crystallographic()

    def is_simply_laced(self):
        """
        Return if ``self`` is simply-laced.
        """
        return self._cartan_type.is_simply_laced()

