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

from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.combinat.root_system.cartan_type import CartanType
from sage.matrix.all import MatrixSpace
from sage.rings.all import ZZ
from sage.rings.universal_cyclotomic_field.universal_cyclotomic_field import UniversalCyclotomicField
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.sage_object import SageObject
from sage.combinat.root_system.cartan_type import CartanType

class CoxeterType(object):
    """
    Abstract class for Coxeter types.
    """
    __metaclass__ = ClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, x):
        """
        Parse input ``x``.

        EXAMPLES::

            sage: CoxeterType(['A',3])
            Coxeter type of ['A', 3]
        """
        if isinstance(x, CoxeterType):
            return x

        try:
            return CoxeterTypeFromCartanType(CartanType(x))
        except (ValueError, TypeError):
            pass

        raise NotImplementedError("Coxeter types not from Cartan types not yet implemented")

    @abstract_method
    def rank(self):
        """
        Return the rank of ``self``.

        This is the number of nodes of the associated Coxeter diagram.

        EXAMPLES::

            sage: CoxeterType(['A', 4]).rank()
            4
            sage: CoxeterType(['A', 7, 2]).rank()
            5
            sage: CoxeterType(['I', 8]).rank()
            2
        """

    @abstract_method
    def index_set(self):
        """
        Return the index set for ``self``.

        This is the list of the nodes of the associated Coxeter diagram.

        EXAMPLES::

            sage: CoxeterType(['A', 3, 1]).index_set()
            (0, 1, 2, 3)
            sage: CoxeterType(['D', 4]).index_set()
            (1, 2, 3, 4)
            sage: CoxeterType(['A', 7, 2]).index_set()
            (0, 1, 2, 3, 4)
            sage: CoxeterType(['A', 7, 2]).index_set()
            (0, 1, 2, 3, 4)
            sage: CoxeterType(['A', 6, 2]).index_set()
            (0, 1, 2, 3)
            sage: CoxeterType(['D', 6, 2]).index_set()
            (0, 1, 2, 3, 4, 5)
            sage: CoxeterType(['E', 6, 1]).index_set()
            (0, 1, 2, 3, 4, 5, 6)
            sage: CoxeterType(['E', 6, 2]).index_set()
            (0, 1, 2, 3, 4)
            sage: CoxeterType(['A', 2, 2]).index_set()
            (0, 1)
            sage: CoxeterType(['G', 2, 1]).index_set()
            (0, 1, 2)
            sage: CoxeterType(['F', 4, 1]).index_set()
            (0, 1, 2, 3, 4)
        """

    @abstract_method
    def coxeter_matrix(self):
        """
        Return the Coxeter matrix associated to ``self``.

        EXAMPLES::

            sage: CoxeterType(['A', 3]).coxeter_matrix()
            [1 3 2]
            [3 1 3]
            [2 3 1]
            sage: CoxeterType(['A', 3, 1]).coxeter_matrix()
            [1 3 2 3]
            [3 1 3 2]
            [2 3 1 3]
            [3 2 3 1]
        """

    @abstract_method
    def coxeter_diagram(self):
        """
        Return the Coxeter diagram associated to ``self``.

        EXAMPLES::

            sage: CoxeterType(['A', 3]).coxeter_diagram()
            Graph on 3 vertices
            sage: CoxeterType(['A', 3, 1]).coxeter_diagram()
            Graph on 4 vertices
        """

    @abstract_method
    def is_finite(self):
        """
        Return whether ``self`` is finite.

        EXAMPLES::

            sage: CoxeterType(['A',4]).is_finite()
            True
            sage: CoxeterType(['A',4, 1]).is_finite()
            False
        """

    @abstract_method
    def is_affine(self):
        """
        Return whether ``self`` is affine.

        EXAMPLES::

            sage: CoxeterType(['A', 3]).is_affine()
            False
            sage: CoxeterType(['A', 3, 1]).is_affine()
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

    @cached_method
    def bilinear_form(self, R=None):
        """
        Return the bilinear form over ``R`` associated to ``self``.

        INPUT:

        - ``R`` -- (default: universal cyclotomic field) a ring used to
          compute the bilinear form

        EXAMPLES::

            sage: CoxeterType(['A', 2, 1]).bilinear_form()
            [   1 -1/2 -1/2]
            [-1/2    1 -1/2]
            [-1/2 -1/2    1]
            sage: CoxeterType(['H', 3]).bilinear_form()
            [                      1                    -1/2                       0]
            [                   -1/2                       1 1/2*E(5)^2 + 1/2*E(5)^3]
            [                      0 1/2*E(5)^2 + 1/2*E(5)^3                       1]
            sage: C = CoxeterMatrix([[1,-1,-1],[-1,1,-1],[-1,-1,1]])
            sage: C.bilinear_form()
            [ 1 -1 -1]
            [-1  1 -1]
            [-1 -1  1]
        """
        if R is None:
            R = UniversalCyclotomicField()
        # Compute the matrix with entries `- \cos( \pi / m_{ij} )`.
        if R is UniversalCyclotomicField():
            val = lambda x: (R.gen(2*x) + ~R.gen(2*x)) / R(-2) if x != -1 else -R.one()
        else:
            from sage.functions.trig import cos
            from sage.symbolic.constants import pi
            val = lambda x: -R(cos(pi / x)) if x != -1 else -R.one()

        n = self.rank()
        MS = MatrixSpace(R, n, sparse=True)
        MC = MS._get_matrix_class()
        mat = self.coxeter_matrix()
        bilinear = MC(MS, entries={(i, j): val(mat[i, j])
                                   for i in range(n) for j in range(n)
                                   if mat[i, j] != 2},
                      coerce=True, copy=True)
        bilinear.set_immutable()
        return bilinear

class CoxeterTypeFromCartanType(SageObject, CoxeterType, UniqueRepresentation):
    """
    A Coxeter type associated to a Cartan type.
    """
    @staticmethod
    def __classcall_private__(cls, cartan_type):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: from sage.combinat.root_system.coxeter_type import CoxeterTypeFromCartanType
            sage: C1 = CoxeterTypeFromCartanType(['A',3])
            sage: C2 = CoxeterTypeFromCartanType(CartanType(['A',3]))
            sage: C1 is C2
            True
        """
        return super(CoxeterTypeFromCartanType, cls).__classcall__(cls,
                         CartanType(cartan_type))

    def __init__(self, cartan_type):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: C = CoxeterType(['A',3])
            sage: TestSuite(C).run()
            sage: C = CoxeterType(['H',4])
            sage: TestSuite(C).run()
            sage: C = CoxeterType(['C',3,1])
            sage: TestSuite(C).run()
        """
        self._cartan_type = cartan_type

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: CoxeterType(['A',3])
            Coxeter type of ['A', 3]
        """
        return "Coxeter type of {}".format(self._cartan_type)

    def coxeter_matrix(self):
        """
        Return the Coxeter matrix associated to ``self``.

        EXAMPLES::

            sage: C = CoxeterType(['H',3])
            sage: C.coxeter_matrix()
            [1 3 2]
            [3 1 5]
            [2 5 1]
        """
        return self._cartan_type.coxeter_matrix()

    def coxeter_diagram(self):
        """
        Return the Coxeter digramh of ``self``.

        EXAMPLES::

            sage: C = CoxeterType(['H',3])
            sage: C.coxeter_diagram().edges()
            [(1, 2, 3), (2, 3, 5)]
        """
        return self._cartan_type.coxeter_diagram()

    def cartan_type(self):
        """
        Return the Cartan type used to construct ``self``.

        EXAMPLES::

            sage: C = CoxeterType(['C',3])
            sage: C.cartan_type()
            ['C', 3]
        """
        return self._cartan_type

    def rank(self):
        """
        Return the rank of ``self``.

        EXAMPLES::

            sage: C = CoxeterType(['I', 16])
            sage: C.rank()
            2
        """
        return self._cartan_type.rank()

    def index_set(self):
        """
        Return the index set of ``self``.

        EXAMPLES::

            sage: C = CoxeterType(['A', 4])
            sage: C.index_set()
            (1, 2, 3, 4)
        """
        return self._cartan_type.index_set()

    def is_finite(self):
        """
        Return if ``self`` is a finite type.

        EXAMPLES::

            sage: C = CoxeterType(['E', 6])
            sage: C.is_finite()
            True
        """
        return self._cartan_type.is_finite()

    def is_affine(self):
        """
        Return if ``self`` is an affine type.

        EXAMPLES::

            sage: C = CoxeterType(['F', 4, 1])
            sage: C.is_affine()
            True
        """
        return self._cartan_type.is_affine()

    def is_crystallographic(self):
        """
        Return if ``self`` is crystallographic.

        EXAMPLES::

            sage: C = CoxeterType(['C', 3])
            sage: C.is_crystallographic()
            True

            sage: C = CoxeterType(['H', 3])
            sage: C.is_crystallographic()
            False
        """
        return self._cartan_type.is_crystallographic()

    def is_simply_laced(self):
        """
        Return if ``self`` is simply-laced.

        EXAMPLES::

            sage: C = CoxeterType(['A', 5])
            sage: C.is_simply_laced()
            True

            sage: C = CoxeterType(['B', 3])
            sage: C.is_simply_laced()
            False
        """
        return self._cartan_type.is_simply_laced()

