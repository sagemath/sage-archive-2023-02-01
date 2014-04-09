r"""
Highest weight crystals
"""

#*****************************************************************************
#       Copyright (C) 2009   Anne Schilling <anne at math.ucdavis.edu>
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
#****************************************************************************

from sage.categories.classical_crystals import ClassicalCrystals
from sage.structure.parent import Parent
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.crystals.letters import CrystalOfLetters
from sage.combinat.crystals.tensor_product import TensorProductOfCrystals, \
    TensorProductOfRegularCrystalsElement
from sage.combinat.crystals.littelmann_path import CrystalOfLSPaths


def HighestWeightCrystal(dominant_weight):
    r"""
    Returns an implementation of the highest weight crystal of highest weight `dominant_weight`.

    This is currently only implemented for crystals of type `E_6` and `E_7`.

    TODO: implement highest weight crystals for classical types `A_n`, `B_n`, `C_n`, `D_n` using tableaux.

    EXAMPLES::

        sage: C=CartanType(['E',6])
        sage: La=C.root_system().weight_lattice().fundamental_weights()
        sage: T = HighestWeightCrystal(La[1])
        sage: T.cardinality()
        27
        sage: T = HighestWeightCrystal(La[6])
        sage: T.cardinality()
        27
        sage: T = HighestWeightCrystal(La[2])
        sage: T.cardinality()
        78
        sage: T = HighestWeightCrystal(La[4])
        sage: T.cardinality()
        2925
        sage: T = HighestWeightCrystal(La[3])
        sage: T.cardinality()
        351
        sage: T = HighestWeightCrystal(La[5])
        sage: T.cardinality()
        351

        sage: C=CartanType(['E',7])
        sage: La=C.root_system().weight_lattice().fundamental_weights()
        sage: T = HighestWeightCrystal(La[1])
        sage: T.cardinality()
        133
        sage: T = HighestWeightCrystal(La[2])
        sage: T.cardinality()
        912
        sage: T = HighestWeightCrystal(La[3])
        sage: T.cardinality()
        8645
        sage: T = HighestWeightCrystal(La[4])
        sage: T.cardinality()
        365750
        sage: T = HighestWeightCrystal(La[5])
        sage: T.cardinality()
        27664
        sage: T = HighestWeightCrystal(La[6])
        sage: T.cardinality()
        1539
        sage: T = HighestWeightCrystal(La[7])
        sage: T.cardinality()
        56

        sage: C = CartanType(['C',2,1])
        sage: La = C.root_system().weight_lattice().fundamental_weights()
        sage: T = HighestWeightCrystal(La[1])
        sage: sorted(T.subcrystal(max_depth=3), key=str)
        [(-Lambda[0] + 3*Lambda[1] - Lambda[2] - delta,),
         (-Lambda[0] + Lambda[1] + Lambda[2] - delta,),
         (-Lambda[1] + 2*Lambda[2] - delta,),
         (2*Lambda[0] - Lambda[1],),
         (Lambda[0] + Lambda[1] - Lambda[2],),
         (Lambda[0] - Lambda[1] + Lambda[2],),
         (Lambda[1],)]
    """
    cartan_type = dominant_weight.parent().cartan_type()
    if cartan_type.is_finite() and cartan_type.type() in ['A','B','C','D']:
        raise NotImplementedError
    elif cartan_type == CartanType(['E',6]):
        return FiniteDimensionalHighestWeightCrystal_TypeE6(dominant_weight)
    elif cartan_type == CartanType(['E',7]):
        return FiniteDimensionalHighestWeightCrystal_TypeE7(dominant_weight)
    elif cartan_type.is_affine():
        return CrystalOfLSPaths(cartan_type,[dominant_weight[i] for i in cartan_type.index_set()])
    else:
        raise NotImplementedError

class FiniteDimensionalHighestWeightCrystal_TypeE(TensorProductOfCrystals):
    """
    Commonalities for all finite dimensional type E highest weight crystals

    Subclasses should setup an attribute column_crystal in their
    __init__ method before calling the __init__ method of this class.
    """

    def __init__(self, dominant_weight):
        """
        EXAMPLES::

            sage: C=CartanType(['E',6])
            sage: La=C.root_system().weight_lattice().fundamental_weights()
            sage: T = HighestWeightCrystal(2*La[2])
            sage: T.cartan_type()
            ['E', 6]
            sage: T.module_generators
            [[[(2, -1), (1,)], [(2, -1), (1,)]]]
            sage: T.cardinality()
            2430
            sage: T = HighestWeightCrystal(La[2])
            sage: T.cardinality()
            78
        """
        self._cartan_type = dominant_weight.parent().cartan_type()
        self._highest_weight = dominant_weight
        assert dominant_weight.is_dominant()
        self.rename("Finite dimensional highest weight crystal of type %s and highest weight %s"%(self._cartan_type, dominant_weight))
        Parent.__init__(self, category = ClassicalCrystals())
        self.module_generators = [self.module_generator()]

    Element = TensorProductOfRegularCrystalsElement

    def module_generator(self):
        """
        This yields the module generator (or highest weight element) of the classical
        crystal of given dominant weight in self.

        EXAMPLES::

            sage: C=CartanType(['E',6])
            sage: La=C.root_system().weight_lattice().fundamental_weights()
            sage: T = HighestWeightCrystal(La[2])
            sage: T.module_generator()
            [[(2, -1), (1,)]]
            sage: T = HighestWeightCrystal(0*La[2])
            sage: T.module_generator()
            []

            sage: C=CartanType(['E',7])
            sage: La=C.root_system().weight_lattice().fundamental_weights()
            sage: T = HighestWeightCrystal(La[1])
            sage: T.module_generator()
            [[(-7, 1), (7,)]]
        """
        dominant_weight = self._highest_weight
        tensor = sum(( [self.column_crystal[i]]*dominant_weight.coefficient(i) for i in dominant_weight.support()), [])
        return self._element_constructor_(*[B.module_generators[0] for B in tensor])

class FiniteDimensionalHighestWeightCrystal_TypeE6(FiniteDimensionalHighestWeightCrystal_TypeE):
    r"""
    Class of finite dimensional highest weight crystals of type `E_6`.

    EXAMPLES::

        sage: C=CartanType(['E',6])
        sage: La=C.root_system().weight_lattice().fundamental_weights()
        sage: T = HighestWeightCrystal(La[2]); T
        Finite dimensional highest weight crystal of type ['E', 6] and highest weight Lambda[2]
        sage: B1 = T.column_crystal[1]; B1
        The crystal of letters for type ['E', 6]
        sage: B6 = T.column_crystal[6]; B6
        The crystal of letters for type ['E', 6] (dual)
        sage: t = T(B6([-1]),B1([-1,3])); t
        [(-1,), (-1, 3)]
        sage: [t.epsilon(i) for i in T.index_set()]
        [2, 0, 0, 0, 0, 0]
        sage: [t.phi(i) for i in T.index_set()]
        [0, 0, 1, 0, 0, 0]
        sage: TestSuite(t).run()
    """

    def __init__(self, dominant_weight):
        """
        EXAMPLES::

            sage: C=CartanType(['E',6])
            sage: La=C.root_system().weight_lattice().fundamental_weights()
            sage: p2=2*La[2]
            sage: p1=La[2]
            sage: p0=0*La[2]
            sage: T = HighestWeightCrystal(0*La[2])
            sage: T.cardinality()
            1
            sage: T = HighestWeightCrystal(La[2])
            sage: T.cardinality()
            78
            sage: T = HighestWeightCrystal(2*La[2])
            sage: T.cardinality()
            2430
        """
        B1 = CrystalOfLetters(['E',6])
        B6 = CrystalOfLetters(['E',6], dual = True)
        self.column_crystal = {1 : B1, 6 : B6,
                               4 : TensorProductOfCrystals(B1,B1,B1,generators=[[B1([-3,4]),B1([-1,3]),B1([1])]]),
                               3 : TensorProductOfCrystals(B1,B1,generators=[[B1([-1,3]),B1([1])]]),
                               5 : TensorProductOfCrystals(B6,B6,generators=[[B6([5,-6]),B6([6])]]),
                               2 : TensorProductOfCrystals(B6,B1,generators=[[B6([2,-1]),B1([1])]])}
        FiniteDimensionalHighestWeightCrystal_TypeE.__init__(self, dominant_weight)


class FiniteDimensionalHighestWeightCrystal_TypeE7(FiniteDimensionalHighestWeightCrystal_TypeE):
    r"""
    Class of finite dimensional highest weight crystals of type `E_7`.

    EXAMPLES::

        sage: C=CartanType(['E',7])
        sage: La=C.root_system().weight_lattice().fundamental_weights()
        sage: T = HighestWeightCrystal(La[1])
        sage: T.cardinality()
        133
        sage: B7 = T.column_crystal[7]; B7
        The crystal of letters for type ['E', 7]
        sage: t = T(B7([-5, 6]), B7([-2, 3])); t
        [(-5, 6), (-2, 3)]
        sage: [t.epsilon(i) for i in T.index_set()]
        [0, 1, 0, 0, 1, 0, 0]
        sage: [t.phi(i) for i in T.index_set()]
        [0, 0, 1, 0, 0, 1, 0]
        sage: TestSuite(t).run()
    """

    def __init__(self, dominant_weight):
        """
        EXAMPLES::

            sage: C=CartanType(['E',7])
            sage: La=C.root_system().weight_lattice().fundamental_weights()
            sage: T = HighestWeightCrystal(0*La[1])
            sage: T.cardinality()
            1
            sage: T = HighestWeightCrystal(La[1])
            sage: T.cardinality()
            133
            sage: T = HighestWeightCrystal(2*La[1])
            sage: T.cardinality()
            7371
        """
        B = CrystalOfLetters(['E',7])
        self.column_crystal = {7 : B,
                               1 : TensorProductOfCrystals(B,B,generators=[[B([-7,1]),B([7])]]),
                               2 : TensorProductOfCrystals(B,B,B,generators=[[B([-1,2]),B([-7,1]),B([7])]]),
                               3 : TensorProductOfCrystals(B,B,B,B,generators=[[B([-2,3]),B([-1,2]),B([-7,1]),B([7])]]),
                               4 : TensorProductOfCrystals(B,B,B,B,generators=[[B([-5,4]),B([-6,5]),B([-7,6]),B([7])]]),
                               5 : TensorProductOfCrystals(B,B,B,generators=[[B([-6,5]),B([-7,6]),B([7])]]),
                               6 : TensorProductOfCrystals(B,B,generators=[[B([-7,6]),B([7])]])}
        FiniteDimensionalHighestWeightCrystal_TypeE.__init__(self, dominant_weight)
