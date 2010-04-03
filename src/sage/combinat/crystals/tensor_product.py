"""
Tensor Products of Crystals
"""
#*****************************************************************************
#       Copyright (C) 2007 Anne Schilling <anne at math.ucdavis.edu>
#                          Nicolas Thiery <nthiery at users.sf.net>
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

import operator
from sage.misc.latex import latex
from sage.misc.cachefunc import cached_method
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets, EnumeratedSets
from sage.structure.element import Element, parent
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.cartesian_product import CartesianProduct
from sage.combinat.combinat import CombinatorialObject
from sage.combinat.partition import Partition
from sage.combinat.tableau import Tableau
from sage.combinat.tableau import Tableau_class
from crystals import CrystalElement, ClassicalCrystal, Crystal
from letters import CrystalOfLetters
from sage.misc.flatten import flatten
from sage.rings.integer import Integer

##############################################################################
# Support classes
##############################################################################

from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation

class TestParent(UniqueRepresentation, Parent):
    def _repr_(self):
        return "A parent for tests"

class ImmutableListWithParent(CombinatorialObject, Element):
    r"""
    A class for lists having a parent

    Specification: any subclass C should implement __init__ which
    accepts the following form C(parent, list = list)

    EXAMPLES: We create an immutable list whose parent is the class
    list::

        sage: from sage.combinat.crystals.tensor_product import ImmutableListWithParent, TestParent
        sage: l = ImmutableListWithParent(TestParent(), [1,2,3])
        sage: l._list
        [1, 2, 3]
        sage: l.parent()
        A parent for tests
        sage: l.sibling([2,1]) == ImmutableListWithParent(TestParent(), [2,1])
        True
        sage: l.reversed()
        [3, 2, 1]
        sage: l.set_index(1,4)
        [1, 4, 3]
    """

    def __init__(self, parent, list):
        """
        EXAMPLES::

            sage: from sage.combinat.crystals.tensor_product import ImmutableListWithParent, TestParent
            sage: l = ImmutableListWithParent(TestParent(), [1,2,3])
            sage: l.parent()
            A parent for tests
            sage: parent(l)
            A parent for tests
            sage: TestSuite(l).run(skip = "_test_category")
        """
        Element.__init__(self, parent)
        CombinatorialObject.__init__(self, list)

    def _repr_(self):
        """
        EXAMPLES::

            sage: from sage.combinat.crystals.tensor_product import ImmutableListWithParent, TestParent
            sage: l = ImmutableListWithParent(TestParent(), [1,2,3])
            sage: l.__repr__()
            '[1, 2, 3]'
        """
        return "%s"%self._list

    def __eq__(self, other):
        """
        EXAMPLES::

            sage: from sage.combinat.crystals.tensor_product import ImmutableListWithParent, TestParent
            sage: l = ImmutableListWithParent(TestParent(), [1,2,3])
            sage: m = ImmutableListWithParent(ZZ, [1,2,3])
            sage: n = ImmutableListWithParent(ZZ, [2,3,4])
            sage: l == l
            True
            sage: l == m
            False
            sage: m == n
            False
        """
        return self.__class__ is other.__class__ and \
               self.parent() == other.parent() and \
               self._list == other._list

    def __ne__(self, other):
        return not self.__eq__(other)

    # Should go in Element? Sets.ElementMethods?
    # How to define them conditionally, only of __lt__ is defined?
    def __le__(self, other):
        if self == other:
            return True
        else:
            return self.__le__(other)

    def __gt__(self, other):
        if parent(self) is not parent(other):
            return NotImplemented
        return other.__lt__(self)

    def __ge__(self, other):
        if self == other:
            return True
        else:
            return self.__le__(other)

    def sibling(self, l):
        """
        Returns an ImmutableListWithParent object whose list is l and whose
        parent is the same as self's parent.

        Note that the implementation of this function makes an assumption
        about the constructor for subclasses.

        EXAMPLES::

            sage: from sage.combinat.crystals.tensor_product import ImmutableListWithParent, TestParent
            sage: l = ImmutableListWithParent(TestParent(), [1,2,3])
            sage: m = l.sibling([2,3,4]); m
            [2, 3, 4]
            sage: m.parent()
            A parent for tests
        """
        return self.__class__(self.parent(), list=l)

    def reversed(self):
        """
        Returns the sibling of self which is obtained by reversing the
        elements of self.

        EXAMPLES::

            sage: from sage.combinat.crystals.tensor_product import ImmutableListWithParent, TestParent
            sage: l = ImmutableListWithParent(TestParent(), [1,2,3])
            sage: l.reversed()
            [3, 2, 1]
        """
        return self.sibling([ i for i in reversed(self._list)])

    def set_index(self, k, value):
        """
        Returns the sibling of self obtained by setting the
        `k^{th}` entry of self to value.

        EXAMPLES::

            sage: from sage.combinat.crystals.tensor_product import ImmutableListWithParent, TestParent
            sage: l = ImmutableListWithParent(TestParent(), [1,2,3])
            sage: l.set_index(0,2)
            [2, 2, 3]
            sage: l.set_index(1,4)
            [1, 4, 3]
            sage: _.parent()
            A parent for tests
        """
        l = [i for i in self._list]
        l[k] = value
        return self.sibling(l)

def TensorProductOfCrystals(*crystals, **options):
    r"""
    Tensor product of crystals.

    Given two crystals `B` and `B'` of the same type,
    one can form the tensor product `B \otimes B'`. As a set
    `B \otimes B'` is the Cartesian product
    `B \times B'`. The crystal operators `f_i` and
    `e_i` act on `b \otimes b' \in B \otimes B'` as
    follows:

    .. math::

      f_i(b \otimes b') = \begin{cases}
      f_i(b) \otimes b' & \text{if $\varepsilon_i(b) \ge \varphi_i(b')$}\\
      b \otimes f_i(b') & \text{otherwise}
      \end{cases}

    and

    .. math::

      e_i(b \otimes b') = \begin{cases}
      b \otimes e_i(b') & \text{if $\varepsilon_i(b) \le \varphi_i(b')$}\\
      e_i(b) \otimes b' & \text{otherwise.}
      \end{cases}

    Note that this is the opposite of Kashiwara's convention for tensor
    products of crystals.

    EXAMPLES:

    We construct the type `A_2`-crystal generated by `2 \otimes 1 \otimes 1`::

        sage: C = CrystalOfLetters(['A',2])
        sage: T = TensorProductOfCrystals(C,C,C,generators=[[C(2),C(1),C(1)]])

    It has `8` elements::

        sage: T.list()
        [[2, 1, 1], [2, 1, 2], [2, 1, 3], [3, 1, 3], [3, 2, 3], [3, 1, 1], [3, 1, 2], [3, 2, 2]]

    One can also check the Cartan type of the crystal::

        sage: T.cartan_type()
        ['A', 2]

    Other examples include crystals of tableaux (which internally are represented as tensor products
    obtained by reading the tableaux columnwise)::

        sage: C = CrystalOfTableaux(['A',3], shape=[1,1,0])
        sage: D = CrystalOfTableaux(['A',3], shape=[1,0,0])
        sage: T = TensorProductOfCrystals(C,D, generators=[[C(rows=[[1], [2]]), D(rows=[[1]])], [C(rows=[[2], [3]]), D(rows=[[1]])]])
        sage: T.cardinality()
        24
        sage: T.check()
        True
        sage: T.module_generators
        [[[[1], [2]], [[1]]], [[[2], [3]], [[1]]]]
        sage: [x.weight() for x in T.module_generators]
        [(2, 1, 0, 0), (1, 1, 1, 0)]

    If no module generators are specified, we obtain the full tensor
    product::

        sage: C=CrystalOfLetters(['A',2])
        sage: T=TensorProductOfCrystals(C,C)
        sage: T.list()
        [[1, 1], [1, 2], [1, 3], [2, 1], [2, 2], [2, 3], [3, 1], [3, 2], [3, 3]]
        sage: T.cardinality()
        9

    For a tensor product of crystals without module generators, the
    default implementation of module_generators contains all elements
    in the tensor product of the crystals. If there is a subset of
    elements in the tensor product that still generates the crystal,
    this needs to be implemented for the specific crystal separately::

        sage: T.module_generators.list()
        [[1, 1], [1, 2], [1, 3], [2, 1], [2, 2], [2, 3], [3, 1], [3, 2], [3, 3]]

    For classical highest weight crystals, it is also possible to list
    all highest weight elements::

        sage: C = CrystalOfLetters(['A',2])
        sage: T = TensorProductOfCrystals(C,C,C,generators=[[C(2),C(1),C(1)],[C(1),C(2),C(1)]])
        sage: T.highest_weight_vectors()
        [[2, 1, 1], [1, 2, 1]]
    """
    crystals = tuple(crystals)
    if "cartan_type" in options:
        cartan_type = CartanType(options["cartan_type"])
    else:
        if len(crystals) == 0:
            raise ValueError, "you need to specify the Cartan type if the tensor product list is empty"
        else:
            cartan_type = crystals[0].cartan_type()

    if "generators" in options:
        generators = tuple(tuple(x) if isinstance(x, list) else x for x in options["generators"])
        if all(isinstance(crystal,ClassicalCrystal) for crystal in crystals):
            return TensorProductOfClassicalCrystalsWithGenerators(crystals, generators=generators, cartan_type = cartan_type)
        else:
            return TensorProductOfCrystalsWithGenerators(crystals, generators=generators, cartan_type = cartan_type)
    else:
        if all(isinstance(crystal,ClassicalCrystal) for crystal in crystals):
            return FullTensorProductOfClassicalCrystals(crystals, cartan_type = cartan_type)
        else:
            return FullTensorProductOfCrystals(crystals, cartan_type = cartan_type)


class CrystalOfWords(Crystal):
    """
    Auxiliary class to provide a call method to create tensor product elements.
    This class is shared with several tensor product classes and is also used in CrystalOfTableaux to allow
    tableaux of different tensor product structures in column-reading (and hence different shapes)
    to be considered elements in the same crystal.
    """
    def _element_constructor_(self, *crystalElements):
        """
        EXAMPLES::

            sage: C = CrystalOfLetters(['A',2])
            sage: T = TensorProductOfCrystals(C,C)
            sage: T(1,1)
            [1, 1]
            sage: _.parent()
            Full tensor product of the crystals [The crystal of letters for type ['A', 2], The crystal of letters for type ['A', 2]]
            sage: T = TensorProductOfCrystals(C,C,C,generators=[[C(2),C(1),C(1)]])
            sage: T(C(2), C(1), C(1))
            [2, 1, 1]
        """
        return self.element_class(self, list(crystalElements))


class TensorProductOfCrystalsWithGenerators(CrystalOfWords):

    def __init__(self, crystals, generators, cartan_type):
        """
        EXAMPLES::

            sage: C = CrystalOfLetters(['A',2])
            sage: T = TensorProductOfCrystals(C,C,C,generators=[[C(2),C(1),C(1)]])
            sage: TestSuite(T).run()
        """
        assert isinstance(crystals, tuple)
        assert isinstance(generators, tuple)
        if all(crystal in FiniteEnumeratedSets() for crystal in crystals):
            category = FiniteEnumeratedSets()
        else:
            category = EnumeratedSets()
        super(TensorProductOfCrystalsWithGenerators, self).__init__(category = category)
        self.rename("The tensor product of the crystals %s"%(crystals,))
        self.crystals = crystals
        self._cartan_type = cartan_type
        self.module_generators = [ self(*x) for x in generators ]

class TensorProductOfClassicalCrystalsWithGenerators(TensorProductOfCrystalsWithGenerators, ClassicalCrystal):
    pass

class FullTensorProductOfCrystals(CrystalOfWords):
    def __init__(self, crystals, **options):
        """
        TESTS::

            sage: from sage.combinat.crystals.tensor_product import FullTensorProductOfCrystals
            sage: C = CrystalOfLetters(['A',2])
            sage: T = TensorProductOfCrystals(C,C)
            sage: isinstance(T, FullTensorProductOfCrystals)
            True
            sage: TestSuite(T).run()
        """
        crystals = list(crystals)
        if all(crystal in FiniteEnumeratedSets() for crystal in crystals):
            category = FiniteEnumeratedSets()
        else:
            category = EnumeratedSets()
        super(FullTensorProductOfCrystals, self).__init__(category = category)
        self.rename("Full tensor product of the crystals %s"%(crystals,))
        self.crystals = crystals
        if options.has_key('cartan_type'):
            self._cartan_type = CartanType(options['cartan_type'])
        else:
            if len(crystals) == 0:
                raise ValueError, "you need to specify the Cartan type if the tensor product list is empty"
            else:
                self._cartan_type = crystals[0].cartan_type()
        self.cartesian_product = CartesianProduct(*self.crystals)
        self.module_generators = self

    def __iter__(self):
        """
        EXAMPLES::

            sage: C = CrystalOfLetters(['A',2])
            sage: T = TensorProductOfCrystals(C,C)
            sage: list(T)
            [[1, 1], [1, 2], [1, 3], [2, 1], [2, 2], [2, 3], [3, 1], [3, 2], [3, 3]]
            sage: _[0].parent()
            Full tensor product of the crystals [The crystal of letters for type ['A', 2], The crystal of letters for type ['A', 2]]
        """
        for x in self.cartesian_product:
            yield self(*x)

#    list = CombinatorialClass._CombinatorialClass__list_from_iterator

    def cardinality(self):
        """
        EXAMPLES::

            sage: C = CrystalOfLetters(['A',2])
            sage: T = TensorProductOfCrystals(C,C)
            sage: T.cardinality()
            9
        """
        return self.cartesian_product.cardinality()


class FullTensorProductOfClassicalCrystals(FullTensorProductOfCrystals, ClassicalCrystal):
    pass

class TensorProductOfCrystalsElement(ImmutableListWithParent, CrystalElement):
    r"""
    A class for elements of tensor products of crystals
    """

    def __lt__(self, other):
        """
        Non elements of the crystal are incomparable with elements of the crystal
        (or should it return NotImplemented?)

        Comparison of two elements of this crystal:
         - different length: incomparable
         - otherwise lexicographicaly, considering self[i] and other[i]
           as incomparable if self[i] < other[i] returns NotImplemented
        """
        if parent(self) is not parent(other):
            return False
        if len(self) != len(other):
            return False
        for i in range(len(self)):
            if (self[i] < other[i]) == True:
                return True
            if (self[i] > other[i]) == True:
                return False
        return False

    def e(self, i):
        """
        Returns the action of `e_i` on self.

        EXAMPLES::

            sage: C = CrystalOfLetters(['A',5])
            sage: T = TensorProductOfCrystals(C,C)
            sage: T(C(1),C(2)).e(1) == T(C(1),C(1))
            True
            sage: T(C(2),C(1)).e(1) == None
            True
            sage: T(C(2),C(2)).e(1) == T(C(1),C(2))
            True
        """
        assert i in self.index_set()
        position = self.positions_of_unmatched_plus(i)
        if position == []:
            return None
        k = position[0]
        return self.set_index(k, self[k].e(i))

    def weight(self):
        """
        Returns the weight of self.

        EXAMPLES::

            sage: C = CrystalOfLetters(['A',3])
            sage: T = TensorProductOfCrystals(C,C)
            sage: T(C(1),C(2)).weight()
            (1, 1, 0, 0)
            sage: T=CrystalOfTableaux(['D',4],shape=[])
            sage: T.list()[0].weight()
            (0, 0, 0, 0)
        """
        return sum((self[j].weight() for j in range(len(self))), self.parent().weight_lattice_realization().zero())

    def f(self, i):
        """
        Returns the action of `f_i` on self.

        EXAMPLES::

            sage: C = CrystalOfLetters(['A',5])
            sage: T = TensorProductOfCrystals(C,C)
            sage: T(C(1),C(1)).f(1)
            [1, 2]
            sage: T(C(1),C(2)).f(1)
            [2, 2]
            sage: T(C(2),C(1)).f(1) is None
            True
        """
        assert i in self.index_set()
        position = self.positions_of_unmatched_minus(i)
        if position == []:
            return None
        k = position[len(position)-1]
        return self.set_index(k, self[k].f(i))

    def phi(self, i):
        """
        EXAMPLES::

            sage: C = CrystalOfLetters(['A',5])
            sage: T = TensorProductOfCrystals(C,C)
            sage: T(C(1),C(1)).phi(1)
            2
            sage: T(C(1),C(2)).phi(1)
            1
            sage: T(C(2),C(1)).phi(1)
            0
        """
        self = self.reversed()
        height = 0
        for j in range(len(self)):
            plus = self[j].epsilon(i)
            minus = self[j].phi(i)
            if height-plus < 0:
                height = minus
            else:
                height = height - plus + minus
        return height

    def epsilon(self, i):
        """
        EXAMPLES::

            sage: C = CrystalOfLetters(['A',5])
            sage: T = TensorProductOfCrystals(C,C)
            sage: T(C(1),C(1)).epsilon(1)
            0
            sage: T(C(1),C(2)).epsilon(1)
            1
            sage: T(C(2),C(1)).epsilon(1)
            0
        """
        height = 0
        for j in range(len(self)):
            minus = self[j].phi(i)
            plus = self[j].epsilon(i)
            if height-minus < 0:
                height = plus
            else:
                height = height - minus + plus
        return height

    def positions_of_unmatched_minus(self, i, dual=False, reverse=False):
        """
        EXAMPLES::

            sage: C = CrystalOfLetters(['A',5])
            sage: T = TensorProductOfCrystals(C,C)
            sage: T(C(2),C(1)).positions_of_unmatched_minus(1)
            []
            sage: T(C(1),C(2)).positions_of_unmatched_minus(1)
            [0]
        """
        unmatched_plus = []
        height = 0
        if reverse == True:
            self = self.reversed()
        if dual == False:
            for j in range(len(self)):
                minus = self[j].phi(i)
                plus = self[j].epsilon(i)
                if height-minus < 0:
                    unmatched_plus.append(j)
                    height = plus
                else:
                    height = height - minus + plus
        else:
            for j in range(len(self)):
                plus = self[j].epsilon(i)
                minus = self[j].phi(i)
                if height-plus < 0:
                    unmatched_plus.append(j)
                    height = minus
                else:
                    height = height - plus + minus
        return unmatched_plus

    def positions_of_unmatched_plus(self, i):
        """
        EXAMPLES::

            sage: C = CrystalOfLetters(['A',5])
            sage: T = TensorProductOfCrystals(C,C)
            sage: T(C(2),C(1)).positions_of_unmatched_plus(1)
            []
            sage: T(C(1),C(2)).positions_of_unmatched_plus(1)
            [1]
        """
        l = self.positions_of_unmatched_minus(i, dual=True, reverse=True)
        l.reverse()
        return [len(self)-1-l[j] for j in range(len(l))]

    def _latex_(self):
        """
        Returns latex code for self.

        EXAMPLES::

            sage: C = CrystalOfLetters(["A",2])
            sage: D = CrystalOfTableaux(["A",2], shape=[2])
            sage: E = TensorProductOfCrystals(C,D)
            sage: E.module_generators[0]._latex_()
            '1\\otimes{\\def\\lr#1{\\multicolumn{1}{|@{\\hspace{.6ex}}c@{\\hspace{.6ex}}|}{\\raisebox{-.3ex}{$#1$}}}\n\\raisebox{-.6ex}{$\\begin{array}[b]{cc}\n\\cline{1-1}\\cline{2-2}\n\\lr{1}&\\lr{1}\\\\\n\\cline{1-1}\\cline{2-2}\n\\end{array}$}\n}'
        """
        return '\otimes'.join(latex(c) for c in self)

CrystalOfWords.Element = TensorProductOfCrystalsElement

class CrystalOfTableaux(CrystalOfWords, ClassicalCrystal):
    r"""
    Crystals of tableaux. Input: a Cartan Type type and "shape", a
    partition of length <= type[1]. Produces a classical crystal with
    the given Cartan Type and highest weight corresponding to the given
    shape.

    If the type is ['D',r] then the shape is permitted to have a
    negative value in the r-th position. Thus if
    shape=`[s_1,s_2,...,s_r]` then s_r may be negative but
    in any case `s1 \ge s2 \ge ... s_{r-1} \ge |s_r|`. This
    crystal is related to `[s_1,\cdots,|s_r|]` by the outer
    automorphism of SO(2r).

    Note that crystals of tableaux are constructed using an embedding
    into tensor products following Kashiwara and Nakashima [Kashiwara,
    Masaki; Nakashima, Toshiki,
    *Crystal graphs for representations of the q-analogue of classical Lie algebras*,
    J. Algebra 165 (1994), no. 2, 295-345.]. Sage's tensor product rule
    for crystals differs from that of Kashiwara and Nakashima by
    reversing the order of the tensor factors. Sage produces the same
    crystals of tableaux as Kashiwara and Nakashima. With Sage's
    convention, the tensor product of crystals is the same as the
    monoid operation on tableaux and hence the plactic monoid.

    EXAMPLES:

    We create the crystal of tableaux for type `A_2`, with
    highest weight given by the partition [2,1,1]::

        sage: Tab = CrystalOfTableaux(['A',3], shape = [2,1,1])

    Here is the list of its elements::

        sage: Tab.list()
        [[[1, 1], [2], [3]], [[1, 2], [2], [3]], [[1, 3], [2], [3]],
         [[1, 4], [2], [3]], [[1, 4], [2], [4]], [[1, 4], [3], [4]],
         [[2, 4], [3], [4]], [[1, 1], [2], [4]], [[1, 2], [2], [4]],
         [[1, 3], [2], [4]], [[1, 3], [3], [4]], [[2, 3], [3], [4]],
         [[1, 1], [3], [4]], [[1, 2], [3], [4]], [[2, 2], [3], [4]]]

    One can get (currently) crude plotting via::

        #        sage: Tab.plot()  # random

    One can get instead get a LaTeX drawing ready to be copy-pasted
    into a LaTeX file::

        #        sage: Tab.latex() # random

    See sage.combinat.crystals.crystals? for general help on using
    crystals

    Internally, a tableau of a given Cartan type is represented as a
    tensor product of letters of the same type. The order in which the
    tensor factors appear is by reading the columns of the tableaux
    left to right, top to bottom (in French notation). As an example::

        sage: T = CrystalOfTableaux(['A',2], shape = [3,2])
        sage: T.module_generators[0]
        [[1, 1, 1], [2, 2]]
        sage: T.module_generators[0]._list
        [2, 1, 2, 1, 1]

    To create a tableau, one can use::

        sage: Tab = CrystalOfTableaux(['A',3], shape = [2,2])
        sage: Tab(rows=[[1,2],[3,4]])
        [[1, 2], [3, 4]]
        sage: Tab(columns=[[3,1],[4,2]])
        [[1, 2], [3, 4]]

    FIXME: do we want to specify the columns increasingly or
    decreasingly That is, should this be Tab(columns = [[1,3],[2,4]])

    TODO: make this fully consistent with Tableau!

    TESTS:

    Base cases::

        sage: T = CrystalOfTableaux(['A',2], shape = [])
        sage: T.list()
        [[]]
        sage: T = CrystalOfTableaux(['C',2], shape = [1])
        sage: T.check()
        True
        sage: T.list()
        [[[1]], [[2]], [[-2]], [[-1]]]
        sage: T = CrystalOfTableaux(['A',2], shapes = [[],[1],[2]])
        sage: T.list()
        [[], [[1]], [[2]], [[3]], [[1, 1]], [[1, 2]], [[2, 2]], [[1, 3]], [[2, 3]], [[3, 3]]]
        sage: T.module_generators
        ([], [[1]], [[1, 1]])
        sage: T = CrystalOfTableaux(['B',2],shape=[3])
        sage: T(rows=[[1,1,0]])
        [[1, 1, 0]]


    Input tests::

        sage: Tab = CrystalOfTableaux(['A',3], shape = [2,2])
        sage: C = Tab.letters
        sage: Tab(rows    = [[1,2],[3,4]])._list == [C(3),C(1),C(4),C(2)]
        True
        sage: Tab(columns = [[3,1],[4,2]])._list == [C(3),C(1),C(4),C(2)]
        True

    And for compatibility with TensorProductOfCrystal we should also
    allow as input the internal list / sequence of elements::

        sage: Tab(list    = [3,1,4,2])._list     == [C(3),C(1),C(4),C(2)]
        True
        sage: Tab(3,1,4,2)._list                 == [C(3),C(1),C(4),C(2)]
        True

    Type D, illustrating that the last parameter in the shape can be
    negative::

        sage: C = CrystalOfTableaux(['D',4],shape=[1,1,1,-1])
        sage: C.cardinality()
        35
        sage: C.check()
        True

    The next example checks whether a given tableau is in fact a valid type C tableau or not:

        sage: T = CrystalOfTableaux(['C',3], shape = [2,2,2])
        sage: t = T(rows=[[1,3],[2,-3],[3,-1]])
        sage: t in T.list()
        True
        sage: t = T(rows=[[2,3],[3,-3],[-3,-2]])
        sage: t in T.list()
        False
    """

    @staticmethod
    def __classcall__(cls, cartan_type, shapes = None, shape = None):
        """
        Normalizes the input arguments to ensure unique representation

        EXAMPLES::

            sage: T1 = CrystalOfTableaux(CartanType(['A',3]), shape  = [2,2])
            sage: T2 = CrystalOfTableaux(['A',3],             shape  = (2,2))
            sage: T3 = CrystalOfTableaux(['A',3],             shapes = ([2,2],))
            sage: T2 is T1, T3 is T1
            (True, True)
        """
        cartan_type = CartanType(cartan_type)
        # standardize shape/shapes input into a tuple of tuples
        assert operator.xor(shape is not None, shapes is not None)
        if shape is not None:
            shapes = (shape,)
        shapes = tuple( tuple(shape) for shape in shapes )
        return super(CrystalOfTableaux, cls).__classcall__(cls, cartan_type, shapes)


    def __init__(self, cartan_type, shapes):
        """
        Construct the crystal of all tableaux of the given shapes

        INPUT:
         - ``cartan_type`` - (data coercible into) a cartan type
         - ``shapes``      - a list (or iterable) of shapes
         - ``shape` `      - a shape

        shapes themselves are lists (or iterable) of integers

        EXAMPLES::

            sage: T = CrystalOfTableaux(['A',3], shape = [2,2])
            sage: TestSuite(T).run()
        """
        super(CrystalOfTableaux, self).__init__(category = FiniteEnumeratedSets())
        self.letters = CrystalOfLetters(cartan_type)
        self.module_generators = tuple(self.module_generator(la) for la in shapes)
        self.rename("The crystal of tableaux of type %s and shape(s) %s"%(cartan_type, list(list(shape) for shape in shapes)))

    def cartan_type(self):
        """
        Returns the Cartan type of the associated crystal

        EXAMPLES::
            sage: T = CrystalOfTableaux(['A',3], shape = [2,2])
            sage: T.cartan_type()
            ['A', 3]
        """
        return self.letters.cartan_type()

    def module_generator(self, shape):
        """
        This yields the module generator (or highest weight element) of a classical
        crystal of given shape. The module generator is the unique tableau with equal
        shape and content.

        EXAMPLE:
            sage: T = CrystalOfTableaux(['D',3], shape = [1,1])
            sage: T.module_generator([1,1])
            [[1], [2]]
        """
        type = self.cartan_type()
        if type[0] == 'D' and len(shape) == type[1] and shape[type[1]-1] < 0:
            invert = True
            shape = shape[:-1]+(-shape[type[1]-1],)
        else:
            invert = False
        p = Partition(shape).conjugate()
        # The column canonical tableau, read by columns
        module_generator = flatten([[p[j]-i for i in range(p[j])] for j in range(len(p))])
        if invert:
            for i in range(type[1]):
                if module_generator[i] == type[1]:
                    module_generator[i] = -type[1]
        return self(list=[self.letters(x) for x in module_generator])

    def _element_constructor_(self, *args, **options):
        """
        Returns a CrystalOfTableauxElement

        EXAMPLES::

            sage: T = CrystalOfTableaux(['A',3], shape = [2,2])
            sage: T(rows=[[1,2],[3,4]])
            [[1, 2], [3, 4]]
            sage: T(columns=[[3,1],[4,2]])
            [[1, 2], [3, 4]]
        """
        return self.element_class(self, *args, **options)

class CrystalOfTableauxElement(TensorProductOfCrystalsElement):
    def __init__(self, parent, *args, **options):
        """
        There are several ways to input tableaux, by rows,
        by columns, as the list of column elements, or as a sequence of numbers
        in column reading.

        EXAMPLES::

            sage: T = CrystalOfTableaux(['A',3], shape = [2,2])
            sage: t = T(rows=[[1,2],[3,4]])
            sage: t
            [[1, 2], [3, 4]]
            sage: TestSuite(t).run()

            sage: t = T(columns=[[3,1],[4,2]])
            sage: t
            [[1, 2], [3, 4]]
            sage: TestSuite(t).run()

            sage: t = T(list=[3,1,4,2])
            sage: t
            [[1, 2], [3, 4]]

            sage: t = T(3,1,4,2)
            sage: t
            [[1, 2], [3, 4]]

        Currently inputting the empty tableau as an empty sequence is broken due to a bug in
        the generic __call__ method (see trac ticket #8648)

        EXAMPLES::

            sage: T = CrystalOfTableaux(['A',3], shape=[])
            sage: t = T()
            sage: t._list
            [0]
        """
        if len(args) == 1:
            if isinstance(args[0], Tableau_class):
                options['rows'] = args[0]
        if options.has_key('list'):
            list = options['list']
        elif options.has_key('rows'):
            rows=options['rows']
#            list=Tableau(rows).to_word_by_column()
            rows=Tableau(rows).conjugate()
            list=[]
            for col in rows:
                col.reverse()
                list+=col
        elif options.has_key('columns'):
            columns=options['columns']
            list=[]
            for col in columns:
                list+=col
        else:
            list = [i for i in args]
        if list != [] and type(list[0]) is Integer:
            list=[parent.letters(x) for x in list]
        TensorProductOfCrystalsElement.__init__(self, parent, list)

    def __repr__(self):
        """
        EXAMPLES::

            sage: T = CrystalOfTableaux(['A',3], shape = [2,2])
            sage: t = T(rows=[[1,2],[3,4]])
            sage: t.__repr__()
            '[[1, 2], [3, 4]]'
        """
        return repr(self.to_tableau())

    def pp(self):
        """
        EXAMPLES::

            sage: T = CrystalOfTableaux(['A',3], shape = [2,2])
            sage: t = T(rows=[[1,2],[3,4]])
            sage: t.pp()
            1  2
            3  4
        """
        return self.to_tableau().pp()

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: T = CrystalOfTableaux(['A',3], shape = [4,2])
            sage: t = T(rows=[[1,1,2,3],[2,3]])
            sage: latex(t) # indirect doctest
            {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
            \raisebox{-.6ex}{$\begin{array}[b]{cccc}
            \cline{1-1}\cline{2-2}\cline{3-3}\cline{4-4}
            \lr{1}&\lr{1}&\lr{2}&\lr{3}\\
            \cline{1-1}\cline{2-2}\cline{3-3}\cline{4-4}
            \lr{2}&\lr{3}\\
            \cline{1-1}\cline{2-2}
            \end{array}$}
            }
        """
        return latex(self.to_tableau())

    @cached_method
    def to_tableau(self):
        """
        Returns the Tableau object corresponding to self.

        EXAMPLES::

            sage: T = CrystalOfTableaux(['A',3], shape = [2,2])
            sage: t = T(rows=[[1,2],[3,4]]).to_tableau(); t
            [[1, 2], [3, 4]]
            sage: type(t)
            <class 'sage.combinat.tableau.Tableau_class'>
            sage: type(t[0][0])
            <type 'sage.rings.integer.Integer'>
            sage: T = CrystalOfTableaux(['D',3], shape = [1,1])
            sage: t=T(rows=[[-3],[3]]).to_tableau(); t
            [[-3], [3]]
            sage: t=T(rows=[[3],[-3]]).to_tableau(); t
            [[3], [-3]]
            sage: T = CrystalOfTableaux(['B',2], shape = [1,1])
            sage: t = T(rows=[[0],[0]]).to_tableau(); t
            [[0], [0]]
        """
        if self._list == []:
            return Tableau([])
        tab = [ [self[0].value] ]
        for i in range(1,len(self)):
            if self[i-1] < self[i] or (self[i-1].value != 0 and self[i-1] == self[i]):
                tab.append([self[i].value])
            else:
                l = len(tab)-1
                tab[l].append(self[i].value)
        for x in tab:
            x.reverse()
        return Tableau(tab).conjugate()

CrystalOfTableaux.Element = CrystalOfTableauxElement
