r"""
Tensor Product of Kirillov-Reshetikhin Tableaux

A tensor product of :class:`KirillovReshetikhinTableaux` which are tableaux of
`r` rows and `s` columns which naturally arise in the bijection between rigged
configurations and tableaux and which are in bijection with the elements of the
Kirillov-Reshetikhin crystal `B^{r,s}`, see :class:`KirillovReshetikhinCrystal`.
They do not have to satisfy the semistandard row or column
restrictions. These tensor products are the result from the bijection from
rigged configurations [RigConBijection]_.

AUTHORS:

- Travis Scrimshaw (2010-09-26): Initial version

EXAMPLES:

Type `A_n^{(1)}` examples::

    sage: KRT = TensorProductOfKirillovReshetikhinTableaux(['A',3,1], [[3,1], [2,1]])
    sage: KRT
    Tensor product of Kirillov-Reshetikhin tableaux of type ['A', 3, 1] and factor(s) ((3, 1), (2, 1))
    sage: KRT.cardinality()
    24
    sage: KRT = TensorProductOfKirillovReshetikhinTableaux(['A',3,1], [[1,1], [2,1], [3,1]])
    sage: KRT
    Tensor product of Kirillov-Reshetikhin tableaux of type ['A', 3, 1] and factor(s) ((1, 1), (2, 1), (3, 1))
    sage: len(KRT.module_generators)
    5
    sage: KRT = TensorProductOfKirillovReshetikhinTableaux(['A',3,1], [[1,1], [2,1], [3,1]])
    sage: KRT.cardinality()
    96

Type `D_n^{(1)}` examples::

    sage: KRT = TensorProductOfKirillovReshetikhinTableaux(['D', 4, 1], [[1, 1], [2, 1], [1, 1]])
    sage: KRT
    Tensor product of Kirillov-Reshetikhin tableaux of type ['D', 4, 1] and factor(s) ((1, 1), (2, 1), (1, 1))
    sage: T = KRT(pathlist=[[1], [-2, 2], [1]])
    sage: T
    [[1]] (X) [[2], [-2]] (X) [[1]]
    sage: T2 = KRT(pathlist=[[1], [2, -2], [1]])
    sage: T2
    [[1]] (X) [[-2], [2]] (X) [[1]]
    sage: T == T2
    False
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

from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.structure.unique_representation import UniqueRepresentation

from sage.combinat.crystals.tensor_product import FullTensorProductOfRegularCrystals
from sage.combinat.crystals.letters import CrystalOfLetters
from sage.combinat.root_system.cartan_type import CartanType

from sage.combinat.rigged_configurations.tensor_product_kr_tableaux_element import TensorProductOfKirillovReshetikhinTableauxElement
from sage.combinat.rigged_configurations.kr_tableaux import KirillovReshetikhinTableaux, \
  KirillovReshetikhinTableauxElement

from sage.rings.integer import Integer

class HighestWeightTensorKRT(UniqueRepresentation):
    """
    Class so we do not have to build the module generators for
    :class:`TensorProductOfKirillovReshetikhinTableaux` at initialization.

    .. WARNING::

        This class is for internal use only!
    """
    def __init__(self, tp_krt):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: KRT = TensorProductOfKirillovReshetikhinTableaux(['D',4,1], [[2,2]])
            sage: from sage.combinat.rigged_configurations.tensor_product_kr_tableaux import HighestWeightTensorKRT
            sage: hw = HighestWeightTensorKRT(KRT)
            sage: hw2 = HighestWeightTensorKRT(KRT)
            sage: hw is hw2
            True
        """
        self.tp_krt = tp_krt
        self._cache = None

    def __getitem__(self, i):
        """
        Return the `i`-th highest weight element in the cache.

        TESTS::

            sage: KRT = TensorProductOfKirillovReshetikhinTableaux(['A',3,1], [[3,1], [2,1]])
            sage: KRT.module_generators[0]
            [[1], [2], [3]] (X) [[1], [2]]
        """
        if self._cache is None:
            self._cache = [x.to_tensor_product_of_Kirillov_Reshetikhin_tableaux()
                           for x in self.tp_krt.rigged_configurations().module_generators]
        return self._cache[i]

    def __iter__(self):
        """
        Iterate over the highest weight elements.

        EXAMPLES::

            sage: KRT = TensorProductOfKirillovReshetikhinTableaux(['D',4,1], [[2,1]])
            sage: from sage.combinat.rigged_configurations.tensor_product_kr_tableaux import HighestWeightTensorKRT
            sage: for x in HighestWeightTensorKRT(KRT): x
            ...
            [[1], [2]]
            [[1], [-1]]
        """
        if self._cache is None:
            self._cache = [x.to_tensor_product_of_Kirillov_Reshetikhin_tableaux()
                           for x in self.tp_krt.rigged_configurations().module_generators]
        for x in self._cache:
            yield x

    def __repr__(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: KRT = TensorProductOfKirillovReshetikhinTableaux(['D',4,1], [[2,1]])
            sage: from sage.combinat.rigged_configurations.tensor_product_kr_tableaux import HighestWeightTensorKRT
            sage: HighestWeightTensorKRT(KRT) # indirect doctest
            Highest weight elements of Tensor product of Kirillov-Reshetikhin tableaux of type ['D', 4, 1] and factor(s) ((2, 1),)
        """
        return "Highest weight elements of %s"%self.tp_krt

    @cached_method
    def cardinality(self):
        """
        Return the cardinality of ``self`` which is the number of highest
        weight elements.

        EXAMPLES::

            sage: KRT = TensorProductOfKirillovReshetikhinTableaux(['D',4,1], [[2,2]])
            sage: from sage.combinat.rigged_configurations.tensor_product_kr_tableaux import HighestWeightTensorKRT
            sage: HW = HighestWeightTensorKRT(KRT)
            sage: HW.cardinality()
            3
            sage: len(HW)
            3
            sage: len(KRT.module_generators)
            3
        """
        count = 0
        for x in self:
            count += 1
        return Integer(count)

    __len__ = cardinality

class TensorProductOfKirillovReshetikhinTableaux(FullTensorProductOfRegularCrystals):
    r"""
    A tensor product of :class:`KirillovReshetikhinTableaux`.

    Through the bijection with rigged configurations, the tableaux that are
    produced in the Kirillov-Reshetikhin model for type `D_n^{(1)}` are all of
    rectangular shapes and do not necessarily obey the usual strict increase in
    columns and weak increase in rows. The relation between the two tableaux
    models is given by a filling map.

    .. NOTE::

        The tableaux for all non-simply-laced provably holds if the bijection
        with :class:`rigged configurations <RiggedConfigurations>` holds.
        Therefore this is only proven for all factors`B^{r,1}` or all factors
        `B^{1,s}`, and in general for types `A_n^{(1)}` and `D_n^{(1)}`.

    For more information see [OSS2011]_ and
    :class:`KirillovReshetikhinTableaux`.

    REFERENCES:

    .. [OSS2011] Masato Okado, Reiho Sakamoto, Anne Schilling
       Affine crystal structure on rigged configurations of type `D_n^{(1)}`
       J. Algebraic Combinatorics, to appear, doi:10.1007/s10801-012-0383-z (arXiv:1109.3523 [math.QA])

    For more information on KR crystals, see
    :mod:`sage.combinat.crystals.kirillov_reshetikhin`.

    INPUT:

    - ``cartan_type``    -- The Cartan type

    - ``B``              -- An (ordered) list of dimensions

    The dimensions (i.e. `B`) is a list whose entries are lists of the
    form `[r, s]` which correspond to a tableau with `r` rows and `s`
    columns (or in type `A_n^{(1)}` of shape ``[r]*s``) and corresponds
    to a Kirillov-Reshetikhin crystal `B^{r,s}`.

    EXAMPLES:

    We can go between tensor products of KR crystals and rigged
    configurations::

        sage: KRT = TensorProductOfKirillovReshetikhinTableaux(['A',3,1], [[3,1],[2,2]])
        sage: tp_krt = KRT(pathlist=[[3,2,1],[3,2,3,2]]); tp_krt
        [[1], [2], [3]] (X) [[2, 2], [3, 3]]
        sage: RC = RiggedConfigurations(['A',3,1], [[3,1],[2,2]])
        sage: rc_elt = tp_krt.to_rigged_configuration(); rc_elt
        <BLANKLINE>
        -2[ ][ ]-2
        <BLANKLINE>
        0[ ][ ]0
        <BLANKLINE>
        (/)
        <BLANKLINE>
        sage: tp_krc = tp_krt.to_tensor_product_of_Kirillov_Reshetikhin_crystals(); tp_krc
        [[[1], [2], [3]], [[2, 2], [3, 3]]]
        sage: KRT(tp_krc) == tp_krt
        True
        sage: rc_elt == tp_krt.to_rigged_configuration()
        True
        sage: KR1 = KirillovReshetikhinCrystal(['A',3,1], 3,1)
        sage: KR2 = KirillovReshetikhinCrystal(['A',3,1], 2,2)
        sage: T = TensorProductOfCrystals(KR1, KR2)
        sage: t = T(KR1(3,2,1), KR2(3,2,3,2))
        sage: KRT(t) == tp_krt
        True
        sage: t == tp_krc
        True

    We can get the highest weight elements by using the attribute
    ``module_generators``::

        sage: KRT = TensorProductOfKirillovReshetikhinTableaux(['A',3,1], [[3,1], [2,1]])
        sage: list(KRT.module_generators)
        [[[1], [2], [3]] (X) [[1], [2]], [[1], [3], [4]] (X) [[1], [2]]] 
    """
    @staticmethod
    def __classcall_private__(cls, cartan_type, B):
        """
        Normalize the input arguments to ensure unique representation.

        EXAMPLES::

            sage: T1 = TensorProductOfKirillovReshetikhinTableaux(CartanType(['A',3,1]), [[2,2]])
            sage: T2 = TensorProductOfKirillovReshetikhinTableaux(['A',3,1], [(2,2)])
            sage: T3 = TensorProductOfKirillovReshetikhinTableaux(['A',3,1], ((2,2),))
            sage: T2 is T1, T3 is T1
            (True, True)
        """
        cartan_type = CartanType(cartan_type)
        if not cartan_type.is_affine():
            raise ValueError("The Cartan type must be affine")

        # Standardize B input into a tuple of tuples
        B = tuple(tuple(dim) for dim in B)
        return super(TensorProductOfKirillovReshetikhinTableaux, cls).__classcall__(cls, cartan_type, B)

    def __init__(self, cartan_type, B):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: KRT = TensorProductOfKirillovReshetikhinTableaux(['A',3,1], [[3,1],[2,2]]); KRT
            Tensor product of Kirillov-Reshetikhin tableaux of type ['A', 3, 1] and factor(s) ((3, 1), (2, 2))
            sage: TestSuite(KRT).run() # long time
            sage: KRT = TensorProductOfKirillovReshetikhinTableaux(['D',4,1], [[2,2]])
            sage: TestSuite(KRT).run() # long time
            sage: KRT = TensorProductOfKirillovReshetikhinTableaux(['D',4,1], [[3,1]])
            sage: TestSuite(KRT).run() # long time
            sage: KRT = TensorProductOfKirillovReshetikhinTableaux(['D',4,1], [[4,3]])
            sage: TestSuite(KRT).run() # long time
        """
        self.dims = B
        self.letters = CrystalOfLetters(cartan_type.classical())
        tensor_prod = []
        for rect_dims in B:
            tensor_prod.append(KirillovReshetikhinTableaux(cartan_type, rect_dims[0], rect_dims[1]))
        FullTensorProductOfRegularCrystals.__init__(self, tuple(tensor_prod), cartan_type=cartan_type)
        # This is needed to override the module_generators set in FullTensorProductOfRegularCrystals
        self.module_generators = HighestWeightTensorKRT(self)
        self.rename("Tensor product of Kirillov-Reshetikhin tableaux of type %s and factor(s) %s"%(\
          cartan_type, B))

    def __iter__(self):
        """
        Returns the iterator of ``self``.

        EXAMPLES::

            sage: KRT = TensorProductOfKirillovReshetikhinTableaux(['A', 3, 1], [[2,1], [1,1]])
            sage: g = KRT.__iter__()
            sage: g.next()
            [[2], [3]] (X) [[1]]
            sage: g.next()
            [[2], [4]] (X) [[1]]
        """
        index_set = self._cartan_type.classical().index_set()
        from sage.combinat.backtrack import TransitiveIdeal
        return TransitiveIdeal(lambda x: [x.f(i) for i in index_set],
                               self.module_generators).__iter__()

    def _test_bijection(self, **options):
        r"""
        Test function to make sure that the bijection between rigged
        configurations and Kirillov-Reshetikhin tableaux is correct.

        EXAMPLES::

            sage: KRT = TensorProductOfKirillovReshetikhinTableaux(['A', 4, 1], [[3,2],[4,1]])
            sage: KRT._test_bijection()
        """
        tester = self._tester(**options)
        rejects = []
        for x in self:
            y = x.to_rigged_configuration()
            z = y.to_tensor_product_of_Kirillov_Reshetikhin_tableaux()
            if z != x:
                rejects.append((x, z))

        tester.assertTrue(len(rejects) == 0, "Bijection is not correct: %s"%rejects)
        if len(rejects) != 0:
            return rejects

    def _element_constructor_(self, *path, **options):
        r"""
        Construct a TensorProductOfKRTableauxElement.

        Typically the user will call this with the option **pathlist** which
        will receive a list and coerce it into a path.

        EXAMPLES::

            sage: KRT = TensorProductOfKirillovReshetikhinTableaux(['A',3,1], [[3,1], [2,1]])
            sage: KRT(pathlist=[[4, 2, 1], [2, 1]]) # indirect doctest
            [[1], [2], [4]] (X) [[1], [2]]
        """
        if isinstance(path[0], KirillovReshetikhinTableauxElement):
            return self.element_class(self, path)
        if isinstance(path[0], TensorProductOfKirillovReshetikhinTableauxElement):
            return path[0]

        from sage.combinat.crystals.kirillov_reshetikhin import KirillovReshetikhinGenericCrystalElement
        if isinstance(path[0], KirillovReshetikhinGenericCrystalElement):
            return self.element_class(self, [x.to_Kirillov_Reshetikhin_tableau() for x in path])

        from sage.combinat.crystals.tensor_product import TensorProductOfRegularCrystalsElement
        if isinstance(path[0], TensorProductOfRegularCrystalsElement) and \
          isinstance(path[0][0], KirillovReshetikhinGenericCrystalElement):
            return self.element_class(self, [x.to_Kirillov_Reshetikhin_tableau() for x in path[0]])

        from sage.combinat.rigged_configurations.rigged_configuration_element import RiggedConfigurationElement
        if isinstance(path[0], RiggedConfigurationElement):
            if self.rigged_configurations() != path[0].parent():
                raise ValueError("incorrect bijection image")
            return path[0].to_tensor_product_of_Kirillov_Reshetikhin_tableaux()

        return self.element_class(self, list(path), **options)

    @cached_method
    def _module_generators_brute_force(self):
        """
        Return the module generators of ``self`` by brute force searching
        through all elements of ``self`` as a Cartesian product.

        EXAMPLES::

            sage: KRT = TensorProductOfKirillovReshetikhinTableaux(['A',3,1], [[1,3], [2,1]])
            sage: tuple(KRT.module_generators)
            ([[1, 1, 1]] (X) [[1], [2]], [[1, 1, 3]] (X) [[1], [2]])
            sage: KRT._module_generators_brute_force()
            ([[1, 1, 1]] (X) [[1], [2]], [[1, 1, 3]] (X) [[1], [2]])
        """
        index_set = self.cartan_type().classical().index_set()
        return tuple(x for x in FullTensorProductOfRegularCrystals.__iter__(self)
                                if x.is_highest_weight(index_set))

    @cached_method
    def rigged_configurations(self):
        """
        Return the corresponding set of rigged configurations.

        EXAMPLES::

            sage: KRT = TensorProductOfKirillovReshetikhinTableaux(['A',3,1], [[1,3], [2,1]])
            sage: KRT.rigged_configurations()
            Rigged configurations of type ['A', 3, 1] and factor(s) ((1, 3), (2, 1))
        """
        from sage.combinat.rigged_configurations.rigged_configurations import RiggedConfigurations
        return RiggedConfigurations(self.cartan_type(), self.dims)

    @cached_method
    def tensor_product_of_Kirillov_Reshetikhin_crystals(self):
        """
        Return the corresponding tensor product of Kirillov-Reshetikhin
        crystals.

        EXAMPLES::

            sage: KRT = TensorProductOfKirillovReshetikhinTableaux(['A',3,1], [[3,1],[2,2]])
            sage: KRT.tensor_product_of_Kirillov_Reshetikhin_crystals()
            Full tensor product of the crystals [Kirillov-Reshetikhin crystal of type ['A', 3, 1] with (r,s)=(3,1),
            Kirillov-Reshetikhin crystal of type ['A', 3, 1] with (r,s)=(2,2)]

        TESTS::

            sage: KRT = TensorProductOfKirillovReshetikhinTableaux(['D', 4, 1], [[4,1], [3,3]])
            sage: KR1 = KirillovReshetikhinCrystal(['D', 4, 1], 4, 1)
            sage: KR2 = KirillovReshetikhinCrystal(['D', 4, 1], 3, 3)
            sage: T = TensorProductOfCrystals(KR1, KR2)
            sage: T == KRT.tensor_product_of_Kirillov_Reshetikhin_crystals()
            True
            sage: T is KRT.tensor_product_of_Kirillov_Reshetikhin_crystals()
            True
        """
        return FullTensorProductOfRegularCrystals(tuple(x.Kirillov_Reshetikhin_crystal() for x in self.crystals),
                                           cartan_type=self.cartan_type())

TensorProductOfKirillovReshetikhinTableaux.Element = TensorProductOfKirillovReshetikhinTableauxElement

