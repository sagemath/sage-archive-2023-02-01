r"""
Tensor product of Kirillov-Reshetikhin tableaux

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
    Tensor product of Kirillov-Reshetikhin tableaux of type ['A', 3, 1] and tableau shape(s) [[1, 1, 1], [1, 1]]
    sage: KRT.cardinality()
    24
    sage: HW = HighestWeightTensorProductOfKirillovReshetikhinTableaux(['A',3,1], [[1,1], [2,1], [3,1]])
    sage: HW
    Highest weight tensor product of Kirillov-Reshetikhin tableaux of type ['A', 3, 1] and tableau shape(s) [[1], [1, 1], [1, 1, 1]]
    sage: HW.cardinality()
    5
    sage: len(HW)
    5
    sage: KRT = TensorProductOfKirillovReshetikhinTableaux(['A',3,1], [[1,1], [2,1], [3,1]])
    sage: KRT.cardinality()
    96

Type `D_n^{(1)}` examples::

    sage: KRT = TensorProductOfKirillovReshetikhinTableaux(['D', 4, 1], [[1, 1], [2, 1], [1, 1]])
    sage: KRT
    Tensor product of Kirillov-Reshetikhin tableaux of type ['D', 4, 1] and tableau shape(s) [[1], [1, 1], [1]]
    sage: HW = HighestWeightTensorProductOfKirillovReshetikhinTableaux(['D', 4, 1], [[1, 1], [2, 1], [1, 1]])
    sage: T = HW(pathlist=[[1], [-2, 2], [1]])
    sage: T
    [[1]] (X) [[2], [-2]] (X) [[1]]
    sage: T2 = HW(pathlist=[[1], [2, -2], [1]])
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

from sage.combinat.crystals.tensor_product import FullTensorProductOfRegularCrystals
from sage.combinat.crystals.letters import CrystalOfLetters
from sage.combinat.root_system.cartan_type import CartanType

from sage.combinat.rigged_configurations.tensor_product_kr_tableaux_element import TensorProductOfKirillovReshetikhinTableauxElement
from sage.combinat.rigged_configurations.kr_tableaux import KirillovReshetikhinTableaux

from sage.rings.integer import Integer

class AbstractTensorProductOfKRTableaux(FullTensorProductOfRegularCrystals):
    r"""
    Abstract class for all of tensor product of KR tableaux of a given Cartan type.

    See :class:`TensorProductOfKirillovReshetikhinTableaux`. This class should
    never be created directly.
    """

    def __init__(self, cartan_type, B, biject_class):
        r"""
        Construct a tensor product of KR tableaux.

        INPUT:

        - ``cartan_type``    -- The crystal type and n value
        - ``B``              -- An (ordered) list of dimensions
        - ``biject_class``   -- The class the bijection creates

        The dimensions (i.e. `B`) is a list whose entries are lists of the
        form `[r, s]` which correspond to a tableau with `r` rows and `s`
        columns (or of shape `[r]*s`) and corresponds to a
        Kirillov-Reshetikhin crystal `B^{r,s}`.

        TESTS::

            sage: HW = HighestWeightTensorProductOfKirillovReshetikhinTableaux(['A',3,1], [[3,1], [2,2]]); HW # indirect doctest
            Highest weight tensor product of Kirillov-Reshetikhin tableaux of type ['A', 3, 1] and tableau shape(s) [[1, 1, 1], [2, 2]]
        """
        assert cartan_type.is_affine()

        self.affine_ct = cartan_type
        self.dims = B
        self.letters = CrystalOfLetters(cartan_type)
        self._bijection_class = biject_class
        tensorProd = []
        for rectDims in B:
            tensorProd.append(KirillovReshetikhinTableaux(
              self.letters.cartan_type(),
              rectDims[0], rectDims[1]))
        FullTensorProductOfRegularCrystals.__init__(self, tuple(tensorProd))

    def _highest_weight_iter(self):
        r"""
        Iterate through all of the highest weight tensor product of Kirillov-Reshetikhin tableaux.

        EXAMPLES::

            sage: HW = HighestWeightTensorProductOfKirillovReshetikhinTableaux(['A',3,1], [[3,1], [2,1]])
            sage: list(HW) # indirect doctest
            [[[1], [2], [3]] (X) [[1], [2]], [[1], [3], [4]] (X) [[1], [2]]]
            sage: HW = HighestWeightTensorProductOfKirillovReshetikhinTableaux(['D', 4, 1], [[2,1]])
            sage: for x in HW: x # indirect doctest
            ...
            [[1], [2]]
            [[1], [-1]]
        """
        # This is a hack solution since during construction, the bijection will
        #   (attempt to) create a new KRT object which hasn't been fully created
        #   and stored in the UniqueRepresentation's cache. So this will be
        #   called again, causing the cycle to repeat. This hack just passes
        #   our self as an optional argument to hide it from the end-user and
        #   so we don't try to create a new KRT object.
        from sage.combinat.rigged_configurations.rigged_configurations import HighestWeightRiggedConfigurations
        for x in HighestWeightRiggedConfigurations(self.affine_ct, self.dims):
            yield x.to_tensor_product_of_Kirillov_Reshetikhin_tableaux(KRT_init_hack=self)

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
        from sage.combinat.crystals.kirillov_reshetikhin import KirillovReshetikhinGenericCrystalElement
        if isinstance(path[0], KirillovReshetikhinGenericCrystalElement):
            return self.element_class(self, *[x.to_Kirillov_Reshetikhin_tableau() for x in path])

        from sage.combinat.crystals.tensor_product import TensorProductOfRegularCrystalsElement
        if isinstance(path[0], TensorProductOfRegularCrystalsElement) and \
          isinstance(path[0][0], KirillovReshetikhinGenericCrystalElement):
            return self.element_class(self, *[x.to_Kirillov_Reshetikhin_tableau() for x in path[0]])

        from sage.combinat.rigged_configurations.rigged_configuration_element import RiggedConfigurationElement
        if isinstance(path[0], RiggedConfigurationElement):
            if self.rigged_configurations() != path[0].parent():
                raise ValueError("Incorrect bijection image.")
            return path[0].to_tensor_product_of_Kirillov_Reshetikhin_tableaux()

        return self.element_class(self, *path, **options)

    def _convert_to_letters(self, index, tableauList):
        """
        Convert the entries of the list to a list of letters.

        This is a helper function to convert the list of ints to letters since
        we do not convert an int to an Integer at compile time.

        TESTS::

            sage: KRT = TensorProductOfKirillovReshetikhinTableaux(['A',3,1], [[1,3]])
            sage: L = KRT._convert_to_letters(0, [3, 2, 2]); L
            [3, 2, 2]
            sage: type(L[0])
            <type 'sage.combinat.crystals.letters.Crystal_of_letters_type_A_element'>
            sage: L[0].value
            3
        """
        return([self.letters(x) for x in tableauList])

    def rigged_configurations(self):
        """
        Return the corresponding set of rigged configurations.

        EXAMPLES::

            sage: KRT = TensorProductOfKirillovReshetikhinTableaux(['A',3,1], [[1,3], [2,1]])
            sage: KRT.rigged_configurations()
            Rigged configurations of type ['A', 3, 1] and factors ((1, 3), (2, 1))
        """
        return self._bijection_class(self.affine_ct, self.dims)

    def list(self):
        r"""
        Create a list of the elements by using the iterator.

        TESTS::

            sage: HW = HighestWeightTensorProductOfKirillovReshetikhinTableaux(['A',3,1], [[3,1], [2,1]])
            sage: HW.list()
            [[[1], [2], [3]] (X) [[1], [2]], [[1], [3], [4]] (X) [[1], [2]]]
        """
        # This is needed to overwrite the list method from the FiniteCrystals
        #   category which generates the list via f_a applications.
        return [x for x in self]

class HighestWeightTensorProductOfKirillovReshetikhinTableaux(AbstractTensorProductOfKRTableaux):
    r"""
    Container class of all highest weight tensor product of KR tableaux.

    A tensor product of KR tableaux is highest weight if the action of `e_i`
    for `i \in I \setminus \{0\}` are all undefined.

    For more on tensor product of Kirillov-Reshetikhin tableaux, see
    :class:`TensorProductOfKirillovReshetikhinTableaux`.
    """

    @staticmethod
    def __classcall_private__(cls, cartan_type, B):
        """
        Normalizes the input arguments to ensure unique representation.

        EXAMPLES::

            sage: T1 = HighestWeightTensorProductOfKirillovReshetikhinTableaux(CartanType(['A',3,1]), [[2,2]])
            sage: T2 = HighestWeightTensorProductOfKirillovReshetikhinTableaux(['A',3,1], [(2,2)])
            sage: T3 = HighestWeightTensorProductOfKirillovReshetikhinTableaux(['A',3,1], ((2,2),))
            sage: T2 is T1, T3 is T1
            (True, True)
        """
        cartan_type = CartanType(cartan_type)
        # Standardize B input into a tuple of tuples
        assert B is not None
        B = tuple(tuple(dim) for dim in B)
        return super(HighestWeightTensorProductOfKirillovReshetikhinTableaux, cls).__classcall__(cls, cartan_type, B)

    def __init__(self, cartan_type, B):
        r"""
        Initialize the class.

        INPUT:

        - ``cartan_type``   -- The crystal type and n value
        - ``B``             -- An (ordered) list of dimensions.

        The dimensions (i.e. `B`) is a list whose entries are lists of the
        form `[r, c]` which correspond to a tableau with `r` rows and `c`
        columns (or of shape `[r]*c`).

        EXAMPLES::

            sage: HW = HighestWeightTensorProductOfKirillovReshetikhinTableaux(['A',3,1], [[3,1], [2,2]]); HW
            Highest weight tensor product of Kirillov-Reshetikhin tableaux of type ['A', 3, 1] and tableau shape(s) [[1, 1, 1], [2, 2]]
            sage: TestSuite(HW).run()

        TESTS:

        That `__len__()` returns the correct value::

            sage: HW = HighestWeightTensorProductOfKirillovReshetikhinTableaux(['A',3,1], [[3,1], [2,2]])
            sage: len(HW)
            2
            sage: HW.cardinality() == len(HW)
            True
        """
        from rigged_configurations import HighestWeightRiggedConfigurations
        AbstractTensorProductOfKRTableaux.__init__(self, cartan_type, B, HighestWeightRiggedConfigurations)
        self.rename("Highest weight tensor product of Kirillov-Reshetikhin tableaux of type %s and tableau shape(s) %s" % (\
          cartan_type, list([rectDims[1]] * rectDims[0] for rectDims in B)))

    __iter__ = AbstractTensorProductOfKRTableaux._highest_weight_iter

#        Old iterator doc-string
#        For type `A`, this works via a backtracing algorithm, however for all
#        other types, it currently just iterates through all of them, checking
#        to see if they are highest weight.

    # TODO Currently this will only be constructed for type A_n, so it needs
    #   to be generalized to all types.
    # TODO Make it work for multiple columns
    def _highest_weight_iterator_type_A(self):
        r"""
        Iterate over the highest weight tensor product of Kirillov-Reshetikhin
        tableaux for type `A`.

        This uses a simple backtracing algorithm to determine all of the
        highest weight crystals and will be significantly faster than iterating
        over all crystals and checking if they are highest weight.

        TESTS::

            sage: HW = HighestWeightTensorProductOfKirillovReshetikhinTableaux(['A',3,1], [[3,1], [2,1]])
            sage: for x in HW._highest_weight_iterator_type_A(): x
            ...
            [[1], [2], [3]] (X) [[1], [2]]
            [[1], [3], [4]] (X) [[1], [2]]
        """
        pathLen = len(self.crystals) # This corresponds to the number of tableaux
        curPath = [] # Stupid python references...
        for i in range(pathLen):
            curPath.append([])
        maxValue = self._cartan_type.n + 1
        # Specifies which tableau in the tensor product we are working in
        tensorIndex = pathLen - 1
        value = 1 # The value we are trying to place at the current location
        # Count of all the values to verify we have a highest weight crystal
        valueCount = [0] * maxValue

        while True: # There is a return statement called if tensorIndex >= pathLen
            if value <= maxValue:
                if (value == 1 or valueCount[value - 2] >= valueCount[value - 1] + 1) \
                  and (value == maxValue or
                  valueCount[value - 1] + 1 >= valueCount[value]):
                    # If it still makes it a highest weight crystal
                    # Prepend the value to our current location
                    curPath[tensorIndex].insert(0, value) # Assuming type A_n
                    valueCount[value - 1] += 1 # value-1 to convert to comp. indices

                    # Assuming type A_n in here as well
                    if len(curPath[tensorIndex]) == self.dims[tensorIndex][0]:
                        # If the current tableau is filled
                        if tensorIndex == 0:
                            # If it is the last tableau in the path, then add it to the
                            #   return list and then pop it off and increment the value
                            # Coerce into a list of crystal of tableaux
                            yield self(*[self.crystals[i](self._convert_to_letters(\
                              i, tab)) for i, tab in enumerate(curPath)])
                            valueCount[curPath[tensorIndex].pop(0) - 1] -= 1
                            value += 1
                        else:
                            # Otherwise we just continue on to the next tableau
                            tensorIndex -= 1
                            value = 1
                    else:
                        # Otherwise the tableau still has more values to add and since
                        #   is an SSYT, it must be column strict inequality (note we
                        #   are subtly again assuming its type A_n)
                        value += 1
                else:
                    # Otherwise the current value will not work
                    value += 1
            else:
                # Otherwise its time to backtrace since we've exhausted all
                #   possible values for this position

                # If the current tableau is empty, back up to the previous tableau
                if len(curPath[tensorIndex]) == 0:
                    tensorIndex += 1
                    if tensorIndex == pathLen:
                        # Nothing more to do since there are no more tableau
                        return

                # Get the previous value (assuming type A_n) and increment
                value = curPath[tensorIndex].pop(0)
                valueCount[value - 1] -= 1
                value += 1

    @cached_method
    def cardinality(self):
        r"""
        Return the number of highest weight tensor product of Kirillov-Reshetikhin tableaux.

        EXAMPLES::

            sage: HW = HighestWeightTensorProductOfKirillovReshetikhinTableaux(['A', 4, 1], [[2, 1]])
            sage: HW.cardinality()
            1
        """
        count = Integer(0)
        for x in self.module_generators:
            count += 1
        return(count)

    # To override the default crystal __len__ which is the cardinality
    # for the whole crystal.
    __len__ = cardinality

    @lazy_attribute
    def module_generators(self):
        r"""
        Module generators for this tensor product of KR tableaux.

        EXAMPLES::

            sage: HWKR = HighestWeightTensorProductOfKirillovReshetikhinTableaux(['A', 4, 1], [[2,2]])
            sage: for x in HWKR.module_generators: x
            ...
            [[1, 1], [2, 2]]
        """
        # This was separated out for speed; only needed for __iter__
        # lazy_attribute does not inherit
        return [x for x in self._highest_weight_iter()]

HighestWeightTensorProductOfKirillovReshetikhinTableaux.Element = TensorProductOfKirillovReshetikhinTableauxElement

class TensorProductOfKirillovReshetikhinTableaux(AbstractTensorProductOfKRTableaux):
    r"""
    A tensor product of :class:`KirillovReshetikhinTableaux`.

    Through the bijection with rigged configurations, the tableaux that are
    produced in the Kirillov-Reshetikhin model for type `D_n^{(1)}` are all of
    rectangular shapes and do not necessarily obey the usual strict increase in
    columns and weak increase in rows. The relation between the two tableaux
    models is given by a filling map.

    For more information see [OSS2011]_ and
    :class:`KirillovReshetikhinTableaux`.

    REFERENCES:

    .. [OSS2011] Masato Okado, Reiho Sakamoto, Anne Schilling
       Affine crystal structure on rigged configurations of type `D_n^{(1)}`
       J. Algebraic Combinatorics, to appear, doi:10.1007/s10801-012-0383-z (arXiv:1109.3523 [math.QA])

    For more information on KR crystals, see
    :mod:`sage.combinat.crystals.kirillov_reshetikhin`.

    INPUT:

    - ``cartan_type`` -- An affine Cartan type
    - ``B`` -- An (ordered) list of dimensions.

    The dimensions (i.e. ``B``) is a list whose entries are lists of the
    form ``[r, s]`` which correspond to Kirillov-Reshetikhin tableaux with
    ``r`` rows and ``s`` columns.

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
        # Standardize B input into a tuple of tuples
        assert B is not None
        B = tuple(tuple(dim) for dim in B)
        return super(TensorProductOfKirillovReshetikhinTableaux, cls).__classcall__(cls, cartan_type, B)

    def __init__(self, cartan_type, B):
        r"""
        Initialize the TensorProductOfKirillovReshetikhinTableaux class.

        INPUT:

        - ``cartan_type``   -- The crystal type and n value
        - ``B``             -- An (ordered) list of dimensions.

        The dimensions (i.e. `B`) is a list whose entries are lists of the
        form `[r, c]` which correspond to a tableau with `r` rows and `c`
        columns (or of shape `[r]*c`).

        EXAMPLES::

            sage: KRT = TensorProductOfKirillovReshetikhinTableaux(['A',3,1], [[3,1],[2,2]]); KRT
            Tensor product of Kirillov-Reshetikhin tableaux of type ['A', 3, 1] and tableau shape(s) [[1, 1, 1], [2, 2]]
            sage: TestSuite(KRT).run() # long time
            sage: KRT = TensorProductOfKirillovReshetikhinTableaux(['D',4,1], [[2,2]])
            sage: TestSuite(KRT).run() # long time
            sage: KRT = TensorProductOfKirillovReshetikhinTableaux(['D',4,1], [[3,1]])
            sage: TestSuite(KRT).run() # long time
            sage: KRT = TensorProductOfKirillovReshetikhinTableaux(['D',4,1], [[4,3]])
            sage: TestSuite(KRT).run() # long time
        """
        from rigged_configurations import RiggedConfigurations
        AbstractTensorProductOfKRTableaux.__init__(self, cartan_type, B, RiggedConfigurations)
        self.rename("Tensor product of Kirillov-Reshetikhin tableaux of type %s and tableau shape(s) %s" % (\
          cartan_type, list([rectDims[1]] * rectDims[0] for rectDims in B)))
        self.module_generators = HighestWeightTensorProductOfKirillovReshetikhinTableaux(cartan_type, B)

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
                                           cartan_type=self.affine_ct)

TensorProductOfKirillovReshetikhinTableaux.Element = TensorProductOfKirillovReshetikhinTableauxElement

