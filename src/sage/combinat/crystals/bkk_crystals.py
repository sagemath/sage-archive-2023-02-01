"""
Benkart-Kang-Kashiwara crystals for the general-linear Lie superalgebra
"""

#*****************************************************************************
#       Copyright (C) 2017 Franco Saliola <saliola@gmail.com>
#                     2017 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

#######################################################################
#                   tensor product of BKK crystals                    #
#######################################################################

from sage.categories.cartesian_product import cartesian_product
from sage.categories.enumerated_sets import EnumeratedSets
from sage.combinat.partition import _Partitions
from sage.misc.cachefunc import cached_method
from sage.structure.element_wrapper import ElementWrapper
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.lazy_attribute import lazy_attribute

from sage.structure.parent import Parent
from sage.combinat.tableau import Tableaux, SemistandardTableaux
from sage.combinat.skew_tableau import SkewTableau, SkewTableaux, SemistandardSkewTableaux
from sage.combinat.root_system.cartan_type import CartanType
from sage.structure.unique_representation import UniqueRepresentation

from sage.combinat.crystals.letters import CrystalOfBKKLetters
from sage.combinat.crystals.letters import CrystalOfBKKLetters as BKKOneBoxCrystal # Backwards compatibility placeholder
from sage.combinat.crystals.tensor_product import TensorProductOfCrystals, CrystalOfWords
from sage.combinat.crystals.tensor_product_element import (TensorProductOfCrystalsElement,
        TensorProductOfSuperCrystalsElement, CrystalOfBKKTableauxElement)
from sage.categories.regular_supercrystals import RegularSuperCrystals

class TensorProductOfSuperCrystalsElement_old(TensorProductOfCrystalsElement):

    def f(self, i):
        # FIXME: generalize this to deal with more than two tensor factors
        (b1, b2) = self

        # Proposition 2.8.ii.a
        if i < 0:
            if b1.phi(i) > b2.epsilon(i):
                x = b1.f(i)
                if x is not None:
                    return self.parent()(x, b2)
            else:
                y = b2.f(i)
                if y is not None:
                    return self.parent()(b1, y)

        # Proposition 2.8.ii.b
        elif 0 < i:
            if b2.phi(i) > b1.epsilon(i):
                y = b2.f(i)
                if y is not None:
                    return self.parent()(b1, y)
            else:
                x = b1.f(i)
                if x is not None:
                    return self.parent()(x, b2)

        # Proposition 2.8.ii.c
        elif i == 0:
            if b1.phi(i) + b1.epsilon(i):
                x = b1.f(i)
                if x is not None:
                    return self.parent()(x, b2)
            else:
                y = b2.f(i)
                if y is not None:
                    return self.parent()(b1, y)

    def e(self, i):
        (b1, b2) = self

        # Proposition 2.8.ii.a
        if i < 0:
            if b1.phi(i) >= b2.epsilon(i):
                x = b1.e(i)
                if x is not None:
                    return self.parent()(x, b2)
            else:
                y = b2.e(i)
                if y is not None:
                    return self.parent()(b1, y)

        # Proposition 2.8.ii.b
        elif 0 < i:
            if b2.phi(i) >= b1.epsilon(i):
                y = b2.e(i)
                if y is not None:
                    return self.parent()(b1, y)
            else:
                x = b1.e(i)
                if x is not None:
                    return self.parent()(x, b2)

        # Proposition 2.8.ii.c
        elif i == 0:
            if b1.phi(i) + b1.epsilon(i):
                x = b1.e(i)
                if x is not None:
                    return self.parent()(x, b2)
            else:
                y = b2.e(i)
                if y is not None:
                    return self.parent()(b1, y)


class TensorProductOfSuperCrystals(TensorProductOfCrystals):
    r"""
    Tensor product of super crystals.

    EXAMPLES::

        sage: L = crystals.Letters(['A', [1,1]])
        sage: T = tensor([L,L,L])
        sage: T.cardinality()
        64

    TESTS:

    In the BKK paper, they point out that there are elements that are highest
    weight vectors, but not "genuine". Here we verify their examples.

    First some import statements::

        sage: from bkk_crystals import BKKOneBoxCrystal
        sage: from bkk_tableaux import japanese_reading_word

    Construct the 6-fold tensor product of the BKK crystal for shape [1]::

        sage: c = BKKOneBoxCrystal(2, 2)
        sage: # TODO: extend tensor product defintion so that tensor([c]*6) works
        sage: c6 = reduce(lambda x, y: tensor((x,y)), [c]*6)

    Here is a function that finds the element of ``c6`` whose reading word is
    the Japanese reading word of the given tableau::

        sage: def find_element_by_japanese_reading_word(t):
        ....:     w = japanese_reading_word(t)
        ....:     s = str(reduce(lambda x, y : [x, y], w))
        ....:     for x in c6:
        ....:         if str(x) == s:
        ....:             return x

        sage: t = Tableau([ [-2,-2,-2], [-1,-1], [1] ])
        sage: ascii_art(t)
         -2 -2 -2
         -1 -1
          1
        sage: japanese_reading_word(t)
        word: -2,-2,-1,-2,-1,1
        sage: find_element_by_japanese_reading_word(t)
        [[[[[-2, -2], -1], -2], -1], 1]

    The following three elements are all highest weight vectors,
    and they all belong to the same connected component of ``c6``::

        sage: H0 = Tableau([ [-2,-2,-2], [-1,-1], [1] ])
        sage: x0 = find_element_by_japanese_reading_word(H0)
        sage: x0.is_highest_weight()
        True
        sage: x0.weight()
        (3, 2, 1, 0)

        sage: H1 = Tableau([ [-2,-2,-2], [-1,2], [1] ])
        sage: x1 = find_element_by_japanese_reading_word(H1)
        sage: x1.is_highest_weight()
        True
        sage: x1.weight()
        (3, 1, 1, 1)

        sage: H2 = Tableau([ [-2,-2,2], [-1,-1], [1] ])
        sage: x2 = find_element_by_japanese_reading_word(H2)
        sage: x2.is_highest_weight()
        True
        sage: x2.weight()
        (2, 2, 1, 1)

        sage: for cc in c6.connected_components():
        ....:    if x0 in cc:
        ....:        break
        sage: (x0 in cc) and (x1 in cc) and (x2 in cc)
        True

    Similarly, the following three elements are all lowest weight vectors,
    and they all belong to the same connected component of ``c6``::

        sage: L0 = Tableau([ [-1,1,2], [1,2], [2] ])
        sage: y0 = find_element_by_japanese_reading_word(L0)
        sage: y0.is_lowest_weight()
        True
        sage: y0.weight()
        (0, 1, 2, 3)

        sage: L1 = Tableau([ [-2,1,2], [-1,2], [2] ])
        sage: y1 = find_element_by_japanese_reading_word(L1)
        sage: y1.is_lowest_weight()
        True
        sage: y1.weight()
        (1, 1, 1, 3)

        sage: L2 = Tableau([ [-2,1,2], [-1,2], [1] ])
        sage: y2 = find_element_by_japanese_reading_word(L2)
        sage: y2.is_lowest_weight()
        True
        sage: y2.weight()
        (1, 1, 2, 2)

        sage: for cc in c6.connected_components():
        ....:    if y0 in cc:
        ....:        break
        sage: (y0 in cc) and (y1 in cc) and (y2 in cc)
        True

    """
    def __init__(self, crystals, old=False):
        assert isinstance(crystals, tuple)

        if old:
            self.Element = TensorProductOfSuperCrystalsElement_old
            self._old = old

        self._cartan_type = crystals[0]._cartan_type
        Parent.__init__(self, category=RegularSuperCrystals())
        self.crystals = crystals

        self.module_generators = self

    def __iter__(self):
        for x in cartesian_product([c.list() for c in self.crystals]):
            yield self(*x)

    def _element_constructor_(self, *s):
        # Hack to build the elements (recursively)
        if len(s) == 1 and isinstance(s[0], list):
            s = s[0]
            return self.element_class(self, [C(s[i]) for i,C in enumerate(self.crystals)])
        return TensorProductOfCrystals._element_constructor_(self, *s)

    class Element(TensorProductOfSuperCrystalsElement):
        pass

#####################################################################
## Tableaux

class BKKTableaux_old(UniqueRepresentation, Parent):
    r"""
    Semistandard tableaux defined by Benkart-Kang-Kashiwara.

    These are fillings of a skew Young diagram with entries
    from the alphabet:

        -m, ..., -2, -1, 1, 2, ..., n

    subject to the following two constraints:

    - entries in each row are increasing, allowing repetition of -m, ..., -1,
      but not of 1, 2, ..., n;

    - entries in each column are increasing, allowing repetition of 1, ..., n,
      but not of -m, ..., -1.
    """

    # FIXME: Hack around the problem. Should be a proper element.
    class Element(SkewTableau):
        def e(self, i):
            read_word = japanese_reading_order(self)
            P = self.parent()
            C = BKKOneBoxCrystal(P._m, P._n)
            T = reduce(lambda x, y: tensor((x,y)), [C]*len(read_word))
            x = reduce(lambda x, y: [x, y], [C(self[r][c]) for r,c in read_word])
            x = T(x).e(i)
            if x is None:
                return None
            ret = [[None]*len(row) for row in self]
            for r,c in reversed(read_word[1:]):
                ret[r][c] = x[1].value
                x = x[0]
            r,c = read_word[0]
            ret[r][c] = x.value
            return self.__class__(P, ret)

        def f(self, i):
            read_word = japanese_reading_order(self)
            P = self.parent()
            C = BKKOneBoxCrystal(P._m, P._n)
            T = reduce(lambda x, y: tensor((x,y)), [C]*len(read_word))
            x = reduce(lambda x, y: [x, y], [C(self[r][c]) for r,c in read_word])
            x = T(x).f(i)
            if x is None:
                return None
            ret = [[None]*len(row) for row in self]
            for r,c in reversed(read_word[1:]):
                ret[r][c] = x[1].value
                x = x[0]
            r,c = read_word[0]
            ret[r][c] = x.value
            return self.__class__(P, ret)

        def weight(self):
            C = BKKOneBoxCrystal(self.parent()._m, self.parent()._n)
            elements = list(C)
            from sage.modules.free_module_element import vector
            v = vector([0]*len(elements))
            for row in self:
                for val in row:
                    i = elements.index(C(val))
                    v[i] += 1
            return v

    @staticmethod
    def __classcall_private__(cls, m, n, shape):
        return super(BKKTableaux_old, cls).__classcall__(cls, m, n, SkewPartition([shape, []]))

    def __init__(self, m, n, shape):
        r"""
        EXAMPLES::

            sage: from bkk_tableaux import BKKTableaux
            sage: T = BKKTableaux(2, 2, [2,1])
            sage: T
            BKK tableaux of skew shape [[2, 1], []] and entries from (-2, -1, 1, 2)

            sage: TestSuite(T).run(skip=["_test_elements", "_test_pickling"])

        """
        self._shape = shape
        self._m = m
        self._n = n
        if self._shape[1] != []:
            raise NotImplementedError
        tab = [[i-self._m]*self._shape[0][i] for i in range(min(len(self._shape[0]), self._m))]
        for ell in self._shape[0][self._m:]:
            if ell > self._n:
                raise ValueError("Invalid shape")
            tab.append(list(range(1,ell+1)))
        self.module_generators = (self.element_class(self, SkewTableau(tab)),)
        Parent.__init__(self, category=RegularSuperCrystals())

    def alphabet(self):
        alphabet = range(-self._m, 0) + range(1, self._n + 1)
        return tuple(alphabet)

    def _repr_(self):
        return "BKK tableaux of skew shape {} and entries from {}".format(self.shape(), self.alphabet())

    def __contains__(self, t):
        # check that t is a tableaux of the correct external shape
        if not isinstance(t, SkewTableau) and t.shape() != self._external_shape:
            return False
        return self._check_verifies_bkk_conditions(t)

    def _check_verifies_bkk_conditions(self, t):
        r"""
        EXAMPLES::

            sage: from bkk_tableaux import BKKTableaux
            sage: T = SkewTableau([[1],[1],[1]])
            sage: B = BKKTableaux([1,1,1], 2, 2)
            sage: B._check_verifies_bkk_conditions(T)
            True
        """

        alphabet_or_None = self.alphabet() + (None,)

        # entries in each row belong to -m, ..., -1, 1, 2, ..., n
        for row in t:
            for entry in row:
                if entry not in alphabet_or_None:
                    return False

        # entries in each row are increasing,
        # allowing repetition of -m, ..., -1,
        # but not of 1, 2, ..., n
        for row in t:
            for i in range(len(row)-1):
                if row[i] is not None:
                    if row[i] > row[i+1]:
                        return False
                    elif row[i] == row[i+1] and row[i] > 0:
                        return False

        # entries in each column are increasing,
        # allowing repetition of 1, 2, ..., n
        # but not of -m, ..., -1
        tc = t.conjugate()
        for row in tc:
            for i in range(len(row)-1):
                if row[i] is not None:
                    if row[i] > row[i+1]:
                        return False
                    elif row[i] == row[i+1] and row[i] < 0:
                        return False

        return True

    def shape(self):
        r"""
        EXAMPLES::

            sage: from bkk_tableaux import BKKTableaux
            sage: B = BKKTableaux([1,1,1], 2, 2)
            sage: B.shape()
            [1, 1, 1] / []
            sage: B = BKKTableaux(([3,1,1], [2,1]), 2, 2)
            sage: B.shape()
            [3, 1, 1] / [2, 1]
        """
        return self._shape

class CrystalOfBKKTableaux(CrystalOfWords):
    """
    Crystal of tableaux for type `A(m|n)`.

    EXAMPLES::

        sage: from sage.combinat.crystals.bkk_crystals import CrystalOfBKKTableaux
        sage: T = CrystalOfBKKTableaux(['A', [1,1]], [2,1])
    """
    @staticmethod
    def __classcall_private__(cls, ct, shape):
        ct = CartanType(ct)
        shape = _Partitions(shape)
        if len(shape) > ct.m + 1 and shape[ct.m] > ct.n + 1:
            raise ValueError("invalid hook shape")
        return super(CrystalOfBKKTableaux, cls).__classcall__(cls, ct, shape)

    def __init__(self, ct, shape):
        r"""
        EXAMPLES::

            sage: from sage.combinat.crystals.bkk_crystals import CrystalOfBKKTableaux
            sage: T = CrystalOfBKKTableaux(['A', [1,1]], [2,1])
            sage: T
            Crystal of BKK tableaux of skew shape [2, 1] of gl(2|2)

            sage: TestSuite(T).run()
        """
        self._shape = shape
        self._cartan_type = ct
        m = ct.m + 1
        n = ct.n + 1
        C = CrystalOfBKKLetters(ct)
        tr = shape.conjugate()
        mg = []
        for i,col_len in enumerate(tr):
            for j in range(col_len - m):
                mg.append(C(i+1))
            for j in range(max(0, m - col_len), m):
                mg.append(C(-j-1))
        mg = list(reversed(mg))
        Parent.__init__(self, category=RegularSuperCrystals())
        self.module_generators = (self.element_class(self, mg),)

    def _repr_(self):
        m = self._cartan_type.m + 1
        n = self._cartan_type.n + 1
        return "Crystal of BKK tableaux of skew shape {} of gl({}|{})".format(self.shape(), m, n)

    def shape(self):
        r"""
        EXAMPLES::

            sage: from bkk_tableaux import BKKTableaux
            sage: B = BKKTableaux([1,1,1], 2, 2)
            sage: B.shape()
            [1, 1, 1] / []
            sage: B = BKKTableaux(([3,1,1], [2,1]), 2, 2)
            sage: B.shape()
            [3, 1, 1] / [2, 1]
        """
        return self._shape

    class Element(CrystalOfBKKTableauxElement):
        pass

#######################################################################
#                          utility functions                          #
#######################################################################

def japanese_reading_word(t):
    r"""
    A Japanese reading proceeds down columns from top to bottom and from right
    to left.

    EXAMPLES::

        sage: from bkk_tableaux import japanese_reading_word
        sage: s = SkewTableau([[None, -3, -2], [-4, -1, 1], [1, 3], [2]])
        sage: ascii_art(s)
          . -3 -2
         -4 -1  1
          1  3
          2
        sage: japanese_reading_word(s)
        word: -2,1,-3,-1,3,-4,1,2
    """
    from sage.combinat.words.word import Word
    columns = list(t.conjugate())
    w = []
    for column in columns[::-1]:
        for entry in column:
            if entry is not None:
                w.append(entry)
    return Word(w)

def japanese_reading_order(t):
    r"""
    Return the coordinates of the cells in a Japanese reading of the tableau.

    A Japanese reading proceeds down columns from top to bottom and from right
    to left. 

    EXAMPLES::

        sage: from bkk_tableaux import japanese_reading_order
        sage: s = SkewTableau([[None, -3, -2], [-4, -1, 1], [1, 3], [2]])
        sage: ascii_art(s)
          . -3 -2
         -4 -1  1
          1  3
          2
        sage: japanese_reading_order(s)
        [(0, 2), (1, 2), (0, 1), (1, 1), (2, 1), (1, 0), (2, 0), (3, 0)]
    """
    columns = list(t.conjugate())
    order = []
    for j in range(len(columns))[::-1]:
        for i in range(len(columns[j])):
            if columns[j][i] is not None:
                order.append((i,j))
    return order

