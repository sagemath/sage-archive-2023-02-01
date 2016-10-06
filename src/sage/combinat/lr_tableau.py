r"""
Littlewood-Richardson tableaux

A semistandard tableau is Littlewood-Richardson with respect to
the sequence of partitions `(\mu^{(1)},\ldots,\mu^{(k)})` if,
when restricted to each alphabet `\{|\mu^{(1)}|+\cdots+|\mu^{(i-1)}|+1,\ldots,
|\mu^{(1)}|+\cdots+|\mu^{(i)}|-1\}`, is Yamanouchi.
s
Authors:

- Maria Gillespie and Anne Schilling (2016): initial version
"""
#*****************************************************************************
#       Copyright (C) 2016 Maria Gillespie
#                          Anne Schilling <anne at math.ucdavis.edu>
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

from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.structure.parent import Parent
from sage.structure.list_clone import ClonableList
from sage.combinat.tableau import SemistandardTableau, SemistandardTableaux
from sage.combinat.partition import Partition

class LittlewoodRichardsonTableau(SemistandardTableau):

    @staticmethod
    def __classcall_private__(cls, t, weight):
        r"""
        Implements the shortcut ``LittlewoodRichardsonTableau(t, weight)`` to
        ``LittlewoodRichardsonTableaux(shape , weight)(t)``
        where ``shape`` is the shape of the tableau.

        TESTS::

            sage: from sage.combinat.lr_tableau import LittlewoodRichardsonTableaux
            sage: LR = LittlewoodRichardsonTableaux([3,2,1],[[2,1],[2,1]])
            sage: t = LR([[1, 1, 3], [2, 3], [4]])
            sage: t.check()
            sage: type(t)
            <class 'sage.combinat.lr_tableau.LittlewoodRichardsonTableaux_with_category.element_class'>
            sage: TestSuite(t).run()
            sage: from sage.combinat.lr_tableau import LittlewoodRichardsonTableau
            sage: LittlewoodRichardsonTableau([[1,1,3],[2,3],[4]], [[2,1],[2,1]])
            [[1, 1, 3], [2, 3], [4]]
        """
        if isinstance(t, cls):
            return t
        tab = SemistandardTableau(list(t))
        shape = tab.shape()
        return LittlewoodRichardsonTableaux(shape, weight)(t)

    def __init__(self, parent, t):
        r"""
        Initialization of Littlewood-Richardson tableau ``t``.

        INPUT:

        - ``t`` -- Littlewood-Richardson tableau; the input is supposed to be a list
          of lists specifying the rows of the tableau.

        TESTS::

            sage: from sage.combinat.lr_tableau import LittlewoodRichardsonTableaux
            sage: LR = LittlewoodRichardsonTableaux([3,2,1],[[2,1],[2,1]])
            sage: t = LR([[1, 1, 3], [2, 3], [4]])
            sage: from sage.combinat.lr_tableau import LittlewoodRichardsonTableau
            sage: s = LittlewoodRichardsonTableau([[1,1,3],[2,3],[4]], [[2,1],[2,1]])
            sage: s == t
            True
            sage: type(t)
            <class 'sage.combinat.lr_tableau.LittlewoodRichardsonTableaux_with_category.element_class'>
            sage: t.parent()
            Littlewood-Richardson Tableaux of shape [3, 2, 1] and weight ([2, 1], [2, 1])
            sage: TestSuite(t).run()
        """
        self._shape = parent._shape
        self._weight = parent._weight
        super(LittlewoodRichardsonTableau, self).__init__(parent, t)

    def check(self):
        r"""
        Check that ``self`` is a valid Littelwood-Richardson tableau.

        EXAMPLES::

            sage: from sage.combinat.lr_tableau import LittlewoodRichardsonTableau
            sage: t = LittlewoodRichardsonTableau([[1,1,3],[2,3],[4]], [[2,1],[2,1]])
            sage: t.check()

        TESTS::

            sage: from sage.combinat.lr_tableau import LittlewoodRichardsonTableaux
            sage: LR = LittlewoodRichardsonTableaux([3,2,1],[[2,1],[2,1]])
            sage: t = LR([[1, 1, 2], [3, 3], [4]])
            Traceback (most recent call last):
            ...
            ValueError: This is not a proper Littlewood-Richardson tableau of the correct weight
            sage: t = LR([[1, 1, 3, 3], [2, 3], [4]])
            Traceback (most recent call last):
            ...
            ValueError: The weight of the parent does not agree with the weight of the tableau!
            sage: t = LR([[1, 1, 3], [3, 3], [4]])
            Traceback (most recent call last):
            ...
            ValueError: The weight of the parent does not agree with the weight of the tableau!
        """
        t = SemistandardTableau(list(self))
        if not [i for a in self.parent()._weight for i in a] == t.weight():
            raise ValueError("The weight of the parent does not agree "
                             "with the weight of the tableau!")
        if not t.shape() == self.parent()._shape:
            raise ValueError("The shape of the parent does not agree "
                             "with the shape of the tableau!")
        heights = [a.length() for a in self._weight]
        if not is_littlewood_richardson(self,heights):
            raise ValueError("This is not a proper Littlewood-Richardson tableau of the correct weight")

class LittlewoodRichardsonTableaux(SemistandardTableaux):

    @staticmethod
    def __classcall_private__(cls, shape, weight):
        r"""
        Straighten arguments before unique representation.

        TESTS::

            sage: from sage.combinat.lr_tableau import LittlewoodRichardsonTableaux
            sage: LR = LittlewoodRichardsonTableaux([3,2,1],[[2,1],[2,1]])
            sage: TestSuite(LR).run()
        """
        shape = Partition(shape)
        weight = tuple(Partition(a) for a in weight)
        if shape.size() != sum(a.size() for a in weight):
            raise ValueError("The sizes of shapes and sequence of weights do not match")
        return super(LittlewoodRichardsonTableaux, cls).__classcall__(cls, shape, weight)

    def __init__(self, shape, weight):
        r"""
        Initializes the parent class of Littlewood-Richardson tableaux.

        INPUT:

        - ``shape`` -- the shape of the Littlewood-Richardson tableaux
        - ``weight`` -- the weight is a sequence of partitions

        TESTS::

            sage: from sage.combinat.lr_tableau import LittlewoodRichardsonTableaux
            sage: LR = LittlewoodRichardsonTableaux([3,2,1],[[2,1],[2,1]])
            sage: TestSuite(LR).run()
        """
        self._shape = shape
        self._weight = weight
        super(LittlewoodRichardsonTableaux, self).__init__(category = FiniteEnumeratedSets())

    def _repr_(self):
        """
        TESTS::

            sage: from sage.combinat.lr_tableau import LittlewoodRichardsonTableaux
            sage: LittlewoodRichardsonTableaux([3,2,1],[[2,1],[2,1]])
            Littlewood-Richardson Tableaux of shape [3, 2, 1] and weight ([2, 1], [2, 1])
        """
        return "Littlewood-Richardson Tableaux of shape %s and weight %s"%(self._shape, self._weight)

    def __iter__(self):
        r"""
        TESTS::

            sage: from sage.combinat.lr_tableau import LittlewoodRichardsonTableaux
            sage: LR = LittlewoodRichardsonTableaux([3,2,1],[[2,1],[2,1]])
            sage: LR.list()
            [[[1, 1, 3], [2, 3], [4]], [[1, 1, 3], [2, 4], [3]]]
        """
        s = [[i for i in a] for a in self._weight]
        mu = sum((a for a in s),[])
        T = SemistandardTableaux(shape=self._shape,eval=mu)
        heights = [a.length() for a in self._weight]
        for t in T:
            if is_littlewood_richardson(t,heights):
                yield t

    Element = LittlewoodRichardsonTableau

#### common or global functions related to LR tableaux

def is_littlewood_richardson(t, heights):
    """
    Return whether semistandard tableau `t` is Littleword-Richardson with respect to `heights`.

    A tableau is Littlewood-Richardson with respect to `heights = (h_1,h_2,\ldots)`
    if each subtableau with respect to the alphabets `\{1,2,\ldots,h_1\}`,`\{h_1+1,\ldots, h_1+h_2\}`,
    etc. is Yamanouchi.

    EXAMPLES::

        sage: from sage.combinat.lr_tableau import is_littlewood_richardson
        sage: t = Tableau([[1,1,2,3,4],[2,3,3],[3]])
        sage: is_littlewood_richardson(t,[2,2])
        False
        sage: t = Tableau([[1,1,3],[2,3],[4,4]])
        sage: is_littlewood_richardson(t,[2,2])
        True
    """
    from sage.combinat.words.word import Word
    partial = [sum(heights[i] for i in range(j)) for j in range(len(heights)+1)]
    w = t.to_word()
    for i in range(len(heights)):
        alphabet = range(partial[i]+1,partial[i+1]+1)
        subword = Word([j for j in w if j in alphabet])
        if not subword.is_yamanouchi():
            return False
    return True
