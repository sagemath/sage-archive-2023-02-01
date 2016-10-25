r"""
Littlewood-Richardson tableaux

A semistandard tableau is Littlewood-Richardson with respect to
the sequence of partitions `(\mu^{(1)},\ldots,\mu^{(k)})` if,
when restricted to each alphabet `\{|\mu^{(1)}|+\cdots+|\mu^{(i-1)}|+1,\ldots,
|\mu^{(1)}|+\cdots+|\mu^{(i)}|-1\}`, is Yamanouchi.

AUTHORS:

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
from sage.libs.symmetrica.all import kostka_tab

class LittlewoodRichardsonTableau(SemistandardTableau):
    r"""
    A semistandard tableau is Littlewood-Richardson with respect to
    the sequence of partitions `(\mu^{(1)}, \ldots, \mu^{(k)})` if,
    when restricted to each alphabet `\{|\mu^{(1)}|+\cdots+|\mu^{(i-1)}|+1,
    \ldots, |\mu^{(1)}|+\cdots+|\mu^{(i)}|-1\}`, is Yamanouchi.

    INPUT:

    - ``t`` -- Littlewood-Richardson tableau; the input is supposed to be
      a list of lists specifying the rows of the tableau

    EXAMPLES::

        sage: from sage.combinat.lr_tableau import LittlewoodRichardsonTableau
        sage: LittlewoodRichardsonTableau([[1,1,3],[2,3],[4]], [[2,1],[2,1]])
        [[1, 1, 3], [2, 3], [4]]
    """
    @staticmethod
    def __classcall_private__(cls, t, weight):
        r"""
        Implements the shortcut ``LittlewoodRichardsonTableau(t, weight)`` to
        ``LittlewoodRichardsonTableaux(shape , weight)(t)``
        where ``shape`` is the shape of the tableau.

        TESTS::

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
        Initialize ``self``.

        TESTS::

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
        super(LittlewoodRichardsonTableau, self).__init__(parent, list(t))

    def check(self):
        r"""
        Check that ``self`` is a valid Littlewood-Richardson tableau.

        EXAMPLES::

            sage: from sage.combinat.lr_tableau import LittlewoodRichardsonTableau
            sage: t = LittlewoodRichardsonTableau([[1,1,3],[2,3],[4]], [[2,1],[2,1]])
            sage: t.check()

        TESTS::

            sage: LR = LittlewoodRichardsonTableaux([3,2,1],[[2,1],[2,1]])
            sage: LR([[1, 1, 2], [3, 3], [4]])
            Traceback (most recent call last):
            ...
            ValueError: [[1, 1, 2], [3, 3], [4]] is not an element of
             Littlewood-Richardson Tableaux of shape [3, 2, 1] and weight ([2, 1], [2, 1]).
            sage: LR([[1, 1, 2, 3], [3], [4]])
            Traceback (most recent call last):
            ...
            ValueError: [[1, 1, 2, 3], [3], [4]] is not an element of
             Littlewood-Richardson Tableaux of shape [3, 2, 1] and weight ([2, 1], [2, 1]).
            sage: LR([[1, 1, 3], [3, 3], [4]])
            Traceback (most recent call last):
            ...
            ValueError: weight of the parent does not agree with the weight of the tableau
        """
        super(LittlewoodRichardsonTableau, self).check()
        if not [i for a in self.parent()._weight for i in a] == self.weight():
            raise ValueError("weight of the parent does not agree "
                             "with the weight of the tableau")
        if not self.shape() == self.parent()._shape:
            raise ValueError("shape of the parent does not agree "
                             "with the shape of the tableau")

class LittlewoodRichardsonTableaux(SemistandardTableaux):
    r"""
    Littlewood-Richardson tableaux.

    A semistandard tableau `t` is *Littlewood-Richardson* with respect to
    the sequence of partitions `(\mu^{(1)},\ldots,\mu^{(k)})` (called the
    weight) if `t` is Yamanouchi when restricted to each alphabet
    `\{|\mu^{(1)}|+\cdots+|\mu^{(i-1)}|+1, \ldots,
    |\mu^{(1)}|+\cdots+|\mu^{(i)}|-1\}`.

    INPUT:

    - ``shape`` -- the shape of the Littlewood-Richardson tableaux
    - ``weight`` -- the weight is a sequence of partitions

    EXAMPLES::

        sage: LittlewoodRichardsonTableaux([3,2,1],[[2,1],[2,1]])
        Littlewood-Richardson Tableaux of shape [3, 2, 1] and weight ([2, 1], [2, 1])
    """
    @staticmethod
    def __classcall_private__(cls, shape, weight):
        r"""
        Straighten arguments before unique representation.

        TESTS::

            sage: LR = LittlewoodRichardsonTableaux([3,2,1],[[2,1],[2,1]])
            sage: TestSuite(LR).run()
            sage: LittlewoodRichardsonTableaux([3,2,1],[[2,1]])
            Traceback (most recent call last):
            ...
            ValueError: the sizes of shapes and sequence of weights do not match
        """
        shape = Partition(shape)
        weight = tuple(Partition(a) for a in weight)
        if shape.size() != sum(a.size() for a in weight):
            raise ValueError("the sizes of shapes and sequence of weights do not match")
        return super(LittlewoodRichardsonTableaux, cls).__classcall__(cls, shape, weight)

    def __init__(self, shape, weight):
        r"""
        Initializes the parent class of Littlewood-Richardson tableaux.

        INPUT:

        - ``shape`` -- the shape of the Littlewood-Richardson tableaux
        - ``weight`` -- the weight is a sequence of partitions

        TESTS::

            sage: LR = LittlewoodRichardsonTableaux([3,2,1],[[2,1],[2,1]])
            sage: TestSuite(LR).run()
        """
        self._shape = shape
        self._weight = weight
        self._heights = [a.length() for a in self._weight]
        super(LittlewoodRichardsonTableaux, self).__init__(category=FiniteEnumeratedSets())

    def _repr_(self):
        """
        TESTS::

            sage: LittlewoodRichardsonTableaux([3,2,1],[[2,1],[2,1]])
            Littlewood-Richardson Tableaux of shape [3, 2, 1] and weight ([2, 1], [2, 1])
        """
        return "Littlewood-Richardson Tableaux of shape %s and weight %s"%(self._shape, self._weight)

    def __iter__(self):
        r"""
        TESTS::

            sage: LR = LittlewoodRichardsonTableaux([3,2,1], [[2,1],[2,1]])
            sage: LR.list()
            [[[1, 1, 3], [2, 3], [4]], [[1, 1, 3], [2, 4], [3]]]
        """
        exp = [i for a in self._weight for i in a]
        for t in kostka_tab(self._shape, exp):
            if is_littlewood_richardson(t, self._heights):
                yield self.element_class(self, t)

    def __contains__(self, t):
        """
        Check if ``t`` is contained in ``self``.

        TESTS::

            sage: LR = LittlewoodRichardsonTableaux([3,2,1], [[2,1],[2,1]])
            sage: SST = SemistandardTableaux([3,2,1], [2,1,2,1])
            sage: [t for t in SST if t in LR]
            [[[1, 1, 3], [2, 3], [4]], [[1, 1, 3], [2, 4], [3]]]
            sage: [t for t in SST if t in LR] == LR.list()
            True

            sage: LR = LittlewoodRichardsonTableaux([3,2,1], [[2,1],[2,1]])
            sage: T = [[1,1,3], [2,3], [4]]
            sage: T in LR
            True
        """
        return (SemistandardTableaux.__contains__(self, t)
                and is_littlewood_richardson(t, self._heights))

    Element = LittlewoodRichardsonTableau

#### common or global functions related to LR tableaux

def is_littlewood_richardson(t, heights):
    """
    Return whether semistandard tableau ``t`` is Littleword-Richardson
    with respect to ``heights``.

    A tableau is Littlewood-Richardson with respect to ``heights`` given by
    `(h_1,h_2,\ldots)` if each subtableau with respect to the alphabets
    `\{1, 2, \ldots, h_1\}`,`\{h_1+1, \ldots, h_1+h_2\}`,
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
    try:
        w = t.to_word()
    except AttributeError:  # Not an instance of Tableau
        w = sum(reversed(t), [])
    for i in range(len(heights)):
        alphabet = set(range(partial[i]+1, partial[i+1]+1))
        subword = Word([j for j in w if j in alphabet])
        if not subword.is_yamanouchi():
            return False
    return True

# jake's additions here

from sage.libs.lrcalc.lrcalc import lrskew
import itertools
from sage.combinat.skew_tableau import SkewTableau

def LRtabs_single(outer, weight, containing_inner=None, shift=0):
    """
    An iterator for the Yamanouchi skew tableaux of specified
    outer shape ``outer`` and weight, but with arbitrary inner shape.

    The optional argument ``containing_inner`` specifies a minimal inner
    shape, i.e. the output tableaux must fit inside the skew shape
    ``outer``/``containing_inner``, (though they do not need to fill
    it up).

    The optional argument ``shift``, a positive integer, increases
    all the entries of the tableau by that amount.
    """


    # TODO: fix lrcalc to check this condition!
    # (this appears to be a bug in lrcalc)
    if not Partition(outer).contains(weight):
        return

    # The highest-weight standard tableau of shape w (for use below)
    w = SkewTableau([[i+1]*val for i, val in enumerate(weight)])
    w = w.standardization()

    # Step 1. From lrcalc, get an iterator for the LR tableaux of
    # skew shape = outer/weight, and arbitrary content
    related_tabs = lrskew(outer, weight)

    # Step 2. Convert them to the desired tableaux using JDT
    for t in related_tabs:
        # Shuffle t with w (returns (rect(t),std))
        # Note: std now has shape outer/(something) and content = weight
        std = w.shuffle(t.standardization())[1] # shuffle returns (rect(t),std)

        # Step 3. Check the inner shape (if necessary)
        if (containing_inner is not None) and not(std.inner_shape().contains(containing_inner)):
            continue

        # Step 4. Destandardize (shifting if necessary)
        yield _destandardize(std, weight, shift)

def LRtabs_multi(outer, weights, inner=[]):
    """
    An iterator for the LR tableaux of skew shape outer/inner,
    which are unions of Yamanouchi sub-tableaux of the specified
    weights.

    Example usage:
        sage: LRtabs_multi([7, 4, 2], [[2, 1], [3, 2], [2]], [2, 1])
        ...
        Returns: the 25 tableaux of shape [7, 4, 2] / [2, 1], such that
            the subtableau of {1, 2}'s is LR of content [2, 1]
            the subtableau of {3, 4}'s is LR of content [3, 2]
            the subtableau of    {5}'s is LR of content [2]

    This should be equivalent to calling
        ``LittlewoodRichardsonTableax(outer, [inner]+weights)``
    except that the inner shape is left blank.

    The algorithm proceeds by starting with an empty tableau of
    shape = outer, then successively attaching, along the inner
    edge, an LR tableau of the next required weight (working backwards).
    """
    shift = sum([len(w) for w in weights])
    start = SkewTableau([[None]*l for l in outer])
    return _LRtabs_multi_recursive(start, list(reversed(weights)), inner, shift)

# This function builds up the tableaux by starting at the outer edge,
# then successively attaching, along the inner edge, an LR tableau 
# of the next required weight
def _LRtabs_multi_recursive(partial_tab, weights, inner, shift):
    if weights == []:
        if partial_tab.inner_shape() == inner:
            yield partial_tab
        return

    w = weights[0]
    for t in LRtabs_single(partial_tab.inner_shape(), w, inner, shift-len(w)):
        for res in _LRtabs_multi_recursive(_tableauJoin(t, partial_tab),
                                           weights[1:],
                                           inner,
                                           shift-len(w)):
            yield res

def _tableauJoin(t1, t2):
    """
    A helper function:
    Concatenate the rows of t1 and t2, dropping any None's from t2.
    (Intended for the case where t1.outer_shape() == t2.inner_shape().)
    """
    return SkewTableau([[e1 for e1 in row1]+[e2 for e2 in row2 if e2 is not None] for (row1, row2) in itertools.izip_longest(t1, t2, fillvalue=[])])

def _destandardize(T, beta, shift=0):
    """
    A helper function:
    Convert a standard tableau back to a semistandard one
    of specified content ``beta`` and entries in the alphabet
    {``shift``+1, ..., ``shift``+``len(beta)``}.

    NOTE:
    There is no error-checking! It just changes the entries
    1, ..., beta[0] to 1's, then beta[0]+1, ..., beta[1] to 2's, and so on.
    """
    S = [[e for e in row] for row in T]
    conversion = [i+1 for i, val in enumerate(beta) for _ in range(val)]
    for (r, c) in T.cells():
        S[r][c] = conversion[T[r][c]-1]+shift
    return SkewTableau(S)
