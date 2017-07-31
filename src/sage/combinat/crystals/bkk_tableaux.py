from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.structure.parent import Parent
from sage.combinat.tableau import Tableaux, SemistandardTableaux
from sage.combinat.skew_tableau import SkewTableau, SkewTableaux, SemistandardSkewTableaux
from sage.structure.unique_representation import UniqueRepresentation
from sage.combinat.crystals.bkk_crystals import BKKCrystalCategory

class BKKTableaux(UniqueRepresentation, Parent):
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

    # TODO: maybe we should create our own dictionary of options? For now we
    # use that of Tableaux.
    options = SkewTableaux.options

    @staticmethod
    def __classcall_private__(cls, shape, m, n):
        from sage.combinat.skew_partition import SkewPartition
        try:
            shape = SkewPartition(shape)
        except ValueError:
            shape = SkewPartition((shape, []))
        return super(BKKTableaux, cls).__classcall__(cls, shape, m, n)

    def __init__(self, shape, m, n):
        r"""
        EXAMPLES::

            sage: from bkk_tableaux import BKKTableaux
            sage: T = BKKTableaux([2,1], 2, 2)
            sage: T
            BKK tableaux of skew shape [[2, 1], []] and entries from (-2, -1, 1, 2)

            sage: TestSuite(T).run(skip=["_test_elements", "_test_pickling"])

        """
        self._shape = shape
        self._m = m
        self._n = n
        self._index_set = tuple(range(-m + 1, n))
        if self._shape[1] != []:
            raise NotImplementedError
        tab = [[i-self._m]*self._shape[0][i] for i in range(min(len(self._shape[0]), self._m))]
        for ell in self._shape[0][self._m:]:
            if ell > self._n:
                raise ValueError("Invalid shape")
            tab.append(list(range(1,ell+1)))
        self.module_generators = (self.element_class(self, SkewTableau(tab)),)
        Parent.__init__(self, category=BKKCrystalCategory())

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

    def __iter__old(self):
        r"""

        EXAMPLES::

            sage: from bkk_tableaux import BKKTableaux
            sage: BKKTableaux([3], 3, 3).cardinality()
            38
            sage: ascii_art(BKKTableaux([3], 0, 4).list())
            [   1  2  3,   1  2  4,   1  3  4,   2  3  4 ]
            sage: ascii_art(BKKTableaux([3], 2, 0).list())
            [  -2 -2 -2,  -2 -2 -1,  -2 -1 -1,  -1 -1 -1 ]
            sage: ascii_art(BKKTableaux([2,1,1], 3, 0).list())
            [  -3 -3   -3 -2   -3 -1 ]
            [  -2      -2      -2    ]
            [  -1   ,  -1   ,  -1    ]
            sage: ascii_art(BKKTableaux([1,1,1], 0, 3).list())
            [   1    1    1    1    1    1    2    2    2    3 ]
            [   2    1    1    1    2    3    2    2    3    3 ]
            [   3,   1,   2,   3,   2,   3,   2,   3,   3,   3 ]

            sage: B = BKKTableaux(([3,1,1], [2,1]), 2, 2)

        """
        alphabet_dict = dict(enumerate(self.alphabet(), 1))
        max_entry = len(alphabet_dict)
        alphabet_dict[None] = None
        shape = self.shape()
        seen = set()
        for t in SemistandardSkewTableaux(shape, max_entry=max_entry):
            s = SkewTableau([[alphabet_dict[entry] for entry in row] for row in t])
            if self._check_verifies_bkk_conditions(s):
                if s not in seen:
                    seen.add(s)
                    yield self.element_class(self, s)
        # FIXME: this is over generating!
        for t in SemistandardSkewTableaux(shape.conjugate(), max_entry=max_entry):
            s = SkewTableau([[alphabet_dict[entry] for entry in row] for row in t]).conjugate()
            if self._check_verifies_bkk_conditions(s):
                if s not in seen:
                    seen.add(s)
                    yield self.element_class(self, s)

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

def arabic_reading_word(t):
    r"""
    An Arabic reading proceeds across rows from right to left and from top to
    bottom.

    EXAMPLES::

        sage: from bkk_tableaux import arabic_reading_word
        sage: s = SkewTableau([[None, -3, -2], [-4, -1, 1], [1, 3], [2]])
        sage: ascii_art(s)
          . -3 -2
         -4 -1  1
          1  3
          2
        sage: arabic_reading_word(s)
        word: -2,-3,1,-1,-4,3,1,2
    """
    from sage.combinat.words.word import Word
    w = []
    for row in t:
        for entry in row[::-1]:
            if entry is not None:
                w.append(entry)
    return Word(w)

def skew_partition_iterator(n):
    from sage.combinat.skew_partition import SkewPartition
    from sage.combinat.partition import Partitions
    for partition in Partitions(n):
        for mu in Partitions():
            if mu.size() == n+1:
                break
            try:
                sk = SkewPartition([partition, mu])
            except ValueError:
                continue
            yield sk


