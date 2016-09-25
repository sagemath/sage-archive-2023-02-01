"""
TODO:

* implement domino insertion
* implement Young-Fibonacci
* guess shape from labels
* move SkewTableau_from_partitions to where it belongs

"""

def SkewTableau_from_partitions(l):
    """
    Return the tableau corresponding to the list of partitions.

    sage: from sage.combinat.growth import SkewTableau_from_partitions
    sage: l = [[2], [2, 1], [3, 1], [4, 3, 2, 1]]
    sage: SkewTableau_from_partitions(l).pp()
    .  .  2  3
    1  3  3
    3  3
    3

    """
    shape = l[-1]
    T = [[None for _ in range(r)] for r in shape]
    for i in range(1,len(l)):
        la = l[i]
        mu = l[i-1]
        mu += [0]*(len(la) - len(mu))

        for r in range(len(la)):
            for c in range(mu[r], la[r]):
                T[r][c] = i
    return SkewTableau(T)

from sage.structure.sage_object import SageObject
from sage.combinat.words.word import Word
from sage.combinat.partition import Partition, Partitions
from sage.combinat.skew_partition import SkewPartition
from sage.combinat.skew_tableau import SkewTableau
from copy import copy

def Robinson_Schensted_Knuth_forward(shape3, shape2, shape1, content):
    """see Krattenthaler, "Growth diagrams, and increasing and
    decreasing chains in fillings of Ferrers shapes", page 15,
    (F^1 0)--(F^1 2)

    INPUT:

    - three partitions from a cell in a growth diagram, labelled as

      shape2 shape1
      shape3

    and the content of the cell.

    OUTPUT:

    - the fourth partition shape4 according to
      Robinson-Schensted-Knuth.

    """
    carry = content
    shape4 = []
    while True:
        if shape1 == []:
            row1 = 0
        else:
            row1 = shape1[0]
        if shape3 == []:
            row3 = 0
        else:
            row3 = shape3[0]
        newPart = max(row1, row3) + carry
        if newPart == 0:
            return Partition(shape4[::-1])
        else:
            shape4 = [newPart] + shape4
            if shape2 == []:
                carry = min(row1, row3)
            else:
                carry = min(row1, row3) - shape2[0]
            shape1 = shape1[1:]
            shape2 = shape2[1:]
            shape3 = shape3[1:]

def Robinson_Schensted_Knuth_backward(shape3, shape4, shape1):
    """
    see Krattenthaler, "Growth diagrams, and increasing and
    decreasing chains in fillings of Ferrers shapes", page 15,
    (B^1 0)--(B^1 2)

    INPUT:
    the three partitions from a cell in a growth diagram, labelled as

           shape1
    shape3 shape4

    OUTPUT:
    a pair (shape2, content) consisting of the shape of the fourth
    partition acording to Robinson-Schensted-Knuth and the content of
    the cell.

    """
    carry = 0
    i = len(shape4)
    shape2 = []
    while i > 0:
        if len(shape1) < i:
            row1 = 0
        else:
            row1 = shape1[i-1]
        if len(shape3) < i:
            row3 = 0
        else:
            row3 = shape3[i-1]
        shape2 = [min(row1, row3) - carry] + shape2
        carry = shape4[i-1] - max(row1, row3)
        i = i-1
    return (Partition(shape2), carry)

def Burge_forward(shape3, shape2, shape1, content):
    """
    see Krattenthaler, "Growth diagrams, and increasing and
    decreasing chains in fillings of Ferrers shapes", page 15,
    (F^4 0)--(F^4 2)

    INPUT:
    three partitions from a cell in a growth diagram, labelled as

    shape2 shape1
    shape3

    and the content of the cell.

    OUTPUT:
    the fourth partition shape4 according to Robinson-Schensted-Knuth.

    """
    carry = content
    shape4 = []
    while True:
        if shape1 == []:
            row1 = 0
        else:
            row1 = shape1[0]
        if shape2 == []:
            row2 = 0
        else:
            row2 = shape2[0]
        if shape3 == []:
            row3 = 0
        else:
            row3 = shape3[0]
        newPart = max(row1, row3) + min(ZZ(row1 == row2 == row3), carry)
        if newPart == 0:
            return Partition(shape4[::-1])
        else:
            shape4 = [newPart] + shape4
            carry = carry - min(ZZ(row1 == row2 == row3), carry) + min(row1, row3) - row2
            shape1 = shape1[1:]
            shape2 = shape2[1:]
            shape3 = shape3[1:]

def Burge_backward(shape3, shape4, shape1):
    """
    see Krattenthaler, "Growth diagrams, and increasing and
    decreasing chains in fillings of Ferrers shapes", page 15, (B^4
    0)--(B^4 2).  There is a typo in the computation of carry in (B^4
    2) , \rho must be replaced by \lambda.

    INPUT:
    the three partitions from a cell in a growth diagram, labelled as

           shape1
    shape3 shape4

    OUTPUT:
    a pair (shape2, content) consisting of the shape of the fourth
    partition acording to Robinson-Schensted-Knuth and the content of
    the cell.

    """
    carry = 0
    i = len(shape4)
    shape2 = []
    while i > 0:
        if len(shape1) < i:
            row1 = 0
        else:
            row1 = shape1[i-1]
        if len(shape2) < i:
            row2 = 0
        else:
            row2 = shape2[i-1]
        if len(shape3) < i:
            row3 = 0
        else:
            row3 = shape3[i-1]
        shape2 = [min(row1, row3) - min(ZZ(row1 == row2 == row3), carry)] + shape2
        carry = carry - min(ZZ(row1 == row2 == row3), carry) + row4 - max(row1, row2)
        i = i-1
    return (Partition(shape2), carry)

def BinWord_forward(shape3, shape2, shape1, content):
    """See Fomin, "Schensted algorithms for dual graded graphs", Lemma
    4.6.1, page 40.

    INPUT: three binary words from a cell in a growth diagram,
    labelled as

    shape2 shape1     t x
    shape3            y

    and the content of the cell.

    OUTPUT:
    the fourth binary word shape4 according to Viennot's bijection.

    EXAMPLES:

    sage: G = GrowthDiagramBinWord([4,1,8,3,6,5,2,7,9]); G
    0  1  0  0  0  0  0  0  0
    0  0  0  0  0  0  1  0  0
    0  0  0  1  0  0  0  0  0
    1  0  0  0  0  0  0  0  0
    0  0  0  0  0  1  0  0  0
    0  0  0  0  1  0  0  0  0
    0  0  0  0  0  0  0  1  0
    0  0  1  0  0  0  0  0  0
    0  0  0  0  0  0  0  0  1
    sage: G._top_labels[9]
    word: 101010011

    """
    if shape1 == shape2 == shape3:
        if content == 0:
            r = shape1
        elif content == 1:
            r = Word(list(shape3) + [1], alphabet=[0,1])
        else:
            raise NotImplementedError
    elif shape1 == shape2:
        r = shape3
    elif shape3 == shape2:
        r = shape1
    else:
        if shape1 != shape3:
            r = Word(list(shape3) + [shape1[-1]], alphabet=[0,1])
        elif shape3 != shape2:
            r = Word(list(shape3) + [0], alphabet=[0,1])
        else:
            raise NotImplementedError
    return r

def BinWord_backward(y, z, x):
    """See Fomin, "Schensted algorithms for dual graded graphs", Lemma
    4.6.1, page 40.

    INPUT: three binary words from a cell in a growth diagram,
    labelled as

      x
    y z

    OUTPUT: a pair (t, content) consisting of the shape of the
    fourth word acording to Viennot and the content of the cell.

    """
    if x == y == z:
        return (x, 0)
    elif x == z != y:
        return (y, 0)
    elif x != z == y:
        return (x, 0)
    else:
        if x != y or (len(z) > 0 and z[-1] == 0):
            return (x[:-1], 0)
        elif x == y and len(z) > 0 and z[-1] == 1:
            return (x, 1)

class GrowthDiagram(SageObject):
    def __init__(self,
                 filling = None,
                 shape = None,
                 labels = None,
                 zero = None,
                 forward_rule = None,
                 backward_rule = None):
        """Initialise a generalized Schensted growth diagram in the sense of
        Fomin.

        An instance of the class is a growth diagram consisting of a
        filling and labels on the boundary.  To initialise it, it is
        necessary to provide either a filling and labels for the
        border on the side of the origin (forward growth), or labels
        for the boundary opposite of the origin (backward growth).

        Currently, this class only provides methods for dual graded
        graphs without multiple edges.

        Coordinates are of the form (col, row) where the origin is in
        the upper left, to be consistent with permutation matrices
        and skew tableaux.  This is different from Fomin's
        convention, who uses a Cartesian coordinate system.

        INPUT::

        - `filling` is `None`, if the growth diagram should be
          determined by applying the backward rules to the `labels`
          decorating the boundary opposite of the origin of the
          `shape`.  Otherwise it can be a dictionary with keys being
          coordinates and integer values, a sequence of sequences of
          integers (including matrices), or a word with integer
          letters (including permutations).

        - `shape` is a (possibly skew) partition or `None`.  In the
          latter case it is determined as the bounding rectangle
          (including the origin) of `filling`, or, if possible,
          inferred from `labels`.

        - `labels` is `None` or a list.  If it is a list, its length
          should be the half-perimeter of `shape`.  If `filling` is
          `None`, its elements are the labels of the boundary
          opposite of the origin.  Otherwise its elements are the
          labels on the boundary on the side of the origin.  If
          `labels` is `None` (in which case `filling` must not be
          `None`) the value of the parameter `zero` is used to
          initialise `labels`.

        EXAMPLES:

        sage: w = [3,3,2,4,1]; G = GrowthDiagramRSK(w)
        sage: from sage.combinat.growth import SkewTableau_from_partitions
        sage: [SkewTableau_from_partitions(G._top_labels[len(w):][::-1]), SkewTableau_from_partitions(G._top_labels[:len(w)+1])]
        [[[1, 3, 4], [2], [3]], [[1, 2, 4], [3], [5]]]
        sage: RSK(w)
        [[[1, 3, 4], [2], [3]], [[1, 2, 4], [3], [5]]]

        sage: p = Tableau([[1,2,2],[2]]); q = Tableau([[1,3,3],[2]])
        sage: gp = RSK_inverse(p, q); gp
        [[1, 2, 3, 3], [2, 1, 2, 2]]
        sage: GrowthDiagramRSK(labels = q.to_chain() + p.to_chain()[1::-1])
        0  1  0
        1  0  2

        """
        self._forward_rule = forward_rule
        self._backward_rule = backward_rule

        if filling is None:
            if shape is None:
                shape = self._shape_from_labels(labels)
                
            self._lambda, self._mu = self._init_shape_from_various_input(shape)
            self._top_labels = labels
            self._check_labels(self._top_labels)
            self._shrink()
        else:
            self._filling = self._init_filling_from_various_input(filling)
            self._lambda, self._mu = self._init_shape_from_various_input(shape)
            self._bot_labels = self._init_labels_forward_from_various_input(labels, zero)
            self._check_labels(self._bot_labels)
            self._grow()

    def _shape_from_labels(self, labels):
        """Return the shape of the growth diagram.

        This has to be implemented in the subclass, because it
        depends on the underlying graded graphs.

        The shape can be determined from the labels if the size of
        each label differs from the size of its successor.

        Otherwise raise an error.

        """
        raise NotImplementedError

    def __iter__(self):
        """
        Return the rows of the filling.

        TESTS:

        sage: G = GrowthDiagramRSK({(0,1):1, (1,0):1}, SkewPartition([[2,1],[1]]))
        sage: list(G)
        [[None, 1], [1]]

        sage: pi = Permutation([2,3,1,6,4,5])
        sage: G = GrowthDiagramRSK(pi.to_matrix())
        sage: list(G)
        [[0, 0, 1, 0, 0, 0],
         [1, 0, 0, 0, 0, 0],
         [0, 1, 0, 0, 0, 0],
         [0, 0, 0, 0, 1, 0],
         [0, 0, 0, 0, 0, 1],
         [0, 0, 0, 1, 0, 0]]
        """
        return ([None]*self._mu[r] + [self._filling.get((self._mu[r]+j,r), 0)
                                      for j in range(self._lambda[r]-self._mu[r])]
                for r in range(len(self._lambda)))

    def __repr__(self):
        """

        TESTS:

        sage: GrowthDiagramRSK({(0,1):1, (1,0):1}, SkewPartition([[2,1],[1]]))
        .  1
        1

        sage: GrowthDiagramRSK({(0,1):1, (2,0):1}, SkewPartition([[3,1],[1]]))
        .  0  1
        1
        """
        return SkewTableau(expr=[self._mu,
                                 [[self._filling.get((self._mu[r]+j,r), 0)
                                   for j in range(self._lambda[r]-self._mu[r])]
                                  for r in range(len(self._lambda))][::-1]])._repr_diagram()

    def _check_labels(self, labels):
        """Check sanity of the parameter `labels`.

        Assumes that `self._lambda` is already set.

        sage: GrowthDiagramRSK(shape=[1], labels=[[], [1]])
        Traceback (most recent call last):
        ...
        ValueError: the number of labels is 2, but for this shape we need 3.
        """
        if len(labels) != self._lambda[0]+len(self._lambda)+1:
            raise ValueError("the number of labels is %s, but for this shape we need %s."
                             %(len(labels), self._lambda[0]+len(self._lambda)+1))


    def _init_labels_forward_from_various_input(self, labels, zero):
        """Return a list of labels decorating the boundary near the origin.

        Assumes that `self._lambda` is already set.

        TESTS:

        `labels` is None:

        sage: filling = []
        sage: shape = SkewPartition([[4,2,1,1],[2,1,1]])
        sage: G = GrowthDiagramRSK(filling, shape)
        sage: G._bot_labels
        [[], [], [], [], [], [], [], [], []]

        `labels` is a list of partitions:

        sage: filling = []
        sage: labels = [[],[1],[],[1],[]]
        sage: shape = SkewPartition([[2,1],[1]])
        sage: G = GrowthDiagramRSK(filling=filling, shape=shape, labels=labels)
        sage: G._bot_labels
        [[], [1], [], [1], []]

        """
        if labels is None:
            return [zero for i in range(self._lambda[0]+len(self._lambda)+1)]
        else:
            return labels

    def _init_shape_from_various_input(self, shape):
        """Return a pair of partitions describing the region of the growth
        diagram.

        Assumes that `self._filling` is already set.

        TESTS:

        `shape` is a skew partition:

        sage: filling = []
        sage: shape = SkewPartition([[4,2,1,1],[2,1,1]])
        sage: G = GrowthDiagramRSK(filling, shape)
        sage: G._lambda, G._mu
        ([4, 2, 1, 1], [2, 1, 1, 0])

        `shape` is a partition:

        sage: filling = []
        sage: shape = Partition([3,2,1,1])
        sage: G = GrowthDiagramRSK(filling, shape)
        sage: G._lambda, G._mu
        ([3, 2, 1, 1], [0, 0, 0, 0])

        `shape` is a list:

        sage: pi = Permutation([2,3,1,6,4,5])
        sage: G = GrowthDiagramRSK(pi.to_matrix(), [len(pi)]*6)
        sage: G._lambda, G._mu
        ([6, 6, 6, 6, 6, 6], [0, 0, 0, 0, 0, 0])

        `shape` is None:

        sage: filling = {(1,3):1, (2,1):1, (1,1):1}
        sage: G = GrowthDiagramRSK(filling)
        sage: G._lambda, G._mu
        ([3, 3, 3, 3], [0, 0, 0, 0])

        """
        if shape is None:
            assert self._filling is not None, "The shape of the diagram must be given."
            # find bounding rectangle of `filling`
            max_row = max(i for i, _ in self._filling)+1
            max_col = max(j for _, j in self._filling)+1
            return ([max_row]*max_col, [0]*max_col)
        else:
            try:
                shape = Partition(shape)
                return (list(shape), [0]*len(shape))
            except ValueError:
                try:
                    shape = SkewPartition(shape)
                    return (list(shape[0]),
                            list(shape[1]) + [0]*(len(shape[0])-len(shape[1])))
                except ValueError:
                    raise ValueError("Cannot make sense of shape %s" %shape)


    def _init_filling_from_various_input(self, filling):
        """Return a dict `F`, such that `F[(i,j)]` is the element in row `i`
        and column `j`.

        TESTS:

        `filling` is a dict of coordinates:

        sage: pi = Permutation([2,3,1,6,4,5])
        sage: G = GrowthDiagramRSK({(i,pi[i]-1):1 for i in range(len(pi))})
        sage: G._filling
        {(0, 1): 1, (1, 2): 1, (2, 0): 1, (3, 5): 1, (4, 3): 1, (5, 4): 1}

        `filling` is a dict of dicts:

        sage: G = GrowthDiagramRSK({i:{pi[i]-1:1} for i in range(len(pi))})
        sage: G._filling
        {(0, 1): 1, (1, 2): 1, (2, 0): 1, (3, 5): 1, (4, 3): 1, (5, 4): 1}

        `filling` is a matrix:

        sage: G = GrowthDiagramRSK(pi.to_matrix())
        sage: G._filling
        {(0, 1): 1, (1, 2): 1, (2, 0): 1, (3, 5): 1, (4, 3): 1, (5, 4): 1}

        `filling` is a permutation:

        sage: G = GrowthDiagramRSK(pi)
        sage: G._filling
        {(0, 1): 1, (1, 2): 1, (2, 0): 1, (3, 5): 1, (4, 3): 1, (5, 4): 1}

        `filling` is a word:

        sage: G = GrowthDiagramRSK([3,1,4,1,5])
        sage: G._filling
        {(0, 2): 1, (1, 0): 1, (2, 3): 1, (3, 0): 1, (4, 4): 1}

        """
        if isinstance(filling, dict):
            v = filling.itervalues().next()
            if isinstance(v, dict):
                # it is a dict of dicts
                F = dict()
                for (i, row) in filling.iteritems():
                    for (j, v) in row.iteritems():
                        if v != 0:
                            F[(i,j)] = int(v)
                return F
            else:
                # it is dict of coordinates
                return filling

        else:
            # it is a sequence
            F = dict()
            try:
                # it is a sequence of sequences
                for i, row in enumerate(filling):
                    for j, v in enumerate(row):
                        if v != 0:
                            F[(j,i)] = int(v)
            except TypeError:
                # it is a word
                for i, l in enumerate(filling):
                    F[(i, l-1)] = 1

            return F

    def _grow(self):
        """Compute the labels on the boundary opposite of the origin, given
        the filling.

        TESTS:

        sage: pi = Permutation([1])
        sage: G = GrowthDiagramRSK(pi.to_matrix())
        sage: G._top_labels
        [[], [1], []]

        sage: pi = Permutation([1,2])
        sage: G = GrowthDiagramRSK(pi.to_matrix())
        sage: G._top_labels
        [[], [1], [2], [1], []]

        sage: pi = Permutation([2,1])
        sage: G = GrowthDiagramRSK(pi.to_matrix())
        sage: G._top_labels
        [[], [1], [1, 1], [1], []]

        sage: pi = Permutation([2,3,6,1,4,5])
        sage: G = GrowthDiagramRSK(pi.to_matrix())
        sage: G._top_labels
        [[],
         [1],
         [2],
         [3],
         [3, 1],
         [3, 2],
         [4, 2],
         [4, 1],
         [3, 1],
         [2, 1],
         [1, 1],
         [1],
         []]

        sage: G = GrowthDiagramRSK({(0,1):1, (1,0):1}, SkewPartition([[2,1],[1]]))
        sage: G._top_labels
        [[], [1], [], [1], []]

        sage: G = GrowthDiagramRSK({(1,1):1}, SkewPartition([[2,2],[1]]), labels=[[],[],[1],[],[]])
        sage: G._top_labels
        [[], [1], [2], [1], []]
        """
        labels = copy(self._bot_labels)
        l = len(self._lambda)
        for r in range(l):
            for c in range(self._mu[r]+l-r, self._lambda[r]+l-r):
                j = r
                i = c-l+r
                labels[c] = self._forward_rule(labels[c-1],
                                               labels[c],
                                               labels[c+1],
                                               self._filling.get((i,j), 0))

        self._top_labels = labels

    def _shrink(self):
        """TESTS:

        sage: filling = [[0,0,1,0,0,0,0], [0,1,0,0,0,0,0], [1,0,0,0,0,0,0], [0,0,0,1,0,0,0], [0,0,0,0,0,0,1], [0,0,0,0,0,1,0], [0,0,0,0,1,0,0]]
        sage: G = GrowthDiagramRSK(filling)
        sage: list(GrowthDiagramRSK(labels=G._top_labels)) == filling
        True

        sage: labels = [[], [1], []]
        sage: G = GrowthDiagramRSK(labels=labels)
        sage: G._filling
        {(0, 0): 1}
        sage: G._bot_labels
        [[], [], []]

        sage: labels = [[], [1], [2], [2,1], [1,1], [1], []]
        sage: G = GrowthDiagramRSK(labels=labels)
        sage: G._filling
        {(0, 1): 1, (1, 2): 1, (2, 0): 1}
        sage: G._bot_labels
        [[], [], [], [], [], [], []]

        sage: labels = [[], [1], [2], [3], [3, 1], [3, 2], [4, 2], [4, 1], [3, 1], [2, 1], [1, 1], [1], []]
        sage: G = GrowthDiagramRSK(labels=labels)
        sage: G._filling
        {(0, 1): 1, (1, 2): 1, (2, 5): 1, (3, 0): 1, (4, 3): 1, (5, 4): 1}

        sage: labels = [[],[1],[1],[2],[2],[2,1],[2]]
        sage: G = GrowthDiagramRSK(labels=labels)
        Traceback (most recent call last):
        ...
        ValueError: Can only determine the shape of the growth diagram if sizes of successive partitions differ.
        sage: G = GrowthDiagramRSK(shape=[3,2,1], labels=labels)
        sage: G._filling
        {(1, 0): 1}
        sage: G._bot_labels
        [[], [], [], [], [1], [1], [2]]

        sage: labels = [[], [1],[1],[2],[2],[2,1],[2],[2,1],[1,1],[2,1],[1,1]]
        sage: G = GrowthDiagramRSK(shape=[5,4,3,2,1], labels=labels)
        sage: G._filling
        {(1, 2): 1, (2, 1): 1, (4, 0): 1}
        sage: G._bot_labels
        [[], [], [], [], [], [], [1], [1], [1], [1, 1], [1, 1]]

        sage: labels = [[], [1],[1],[2],[2],[2,1],[2],[2,1],[1,1],[2,1],[1,1]]
        sage: G = GrowthDiagramRSK(shape=SkewPartition([[5,4,3,2,1],[3,2,1]]), labels=labels)
        sage: G._filling
        {(1, 2): 1, (2, 1): 1, (4, 0): 1}
        sage: G._bot_labels
        [[], [], [], [1], [1], [1], [1], [1], [1], [1, 1], [1, 1]]
        """
        F = dict()
        labels = copy(self._top_labels)
        l = len(self._lambda)
        for r in range(l):
            for c in range(self._lambda[l-r-1]+r, self._mu[l-r-1]+r, -1):
                j = l-r-1
                i = c-r-1
                labels[c], v = self._backward_rule(labels[c-1],
                                                   labels[c],
                                                   labels[c+1])
                if v != 0:
                    F[(i,j)] = v

        self._bot_labels = labels
        self._filling = F


class GrowthDiagramBinWord(GrowthDiagram):
    def __init__(self,
                 filling = None,
                 shape = None,
                 labels = None):
        super(GrowthDiagramBinWord, self).__init__(filling = filling,
                                                   shape = shape,
                                                   labels = labels,
                                                   zero = Word([], alphabet=[0,1]),
                                                   forward_rule=BinWord_forward,
                                                   backward_rule=BinWord_backward)

    def _shape_from_labels(self, labels):
        def right_left(la, mu):
            if len(la) < len(mu):
                return 1
            elif len(la) > len(mu):
                return 0
            else:
                raise ValueError("Can only determine the shape of the growth diagram if sizes of successive words differ.")
        return Partitions().from_zero_one([right_left(labels[i], labels[i+1]) for i in range(len(labels)-1)])

        
class GrowthDiagramRSK(GrowthDiagram):
    def __init__(self,
                 filling = None,
                 shape = None,
                 labels = None):
        if labels is not None:
            labels = [Partition(la) for la in labels]
        super(GrowthDiagramRSK, self).__init__(filling = filling,
                                               shape = shape,
                                               labels = labels,
                                               zero = Partition([]),
                                               forward_rule = Robinson_Schensted_Knuth_forward,
                                               backward_rule = Robinson_Schensted_Knuth_backward)

    def _shape_from_labels(self, labels):
        def right_left(la, mu):
            if la.size() < mu.size():
                return 1
            elif la.size() > mu.size():
                return 0
            else:
                raise ValueError("Can only determine the shape of the growth diagram if sizes of successive partitions differ.")
        return Partitions().from_zero_one([right_left(labels[i], labels[i+1]) for i in range(len(labels)-1)])
        
#        Domino(shape1, shape2, shape3, content) ==
#
#            if not  member?(content,[0,1,2]$List(NonNegativeInteger)) then
#                error "Domino: The content of the filling must be in {0,1,2}"
#
#            (content = 1) =>
#                l1 := copy(shape1::List(Integer))
#                if empty?(l1) then return partition([2])$Partition
#                else
#                    l1.1:=l1.1 + 2
#                    return partition(l1)$Partition
#
#            (content = 2) =>
#                  l1 := copy(shape1::List(Integer))
#                  if empty?(l1) then return partition([1,1])$Partition
#                  else if #l1=1 then
#                      l1.1:=l1.1 + 1
#                      return partition(append(l1, [1]))$Partition
#                  else
#                    l1.1:=l1.1 + 1
#                    l1.2:=l1.2 + 1
#                    return partition(l1)$Partition
#
#            (content=0 and (shape2=shape1 or shape2=shape3)) =>
#                return union(shape1,shape3)
#
#            (content=0 and shape2~=shape1 and shape1=shape3) =>
#                l1 := copy(shape1::List(Integer))
#                l2 := shape2::List(Integer)
#                i := 1
#                j :=#l2
#                while i<=j  and l1.i = l2.i repeat
#                    i := i+1
#                if j<i then
#                    if l1.i=2 then return partition(append(l1, [2]))$Partition
#                    else
#                        l1.i     :=  l1.i + 1
#                        l1.(i+1) := l1.(i+1) + 1
#                        return partition(l1)$Partition
#                else
#                    if  l1.i=l2.i + 2 then return partition(append(l1, [2]))$Partition
#                    else
#                        l1.i:= l1.i + 1
#                        l1.(i+1) := l1.(i+1) + 1
#                        return partition(l1)$Partition
#
#
#            (content=0 and shape2~=shape1 and shape1~=shape3) =>
#                l1 := copy(shape1::List(Integer))
#                l2 := shape2::List(Integer)
#                l3 := copy(shape3::List(Integer))
#                i := 1
#                j :=#l2
#                shape := union(union(shape2,shape1),shape3)
#                while i<=j  and l1.i = l2.i repeat
#                    i := i+1
#                if  intersect(shape1,shape3) ~= shape2 then
#                    list  := shape::List(Integer)
#                    list.(i+1) := list.(i+1)+1
#                    return (partition(list)$Partition)
#                else return shape
