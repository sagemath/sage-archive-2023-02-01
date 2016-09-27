r"""
Growth diagrams and dual graded graphs.

AUTHORS:

- Martin Rubey (2016-09): Initial version


.. TODO::

    - when shape is given, check that it is compatible with filling or labels
    - implement domino insertion
    - implement Young-Fibonacci

REFERENCES:

.. [Fom1995] Sergey V. Fomin
   *Schensted algorithms for dual graded graphs*.
   Journal of Algebraic Combinatorics Volume 4, Number 1 (1995), pp. 5-45

.. [Kra2006] Christian Krattenthaler
   *Growth diagrams, and increasing and decreasing chains in fillings of Ferrers shapes*.
   Advances in Applied Mathematics Volume 37, Number 3 (2006), pp. 404-431

.. [Lam2004] Thomas Lam
   *Growth diagrams, domino insertion and sign-imbalance*.
   Journal of Combinatorial Theory, Series A Volume 107, Number 1 (2004), pp. 87-115
"""
from sage.structure.sage_object import SageObject
from sage.combinat.words.word import Word
from sage.combinat.partition import Partition, Partitions
from sage.combinat.skew_partition import SkewPartition
from sage.combinat.skew_tableau import SkewTableau
from copy import copy

"""
local rules
-----------
RSK
^^^
"""
def Robinson_Schensted_Knuth_forward(shape3, shape2, shape1, content):
    """Return the output shape given three shapes and the content.

    See [Kra2006]_ ``(F^1 0)-(F^1 2)``.

    INPUT:

    - ``shape3, shape2, shape1`` -- three partitions from a cell in a
      growth diagram, labelled as::

          shape2 shape1
          shape3

    - ``content`` -- a non-negative integer, the content of the cell.

    OUTPUT:

    The fourth partition according to the Robinson-Schensted-Knuth
    correspondence.

    EXAMPLES::

        sage: pi = Permutation([2,3,6,1,4,5])
        sage: G = GrowthDiagramRSK(pi)
        sage: G._out_labels
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
        sage: P, Q = RSK(pi); Q.to_chain()[:-1] + P.to_chain()[::-1] == G._out_labels
        True

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
            # returning this as a Partition costs a lot of time
            return shape4[::-1]
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
    """Return the content and the input shape.

    See [Kra2006]_ ``(B^1 0)-(B^1 2)``.

    INPUT:

    - ``shape3, shape4, shape1`` -- three partitions from a cell in a
    growth diagram, labelled as::

               shape1
        shape3 shape4

    OUTPUT:

    A pair ``(shape2, content)`` consisting of the shape of the
    fourth word acording to the Robinson-Schensted-Knuth
    correspondence and the content of the cell.

    TEST::

        sage: w = [4,1,8,3,6,5,2,7,9]; G = GrowthDiagramRSK(w);
        sage: G._out_labels
        [[],
         [1],
         [1, 1],
         [2, 1],
         [2, 2],
         [3, 2],
         [3, 2, 1],
         [3, 2, 1, 1],
         [4, 2, 1, 1],
         [5, 2, 1, 1],
         [4, 2, 1, 1],
         [4, 2, 1],
         [3, 2, 1],
         [3, 1, 1],
         [2, 1, 1],
         [2, 1],
         [2],
         [1],
         []]
        sage: GrowthDiagramRSK(labels=G._out_labels).to_word() == w
        True

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

"""
Burge
^^^^^
"""
def Burge_forward(shape3, shape2, shape1, content):
    """Return the output shape given three shapes and the content.

    See [Kra2006]_ ``(F^4 0)-(F^4 2)``.

    INPUT:

    - ``shape3, shape2, shape1`` -- three  from a cell in a growth diagram,
      labelled as::

          shape2 shape1
          shape3

    - ``content`` -- a non-negative integer, the content of the cell.

    OUTPUT:

    The fourth partition according to the Burge correspondence.

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
    """Return the content and the input shape.

    See [Kra2006]_ ``(B^4 0)-(B^4 2)``.  There is a typo in the
    computation of carry in ``(B^4 2)``, `\rho` must be replaced by
    `\lambda`.

    INPUT:

    - ``shape3, shape4, shape1`` -- three partitions from a cell in a
      growth diagram, labelled as::

               shape1
        shape3 shape4

    OUTPUT:

    A pair ``(t, content)`` consisting of the shape of the fourth
    partition acording to the Burge correspondence and the content of
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

"""
BinWord - BinTree
^^^^^^^^^^^^^^^^^
"""
def BinWord_forward(shape3, shape2, shape1, content):
    """Return the output shape given three shapes and the content.

    See [Fom1995]_ Lemma 4.6.1, page 40.

    INPUT:

    - ``shape3, shape2, shape1`` -- three binary words from a cell in a growth diagram,
      labelled as::

          shape2 shape1
          shape3

    - ``content`` -- 0 or 1, the content of the cell.

    OUTPUT:

    The fourth binary word shape4 according to Viennot's bijection.

    EXAMPLES::

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
        sage: G._out_labels[9]
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
    """Return the content and the input shape.

    See [Fom1995]_ Lemma 4.6.1, page 40.

    - ``y, z, x`` -- three binary words from a cell in a growth diagram,
      labelled as::

            x
          y z

    OUTPUT:

    A pair ``(t, content)`` consisting of the shape of the fourth
    word acording to Viennot's bijection and the content of the cell.

    TEST::

        sage: w = [4,1,8,3,6,5,2,7,9]; G = GrowthDiagramBinWord(w);
        sage: G._out_labels
        [word: ,
          word: 1,
          word: 10,
          word: 101,
          word: 1010,
          word: 10101,
          word: 101010,
          word: 1010100,
          word: 10101001,
          word: 101010011,
          word: 10101001,
          word: 1001001,
          word: 100100,
          word: 10010,
          word: 1000,
          word: 110,
          word: 11,
          word: 1,
          word: ]
        sage: GrowthDiagramBinWord(labels=G._out_labels).to_word() == w
        True

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

"""
Domino Robinson-Schensted
^^^^^^^^^^^^^^^^^^^^^^^^^
"""
def Domino_forward(shape3, shape2, shape1, content):
    """Return the output shape given three shapes and the content.

    See [Lam2004]_ Section 3.1.

    INPUT:

    - ``shape3, shape2, shape1`` -- three partitions from a cell in a
      growth diagram, labelled as::

          shape2 shape1
          shape3

    - ``content`` -- -1, 0 or 1, the content of the cell.

    OUTPUT:

    The fourth partition according to Domino insertion.

    EXAMPLES::

        sage: from sage.combinat.growth import Domino_forward

        Rule 1:
        sage: Domino_forward([], [], [], 1)
        [2]

        sage: Domino_forward([1,1], [1,1], [1,1], 1)
        [3, 1]

        Rule 2:
        sage: Domino_forward([1,1], [1,1], [1,1], -1)
        [2, 2]

        Rule 3:
        sage: Domino_forward([1,1], [1,1], [2,2], 0)
        [2, 2]

        Rule 4:
        sage: Domino_forward([2,2,2], [2,2], [3,3], 0)
        [3, 3, 2]

        sage: Domino_forward([2], [], [1,1], 0)
        [2, 2]

        sage: Domino_forward([1,1], [], [2], 0)
        [2, 2]

        sage: Domino_forward([2], [], [2], 0)
        [2, 2]

        sage: Domino_forward([4], [2], [4], 0)
        [4, 2]

        sage: Domino_forward([1,1,1,1], [1,1], [1,1,1,1], 0)
        [2, 2, 1, 1]

    """
    if content not in [0,1,-1]:
        raise ValueError("Domino: The content of the filling must be in {-1,0,1}")

    if content == 1:
        if shape2 == []:
            shape4 = [2]
        else:
            shape4 = [shape2[0] + 2] + shape2[1:]

    elif content == -1:
        if shape2 == []:
            shape4 = [1,1]
        elif len(shape2) == 1:
            shape4 = [shape2[0] + 1] + [1]
        else:
            shape4 = [shape2[0] + 1, shape2[1] + 1] + shape2[2:]

    elif content == 0 and (shape2 == shape1 or shape2 == shape3):
        if len(shape1) > len(shape3):
            shape4 = shape1
        elif len(shape1) < len(shape3):
            shape4 = shape3
        else:
            shape4 = [max(p,q) for (p,q) in zip(shape1, shape3)]

    else:
        # content == 0 and shape2 differs from shape1 and shape3 by
        # dominos gamma1 and gamma3

        # the following is certainly very slow
        gamma3 = set(SkewPartition([shape3, shape2]).cells())
        gamma1 = set(SkewPartition([shape1, shape2]).cells())

        diff = gamma1.intersection(gamma3)
        cell1, cell2 = gamma3
        shape4 = copy(shape1)

        if len(diff) == 0:
            # add gamma3 to shape1
            if len(shape4) <= cell1[0]:
                shape4 += [1]
            else:
                shape4[cell1[0]] += 1
            if len(shape4) <= cell2[0]:
                shape4 += [1]
            else:
                shape4[cell2[0]] += 1

        elif len(diff) == 1:
            # diff is a single cell
            (k,l) = diff.pop()
            # add (k+1, l+1) to shape1
            # either (k, l+1) or (k+1, l) must also be added
            if shape4[k] <= l+1:
                shape4[k] += 1
                shape4[k+1] += 1
            else:
                if len(shape4) <= k+1:
                    shape4 += [2]
                else:
                    shape4[k+1] += 2

        # diff has size 2, that is shape1 == shape3

        elif cell1[0] == cell2[0]:
            # a horizontal domino
            if len(shape4) <= cell1[0]+1:
                shape4 += [2]
            else:
                shape4[cell1[0]+1] += 2
        else:
            # a vertical domino
            # find first row shorter than cell1[1]
            for r, p in enumerate(shape4):
                if p <= cell1[1]:
                    shape4[r] += 1
                    shape4[r+1] += 1
            else:
                shape4[0] += 1
                shape4[1] += 1

    return shape4

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

        INPUT:

        - `filling` is `None`, if the growth diagram should be
          determined by applying the backward rules to the `labels`
          decorating the boundary opposite of the origin of the
          `shape`.  Otherwise it can be a dictionary with keys being
          coordinates and integer values, a sequence of sequences of
          integers (including matrices), or a word with integer
          letters (including permutations).

        - `shape` is a (possibly skew) partition or `None`.  In the
          latter case it is determined as the Ferrers shape given by
          `filling`, if that is a sequence of sequences, the bounding
          rectangle (including the origin) of `filling`, or, if
          possible, inferred from `labels`, if `filling` is `None`.

        - `labels` is `None` or a list.  If it is a list, its length
          should be the half-perimeter of `shape`.  If `filling` is
          `None`, its elements are the labels of the boundary
          opposite of the origin.  Otherwise its elements are the
          labels on the boundary on the side of the origin.  If
          `labels` is `None` (in which case `filling` must not be
          `None`) the value of the parameter `zero` is used to
          initialise `labels`.

        EXAMPLES::

            sage: w = [3,3,2,4,1]; G = GrowthDiagramRSK(w)
            sage: [G.P_symbol(), G.Q_symbol()]
            [[[1, 3, 4], [2], [3]], [[1, 2, 4], [3], [5]]]
            sage: RSK(w)
            [[[1, 3, 4], [2], [3]], [[1, 2, 4], [3], [5]]]

        """
        self._forward_rule = forward_rule
        self._backward_rule = backward_rule

        if filling is None:
            if shape is None:
                shape = self._shape_from_labels(labels)

            self._lambda, self._mu = self._init_shape_from_various_input(shape)
            self._out_labels = labels
            self._check_labels(self._out_labels)
            self._shrink()
        else:
            self._filling, (self._lambda, self._mu) = self._init_filling_and_shape_from_various_input(filling, shape)
            self._in_labels = self._init_labels_forward_from_various_input(labels, zero)
            self._check_labels(self._in_labels)
            self._grow()

    def _shape_from_labels(self, labels):
        """Determine the shape of the growth diagram given a list of labels
        during initialisation.

        This has to be implemented in the subclass, because it
        depends on the underlying graded graphs.

        The shape can be determined from the labels if the size of
        each label differs from the size of its successor.

        Otherwise raise an error.

        """
        raise NotImplementedError

    def conjugate(self):
        """Return the growth diagram with the filling reflected on the main
        diagonal.

        When the filling is a permutation, the conjugate filling corresponds to its inverse.

        EXAMPLES::

            sage: G = GrowthDiagramRSK([[0,1,0], [1,0,2]])
            sage: Gc = G.conjugate()
            sage: (Gc.P_symbol(), Gc.Q_symbol()) == (G.Q_symbol(), G.P_symbol())
            True
        """
        F = {(j,i): v for (i,j),v in self._filling.iteritems()}
        return self.parent()(filling = F)

    def rotate(self):
        """Return the growth diagram with the filling rotated by 180 degrees.

        For RSK-growth diagrams and rectangular fillings, this
        corresponds to evacutation of the P- and the Q-symbol.

        EXAMPLES::

            sage: G = GrowthDiagramRSK([[0,1,0], [1,0,2]])
            sage: Gc = G.rotate()
            sage: ascii_art(Gc.P_symbol(), Gc.Q_symbol())
            1  1  1  1  1  2
            2        3

            sage: ascii_art(SemistandardTableau(Tableau(G.P_symbol())).evacuation(), SemistandardTableau(Tableau(G.Q_symbol())).evacuation())
            1  1  1  1  1  2
            2        3

        """
        max_row = max(i for i, _ in self._filling)
        max_col = max(j for _, j in self._filling)
        F = {(max_row-i,max_col-j): v for (i,j),v in self._filling.iteritems()}
        return self.parent()(filling = F)

    def shape(self):
        """
        Return the shape of the growth diagram as a skew partition.

        EXAMPLES::

            sage: GrowthDiagramRSK([1]).shape()
            [1] / []

        """
        return SkewPartition([self._lambda, self._mu])

    def P_symbol(self):
        """Return the labels along the vertical boundary of a rectangular
        growth diagram as a skew tableau.

        EXAMPLES::

            sage: G = GrowthDiagramRSK([[0,1,0], [1,0,2]])
            sage: G.P_symbol().pp()
            1  2  2
            2
        """
        if self.is_rectangular():
            return SkewTableau(chain = self._out_labels[self._lambda[0]:][::-1])
        else:
            raise ValueError("The P symbol is only defined for rectangular shapes.")

    def Q_symbol(self):
        """ Return the labels along the horizontal boundary of a rectangular
        growth diagram as a skew tableau.

        EXAMPLES::

            sage: G = GrowthDiagramRSK([[0,1,0], [1,0,2]])
            sage: G.Q_symbol().pp()
            1  3  3
            2
        """
        if self.is_rectangular():
            return SkewTableau(chain = self._out_labels[:self._lambda[0]+1])
        else:
            raise ValueError("The Q symbol is only defined for rectangular shapes.")

    def is_rectangular(self):
        """
        Return ``True`` if the shape of the growth diagram is rectangular.

        EXAMPLES::

            sage: GrowthDiagramRSK([2,3,1]).is_rectangular()
            True
            sage: GrowthDiagramRSK([[1,0,1],[0,1]]).is_rectangular()
            False
        """
        return all(x == 0 for x in self._mu) and all(x == self._lambda[0] for x in self._lambda)

    def to_word(self):
        """Return the filling as a word, if the shape is rectangular and
        there is at most one nonzero entry in each column, which must be 1.

        EXAMPLES::

            sage: w = [3,3,2,4,1]; G = GrowthDiagramRSK(w)
            sage: G
            0  0  0  0  1
            0  0  1  0  0
            1  1  0  0  0
            0  0  0  1  0
            sage: G.to_word()
            [3, 3, 2, 4, 1]

        """
        if self.is_rectangular():
            w = [0]*self._lambda[0]
            for ((i,j), v) in self._filling.iteritems():
                if v != 0:
                    if v == 1:
                        if w[i] == 0:
                            w[i] = j+1
                        else:
                            raise ValueError("Can only convert fillings with at most one entry per column to words.")
                    else:
                        raise ValueError("Can only convert 0-1 fillings to words.  Try `to_biword`.")
            return w
        else:
            raise ValueError("Can only convert fillings of rectangular shapes to words.")

    def to_biword(self):
        """Return the filling as a biword, if the shape is rectangular.

        EXAMPLES::

            sage: P = Tableau([[1,2,2],[2]]); Q = Tableau([[1,3,3],[2]])
            sage: bw = RSK_inverse(P, Q); bw
            [[1, 2, 3, 3], [2, 1, 2, 2]]
            sage: G = GrowthDiagramRSK(labels = Q.to_chain()[:-1] + P.to_chain()[::-1]); G
            0  1  0
            1  0  2

            sage: P = SemistandardTableau([[1, 1, 2], [2]])
            sage: Q = SemistandardTableau([[1, 2, 2], [2]])
            sage: G = GrowthDiagramRSK(labels = Q.to_chain()[:-1] + P.to_chain()[::-1]); G
            0  2
            1  1
            sage: G.to_biword()
            ([1, 2, 2, 2], [2, 1, 1, 2])
            sage: RSK([1, 2, 2, 2], [2, 1, 1, 2])
            [[[1, 1, 2], [2]], [[1, 2, 2], [2]]]

        """
        if self.is_rectangular():
            w1 = []
            w2 = []
            for ((i,j), v) in sorted(self._filling.iteritems()):
                if v >= 0:
                    w1.extend([i+1]*v)
                    w2.extend([j+1]*v)
                else:
                    raise ValueError("Can only convert fillings with nonnegative entries to words.")
            return (w1, w2)
        else:
            raise ValueError("Can only convert fillings of rectangular shapes to words.")

    def __iter__(self):
        """
        Return the rows of the filling.

        TESTS::

            sage: G = GrowthDiagramRSK({(0,1):1, (1,0):1}, SkewPartition([[2,1],[1]]))
            sage: list(G)
            [[None, 1], [1]]

            sage: pi = Permutation([2,3,1,6,4,5])
            sage: G = GrowthDiagramRSK(pi)
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
        """Print the filling of the growth diagram as a skew tableau.

        TESTS::

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

        TEST::

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

        `labels` is None::

            sage: filling = []
            sage: shape = SkewPartition([[4,2,1,1],[2,1,1]])
            sage: G = GrowthDiagramRSK(filling, shape)
            sage: G._in_labels
            [[], [], [], [], [], [], [], [], []]

        `labels` is a list of partitions::

            sage: filling = []
            sage: labels = [[],[1],[],[1],[]]
            sage: shape = SkewPartition([[2,1],[1]])
            sage: G = GrowthDiagramRSK(filling=filling, shape=shape, labels=labels)
            sage: G._in_labels
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

        `shape` is a skew partition::

            sage: filling = []
            sage: shape = SkewPartition([[4,2,1,1],[2,1,1]])
            sage: G = GrowthDiagramRSK(filling, shape)
            sage: G._lambda, G._mu
            ([4, 2, 1, 1], [2, 1, 1, 0])

        `shape` is a partition::

            sage: filling = []
            sage: shape = Partition([3,2,1,1])
            sage: G = GrowthDiagramRSK(filling, shape)
            sage: G._lambda, G._mu
            ([3, 2, 1, 1], [0, 0, 0, 0])

        """
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


    def _init_filling_and_shape_from_various_input(self, filling, shape):
        """Return a dict `F`, such that `F[(i,j)]` is the element in row `i`
        and column `j`, and a pair of partitions describing the
        region of the growth diagram.

        TESTS:

        `filling` is a dict of coordinates::

            sage: pi = Permutation([2,3,1,6,4,5])
            sage: G = GrowthDiagramRSK({(i,pi[i]-1):1 for i in range(len(pi))})
            sage: G._filling
            {(0, 1): 1, (1, 2): 1, (2, 0): 1, (3, 5): 1, (4, 3): 1, (5, 4): 1}
            sage: G.shape()
            [6, 6, 6, 6, 6, 6] / []

        `filling` is a dict of dicts::

            sage: G = GrowthDiagramRSK({i:{pi[i]-1:1} for i in range(len(pi))})
            sage: G._filling
            {(0, 1): 1, (1, 2): 1, (2, 0): 1, (3, 5): 1, (4, 3): 1, (5, 4): 1}
            sage: G.shape()
            [6, 6, 6, 6, 6, 6] / []

        `filling` is a matrix::

            sage: G = GrowthDiagramRSK(pi.to_matrix())
            sage: G._filling
            {(0, 1): 1, (1, 2): 1, (2, 0): 1, (3, 5): 1, (4, 3): 1, (5, 4): 1}
            sage: G.shape()
            [6, 6, 6, 6, 6, 6] / []

        `filling` is a permutation::

            sage: G = GrowthDiagramRSK(pi)
            sage: G._filling
            {(0, 1): 1, (1, 2): 1, (2, 0): 1, (3, 5): 1, (4, 3): 1, (5, 4): 1}
            sage: G.shape()
            [6, 6, 6, 6, 6, 6] / []

        `filling` is a list::

            sage: G = GrowthDiagramRSK([3,1,4,1,5])
            sage: G._filling
            {(0, 2): 1, (1, 0): 1, (2, 3): 1, (3, 0): 1, (4, 4): 1}
            sage: G.shape()
            [5, 5, 5, 5, 5] / []

        `filling` is a list of lists::

            sage: G = GrowthDiagramRSK([[1,0,1],[0,1]])
            sage: G._filling
            {(0, 0): 1, (1, 1): 1, (2, 0): 1}
            sage: G.shape()
            [3, 2] / []

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
            else:
                # it is dict of coordinates
                F = filling

        else:
            # it is a sequence
            F = dict()
            try:
                # it is a sequence of sequences
                for i, row in enumerate(filling):
                    for j, v in enumerate(row):
                        if v != 0:
                            F[(j,i)] = int(v)
                if shape is None:
                    shape = [len(row) for row in filling]

            except TypeError:
                # it is a word
                for i, l in enumerate(filling):
                    F[(i, l-1)] = 1

        if shape is None:
            # find bounding rectangle of `filling`
            max_row = max(i for i, _ in F)+1
            max_col = max(j for _, j in F)+1
            shape = [max_row]*max_col


        return (F, self._init_shape_from_various_input(shape))

    def _grow(self):
        """Compute the labels on the boundary opposite of the origin, given
        the filling.

        TESTS::

            sage: pi = Permutation([1])
            sage: G = GrowthDiagramRSK(pi)
            sage: G._out_labels
            [[], [1], []]

            sage: pi = Permutation([1,2])
            sage: G = GrowthDiagramRSK(pi)
            sage: G._out_labels
            [[], [1], [2], [1], []]

            sage: pi = Permutation([2,1])
            sage: G = GrowthDiagramRSK(pi)
            sage: G._out_labels
            [[], [1], [1, 1], [1], []]

            sage: G = GrowthDiagramRSK({(0,1):1, (1,0):1}, SkewPartition([[2,1],[1]]))
            sage: G._out_labels
            [[], [1], [], [1], []]

            sage: G = GrowthDiagramRSK({(1,1):1}, SkewPartition([[2,2],[1]]), labels=[[],[],[1],[],[]])
            sage: G._out_labels
            [[], [1], [2], [1], []]
        """
        assert self._forward_rule is not None, "For computing a growth diagram, a forward rule has to be specified."
        labels = copy(self._in_labels)
        l = len(self._lambda)
        for r in range(l):
            for c in range(self._mu[r]+l-r, self._lambda[r]+l-r):
                j = r
                i = c-l+r
                labels[c] = self._forward_rule(labels[c-1],
                                               labels[c],
                                               labels[c+1],
                                               self._filling.get((i,j), 0))

        self._out_labels = labels

    def _shrink(self):
        """Compute the labels on the boundary near the origin, and the filling.

        TESTS::

            sage: filling = [[0,0,1,0,0,0,0], [0,1,0,0,0,0,0], [1,0,0,0,0,0,0], [0,0,0,1,0,0,0], [0,0,0,0,0,0,1], [0,0,0,0,0,1,0], [0,0,0,0,1,0,0]]
            sage: G = GrowthDiagramRSK(filling)
            sage: list(GrowthDiagramRSK(labels=G._out_labels)) == filling
            True

            sage: labels = [[], [1], []]
            sage: G = GrowthDiagramRSK(labels=labels)
            sage: G._filling
            {(0, 0): 1}
            sage: G._in_labels
            [[], [], []]

            sage: labels = [[], [1], [2], [2,1], [1,1], [1], []]
            sage: G = GrowthDiagramRSK(labels=labels)
            sage: G._filling
            {(0, 1): 1, (1, 2): 1, (2, 0): 1}
            sage: G._in_labels
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
            sage: G._in_labels
            [[], [], [], [], [1], [1], [2]]

            sage: labels = [[], [1],[1],[2],[2],[2,1],[2],[2,1],[1,1],[2,1],[1,1]]
            sage: G = GrowthDiagramRSK(shape=[5,4,3,2,1], labels=labels)
            sage: G._filling
            {(1, 2): 1, (2, 1): 1, (4, 0): 1}
            sage: G._in_labels
            [[], [], [], [], [], [], [1], [1], [1], [1, 1], [1, 1]]

            sage: labels = [[], [1],[1],[2],[2],[2,1],[2],[2,1],[1,1],[2,1],[1,1]]
            sage: G = GrowthDiagramRSK(shape=SkewPartition([[5,4,3,2,1],[3,2,1]]), labels=labels)
            sage: G._filling
            {(1, 2): 1, (2, 1): 1, (4, 0): 1}
            sage: G._in_labels
            [[], [], [], [1], [1], [1], [1], [1], [1], [1, 1], [1, 1]]

        """
        assert self._backward_rule is not None, "For computing a filling, a backward rule has to be specified."
        F = dict()
        labels = copy(self._out_labels)
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

        self._in_labels = labels
        self._filling = F

class GrowthDiagramBinWord(GrowthDiagram):
    """A class modelling a Schensted-like correspondence for binary
    words.
    """
    def __init__(self,
                 filling = None,
                 shape = None,
                 labels = None):
        # TODO: should check that the filling is standard
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
    """
    A class modelling Robinson-Schensted-Knuth insertion.
    """
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

class GrowthDiagramDomino(GrowthDiagram):
    """
    A class modelling domino insertion.

    EXAMPLES::

        sage: G = GrowthDiagramDomino([[1]]); G
        1
        sage: G._out_labels
        [[], [2], []]

        sage: G = GrowthDiagramDomino([[-1]]); G
        -1
        sage: G._out_labels
        [[], [1, 1], []]


        sage: G = GrowthDiagramDomino([[0,0,0,-1],[0,0,1,0],[-1,0,0,0],[0,1,0,0]]); G
         0  0  0 -1
         0  0  1  0
        -1  0  0  0
         0  1  0  0
        sage: G._out_labels
        [[], [1, 1], [3, 1], [3, 3], [3, 3, 2], [2, 2, 2], [2, 2], [1, 1], []]

    """
    def __init__(self,
                 filling = None,
                 shape = None,
                 labels = None):
        if labels is not None:
            labels = [Partition(la) for la in labels]
        super(GrowthDiagramDomino, self).__init__(filling = filling,
                                                  shape = shape,
                                                  labels = labels,
                                                  zero = Partition([]),
                                                  forward_rule = Domino_forward,
                                                  backward_rule = None)

    def _shape_from_labels(self, labels):
        def right_left(la, mu):
            if la.size() < mu.size():
                return 1
            elif la.size() > mu.size():
                return 0
            else:
                raise ValueError("Can only determine the shape of the growth diagram if sizes of successive partitions differ.")
        return Partitions().from_zero_one([right_left(labels[i], labels[i+1]) for i in range(len(labels)-1)])
