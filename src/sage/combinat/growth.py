r"""
Growth diagrams and dual graded graphs

AUTHORS:

- Martin Rubey (2016-09): Initial version

.. TODO::

    - when shape is given, check that it is compatible with filling or labels
    - implement backward rules for :class:`GrowthDiagramDomino` and :class:`GrowthDiagramSylvester`
    - optimise rules, mainly for :class:`GrowthDiagramRSK` and :class:`GrowthDiagramBurge`
    - make semistandard extension generic

Growth diagrams, invented by Sergey Fomin [Fom1995]_, provide a vast
generalisation of the Robinson-Schensted-Knuth correspondence between
matrices with non-negative integer entries and pairs of semistandard
Young tableaux of the same shape.

The main fact is that many correspondences similar to RSK can be
defined by providing a so-called 'forward' rule: a function whose
input are three vertices x, y and t of a certain directed graph (in
the case of Robinson-Schensted: the directed graph corresponding to
Young's lattice) and an integer (in the case of Robinson-Schensted:
zero or one), and whose output is a fourth vertex z.  This rule
should be invertible in the following sense: there is a so-called
'backward' rule that recovers the integer and t given z, x and y.

The classical Robinson-Schensted-Knuth correspondence is provided by
:class:`GrowthDiagramRSK`.  Note that a growth diagram is printed
with matrix coordinates, the origin being in the top-left corner::

    sage: w = [2,3,6,1,4,5]; G = GrowthDiagramRSK(w); G
    0  0  0  1  0  0
    1  0  0  0  0  0
    0  1  0  0  0  0
    0  0  0  0  1  0
    0  0  0  0  0  1
    0  0  1  0  0  0

The 'forward' rule just mentioned assigns partitions to the corners
of each of the 36 cells of this matrix - with the exception of the
corners on the left and top boundary, which are initialized with the
empty partition.  The partitions along the boundary opposite of the
origin are obtained by using the method
:meth:`GrowthDiagramRSK.out_labels()`::

    sage: G.out_labels()
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

However, in the case of a rectangular filling, it is more practical
to split this sequence of labels in two.  We then obtain the `P` and
`Q` symbols::

    sage: [G.P_symbol(), G.Q_symbol()]
    [[[1, 3, 4, 5], [2, 6]], [[1, 2, 3, 6], [4, 5]]]
    sage: RSK(w)
    [[[1, 3, 4, 5], [2, 6]], [[1, 2, 3, 6], [4, 5]]]

A great advantage of growth diagrams is that we immediately have
access also to the skew version of the RSK-correspondence, by
providing different initialisation for the labels near the origin.
We reproduce the original example of Bruce Sagan and Richard Stanley,
see also Tom Roby's thesis [Rob1991]_.  We can represent the
generalised permutation::

    1 2 4
    4 2 3

as a dictionary of coordinates, subtracting `1` from all entries
because lists in SageMath are zero-based::

    sage: w = {(1-1,4-1):1, (2-1,2-1):1, (4-1,3-1):1}
    sage: T = SkewTableau([[None, None],[None,5],[1]])
    sage: U = SkewTableau([[None, None],[None,3],[5]])
    sage: G = GrowthDiagramRSK(filling = w, shape = [5]*5, labels = T.to_chain()[::-1]+U.to_chain()[1:]); G
    0  0  0  0  0
    0  1  0  0  0
    0  0  0  1  0
    1  0  0  0  0
    0  0  0  0  0
    sage: G.P_symbol(), G.Q_symbol()
    ([[None, None, 2, 3], [None, None], [None, 4], [1], [5]],
     [[None, None, 1, 4], [None, None], [None, 2], [3], [5]])

Moreover, we are not forced to use rectangular fillings.  For
example, consider the Stanley-Sundaram correspondence between (skew)
oscillating tableaux and (partial) perfect matchings.  Again, from
Tom Roby's thesis::

    sage: o = [[2,1],[2,2],[3,2],[4,2],[4,1],[4,1,1],[3,1,1],[3,1],[3,2],[3,1],[2,1]]
    sage: l = [None]*(2*len(o)-1)
    sage: l[::2] = [Partition(la) for la in o]
    sage: l[1::2] = [l[2*i] if l[2*i].size() < l[2*i+2].size() else l[2*i+2] for i in range(len(o)-1)]
    sage: G = GrowthDiagramRSK(labels=l[1:-1], shape=[i for i in range(len(o)-2,0,-1)]); G
    0  0  0  0  0  0  0  1  0
    0  1  0  0  0  0  0  0
    0  0  0  0  0  0  0
    0  0  0  0  0  0
    0  0  1  0  0
    0  0  0  0
    0  0  0
    0  0
    0
    sage: ascii_art(SkewTableau(chain=G.in_labels()[len(o)-2:]), SkewTableau(chain=G.in_labels()[:len(o)-1][::-1]))
    .  1  .  7
    5     4

As mentioned at the beginning, the Robinson-Schensted-Knuth
correspondence is just a special case of growth diagrams.  In
particular, we have implemented local rules for the variation of RSK
originally due to Burge (:class:`GrowthDiagramBurge`), a
correspondence producing binary words originally due to Viennot
(:class:`GrowthDiagramBinWord`), and a correspondence producing
domino tableaux (:class:`GrowthDiagramDomino`) originally due to
Barbasch and Vogan.

REFERENCES:

.. [Fom1995] Sergey V. Fomin.
   *Schensted algorithms for dual graded graphs*.
   Journal of Algebraic Combinatorics Volume 4, Number 1 (1995), pp. 5-45

.. [Rob1991] Tom Roby.
   *Applications and extensions of Fomin's generalization of the Robinson-Schensted correspondence to differential posets*.
   M.I.T., Cambridge, Massachusetts

.. [Kra2006] Christian Krattenthaler.
   *Growth diagrams, and increasing and decreasing chains in fillings of Ferrers shapes*.
   Advances in Applied Mathematics Volume 37, Number 3 (2006), pp. 404-431

.. [Lam2004] Thomas Lam.
   *Growth diagrams, domino insertion and sign-imbalance*.
   Journal of Combinatorial Theory, Series A Volume 107, Number 1 (2004), pp. 87-115

.. [LamShi2007] Thomas Lam and Mark Shimozono.
   *Dual graded graphs for Kac-Moody algebras*.
   Algebra & Number Theory 1.4 (2007): pp. 451-488

.. [HivNze] Florent Hivert and Janvier Nzeutchap.
   *Dual Graded Graphs in Combinatorial Hopf Algebras*.
   https://www.lri.fr/~hivert/PAPER/commCombHopfAlg.pdf

.. [Vie1983] Xavier G. Viennot.
   *Maximal chains of subwords and up-down sequences of permutations*.
   Journal of Combinatorial Theory, Series A Volume 34, (1983), pp. 1-14

.. [Nze2007] Janvier Nzeutchap.
   *Binary Search Tree insertion, the Hypoplactic insertion, and Dual Graded Graphs*.
   :arXiv:`0705.2689` (2007)

"""
from sage.structure.sage_object import SageObject
from sage.combinat.words.word import Word
from sage.combinat.binary_tree import BinaryTree
from sage.combinat.partition import Partition, Partitions
from sage.combinat.skew_partition import SkewPartition
from sage.combinat.skew_tableau import SkewTableau
from copy import copy

class GrowthDiagram(SageObject):
    r"""
    The base class all variants of growth diagrams inherit from.

    Inheriting classes should provide an ``__init__`` method that
    checks that ``labels``, when provided, are of the correct
    type, and sets the following attributes:

    - ``self._zero``, the the zero element of the vertices of the
      graphs,

    - ``self._rank_function``, the rank function of the dual
      graded graphs,

    - ``self._covers_1``, ``self._covers_2``, functions taking
      two vertices as arguments and return True if the first
      covers the second in the respective graded graph.

    It should then call the ``__init__`` method of this class.

    EXAMPLES::

        sage: w = [3,3,2,4,1]; G = GrowthDiagramRSK(w)
        sage: [G.P_symbol(), G.Q_symbol()]
        [[[1, 3, 4], [2], [3]], [[1, 2, 4], [3], [5]]]
        sage: RSK(w)
        [[[1, 3, 4], [2], [3]], [[1, 2, 4], [3], [5]]]

    TESTS::

        sage: G = GrowthDiagramRSK()
        Traceback (most recent call last):
        ...
        ValueError: please provide a filling or a sequence of labels.
    """
    def __init__(self,
                 filling = None,
                 shape = None,
                 labels = None):
        r"""
        Initialise a generalized Schensted growth diagram in the sense of
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

        - ``filling`` is ``None``, if the growth diagram should be
          determined by applying the backward rules to the ``labels``
          decorating the boundary opposite of the origin of the
          ``shape``.  Otherwise it can be a dictionary with keys
          being coordinates and integer values, a sequence of
          sequences of integers (including matrices), or a word with
          integer letters (including permutations).

        - ``shape`` is a (possibly skew) partition or ``None``.  In
          the latter case it is determined as the Ferrers shape given
          by ``filling``, if that is a sequence of sequences, the
          bounding rectangle (including the origin) of ``filling``,
          or, if possible, inferred from ``labels``, if ``filling``
          is ``None``.

        - ``labels`` is ``None`` or a list.  If it is a list, its
          length should be the half-perimeter of ``shape``.  If
          ``filling`` is ``None``, its elements are the labels of the
          boundary opposite of the origin.  Otherwise its elements
          are the labels on the boundary on the side of the origin.
          If ``labels`` is ``None`` (in which case ``filling`` must
          not be ``None``) the value of ``self._zero`` is used to
          initialise ``labels``.
        """
        try:
            self._covers_1((self._zero, self._zero))
        except AttributeError:
            self._covers_1 = lambda (a,b): True
        try:
            self._covers_2((self._zero, self._zero))
        except AttributeError:
            self._covers_2 = lambda (a,b): True

        if filling is None:
            if labels is None:
                raise ValueError("please provide a filling or a sequence of labels.")

            if shape is None:
                shape = self._shape_from_labels(labels)

            self._lambda, self._mu = self._init_shape_from_various_input(shape)
            self._out_labels = labels
            self._check_labels(self._out_labels)
            self._shrink()
        else:
            self._filling, (self._lambda, self._mu) = self._init_filling_and_shape_from_various_input(filling, shape)
            self._in_labels = self._init_labels_forward_from_various_input(labels)
            self._check_labels(self._in_labels)
            self._grow()

    def filling(self):
        r"""
        Return the filling of the diagram as a dictionary.

        EXAMPLES::

            sage: G = GrowthDiagramRSK([[0,1,0], [1,0,2]])
            sage: G.filling()
            {(0, 1): 1, (1, 0): 1, (2, 1): 2}
        """
        return self._filling

    def conjugate(self):
        r"""
        Return the conjugate growth diagram of ``self``. This is
        the growth diagram with the filling reflected over the
        main diagonal.

        When the filling is a permutation, the conjugate filling
        corresponds to its inverse.

        EXAMPLES::

            sage: G = GrowthDiagramRSK([[0,1,0], [1,0,2]])
            sage: Gc = G.conjugate()
            sage: (Gc.P_symbol(), Gc.Q_symbol()) == (G.Q_symbol(), G.P_symbol())
            True
        """
        F = {(j,i): v for (i,j),v in self._filling.items()}
        return self.parent()(filling = F)

    def rotate(self):
        r"""
        Return the growth diagram with the filling rotated by 180 degrees.

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
        F = {(max_row-i,max_col-j): v for (i,j),v in self._filling.items()}
        return self.parent()(filling = F)

    def shape(self):
        r"""
        Return the shape of the growth diagram as a skew partition.

        .. WARNING::

            In the literature the label on the corner opposite of the
            origin of a rectangular filling is often called the shape
            of the filling.  This method returns the shape of the
            region instead.

        EXAMPLES::

            sage: GrowthDiagramRSK([1]).shape()
            [1] / []

        """
        return SkewPartition([self._lambda, self._mu])

    def out_labels(self):
        r"""
        Return the labels along the boundary opposite of the origin.

        EXAMPLES::

            sage: G = GrowthDiagramRSK([[0,1,0], [1,0,2]])
            sage: G.out_labels()
            [[], [1], [1, 1], [3, 1], [1], []]
        """
        return self._out_labels

    def in_labels(self):
        r"""
        Return the labels along the boundary on the side of the origin.

        EXAMPLES::

            sage: G = GrowthDiagramRSK(labels=[[2,2],[3,2],[3,3],[3,2]]); G
            1 0
            sage: G.in_labels()
            [[2, 2], [2, 2], [2, 2], [3, 2]]
        """
        return self._in_labels

    def P_chain(self):
        r"""
        Return the labels along the vertical boundary of a rectangular
        growth diagram.

        EXAMPLES::

            sage: G = GrowthDiagramBinWord([4, 1, 2, 3])
            sage: G.P_chain()
            [word: , word: 1, word: 11, word: 111, word: 1011]
        """
        if self.is_rectangular():
            return self._out_labels[self._lambda[0]:][::-1]
        else:
            raise ValueError("The P symbol is only defined for rectangular shapes.")

    def Q_chain(self):
        r"""
        Return the labels along the horizontal boundary of a rectangular
        growth diagram.

        EXAMPLES::

            sage: G = GrowthDiagramBinWord([[0,1,0,0], [0,0,1,0], [0,0,0,1], [1,0,0,0]])
            sage: G.Q_chain()
            [word: , word: 1, word: 10, word: 101, word: 1011]
        """
        if self.is_rectangular():
            return self._out_labels[:self._lambda[0]+1]
        else:
            raise ValueError("The Q symbol is only defined for rectangular shapes.")

    def is_rectangular(self):
        r"""
        Return ``True`` if the shape of the growth diagram is rectangular.

        EXAMPLES::

            sage: GrowthDiagramRSK([2,3,1]).is_rectangular()
            True
            sage: GrowthDiagramRSK([[1,0,1],[0,1]]).is_rectangular()
            False
        """
        return all(x == 0 for x in self._mu) and all(x == self._lambda[0] for x in self._lambda)

    def to_word(self):
        r"""
        Return the filling as a word, if the shape is rectangular and
        there is at most one nonzero entry in each column, which must
        be 1.

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
            for ((i,j), v) in self._filling.items():
                if v != 0:
                    if v == 1:
                        if w[i] == 0:
                            w[i] = j+1
                        else:
                            raise ValueError("Can only convert fillings with at most one entry per column to words.")
                    else:
                        raise ValueError("Can only convert 0-1 fillings to words.  Try ``to_biword``.")
            return w
        else:
            raise ValueError("Can only convert fillings of rectangular shapes to words.")

    def to_biword(self):
        r"""
        Return the filling as a biword, if the shape is rectangular.

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
            for ((i,j), v) in sorted(self._filling.items()):
                if v >= 0:
                    w1.extend([i+1]*v)
                    w2.extend([j+1]*v)
                else:
                    raise ValueError("Can only convert fillings with nonnegative entries to words.")
            return (w1, w2)
        else:
            raise ValueError("Can only convert fillings of rectangular shapes to words.")

    def __iter__(self):
        r"""
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
        r"""
        Print the filling of the growth diagram as a skew tableau.

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

    def __eq__(self, other):
        r"""
        Return ``True`` if the growth diagram ``other`` has the same
        shape and the same filling as ``self``.

        EXAMPLES:

        Equality ignores zeros in fillings::

            sage: G1 = GrowthDiagramRSK({(0, 1): 1, (1, 0): 1})
            sage: G2 = GrowthDiagramRSK({(0, 0): 0, (0, 1): 1, (1, 0): 1})
            sage: G1 == G2
            True

        Growth diagrams with different shapes are different::

            sage: G1 = GrowthDiagramRSK([[0,1,0],[1,0]])
            sage: G2 = GrowthDiagramRSK([[0,1,0],[1]])
            sage: G1 == G2
            False

        Growth diagrams with different rules are different::

            sage: G1 = GrowthDiagramRSK({(0, 1): 1, (1, 0): 1})
            sage: G2 = GrowthDiagramBinWord({(0, 1): 1, (1, 0): 1})
            sage: G1 == G2
            False
        """
        return (self.parent() == other.parent() and
                self._lambda == other._lambda and
                self._mu == other._mu and
                self._filling == other._filling)

    def _half_perimeter(self):
        r"""
        Return half the perimeter of the shape of the growth diagram.

        Assumes that ``self._lambda`` is already set.

        TESTS::

            sage: G = GrowthDiagramRSK({(0,1):1, (2,0):1}, SkewPartition([[3,1],[1]])); G
            .  0  1
            1
            sage: G._half_perimeter()
            6

        """
        if len(self._lambda) == 0:
            return 1
        else:
            return self._lambda[0]+len(self._lambda)+1

    def _shape_from_labels(self, labels):
        r"""
        Determine the shape of the growth diagram given a list of labels
        during initialisation.

        The shape can be determined from the labels if the size of
        each label differs from the size of its successor.

        Assumes that ``self._rank_function`` is set.

        Otherwise raise an error.

        TESTS::

            sage: labels = [[],[1],[],[1],[]]
            sage: G = GrowthDiagramRSK(labels=labels); G
            0 1
            1
            sage: G._shape_from_labels(G.out_labels())
            [2, 1]

        """
        def right_left(la, mu):
            if self._rank_function(la) < self._rank_function(mu):
                assert self._covers_1((mu, la)), "%s has smaller rank than %s but isn't covered by it!" %(la, mu)
                return 1
            elif self._rank_function(la) > self._rank_function(mu):
                assert self._covers_2((la, mu)), "%s has smaller rank than %s but isn't covered by it!" %(mu, la)
                return 0
            else:
                raise ValueError("Can only determine the shape of the growth diagram if ranks of successive labels differ.")
        return Partitions().from_zero_one([right_left(labels[i], labels[i+1]) for i in range(len(labels)-1)])

    def _check_labels(self, labels):
        r"""
        Check sanity of the parameter ``labels``.

        Assumes that ``self._lambda`` is already set.

        TEST::

            sage: GrowthDiagramRSK(shape=[1], labels=[[], [1]])                 # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: the number of labels is 2, but for this shape we need 3.

            sage: GrowthDiagramRSK(labels=[[], [1], [2], [2,1]])                # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: the number of labels is 4, but for this shape we need 1.

        .. TODO::

            Can we do something more sensible when the chain of labels is strictly increasing?
        """
        half_perimeter = self._half_perimeter()
        if len(labels) != half_perimeter:
            raise ValueError("the number of labels is %s, but for this shape we need %s."
                             %(len(labels), half_perimeter))

    def _init_labels_forward_from_various_input(self, labels):
        r"""
        Return a list of labels decorating the boundary near the origin.

        Assumes that ``self._lambda`` is already set.

        TESTS:

        ``labels`` is None::

            sage: filling = []
            sage: shape = SkewPartition([[4,2,1,1],[2,1,1]])
            sage: G = GrowthDiagramRSK(filling, shape)
            sage: G.in_labels()                                                 # indirect doctest
            [[], [], [], [], [], [], [], [], []]

        ``labels`` is a list of partitions::

            sage: filling = []
            sage: labels = [[],[1],[],[1],[]]
            sage: shape = SkewPartition([[2,1],[1]])
            sage: G = GrowthDiagramRSK(filling=filling, shape=shape, labels=labels)
            sage: G.in_labels()                                                 # indirect doctest
            [[], [1], [], [1], []]

        """
        if labels is None:
            return [self._zero for i in range(self._half_perimeter())]
        else:
            return labels

    def _init_shape_from_various_input(self, shape):
        r"""
        Return a pair of partitions describing the region of the growth
        diagram.

        Assumes that ``self._filling`` is already set.

        TESTS:

        ``shape`` is a skew partition::

            sage: filling = []
            sage: shape = SkewPartition([[4,2,1,1],[2,1,1]])
            sage: G = GrowthDiagramRSK(filling, shape)
            sage: G._lambda, G._mu                                              # indirect doctest
            ([4, 2, 1, 1], [2, 1, 1, 0])

        ``shape`` is a partition::

            sage: filling = []
            sage: shape = Partition([3,2,1,1])
            sage: G = GrowthDiagramRSK(filling, shape)
            sage: G._lambda, G._mu                                              # indirect doctest
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
        r"""
        Return a dict ``F``, such that ``F[(i,j)]`` is the element in row ``i``
        and column ``j``, and a pair of partitions describing the
        region of the growth diagram.

        TESTS:

        ``filling`` is a dict of coordinates::

            sage: pi = Permutation([2,3,1,6,4,5])
            sage: G = GrowthDiagramRSK({(i,pi[i]-1):1 for i in range(len(pi))})
            sage: G._filling                                                    # indirect doctest
            {(0, 1): 1, (1, 2): 1, (2, 0): 1, (3, 5): 1, (4, 3): 1, (5, 4): 1}
            sage: G.shape()
            [6, 6, 6, 6, 6, 6] / []

        ``filling`` is a dict of dicts::

            sage: G = GrowthDiagramRSK({i:{pi[i]-1:1} for i in range(len(pi))})
            sage: G._filling                                                    # indirect doctest
            {(0, 1): 1, (1, 2): 1, (2, 0): 1, (3, 5): 1, (4, 3): 1, (5, 4): 1}
            sage: G.shape()
            [6, 6, 6, 6, 6, 6] / []

        ``filling`` is a matrix::

            sage: G = GrowthDiagramRSK(pi.to_matrix())
            sage: G._filling                                                    # indirect doctest
            {(0, 1): 1, (1, 2): 1, (2, 0): 1, (3, 5): 1, (4, 3): 1, (5, 4): 1}
            sage: G.shape()
            [6, 6, 6, 6, 6, 6] / []

        ``filling`` is a permutation::

            sage: G = GrowthDiagramRSK(pi)
            sage: G._filling                                                    # indirect doctest
            {(0, 1): 1, (1, 2): 1, (2, 0): 1, (3, 5): 1, (4, 3): 1, (5, 4): 1}
            sage: G.shape()
            [6, 6, 6, 6, 6, 6] / []

        ``filling`` is a list::

            sage: G = GrowthDiagramRSK([3,1,4,1,5])
            sage: G._filling                                                    # indirect doctest
            {(0, 2): 1, (1, 0): 1, (2, 3): 1, (3, 0): 1, (4, 4): 1}
            sage: G.shape()
            [5, 5, 5, 5, 5] / []

        ``filling`` is a list of lists::

            sage: G = GrowthDiagramRSK([[1,0,1],[0,1]])
            sage: G._filling                                                    # indirect doctest
            {(0, 0): 1, (1, 1): 1, (2, 0): 1}
            sage: G.shape()
            [3, 2] / []

        ``filling`` is a list of lists and shape is given::

            sage: G = GrowthDiagramRSK([[1,0,1],[0,1]], shape=SkewPartition([[3,2],[1]]))
            sage: G._filling                                                    # indirect doctest
            {(0, 0): 1, (1, 1): 1, (2, 0): 1}
            sage: G.shape()
            [3, 2] / [1]

        ``filling`` is empty and shape is ``None``::

            sage: G = GrowthDiagramRSK({})
            sage: (G.filling(), G.shape())
            ({}, [] / [])

        """
        if isinstance(filling, dict):
            try:
                from six import itervalues
                v = next(itervalues(filling))
                if isinstance(v, dict):
                    # it is a dict of dicts
                    F = dict()
                    for (i, row) in filling.items():
                        for (j, v) in row.items():
                            if v != 0:
                                F[(i,j)] = int(v)
                else:
                    # it is dict of coordinates
                    F = {(i,j): v for ((i,j), v) in filling.items()
                         if v != 0}
            except StopIteration:
                # it is an empty dict of coordinates
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
            if F == {}:
                shape = []
            else:
                # find bounding rectangle of ``filling``
                max_row = max(i for i, _ in F)+1
                max_col = max(j for _, j in F)+1
                shape = [max_row]*max_col

        return (F, self._init_shape_from_various_input(shape))

    def _grow(self):
        r"""
        Compute the labels on the boundary opposite of the origin, given
        the filling.

        TESTS::

            sage: pi = Permutation([1])
            sage: G = GrowthDiagramRSK(pi)
            sage: G._out_labels                                                 # indirect doctest
            [[], [1], []]

            sage: pi = Permutation([1,2])
            sage: G = GrowthDiagramRSK(pi)
            sage: G._out_labels                                                 # indirect doctest
            [[], [1], [2], [1], []]

            sage: pi = Permutation([2,1])
            sage: G = GrowthDiagramRSK(pi)
            sage: G._out_labels                                                 # indirect doctest
            [[], [1], [1, 1], [1], []]

            sage: G = GrowthDiagramRSK({(0,1):1, (1,0):1}, SkewPartition([[2,1],[1]]))
            sage: G._out_labels                                                 # indirect doctest
            [[], [1], [], [1], []]

            sage: G = GrowthDiagramRSK({(1,1):1}, SkewPartition([[2,2],[1]]), labels=[[],[],[1],[],[]])
            sage: G._out_labels                                                 # indirect doctest
            [[], [1], [2], [1], []]

            sage: G = GrowthDiagramBinWord({(1,1):1}, SkewPartition([[2,2],[1]]), labels=[[],[],[1],[],[]])
            sage: G._out_labels                                                 # indirect doctest
            [word: , word: 1, word: 11, word: 1, word: ]
        """
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
        r"""
        Compute the labels on the boundary near the origin, and the filling.

        TESTS::

            sage: filling = [[0,0,1,0,0,0,0], [0,1,0,0,0,0,0], [1,0,0,0,0,0,0], [0,0,0,1,0,0,0], [0,0,0,0,0,0,1], [0,0,0,0,0,1,0], [0,0,0,0,1,0,0]]
            sage: G = GrowthDiagramRSK(filling)
            sage: list(GrowthDiagramRSK(labels=G._out_labels)) == filling
            True

            sage: labels = [[], [1], []]
            sage: G = GrowthDiagramRSK(labels=labels)
            sage: G._filling
            {(0, 0): 1}
            sage: G._in_labels                                                  # indirect doctest
            [[], [], []]

            sage: labels = [[], [1], [2], [2,1], [1,1], [1], []]
            sage: G = GrowthDiagramRSK(labels=labels)
            sage: G._filling
            {(0, 1): 1, (1, 2): 1, (2, 0): 1}
            sage: G._in_labels                                                  # indirect doctest
            [[], [], [], [], [], [], []]

            sage: labels = [[], [1], [2], [3], [3, 1], [3, 2], [4, 2], [4, 1], [3, 1], [2, 1], [1, 1], [1], []]
            sage: G = GrowthDiagramRSK(labels=labels)
            sage: G._filling                                                    # indirect doctest
            {(0, 1): 1, (1, 2): 1, (2, 5): 1, (3, 0): 1, (4, 3): 1, (5, 4): 1}

            sage: labels = [[],[1],[1],[2],[2],[2,1],[2]]
            sage: G = GrowthDiagramRSK(labels=labels)
            Traceback (most recent call last):
            ...
            ValueError: Can only determine the shape of the growth diagram if ranks of successive labels differ.
            sage: G = GrowthDiagramRSK(shape=[3,2,1], labels=labels)
            sage: G._filling                                                    # indirect doctest
            {(1, 0): 1}
            sage: G._in_labels
            [[], [], [], [], [1], [1], [2]]

            sage: labels = [[], [1],[1],[2],[2],[2,1],[2],[2,1],[1,1],[2,1],[1,1]]
            sage: G = GrowthDiagramRSK(shape=[5,4,3,2,1], labels=labels)
            sage: G._filling
            {(1, 2): 1, (2, 1): 1, (4, 0): 1}
            sage: G._in_labels                                                  # indirect doctest
            [[], [], [], [], [], [], [1], [1], [1], [1, 1], [1, 1]]

            sage: labels = [[], [1],[1],[2],[2],[2,1],[2],[2,1],[1,1],[2,1],[1,1]]
            sage: G = GrowthDiagramRSK(shape=SkewPartition([[5,4,3,2,1],[3,2,1]]), labels=labels)
            sage: G._filling
            {(1, 2): 1, (2, 1): 1, (4, 0): 1}
            sage: G._in_labels                                                  # indirect doctest
            [[], [], [], [1], [1], [1], [1], [1], [1], [1, 1], [1, 1]]

        """
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
    r"""
    A class modelling a Schensted-like correspondence for binary
    words.

    EXAMPLES::

        sage: pi = Permutation([4,1,8,3,6,5,2,7,9]); G = GrowthDiagramBinWord(pi); G
        0  1  0  0  0  0  0  0  0
        0  0  0  0  0  0  1  0  0
        0  0  0  1  0  0  0  0  0
        1  0  0  0  0  0  0  0  0
        0  0  0  0  0  1  0  0  0
        0  0  0  0  1  0  0  0  0
        0  0  0  0  0  0  0  1  0
        0  0  1  0  0  0  0  0  0
        0  0  0  0  0  0  0  0  1
        sage: G.out_labels()[9]
        word: 101010011

    The Kleitman Greene invariant is the descent word::

        sage: pi.descents(from_zero=False)
        [1, 3, 5, 6]

    .. automethod:: _forward_rule
    .. automethod:: _backward_rule

    TESTS::

        sage: G = GrowthDiagramBinWord([1,3,2])
        sage: G._zero
        word:
        sage: G = GrowthDiagramBinWord(labels = [[1,1],[1,1,0],[0,1]])
        Traceback (most recent call last):
        ...
        AssertionError: 01 has smaller rank than 110 but isn't covered by it!

        sage: G = GrowthDiagramBinWord(labels = [[1,1],[1,0,1],[0,1]])
        Traceback (most recent call last):
        ...
        AssertionError: 11 has smaller rank than 101 but isn't covered by it!

        sage: pi = Permutations(10).random_element()
        sage: G = GrowthDiagramBinWord(pi)
        sage: list(GrowthDiagramBinWord(labels=G.out_labels())) == list(G)
        True

    Test that the Kleitman Greene invariant is indeed the descent word::

        sage: r=4; all(Word([0 if i in w.descents(from_zero=False) else 1 for i in range(r)]) == GrowthDiagramBinWord(w).out_labels()[r] for w in Permutations(r))
        True

    """
    def __init__(self,
                 filling = None,
                 shape = None,
                 labels = None):
        # TODO: should check that the filling is standard
        if labels is not None:
            labels = [Word(la, alphabet=[0,1]) for la in labels]
        self._zero = Word([], alphabet=[0,1])
        self._rank_function = lambda w: len(w)
        self._covers_1 = lambda (w, v): w[:-1] == v
        self._covers_2 = lambda (w, v): v.is_subword_of(w)
        super(GrowthDiagramBinWord, self).__init__(filling = filling,
                                                   shape = shape,
                                                   labels = labels)
    __init__.__doc__ = GrowthDiagram.__init__.__doc__


    @staticmethod
    def _forward_rule(y, t, x, content):
        r"""
        Return the output shape given three shapes and the content.

        See [Fom1995]_ Lemma 4.6.1, page 40.

        INPUT:

        - ``y, t, x`` -- three binary words from a cell in a growth
          diagram, labelled as::

              t x
              y

        - ``content`` -- 0 or 1, the content of the cell.

        OUTPUT:

        The fourth binary word z according to Viennot's bijection [Vie1983]_.

        TESTS::

            sage: G = GrowthDiagramBinWord

            sage: G._forward_rule([], [], [], 1)
            word: 1

            sage: G._forward_rule([1], [1], [1], 1)
            word: 11

        if ``x != y`` append last letter of ``x`` to ``y``::

            sage: G._forward_rule([1,0], [1], [1,1], 0)
            word: 101

        if ``x == y != t`` append ``0`` to ``y``::

            sage: G._forward_rule([1,1], [1], [1,1], 0)
            word: 110

        """
        if x == t == y:
            if content == 0:
                z = x
            elif content == 1:
                z = Word(list(y) + [1], alphabet=[0,1])
            else:
                raise NotImplementedError
        elif content != 0:
            raise ValueError("For y=%s, t=%s, x=%s, the content should be 0 but is %s" %(y, t, x, content))
        elif x != t == y:
            z = x
        elif x == t != y:
            z = y
        else:
            if x != y:
                z = Word(list(y) + [x[-1]], alphabet=[0,1])
            elif x == y != t:
                z = Word(list(y) + [0], alphabet=[0,1])
            else:
                raise NotImplementedError
        return z

    @staticmethod
    def _backward_rule(y, z, x):
        r"""
        Return the content and the input shape.

        See [Fom1995]_ Lemma 4.6.1, page 40.

        - ``y, z, x`` -- three binary words from a cell in a growth diagram,
          labelled as::

                x
              y z

        OUTPUT:

        A pair ``(t, content)`` consisting of the shape of the fourth
        word and the content of the cell acording to Viennot's
        bijection [Vie1983]_.

        TEST::

            sage: w = [4,1,8,3,6,5,2,7,9]; G = GrowthDiagramBinWord(w);
            sage: GrowthDiagramBinWord(labels=G.out_labels()).to_word() == w    # indirect doctest
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

class GrowthDiagramSylvester(GrowthDiagram):
    r"""
    A class modelling a Schensted-like correspondence for binary
    trees.

    EXAMPLES:

    From [Nze2007]_::

        sage: pi = Permutation([3,5,1,4,2,6]); G = GrowthDiagramSylvester(pi); G
        0  0  1  0  0  0
        0  0  0  0  1  0
        1  0  0  0  0  0
        0  0  0  1  0  0
        0  1  0  0  0  0
        0  0  0  0  0  1
        sage: ascii_art(G.out_labels()[len(pi)])
          __o__
         /     \
        o       o
         \     / \
          o   o   o

    .. automethod:: _forward_rule

    TESTS::

        sage: G = GrowthDiagramSylvester([1,3,2])
        sage: G._zero
        .

    """
    def __init__(self,
                 filling = None,
                 shape = None,
                 labels = None):
        # TODO: should check that the filling is standard
        if labels is not None:
            labels = [BinaryTree(la) for la in labels]
        self._zero = BinaryTree()
        self._rank_function = lambda w: w.node_number()
        super(GrowthDiagramSylvester, self).__init__(filling = filling,
                                                     shape = shape,
                                                     labels = labels)
    __init__.__doc__ = GrowthDiagram.__init__.__doc__

    @staticmethod
    def _forward_rule(y, t, x, content):
        r"""
        Return the output shape given three shapes and the content.

        See [Nze2007]_, page 9.

        INPUT:

        - ``y, t, x`` -- three binary trees from a cell in a growth
          diagram, labelled as::

              t x
              y

        - ``content`` -- 0 or 1, the content of the cell.

        OUTPUT:

        The fourth binary tree z.

        TESTS::

            sage: G = GrowthDiagramSylvester
            sage: B = BinaryTree; E = B(); N = B([]); L = B([[],None]); R = B([None,[]]); T = B([[],[]])

            sage: ascii_art(G._forward_rule(E, E, E, 1))
            o
            sage: ascii_art(G._forward_rule(N, N, N, 1))
            o
             \
              o
            sage: ascii_art(G._forward_rule(L, L, L, 1))
              o
             / \
            o   o
            sage: ascii_art(G._forward_rule(R, R, R, 1))
            o
             \
              o
               \
                o

        if ``x != y``, obtain ``z`` from ``x`` adding a node such
        that deleting the right most gives ``y``::

            sage: ascii_art(G._forward_rule(R, N, L, 0))
              o
             /
            o
             \
              o

            sage: ascii_art(G._forward_rule(L, N, R, 0))
              o
             / \
            o   o

        if ``x == y != t``, obtain ``z`` from ``y`` by adding a node
        as left child to the right most node::

            sage: ascii_art(G._forward_rule(N, E, N, 0))
              o
             /
            o
            sage: ascii_art(G._forward_rule(T, L, T, 0))
              _o_
             /   \
            o     o
                 /
                o

        """
        def successors(b):
            """
            Return all trees obtained from ``b`` by adding a node.
            """
            if b.is_empty():
                yield BinaryTree([])
            else:
                for t in successors(b[0]):
                    yield BinaryTree([t, b[1]])
                for t in successors(b[1]):
                    yield BinaryTree([b[0], t])

        def delete_right_most_node(b):
            """
            Return the tree obtained by deleting the right most node from ``b``.
            """
            if b.is_empty():
                raise ValueError("Cannot delete right most node from empty tree")
            elif b[1].is_empty():
                return b[0]
            else:
                return BinaryTree([b[0], delete_right_most_node(b[1])])

        def union(x, y):
            """
            Return the unique tree obtained by adding a node to ``x`` such
            that deleting the right most node gives ``y``.
            """
            for t in successors(x):
                if delete_right_most_node(t) == y:
                    return t
            raise ValueError("Couldn't find union of %s and %s" %(x,y))

        def add_left_child_to_right_most_node(b):
            """
            Return the tree obtained from ``b`` by adding a node as left
            child to the right most node.
            """
            if b.is_empty():
                raise ValueError("Cannot add left child to empty tree")
            elif b == BinaryTree([]):
                return BinaryTree([[], None])
            else:
                return BinaryTree([b[0], add_left_child_to_right_most_node(b[1])])

        if x == t == y:
            if content == 0:
                z = x
            elif content == 1:
                z = t.over(BinaryTree([]))
            else:
                raise NotImplementedError
        elif content != 0:
            raise ValueError("For y=%s, t=%s, x=%s, the content should be 0 but is %s" %(y, t, x, content))
        elif x != t == y:
            z = x
        elif x == t != y:
            z = y
        else:
            if x != y:
                return union(x, y)
            elif x == y != t:
                z = add_left_child_to_right_most_node(y)
            else:
                raise NotImplementedError
        return z

    @staticmethod
    def _backward_rule(y, z, x):
        r"""
        Return the content and the input shape.

        See [Fom1995]_ Lemma 4.6.1, page 40.

        - ``y, z, x`` -- three binary words from a cell in a growth diagram,
          labelled as::

                x
              y z

        OUTPUT:

        A pair ``(t, content)`` consisting of the shape of the fourth
        word and the content of the cell acording to Viennot's
        bijection [Vie1983]_.

        TEST::

            sage: w = [4,1,8,3,6,5,2,7,9]; G = GrowthDiagramBinWord(w);
            sage: GrowthDiagramBinWord(labels=G.out_labels()).to_word() == w    # indirect doctest
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

class GrowthDiagramYoungFibonacci(GrowthDiagram):
    r"""

    A class modelling a Schensted-like correspondence for
    Young-Fibonacci-tableaux.

    EXAMPLES::

        sage: G = GrowthDiagramYoungFibonacci([4,1,8,3,6,5,2,7,9]); G
        0  1  0  0  0  0  0  0  0
        0  0  0  0  0  0  1  0  0
        0  0  0  1  0  0  0  0  0
        1  0  0  0  0  0  0  0  0
        0  0  0  0  0  1  0  0  0
        0  0  0  0  1  0  0  0  0
        0  0  0  0  0  0  0  1  0
        0  0  1  0  0  0  0  0  0
        0  0  0  0  0  0  0  0  1
        sage: G.out_labels()[9]
        word: 122121

    The Kleitman Greene invariant is: take the last letter and the
    largest letter of the permutation and remove them.  If they
    coincide write 1, otherwise write 2.

    .. automethod:: _forward_rule
    .. automethod:: _backward_rule

    TESTS::

        sage: G = GrowthDiagramYoungFibonacci([1,3,2])
        sage: G._zero
        word:
        sage: G = GrowthDiagramYoungFibonacci(labels = [[1],[1,0],[1]])
        Traceback (most recent call last):
        ...
        ValueError: 0 not in alphabet!

        sage: G = GrowthDiagramYoungFibonacci(labels = [[1,1],[1,2]])
        Traceback (most recent call last):
        ...
        AssertionError: 11 has smaller rank than 12 but isn't covered by it!

        sage: G = GrowthDiagramYoungFibonacci(labels = [[1,2],[1,1]])
        Traceback (most recent call last):
        ...
        AssertionError: 11 has smaller rank than 12 but isn't covered by it!

        sage: pi = Permutations(10).random_element()
        sage: G = GrowthDiagramYoungFibonacci(pi)
        sage: list(GrowthDiagramYoungFibonacci(labels=G.out_labels())) == list(G)
        True

    """
    def __init__(self,
                 filling = None,
                 shape = None,
                 labels = None):
        # TODO: should check that the filling is standard
        if labels is not None:
            labels = [Word(la, alphabet=[1,2]) for la in labels]
        self._zero = Word([], alphabet=[1,2])
        self._rank_function = lambda w: sum(w)

        def covers(c):
            for i in range(len(c)+1):
                d = list(c)
                d.insert(i, 1)
                yield Word(d, alphabet=[1,2])
                if i < len(c) and c[i] == 1:
                    d = list(c)
                    d[i] = 2
                    yield Word(d, alphabet=[1,2])
                    break

        self._covers_1 = self._covers_2 = lambda (w, v): w in covers(v)

        super(GrowthDiagramYoungFibonacci, self).__init__(filling = filling,
                                                          shape = shape,
                                                          labels = labels)
    __init__.__doc__ = GrowthDiagram.__init__.__doc__

    @staticmethod
    def _forward_rule(shape3, shape2, shape1, content):
        r"""
        Return the output shape given three shapes and the content.

        See [Fom1995]_ Lemma 4.4.1, page 35.

        INPUT:

        - ``shape3, shape2, shape1`` -- three Fibonacci words from a
          cell in a growth diagram, labelled as::

              shape2 shape1
              shape3

        - ``content`` -- 0 or 1, the content of the cell.

        OUTPUT:

        The fourth Fibonacci word.

        TESTS::

            sage: G = GrowthDiagramYoungFibonacci

            sage: G._forward_rule([], [], [], 1)
            word: 1

            sage: G._forward_rule([1], [1], [1], 1)
            word: 11

            sage: G._forward_rule([1,2], [1], [1,1], 0)
            word: 21

            sage: G._forward_rule([1,1], [1], [1,1], 0)
            word: 21

        """
        if shape1 == shape2 == shape3:
            if content == 0:
                r = shape1
            elif content == 1:
                r = Word([1]+list(shape3), alphabet=[1,2])
            else:
                raise NotImplementedError
        elif content != 0:
            raise ValueError("For shape3=%s, shape2=%s, shape1=%s, the content should be 0 but is %s" %(shape3, shape2, shape1, content))
        elif shape1 == shape2:
            r = shape3
        elif shape3 == shape2:
            r = shape1
        else:
            if shape1 != shape2 != shape3:
                r = Word([2] + list(shape2), alphabet=[1,2])
            else:
                raise NotImplementedError("For shape3=%s, shape2=%s, shape1=%s, content %s we have no rule." %(shape3, shape2, shape1, content))
        return r

    @staticmethod
    def _backward_rule(y, z, x):
        r"""
        Return the content and the input shape.

        See [Fom1995]_ Lemma 4.4.1, page 35.

        - ``y, z, x`` -- three Fibonacci words from a cell in a
          growth diagram, labelled as::

                x
              y z

        OUTPUT:

        A pair ``(t, content)`` consisting of the shape of the fourth
        word and the content of the cell.

        TEST::

            sage: w = [4,1,8,3,6,5,2,7,9]; G = GrowthDiagramYoungFibonacci(w);
            sage: GrowthDiagramYoungFibonacci(labels=G.out_labels()).to_word() == w    # indirect doctest
            True

        """
        if x == y == z:
            return (x, 0)
        elif x == z != y:
            return (y, 0)
        elif x != z == y:
            return (x, 0)
        else:
            if z[0] == 1:
                return (z[1:], 1)
            elif z[0] == 2:
                return (z[1:], 0)


class GrowthDiagramOnPartitions(GrowthDiagram):
    r"""
    A class for growth diagrams on Young's lattice on integer
    partitions graded by size.  Since we do not use the covering
    relation, but only check partition containment, this class can
    also be used as a base class for growth diagrams modelling domino
    insertion, :class:`GrowthDiagramDomino`.

    TESTS::

        sage: G = GrowthDiagramBurge([])                                    # indirect doctest
        sage: G._zero
        []
        sage: G = GrowthDiagramBurge(labels = [[1],[1]])                    # indirect doctest
        Traceback (most recent call last):
        ...
        ValueError: Can only determine the shape of the growth diagram if ranks of successive labels differ.

    """
    def __init__(self,
                 filling = None,
                 shape = None,
                 labels = None):
        if labels is not None:
            labels = [Partition(la) for la in labels]
        self._zero = Partition([])
        self._rank_function = lambda p: p.size()
        super(GrowthDiagramOnPartitions, self).__init__(filling = filling,
                                                        shape = shape,
                                                        labels = labels)
    __init__.__doc__ = GrowthDiagram.__init__.__doc__

    def P_symbol(self):
        r"""
        Return the labels along the vertical boundary of a rectangular
        growth diagram as a skew tableau.

        EXAMPLES::

            sage: G = GrowthDiagramRSK([[0,1,0], [1,0,2]])
            sage: G.P_symbol().pp()
            1  2  2
            2
        """
        return SkewTableau(chain = self.P_chain())

    def Q_symbol(self):
        r"""
        Return the labels along the horizontal boundary of a rectangular
        growth diagram as a skew tableau.

        EXAMPLES::

            sage: G = GrowthDiagramRSK([[0,1,0], [1,0,2]])
            sage: G.Q_symbol().pp()
            1  3  3
            2
        """
        return SkewTableau(chain = self.Q_chain())

class GrowthDiagramRSK(GrowthDiagramOnPartitions):
    r"""
    A class modelling Robinson-Schensted-Knuth insertion.

    EXAMPLES::

        sage: pi = Permutation([2,3,6,1,4,5])
        sage: G = GrowthDiagramRSK(pi)
        sage: G.P_symbol(), G.Q_symbol()
        ([[1, 3, 4, 5], [2, 6]], [[1, 2, 3, 6], [4, 5]])
        sage: RSK(pi)
        [[[1, 3, 4, 5], [2, 6]], [[1, 2, 3, 6], [4, 5]]]

    .. automethod:: _forward_rule
    .. automethod:: _backward_rule

    """
    @staticmethod
    def _forward_rule(shape3, shape2, shape1, content):
        r"""
        Return the output shape given three shapes and the content.

        See [Kra2006]_ `(F^1 0)-(F^1 2)`.

        INPUT:

        - ``shape3, shape2, shape1`` -- three partitions from a cell in a
          growth diagram, labelled as::

              shape2 shape1
              shape3

        - ``content`` -- a non-negative integer, the content of the cell.

        OUTPUT:

        The fourth partition according to the Robinson-Schensted-Knuth
        correspondence.


        TESTS::

            sage: G = GrowthDiagramRSK
            sage: G._forward_rule([2,1],[2,1],[2,1],1)
            [3, 1]

            sage: G._forward_rule([1],[],[2],2)
            [4, 1]
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

    @staticmethod
    def _backward_rule(shape3, shape4, shape1):
        r"""
        Return the content and the input shape.

        See [Kra2006]_ `(B^1 0)-(B^1 2)`.

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

            sage: w = [4,1,8,3,6,5,2,7,9]; G = GrowthDiagramRSK(w)
            sage: GrowthDiagramRSK(labels=G._out_labels).to_word() == w         # indirect doctest
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

class GrowthDiagramBurge(GrowthDiagramOnPartitions):
    r"""
    A class modelling Burge insertion.

    EXAMPLES::

        sage: GrowthDiagramBurge(labels=[[],[1,1,1],[2,1,1,1],[2,1,1],[2,1],[1,1],[]])
        1  1
        0  1
        1  0
        1  0

    .. automethod:: _forward_rule
    .. automethod:: _backward_rule
    """
    @staticmethod
    def _forward_rule(shape3, shape2, shape1, content):
        r"""
        Return the output shape given three shapes and the content.

        See [Kra2006]_ `(F^4 0)-(F^4 2)`.

        INPUT:

        - ``shape3, shape2, shape1`` -- three  from a cell in a growth diagram,
          labelled as::

              shape2 shape1
              shape3

        - ``content`` -- a non-negative integer, the content of the cell.

        OUTPUT:

        The fourth partition according to the Burge correspondence.

        TESTS::

            sage: G = GrowthDiagramBurge
            sage: G._forward_rule([2,1],[2,1],[2,1],1)
            [3, 1]

            sage: G._forward_rule([1],[],[2],2)
            [2, 1, 1, 1]

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
            newPart = max(row1, row3) + min(int(row1 == row2 == row3), carry)
            if newPart == 0:
                return Partition(shape4[::-1])
            else:
                shape4 = [newPart] + shape4
                carry = carry - min(int(row1 == row2 == row3), carry) + min(row1, row3) - row2
                shape1 = shape1[1:]
                shape2 = shape2[1:]
                shape3 = shape3[1:]

    @staticmethod
    def _backward_rule(shape3, shape4, shape1):
        r"""
        Return the content and the input shape.

        See [Kra2006]_ `(B^4 0)-(B^4 2)`.  There is a typo in the
        computation of carry in `(B^4 2)` in the arXiv version of the
        article, `\rho` must be replaced by `\lambda`.

        INPUT:

        - ``shape3, shape4, shape1`` -- three partitions from a cell in a
          growth diagram, labelled as::

                   shape1         mu
            shape3 shape4      nu lambda

        OUTPUT:

        A pair ``(t, content)`` consisting of the shape of the fourth
        partition acording to the Burge correspondence and the content of
        the cell.

        TEST::

            sage: w = [4,1,8,3,6,5,2,7,9]; G = GrowthDiagramBurge(w)
            sage: GrowthDiagramBurge(labels=G._out_labels).to_word() == w       # indirect doctest
            True

            sage: G = GrowthDiagramBurge
            sage: G._backward_rule([1,1,1],[2,1,1,1],[2,1,1])
            ([1, 1], 0)

        """
        carry = 0
        shape2 = []
        i = len(shape4)
        while i > 0:
            mu_i = 0 if len(shape1) < i else shape1[i-1]
            la_i = 0 if len(shape4) < i else shape4[i-1]
            nu_i = 0 if len(shape3) < i else shape3[i-1]

            shape2 = [min(mu_i, nu_i) - min(int(mu_i == nu_i == la_i), carry)] + shape2
            carry = carry - min(int(mu_i == nu_i == la_i), carry) + la_i - max(mu_i, nu_i)
            i = i-1
        return (Partition(shape2), carry)

class GrowthDiagramDomino(GrowthDiagramOnPartitions):
    r"""
    A class modelling domino insertion.

    EXAMPLES:

    Figure 3 in [Lam2004]_::

        sage: G = GrowthDiagramDomino([[0,0,0,-1],[0,0,1,0],[-1,0,0,0],[0,1,0,0]]); G
         0  0  0 -1
         0  0  1  0
        -1  0  0  0
         0  1  0  0

        sage: ascii_art(G.P_symbol(), G.Q_symbol())
        1  2  4  1  2  2
        1  2  4  1  3  3
        3  3     4  4

    .. automethod:: _forward_rule
    """
    @staticmethod
    def _forward_rule(shape3, shape2, shape1, content):
        r"""
        Return the output shape given three shapes and the content.

        See [Lam2004]_ Section 3.1.

        INPUT:

        - ``shape3, shape2, shape1`` -- three partitions from a cell in a
          growth diagram, labelled as::

              shape2 shape1
              shape3

        - ``content`` -- -1, 0 or 1, the content of the cell.

        OUTPUT:

        The fourth partition according to domino insertion.

        TESTS::

            sage: G = GrowthDiagramDomino

        Rule 1::

            sage: G._forward_rule([], [], [], 1)
            [2]

            sage: G._forward_rule([1,1], [1,1], [1,1], 1)
            [3, 1]

        Rule 2::

            sage: G._forward_rule([1,1], [1,1], [1,1], -1)
            [2, 2]

        Rule 3::

            sage: G._forward_rule([1,1], [1,1], [2,2], 0)
            [2, 2]

        Rule 4::

            sage: G._forward_rule([2,2,2], [2,2], [3,3], 0)
            [3, 3, 2]

            sage: G._forward_rule([2], [], [1,1], 0)
            [2, 2]

            sage: G._forward_rule([1,1], [], [2], 0)
            [2, 2]

            sage: G._forward_rule([2], [], [2], 0)
            [2, 2]

            sage: G._forward_rule([4], [2], [4], 0)
            [4, 2]

            sage: G._forward_rule([1,1,1,1], [1,1], [1,1,1,1], 0)
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
