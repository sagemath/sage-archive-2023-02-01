r"""
Combinatorial diagrams

A combinatorial diagram is a collection of cells `(i,j)` indexed by pairs of
natural numbers. The positions are indexed by rows and columns. For example,
a Ferrer's diagram is a diagram obtained from a partition
`\lambda = (\lambda_0, \lambda_1, \ldots, \lambda_\ell)` where the cells are
in rows `i` for `0 \leq i \leq \ell` and the cells in row `i` consist of
`(i,j)` for `0 \leq j < \lambda_i`. In English notation, the indices are read
from top left to bottom right as in a matrix.

Indexing conventions are the same as
:class:`~sage.combinat.partition.Partition`.

EXAMPLES:

To create an arbirtrary diagram, pass a list of all cells::

    sage: from sage.combinat.diagram import Diagram
    sage: cells = [(0,0), (0,1), (1,0), (1,1), (4,4), (4,5), (4,6), (5,4), (7, 6)]
    sage: D = Diagram(cells); D
    [(0, 0), (0, 1), (1, 0), (1, 1), (4, 4), (4, 5), (4, 6), (5, 4), (7, 6)]

We can visualize the diagram by printing ``O``'s and ``.``'s. ``O``'s are
present in the cells which are present in the diagram and a ``.`` represents
the absence of a cell in the diagram::

    sage: D.pp()
    O O . . . . .
    O O . . . . .
    . . . . . . .
    . . . . . . .
    . . . . O O O
    . . . . O . .
    . . . . . . .
    . . . . . . O

We can also check if certain cells are contained in a given diagram::

    sage: (1, 0) in D
    True
    sage: (2, 2) in D
    False

If you know that there are entire empty rows or columns at the end of the
diagram, you can manually pass them with keyword arguments ``n_rows=`` or
``n_cols=``::

    sage: Diagram([(0,0), (0,3), (2,2), (2,4)]).pp()
    O . . O .
    . . . . .
    . . O . O
    sage: Diagram([(0,0), (0,3), (2,2), (2,4)], n_rows=6, n_cols=6).pp()
    O . . O . .
    . . . . . .
    . . O . O .
    . . . . . .
    . . . . . .
    . . . . . .


A ``Diagram`` is an element of the parent class ``Diagrams``, so we can also
construct a diagram by calling an instance of ``Diagrams``::

    sage: from sage.combinat.diagram import Diagrams
    sage: Dgms = Diagrams()
    sage: D in Dgms
    True
    sage: Dgms([(0,1),(3,3)]).pp()
    . O . .
    . . . .
    . . . .
    . . . O

There are two other specific types of diagrams which are implemented in Sage,
namely northwest diagrams (:class:`NorthwestDiagram`) and Rothe diagrams
(:func:`RotheDiagram`, a special kind of northwest diagram).

A diagram is a
*northwest diagram* if it satsifies the property that: the presence of two
cells `(i_1, j_1)` and `(i_2, j_2)` in a diagram `D` implies the presence of
the cell `(\min(i_1, i_2), \min(j_1, j_2))`::

    sage: from sage.combinat.diagram import NorthwestDiagram
    sage: N = NorthwestDiagram([(0,0), (0, 10), (5,0)]); N.pp()
    O . . . . . . . . . O
    . . . . . . . . . . .
    . . . . . . . . . . .
    . . . . . . . . . . .
    . . . . . . . . . . .
    O . . . . . . . . . .

Note that checking whether or not the northwest property is satisfied is
automatically checked. The diagram found by adding the cell `(1,1)` to the
diagram above is *not* a northwest diagram. The cell `(1,0)` should be
present due to the presence of `(5,0)` and `(1,1)`::

    sage: Diagram([(0, 0), (0, 10), (5, 0), (1, 1)]).pp()
    O . . . . . . . . . O
    . O . . . . . . . . .
    . . . . . . . . . . .
    . . . . . . . . . . .
    . . . . . . . . . . .
    O . . . . . . . . . .
    sage: NorthwestDiagram([(0, 0), (0, 10), (5, 0), (1, 1)])
    Traceback (most recent call last):
    ...
    AssertionError

However, this behavior can be turned off if you are confident that you are
providing a northwest diagram::

    sage: N = NorthwestDiagram([(0, 0), (0, 10), (5, 0),
    ....:                      (1, 1), (0, 1), (1, 0)],
    ....:                      check=False)
    sage: N.pp()
    O O . . . . . . . . O
    O O . . . . . . . . .
    . . . . . . . . . . .
    . . . . . . . . . . .
    . . . . . . . . . . .
    O . . . . . . . . . .

Note that arbitrary diagrams which happen to be northwest diagrams only live
in the parent of :class:`Diagrams`::

    sage: D = Diagram([(0, 0), (0, 10), (5, 0), (1, 1), (0, 1), (1, 0)])
    sage: D.pp()
    O O . . . . . . . . O
    O O . . . . . . . . .
    . . . . . . . . . . .
    . . . . . . . . . . .
    . . . . . . . . . . .
    O . . . . . . . . . .
    sage: from sage.combinat.diagram import NorthwestDiagrams
    sage: D in NorthwestDiagrams()
    False

One particular way of constructing a northwest diagram from a permutation is
by constructing its Rothe diagram. Formally, if `\omega` is a permutation,
then the Rothe diagram `D(\omega)` is the diagram whose cells are

.. MATH::

    D(\omega) = \{(\omega_j, i) : i<j,\, \omega_i > \omega_j \}.

Informally, one can construct the Rothe diagram by starting with all `n^2`
possible cells, and then deleting the cells `(i, \omega(i))` as well as all
cells to the right and below. (These are the "death rays" of []_.) To compute
a Rothe diagram in Sage, start with a
:class:`~sage.combinat.permutation.Permutation`::

    sage: w = Permutations(8)([2,5,4,1,3,6,7,8])

Then call :func:`RotheDiagram` on ``w``::

    sage: from sage.combinat.diagram import RotheDiagram
    sage: RotheDiagram(w).pp()
    O . . . . . . .
    O . O O . . . .
    O . O . . . . .
    . . . . . . . .
    . . . . . . . .
    . . . . . . . .
    . . . . . . . .
    . . . . . . . .

For a fixed northwest diagram `D`, we say that a Young tableau `T` is
`D`-peelable if:

- the row indices of the cells in the first column of `D` are
the entries in an initial segment in the first column of `T` and

- the tableau `Q` obtained by removing those cells from `T` and playing jeu de
taquin is `D-C`-peelable, where `D-C` is the diagram formed by forgetting the
first column of `D`.

Reiner and Shimozono [RS1995]_ showed that the number `\operatorname{red}(w)`
of reduced words of a permutation `w` may be computed using the peelable
tableaux of the Rothe diagram `D(w)`. Explicitly,

.. MATH::

    \operatorname{red}(w) = \sum_{T} f_{\operatorname{shape} T}

where the sum runs over the `D(w)`-peelable tableaux `T` and `f_\lambda` is the
number of standard Young tableaux of shape `\lambda` (which may be computed
using the hook-length formula).

We can compute the `D`-peelable diagrams for a northwest diagram `D`::

    sage: cells = [(0,0), (0,1), (0,2), (1,0), (2,0), (2,2), (2,4),
    ....:          (4,0), (4,2)]
    sage: D = NorthwestDiagram(cells); D.pp()
    O O O . .
    O . . . .
    O . O . O
    . . . . .
    O . O . .
    sage: D.peelable_tableaux()
    {[[1, 1, 1], [2, 3, 3], [3, 5], [5]],
    [[1, 1, 1, 3], [2, 3], [3, 5], [5]]}

AUTHORS:

- Trevor K. Karn (2022-08-01): initial version
"""

# ****************************************************************************
#       Copyright (C) 2022 Trevor K. Karn <karnx018 (at) umn.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories.sets_cat import Sets
from sage.sets.non_negative_integers import NonNegativeIntegers as NN
from sage.combinat.partition import Partition
from sage.combinat.permutation import Permutations
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.list_clone import ClonableArray
from sage.structure.parent import Parent
from sage.combinat.tableau import Tableau
from sage.combinat.skew_tableau import SkewTableaux
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass


class Diagram(ClonableArray, metaclass=InheritComparisonClasscallMetaclass):
    r"""
    A class to model arbitrary combinatorial diagrams.

    EXAMPLES::

        sage: from sage.combinat.diagram import Diagram
        sage: D = Diagram([(0,0), (0,3), (2,2), (2,4)]); D
        [(0, 0), (0, 3), (2, 2), (2, 4)]

    TESTS::

        sage: from sage.combinat.diagram import Diagrams
        sage: D = Diagrams().an_element()
        sage: TestSuite(D).run()
    """
    @staticmethod
    def __classcall_private__(self, cells, n_rows=None, n_cols=None, check=False):
        r"""
        Normalize the input so that it lives in the correct parent.

        EXAMPLES::

            sage: from sage.combinat.diagram import Diagram
            sage: D = Diagram([(0,0), (0,3), (2,2), (2,4)])
            sage: D.parent()
            Combinatorial diagrams
        """
        return Diagrams()(cells, n_rows, n_cols, check)

    def __init__(self, parent, cells, n_rows=None, n_cols=None, check=False):
        r"""
        EXAMPLES::

            sage: from sage.combinat.diagram import Diagram
            sage: D1 = Diagram([(0,2),(0,3),(1,1),(3,2)])
            sage: D1.cells()
            [(0, 2), (0, 3), (1, 1), (3, 2)]
            sage: D1.n_rows()
            4
            sage: D1.n_cols()
            4

        We can specify the number of rows and columns explicitly,
        in case they are supposed to be empty::

            sage: D2 = Diagram([(0,2),(0,3),(1,1),(3,2)], n_cols=5)
            sage: D2.cells()
            [(0, 2), (0, 3), (1, 1), (3, 2)]
            sage: D2.n_cols()
            5
            sage: D2.pp()
            . . O O .
            . O . . .
            . . . . .
            . . O . .
        """
        self._cells = {c: True for c in cells}
        if n_rows is not None:
            if n_rows < max(c[0] for c in self._cells) + 1:
                raise ValueError('n_rows is too small')
            self._n_rows = n_rows
        else:
            self._n_rows = max(c[0] for c in self._cells) + 1
        if n_cols is not None:
            if n_cols < max(c[1] for c in self._cells) + 1:
                raise ValueError('n_cols is too small')
            self._n_cols = n_cols
        else:
            self._n_cols = max(c[1] for c in self._cells) + 1

        self._n_nonempty_rows = len(set(i for i, j in self._cells))
        self._n_nonempty_cols = len(set(j for i, j in self._cells))

        ClonableArray.__init__(self, parent, cells, check=False)

    def _hash_(self):
        r"""
        TESTS::

            sage: from sage.combinat.diagram import Diagram
            sage: D = Diagram([(0,2),(0,3),(1,1),(3,2)])
            sage: hash(D)
            5125392318921082108
        """
        return hash(tuple(sorted(self._cells)))

    def __contains__(self, other):
        r"""
        EXAMPLES::

            sage: from sage.combinat.diagram import Diagram
            sage: D = Diagram([(0,2),(0,3),(1,1),(3,2)])
            sage: (1, 2) in D
            False
            sage: (0, 2) in D
            True
            sage: (2, 1) in D
            False
            sage: (3, 2) in D
            True
        """
        return other in self._cells

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: from sage.combinat.diagram import Diagram
            sage: D = Diagram([(0,2),(0,3),(1,1),(3,2)]); D
            [(0, 2), (0, 3), (1, 1), (3, 2)]
        """
        return str(list(self._cells))

    def pp(self):
        r"""
        Return a visualization of the diagram. Cells which are present in the
        diagram are filled with a ``O``. Cells which are not present in the
        diagram are filled with a ``.``.

        EXAMPLES::

            sage: from sage.combinat.diagram import Diagram
            sage: Diagram([(0,0), (0,3), (2,2), (2,4)]).pp()
            O . . O .
            . . . . .
            . . O . O
            sage: Diagram([(0,0), (0,3), (2,2), (2,4)], n_rows=6, n_cols=6).pp()
            O . . O . .
            . . . . . .
            . . O . O .
            . . . . . .
            . . . . . .
            . . . . . .
        """
        output_str = ''

        for i in range(self._n_rows):
            for j in range(self._n_cols):
                if (i, j) in self:
                    output_str += 'O '
                else:
                    output_str += '. '
            output_str += '\n'

        print(output_str)

    def n_rows(self):
        r"""
        Return the total number of rows of the cell, including those which do
        not have any cells.

        EXAMPLES:

        The following example has three rows which are filled, but they
        are contained in rows 0 to 3 (for a total of four)::

            sage: from sage.combinat.diagram import Diagram
            sage: D1 = Diagram([(0,2),(0,3),(1,1),(3,2)])
            sage: D1.n_rows()
            4

        We can also include empty rows at the end::

            sage: from sage.combinat.diagram import Diagram
            sage: D = Diagram([(0,2),(0,3),(1,1),(3,2)], n_rows=6)
            sage: D.n_rows()
            6
            sage: D.pp()
            . . O O
            . O . .
            . . . .
            . . O .
            . . . .
            . . . .
        """
        return self._n_rows

    def n_cols(self):
        r"""
        Return the total number of rows of the cell, including those which do
        not have any cells.

        EXAMPLES:

        The following example has three columns which are filled, but they
        are contained in rows 0 to 3 (for a total of four)::

            sage: from sage.combinat.diagram import Diagram
            sage: D = Diagram([(0,2),(0,3),(1,1),(3,2)])
            sage: D.n_cols()
            4

        We can also include empty columns at the end::

            sage: from sage.combinat.diagram import Diagram
            sage: D = Diagram([(0,2),(0,3),(1,1),(3,2)], n_cols=6)
            sage: D.n_cols()
            6
            sage: D.pp()
            . . O O . .
            . O . . . .
            . . . . . .
            . . O . . .
        """

        return self._n_cols

    def cells(self):
        r"""
        Return a ``list`` of the cells contained in the diagram ``self``.

        EXAMPLES::

            sage: from sage.combinat.diagram import Diagram
            sage: D1 = Diagram([(0,2),(0,3),(1,1),(3,2)])
            sage: D1.cells()
            [(0, 2), (0, 3), (1, 1), (3, 2)]
        """
        return list(self._cells.keys())

    def n_cells(self):
        r"""
        Return the total number of cells contained in the diagram ``self``.

        EXAMPLES::

            sage: from sage.combinat.diagram import Diagram
            sage: D1 = Diagram([(0,2),(0,3),(1,1),(3,2)])
            sage: D1.n_cells()
            4
        """
        return len(self._cells)

    def check(self):
        r"""
        Check that this is a valid diagram by checking that it is an iterable
        of length-two tuples. (The fact that each tuple is length-two is
        implicitly checked during creation of the diagram).

        .. WARNING::

            This method is required for ``Diagram`` to be a subclass of
            :class:`~sage.structure.list_clone.ClonableArray`, however the check
            is *not* automatically performed upon creation of the element.

        EXAMPLES::

            sage: from sage.combinat.diagram import Diagram
            sage: D = Diagram([(0,0), (0,3), (2,2), (2,4)])
            sage: D.check()

        In the next two examples, the diagram ``D`` is initialized but with
        bad information that is not detected untill we call ``D.check()``.
        The first example fails because one cells is indexed by negative
        integers::

            sage: D = Diagram([(0,0), (0,-3), (2,2), (2,4)])
            sage: D.check()
            Traceback (most recent call last):
            ...
            AssertionError

        The next example fails because one cell is indexed by rational
        numbers::

            sage: D = Diagram([(0,0), (0,3), (2/3,2), (2,4)])
            sage: D.check()
            Traceback (most recent call last):
            ...
            AssertionError
        """
        from sage.sets.non_negative_integers import NonNegativeIntegers
        NN = NonNegativeIntegers()
        assert all(all(list(i in NN for i in c)) for c in self._cells.keys())


class Diagrams(UniqueRepresentation, Parent):
    r"""
    The class of combinatorial diagrams.

    A *combinatorial diagram* is a set of cells indexed by pairs of natural
    numbers. Calling an instance of :class:`Diagrams` is one way to construct
    diagrams.

    EXAMPLES::

        sage: from sage.combinat.diagram import Diagrams
        sage: Dgms = Diagrams()
        sage: D = Dgms([(0,0), (0,3), (2,2), (2,4)])
        sage: D.parent()
        Combinatorial diagrams

    """
    def __init__(self):
        r"""
        EXAMPLES::

            sage: from sage.combinat.diagram import Diagrams
            sage: Dgms = Diagrams(); Dgms
            Combinatorial diagrams

        TESTS::

            sage: TestSuite(Dgms).run()
        """

        Parent.__init__(self, category=Sets())

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: from sage.combinat.diagram import Diagrams
            sage: Dgms = Diagrams(); Dgms
            Combinatorial diagrams
        """
        return 'Combinatorial diagrams'

    def _element_constructor_(self, cells, n_rows=None, n_cols=None, check=False):
        r"""
        EXAMPLES::

            sage: from sage.combinat.diagram import Diagrams
            sage: Dgms = Diagrams()
            sage: Dgms([(0,1),(2,2)]).pp()
            . O .
            . . .
            . . O

        TESTS::

            sage: TestSuite(Dgms).run()
        """
        return self.element_class(self, cells, n_rows, n_cols, check)

    def _an_element_(self):
        r"""
        EXAMPLES::

            sage: from sage.combinat.diagram import Diagrams
            sage: Dgms = Diagrams()
            sage: D = Dgms.an_element(); D
            [(0, 2), (1, 1), (2, 3)]
            sage: D.pp()
            . . O .
            . O . .
            . . . O
        """
        return self([(0, 2), (1, 1), (2, 3)])

    Element = Diagram


####################
# Northwest diagrams
####################

class NorthwestDiagram(Diagram, metaclass=InheritComparisonClasscallMetaclass):
    r"""
    A diagram is a set of cells indexed by natural numbers. Such a diagram
    has the *northwest property* if the presence of cells `(i1, j1)` and
    `(i2, j2)` implies the presence of the cell
    `(\min(i1, i2), \min(j1, j2))`. Diagrams with the northwest property are
    called *northwest diagrams*.

    For general diagrams see :class:`Diagram`.

    EXAMPLES::

        sage: from sage.combinat.diagram import NorthwestDiagram
        sage: N = NorthwestDiagram([(0,0), (0, 2), (2,0)])

    To visualize them, use the ``.pp()`` method::

        sage: N.pp()
        O . O
        . . .
        O . .

    One can also create northwest diagrams from a partition::

        sage: mu = Partition([5,4,3,2,1])
        sage: mu.pp()
        *****
        ****
        ***
        **
        *
        sage: D = NorthwestDiagram(mu.cells()).pp()
        O O O O O
        O O O O .
        O O O . .
        O O . . .
        O . . . .
    """
    @staticmethod
    def __classcall_private__(self, cells, n_rows=None, n_cols=None, check=True):
        """
        Normalize input to ensure a correct parent.

        EXAMPLES::

            sage: from sage.combinat.diagram import NorthwestDiagram, NorthwestDiagrams
            sage: N1 = NorthwestDiagram([(0,1), (0,2)])
            sage: N2 = NorthwestDiagram([(0,1), (0,3)])
            sage: N1.parent() is N2.parent()
            True
            sage: N3 = NorthwestDiagrams()([(0,1), (0,2)])
            sage: N3.parent() is NorthwestDiagrams()
            True
            sage: N1.parent() is NorthwestDiagrams()
            True
        """
        # TODO: Assert that cells is sorted in lex order to speed up lookup.
        return NorthwestDiagrams()(cells, n_rows, n_cols, check)

    def check(self):
        r"""
        A diagram has the northwest property if the presence of cells
        `(i1, j1)` and `(i2, j2)` implies the presence of the cell
        `(min(i1, i2), min(j1, j2))`. This method checks if the northwest
        property is satisfied for ``self``

        EXAMPLES::

            sage: from sage.combinat.diagram import NorthwestDiagram
            sage: N = NorthwestDiagram([(0,0), (0,3), (3,0)])
            sage: N.check()

        Here is a non-example::

            sage: notN = NorthwestDiagram([(0,1), (1,0)])  #.check() is implicit
            Traceback (most recent call last):
            ...
            AssertionError
        """
        from itertools import combinations
        assert all((min(i1, i2), min(j1, j2)) in self
                   for (i1, j1), (i2, j2) in combinations(self._cells, 2))

    def peelable_tableaux(self):
        r"""
        Given a northwest diagram `D`, a tableau `T` is said to be
        `D`-peelable if...

        EXAMPLES:

        If the diagram is only one column, there is only one peelable tableau::

            sage: from sage.combinat.diagram import NorthwestDiagram
            sage: NWD = NorthwestDiagram([(0,0), (2,0)])
            sage: NWD.peelable_tableaux()
            {[[1], [3]]}

        From [...]_ we know that there is only one peelable tableau for the
        Rothe diagram of the permutation (in one line notation) `251643`::

            sage: D = NorthwestDiagram([(1, 2), (1, 3), (3, 2), (3, 3), (4, 2)])
            sage: D.pp()
            . . . .
            . . O O
            . . . .
            . . O O
            . . O .

            sage: D.peelable_tableaux()
            {[[2, 2], [4, 4], [5]]}

        Here are all the intermediate steps to compute the peelables for the
        Rothe diagram of (in one-line notation) `64817235`. They are listed from
        deepest in the recursion to the final step. The recursion has depth five
        in this case so we will label the intermediate tableaux by `D_i` where
        `i` is the step in the recursion at which they appear.

        Start with the one that has a single column::

            sage: D5 = NorthwestDiagram([(2,0)]); D5.pp()
            .
            .
            O
            sage: D5.peelable_tableaux()
            {[[3]]}

        Now we know all of the `D_5` peelables, so we can compute the `D_4`
        peelables.

            sage: D4 = NorthwestDiagram([(0, 0), (2,0), (4, 0), (2, 2)])
            sage: D4.pp()
            O . .
            . . .
            O . O
            . . .
            O . .

            sage: D4.peelable_tableaux()
            {[[1, 3], [3], [5]]}

        There is only one `D_4` peelable, so we can compute the `D_3`
        peelables::

            sage: D3 = NorthwestDiagram([(0,0), (0,1), (2, 1), (2, 3), (4,1)])
            sage: D3.pp()
            O O . .
            . . . .
            . O . O
            . . . .
            . O . .

            sage: D3.peelable_tableaux()
            {[[1, 1], [3, 3], [5]], [[1, 1, 3], [3], [5]]}

        Now compute the `D_2` peelables::

            sage: cells = [(0,0), (0,1), (0,2), (1,0), (2,0), (2,2), (2,4),
            ....:          (4,0), (4,2)]
            sage: D2 = NorthwestDiagram(cells); D2.pp()
            O O O . .
            O . . . .
            O . O . O
            . . . . .
            O . O . .

            sage: D2.peelable_tableaux()
            {[[1, 1, 1], [2, 3, 3], [3, 5], [5]],
            [[1, 1, 1, 3], [2, 3], [3, 5], [5]]}

        And the `D_1` peelables::

            sage: cells = [(0,0), (0,1), (0,2), (0,3), (1,0), (1,1), (2,0),
            ....:          (2,1), (2,3), (2,5), (4,0), (4,1), (4,3)]
            sage: D1 = NorthwestDiagram(cells); D1.pp()
            O O O O . .
            O O . . . .
            O O . O . O
            . . . . . .
            O O . O . .

            sage: D1.peelable_tableaux()
            {[[1, 1, 1, 1], [2, 2, 3, 3], [3, 3, 5], [5, 5]],
            [[1, 1, 1, 1, 3], [2, 2, 3], [3, 3, 5], [5, 5]]}

        Which we can use to get the `D` peelables::

            sage: cells = [(0,0), (0,1), (0,2), (0,3), (0,4),
            ....:          (1,0), (1,1), (1,2),
            ....:          (2,0), (2,1), (2,2), (2,4), (2,6),
            ....:                 (4,1), (4,2), (4,4)]
            sage: D = NorthwestDiagram(cells); D.pp()
            O O O O O . .
            O O O . . . .
            O O O . O . O
            . . . . . . .
            . O O . O . .
            sage: D.peelable_tableaux()
            {[[1, 1, 1, 1, 1], [2, 2, 2, 3, 3], [3, 3, 3], [5, 5, 5]],
             [[1, 1, 1, 1, 1], [2, 2, 2, 3, 3], [3, 3, 3, 5], [5, 5]],
             [[1, 1, 1, 1, 1, 3], [2, 2, 2, 3], [3, 3, 3], [5, 5, 5]],
             [[1, 1, 1, 1, 1, 3], [2, 2, 2, 3], [3, 3, 3, 5], [5, 5]]}

        .. ALGORITHM::

            This implementation uses the algorithm suggested in remark 25
            of [...]_.
        """
        # TODO: There is a condition on the first column (if the rows in Dhat
        # are a subset of the rows in the first column) which simplifies the
        # description without performing JDT, so we should implement that

        # if there is a single column in the diagram then there is only
        # one posslbe peelable tableau.
        if self._n_nonempty_cols == 1:
            return {Tableau([[i+1] for i, j in self._cells])}

        first_col = min(j for i, j in self._cells)

        # TODO: The next two lines of code could be optimized by only
        # looping over self._cells once (rather than two separate times)
        # get the diagram without the first column
        Dhat = NorthwestDiagram([c for c in self._cells if c[1] != first_col])

        # get the values from the first column
        new_vals = sorted(i + 1 for i, j in self._cells if j == first_col)

        k = self.n_cells() - Dhat.n_cells()

        peelables = set()

        for Q in Dhat.peelable_tableaux():
            # get the vertical strips
            mu = Q.shape()
            vertical_strips = mu.add_vertical_border_strip(k)
            for s in vertical_strips:
                sQ = SkewTableaux()(Q)  # sQ is skew - get it?
                new_cells = (s/mu).cells()
                # perform the jeu de taquin slides
                for c in new_cells:
                    sQ = sQ.backward_slide(c)
                # create the new tableau by filling the columns
                sQ_new = sQ.to_list()
                for n in range(k):
                    sQ_new[n][0] = new_vals[n]

                T = Tableau(sQ_new)
                if T.is_column_strict():
                    peelables.add(T)

        return peelables


class NorthwestDiagrams(Diagrams):
    r"""
    The class of northwest diagrams.

    A diagram has the *northwest property* if the presence of cells
    `(i1, j1)` and `(i2, j2)` implies the presence of the cell
    `(min(i1, i2), min(j1, j2))`. For the class of general diagrams, see
    :class:`Diagrams`.

    One may create an instance of a northwest diagram either by directly
    calling :class:`NorthwestDiagram` or by creating an instance of the
    parent class `NorthwestDiagrams` (with an `s`) and calling it on a `list`
    of tuples consisting of the coordinates of the diagram.

    EXAMPLES::

    Additionally, there are natural constructions of a northwest diagram
    given the data of a permutation (Rothe diagrams are the protypical example
    of northwest diagrams), or the data of a partition of an integer, or a
    skew partition.

    The Rothe diagram `D(\omega)` of a permutation `\omega` is specified by
    the cells

    .. MATH::

        \{ (i, j) : i < \omega^{-1}(j) \text{ and } j < \omega(i)}

    We can construct one by calling :meth:`rothe_diagram` method on the parent
    class :class:`NorthwestDiagrams`.

    EXAMPLES::

        sage: from sage.combinat.diagram import NorthwestDiagrams
        sage: NWDgms = NorthwestDiagrams

    To turn a Ferrers diagram into a northwest diagram, we may call the
    :meth:`ferrers_diagram` method. This will return a Ferrer's diagram in the
    parent of all northwest diagrams. For many use-cases it is probably better
    to get Ferrer's diagrams by the corresponding method on partitons, namely
    :meth:`sage.combinat.partitions.Partitions.ferrers_diagram`.

    It is also possible to turn a Ferrers diagram of a skew partition into a
    northwest diagram, altough it is more subtle than just using the skew
    diagram itself.
    """

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: from sage.combinat.diagram import NorthwestDiagrams
            sage: NWDgms = NorthwestDiagrams(); NWDgms
            Combinatorial northwest diagrams
        """
        return 'Combinatorial northwest diagrams'

    def _an_element_(self):
        r"""
        EXAMPLES::

            sage: from sage.combinat.diagram import NorthwestDiagrams
            sage: NWDgms = NorthwestDiagrams()
            sage: NWD = NWDgms.an_element(); NWD
            [(0, 1), (0, 2), (1, 1), (2, 3)]
            sage: NWD.pp()
            . O O .
            . O . .
            . . . O
            sage: NWD.parent() is NWDgms
            True
        """
        return self([(0, 1), (0, 2), (1, 1), (2, 3)])

    Element = NorthwestDiagram


def RotheDiagram(w):
    r"""
    A constructor to build the Rothe diagram of a permutation w as an element
    in the parent of :class:`NorthwestDiagrams`.

    EXAMPLES::

        sage: w = Permutations(9)([1, 7, 4, 5, 9, 3, 2, 8, 6])
        sage: from sage.combinat.diagram import RotheDiagram
        sage: D = RotheDiagram(w); D.pp()
        . . . . . . . . .
        . O O O O O . . .
        . O O . . . . . .
        . O O . . . . . .
        . O O . . O . O .
        . O . . . . . . .
        . . . . . . . . .
        . . . . . O . . .
        . . . . . . . . .


    Currently, only elements of the parent
    :class:`sage.combinat.permutations.Permutations` are supported. In
    particular, elements of permutation groups are not supported::

        sage: w = SymmetricGroup(9).an_element()
        sage: RotheDiagram(w)
        Traceback (most recent call last):
        ...
        ValueError: w must be a Permutation
    """

    from sage.misc.mrange import cartesian_product_iterator

    if w not in Permutations():
        raise ValueError('w must be a Permutation')

    N = w.size()

    cells = [c for c in cartesian_product_iterator((range(N), range(N)))
             if c[0]+1 < w.inverse()(c[1]+1) and c[1]+1 < w(c[0]+1)]

    return NorthwestDiagram(cells, n_rows=N, n_cols=N, check=False)
