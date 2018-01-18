# -*- coding: utf-8 -*-
"""
Shifted primed tableaux

AUTHORS:

- Kirill Paramonov (2017-08-18): initial implementation
"""

#*****************************************************************************
#       Copyright (C) 2017 Kirill Paramonov <kbparamonov at ucdavis.edu>,
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import print_function, absolute_import
from six import add_metaclass

from sage.combinat.partition import Partition, Partitions, _Partitions, OrderedPartitions
from sage.combinat.tableau import Tableaux
from sage.combinat.skew_partition import SkewPartition
from sage.combinat.integer_vector import IntegerVectors
from sage.rings.integer import Integer
from sage.rings.rational_field import QQ
from sage.rings.rational import Rational
from sage.rings.integer_ring import ZZ

from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.misc.lazy_attribute import lazy_attribute

from sage.structure.list_clone import ClonableArray
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation

from sage.categories.regular_crystals import RegularCrystals
from sage.categories.classical_crystals import ClassicalCrystals
from sage.categories.sets_cat import Sets
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.combinat.root_system.cartan_type import CartanType


@add_metaclass(InheritComparisonClasscallMetaclass)
class ShiftedPrimedTableau(ClonableArray):
    r"""
    A shifted primed tableau with primed elements stored as half-integers.

    A primed tableau is a tableau of shifted shape in the alphabet
    `X' = \{1' < 1 < 2' < 2 < \cdots < n' < n\}` such that

        1. The entries are weakly increasing along rows and columns;
        2. A row cannot have two repeated primed elements, and a column
           cannot have two repeated non-primed elements;
        3. There are only non-primed elements on the main diagonal.

    EXAMPLES::

        sage: T = ShiftedPrimedTableaux([4,2])
        sage: T([[1,"2'","3'",3],[2,"3'"]])[1]
        (2, 3')
        sage: t = ShiftedPrimedTableau([[1,"2p",2.5,3],[2,2.5]])
        sage: t[1]
        (2, 3')
        sage: ShiftedPrimedTableau([["2p",2,3],["2p","3p"],[2]], skew=[2,1])
        [(None, None, 2', 2, 3), (None, 2', 3'), (2,)]

    TESTS::

        sage: t = ShiftedPrimedTableau([[1,2,2.5,3],[2,2.5]])
        Traceback (most recent call last):
        ...
        ValueError: [[1, 2, 2.50000000000000, 3], [2, 2.50000000000000]]
         is not an element of Shifted Primed Tableaux
    """
    @staticmethod
    def __classcall_private__(cls, T, skew=None):
        r"""
        Ensure that a shifted tableau is only ever constructed as an
        ``element_class`` call of an appropriate parent.

        EXAMPLES::

            sage: data = [[1,"2'","2",3],[2,"3'"]]
            sage: t = ShiftedPrimedTableau(data)
            sage: T = ShiftedPrimedTableaux(shape=[4,2],weight=(1,3,2))
            sage: t == T(data)
            True
            sage: S = ShiftedPrimedTableaux(shape=[4,2])
            sage: t == S(data)
            True
            sage: t = ShiftedPrimedTableau([["2p",2,3],["2p"]],skew=[2,1])
            sage: t.parent()
            Shifted Primed Tableaux skewed by [2, 1]
        """
        if (isinstance(T, cls) and T._skew == skew):
            return T
        return ShiftedPrimedTableaux(skew=skew)(T)

    def __init__(self, parent, T, skew=None, check=True, preprocessed=False):
        r"""
        Initialize a shifted tableau.

        TESTS::

            sage: s = ShiftedPrimedTableau([[1,"2'","3'",3],[2,"3'"]])
            sage: t = ShiftedPrimedTableaux([4,2])([[1,"2p","3p",3],
            ....: [2,"3p"]])
            sage: s==t
            True
            sage: t.parent()
            Shifted Primed Tableaux of shape [4, 2]
            sage: s.parent()
            Shifted Primed Tableaux
            sage: r = ShiftedPrimedTableaux([4, 2])(s); r.parent()
            Shifted Primed Tableaux of shape [4, 2]
            sage: s is t # identical shifted tableaux are distinct objects
            False

        A shifted primed tableau is shallowly immutable, the rows are
        represented as tuples.

        ::

            sage: t = ShiftedPrimedTableau([[1,"2p","3p",3],[2,"3p"]])
            sage: t[0][1] = 3
            Traceback (most recent call last):
            ...
            TypeError: 'tuple' object does not support item assignment
        """
        if not preprocessed:
            T = self._preprocess_(T, skew=skew)
        self._skew = skew
        ClonableArray.__init__(self, parent, T, check=check)

    @staticmethod
    def _preprocess_(T, skew=None):
        """
        Preprocessing list ``T`` to initialize the tableau.
        """
        if isinstance(T, ShiftedPrimedTableau):
            return T

        # Preprocessing list t for primes and other symbols
        T = [[ShiftedPrimedTableauEntry(entry) for entry in row] for row in T]

        if skew is not None:
            T = ([(None,)*skew[i] + tuple(T[i])
                  for i in range(len(skew))]
                 + [tuple(T[i]) for i in range(len(skew), len(T))])
        else:
            T = [tuple(row) for row in T]
        return T

    def check(self):
        """
        Check that ``self`` is a valid primed tableau.

        EXAMPLES::

            sage: T = ShiftedPrimedTableaux([4,2])
            sage: t = T([[1,'2p',2,2],[2,'3p']])
            sage: t.check()
            sage: s = ShiftedPrimedTableau([["2p",2,3],["2p"],[2]],skew=[2,1])
            sage: s.check()
            sage: t = T([['1p','2p',2,2],[2,'3p']])
            Traceback (most recent call last):
            ....
            ValueError: [['1p', '2p', 2, 2], [2, '3p']] is not an element of
            Shifted Primed Tableaux of shape [4, 2]
        """
        if not self.parent()._contains_tableau_(self):
            raise ValueError("{} is not an element of Shifted Primed Tableaux"
                             .format(self))

    def __eq__(self, other):
        """
        Check whether ``self`` is equal to ``other``.

        INPUT:

        - ``other`` -- the element that ``self`` is compared to

        OUTPUT: Boolean

        EXAMPLES::

            sage: t = ShiftedPrimedTableau([[1,"2p"]])
            sage: t == ShiftedPrimedTableaux([2])([[1,1.5]])
            True
        """
        if isinstance(other, ShiftedPrimedTableau):
            return (self._skew == other._skew and list(self) == list(other))
        try:
            Tab = ShiftedPrimedTableau(other)
        except (ValueError, TypeError):
            return False
        return (self._skew == Tab._skew and list(self) == list(Tab))

    def __ne__(self, other):
        """
        Check whether ``self`` is not equal to ``other``.

        INPUT:

        - ``other`` -- the element that ``self`` is compared to

        OUTPUT: Boolean

        EXAMPLES::

            sage: t = ShiftedPrimedTableau([[1,"2p"]])
            sage: t != ShiftedPrimedTableaux([2])([[1,1]])
            True
        """
        if isinstance(other, ShiftedPrimedTableau):
            return (self._skew != other._skew or list(self) != list(other))
        try:
            Tab = ShiftedPrimedTableau(other)
        except (ValueError, TypeError):
            return True
        return (self._skew != Tab._skew or list(self) != list(Tab))

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: t = ShiftedPrimedTableau([[1,'2p',2,2],[2,'3p']])
            sage: t
            [(1, 2', 2, 2), (2, 3')]
            sage: ShiftedPrimedTableau([["2p",2,3],["2p"]],skew=[2,1])
            [(None, None, 2', 2, 3), (None, 2')]
            """
        return self.parent().options._dispatch(self, '_repr_', 'display')

    def _repr_list(self):
        """
        Return a string representation of ``self`` as a list of tuples.

        EXAMPLES::

            sage: print(ShiftedPrimedTableau([['2p',3],[2,2]], skew=[2]))
            [(None, None, 2', 3), (2, 2)]
            """
        return (repr([row for row in self]))

    def _repr_tab(self):
        """
        Return a nested list of strings representing the elements.

        TESTS::

            sage: t = ShiftedPrimedTableau([[1,'2p',2,2],[2,'3p']])
            sage: t._repr_tab()
            [[' 1 ', " 2'", ' 2 ', ' 2 '], [' 2 ', " 3'"]]
            sage: s = ShiftedPrimedTableau([["2p",2,3],["2p"]],skew=[2,1])
            sage: s._repr_tab()
            [[' . ', ' . ', " 2'", ' 2 ', ' 3 '], [' . ', " 2'"]]
        """
        max_len = len(str(self.max_entry()))+2
        repr_tab = []
        for row in self:
            repr_row = []
            for entry in row:
                if entry is None:
                    repr_row.append('. '.rjust(max_len))
                elif entry.is_primed():
                    repr_row.append(repr(entry).rjust(max_len))
                elif entry.is_unprimed():
                    repr_row.append(repr(entry).rjust(max_len-1)+" ")
            repr_tab.append(repr_row)
        return repr_tab

    def _repr_diagram(self):
        """
        Return a string representation of ``self`` as an array.

        EXAMPLES::

            sage: t = ShiftedPrimedTableau([[1,'2p',2,2],[2,'3p']])
            sage: print(t._repr_diagram())
             1  2' 2  2
                2  3'
            sage: t = ShiftedPrimedTableau([[10,'11p',11,11],[11,'12']])
            sage: print(t._repr_diagram())
             10  11' 11  11
                 11  12
            sage: s = ShiftedPrimedTableau([["2p",2,3],["2p"]],skew=[2,1])
            sage: print(s._repr_diagram())
             .  .  2' 2  3
                .  2'
        """
        max_len = len(str(self.max_entry()))+2
        return "\n".join([" "*max_len*i + "".join(_)
                          for i, _ in enumerate(self._repr_tab())])

    _repr_compact = _repr_diagram

    def _ascii_art_(self):
        """
        Return ASCII representation of ``self``.

        EXAMPLES::

            sage: ascii_art(ShiftedPrimedTableau([[1,'2p',2,2],[2,'3p']]))
            +---+---+---+---+
            | 1 | 2'| 2 | 2 |
            +---+---+---+---+
                | 2 | 3'|
                +---+---+
            sage: s = ShiftedPrimedTableau([["2p",2,3],["2p"]],skew=[2,1])
            sage: ascii_art(s)
            +---+---+---+---+---+
            | . | . | 2'| 2 | 3 |
            +---+---+---+---+---+
                | . | 2'|
                +---+---+

        TESTS::

            sage: ascii_art(ShiftedPrimedTableau([]))
            ++
            ++
        """
        from sage.typeset.ascii_art import AsciiArt
        return AsciiArt(self._ascii_art_table(unicode=False).splitlines())

    def _unicode_art_(self):
        """
        Return a Unicode representation of ``self``.

        EXAMPLES::

            sage: unicode_art(ShiftedPrimedTableau([[1,'2p',2,2],[2,'3p']]))
            ┌───┬───┬───┬───┐
            │ 1 │ 2'│ 2 │ 2 │
            └───┼───┼───┼───┘
                │ 2 │ 3'│
                └───┴───┘
            sage: s = ShiftedPrimedTableau([["2p",2,3],["2p"]],skew=[2,1])
            sage: unicode_art(s)
            ┌───┬───┬───┬───┬───┐
            │ . │ . │ 2'│ 2 │ 3 │
            └───┼───┼───┼───┴───┘
                │ . │ 2'│
                └───┴───┘

        TESTS::

            sage: unicode_art(ShiftedPrimedTableau([]))
            ┌┐
            └┘
        """
        from sage.typeset.unicode_art import UnicodeArt
        return UnicodeArt(self._ascii_art_table(unicode=True).splitlines())

    def _ascii_art_table(self, unicode=False):
        """
        TESTS::

            sage: t = ShiftedPrimedTableau([[1,'2p',2],[2,'3p']])
            sage: print(t._ascii_art_table(unicode=True))
            ┌───┬───┬───┐
            │ 1 │ 2'│ 2 │
            └───┼───┼───┤
                │ 2 │ 3'│
                └───┴───┘
            sage: print(t._ascii_art_table())
            +---+---+---+
            | 1 | 2'| 2 |
            +---+---+---+
                | 2 | 3'|
                +---+---+
            sage: s = ShiftedPrimedTableau([[1,'2p',2, 23],[2,'30p']])
            sage: print(s._ascii_art_table(unicode=True))
            ┌────┬────┬────┬────┐
            │  1 │  2'│  2 │ 23 │
            └────┼────┼────┼────┘
                 │  2 │ 30'│
                 └────┴────┘
            sage: print(s._ascii_art_table(unicode=False))
            +----+----+----+----+
            |  1 |  2'|  2 | 23 |
            +----+----+----+----+
                 |  2 | 30'|
                 +----+----+
            sage: s = ShiftedPrimedTableau([["2p",2,10],["2p"]],skew=[2,1])
            sage: print(s._ascii_art_table(unicode=True))
            ┌────┬────┬────┬────┬────┐
            │  . │  . │  2'│  2 │ 10 │
            └────┼────┼────┼────┴────┘
                 │  . │  2'│
                 └────┴────┘
        """
        if unicode:
            import unicodedata
            v = unicodedata.lookup('BOX DRAWINGS LIGHT VERTICAL')
            h = unicodedata.lookup('BOX DRAWINGS LIGHT HORIZONTAL')
            dl = unicodedata.lookup('BOX DRAWINGS LIGHT DOWN AND LEFT')
            dr = unicodedata.lookup('BOX DRAWINGS LIGHT DOWN AND RIGHT')
            ul = unicodedata.lookup('BOX DRAWINGS LIGHT UP AND LEFT')
            ur = unicodedata.lookup('BOX DRAWINGS LIGHT UP AND RIGHT')
            vl = unicodedata.lookup('BOX DRAWINGS LIGHT VERTICAL AND LEFT')
            uh = unicodedata.lookup('BOX DRAWINGS LIGHT UP AND HORIZONTAL')
            dh = unicodedata.lookup('BOX DRAWINGS LIGHT DOWN AND HORIZONTAL')
            vh = unicodedata.lookup(
                'BOX DRAWINGS LIGHT VERTICAL AND HORIZONTAL')
        else:
            v = '|'
            h = '-'
            dl = dr = ul = ur = vl = uh = dh = vh = '+'

        if self.shape() == []:
            return dr + dl + '\n' + ur + ul

        # Get the widths of the columns
        str_tab = self._repr_tab()
        width = len(str_tab[0][0])
        str_list = [dr + (h*width + dh)*(len(str_tab[0])-1) + h*width + dl]
        for nrow, row in enumerate(str_tab):
            l1 = " " * (width+1) * nrow
            l2 = " " * (width+1) * nrow
            n = len(str_tab[nrow+1]) if nrow+1 < len(str_tab) else -1
            for i, e in enumerate(row):
                if i == 0:
                    l1 += ur + h*width
                elif i <= n+1:
                    l1 += vh + h*width
                else:
                    l1 += uh + h*width
                if unicode:
                    l2 += u"{}{:^{width}}".format(v, e, width=width)
                else:
                    l2 += "{}{:^{width}}".format(v, e, width=width)
            if i <= n:
                l1 += vl
            else:
                l1 += ul
            l2 += v
            str_list.append(l2)
            str_list.append(l1)
        return "\n".join(str_list)

    def pp(self):
        """
        Pretty print ``self``.

        EXAMPLES::

            sage: t = ShiftedPrimedTableau([[1,'2p',2,2],[2,'3p']])
            sage: t.pp()
            1  2' 2  2
               2  3'
            sage: t = ShiftedPrimedTableau([[10,'11p',11,11],[11,'12']])
            sage: t.pp()
            10  11' 11  11
                11  12
            sage: s = ShiftedPrimedTableau([['2p',2,3],['2p']],skew=[2,1])
            sage: s.pp()
             .  .  2' 2  3
                .  2'
        """
        print(self._repr_diagram())

    def _latex_(self):
        r"""
        Return LaTex code for ``self``.

        EXAMPLES::

            sage: T = ShiftedPrimedTableaux([4,2])
            sage: latex(T([[1,"2p",2,"3p"],[2,3]]))
            {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
            \raisebox{-.6ex}{$\begin{array}[b]{*{4}c}\cline{1-4}
            \lr{ 1 }&\lr{ 2'}&\lr{ 2 }&\lr{ 3'}\\\cline{1-4}
            &\lr{ 2 }&\lr{ 3 }\\\cline{2-3}
            \end{array}$}
            }
        """
        from sage.combinat.output import tex_from_array
        L = [[None]*i + row for i, row in enumerate(self._repr_tab())]
        return tex_from_array(L)

    def max_entry(self):
        """
        Return the maximum entry in the primed tableaux ``self``, rounded up.

        EXAMPLES::

            sage: Tab = ShiftedPrimedTableau([(1,1,'2p','3p'),(2,2)])
            sage: Tab.max_entry()
            3
        """
        if self == []:
            return 0
        else:
            flat = [entry.unprime() for row in self
                    for entry in row if entry is not None]
            return max(flat)

    def shape(self):
        """
        Return the shape of the underlying partition of ``self`` in list
        format.

        EXAMPLES::

            sage: t = ShiftedPrimedTableau([[1,'2p',2,2],[2,'3p']])
            sage: t.shape()
            [4, 2]
            sage: s = ShiftedPrimedTableau([["2p",2,3],["2p"]],skew=[2,1])
            sage: s.shape()
            [5, 2] / [2, 1]
            """
        if self._skew is None:
            return Partition([len(row) for row in self])
        return SkewPartition(([len(row) for row in self], self._skew))

    def weight(self):
        r"""
        Return the weight of ``self``.

        The weight of a shifted primed tableau is defined to be the vector
        with `i`-th component equal to the number of entries i and i' in the
        tableau.

        EXAMPLES::

           sage: t = ShiftedPrimedTableau([[1,'2p',2,2],[2,'3p']])
           sage: t.weight()
           (1, 4, 1)
        """
        flat = [entry.unprime() for row in self for entry in row]
        if flat == []:
            max_ind = 0
        else:
            max_ind = max(flat)
        weight = tuple([flat.count(i+1) for i in range(max_ind)])
        return weight


class CrystalElementShiftedPrimedTableau(ShiftedPrimedTableau):
    """
    Class for elements of ``crystals.ShiftedPrimedTableau``.
    """
    def _to_matrix(self):
        """
        Return a 2-dimensional array representation of a shifted tableau.

        EXAMPLES::

            sage: t = ShiftedPrimedTableaux([4,2,1])([[1,'2p',2,2],[2,'3p'],[3]])
            sage: mat = t._to_matrix()
            sage: mat
            [[1, 2', 2, 2], [None, 2, 3', None], [None, None, 3, None]]
        """
        m = self.shape()[0]
        return [[None]*i + list(row) + [None]*(m-i-len(row))
                for i, row in enumerate(self)]

    def _reading_word_with_positions(self):
        """
        Return the reading word of ``self`` together with positions of the
        corresponding letters in ``self``.

        The reading word of a shifted primed tableau is constructed
        as follows:

        1. List all primed entries in the tableau, column by
           column, in decreasing order within each column, moving
           from the rightmost column to the left, and with all
           the primes removed (i.e. all entries are increased by
            half a unit).

        2. Then list all unprimed entries, row by row, in
           increasing order within each row, moving from the
           bottommost row to the top.

        EXAMPLES::

            sage: t = ShiftedPrimedTableaux([4,2])([[1,'2p',2,2],[2,'3p']])
            sage: t._reading_word_with_positions()
            [((1, 2), 3), ((0, 1), 2), ((1, 1), 2), ((0, 0), 1),
            ((0, 2), 2), ((0, 3), 2)]
        """
        mat = self._to_matrix()
        ndim, mdim = len(mat), len(mat[0])
        list_with_positions = []
        for j in reversed(range(mdim)):
            for i in range(ndim):
                x = mat[i][j]
                if x is not None and x.is_primed():
                    list_with_positions.append(((i, j), x.unprime()))
        for i in reversed(range(ndim)):
            for j in range(mdim):
                x = mat[i][j]
                if x is not None and x.is_unprimed():
                    list_with_positions.append(((i, j), x.unprime()))
        return list_with_positions

    def reading_word(self):
        """
        Return the reading word of ``self``.

        The reading word of a shifted primed tableau is constructed
        as follows:

        1. List all primed entries in the tableau, column by
           column, in decreasing order within each column, moving
           from the rightmost column to the left, and with all
           the primes removed (i.e. all entries are increased by
            half a unit).

        2. Then list all unprimed entries, row by row, in
           increasing order within each row, moving from the
           bottommost row to the top.

        EXAMPLES::

            sage: t = ShiftedPrimedTableaux([4,2])([[1,'2p',2,2],[2,'3p']])
            sage: t.reading_word()
            [3, 2, 2, 1, 2, 2]
        """
        if self._skew is not None:
            raise NotImplementedError('skew tableau must be empty')
        return [tup[1] for tup in self._reading_word_with_positions()]

    def f(self, ind):
        r"""
        Compute the action of the crystal operator `f_i` on a shifted primed
        tableau using cases from the paper [HPS2017]_.

        INPUT:

        - ``ind`` -- element in the index set of the crystal

        OUTPUT:

        Primed tableau or ``None``.

        EXAMPLES::

            sage: t = ShiftedPrimedTableaux([5,4,2])([[1,1,1,1,'3p'],[2,2,2,'3p'],[3,3]])
            sage: t.pp()
            1  1  1  1  3'
               2  2  2  3'
                  3  3
            sage: s = t.f(2)
            sage: print(s)
            None
            sage: t = ShiftedPrimedTableaux([5,4,2])([[1,1,1,'2p','3p'],[2,2,3,3],[3,4]])
            sage: t.pp()
            1  1  1  2' 3'
               2  2  3  3
                  3  4
            sage: s = t.f(2)
            sage: s.pp()
            1  1  1  2' 3'
               2  3' 3  3
                  3  4

        """
        T = self._to_matrix()

        read_word = self._reading_word_with_positions()
        read_word = [num
                     for num in read_word
                     if num[1] == ind or num[1] == ind+1]

        element_to_change = None
        count = 0

        for element in read_word:
            if element[1] == ind+1:
                count += 1
            elif count == 0:
                element_to_change = element
            else:
                count -= 1

        if element_to_change is None:
            return None

        (r, c), elt = element_to_change

        if T[r][c].is_primed():
            T = [[elt+.5 if elt is not None else elt
                  for elt in row] for row in T]
            T = map(list, zip(*T))
            r, c = c, r
        h, l = len(T), len(T[0])

        if (c+1 == l or T[r][c+1] is None or T[r][c+1] >= ind+1):
            (tp_r, tp_c) = (r, c)
            while True:
                if (tp_r+1 == h or T[tp_r+1][tp_c] is None or T[tp_r+1][tp_c] > ind+1):
                    break
                if tp_r <= tp_c and T[tp_r+1][tp_r+1] == ind+1:
                    tp_r += 1
                    tp_c = tp_r
                    break
                if (ind+.5 not in T[tp_r+1]):
                    break
                tp_r += 1
                tp_c = T[tp_r].index(ind+.5)

            if tp_r == r:
                T[r][c] = T[r][c]+1
            elif tp_r == tp_c:
                T[r][c] = T[r][c]+.5
            else:
                T[r][c] = T[r][c]+.5
                T[tp_r][tp_c] = T[tp_r][tp_c]+.5

        elif T[r][c+1] == ind+.5:
            T[r][c+1] = T[r][c+1]+.5
            T[r][c] = T[r][c]+.5

        if r > c:
            T = [[elt-.5 if elt is not None else elt
                  for elt in row] for row in T]
            T = map(list, zip(*T))

        T = [tuple(elt for elt in row if elt is not None) for row in T]
        return type(self)(self.parent(), T, check=False, preprocessed=True)

    def e(self, ind):
        r"""
        Compute the action of the crystal operator `e_i` on a shifted primed
        tableau using cases from the paper [HPS2017]_.

        INPUT:

        - ``ind`` -- an element in the index set of the crystal.

        OUTPUT:

        Primed tableau or ``None``.

        EXAMPLES::

            sage: t = ShiftedPrimedTableaux([5,4,2])([[1,1,1,'2p','3p'],
            ....: [2,'3p',3,3],[3,4]])
            sage: t.pp()
            1  1  1  2' 3'
               2  3' 3  3
                  3  4
            sage: s = t.e(2)
            sage: s.pp()
            1  1  1  2' 3'
               2  2  3  3
                  3  4
            sage: t == s.f(2)
            True

        """
        T = self._to_matrix()

        read_word = self._reading_word_with_positions()
        read_word = [num for num in read_word
                     if num[1] == ind or num[1] == ind+1]

        element_to_change = None
        count = 0

        for element in read_word[::-1]:
            if element[1] == ind:
                count += 1
            elif count == 0:
                element_to_change = element
            else:
                count -= 1

        if element_to_change is None:
            return None
        (r, c), elt = element_to_change

        if T[r][c].is_primed():
            T = [[elt+.5 if elt is not None else elt
                  for elt in row] for row in T]
            T = map(list, zip(*T))
            r, c = c, r

        if (c == 0 or T[r][c-1] is None or T[r][c-1] <= ind):
            (tp_r, tp_c) = (r, c)
            while True:
                if (tp_r == 0 or T[tp_r-1][tp_c] is None or T[tp_r-1][tp_c] < ind):
                    break
                if (ind+.5 not in T[tp_r-1]):
                    break
                tp_r -= 1
                tp_c = T[tp_r].index(ind+.5)

            if tp_r == r:
                T[r][c] = T[r][c]-1.0
            elif tp_r == tp_c:
                T[r][c] = T[r][c]-.5
            else:
                T[r][c] = T[r][c]-.5
                T[tp_r][tp_c] = T[tp_r][tp_c]-.5

        elif T[r][c-1] == ind+.5:
            T[r][c-1] = T[r][c-1]-.5
            T[r][c] = T[r][c]-.5
        if r > c:
            T = [[elt-.5 if elt is not None else elt
                  for elt in row] for row in T]
            T = map(list, zip(*T))

        T = [tuple(elt for elt in row if elt is not None) for row in T]
        return type(self)(self.parent(), T, check=False, preprocessed=True)

    def is_highest_weight(self, index_set=None):
        r"""
        Return whether ``self`` is a highest weight element of the crystal.

        An element is highest weight if it vanishes under all crystal operators
        `e_i`.

        EXAMPLES::

            sage: t = ShiftedPrimedTableaux([5,4,2])([(1.0, 1.0, 1.0, 1.0, 1.0),
            ....: (2.0, 2.0, 2.0, 2.5), (3.0, 3.0)])
            sage: t.is_highest_weight()
            True
            sage: s = ShiftedPrimedTableaux([5,4])([(1.0, 1.0, 1.0, 1.0, 1.0),
            ....: (2.0, 2.0, 2.5, 3.0)])
            sage: s.is_highest_weight(index_set = [1])
            True
        """
        read_w = self.reading_word()
        max_entry = max(read_w)
        count = {i: 0 for i in range(max_entry+1)}
        if index_set is None:
            index_set = self.parent().index_set()
        for l in read_w[::-1]:
            count[l] += 1
            if (l-1 in index_set) and (l > 1) and (count[l] > count[l-1]):
                return False
        return True

    def weight(self):
        r"""
        Return the weight of ``self``.

        The weight of a shifted primed tableau is defined to be the vector
        with `i`-th component equal to the number of entries i and i' in the
        tableau.

        EXAMPLES::

           sage: t = ShiftedPrimedTableau([[1,'2p',2,2],[2,'3p']])
           sage: t.weight()
           (1, 4, 1)
        """
        flat = [entry.unprime() for row in self for entry in row]
        if flat == []:
            max_ind = 0
        else:
            max_ind = max(flat)
        weight = tuple([flat.count(i+1) for i in range(max_ind)])
        return self.parent().weight_lattice_realization()(weight)


class ShiftedPrimedTableauEntry(Rational):
    """
    The class of entries in shifted primed tableaux.
    """
    def __init__(self, entry):
        """
        Normalize the entry to be a half-integer in ``QQ``.

        TEST::
            sage: ShiftedPrimedTableau([[1,"2p"]])[0][1]
            2'
            sage: ShiftedPrimedTableau([[1,"2'"]])[0][1]
            2'
            sage: ShiftedPrimedTableau([[1,1.5]])[0][1]
            2'
        """
        if isinstance(entry, str):
            if (entry[-1] == "'" or entry[-1] == "p") and entry[:-1].isdigit() is True:
                # Check if an element has "'" or "p" at the end
                entry = QQ(entry[:-1]) - .5
        try:
            if (QQ(entry)+.5 in ZZ) or (QQ(entry) in ZZ):
                # Check if an element is a half-integer
                entry = QQ(entry)
            else:
                raise ValueError("all numbers must be half-integers")
        except (TypeError, ValueError):
            raise ValueError("primed elements have wrong format")
        Rational.__init__(self, entry)

    def __repr__(self):
        """
        Represent ``self`` as primed or unprimed integer.

        TEST::
            sage: ShiftedPrimedTableau([[1,"2p"]])[0][1]
            2'
        """
        if self.is_primed():
            return Integer(self.unprime()).__repr__() + "'"
        else:
            return Integer(self.unprime()).__repr__()

    def __add__(self, other):
        """
        Sum of ``self`` with rational number ``other``.

        TEST::
            sage: a = ShiftedPrimedTableau([[1,"2p"]])[0][1]
            sage: b = a+.5
            sage: type(b)
            <class 'sage.combinat.shifted_primed_tableau.ShiftedPrimedTableauEntry'>
        """
        return ShiftedPrimedTableauEntry(Rational.__add__(self, float(other)))

    def __sub__(self, other):
        """
        Sum of ``self`` with rational number ``other``.

        TEST::
            sage: a = ShiftedPrimedTableau([[1,"2p"]])[0][1]
            sage: b = a-.5
            sage: type(b)
            <class 'sage.combinat.shifted_primed_tableau.ShiftedPrimedTableauEntry'>
        """
        return ShiftedPrimedTableauEntry(Rational.__sub__(self, float(other)))


    def is_unprimed(self):
        """
        Checks if ``self`` is an unprimed element.

        TEST::
            sage: a = ShiftedPrimedTableau([[1,"2p"]])[0][1]
            sage: a.is_unprimed()
            False
        """
        return (self in ZZ)

    def is_primed(self):
        """
        Checks if ``self`` is a primed element.

        TEST::
            sage: a = ShiftedPrimedTableau([[1,"2p"]])[0][1]
            sage: a.is_primed()
            True
        """
        return (self not in ZZ)

    def unprime(self):
        """
        Unprime ``self`` if it was a prime element.

        TEST::
            sage: a = ShiftedPrimedTableau([[1,"2p"]])[0][1]
            sage: a.unprime()
            2
        """
        if self.is_primed():
            return Integer(self + .5)
        else:
            return Integer(self)


class ShiftedPrimedTableaux(UniqueRepresentation, Parent):
    r"""
    Returns the combinatorial class of shifted primed tableaux subject
    ro the constraints given by the arguments.

    A primed tableau is a tableau of shifted shape on the alphabet
    `X' = \{1' < 1 < 2' < 2 < \cdots < n' < n\}` such that

    1. the entries are weakly increasing along rows and columns

    2. a row cannot have two repeated primed entries, and a column
       cannot have two repeated non-primed entries

    3. there are only non-primed entries along the main diagonal

    INPUT:

    Valid optional keywords:

    - ``shape`` -- the (outer skew) shape of tableaux

    - ``weight`` -- the weight of tableaux

    - ``max_entry`` -- the maximum entry of tableaux

    - ``skew`` -- the inner skew shape of tableaux

    The weight of a tableau is defined to be the vector with `i`-th
    component equal to the number of entries `i` and `i'` in the tableau.
    The sum of the coordinates in the weight vector must be equal to the
    number of entries in the partition.

    The ``shape`` and ``skew`` must be strictly decreasing partitions.

    EXAMPLES::

        sage: SPT = ShiftedPrimedTableaux(weight=(1,2,2), shape=[3,2]); SPT
        Shifted Primed Tableaux of weight (1, 2, 2) and shape [3, 2]
        sage: SPT.list()
        [[(1, 2, 2), (3, 3)],
         [(1, 2', 3'), (2, 3)],
         [(1, 2', 3'), (2, 3')],
         [(1, 2', 2), (3, 3)]]
        sage: SPT = ShiftedPrimedTableaux(weight=(1,2)); SPT
        Shifted Primed Tableaux of weight (1, 2)
        sage: list(SPT)
        [[(1, 2, 2)], [(1, 2', 2)], [(1, 2'), (2,)]]
        sage: SPT = ShiftedPrimedTableaux([3,2], max_entry = 2); SPT
        Shifted Primed Tableaux of shape [3, 2] and maximum entry 2
        sage: list(SPT)
        [[(1, 1, 1), (2, 2)], [(1, 1, 2'), (2, 2)]]

    TESTS::

        sage: [(1,'2p',2,2),(2,'3p')] in ShiftedPrimedTableaux()
        True
        sage: [(1,1),(2,2)] in ShiftedPrimedTableaux()
        False
        sage: [] in ShiftedPrimedTableaux()
        True

    .. SEEALSO::

        - :class:`ShiftedPrimedTableau`
    """
    Element = ShiftedPrimedTableau
    options = Tableaux.options

    @staticmethod
    def __classcall_private__(cls, shape=None, weight=None,
                              max_entry=None, skew=None):
        r"""
        Normalize and process input to return the correct parent and ensure a
        unique representation. See the documentation for
        :class:`ShiftedPrimedTableaux` for more information.

        TESTS::

            sage: ShiftedPrimedTableaux([])
            Shifted Primed Tableaux of shape []
            sage: ShiftedPrimedTableaux(3)
            Traceback (most recent call last):
            ...
            ValueError: invalid shape argument
            sage: ShiftedPrimedTableaux(weight=(2,2,2), shape=[3,2])
            Traceback (most recent call last):
            ...
            ValueError: weight and shape are incompatible
            sage: ShiftedPrimedTableaux([[1]])
            Traceback (most recent call last):
            ...
            ValueError: invalid shape argument
            sage: ShiftedPrimedTableaux(weight=(2,2,2), max_entry=2)
            Traceback (most recent call last):
            ...
            ValueError: maximum entry is incompatible with the weight
            sage: ShiftedPrimedTableaux(shape=[4,1],skew=[3,2])
            Traceback (most recent call last):
            ...
            ValueError: skew shape must be inside the given tableau shape
        """
        if skew is not None:
            try:
                skew = Partition(skew)
            except ValueError:
                raise ValueError('invalid skew argument')
            if not all(skew[i] > skew[i+1] for i in range(len(skew)-1)):
                raise ValueError('skew shape must be a strict partition')

        if weight is not None:
            weight = tuple(weight)

        if shape is not None:
            if isinstance(shape, SkewPartition):
                skew = shape.inner()
                shape = shape.outer()
            try:
                shape = Partition(shape)
            except (ValueError, TypeError):
                raise ValueError('invalid shape argument')

            if not all(shape[i] > shape[i+1] for i in range(len(shape)-1)):
                raise ValueError("shape {} is not a strict partition".format(shape))

            if (skew is not None and not all(skew[i] <= shape[i]
                                             for i in range(len(skew)))):
                raise ValueError('skew shape must be inside the given tableau shape')

        if weight is not None:
            while weight[-1] == 0:
                weight = weight[:-1]

        if max_entry is not None and weight is not None:
            if len(weight) > max_entry:
                raise ValueError("maximum entry is incompatible with the weight")

        if shape is None:
            if weight is None:
                if max_entry is not None:
                    raise ValueError("specify shape or weight argument")
                return ShiftedPrimedTableaux_all(skew=skew)
            else:
                return ShiftedPrimedTableaux_weight(weight, skew=skew)
        else:
            if weight is None:
                return ShiftedPrimedTableaux_shape(shape, max_entry=max_entry, skew=skew)

            if (skew is not None and sum(shape) - sum(skew) != sum(weight)
                    or skew is None and sum(shape) != sum(weight)):
                raise ValueError("weight and shape are incompatible")

            return ShiftedPrimedTableaux_weight_shape(weight, shape, skew=skew)

    def __init__(self, skew=None):
        """
        Initialization of the parent class with given skew shape.
        """
        self._skew = skew

    def _element_constructor_(self, T):
        """
        Construct an object from ``T`` as an element of shifted primed
        tableaux, if possible.

        INPUT:

        - ``T`` -- data which can be interpreted as a primed tableau

        OUTPUT:

        - the corresponding primed tableau object

        """
        try:
            return self.element_class(self, T, skew=self._skew)
        except ValueError:
            raise ValueError("{} is not an element of {}".format(T, self))

    def _contains_tableau_(self, T):
        """
        Check if ``self`` contains preprocessed tableau ``T``.
        """
        if not all(len(T[i]) > len(T[i+1]) for i in range(len(T)-1)):
            return False
        if self._skew is not None:
            skew = self._skew + [0]*(len(T)-len(self._skew))
        else:
            skew = [0]*len(T)
        for i, row in enumerate(T):
            if i > 0:
                if not all(val > T[i-1][j+1]
                           for j, val in enumerate(row)
                           if j+1 >= skew[i-1]
                           if int(val) == val):
                    return False
                if not all(val >= T[i-1][j+1]
                           for j, val in enumerate(row)
                           if j+1 >= skew[i-1]
                           if int(val) != val):
                    return False
            if not all(row[j] <= row[j+1]
                       for j in range(skew[i], len(row)-1)
                       if int(row[j]) == row[j]):
                return False
            if not all(row[j] < row[j+1]
                       for j in range(skew[i], len(row)-1)
                       if int(row[j]) != row[j]):
                return False
        if not all(int(row[0]) == row[0]
                   for i, row in enumerate(T)
                   if skew[i] == 0):
            return False
        return True


class ShiftedPrimedTableaux_all(ShiftedPrimedTableaux):
    """
    The class of all shifted primed tableaux.
    """
    def __init__(self, skew=None):
        """
        Initialize the class of all shifted tableaux.

        TESTS::

            sage: [[1,1.5],[2]] in ShiftedPrimedTableaux()
            True
            sage: [[1,1.5],[1.5]] in ShiftedPrimedTableaux()
            False
            sage: [[1,1],[1]] in ShiftedPrimedTableaux()
            False
            sage: [[1,1],[2,2]] in ShiftedPrimedTableaux()
            False
        """
        Parent.__init__(self, category=InfiniteEnumeratedSets())
        ShiftedPrimedTableaux.__init__(self, skew=skew)
        self._skew = skew

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TESTS::

            sage: ShiftedPrimedTableaux()
            Shifted Primed Tableaux
        """
        if self._skew is None:
            return "Shifted Primed Tableaux"
        return "Shifted Primed Tableaux skewed by {}".format(self._skew)

    def _element_constructor_(self, T):
        """
        Construct an object from ``T`` as an element of ``self``, if
        possible.

        INPUT:

        - ``T`` -- data which can be interpreted as a tableau

        OUTPUT:

        - the corresponding tableau object

        TESTS::

            sage: Tab=ShiftedPrimedTableaux()([[1,1,"2p"]]); Tab
            [(1, 1, 2')]
            sage: Tab.parent()
            Shifted Primed Tableaux
            sage: Tab=ShiftedPrimedTableaux()([[1,1,2],[2,2]])
            Traceback (most recent call last):
            ...
            ValueError: [[1, 1, 2], [2, 2]] is not an element of
            Shifted Primed Tableaux

        """
        return (super(ShiftedPrimedTableaux_all, self)._element_constructor_(T))

    def __iter__(self):
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: Tabs = ShiftedPrimedTableaux()
            sage: Tabs[:5]
            [[], [(1,)], [(2,)], [(1, 2)], [(1, 2')]]
        """
        if self._skew is not None:
            raise NotImplementedError('skew tableau must be empty')
        yield self.element_class(self, [])

        max_entry = 1
        while(True):
            for size in range(1, max_entry+1):
                for shape in Partitions(size, max_slope=-1):
                    for weight in OrderedPartitions(size+max_entry-1,
                                                    k=max_entry):
                        weight = [weight[i]-1 for i in range(max_entry)]
                        weight[-1] += 1
                        for tab in ShiftedPrimedTableaux(shape=shape,
                                                         weight=weight):
                            yield self.element_class(self, tab, check=False,
                                                     preprocessed=True)
            max_entry += 1


class ShiftedPrimedTableaux_shape(ShiftedPrimedTableaux):
    r"""
    Shifted primed tableaux of a fixed shape.

    Shifted primed tableaux admit a type `A_n` classical crystal structure
    with highest weights corresponding to a given shape.

    The list of module generators consists of all elements of the
    crystal with nonincreasing weight entries.

    The crystal is constructed following operations described in [HPS2017]_.

    EXAMPLES::

        sage: ShiftedPrimedTableaux([4,3,1], max_entry=4)
        Shifted Primed Tableaux of shape [4, 3, 1] and maximum entry 4
        sage: ShiftedPrimedTableaux([4,3,1], max_entry=4).cardinality()
        384

    We compute some of the crystal structure::

        sage: SPTC = crystals.ShiftedPrimedTableaux([3,2], 3)
        sage: T = SPTC.module_generators[-1]
        sage: T
        [(1, 1, 2'), (2, 3')]
        sage: T.f(2)
        [(1, 1, 3'), (2, 3')]
        sage: len(SPTC.module_generators)
        7
        sage: SPTC[0]
        [(1, 1, 1), (2, 2)]
        sage: SPTC.cardinality()
        24
    """
    @staticmethod
    def __classcall_private__(cls, shape, max_entry=None, skew=None):
        """
        Normalize the attributes for the class.

        TEST::

            sage: SPT = ShiftedPrimedTableaux(shape=[2,1])
            sage: SPT._shape.parent()
            Partitions
        """
        shape = _Partitions(shape)
        return (super(ShiftedPrimedTableaux_shape, cls)
                .__classcall__(cls, shape=shape, max_entry=max_entry, skew=skew))

    def __init__(self, shape, max_entry, skew):
        """
        Initialize the class of shifted primed tableaux of a given shape.

        If ``max_elt`` is specified, a finite set with entries smaller
        or equal to ``max_elt``.

        TESTS::

            sage: TestSuite(ShiftedPrimedTableaux([4,2,1], max_entry=4)).run()
        """
        # Determine the correct category
        if max_entry is None:
            if skew is None:
                category = RegularCrystals().Infinite()
                self._cartan_type = CartanType(['A+oo'])
                self.Element = CrystalElementShiftedPrimedTableau
            else:
                category = Sets().Infinite()
        else:
            if skew is None:
                category = ClassicalCrystals()
                self._cartan_type = CartanType(['A', max_entry-1])
                self.Element = CrystalElementShiftedPrimedTableau
            else:
                category = Sets().Finite()

        ShiftedPrimedTableaux.__init__(self, skew=skew)
        Parent.__init__(self, category=category)
        self._max_entry = max_entry
        self._skew = skew
        if skew is None:
            self._shape = Partition(shape)
        else:
            self._shape = SkewPartition((shape, skew))

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TESTS::

            sage: ShiftedPrimedTableaux([3,2,1])
            Shifted Primed Tableaux of shape [3, 2, 1]
        """
        base = "Shifted Primed Tableaux of shape " + self._shape._repr_()
        if self._max_entry is not None:
            base += " and maximum entry {}".format(self._max_entry)
        return base

    def _contains_tableau_(self, T):
        """
        TESTS::

           sage: [[1,'2p',2,2],[2,'3p']] in ShiftedPrimedTableaux([4,2],
           ....: max_entry=4)
           True
           sage: [[1,'2p',2],[2,'3p']] in ShiftedPrimedTableaux([4,2])
           False
        """
        if not super(ShiftedPrimedTableaux_shape, self)._contains_tableau_(T):
            return False

        shape = [len(row) for row in T]
        skew = [row.count(None) for row in T]
        if sum(skew) == 0:
            shape = Partition(shape)
        else:
            shape = SkewPartition((shape, skew))
        if self._shape != shape:
            return False

        if self._max_entry is not None:
            flat = [round(item) for sublist in T for item in sublist]
            if flat == []:
                max_entry = 0
            else:
                max_entry = int(max(flat))
            if max_entry > self._max_entry:
                return False

        return True

    def _element_constructor_(self, T):
        """
        Construct an object from ``T`` as an element of ``self``, if
        possible.

        INPUT:

        - ``T`` -- data which can be interpreted as a primed tableau

        OUTPUT:

        - the corresponding primed tableau object

        TESTS::

            sage: tab= ShiftedPrimedTableaux([3])([[1,1,1.5]]); tab
            [(1, 1, 2')]
            sage: tab.parent()
            Shifted Primed Tableaux of shape [3]
            sage: ShiftedPrimedTableaux([3])([[1,1]])
            Traceback (most recent call last):
            ...
            ValueError: [[1, 1]] is not an element of Shifted Primed Tableaux
             of shape [3]
        """
        try:
            return (super(ShiftedPrimedTableaux_shape, self)
                    ._element_constructor_(T))
        except ValueError:
            raise ValueError("{} is not an element of {}".format(T, self))

    @lazy_attribute
    def module_generators(self):
        """
        Return the generators of ``self`` as a crystal.

        TEST::

            sage: SPT = ShiftedPrimedTableaux(shape=[2,1])
            sage: SPT.module_generators
            ([(1, 1), (2,)], [(1, 2), (3,)], [(1, 2'), (3,)])
        """
        if self._skew is not None:
            raise NotImplementedError("only for non-skew shapes")
        list_dw = []
        if self._max_entry is None:
            max_entry = sum(self._shape)
        else:
            max_entry = self._max_entry
        for weight in (Partition(self._shape).dominated_partitions(rows=max_entry)):
            list_dw.extend([self.element_class(self, T, check=False,
                                               preprocessed=True)
                            for T in ShiftedPrimedTableaux(weight=tuple(weight),
                                                           shape=self._shape)])
        return tuple(list_dw)

    def shape(self):
        """
        Return the shape of the shifted tableaux ``self``.

        TESTS::

            sage: ShiftedPrimedTableaux([6,4,3,1]).shape()
            [6, 4, 3, 1]
        """
        return self._shape

    def __iter__(self):
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: Tabs = ShiftedPrimedTableaux([3,2], max_entry=3)
            sage: Tabs[:3]
            [[(1, 1, 1), (2, 2)], [(1, 1, 1), (2, 3)], [(1, 1, 1), (2, 3')]]
            sage: len(list(Tabs))
            24
            sage: Tabs = ShiftedPrimedTableaux([3,2])
            sage: Tabs[:3]
            [[(1, 1, 1), (2, 2)], [(1, 1, 2'), (2, 2)], [(1, 1, 1), (2, 2)]]
        """
        if self._skew is not None:
            raise NotImplementedError('skew tableau must be empty')
        if self._max_entry is not None:
            for weight in OrderedPartitions(sum(self._shape)+self._max_entry,
                                            k=self._max_entry):
                weight_n = tuple([w-1 for w in weight])
                for tab in ShiftedPrimedTableaux(shape=self._shape, weight=weight_n):
                    yield self.element_class(self, tab, check=False,
                                             preprocessed=True)
        else:
            max_entry = 1
            while(True):
                for weight in OrderedPartitions(sum(self._shape)+max_entry, k=max_entry):
                    weight_n = tuple([w-1 for w in weight])
                    for tab in ShiftedPrimedTableaux(shape=self._shape, weight=weight_n):
                        yield self.element_class(self, tab, check=False,
                                                 preprocessed=True)
                max_entry += 1


class ShiftedPrimedTableaux_weight(ShiftedPrimedTableaux):
    """
    Shifted primed tableaux of fixed weight.

    TESTS::

        sage: ShiftedPrimedTableaux(weight=(2,3,1))
        Shifted Primed Tableaux of weight (2, 3, 1)
        sage: ShiftedPrimedTableaux(weight=(2,3,1)).cardinality()
        17
    """

    def __init__(self, weight, skew=None):
        """
        Initialize the class of shifted primed tableaux of a given weight.

        TESTS::

            sage: TestSuite( ShiftedPrimedTableaux(weight=(3,2,1)) ).run()
        """
        ShiftedPrimedTableaux.__init__(self, skew=skew)
        Parent.__init__(self, category=FiniteEnumeratedSets())
        self._weight = weight
        self._skew = skew

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TESTS::

            sage: ShiftedPrimedTableaux(weight=(3,2,1))
            Shifted Primed Tableaux of weight (3, 2, 1)
        """
        if self._skew is None:
            return "Shifted Primed Tableaux of weight {}".format(self._weight)
        return ("Shifted Primed Tableaux of weight {} skewed by {}"
                .format(self._weight, self._skew))

    def _element_constructor_(self, T):
        """
        Construct an object from ``T`` as an element of ``self``, if
        possible.

        INPUT:

        - ``T`` -- data which can be interpreted as a primed tableau

        OUTPUT:

        - the corresponding primed tableau object

        TESTS::

            sage: tab= ShiftedPrimedTableaux(weight=(2,1))([[1,1,1.5]]); tab
            [(1, 1, 2')]
            sage: tab.parent()
            Shifted Primed Tableaux of weight (2, 1)
        """
        try:
            return (super(ShiftedPrimedTableaux_weight, self)
                    ._element_constructor_(T))
        except ValueError:
            raise ValueError("{} is not an element of {}".format(T, self))

    def _contains_tableau_(self, T):
        """
        Check if ``self`` contains preprocessed tableau ``T``.

        TESTS::

            sage: [[1,1.5],[2]] in ShiftedPrimedTableaux(weight=(1,2))
            True
            sage: [[1,1.5],[3]] in ShiftedPrimedTableaux(weight=(1,2))
            False
        """
        if not super(ShiftedPrimedTableaux_weight, self)._contains_tableau_(T):
            return False

        flat = [round(item) for sublist in T for item in sublist]
        if flat == []:
            max_ind = 0
        else:
            max_ind = int(max(flat))
        weight = tuple([flat.count(i+1) for i in range(max_ind)])
        return self._weight == weight

    def __iter__(self):
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: Tabs = ShiftedPrimedTableaux(weight=(2,3))
            sage: Tabs[:4]
            [[(1, 1, 2, 2, 2)],
             [(1, 1, 2', 2, 2)],
             [(1, 1, 2, 2), (2,)],
             [(1, 1, 2', 2), (2,)]]
            sage: len(list(Tabs))
            5
        """
        for shape_ in Partitions(sum(self._weight)):
            if all(shape_[i] > shape_[i+1] for i in range(len(shape_)-1)):
                for tab in ShiftedPrimedTableaux(shape=shape_, weight=self._weight,
                                                 skew=self._skew):
                    yield self.element_class(self, tab, check=False,
                                             preprocessed=True)


class ShiftedPrimedTableaux_weight_shape(ShiftedPrimedTableaux):
    """
    Shifted primed tableaux of the fixed weight and shape.

    TESTS::

        sage: ShiftedPrimedTableaux([4,2,1], weight=(2,3,2))
        Shifted Primed Tableaux of weight (2, 3, 2) and shape [4, 2, 1]
        sage: ShiftedPrimedTableaux([4,2,1], weight=(2,3,2)).cardinality()
        4
    """
    def __init__(self, weight, shape, skew=None):
        """
        Initialize the class of shifted primed tableaux of the given weight
        and shape.

        TESTS::

            sage: TestSuite( ShiftedPrimedTableaux([4,2,1], weight=(3,2,2)) ).run()
        """
        ShiftedPrimedTableaux.__init__(self, skew=skew)
        Parent.__init__(self, category=FiniteEnumeratedSets())
        self._weight = weight
        self._skew = skew
        if skew is None:
            self._shape = Partition(shape)
        else:
            self._shape = SkewPartition((shape, skew))

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TESTS::

            sage: ShiftedPrimedTableaux([3,2,1], weight=(4,2))
            Shifted Primed Tableaux of weight (4, 2) and shape [3, 2, 1]
        """
        return ("Shifted Primed Tableaux of weight {} and shape {}"
                .format(self._weight, self._shape))

    def _contains_tableau_(self, T):
        """
        Check if ``self`` contains preprocessed tableau ``T``.

        TESTS::

            sage: [[1,1.5],[2]] in ShiftedPrimedTableaux([2,1], weight=(1,2))
            True
            sage: [[1,1.5],[3]] in ShiftedPrimedTableaux([2,1], weight=(1,2))
            False
            sage: [[1,1.5,2,3],[3]] in ShiftedPrimedTableaux([3,2], weight=(1,2,2))
            False
        """

        if not super(ShiftedPrimedTableaux_weight_shape, self)._contains_tableau_(T):
            return False

        flat = [round(item) for sublist in T for item in sublist]
        if flat == []:
            max_ind = 0
        else:
            max_ind = int(max(flat))
        weight = tuple([flat.count(i+1) for i in range(max_ind)])
        if self._weight!=weight:
            return False

        shape = [len(row) for row in T]
        skew = [row.count(None) for row in T]
        if sum(skew) == 0:
            shape = Partition(shape)
        else:
            shape = SkewPartition((shape, skew))
        if self._shape != shape:
            return False

        return True

    def _element_constructor_(self, T):
        """
        Construct an object from ``T`` as an element of ``self``, if
        possible.

        TESTS::

            sage: tab = ShiftedPrimedTableaux([3], weight=(2,1))([[1,1,1.5]]); tab
            [(1, 1, 2')]
            sage: tab.parent()
            Shifted Primed Tableaux of weight (2, 1) and shape [3]
            sage: ShiftedPrimedTableaux([3], weight=(2,1))([[1,1]])
            Traceback (most recent call last):
            ...
            ValueError: [[1, 1]] is not an element of Shifted Primed Tableaux
            of weight (2, 1) and shape [3]
        """
        try:
            return (super(ShiftedPrimedTableaux_weight_shape, self)
                    ._element_constructor_(T))
        except ValueError:
            raise ValueError("{} is not an element of {}".format(T, self))

    def __iter__(self):
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: Tabs = ShiftedPrimedTableaux([3,2], weight=(1,2,2))
            sage: Tabs[:4]
            [[(1, 2, 2), (3, 3)],
             [(1, 2', 3'), (2, 3)],
             [(1, 2', 3'), (2, 3')],
             [(1, 2', 2), (3, 3)]]
            sage: len(list(Tabs))
            4

        TEST::

            sage: Tabs = ShiftedPrimedTableaux([3,2], weight=(1,4))
            sage: list(Tabs)
            []
        """
        if self._skew is not None:
            raise NotImplementedError('skew tableau must be empty')

        if not self._shape.dominates(
              Partition(sorted(list(self._weight), reverse=True))):
            return
        full_shape = self._shape
        sub_tab = []
        tab_list_new = [[]]
        for i, w in enumerate(self._weight):
            tab_list_old = tab_list_new
            tab_list_new = []
            for sub_tab in tab_list_old:
                sub_shape = [len(sub_tab[r]) for r in range(len(sub_tab))]
                for strip in _add_strip(sub_shape, full_shape, w):
                    l = int(len(strip)/2)
                    if len(sub_shape) < len(full_shape):
                        new_tab = [sub_tab[r]
                                   + [float(i+.5)]*int(strip[r])
                                   + [float(i+1)]*int(strip[-r-1])
                                   for r in range(l-1)]
                        if strip[l] != 0:
                            new_tab.append([float(i+1)]*int(strip[l]))
                    else:
                        new_tab = [sub_tab[r]
                                   + [float(i+.5)]*int(strip[r])
                                   + [float(i+1)]*int(strip[-r-1])
                                   for r in range(l)]
                    tab_list_new.append(new_tab)
        for tab in tab_list_new:
            yield self.element_class(self, tab)


####################
# Helper functions #
####################


def _add_strip(sub_tab, full_tab, length):
    """
    Helper function used in the algorithm to generate all shifted primed
    tableaux of the fixed weight and shape.

    TESTS::

        sage: list(ShiftedPrimedTableaux([3,1], weight=(2,2)))  # indirect doctest
        [[(1, 1, 2), (2,)], [(1, 1, 2'), (2,)]]
    """
    if sum(sub_tab)+length > sum(full_tab):
        raise ValueError("strip does not fit")

    if len(sub_tab) == 0:
        cliff_list = []
    else:
        cliff_list = [int(sub_tab[0] != full_tab[0])]

    for row in range(1, len(sub_tab)):
        if sub_tab[row] == full_tab[row]:
            cliff_list.append(0)
        elif sub_tab[row-1]-1 == sub_tab[row]:
            cliff_list[-1] += 1
        else:
            cliff_list.append(1)

    if len(sub_tab) < len(full_tab):
        cliff_list.append(0)

    for primes_num in range(min(sum(cliff_list), length)+1):
        for primed_list in IntegerVectors(n=primes_num, k=len(cliff_list),
                                          outer=cliff_list):
            row = 0
            primed_strip = list()
            for i, cliff in enumerate(cliff_list):
                if cliff == 0:
                    row += 1
                    primed_strip.append(0)
                    pass
                primed_strip.extend([int(primed_list[i] > j)
                                     for j in range(cliff)])
                row += cliff
            plat_list = list()

            if len(sub_tab) < len(full_tab) and len(sub_tab) != 0:
                plat_list.append(min(sub_tab[-1] + primed_strip[-2] - 1,
                                     full_tab[len(sub_tab)]))
            for row in range(1, len(sub_tab))[::-1]:
                plat_list.append(
                    min(sub_tab[row-1]+primed_strip[row-1]-1, full_tab[row])
                    - sub_tab[row] - primed_strip[row])
            if len(sub_tab) > 0:
                plat_list.append(full_tab[0] - sub_tab[0] - primed_strip[0])
            else:
                plat_list.append(full_tab[0])

            if sum(plat_list) < length - primes_num:
                pass
            for non_primed_strip in IntegerVectors(n=length-primes_num,
                                                   k=len(plat_list),
                                                   outer=plat_list):
                yield (list(primed_strip) + list(non_primed_strip))

