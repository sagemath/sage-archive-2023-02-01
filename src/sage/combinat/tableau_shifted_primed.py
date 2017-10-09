from six import add_metaclass
import numpy as np

from sage.combinat.partition import Partition, Partitions, OrderedPartitions
from sage.combinat.integer_vector import IntegerVectors
from sage.rings.integer import Integer

from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass

from sage.structure.list_clone import ClonableArray
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation

# Imports for the crystal

from sage.categories.classical_crystals import ClassicalCrystals
from sage.combinat.root_system.cartan_type import CartanType


@add_metaclass(InheritComparisonClasscallMetaclass)
class ShiftedPrimedTableau(ClonableArray):
    r"""
    A shifted primed tableau with primed elements stored as half-integers.

    A primed tableau is a shifted tableau on the alphabet
    X' = {1' < 1 < 2' < 2 <...< n' < n} such that

        1. the entries are weakly increasing along rows and columns
        2. a row can't have two repeated primed elements, and a column
           can't have two repeated non-primed elements
        3. there are only non-primed elements on the main diagonal

    EXAMPLES::

        sage: T = ShiftedPrimedTableaux([4,2])
        sage: T([[1,"2'","3'",3],[2,"3'"]])[1]
        (2.0, 2.5)
        sage: t = ShiftedPrimedTableau([[1,"2p",2.5,3],[0,2,2.5]])
        sage: t[1]
        (2.0, 2.5)

    TEST::

        sage: t = ShiftedPrimedTableau([[1,2,2.5,3],[0,2,2.5]])
        Traceback (most recent call last):
        ...
        ValueError: [[1, 2, 2.50000000000000, 3], [0, 2, 2.50000000000000]]
        is not an element of Shifted Primed Tableaux
    """

    @staticmethod
    def __classcall_private__(cls, T):
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

        """

        if isinstance(T, cls):
            return T

        return ShiftedPrimedTableaux()(T)

    def __init__(self, parent, T):
        r"""
        Initialize a shifted tableau.

        TESTS::

            sage: s = ShiftedPrimedTableau([[1,"2'","3'",3],[2,"3'"]])
            sage: t = ShiftedPrimedTableaux([4,2])([[1,"2p","3p",3],
            ....: [0, 2,"3p"]])
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

        if isinstance(T, ShiftedPrimedTableau):
            ClonableArray.__init__(self, parent, T)
            return

        if not isinstance(T, list):
            try:
                t_ = list(T)
            except TypeError:
                t_ = [T]
        else:
            t_ = T

        if not all(isinstance(row, (list, tuple, np.ndarray)) for row in t_):
            t_ = [t_]

        t = []
        # Preprocessing list t for primes and other symbols
        for row in t_:
            row_new = []
            for element in row:

                if isinstance(element, str):
                    if element[-1] == "'" and element[:-1].isdigit() is True:
                        # Check if an element has "'" at the end
                        row_new.append(float(float(element[:-1]) - .5))
                        continue
                    if element[-1] == "p" and element[:-1].isdigit() is True:
                        # Check if an element has "p" at the end
                        row_new.append(float(float(element[:-1]) - .5))
                        continue
                try:
                    if int(float(element)*2) == float(element)*2:
                        # Check if an element is a half-integer
                        row_new.append(float(element))
                        continue
                    else:
                        raise ValueError("all numbers must be half-integers")

                except (TypeError, ValueError):
                    raise ValueError("primed elements have wrong format")

            t.append(row_new)

        # Accounting for zeros at the beginning and at the end of a row'
        i = 0
        while i < len(t):
            row = t[i]
            try:
                while row[0] == 0:
                    row.pop(0)
                while row[-1] == 0:
                    row.pop(-1)
            except IndexError:
                t.remove(t[i])
                continue
            t[i] = row
            i += 1

        t = [tuple(_) for _ in t]
        ClonableArray.__init__(self, parent, t)

    def __eq__(self, other):
        """
        Check whether ``self`` is equal to ``other``.

        INPUT:

        - ``other`` -- the element that ``self`` is compared to

        OUTPUT:

        Boolean.

        EXAMPLE::

            sage: t = ShiftedPrimedTableau([[1,"2p"]])
            sage: t == ShiftedPrimedTableaux([2])([1,1.5])
            True
        """
        try:
            Tab = ShiftedPrimedTableau(other)
        except ValueError:
            return False
        return (list(self) == list(Tab))

    def _to_matrix(self):
        """
        Return a 2-dimensional numpy.array representation of a shifted tableau

        EXAMPLE::

            sage: t = ShiftedPrimedTableau([[1,'2p',2,2],[2,'3p'],[3]])
            sage: mat = t._to_matrix()
            sage: mat
            array([[ 1. ,  1.5,  2. ,  2. ],
                   [ 0. ,  2. ,  2.5,  0. ],
                   [ 0. ,  0. ,  3. ,  0. ]])
            sage: t == ShiftedPrimedTableau(mat)
            True
        """
        array = []
        m = len(self[0])
        for i in range(len(self)):
            array.append([0]*i + list(self[i]) + [0]*(m-i-len(self[i])))
        array = np.array(array, dtype='float')
        return array

    def check(self):
        """
        Check that ``self`` is a valid primed tableaux.

        EXAMPLE::

            sage: T = ShiftedPrimedTableaux([4,2])
            sage: t = T([[1,'2p',2,2],[2,'3p']])
            sage: t.check()
        """
        if not all(len(self[i]) > len(self[i+1]) for i in range(len(self)-1)):
            raise ValueError('shape must be a strict partition')
        for i, row in enumerate(self):
            if i > 0:
                if not all(val > self[i-1][j+1]
                           for j, val in enumerate(row) if int(val) == val):
                    raise ValueError(
                        'column is not strictly increasing in non-primes')
                if not all(val >= self[i-1][j+1]
                           for j, val in enumerate(row) if int(val) != val):
                    raise ValueError(
                        'column is not weakly increasing in primes')
            if not all(row[j] <= row[j+1]
                       for j in range(len(row)-1) if int(row[j]) == row[j]):
                raise ValueError('row is not weakly increasing in non-primes')
            if not all(row[j] < row[j+1]
                       for j in range(len(row)-1) if int(row[j]) != row[j]):
                raise ValueError('row is not strictly increasing in primes')
        if not all(int(row[0]) == row[0] for row in self):
            raise ValueError('diagonal elements must be non-primes')

    def _repr_(self):
        """
        Represent Shifted Primed Tableau as a list of rows,
        rows are represented as tuples of half-integers.

        EXAMPLE::

            sage: t = ShiftedPrimedTableau([[1,'2p',2,2],[2,'3p']])
            sage: t
            [(1.0, 1.5, 2.0, 2.0), (2.0, 2.5)]
        """
        return repr([tuple(_) for _ in self])

    def _repr_tab(self):
        """
        Return a nested list of strings representing the elements.

        TEST::

            sage: t = ShiftedPrimedTableau([[1,'2p',2,2],[2,'3p']])
            sage: t._repr_tab()
            [[' 1 ', " 2'", ' 2 ', ' 2 '], [' 2 ', " 3'"]]
        """
        max_len = len(str(self.max_element()))+1
        string_tab = []
        for i, row in enumerate(self):
            string_row = []
            for val in row:
                if int(val) == val:
                    string_row.append(' '*(max_len-len(str(int(val))))
                                      + str(int(val)) + ' ')
                else:
                    string_row.append(' '*(max_len-len(str(int(val+.5))))
                                      + str(int(val+.5)) + "'")
            string_tab.append(string_row)
        return string_tab

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

        """
        max_len = len(str(self.max_element()))+2
        string_list = ""
        for i, row in enumerate(self):
            string_list += ' '
            string_list += ' ' * i * (max_len)
            for val in row:
                if int(val) == val:
                    string_list += (str(int(val))
                                    + ' ' * (max_len-len(str(int(val)))))
                else:
                    string_list += (str(int(val+.5)) + "'"
                                    + ' '*(max_len-1 - len(str(int(val+.5)))))
            string_list += '\n'
        string_list = string_list[:-2]
        return string_list

    def _ascii_art_(self):
        """
        Return ASCII representaion of a tableau.

        EXAMPLE::

            sage: ascii_art(ShiftedPrimedTableau([[1,'2p',2,2],[2,'3p']]))
            +---+---+---+---+
            | 1 | 2'| 2 | 2 |
            +---+---+---+---+
                | 2 | 3'|
                +---+---+

        TEST::

            sage: ascii_art(ShiftedPrimedTableau([]))
            ++
            ++
        """
        from sage.typeset.ascii_art import AsciiArt
        return AsciiArt(self._ascii_art_table(unicode=False).splitlines())

    def _unicode_art_(self):
        """
        Return a Unicode representation of a tableau.

        EXAMPLE::

            sage: unicode_art(ShiftedPrimedTableau([[1,'2p',2,2],[2,'3p']]))
            ┌───┬───┬───┬───┐
            │ 1 │ 2'│ 2 │ 2 │
            └───┼───┼───┼───┘
                │ 2 │ 3'│
                └───┴───┘

        TEST::
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
        Print out a nice version of ``self``.

        EXAMPLES::

            sage: t = ShiftedPrimedTableau([[1,'2p',2,2],[2,'3p']])
            sage: t.pp()
            1  2' 2  2
               2  3'
            sage: t = ShiftedPrimedTableau([[10,'11p',11,11],[11,'12']])
            sage: t.pp()
            10  11' 11  11
                11  12
        """
        print(self._repr_diagram())

    def _latex_(self):
        r"""
        Return LaTex code for ``self``.

        EXAMPLE::

            sage: T = ShiftedPrimedTableaux([4,2])
            sage: latex(T([[1,"2p",2,"3p"],[2,3]]))
            {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
            \raisebox{-.6ex}{$\begin{array}[b]{*{4}c}\cline{1-4}
            \lr{1}&\lr{2'}&\lr{2}&\lr{3'}\\\cline{1-4}
            &\lr{2}&\lr{3}\\\cline{2-3}
            \end{array}$}
            }
        """
        from sage.combinat.output import tex_from_array
        L = list()
        for i, row in enumerate(self):
            num_list = [None]*i
            for let in row:
                if int(let) == let:
                    num_list.append(int(let))
                else:
                    num_list.append(str(int(let+0.5))+"'")
            L.append(num_list)
        return tex_from_array(L)

    def max_element(self):
        """
        Return the maximum element in the primed tableaux ``self``, rounded up.

        EXAMPLE::

            sage: Tab = ShiftedPrimedTableau([(1,1,1.5,2.5),(2,2)])
            sage: Tab.max_element()
            3
        """
        if self == []:
            return 0
        else:
            flat = [item for sublist in self for item in sublist]
            return int(round(max(flat)))

    def shape(self):
        """
        Return the shape of the underlying partition of ``self`` in list
        format.

        EXAMPLES::

            sage: t = ShiftedPrimedTableau([[1,'2p',2,2],[2,'3p']])
            sage: t.shape()
            [4, 2]
        """
        return ([len(row) for row in self])

    def __call__(self, *cell):
        """
        Function call of ``self``.

        INPUT:

        - ``cell`` -- a pair of integers, tuple, or list specifying a cell in
          the tableau.

        OUTPUT:

        - the element in the corresponding cell. if the element is primed,
          returns half-integer value.

        EXAMPLES::

            sage: t  = ShiftedPrimedTableau([[1,'2p',2,2],[2,'3p']])
            sage: t(1,0)
            2.0
            sage: t((1,2))
            Traceback (most recent call last):
            ...
            IndexError: invalid cell
            sage: t((1,1))
            2.5
        """

        try:
            i, j = cell
        except ValueError:
            i, j = cell[0]

        try:
            return self[i][j]
        except IndexError:
            raise IndexError("invalid cell")

    def weight(self):
        """
        Return the weight of ``self`` as a tuple.

        The weight of a shifted primed tableau is defined to be the vector
        with i-th component equal to the number of letters i and i' in the
        tableau.

        EXAMPLE::

           sage: t = ShiftedPrimedTableau([[1,'2p',2,2],[2,'3p']])
           sage: t.weight()
           (1, 4, 1)
        """
        flat = [round(item) for sublist in self for item in sublist]
        if flat == []:
            max_ind = 0
        else:
            max_ind = int(max(flat))
        weight = tuple([flat.count(i+1) for i in range(max_ind)])
        return tuple(weight)

    def _reading_word_with_positions(self):
        """
        Return the reading word of ``self`` together with positions of the
        corresponding letters in ``self``.

        The reading word of a shifted primed tableau is constructed as follows:

            1. List all primed letters in the tableau, column by
               column, in decreasing order within each column, moving
               from the rightmost column to the left, and with all
               the primes removed (i.e. all letters are increased by
               half a unit).

            2. Then list all unprimed elements, row by row, in
               increasing order within each row, moving from the
               bottommost row to the top.

        EXAMPLE::

            sage: t = ShiftedPrimedTableau([[1,'2p',2,2],[2,'3p']])
            sage: t._reading_word_with_positions()
            [((1, 2), 3), ((0, 1), 2), ((1, 1), 2), ((0, 0), 1),
            ((0, 2), 2), ((0, 3), 2)]

        """
        mat = self._to_matrix()
        list_with_positions = []
        for (i, j), x in np.ndenumerate(mat[:, ::-1].T):
            if int(x) != x:
                list_with_positions.append(((j, mat.shape[1]-i-1), int(x+0.5)))
        for (i, j), x in np.ndenumerate(mat[::-1, :]):
            if int(x) == x and int(x) != 0:
                list_with_positions.append(((mat.shape[0]-i-1, j), int(x)))
        return list_with_positions

    def reading_word(self):
        """
        Return the reading word of ``self``.

        The reading word of a shifted primed tableau is constructed as follows:

            1. List all primed letters in the tableau, column by
               column, in decreasing order within each column, moving
               from the rightmost column to the left, and with all
               the primes removed (i.e. all letters are increased by
               half a unit).

            2. Then list all unprimed elements, row by row, in
               increasing order within each row, moving from the
               bottommost row to the top.

        EXAMPLE::

            sage: t = ShiftedPrimedTableau([[1,'2p',2,2],[2,'3p']])
            sage: t.reading_word()
            [3, 2, 2, 1, 2, 2]
        """
        return [tup[1] for tup in self._reading_word_with_positions()]

    def f(self, ind):
        """
        Compute the action of the crystal operator `f_i` on a Shifted Primed
        Tableau using cases from the paper [GPS.17].

        INPUT::

        - ``self`` -- shifted primed tableau
        - ``ind`` -- index of the crystal operator `f_i`.

        OUTPUT:

        Primed tableau or 'None'.

        EXAMPLES::

            sage: t = ShiftedPrimedTableau([[1,1,1,1,2.5],[2,2,2,2.5],[3,3]])
            sage: t.pp()
            1  1  1  1  3'
               2  2  2  3'
                  3  3
            sage: s = t.f(2)
            sage: print(s)
            None
            sage: t = ShiftedPrimedTableau([[1,1,1,1.5,2.5],[0,2,2,3,3],
            ....: [0,0,3,4]])
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
        if self is None:
            return None

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

        if T[r, c] == ind - .5:
            T = T.T + .5
            r, c = c, r
        h, l = T.shape

        if (c+1 == l or T[r, c+1] >= ind+1 or T[r, c+1] < 1):

            (tp_r, tp_c) = (r, c)
            while True:

                if (tp_r+1 == h or
                        T[tp_r+1, tp_c] > ind+1 or
                        T[tp_r+1, tp_c] < 1):
                    break
                if tp_r <= tp_c and T[tp_r+1, tp_r+1] == ind+1:
                    tp_r += 1
                    tp_c = tp_r
                    break
                if (ind+.5 not in T[tp_r+1]):
                    break
                tp_r += 1
                tp_c = np.where(T[tp_r] == ind+.5)[0]

            if tp_r == r:
                T[r, c] += 1

            elif tp_r == tp_c:
                T[r, c] += .5

            else:
                T[r, c] += .5
                T[tp_r, tp_c] += .5

        elif T[r, c+1] == ind+.5:
            T[r, c+1] += .5
            T[r, c] += .5

        if r > c:
            T = T.T - .5

        return(ShiftedPrimedTableau(T))

    def e(self, ind):
        """
        Compute the action of the crystal operator `e_i` on a Shifted Primed
        Tableau using cases from the paper [HPS.17].

        INPUT:

        - ``self`` -- shifted primed tableau
        - ``ind`` -- index of the crystal operator `e_i`.

        OUTPUT:

        Primed tableau or 'None'.

        EXAMPLES::

            sage: t = ShiftedPrimedTableau([[1,1,1,1.5,2.5],
            ....: [0,2,"3p",3,3],[0,0,3,4]])
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
        if self is None:
            return None
        T = self._to_matrix()

        read_word = self._reading_word_with_positions()
        read_word = [num
                     for num in read_word
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

        if T[r, c] == ind + .5:
            T = T.T + .5
            r, c = c, r
        h, l = T.shape

        if (c == 0 or T[r, c-1] <= ind or T[r, c-1] < 1):

            (tp_r, tp_c) = (r, c)
            while True:

                if tp_r == 0 or T[tp_r-1, tp_c] < ind or T[tp_r-1, tp_c] < 1:
                    break
                if (ind+.5 not in T[tp_r-1]):
                    break
                tp_r -= 1
                tp_c = np.where(T[tp_r] == ind+.5)[0]

            if tp_r == r:
                T[r, c] -= 1

            elif tp_r == tp_c:
                T[r, c] -= .5

            else:
                T[r, c] -= .5
                T[tp_r, tp_c] -= .5

        elif T[r, c-1] == ind+.5:
            T[r, c-1] -= .5
            T[r, c] -= .5

        if r > c:
            T = T.T - .5

        return(ShiftedPrimedTableau(T))

    def is_highest_weight(self):
        """
        Check wether the shifted primed tableau ``self`` is a highest weight
        element of the crystal.


        An element is highest weight if it vanishes under any crystal operator
        `e_i`.

        EXAMPLES::

            sage: t = ShiftedPrimedTableau([(1.0, 1.0, 1.0, 1.0, 1.0),
            ....: (2.0, 2.0, 2.0, 2.5), (3.0, 3.0)])
            sage: print(t.e(1))
            None
            sage: print(t.e(2))
            None
            sage: t.is_highest_weight()
            True

        """
        read_w = self.reading_word()
        count = {}
        for l in read_w[::-1]:
            try:
                count[l] += 1
            except KeyError:
                count[l] = 1
            if l > 1:
                try:
                    if count[l] > count[l-1]:
                        return False
                except KeyError:
                    return False
        return True


class ShiftedPrimedTableaux(UniqueRepresentation, Parent):
    """
    A factory for the various classes of shifted standard tableaux.

    INPUT:

    - a weight and/or a partition

    OUTPUT:

    - with no argument, the class of all primed tableaux

    - with a list or partition argument, the class of all primed
      tableaux of that shape (infinite set if we don't specify the
      weight or maximum element)

    - with an additional integer argument ``max_element``, the class
      of all primed tableaux of a given shape and a given maximum
      element

    - with a tuple argument, the class of all primed tableaux of that
      weight (finite set)

    A primed tableau is a shifted tableau on the alphabet
    X' = {1' < 1 < 2' < 2 <...< n' < n} such that

        1. the entries are weakly increasing along rows and columns
        2. a row can't have two repeated primed elements, and a column
           can't have two repeated non-primed elements
        3. there are only non-primed elements along the main diagonal

    The weight of a tableau is defined to be the vector with i-th
    component equal to the number of letters i and i' in the tableau.
    The sum of the entries in the weight vector must be equal to the
    number of boxes in the partition.

    None of the Shifted Primed Tableaux classes can be iterated over.

    EXAMPLES::

        sage: SPT = ShiftedPrimedTableaux(weight=(1,2,2), shape=[3,2]); SPT
        Shifted Primed Tableaux of weight (1, 2, 2) and shape [3, 2]
        sage: SPT.list()
        [[(1.0, 2.0, 2.0), (3.0, 3.0)],
         [(1.0, 1.5, 2.5), (2.0, 3.0)],
         [(1.0, 1.5, 2.5), (2.0, 2.5)],
         [(1.0, 1.5, 2.0), (3.0, 3.0)]]

        sage: SPT = ShiftedPrimedTableaux((1,2)); SPT
        Shifted Primed Tableaux of weight (1, 2)
        sage: list(SPT)
        [[(1.0, 2.0, 2.0)], [(1.0, 1.5, 2.0)], [(1.0, 1.5), (2.0,)]]

        sage: SPT = ShiftedPrimedTableaux([3,2], max_element = 2); SPT
        Shifted Primed Tableaux of shape [3, 2] and maximum element 2
        sage: list(SPT)
        [[(1.0, 1.0, 1.0), (2.0, 2.0)], [(1.0, 1.0, 1.5), (2.0, 2.0)]]
        sage: SPT.list_highest_weight()
        [[(1.0, 1.0, 1.0), (2.0, 2.0)]]

    .. SEEALSO::

        - :class:`ShiftedPrimedTableau`
    """

    Element = ShiftedPrimedTableau

    @staticmethod
    def __classcall_private__(cls, *args, **kwargs):
        r"""
        This is a factory class which returns the appropriate parent based on
        arguments.  See the documentation for :class:`ShiftedPrimedTableaux`
        for more information.


        TESTS::

            sage: ShiftedPrimedTableaux([])
            Shifted Primed Tableaux of shape []

            sage: ShiftedPrimedTableaux(3)
            Traceback (most recent call last):
            ...
            ValueError: invalid argument for weight or shape

            sage: ShiftedPrimedTableaux(weight=(2,2,2), shape=[3,2])
            Traceback (most recent call last):
            ...
            ValueError: weight and shape are incompatible

            sage: ShiftedPrimedTableaux([[1]])
            Traceback (most recent call last):
            ...
            ValueError: shape [[1]] is not a partition

            sage: ShiftedPrimedTableaux(weight=(2,2,2), max_element=2)
            Traceback (most recent call last):
            ...
            ValueError: maximum element is incompatible with the weight
        """
        weight = None
        shape = None
        max_element = None

        if 'max_elt' in kwargs:
            max_element = int(kwargs['max_elt'])

        if 'max_element' in kwargs:
            max_element = int(kwargs['max_element'])

        if 'max' in kwargs:
            max_element = int(kwargs['max'])

        if 'size' in kwargs and isinstance(kwargs['size'],
                                           (list, tuple, Partition)):
            shape = list(kwargs['size'])

        if 'shape' in kwargs:
            shape = list(kwargs['shape'])

        if 'weight' in kwargs:
            weight = tuple(kwargs['weight'])

        if args:
            if isinstance(args[0], tuple) and weight is None:
                weight = args[0]
                if len(args) > 1:
                    if ((isinstance(args[1], (list, Partition)) or
                         args[1] is None) and shape is None):
                        shape = args[1]
                    else:
                        raise ValueError(
                            'invalid argument for weight or shape')

            elif (isinstance(args[0], (list, Partition)) and shape is None):
                shape = args[0]
                if len(args) > 1:
                    if ((isinstance(args[1], tuple) or args[1] is None) and
                            weight is None):
                        weight = args[1]
                    else:
                        raise ValueError(
                            'invalid argument for weight or shape')
            else:
                raise ValueError('invalid argument for weight or shape')

        if shape is not None:
            try:
                shape = Partition(shape)
            except ValueError:
                raise ValueError('shape {} is not a partition'.format(shape))

        if weight is not None:
            while weight[-1] == 0:
                weight = weight[:-1]

        if max_element is not None and weight is not None:
            if len(weight) != max_element:
                raise ValueError(
                    "maximum element is incompatible with the weight")

        if max_element is not None and shape is not None:
            if max_element < len(shape):
                raise ValueError(
                    "maximum element is incompatible with the shape")

        if shape is None and weight is None:
            if max_element is not None:
                raise ValueError("specify shape or weight argument")
            return ShiftedPrimedTableaux_all()

        elif weight is None and all(shape[i] > shape[i+1]
                                    for i in range(len(shape)-1)):
            return ShiftedPrimedTableaux_shape(Partition(shape), max_element)

        elif shape is None:
            return ShiftedPrimedTableaux_weight(weight)

        if not all(shape[i] > shape[i+1] for i in range(len(shape)-1)):
            raise ValueError(
                "shape {} is not a strict partition".format(shape))

        if sum(shape) != sum(weight):
            raise ValueError(
                "weight and shape are incompatible")

        return ShiftedPrimedTableaux_weight_shape(weight, shape)

    def __contains__(self, T):
        """
        EXAMPLES::

            sage: [(1,'2p',2,2),(2,'3p')] in ShiftedPrimedTableaux()
            True
            sage: [(1,1),(2,2)]           in ShiftedPrimedTableaux()
            False
            sage: [(1,1),('2p',)]         in ShiftedPrimedTableaux()
            False
            sage: [(1,1,'3p'),(2,'3p')]   in ShiftedPrimedTableaux()
            True
            sage: [(1,1,2),(2,2)]         in ShiftedPrimedTableaux()
            False
            sage: [1,1,1] in ShiftedPrimedTableaux()
            True

        TESTS::

            sage: 1 in ShiftedPrimedTableaux()
            True
            sage: [] in ShiftedPrimedTableaux()
            True
        """
        try:
            self.element_class(self, T)
            return True
        except ValueError:
            return False

    _is_a = __contains__


class ShiftedPrimedTableaux_all(ShiftedPrimedTableaux):
    """
    The class of all Shifted Primed Tableaux.
    """
    Element = ShiftedPrimedTableau

    def __init__(self):
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
        Parent.__init__(self, category=FiniteEnumeratedSets())

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TEST::

            sage: ShiftedPrimedTableaux()
            Shifted Primed Tableaux
        """
        return "Shifted Primed Tableaux"

    def _element_constructor_(self, T):
        """
        Construct an object from ``T`` as an element of ``self``, if
        possible.

        INPUT:

        - ``T`` -- data which can be interpreted as a tableau

        OUTPUT:

        - the corresponding tableau object

        TESTS::

            sage: Tab=ShiftedPrimedTableaux()([1,1,"2p"]); Tab
            [(1.0, 1.0, 1.5)]
            sage: Tab.parent()
            Shifted Primed Tableaux
            sage: Tab=ShiftedPrimedTableaux()([[1,1,2],[2,2]])
            Traceback (most recent call last):
            ...
            ValueError: [[1, 1, 2], [2, 2]] is not an element of
            Shifted Primed Tableaux

        """
        try:
            return self.element_class(self, T)
        except ValueError:
            raise ValueError(
                "{} is not an element of Shifted Primed Tableaux".format(T))

    def __contains__(self, T):
        """
        TESTS::

            sage: [1,1,2] in ShiftedPrimedTableaux()
            True
            sage: (2,1) in ShiftedPrimedTableaux()
            False
        """
        try:
            self.element_class(self, T)
            return True
        except ValueError:
            return False


class ShiftedPrimedTableaux_shape(ShiftedPrimedTableaux):
    """
    Shifted Primed Tableaux of a fixed shape.
    """
    Element = ShiftedPrimedTableau

    def __init__(self, shape, max_elt):
        """
        Initialize the class of Shifted Primed Tableaux of a given shape.

        If ``max_elt`` is specified, a finite set with entries smaller
        or equal to ``max_elt``.

        """
        Parent.__init__(self, category=FiniteEnumeratedSets())
        self._max_elt = max_elt
        self._shape = shape

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TEST::

            sage: ShiftedPrimedTableaux([3,2,1])
            Shifted Primed Tableaux of shape [3, 2, 1]
        """
        if self._max_elt is None:
            return "Shifted Primed Tableaux of shape {}".format(self._shape)
        if self._max_elt is not None:
            return (
                "Shifted Primed Tableaux of shape {} and maximum element {}"
                .format(self._shape, self._max_elt))

    def __contains__(self, T):
        """
        TEST::

           sage: [[1,'2p',2,2],[2,'3p']] in ShiftedPrimedTableaux([4,2])
           True
        """

        try:
            Tab = self.element_class(self, T)
        except ValueError:
            return False

        if self._max_elt is None:
            return (Tab.shape() == self._shape)

        return (Tab.shape() == self._shape and
                Tab.max_element() <= self._max_elt)

    def _element_constructor_(self, T):
        """
        Construct an object from ``T`` as an element of ``self``, if
        possible.

        INPUT:

        - ``T`` -- data which can be interpreted as a primed tableau

        OUTPUT:

        - the corresponding primed tableau object

        TESTS::

            sage: tab= ShiftedPrimedTableaux([3])([1,1,1.5]); tab
            [(1.0, 1.0, 1.5)]
            sage: tab.parent()
            Shifted Primed Tableaux of shape [3]
            sage: ShiftedPrimedTableaux([3])([1,1])
            Traceback (most recent call last):
            ...
            ValueError: [1, 1] is not an element of Shifted Primed Tableaux
            of shape [3]
        """

        try:
            Tab = self.element_class(self, T)
        except ValueError:
            raise ValueError(
                "{} is not an element of Shifted Primed Tableaux".format(T))

        if self._max_elt is None:
            if Tab.shape() == self._shape:
                return Tab
        else:
            if (Tab.shape() == self._shape and
                    Tab.max_element() <= self._max_elt):
                return Tab
        raise ValueError("{} is not an element of {}".format(T, self))

    def shape(self):
        """
        Return the shape of the shifted tableaux ``self``.

        TEST::

            sage: ShiftedPrimedTableaux([6,4,3,1]).shape()
            [6, 4, 3, 1]
        """
        return self._shape

    def __iter__(self):
        """
        Iterate over ``self``, if ``max_element`` is specified.

        EXAMPLE::

            sage: Tabs = ShiftedPrimedTableaux([3,2], max_element=3)
            sage: Tabs[:4]
            [[(1.0, 1.0, 1.0), (2.0, 2.0)],
             [(1.0, 1.0, 1.0), (2.0, 3.0)],
             [(1.0, 1.0, 1.0), (2.0, 2.5)],
             [(1.0, 1.0, 1.0), (3.0, 3.0)]]
            sage: len(list(Tabs))
            24

        TEST::

            sage: Tabs = ShiftedPrimedTableaux([3,2])
            sage: Tabs[:3]
            Traceback (most recent call last):
            ...
            ValueError: set is infinite
        """
        if self._max_elt is None:
            raise ValueError("set is infinite")
        for weight in OrderedPartitions(sum(self._shape)+self._max_elt,
                                        k=self._max_elt):
            weight_n = tuple([w-1 for w in weight])
            for tab in ShiftedPrimedTableaux(shape=self._shape,
                                             weight=weight_n):
                yield (tab)

    def list_decreasing_weight(self):
        """
        List elements of ``self`` with weakly decreasing weight.

        EXAMPLE::

            sage: Tabs = ShiftedPrimedTableaux([2,1])
            sage: Tabs.list_decreasing_weight()
            [[(1.0, 1.0), (2.0,)], [(1.0, 2.0), (3.0,)], [(1.0, 1.5), (3.0,)]]
        """
        list_dw = []
        if self._max_elt is None:
            max_element = sum(self._shape)
        else:
            max_element = self._max_elt
        for weight in Partition(self._shape).dominated_partitions(
                rows=max_element):
            list_dw.extend(list(ShiftedPrimedTableaux(weight=tuple(weight),
                                                      shape=self._shape)))
        return list_dw

    def list_highest_weight(self):
        """
        List elements of ``self`` that are highest weight elements
        in the crystal.

        EXAMPLE::

            sage: Tabs = ShiftedPrimedTableaux([3,1])
            sage: Tabs.list_highest_weight()
            [[(1.0, 1.0, 1.0), (2.0,)],
             [(1.0, 1.0, 1.5), (2.0,)],
             [(1.0, 1.0, 2.5), (2.0,)]]
        """
        return [tab
                for tab in self.list_decreasing_weight()
                if tab.is_highest_weight()]


class ShiftedPrimedTableaux_weight(ShiftedPrimedTableaux):
    """
    Shifted Primed Tableaux of fixed weight.
    """
    Element = ShiftedPrimedTableau

    def __init__(self, weight):
        """
        Initialize the class of Shifted Primed Tableaux of a given weight.

        TEST::

            sage: TestSuite( ShiftedPrimedTableaux((3,2,1)) ).run()
        """
        Parent.__init__(self, category=FiniteEnumeratedSets())
        self._weight = weight

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TEST::

            sage: ShiftedPrimedTableaux((3,2,1))
            Shifted Primed Tableaux of weight (3, 2, 1)
        """
        return "Shifted Primed Tableaux of weight {}".format(self._weight)

    def __contains__(self, T):
        """
        TEST::

           sage: [[1,'2p',2,2],[2,'3p']] in ShiftedPrimedTableaux((1,4,1))
           True
        """
        try:
            Tab = self.element_class(self, T)
        except ValueError:
            return False
        return Tab.weight() == self._weight

    def _element_constructor_(self, T):
        """
        Construct an object from ``T`` as an element of ``self``, if
        possible.

        INPUT:

        - ``T`` -- data which can be interpreted as a primed tableau

        OUTPUT:

        - the corresponding primed tableau object

        TEST::

            sage: tab= ShiftedPrimedTableaux((2,1))([1,1,1.5]); tab
            [(1.0, 1.0, 1.5)]
            sage: tab.parent()
            Shifted Primed Tableaux of weight (2, 1)
        """
        try:
            Tab = self.element_class(self, T)
        except ValueError:
            raise ValueError(
                "{} is not an element of Shifted Primed Tableaux".format(T))
        if Tab.weight() == self._weight:
            return Tab
        raise ValueError("{} is not an element of {}".format(T, self))

    def __iter__(self):
        """
        Iterate over ``self``.

        EXAMPLE::

            sage: Tabs = ShiftedPrimedTableaux((2,3))
            sage: Tabs[:4]
            [[(1.0, 1.0, 2.0, 2.0, 2.0)],
             [(1.0, 1.0, 1.5, 2.0, 2.0)],
             [(1.0, 1.0, 2.0, 2.0), (2.0,)],
             [(1.0, 1.0, 1.5, 2.0), (2.0,)]]
            sage: len(list(Tabs))
            5
        """
        return (tab for shape_ in Partitions(sum(self._weight))
                if all(shape_[i] > shape_[i+1] for i in range(len(shape_)-1))
                for tab in ShiftedPrimedTableaux(shape=shape_,
                                                 weight=self._weight))


class ShiftedPrimedTableaux_weight_shape(ShiftedPrimedTableaux):
    """
    Shifted Primed Tableaux of the fixed weight and shape.
    """
    Element = ShiftedPrimedTableau

    def __init__(self, weight, shape):
        """
        Initialize the class of Shifted Primed Tableaux of the given weight
        and shape.
        """
        Parent.__init__(self, category=FiniteEnumeratedSets())
        self._weight = weight
        self._shape = shape

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TEST::

            sage: ShiftedPrimedTableaux([3,2,1],(4,2))
            Shifted Primed Tableaux of weight (4, 2) and shape [3, 2, 1]
        """
        return (
            "Shifted Primed Tableaux of weight {} and shape {}" .format(
                self._weight, self._shape))

    def __contains__(self, T):
        """
        TEST::

            sage: [[1,'2p',2,2],[2,'3p']] in ShiftedPrimedTableaux(
            ....: [4,2],(1,4,1))
            True
        """
        try:
            Tab = self.element_class(self, T)
        except ValueError:
            return False

        return (Tab.weight() == self._weight
                and Tab.shape() == self._shape)

    def _element_constructor_(self, T):
        """
        Construct an object from ``T`` as an element of ``self``, if
        possible.

        TESTS::

            sage: tab= ShiftedPrimedTableaux([3],(2,1))([1,1,1.5]); tab
            [(1.0, 1.0, 1.5)]
            sage: tab.parent()
            Shifted Primed Tableaux of weight (2, 1) and shape [3]
            sage: ShiftedPrimedTableaux([3],(2,1))([1,1])
            Traceback (most recent call last):
            ...
            ValueError: [1, 1] is not an element of Shifted Primed Tableaux
            of weight (2, 1) and shape [3]
        """
        try:
            Tab = self.element_class(self, T)
        except ValueError:
            raise ValueError(
                "{} is not an element of Shifted Primed Tableaux".format(T))

        if Tab.shape() == self._shape and Tab.weight() == self._weight:
            return Tab

        raise ValueError("{} is not an element of {}".format(T, self))

    def __iter__(self):
        """
        Iterate over ``self``.

        EXAMPLE::

            sage: Tabs = ShiftedPrimedTableaux([3,2], (1,2,2))
            sage: Tabs[:4]
            [[(1.0, 2.0, 2.0), (3.0, 3.0)],
             [(1.0, 1.5, 2.5), (2.0, 3.0)],
             [(1.0, 1.5, 2.5), (2.0, 2.5)],
             [(1.0, 1.5, 2.0), (3.0, 3.0)]]
            sage: len(list(Tabs))
            4
        """
        if not self._shape.dominates(Partition(sorted(list(self._weight),
                                                      key=int,
                                                      reverse=True))):
            return
            yield
        full_shape = self._shape
        sub_tab = []
        tab_list_new = [[]]
        for i, w in enumerate(self._weight):
            tab_list_old = tab_list_new
            tab_list_new = []
            for sub_tab in tab_list_old:
                sub_shape = [len(sub_tab[r]) for r in range(len(sub_tab))]
                for strip in add_strip(sub_shape, full_shape, w):
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
            yield(ShiftedPrimedTableau(tab))


###########
# Crystal #
###########


class ShiftedPrimedTableauxCrystal(UniqueRepresentation, Parent):
    r"""
    The class of crystals generated by Shifted Primed Tableaux of fixed shape.

    INPUT:

    -``n`` or ``rank`` -- a nonnegative integer
    -``shape`` -- a strictly decreasing partition of length at most ``n`` plus
    one

    This constructs a classical crystal of type `A_n` with highest weights
    corresponding to a given shape.

    If ``n`` is not given, the rank of the crystal is assumed to be the size of
    the partition ``shape`` minus one.

    The list of module generators consists of all elements of the crystal with
    nonincreasing weight.

    Crystal is constructed following operations described in [HPS17]_.

    .. SEEALSO::

        :class:`ShiftedPrimedTableaux`
        :class:`ShiftedPrimedTableau`

    EXAMPLES::

        sage: SPTC = crystals.ShiftedPrimedTableaux([3,2], 2)
        sage: T = SPTC.module_generators[-1]
        sage: T
        [(1.0, 1.0, 1.5), (2.0, 2.5)]
        sage: T.f(2)
        [(1.0, 1.0, 2.5), (2.0, 2.5)]
        sage: len(SPTC.module_generators)
        7
        sage: SPTC[0]
        [(1.0, 1.0, 1.0), (2.0, 2.0)]
        sage: SPTC.cardinality()
        24
    """
    @staticmethod
    def __classcall_private__(cls, *args, **kwargs):
        """
        Normalize the input.

        EXAMPLES::

            sage: crystals.ShiftedPrimedTableaux(n=2, shape=[4,2])
            Crystal of Shifted Primed Tableaux of type A_2 of shape (4, 2)
            sage: crystals.ShiftedPrimedTableaux([4,2], rank=2)
            Crystal of Shifted Primed Tableaux of type A_2 of shape (4, 2)
            sage: crystals.ShiftedPrimedTableaux([4,2])
            Crystal of Shifted Primed Tableaux of type A_5 of shape (4, 2)

        TESTS::
            sage: crystals.ShiftedPrimedTableaux()
            Traceback (most recent call last):
            ...
            ValueError: shape argument must be specified
            sage: crystals.ShiftedPrimedTableaux([4,4], 3)
            Traceback (most recent call last):
            ...
            ValueError: shape [4, 4] is not a strict partition
            sage: crystals.ShiftedPrimedTableaux([3,2,1], 1)
            Traceback (most recent call last):
            ...
            ValueError: invalid crystal rank
            sage: crystals.ShiftedPrimedTableaux([3,2,1])
            Crystal of Shifted Primed Tableaux of type A_5 of shape (3, 2, 1)
        """

        shape = None
        n = None
        if 'shape' in kwargs:
            shape = tuple(kwargs['shape'])
        if 'rank' in kwargs:
            n = int(kwargs['rank'])
        if 'n' in kwargs:
            n = int(kwargs['n'])
        if args:
            if isinstance(args[0], (list, tuple, Partition)):
                shape = tuple(args[0])
                if len(args) > 1 and isinstance(args[1], (int, Integer)):
                    n = args[1]
            else:
                if isinstance(args[0], (int, Integer)):
                    n = args[0]
                    if (len(args) > 1 and
                            isinstance(args[1], (list, tuple, Partition))):
                        shape = tuple(args[1])
        if shape is None:
            raise ValueError('shape argument must be specified')
        if n is None:
            n = sum(shape)-1

        if n+1 < len(shape):
            raise ValueError('invalid crystal rank')
        return SPTCrystal(shape=shape, n=n)


class SPTCrystal(ShiftedPrimedTableauxCrystal):

    """
    A factory class generating a classical crystal of Shifted Primed Tableaux.

    INPUT:

    -``n``-- a nonnegative integer
    -``shape``-- a strictly decreasing partition of length at most ``n`` plus
    one
    """
    Element = ShiftedPrimedTableau

    def __init__(self, shape, n):
        """
        Initialize the crystal of Shifted Primed Tableaux.

        TESTS::

            sage: SPTC = crystals.ShiftedPrimedTableaux([4,2], 2)
            sage: SPTC._cartan_type
            ['A', 2]
            sage: len(SPTC.module_generators)
            21
        """
        if shape is None or n is None:
            raise ValueError('shape and n must be specified')
        Parent.__init__(self, category=ClassicalCrystals())
        self.n = n
        self._shape = shape
        self._cartan_type = CartanType(['A', n])
        self.module_generators = ShiftedPrimedTableaux(
            shape=shape,
            max_element=n+1).list_decreasing_weight()

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TEST::

            sage: crystals.ShiftedPrimedTableaux([4,2], 2)
            Crystal of Shifted Primed Tableaux of type A_2 of shape (4, 2)
        """
        return ("Crystal of Shifted Primed Tableaux of type A_%s of shape "
                % (self.n) + str(self._shape))

####################
# Helper functions #
####################


def add_strip(sub_tab, full_tab, length):
    """
    Helper function used in the algorithm to generate all Shifted Primed
    Tableaux of the fixed weight and shape.
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
