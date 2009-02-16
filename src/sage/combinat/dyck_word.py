r"""
Dyck Words

AUTHORS:
    -- Mike Hansen
    -- Dan Drake (2008-05-30): DyckWordBacktracker support
"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
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
#*****************************************************************************
from combinat import CombinatorialClass, CombinatorialObject, catalan_number
from backtrack import GenericBacktracker

open_symbol = 1
close_symbol = 0

def replace_parens(x):
    """
    EXAMPLES:
        sage: from sage.combinat.dyck_word import replace_parens
        sage: replace_parens('(')
        1
        sage: replace_parens(')')
        0
        sage: replace_parens(1)
        Traceback (most recent call last):
        ...
        ValueError
    """
    if x == '(':
        return open_symbol
    elif x == ')':
        return close_symbol
    else:
        raise ValueError


def replace_symbols(x):
    """
    EXAMPLES:
        sage: from sage.combinat.dyck_word import replace_symbols
        sage: replace_symbols(1)
        '('
        sage: replace_symbols(0)
        ')'
        sage: replace_symbols(3)
        Traceback (most recent call last):
        ...
        ValueError
    """
    if x == open_symbol:
        return '('
    elif x == close_symbol:
        return ')'
    else:
        raise ValueError


def DyckWord(dw=None, noncrossing_partition=None):
    """
    Returns a Dyck word object

    EXAMPLES:
        sage: dw = DyckWord([1, 0, 1, 0]); dw
        [1, 0, 1, 0]
        sage: print dw
        ()()
        sage: print dw.height()
        1
        sage: dw.to_noncrossing_partition()
        [[1], [2]]

        sage: DyckWord('()()')
        [1, 0, 1, 0]

        sage: DyckWord(noncrossing_partition=[[1],[2]])
        [1, 0, 1, 0]
    """
    if noncrossing_partition is not None:
        return from_noncrossing_partition(noncrossing_partition)

    elif isinstance(dw, str):
        l = map(replace_parens, dw)
    else:
        l = dw

    if isinstance(l, DyckWord_class):
        return l
    elif l in DyckWords() or is_a_prefix(l):
        return DyckWord_class(l)
    else:
        raise ValueError, "invalid Dyck word"

class DyckWord_class(CombinatorialObject):
    def __str__(self):
        """
        Returns a string consisting of matched parentheses corresponding to
        the Dyck word.

        EXAMPLES:
            sage: print DyckWord([1, 0, 1, 0])
            ()()
            sage: print DyckWord([1, 1, 0, 0])
            (())
        """
        return "".join(map(replace_symbols, [x for x in self]))


    def size(self):
        """
        Returns the size of the Dyck word, which is the number
        of opening parentheses in the Dyck word.

        EXAMPLES:
            sage: DyckWord([1, 0, 1, 0]).size()
            2
            sage: DyckWord([1, 0, 1, 1, 0]).size()
            3
        """
        return len(filter(lambda x: x == open_symbol, self))

    def height(self):
        r"""
        Returns the height of the Dyck word.

        We view the Dyck word as a
        Dyck path from $(0,0)$ to $(n,0)$ in the first quadrant by letting
        '1's represent steps in the direction $(1,1)$ and '0's represent
        steps in the direction $(1,-1)$.

        The height is the maximum $y$-coordinate reached.

        EXAMPLES:
            sage: DyckWord([]).height()
            0
            sage: DyckWord([1,0]).height()
            1
            sage: DyckWord([1, 1, 0, 0]).height()
            2
            sage: DyckWord([1, 1, 0, 1, 0]).height()
            2
            sage: DyckWord([1, 1, 0, 0, 1, 0]).height()
            2
            sage: DyckWord([1, 0, 1, 0]).height()
            1
            sage: DyckWord([1, 1, 0, 0, 1, 1, 1, 0, 0, 0]).height()
            3
        """
        height = 0
        height_max = 0
        for letter in self:
            if letter == open_symbol:
                height += 1
                height_max = max(height, height_max)
            elif letter == close_symbol:
                height -= 1
        return height_max

    def associated_parenthesis(self, pos):
        """
        EXAMPLES:
            sage: DyckWord([1, 0]).associated_parenthesis(0)
            1
            sage: DyckWord([1, 0, 1, 0]).associated_parenthesis(0)
            1
            sage: DyckWord([1, 0, 1, 0]).associated_parenthesis(1)
            0
            sage: DyckWord([1, 0, 1, 0]).associated_parenthesis(2)
            3
            sage: DyckWord([1, 0, 1, 0]).associated_parenthesis(3)
            2
            sage: DyckWord([1, 1, 0]).associated_parenthesis(1)
            2
            sage: DyckWord([1, 1]).associated_parenthesis(0)
        """
        d = 0
        height = 0
        if pos >= len(self):
            raise ValueError, "invalid index"

        if self[pos] == open_symbol:
            d += 1
            height += 1
        elif self[pos] == close_symbol:
            d -= 1
            height -= 1
        else:
            raise ValueError, "unknown symbol %s"%self[pos-1]

        while height != 0:
            pos += d
            if pos < 0 or pos >= len(self):
                return None
            if self[pos] == open_symbol:
                height += 1
            elif self[pos] == close_symbol:
                height -= 1

        return pos

    def to_noncrossing_partition(self):
        """
        Bijection of Biane from Dyck words to non crossing partitions
        Thanks to Mathieu Dutour for describing the bijection.

        EXAMPLES:
            sage: DyckWord([1, 0]).to_noncrossing_partition()
            [[1]]
            sage: DyckWord([1, 1, 0, 0]).to_noncrossing_partition()
            [[1, 2]]
            sage: DyckWord([1, 1, 1, 0, 0, 0]).to_noncrossing_partition()
            [[1, 2, 3]]
            sage: DyckWord([1, 0, 1, 0, 1, 0]).to_noncrossing_partition()
            [[1], [2], [3]]
            sage: DyckWord([1, 1, 0, 1, 0, 0]).to_noncrossing_partition()
            [[2], [1, 3]]
        """
        partition = []
        stack = []
        i = 0
        p = 1

        #Invariants:
        # - self[i] = 0
        # - p is the number of opening parens at position i

        while i < len(self):
            stack.append(p)
            j = i + 1
            while j < len(self) and self[j] == close_symbol:
                j += 1

            #Now j points to the next 1 or past the end of self
            nz = j - (i+1) # the number of )'s between i and j
            if nz > 0:
                # Remove the nz last elements of stack and
                # make a new part in partition
                if nz > len(stack):
                    raise ValueError, "incorrect dyck word"

                partition.append( stack[-nz:] )

                stack = stack[: -nz]
            i = j
            p += 1

        if len(stack) > 0:
            raise ValueError, "incorrect dyck word"

        return partition

    def to_ordered_tree(self):
        """
        TESTS:
            sage: DyckWord([1, 1, 0, 0]).to_ordered_tree()
            Traceback (most recent call last):
            ...
            NotImplementedError: TODO
        """
        raise NotImplementedError, "TODO"

    def to_triangulation(self):
        """
        TESTS:
            sage: DyckWord([1, 1, 0, 0]).to_triangulation()
            Traceback (most recent call last):
            ...
            NotImplementedError: TODO
        """
        raise NotImplementedError, "TODO"

    def peaks(self):
        r"""
        Returns a list of the positions of the peaks of a Dyck word.  A
        peak is 1 followed by a 0.  Note that this does not agree with the
        definition given by Haglund in: The $q,t$-Catalan Numbers and the
        Space of Diagonal Harmonics: With an Appendix on the Combinatorics
        of Macdonald Polynomials - James Haglund, University of
        Pennsylvania, Philadelphia - AMS, 2008, 167 pp.

        EXAMPLES:
            sage: DyckWord([1, 0, 1, 0]).peaks()
            [0, 2]
            sage: DyckWord([1, 1, 0, 0]).peaks()
            [1]
            sage: DyckWord([1,1,0,1,0,1,0,0]).peaks() # Haglund's def gives 2
            [1, 3, 5]

        """
        return [i for i in range(len(self)-1) if self[i] == open_symbol and self[i+1] == close_symbol]

    def to_tableau(self):
        """
        EXAMPLES:
            sage: DyckWord([]).to_tableau()
            []
            sage: DyckWord([1, 0]).to_tableau()
            [[2], [1]]
            sage: DyckWord([1, 1, 0, 0]).to_tableau()
            [[3, 4], [1, 2]]
            sage: DyckWord([1, 0, 1, 0]).to_tableau()
            [[2, 4], [1, 3]]
            sage: DyckWord([1]).to_tableau()
            [[1]]
            sage: DyckWord([1, 0, 1]).to_tableau()
            [[2], [1, 3]]
        """
        open_positions = []
        close_positions = []

        for i in range(len(self)):
            if self[i] == open_symbol:
                open_positions.append(i+1)
            else:
                close_positions.append(i+1)
        return filter(lambda x: x != [],  [ close_positions, open_positions ])

    def a_statistic(self):
        """
        Returns the a-statistic for the Dyck word.

        One can view a balanced Dyck word as a lattice path from $(0,0)$ to
        $(n,n)$ in the first quadrant by letting '1's represent steps in
        the direction $(1,0)$ and '0's represent steps in the direction
        $(0,1)$.  The resulting path will remain weakly above the diagonal
        $y = x$.

        The a-statistic, or area statistic, is the number of complete
        squares in the integer lattice which are below the path and above
        the line $y = x$. The 'half-squares' directly above the line $y=x$
        do not contribute to this statistic.


        EXAMPLES:
            sage: dw = DyckWord([1,0,1,0])
            sage: dw.a_statistic() # 2 half-squares, 0 complete squares
            0

            sage: dw = DyckWord([1,1,1,0,1,1,1,0,0,0,1,1,0,0,1,0,0,0])
            sage: dw.a_statistic()
            19

            sage: DyckWord([1,1,1,1,0,0,0,0]).a_statistic()
            6
            sage: DyckWord([1,1,1,0,1,0,0,0]).a_statistic()
            5
            sage: DyckWord([1,1,1,0,0,1,0,0]).a_statistic()
            4
            sage: DyckWord([1,1,1,0,0,0,1,0]).a_statistic()
            3
            sage: DyckWord([1,0,1,1,0,1,0,0]).a_statistic()
            2
            sage: DyckWord([1,1,0,1,1,0,0,0]).a_statistic()
            4
            sage: DyckWord([1,1,0,0,1,1,0,0]).a_statistic()
            2
            sage: DyckWord([1,0,1,1,1,0,0,0]).a_statistic()
            3
            sage: DyckWord([1,0,1,1,0,0,1,0]).a_statistic()
            1
            sage: DyckWord([1,0,1,0,1,1,0,0]).a_statistic()
            1
            sage: DyckWord([1,1,0,0,1,0,1,0]).a_statistic()
            1
            sage: DyckWord([1,1,0,1,0,0,1,0]).a_statistic()
            2
            sage: DyckWord([1,1,0,1,0,1,0,0]).a_statistic()
            3
            sage: DyckWord([1,0,1,0,1,0,1,0]).a_statistic()
            0
        """
        above = 0
        diagonal = 0
        a = 0
        for move in self:
            if move == 1:
                above += 1
            elif move == 0:
                diagonal += 1
                a += above - diagonal
        return a

    def b_statistic(self):
        r"""
        Returns the b-statistic for the Dyck word.

        One can view a balanced Dyck word as a lattice path from $(0,0)$ to
        $(n,n)$ in the first quadrant by letting '1's represent steps in
        the direction $(0,1)$ and '0's represent steps in the direction
        $(1,0)$.  The resulting path will remain weakly above the diagonal
        $y = x$.

        We describe the b-statistic of such a path in terms of what is
        known as the "bounce path".

        We can think of our bounce path as describing the trail of a billiard
        ball shot West from (n, n), which "bounces" down whenever it
        encounters a vertical step and "bounces" left when it encounters the
        line y = x. The bouncing ball will strike the diagonal at places (0,
        0), (j_1, j_1), (j_2, j_2), ... , (j_r-1, j_r-1), (j_r, j_r) = (n, n).
        We define the b-statistic to be the sum sum_{i=1}^{r-1} j_i.

        EXAMPLES:
            sage: dw = DyckWord([1,1,1,0,1,1,1,0,0,0,1,1,0,0,1,0,0,0])
            sage: dw.b_statistic()
            7

            sage: DyckWord([1,1,1,1,0,0,0,0]).b_statistic()
            0
            sage: DyckWord([1,1,1,0,1,0,0,0]).b_statistic()
            1
            sage: DyckWord([1,1,1,0,0,1,0,0]).b_statistic()
            2
            sage: DyckWord([1,1,1,0,0,0,1,0]).b_statistic()
            3
            sage: DyckWord([1,0,1,1,0,1,0,0]).b_statistic()
            3
            sage: DyckWord([1,1,0,1,1,0,0,0]).b_statistic()
            1
            sage: DyckWord([1,1,0,0,1,1,0,0]).b_statistic()
            2
            sage: DyckWord([1,0,1,1,1,0,0,0]).b_statistic()
            1
            sage: DyckWord([1,0,1,1,0,0,1,0]).b_statistic()
            4
            sage: DyckWord([1,0,1,0,1,1,0,0]).b_statistic()
            3
            sage: DyckWord([1,1,0,0,1,0,1,0]).b_statistic()
            5
            sage: DyckWord([1,1,0,1,0,0,1,0]).b_statistic()
            4
            sage: DyckWord([1,1,0,1,0,1,0,0]).b_statistic()
            2
            sage: DyckWord([1,0,1,0,1,0,1,0]).b_statistic()
            6


        REFERENCES:
            [1] The $q,t$-Catalan Numbers and the Space of Diagonal
            Harmonics: With an Appendix on the Combinatorics of Macdonald
            Polynomials - James Haglund, University of Pennsylvania,
            Philadelphia - AMS, 2008, 167 pp.
        """
        x_pos = len(self)/2
        y_pos = len(self)/2

        b = 0

        mode = "left"
        makeup_steps = 0
        l = self._list[:]
        l.reverse()

        for move in l:
            #print x_pos, y_pos, mode, move
            if mode == "left":
                if move == 0:
                    x_pos -= 1
                elif move == 1:
                    y_pos -= 1
                    if x_pos == y_pos:
                        b += x_pos
                    else:
                        mode = "drop"
            elif mode == "drop":
                if move == 0:
                    makeup_steps += 1
                elif move == 1:
                    y_pos -= 1
                    if x_pos == y_pos:
                        b += x_pos
                        mode = "left"
                        x_pos -= makeup_steps
                        makeup_steps = 0

        return b


def DyckWords(k1=None, k2=None):
    r"""
    Returns the combinatorial class of Dyck words.  A Dyck word is a
    sequence $(w_1, ..., w_n)$ consisting of '1's and '0's, with the
    property that for any $i$ with $1 \le i \le n$, the sequence
    $(w_1,...,w_i)$ contains at least as many $1$s as $0$s.

    A Dyck word is balanced if the total number of '1's is equal to the
    total number of '0's.  The number of balanced Dyck words of length $2k$
    is given by the Catalan number $C_k$.

    EXAMPLES:
      If neither k1 nor k2 are specified, then DyckWords returns
      the combinatorial class of all balanced Dyck words.

        sage: DW = DyckWords(); DW
        Dyck words
        sage: [] in DW
        True
        sage: [1, 0, 1, 0] in DW
        True
        sage: [1, 1, 0] in DW
        False

      If just k1 is specified, then it returns the combinatorial
      class of balanced Dyck words with k1 opening parentheses and k1
      closing parentheses.
        sage: DW2 = DyckWords(2); DW2
        Dyck words with 2 opening parentheses and 2 closing parentheses
        sage: DW2.first()
        [1, 0, 1, 0]
        sage: DW2.last()
        [1, 1, 0, 0]
        sage: DW2.count()
        2
        sage: DyckWords(100).count() == catalan_number(100)
        True

      If k2 is specified in addition to k1, then it returns
      the combinatorial class of Dyck words with k1 opening
      parentheses and k2 closing parentheses.
        sage: DW32 = DyckWords(3,2); DW32
        Dyck words with 3 opening parentheses and 2 closing parentheses
        sage: DW32.list()
        [[1, 0, 1, 0, 1],
         [1, 0, 1, 1, 0],
         [1, 1, 0, 0, 1],
         [1, 1, 0, 1, 0],
         [1, 1, 1, 0, 0]]
    """
    if k1 is None and k2 is None:
        return DyckWords_all()
    else:
        if k1 < 0 or (k2 is not None and k2 < 0):
            raise ValueError, "k1 (= %s) and k2 (= %s) must be nonnegative, with k1 >= k2."%(k1, k2)
        if k2 is not None and k1 < k2:
            raise ValueError, "k1 (= %s) must be >= k2 (= %s)"%(k1, k2)
        if k2 is None:
            return DyckWords_size(k1, k1)
        else:
            return DyckWords_size(k1, k2)

class DyckWords_all(CombinatorialClass):
    def __init__(self):
        """
        TESTS:
            sage: DW = DyckWords()
            sage: DW == loads(dumps(DW))
            True
        """
        pass

    def __repr__(self):
        """
        TESTS:
            sage: repr(DyckWords())
            'Dyck words'
        """
        return "Dyck words"

    def __contains__(self, x):
        """
        TESTS:
            sage: [] in DyckWords()
            True
            sage: [1] in DyckWords()
            False
            sage: [0] in DyckWords()
            False
            sage: [1, 0] in DyckWords()
            True
        """
        if isinstance(x, DyckWord_class):
            return True

        if not isinstance(x, list):
            return False

        if len(x) % 2 != 0:
            return False

        return is_a(x)

    def list(self):
        """
        TESTS:
            sage: DyckWords().list()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError


class DyckWordBacktracker(GenericBacktracker):
    """DyckWordBacktracker: this class is an iterator for all Dyck words
    with n opening parentheses and n - endht closing parentheses using
    the backtracker class. It is used by the DyckWords_size class.

    This is not really meant to be called directly, partially because it
    fails in a couple corner cases: DWB(0) yields [0], not the empty
    word, and DWB(k, k+1) yields something (it shouldn't yield
    anything). This could be fixed with a sanity check in _rec(), but
    then we'd be doing the sanity check *every time* we generate new
    objects; instead, we do *one* sanity check in DyckWords and assume
    here that the sanity check has already been made.

    AUTHOR:
        -- Dan Drake (2008-05-30)
    """
    def __init__(self, k1, k2):
        """
        TESTS:
            sage: from sage.combinat.dyck_word import DyckWordBacktracker
            sage: len(list(DyckWordBacktracker(5, 5)))
            42
            sage: len(list(DyckWordBacktracker(6,4)))
            90
            sage: len(list(DyckWordBacktracker(7,0)))
            1
        """
        GenericBacktracker.__init__(self, [], (0, 0))
        # note that the comments in this class think of our objects as
        # Dyck paths, not words; having k1 opening parens and k2 closing
        # parens corresponds to paths of length k1 + k2 ending at height
        # k1 - k2.
        self.n = k1 + k2
        self.endht = k1 - k2

    def _rec(self, path, state):
        """
        TESTS:
            sage: from sage.combinat.dyck_word import DyckWordBacktracker
            sage: dwb = DyckWordBacktracker(3, 3)
            sage: list(dwb._rec([1,1,0],(3, 2)))
            [([1, 1, 0, 0], (4, 1), False), ([1, 1, 0, 1], (4, 3), False)]
            sage: list(dwb._rec([1,1,0,0],(4, 0)))
            [([1, 1, 0, 0, 1], (5, 1), False)]
            sage: list(DyckWordBacktracker(4, 4)._rec([1,1,1,1],(4, 4)))
            [([1, 1, 1, 1, 0], (5, 3), False)]
        """
        len, ht = state

        if len < self.n - 1:
            # if length is less than n-1, new path won't have length n, so
            # don't yield it, and keep building paths

            # if the path isn't too low and is not touching the x-axis, we can
            # yield a path with a downstep at the end
            if ht > (self.endht - (self.n - len)) and ht > 0:
                yield path + [0], (len + 1, ht - 1), False

            # if the path isn't too high, it can also take an upstep
            if ht < (self.endht + (self.n - len)):
                yield path + [1], (len + 1, ht + 1), False
        else:
            # length is n - 1, so add a single step (up or down,
            # according to current height and endht), don't try to
            # construct more paths, and yield the path
            if ht < self.endht:
                yield path + [1], None, True
            else:
                yield path + [0], None, True


class DyckWords_size(CombinatorialClass):
    def __init__(self, k1, k2=None):
        """
        TESTS:
            sage: DW4 = DyckWords(4)
            sage: DW4 == loads(dumps(DW4))
            True
            sage: DW42 = DyckWords(4,2)
            sage: DW42 == loads(dumps(DW42))
            True

        """
        self.k1 = k1
        self.k2 = k2


    def __repr__(self):
        """
        TESTS:
            sage: repr(DyckWords(4))
            'Dyck words with 4 opening parentheses and 4 closing parentheses'
        """
        return "Dyck words with %s opening parentheses and %s closing parentheses"%(self.k1, self.k2)

    def count(self):
        """
        Returns the number of Dyck words of size n, i.e. the n-th Catalan number.

        EXAMPLES:
            sage: DyckWords(4).count()
            14
            sage: ns = range(9)
            sage: dws = [DyckWords(n) for n in ns]
            sage: all([ dw.count() == len(dw.list()) for dw in dws])
            True
        """
        if self.k2 == self.k1:
            return catalan_number(self.k1)
        else:
            return len(self.list())

    def __contains__(self, x):
        """
        EXAMPLES:
             sage: [1, 0] in DyckWords(1)
             True
             sage: [1, 0] in DyckWords(2)
             False
             sage: [1, 1, 0, 0] in DyckWords(2)
             True
             sage: [1, 0, 0, 1] in DyckWords(2)
             False
             sage: [1, 0, 0, 1] in DyckWords(2,2)
             False
             sage: [1, 0, 1, 0] in DyckWords(2,2)
             True
             sage: [1, 0, 1, 0, 1] in DyckWords(3,2)
             True
             sage: [1, 0, 1, 1, 0] in DyckWords(3,2)
             True
             sage: [1, 0, 1, 1] in DyckWords(3,1)
             True
        """
        return is_a(x, self.k1, self.k2)

    def list(self):
        """
        Returns a list of all the Dyck words with k1 opening and k2
        closing parentheses.

        EXAMPLES:
            sage: DyckWords(0).list()
            [[]]
            sage: DyckWords(1).list()
            [[1, 0]]
            sage: DyckWords(2).list()
            [[1, 0, 1, 0], [1, 1, 0, 0]]
        """
        return list(self.iterator())

    def iterator(self):
        """
        Returns an iterator for Dyck words with k1 opening and k2
        closing parentheses.

        EXAMPLES:
            sage: [ w for w in DyckWords(0) ]
            [[]]
            sage: [ w for w in DyckWords(1) ]
            [[1, 0]]
            sage: [ w for w in DyckWords(2) ]
            [[1, 0, 1, 0], [1, 1, 0, 0]]
            sage: len([ 'x' for _ in DyckWords(5) ])
            42
        """
        if self.k1 == 0:
            yield DyckWord_class([])
        elif self.k2 == 0:
            yield DyckWord_class([ open_symbol for _ in range(self.k1) ])
        else:
            for w in DyckWordBacktracker(self.k1, self.k2):
                yield DyckWord_class(w)





def is_a_prefix(obj, k1 = None, k2 = None):
    """
    If k1 is specificied, then the object must have exactly k1 open
    symbols.  If k2 is also specified, then obj must have exactly
    k2 close symbols.

    EXAMPLES:
        sage: from sage.combinat.dyck_word import is_a_prefix
        sage: is_a_prefix([1,1,0])
        True
        sage: is_a_prefix([0,1,0])
        False
        sage: is_a_prefix([1,1,0],2,1)
        True
        sage: is_a_prefix([1,1,0],1,1)
        False
    """
    if k1 is not None and k2 is None:
        k2 = k1
    if k1 is not None and k1 < k2:
        raise ValueError, "k1 (= %s) must be >= k2 (= %s)"%(k1, k2)


    n_opens = 0
    n_closes = 0

    for p in obj:
        if p == open_symbol:
            n_opens += 1
        elif p == close_symbol:
            n_closes += 1
        else:
            return False

        if n_opens < n_closes:
            return False

    if k1 is None and k2 is None:
        return True
    elif k2 is None:
        return n_opens == k1
    else:
        return n_opens == k1 and n_closes == k2


def is_a(obj, k1 = None, k2 = None):
    """
    If k1 is specificied, then the object must have exactly k1 open
    symbols.  If k2 is also specified, then obj must have exactly
    k2 close symbols.

    EXAMPLES:
        sage: from sage.combinat.dyck_word import is_a
        sage: is_a([1,1,0,0])
        True
        sage: is_a([1,0,1,0])
        True
        sage: is_a([1,1,0,0],2)
        True
        sage: is_a([1,1,0,0],3)
        False
    """
    if k1 is not None and k2 is None:
        k2 = k1
    if k1 is not None and k1 < k2:
        raise ValueError, "k1 (= %s) must be >= k2 (= %s)"%(k1, k2)


    n_opens = 0
    n_closes = 0

    for p in obj:
        if p == open_symbol:
            n_opens += 1
        elif p == close_symbol:
            n_closes += 1
        else:
            return False

        if n_opens < n_closes:
            return False

    if k1 is None and k2 is None:
        return n_opens == n_closes
    elif k2 is None:
        return n_opens == n_closes and n_opens == k1
    else:
        return n_opens == k1 and n_closes == k2

def from_noncrossing_partition(ncp):
    """
    TESTS:
        sage: DyckWord(noncrossing_partition=[[1,2]]) # indirect doctest
        [1, 1, 0, 0]
        sage: DyckWord(noncrossing_partition=[[1],[2]])
        [1, 0, 1, 0]

        sage: dws = DyckWords(5).list()
        sage: ncps = map( lambda x: x.to_noncrossing_partition(), dws)
        sage: dws2 = map( lambda x: DyckWord(noncrossing_partition=x), ncps)
        sage: dws == dws2
        True
    """
    l = [ 0 ] * int( sum( [ len(v) for v in ncp ] ) )
    for v in ncp:
        l[v[-1]-1] = len(v)

    res = []
    for i in l:
        res += [ open_symbol ] + [close_symbol]*int(i)
    return DyckWord(res)

def from_ordered_tree(tree):
    """
    TESTS:
        sage: sage.combinat.dyck_word.from_ordered_tree(1)
        Traceback (most recent call last):
        ...
        NotImplementedError: TODO
    """
    raise NotImplementedError, "TODO"
