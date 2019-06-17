r"""
Robinson-Schensted-Knuth correspondence

AUTHORS:

- Travis Scrimshaw (2012-12-07): Initial version

Introduction
============

The Robinson-Schensted-Knuth (RSK) correspondence is most naturally 
stated as a bijection between generalized permutations (also known 
as two-line arrays, biwords, ...) and pairs of semi-standard Young 
tableaux `(P, Q)` of identical shape.

The basic operation in the RSK correspondence is a row insertion 
`P \leftarrow k`(where `P` is a given semi-standard Young tableau, 
and `k` is an integer). Different insertion algorithms have been 
implemented for the RSK correspondence which can be specify in 
the function.

EXAMPLES:

We can perform RSK and the inverse on a variety of objects::

    sage: p = Tableau([[1,2,2],[2]]); q = Tableau([[1,3,3],[2]])
    sage: gp = RSK_inverse(p, q); gp
    [[1, 2, 3, 3], [2, 1, 2, 2]]
    sage: RSK(*gp)
    [[[1, 2, 2], [2]], [[1, 3, 3], [2]]]
    sage: RSK([2,3,2,1,2,3])
    [[[1, 2, 2, 3], [2], [3]], [[1, 2, 5, 6], [3], [4]]]
    sage: RSK([2,3,2,1,2,3], insertion=RSK.rules.EG)
    [[[1, 2, 3], [2, 3], [3]], [[1, 2, 6], [3, 5], [4]]]
    sage: m = RSK_inverse(p, q, 'matrix'); m
    [0 1]
    [1 0]
    [0 2]
    sage: RSK(m)
    [[[1, 2, 2], [2]], [[1, 3, 3], [2]]]

Insertions currently available
------------------------------
The following insertion algorithms for RSK correspondence are currently 
available:

- RSK (:class:`~sage.combinat.rsk.RuleRSK`)
- Edelman-Greene insertion (:class:`~sage.combinat.rsk.RuleEG`), an algorithm 
  defined in [EG1987]_ Definition 6.20 (where it is referred to as 
  Coxeter-Knuth insertion).
- Hecke RSK algorithm (:class:`~sage.combinat.rsk.RuleHecke`) , defined 
  using the Hecke insertion studied in [BKSTY06]_ (but using rows instead 
  of columns).

Background
----------

Edelman-Greene insertion (:class:`~sage.combinat.rsk.RuleEG`)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For a reduced word of a permutation (i.e., an element of a type-`A` 
Coxeter group), one can use Edelman-Greene insertion, an algorithm 
defined in [EG1987]_ Definition 6.20 (where it is referred to as 
Coxeter-Knuth insertion). The Edelman-Greene insertion is similar to the 
standard row insertion except that if `k_i` and `k_i + 1` both exist in row 
`i`, we *only* set `k_{i+1} = k_i + 1` and continue.

Hecke RSK algorithm (:class:`~sage.combinat.rsk.RuleHecke`)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The Hecke RSK algorithm defined using the Hecke insertion studied in 
[BKSTY06]_ (but using rows instead of columns) proceeds similarly to 
the classical RSK algorithm. However, it is not clear in what generality 
it works; thus, following [BKSTY06]_, we shall assume that our biword 
`p` has top line `(1, 2, \ldots, n)` (or, at least, has its top line 
strictly increasing).

The Hecke RSK algorithm returns a pair of an increasing tableau and a 
set-valued standard tableau. If 
`p = ((j_0, k_0), (j_1, k_1), \ldots, (j_{\ell-1}, k_{\ell-1}))`,
then the algorithm recursively constructs pairs
`(P_0, Q_0), (P_1, Q_1), \ldots, (P_\ell, Q_\ell)` of tableaux.
The construction of `P_{t+1}` and `Q_{t+1}` from `P_t`, `Q_t`,
`j_t` and `k_t` proceeds as follows: Set `i = j_t`, `x = k_t`,
`P = P_t` and `Q = Q_t`. We are going to insert `x` into the
increasing tableau `P` and update the set-valued "recording
tableau" `Q` accordingly. As in the classical RSK algorithm, we
first insert `x` into row `1` of `P`, then into row `2` of the
resulting tableau, and so on, until the construction terminates.
The details are different: Suppose we are inserting `x` into
row `R` of `P`. If (Case 1) there exists an entry `y` in row `R`
such that `x < y`, then let `y` be the minimal such entry. We
replace this entry `y` with `x` if the result is still an
increasing tableau; in either subcase, we then continue
recursively, inserting `y` into the next row of `P`.
If, on the other hand, (Case 2) no such `y` exists, then we
append `x` to the end of `R` if the result is an increasing
tableau (Subcase 2.1), and otherwise (Subcase 2.2) do nothing.
Furthermore, in Subcase 2.1, we add the box that we have just
filled with `x` in `P` to the shape of `Q`, and fill it with
the one-element set `\{i\}`. In Subcase 2.2, we find the
bottommost box of the column containing the rightmost box of
row `R`, and add `i` to the entry of `Q` in this box (this
entry is a set, since `Q` is a set-valued). In either
subcase, we terminate the recursion, and set
`P_{t+1} = P` and `Q_{t+1} = Q`.

Implementing your own insertion rule
------------------------------------
The functions RSK() and RSK_inverse() are written so that it is easy to
implement insertion algorithms you come across in your research.

To implement your own insertion algorithm, you first need to import the 
base class for a rule::

    sage: from sage.combinat.rsk import Rule

Next, you need to implement the forward rule and backward rule for RSK() 
and RSK_inverse respectively. For more information, see 
:class:`~sage.combinat.rsk.Rule`.

TESTS:

Check that it is a correspondence between all types of input and
the input is preserved::

    sage: t1 = Tableau([[1, 2, 5], [3], [4]])
    sage: t2 = Tableau([[1, 2, 3], [4], [5]])
    sage: gp = RSK_inverse(t1, t2); gp
    [[1, 2, 3, 4, 5], [1, 4, 5, 3, 2]]
    sage: w = RSK_inverse(t1, t2, 'word'); w
    word: 14532
    sage: m = RSK_inverse(t1, t2, 'matrix'); m
    [1 0 0 0 0]
    [0 0 0 1 0]
    [0 0 0 0 1]
    [0 0 1 0 0]
    [0 1 0 0 0]
    sage: p = RSK_inverse(t1, t2, 'permutation'); p
    [1, 4, 5, 3, 2]
    sage: t1
    [[1, 2, 5], [3], [4]]
    sage: t2
    [[1, 2, 3], [4], [5]]
    sage: RSK(*gp) == [t1, t2]
    True
    sage: RSK(w) == [t1, t2]
    True
    sage: RSK(m) == [t1, t2]
    True
    sage: RSK(p) == [t1, t2]
    True
    sage: gp
    [[1, 2, 3, 4, 5], [1, 4, 5, 3, 2]]
    sage: w
    word: 14532
    sage: m
    [1 0 0 0 0]
    [0 0 0 1 0]
    [0 0 0 0 1]
    [0 0 1 0 0]
    [0 1 0 0 0]
    sage: p
    [1, 4, 5, 3, 2]

REFERENCES:

.. [Knu1970] Donald E. Knuth.
   *Permutations, matrices, and generalized Young tableaux*.
   Pacific J. Math. Volume 34, Number 3 (1970), pp. 709-727.
   http://projecteuclid.org/euclid.pjm/1102971948

.. [EG1987] Paul Edelman, Curtis Greene.
   *Balanced Tableaux*.
   Advances in Mathematics 63 (1987), pp. 42-99.
   http://www.sciencedirect.com/science/article/pii/0001870887900636

.. [BKSTY06] \A. Buch, A. Kresch, M. Shimozono, H. Tamvakis, and A. Yong.
   *Stable Grothendieck polynomials and* `K`-*theoretic factor sequences*.
   Math. Ann. **340** Issue 2, (2008), pp. 359--382.
   :arxiv:`math/0601514v1`.
"""
#*****************************************************************************
#       Copyright (C) 2012 Travis Scrimshaw <tscrim@ucdavis.edu>
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

from builtins import zip
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.sage_object import SageObject

from sage.structure.element import is_Matrix
from sage.matrix.all import matrix

class Rule(UniqueRepresentation):
    r"""
    Generic base class for a insertion rule for RSK correspondence.
    """
    def to_pair(self, obj1=None, obj2=None):
        r"""
        Returns an iterable two-array in pair form for row insertion.

        EXAMPLES::

            sage: from sage.combinat.rsk import Rule
            sage: Rule().to_pair([1, 2, 2, 2], [2, 1, 1, 2])
            [(1, 2), (2, 1), (2, 1), (2, 2)]
            sage: m = Matrix(ZZ, 3, 2, [0,1,1,0,0,2]) ; m
            [0 1]
            [1 0]
            [0 2]
            sage: Rule().to_pair(m)
            [(1, 2), (2, 1), (3, 2), (3, 2)]
        """
        if obj2 is None:
            try:
                itr = obj1._rsk_iter()
            except AttributeError:
                # If this is (something which looks like) a matrix
                #   then build the generalized permutation
                try:
                    t = []
                    b = []
                    for i, row in enumerate(obj1):
                        for j, mult in enumerate(row):
                            if mult > 0:
                                t.extend([i+1]*mult)
                                b.extend([j+1]*mult)
                    itr = zip(t, b)
                    obj1 = t
                    obj2 = b
                except TypeError:
                    itr = zip(range(1, len(obj1)+1), obj1)
                    obj2 = obj1
        else:
            itr = zip(obj1, obj2)
        return list(itr)

    def forward_rule(self, obj1, obj2, check_standard=False):
        r"""
        Returns a pair of tableau from forward insertion of the pair ``[obj1, obj2]``.
        
        INPUT:
          
        - ``obj1, obj2`` -- Can be one of the following:

          - A word in an ordered alphabet
          - Two lists of equal length representing a generalized permutation
          - Any object which has a method ``_rsk_iter()`` which returns an
            iterator over the object represented as generalized permutation or
            a pair of lists.

        -  ``check_standard`` -- (Default: ``False``) Check if either of the
            resulting tableaux is a standard tableau, and if so, typecast it
            as such

        """
        self._forward_verify_input(obj1, obj2)
        itr = self.to_pair(obj1, obj2)
        p = []       #the "insertion" tableau
        q = []       #the "recording" tableau
        lt = 0
        lb = 0
        for i, j in itr:
            self.insertion(i, j, p, q)
        return self._forward_formatOutput(p, q, check_standard)

    def backward_rule(self, p, q, output):
        r"""
        Returns the generalized permutation from the reverse insertion 
        of a pair of tableaux ``(p, q)``.

        INPUT:

        - ``p``, ``q`` -- Two tableaux of the same shape.

        - ``output`` -- (Default: ``'array'``) if ``q`` is semi-standard:

          - ``'array'`` -- as a two-line array (i.e. generalized permutation or
            biword)
          -  ``'matrix'`` -- as an integer matrix

          and if ``q`` is standard, we can have the output:

          - ``'word'`` -- as a word

          and additionally if ``p`` is standard, we can also have the output:

          - ``'permutation'`` -- as a permutation 
        """
        from sage.combinat.tableau import SemistandardTableaux
        self._backward_verify_input(p, q)
        # Make a copy of p since this is destructive to it
        p_copy = [list(row) for row in p]

        if q.is_standard():
            rev_word = [] # This will be our word in reverse
            d = dict((qij,i) for i, Li in enumerate(q) for qij in Li)
            # d is now a dictionary which assigns to each integer k the
            # number of the row of q containing k.

            for key in sorted(d, reverse=True): # Delete last entry from i-th row of p_copy
                self.rev_insertion(d[key], p_copy, rev_word, True)
            return self._backward_formatOutput(rev_word, None, output, p.is_standard(), q.is_standard())

        upper_row = []
        lower_row = []
        #upper_row and lower_row will be the upper and lower rows of the
        #generalized permutation we get as a result, but both reversed.
        d = {}
        for row, Li in enumerate(q):
            for col, val in enumerate(Li):
                if val in d:
                    d[val][col] = row
                else:
                    d[val] = {col: row}
        #d is now a double family such that for every integers k and j,
        #the value d[k][j] is the row i such that the (i, j)-th cell of
        #q is filled with k.
        for value, row_dict in sorted(d.items(), reverse=True, key=lambda x: x[0]):
            for key in sorted(row_dict, reverse=True):
                self.rev_insertion(row_dict[key], p_copy, lower_row, False)
                upper_row.append(value)
        return self._backward_formatOutput(lower_row, upper_row, output, p.is_standard(), q.is_standard())

    def _forward_formatOutput(self, p=None, q=None, check_standard=False):
        r"""
        Return final output as per additional parameters for forward rule.
        """
        from sage.combinat.tableau import SemistandardTableau, StandardTableau

        if check_standard:
            try:
                P = StandardTableau(p)
            except ValueError:
                P = SemistandardTableau(p)
            try:
                Q = StandardTableau(q)
            except ValueError:
                Q = SemistandardTableau(q)
            return [P, Q]
        return [SemistandardTableau(p), SemistandardTableau(q)]
    
    def _backward_formatOutput(self, lower_row=None, upper_row=None, output='array', p_is_standard=True, q_is_standard=True):
        r"""
        Return final output as per additional parameters for backward rule.
        """
        if q_is_standard:
            if output == 'word':
                from sage.combinat.words.word import Word
                return Word(reversed(lower_row))
            if output == 'matrix':
                return to_matrix(list(range(1, len(lower_row)+1)), list(reversed(lower_row)))
            if output == 'array':
                return [list(range(1, len(lower_row)+1)), list(reversed(lower_row))]
            raise ValueError("invalid output option")

        else:
            if output == 'matrix':
                return to_matrix(list(reversed(upper_row)), list(reversed(lower_row)))
            if output == 'array':
                return [list(reversed(upper_row)), list(reversed(lower_row))]
            if output in ['permutation', 'word']:
                raise TypeError("q must be standard to have a %s as valid output"%output)
            raise ValueError("invalid output option")

    def _forward_verify_input(self, obj1, obj2):
        r"""
        Return exception for invalid input in forward rule.
        """
        pass
    
    def _backward_verify_input(self, p, q):
        r"""
        Return exception for invalid input in backward rule.
        """
        from sage.combinat.tableau import SemistandardTableaux

        if p not in SemistandardTableaux():
            raise ValueError("p(=%s) must be a semistandard tableau"%p)

        if q not in SemistandardTableaux() and not q.is_standard():
            raise ValueError("q(=%s) must be a semistandard tableau"%q)

class RuleRSK(Rule):
    r"""
    A rule modeling the classical Robinson-Schensted-Knuth insertion.

    EXAMPLES::

        sage: RSK([1, 2, 2, 2], [2, 1, 1, 2], insertion=RSK.rules.RSK)
        [[[1, 1, 2], [2]], [[1, 2, 2], [2]]]
        sage: p = Tableau([[1,2,2],[2]]); q = Tableau([[1,3,3],[2]])
        sage: RSK_inverse(p, q, insertion=RSK.rules.RSK)
        [[1, 2, 3, 3], [2, 1, 2, 2]]

    For ``RSK`` and ``RSK_inverse``, ``RuleRSK`` behaves same as 
    :class:`~sage.combinat.rsk.Rule`. It is worth noting that in case of 
    ``RSK_inverse`` with ``output = 'permutation'`` ``RuleRSK`` returns 
    an object of class :class:`~sage.combinat.permutation.Permutation`.
    """

    def insertion(self, i, j, p, q):
        r"""
        Inserts the letter ``(i,j)`` from the bi-word to the current tableaux ``p`` and ``q``.
        """
        from bisect import bisect_right
        for r, qr in zip(p,q):
            if r[-1] > j:
                #Figure out where to insert j into the row r.  The
                #bisect command returns the position of the least
                #element of r greater than j.  We will call it y.
                y_pos = bisect_right(r, j)
                #Switch j and y
                j, r[y_pos] = r[y_pos], j
            else:
                break
        else:
            #We made through all of the rows of p without breaking
            #so we need to add a new row to p and q.
            r = []; p.append(r)
            qr = []; q.append(qr)
        r.append(j)
        qr.append(i) # Values are always inserted to the right

    def rev_insertion(self, i, p_copy, rev_word, q_is_standard):
        r"""
        Reverse bump the right-most letter from the `i^{th}` row of the 
        current tableaux p_copy and appends the removed entry from ``p_copy``
        to the list ``rev_word``.
        """
        from bisect import bisect_left
        if q_is_standard:
            x = p_copy[i].pop() # Always the right-most entry
            for row in reversed(p_copy[:i]):
                y_pos = bisect_left(row,x) - 1
                # switch x and y
                x, row[y_pos] = row[y_pos], x
            rev_word.append(x)
        else:
            x = p_copy[i].pop() # Always the right-most entry
            for row in reversed(p_copy[:i]):
                y = bisect_left(row,x) - 1
                x, row[y] = row[y], x
            rev_word.append(x)

    def _backward_formatOutput(self, lower_row=None, upper_row=None, output='array', p_is_standard=True, q_is_standard=True):
        r"""
        Return final output as per additional parameters for backward rule.
        """
        if q_is_standard and output == 'permutation':
            if not p_is_standard:
                raise TypeError("p must be standard to have a valid permutation as output")
            from sage.combinat.permutation import Permutation
            return Permutation(reversed(lower_row))
        else:
            return super(RuleRSK, self)._backward_formatOutput(lower_row, upper_row, output, p_is_standard, q_is_standard)


class RuleEG(Rule):
    r"""
    A rule modeling Edelman-Greene insertion.

    EXAMPLES::

        sage: RSK([2,3,2,1,2,3], insertion=RSK.rules.EG)
        [[[1, 2, 3], [2, 3], [3]], [[1, 2, 6], [3, 5], [4]]]
        sage: pq = RSK([2,1,2,3,2], insertion=RSK.rules.EG); pq
        [[[1, 2, 3], [2, 3]], [[1, 3, 4], [2, 5]]]
        sage: RSK(RSK_inverse(*pq, insertion=RSK.rules.EG, output='matrix'), insertion=RSK.rules.EG)
        [[[1, 2, 3], [2, 3]], [[1, 3, 4], [2, 5]]]
        sage: RSK_inverse(*pq, insertion=RSK.rules.EG)
        [[1, 2, 3, 4, 5], [2, 1, 2, 3, 2]]

    For ``RSK``, ``RuleEG`` provides a bijection from reduced words of 
    permutations/elements of a type-`A` Coxeter group to a pair of 
    semi-standard tableaux tableaux ([EG1987]_ Definition 2.1) of the same shape.

    For ``RSK_inverse``, ``RuleEG`` provides a bijection from a pair of 
    same shaped tableaux to reduced words of a generalized permutation. For 
    ``output = 'permutation'`` RuleEG returns the smallest permutation satisfying 
    the resulting reduced word.

    """

    def insertion(self, i, j, p, q):
        r"""
        Inserts the letter ``(i,j)`` from the bi-word to the current tableaux ``p`` and ``q``.
        """
        from bisect import bisect_right
        for r, qr in zip(p,q):
            if r[-1] > j:
                #Figure out where to insert j into the row r.  The
                #bisect command returns the position of the least
                #element of r greater than j.  We will call it y.
                y_pos = bisect_right(r, j)
                if r[y_pos] == j + 1 and y_pos > 0 and j == r[y_pos - 1]:
                    #Special bump: Nothing to do ejcept increment j by 1
                    j += 1
                else:
                    #Switch j and y
                    j, r[y_pos] = r[y_pos], j
            else:
                break
        else:
            #We made through all of the rows of p without breaking
            #so we need to add a new row to p and q.
            r = []; p.append(r)
            qr = []; q.append(qr)
        r.append(j)
        qr.append(i) # Values are always inserted to the right

    def rev_insertion(self, i, p_copy, rev_word, q_is_standard):
        r"""
        Reverse bump the right-most letter from the `i^{th}` row of the 
        current tableaux p_copy and appends the removed entry from ``p_copy``
        to the list ``rev_word``.
        
        TESTS:
        
        For non-standard p, q::

            sage: RSK_inverse(*RSK([1, 2, 3, 2, 1], insertion='EG'), insertion='EG')
            [[1, 2, 3, 4, 5], [1, 2, 3, 2, 1]]
            sage: RSK_inverse(*RSK([1, 1, 1, 2], [1, 2, 3, 4], insertion=RSK.rules.EG), insertion=RSK.rules.EG)
            [[1, 1, 1, 2], [1, 2, 3, 4]]
            sage: RSK_inverse(*RSK([1, 2, 3, 3], [2, 1, 2, 2], insertion='EG'), insertion='EG')
            [[1, 2, 3, 3], [2, 1, 2, 2]]
            
        Since the column reading of the insertion tableau from Edelman-Greene insertion 
        gives one of reduced words for the original permutation, we can also check for that

            sage: f = lambda p: reversed([x for row in reversed(p) for x in row])
            sage: g = lambda p: RSK(p.reduced_word(), insertion=RSK.rules.EG)[0]
            sage: all(p == Permutations(n).from_reduced_word(f(g(p))) for n in range(8) for p in Permutations(n))
            True
        """
        from bisect import bisect_left
        if q_is_standard:
            x = p_copy[i].pop() # Always the right-most entry
            for row in reversed(p_copy[:i]):
                y_pos = bisect_left(row,x) - 1
                if row[y_pos] == x - 1 and y_pos < len(row)-1 and row[y_pos+1] == x:
                    # Nothing to do except decrement x by 1.
                    # (Case 1 on p. 74 of Edelman-Greene [EG1987]_.)
                    x -= 1
                else:
                    # switch x and y
                    x, row[y_pos] = row[y_pos], x
            rev_word.append(x)
        else:
            x = p_copy[i].pop() # Always the right-most entry
            for row in reversed(p_copy[:i]):
                y = bisect_left(row,x) - 1
                if row[y] == x - 1 and y < len(row)-1 and row[y+1] == x:
                    x-=1
                else:
                    x, row[y] = row[y], x
            rev_word.append(x)

    def _backward_formatOutput(self, lower_row=None, upper_row=None, output='array', p_is_standard=True, q_is_standard=True):
        r"""
        Return final output as per additional parameters for backward rule.
        """
        if q_is_standard and output == 'permutation':
            n = 0
            if list(lower_row):
                n = max(list(lower_row)) + 1
            from sage.combinat.permutation import Permutations
            return Permutations(n).from_reduced_word(list(lower_row))
        else:
            return super(RuleEG, self)._backward_formatOutput(lower_row, upper_row, output, p_is_standard, q_is_standard)


class RuleHecke(Rule):
    r"""
    A rule modeling the Hecke insertion algorithm.

    EXAMPLES::

        sage: w = [5, 4, 1, 3, 4, 2, 5, 1, 2, 1, 4, 2, 4]
        sage: RSK(w, insertion=RSK.rules.Hecke)
        [[[1, 2, 4, 5], [2, 4, 5], [3, 5], [4], [5]],
         [[(1,), (4,), (5,), (7,)],
          [(2,), (9,), (11, 13)],
          [(3,), (12,)],
          [(6,)],
          [(8, 10)]]]
        sage: w = [5, 4, 1, 3, 4, 2, 5, 1, 2, 1, 4, 2, 4]
        sage: P,Q = RSK(w, insertion=RSK.rules.Hecke)
        sage: wp = RSK_inverse(P, Q, insertion=RSK.rules.Hecke, output='list'); wp
        [5, 4, 1, 3, 4, 2, 5, 1, 2, 1, 4, 2, 4]
        sage: wp == w
        True
    """
    def forward_rule(self, obj1, obj2, check_standard=False):
        r"""
        Return the Hecke insertion of the pair ``[obj1, obj2]``.
        
        INPUT:
        
        - ``obj1, obj2`` -- Can be one of the following:

          - A word in an ordered alphabet
          - An integer matrix
          - Two lists of equal length representing a generalized permutation
          - Any object which has a method ``_rsk_iter()`` which returns an
            iterator over the object represented as generalized permutation or
            a pair of lists.

        -  ``check_standard`` -- (Default: ``False``) Check if either of the
            resulting tableaux is a standard tableau, and if so, typecast it
            as such

        """
        from sage.combinat.tableau import SemistandardTableau, Tableau
        if obj2 is None:    
            obj2 = obj1
            obj1 = list(range(1, len(obj1) + 1))

        p = []       #the "insertion" tableau
        q = []       #the "recording" tableau

        for i, j in zip(obj1, obj2):
            self.insertion(i, j, p, q)
        
        return [SemistandardTableau(p), Tableau(q)]

    def backward_rule(self, p, q, output):
        r"""
        Return the reverse Hecke insertion of ``(p, q)``.

        INPUT:  
        
        - ``p``, ``q`` -- Two tableaux of the same shape.

        -  ``output`` -- (Default: ``'array'``) if ``q`` is semi-standard:

          - ``'array'`` -- as a two-line array (i.e. generalized permutation or
            biword)

          and if ``q`` is standard, we can have the output:

          - ``'word'`` -- as a word
          - ``list`` -- as a list
        """
        if p.shape() != q.shape():
            raise ValueError("p(=%s) and q(=%s) must have the same shape"%(p, q))
        from sage.combinat.tableau import SemistandardTableaux
        if p not in SemistandardTableaux():
            raise ValueError("p(=%s) must be a semistandard tableau"%p)

        # Make a copy of p and q since this is destructive to it
        p_copy = [list(row) for row in p]
        q_copy = [[list(v) for v in row] for row in q]
        # We shall work on these copies of p and q. Notice that p might get
        # some empty rows in the process; we do not bother pruning them, as
        # they do not matter.

        #upper_row and lower_row will be the upper and lower rows of the
        #generalized permutation we get as a result, but both reversed.
        upper_row = []
        lower_row = []
        d = {}
        for ri, row in enumerate(q):
            for ci, entry in enumerate(row):
                for val in entry:
                    if val in d:
                        d[val][ci] = ri
                    else:
                        d[val] = {ci: ri}
        #d is now a double family such that for every integers k and j,
        #the value d[k][j] is the row i such that the (i, j)-th cell of
        #q is filled with k.
        for value, row_dict in sorted(d.items(), key=lambda x: -x[0]):
            for i in sorted(row_dict.values(), reverse=True):
                # These are always the right-most entry
                should_be_value = q_copy[i][-1].pop()
                assert value == should_be_value
                self.rev_insertion(i, p_copy, q_copy, lower_row)
                upper_row.append(value)
        return self._backward_formatOutput(lower_row, upper_row, output)

    def insertion(self, i, j, p, q):
        r"""
        Inserts the letter ``(i,j)`` from the bi-word to the current tableaux ``p`` and ``q``.
        """
        from bisect import bisect_right
        for ir,r in enumerate(p):
            if r[-1] > j:
                #Figure out where to insert j into the row r.  The
                #bisect command returns the position of the least
                #element of r greater than j.  We will call it y.
                y_pos = bisect_right(r, j)
                y = r[y_pos]
                # Check to see if we can swap j for y
                if (y_pos == 0 or r[y_pos-1] < j) and (ir == 0 or p[ir-1][y_pos] < j):
                    r[y_pos] = j
                j = y
            else:
                # We must have len(p[ir-1]) > len(r), since j is coming
                # from the previous row.
                if r[-1] < j and (ir == 0 or p[ir-1][len(r)] < j):
                    # We can add a boj to the row
                    r.append(j)
                    q[ir].append((i,)) # Values are always inserted to the right
                else:
                    # We must append i to the bottom of this column
                    l = len(r) - 1
                    while ir < len(q) and len(q[ir]) > l:
                        ir += 1
                    q[ir-1][-1] = q[ir-1][-1] + (i,)
                break
        else:
            #We made through all of the rows of p without breaking
            #so we need to add a new row to p and q.
            p.append([j])
            q.append([(i,)])

    def rev_insertion(self, i, p_copy, q_copy, rev_word):
        r"""
        Reverse bump the right-most letter from the `i^{th}` row of the 
        current tableaux p_copy and appends the removed entry from ``p_copy``
        to the list ``rev_word``.
        """
        
        from bisect import bisect_left
        if not q_copy[i][-1]:
            # That is, if value was alone in cell q_copy[i][-1].
            q_copy[i].pop()
            x = p_copy[i].pop()
        else:
            x = p_copy[i][-1]
        while i > 0:
            i -= 1
            row = p_copy[i]
            y_pos = bisect_left(row,x) - 1
            y = row[y_pos]
            # Check to see if we can swap x for y
            if ((y_pos == len(row) - 1 or x < row[y_pos+1])
                and (i == len(p_copy) - 1 or len(p_copy[i+1]) <= y_pos
                     or x < p_copy[i+1][y_pos])):
                row[y_pos] = x
            x = y
        rev_word.append(x)

    def _backward_formatOutput(self, lower_row=None, upper_row=None, output='array'):
        r"""
        Return final output as per additional parameters for backward rule.
        """
        if output == 'array':
            return [list(reversed(upper_row)), list(reversed(lower_row))]
        is_standard = (upper_row == list(range(len(upper_row), 0, -1)))
        if output == 'word':
            if not is_standard:
                raise TypeError("q must be standard to have a %s as valid output"%output)
            from sage.combinat.words.word import Word
            return Word(reversed(lower_row))
        if output == 'list':
            if not is_standard:
                raise TypeError("q must be standard to have a %s as valid output"%output)
            return list(reversed(lower_row))
        raise ValueError("invalid output option")

class RuleDualRSK(Rule):
    r"""
    A rule modeling the dual RSK insertion algorithm.
    """
    def to_pair(self, obj1=None, obj2=None):
        if obj2 is None:
            try:
                for i, row in enumerate(obj1):
                    for j, mult in enumerate(row):
                        if mult > 1:
                            raise ValueError("the matrix is not a {0, 1}-matrix")
            except TypeError:
                return super(RuleDualRSK, self).to_pair(obj1, obj2)
        else:
            lt = 0
            lb = 0
            for t,b in zip(obj1, obj2):
                if lt == t and b == lb:
                    raise ValueError("invalid strict biword")
                lt = t
                lb = b

        return super(RuleDualRSK, self).to_pair(obj1, obj2)
        
    def insertion(self, i, j, p, q):
        r"""
        Inserts the letter ``(i,j)`` from the bi-word to the current tableaux ``p`` and ``q``.
        """
        from bisect import bisect_left
        for r, qr in zip(p,q):
            if r[-1] >= j:
                #Figure out where to insert j into the row r.  The
                #bisect command returns the position of the least
                #element of r greater than j.  We will call it y.
                y_pos = bisect_left(r, j)
                #Switch j and y
                j, r[y_pos] = r[y_pos], j
            else:
                break
        else:
            #We made through all of the rows of p without breaking
            #so we need to add a new row to p and q.
            r = []; p.append(r)
            qr = []; q.append(qr)
        r.append(j)
        qr.append(i) # Values are always inserted to the right

    def rev_insertion(self, i, p_copy, rev_word, q_is_standard):
        r"""
        Reverse bump the right-most letter from the `i^{th}` row of the 
        current tableaux p_copy and appends the removed entry from ``p_copy``
        to the list ``rev_word``.
        """
        from bisect import bisect_right
        if q_is_standard:
            x = p_copy[i].pop() # Always the right-most entry
            for row in reversed(p_copy[:i]):
                y_pos = bisect_right(row,x) - 1
                # switch x and y
                x, row[y_pos] = row[y_pos], x
            rev_word.append(x)
        else:
            x = p_copy[i].pop() # Always the right-most entry
            for row in reversed(p_copy[:i]):
                y = bisect_right(row,x) - 1
                x, row[y] = row[y], x
            rev_word.append(x)

    def _backward_formatOutput(self, lower_row=None, upper_row=None, output='array', p_is_standard=True, q_is_standard=True):
        r"""
        Return final output as per additional parameters for backward rule.
        """
        if q_is_standard and output == 'permutation':
            if not p_is_standard:
                raise TypeError("p must be standard to have a valid permutation as output")
            from sage.combinat.permutation import Permutation
            return Permutation(reversed(lower_row))
        else:
            return super(RuleDualRSK, self)._backward_formatOutput(lower_row, upper_row, output, p_is_standard, q_is_standard)

    def _forward_formatOutput(self, p=None, q=None, check_standard=False):
        r"""
        Return final output as per additional parameters for forward rule.
        """
        from sage.combinat.tableau import Tableau, StandardTableau

        if check_standard:
            try:
                P = StandardTableau(p)
            except ValueError:
                P = Tableau(p)
            try:
                Q = StandardTableau(q)
            except ValueError:
                Q = Tableau(q)
            return [P, Q]
        return [Tableau(p), Tableau(q)]

    def _backward_verify_input(self, p, q):
        r"""
        Return exception for invalid input in backward rule.
        """
        from sage.combinat.tableau import SemistandardTableaux

        if q not in SemistandardTableaux() and not q.is_standard():
            raise ValueError("q(=%s) must be a semistandard tableau"%q)


class InsertionRules(object):
    r"""
    Catalog of rules for growth diagrams.
    """
    RSK = RuleRSK
    EG = RuleEG
    Hecke = RuleHecke
    dualRSK = RuleDualRSK

#####################################################################

def RSK(obj1=None, obj2=None, insertion=InsertionRules.RSK, check_standard=False, **options):
    r"""
    Perform the Robinson-Schensted-Knuth (RSK) correspondence.

    The Robinson-Schensted-Knuth (RSK) correspondence (also known
    as the RSK algorithm) is most naturally stated as a bijection
    between generalized permutations (also known as two-line arrays,
    biwords, ...) and pairs of semi-standard Young tableaux `(P, Q)`
    of identical shape. The tableau `P` is known as the insertion
    tableau, and `Q` is known as the recording tableau.

    The basic operation is known as row insertion `P \leftarrow k`
    (where `P` is a given semi-standard Young tableau, and `k` is an
    integer). Row insertion is a recursive algorithm which starts by
    setting `k_0 = k`, and in its `i`-th step inserts the number `k_i`
    into the `i`-th row of `P` (we start counting the rows at `0`) by
    replacing the first integer greater than `k_i` in the row by `k_i`
    and defines `k_{i+1}` as the integer that has been replaced. If no
    integer greater than `k_i` exists in the `i`-th row, then `k_i` is
    simply appended to the row and the algorithm terminates at this
    point.

    Now the RSK algorithm, applied to a generalized permutation
    `p = ((j_0, k_0), (j_1, k_1), \ldots, (j_{\ell-1}, k_{\ell-1}))`
    (encoded as a lexicographically sorted list of pairs) starts by
    initializing two semi-standard tableaux `P_0` and `Q_0` as empty
    tableaux. For each nonnegative integer `t` starting at `0`, take
    the pair `(j_t, k_t)` from `p` and set
    `P_{t+1} = P_t \leftarrow k_t`, and define `Q_{t+1}` by adding a
    new box filled with `j_t` to the tableau `Q_t` at the same
    location the row insertion on `P_t` ended (that is to say, adding
    a new box with entry `j_t` such that `P_{t+1}` and `Q_{t+1}` have
    the same shape). The iterative process stops when `t` reaches the
    size of `p`, and the pair `(P_t, Q_t)` at this point is the image
    of `p` under the Robinson-Schensted-Knuth correspondence.

    This correspondence has been introduced in [Knu1970]_, where it has
    been referred to as "Construction A".

    For more information, see Chapter 7 in [Sta-EC2]_.

    We also note that integer matrices are in bijection with generalized
    permutations. Furthermore, we can convert any word `w` (and, in
    particular, any permutation) to a generalized permutation by
    considering the top line to be `(1, 2, \ldots, n)` where `n` is the
    length of `w`.

    The optional argument ``insertion`` allows to specify an alternative
    insertion procedure to be used instead of the standard
    Robinson-Schensted-Knuth insertion.
    
    INPUT:

    - ``obj1, obj2`` -- Can be one of the following:

      - A word in an ordered alphabet
      - An integer matrix
      - Two lists of equal length representing a generalized permutation
      - Any object which has a method ``_rsk_iter()`` which returns an
        iterator over the object represented as generalized permutation or
        a pair of lists.

    - ``insertion`` -- (Default: ``RSK``) The following types of insertion
      are currently supported:

      - ``RSK`` -- Robinson-Schensted-Knuth (:class:`~sage.combinat.rsk.RuleRSK`)
      - ``EG`` -- Edelman-Greene (only for reduced words of
        permutations/elements of a type-`A` Coxeter group)
        (:class:`~sage.combinat.rsk.RuleEG`)
      - ``Hecke`` -- Hecke insertion (only guaranteed for
        generalized permutations whose top row is strictly increasing)
        (:class:`~sage.combinat.rsk.RuleHecke`)

    - ``check_standard`` -- (Default: ``False``) Check if either of the
      resulting tableaux is a standard tableau, and if so, typecast it
      as such

    For precise information see the particular Rule class.

    EXAMPLES:

    If we only give one line, we treat the top line as being
    `(1, 2, \ldots, n)`::

        sage: RSK([3,3,2,4,1])
        [[[1, 3, 4], [2], [3]], [[1, 2, 4], [3], [5]]]
        sage: RSK(Word([3,3,2,4,1]))
        [[[1, 3, 4], [2], [3]], [[1, 2, 4], [3], [5]]]
        sage: RSK(Word([2,3,3,2,1,3,2,3]))
        [[[1, 2, 2, 3, 3], [2, 3], [3]], [[1, 2, 3, 6, 8], [4, 7], [5]]]

    With a generalized permutation::

        sage: RSK([1, 2, 2, 2], [2, 1, 1, 2])
        [[[1, 1, 2], [2]], [[1, 2, 2], [2]]]
        sage: RSK(Word([1,1,3,4,4]), [1,4,2,1,3])
        [[[1, 1, 3], [2], [4]], [[1, 1, 4], [3], [4]]]
        sage: RSK([1,3,3,4,4], Word([6,2,2,1,7]))
        [[[1, 2, 7], [2], [6]], [[1, 3, 4], [3], [4]]]

    If we give it a matrix::

        sage: RSK(matrix([[0,1],[2,1]]))
        [[[1, 1, 2], [2]], [[1, 2, 2], [2]]]

    We can also give it something looking like a matrix::

        sage: RSK([[0,1],[2,1]])
        [[[1, 1, 2], [2]], [[1, 2, 2], [2]]]

    As an example of Hecke insertion we construct Example 2.1
    in :arxiv:`0801.1319v2`::

        sage: w = [5, 4, 1, 3, 4, 2, 5, 1, 2, 1, 4, 2, 4]
        sage: RSK(w, insertion=RSK.rules.Hecke)
        [[[1, 2, 4, 5], [2, 4, 5], [3, 5], [4], [5]],
         [[(1,), (4,), (5,), (7,)],
          [(2,), (9,), (11, 13)],
          [(3,), (12,)],
          [(6,)],
          [(8, 10)]]]

    Using dual RSK insertion, where given only one line,
    we treat the top line as being `(1, 2, \ldots, n)`::

        sage: RSK([3,3,2,4,1], insertion=RSK.rules.dualRSK)
        [[[1, 4], [2], [3], [3]], [[1, 4], [2], [3], [5]]]
        sage: RSK(Word([3,3,2,4,1]), insertion=RSK.rules.dualRSK)
        [[[1, 4], [2], [3], [3]], [[1, 4], [2], [3], [5]]]
        sage: RSK(Word([2,3,3,2,1,3,2,3]), insertion=RSK.rules.dualRSK)
        [[[1, 2, 3], [2, 3], [2, 3], [3]], [[1, 2, 8], [3, 6], [4, 7], [5]]]

    Using dual RSK insertion with a strict biword::

        sage: RSK([1,1,2,4,4,5],[2,4,1,1,3,2], insertion=RSK.rules.dualRSK)
        [[[1, 2], [1, 3], [2, 4]], [[1, 1], [2, 4], [4, 5]]]
        sage: RSK([1,1,2,3,3,4,5],[1,3,2,1,3,3,2], insertion=RSK.rules.dualRSK)
        [[[1, 2, 3], [1, 2], [3], [3]], [[1, 1, 3], [2, 4], [3], [5]]]
        sage: RSK([1, 2, 2, 2], [2, 1, 2, 4], insertion=RSK.rules.dualRSK)
        [[[1, 2, 4], [2]], [[1, 2, 2], [2]]]
        sage: RSK(Word([1,1,3,4,4]), [1,4,2,1,3], insertion=RSK.rules.dualRSK)
        [[[1, 2, 3], [1], [4]], [[1, 1, 4], [3], [4]]]
        sage: RSK([1,3,3,4,4], Word([6,1,2,1,7]), insertion=RSK.rules.dualRSK)
        [[[1, 2, 7], [1], [6]], [[1, 3, 4], [3], [4]]]

    Using dual RSK insertion with a `\{0, 1\}`-matrix::

        sage: RSK(matrix([[0,1],[1,1]]), insertion=RSK.rules.dualRSK)
        [[[1, 2], [2]], [[1, 2], [2]]]

    We can also give it something looking like a matrix::

        sage: RSK([[0,1],[1,1]], insertion=RSK.rules.dualRSK)
        [[[1, 2], [2]], [[1, 2], [2]]]

    There is also :func:`~sage.combinat.rsk.RSK_inverse` which performs
    the inverse of the bijection on a pair of semistandard tableaux. We
    note that the inverse function takes 2 separate tableaux as inputs, so
    to compose with :func:`~sage.combinat.rsk.RSK`, we need to use the
    python ``*`` on the output::

        sage: RSK_inverse(*RSK([1, 2, 2, 2], [2, 1, 1, 2]))
        [[1, 2, 2, 2], [2, 1, 1, 2]]
        sage: P,Q = RSK([1, 2, 2, 2], [2, 1, 1, 2])
        sage: RSK_inverse(P, Q)
        [[1, 2, 2, 2], [2, 1, 1, 2]]

        sage: RSK_inverse(*RSK([1, 2, 2, 2], [2, 1, 2, 3], 
        ....:           insertion=RSK.rules.dualRSK), insertion=RSK.rules.dualRSK)
        [[1, 2, 2, 2], [2, 1, 2, 3]]
        sage: P,Q = RSK([1, 2, 2, 2], [2, 1, 2, 3], insertion=RSK.rules.dualRSK)
        sage: RSK_inverse(P, Q, insertion=RSK.rules.dualRSK)
        [[1, 2, 2, 2], [2, 1, 2, 3]]

    TESTS:

    Empty objects::

        sage: RSK(Permutation([]))
        [[], []]
        sage: RSK(Word([]))
        [[], []]
        sage: RSK(matrix([[]]))
        [[], []]
        sage: RSK([], [])
        [[], []]
        sage: RSK([[]])
        [[], []]
        sage: RSK(Word([]), insertion=RSK.rules.EG)
        [[], []]
        sage: RSK(Word([]), insertion=RSK.rules.Hecke)
        [[], []]

    Empty objects for dual RSK::

        sage: RSK(Permutation([]), insertion=RSK.rules.dualRSK)
        [[], []]
        sage: RSK(Word([]), insertion=RSK.rules.dualRSK)
        [[], []]
        sage: RSK(matrix([[]]), insertion=RSK.rules.dualRSK)
        [[], []]
        sage: RSK([], [], insertion=RSK.rules.dualRSK)
        [[], []]
        sage: RSK([[]], insertion=RSK.rules.dualRSK)
        [[], []]
    """
    from sage.combinat.tableau import SemistandardTableau, StandardTableau, Tableau 

    if isinstance(insertion, str):
        if insertion == 'RSK':
            insertion = RSK.rules.RSK
        elif insertion == 'EG':
            insertion = RSK.rules.EG
        elif insertion == 'hecke':
            insertion = RSK.rules.Hecke
        else:
            raise ValueError("invalid input")

    rule = insertion()
    if not isinstance(rule, Rule):
        raise TypeError("the insertion must be an instance of Rule")

    if obj1 is None and obj2 is None:
        if 'matrix' in options:
            obj1 = matrix(options['matrix'])
        else:
            raise ValueError("invalid input")

    if is_Matrix(obj1):
        obj1 = obj1.rows()
    if len(obj1) == 0:
        return [StandardTableau([]), StandardTableau([])]

    if obj2 is not None:
        if len(obj1) != len(obj2):
            raise ValueError("the two arrays must be the same length")
        # Check it is a generalized permutation
        lt = 0
        lb = 0
        for t,b in zip(obj1, obj2):
            if t < lt or (lt == t and b < lb):
                raise ValueError("invalid generalized permutation")
            lt = t
            lb = b

    output = rule.forward_rule(obj1, obj2, check_standard)
    return output

robinson_schensted_knuth = RSK

RSK.rules = InsertionRules

def RSK_inverse(p, q, output='array', insertion=InsertionRules.RSK):
    r"""
    Returns the generalized permutation corresponding to the pair of
    tableaux `(p,q)` under the inverse of the Robinson-Schensted-Knuth
    correspondence.

    For more information on the bijection, see :func:`RSK`.

    INPUT:

    - ``p``, ``q`` -- Two semi-standard tableaux of the same shape, or
      (in the case when Hecke insertion is used) an increasing tableau and
      a set-valued tableau of the same shape (see the note below for the
      format of the set-valued tableau)

    - ``output`` -- (Default: ``'array'``) if ``q`` is semi-standard:

      - ``'array'`` -- as a two-line array (i.e. generalized permutation or
        biword)
      -  ``'matrix'`` -- as an integer matrix

      and if ``q`` is standard, we can have the output:

      - ``'word'`` -- as a word

      and additionally if ``p`` is standard, we can also have the output:

      - ``'permutation'`` -- as a permutation

    - ``insertion`` -- (Default: ``RSK``) The insertion algorithm used in the
      bijection. Currently the following are supported:

      - ``RSK`` -- Robinson-Schensted-Knuth insertion (:class:`~sage.combinat.rsk.RuleRSK`)
      - ``EG`` -- Edelman-Greene insertion (:class:`~sage.combinat.rsk.RuleEG`)
      - ``Hecke`` -- Hecke insertion (:class:`~sage.combinat.rsk.RuleHecke`)

    For precise information see the particular Rule class.
    
    .. NOTE::

        In the case of Hecke insertion, the input variable ``q`` should
        be a set-valued tableau, encoded as a tableau whose entries are
        strictly increasing tuples of positive integers. Each such tuple
        encodes the set of its entries.

    EXAMPLES:

    If both ``p`` and ``q`` are standard::

        sage: t1 = Tableau([[1, 2, 5], [3], [4]])
        sage: t2 = Tableau([[1, 2, 3], [4], [5]])
        sage: RSK_inverse(t1, t2)
        [[1, 2, 3, 4, 5], [1, 4, 5, 3, 2]]
        sage: RSK_inverse(t1, t2, 'word')
        word: 14532
        sage: RSK_inverse(t1, t2, 'matrix')
        [1 0 0 0 0]
        [0 0 0 1 0]
        [0 0 0 0 1]
        [0 0 1 0 0]
        [0 1 0 0 0]
        sage: RSK_inverse(t1, t2, 'permutation')
        [1, 4, 5, 3, 2]
        sage: RSK_inverse(t1, t1, 'permutation')
        [1, 4, 3, 2, 5]
        sage: RSK_inverse(t2, t2, 'permutation')
        [1, 2, 5, 4, 3]
        sage: RSK_inverse(t2, t1, 'permutation')
        [1, 5, 4, 2, 3]

    If the first tableau is semistandard::

        sage: p = Tableau([[1,2,2],[3]]); q = Tableau([[1,2,4],[3]])
        sage: ret = RSK_inverse(p, q); ret
        [[1, 2, 3, 4], [1, 3, 2, 2]]
        sage: RSK_inverse(p, q, 'word')
        word: 1322

    In general::

        sage: p = Tableau([[1,2,2],[2]]); q = Tableau([[1,3,3],[2]])
        sage: RSK_inverse(p, q)
        [[1, 2, 3, 3], [2, 1, 2, 2]]
        sage: RSK_inverse(p, q, 'matrix')
        [0 1]
        [1 0]
        [0 2]

    If both ``p`` and ``q`` are standard, the dual RSK insertion
    behaves identically to the usual RSK insertion::

        sage: t1 = Tableau([[1, 2, 5], [3], [4]])
        sage: t2 = Tableau([[1, 2, 3], [4], [5]])
        sage: RSK_inverse(t1, t2, insertion=RSK.rules.dualRSK)
        [[1, 2, 3, 4, 5], [1, 4, 5, 3, 2]]
        sage: RSK_inverse(t1, t2, 'word', insertion=RSK.rules.dualRSK)
        word: 14532
        sage: RSK_inverse(t1, t2, 'matrix', insertion=RSK.rules.dualRSK)
        [1 0 0 0 0]
        [0 0 0 1 0]
        [0 0 0 0 1]
        [0 0 1 0 0]
        [0 1 0 0 0]
        sage: RSK_inverse(t1, t2, 'permutation', insertion=RSK.rules.dualRSK)
        [1, 4, 5, 3, 2]
        sage: RSK_inverse(t1, t1, 'permutation', insertion=RSK.rules.dualRSK)
        [1, 4, 3, 2, 5]
        sage: RSK_inverse(t2, t2, 'permutation', insertion=RSK.rules.dualRSK)
        [1, 2, 5, 4, 3]
        sage: RSK_inverse(t2, t1, 'permutation', insertion=RSK.rules.dualRSK)
        [1, 5, 4, 2, 3]

    For dual RSK, the first tableau is merely transpose semistandard::

        sage: p = Tableau([[1,2,2],[1]]); q = Tableau([[1,2,4],[3]])
        sage: ret = RSK_inverse(p, q, insertion=RSK.rules.dualRSK); ret
        [[1, 2, 3, 4], [1, 2, 1, 2]]
        sage: RSK_inverse(p, q, 'word', insertion=RSK.rules.dualRSK)
        word: 1212

    In general for dual RSK::

        sage: p = Tableau([[1,1,2],[1]]); q = Tableau([[1,3,3],[2]])
        sage: RSK_inverse(p, q, insertion=RSK.rules.dualRSK)
        [[1, 2, 3, 3], [1, 1, 1, 2]]
        sage: RSK_inverse(p, q, 'matrix', insertion=RSK.rules.dualRSK)
        [1 0]
        [1 0]
        [1 1]

    Using Hecke insertion::

        sage: w = [5, 4, 1, 3, 4, 2, 5, 1, 2, 1, 4, 2, 4]
        sage: pq = RSK(w, insertion=RSK.rules.Hecke)
        sage: RSK_inverse(*pq, insertion=RSK.rules.Hecke, output='list')
        [5, 4, 1, 3, 4, 2, 5, 1, 2, 1, 4, 2, 4]

    .. NOTE::

        The constructor of ``Tableau`` accepts not only semistandard
        tableaux, but also arbitrary lists that are fillings of a
        partition diagram. (And such lists are used, e.g., for the
        set-valued tableau ``q`` that is passed to
        ``RSK_inverse(p, q, insertion='hecke')``.)
        The user is responsible for ensuring that the tableaux passed to
        ``RSK_inverse`` are of the right types (semistandard, standard,
        increasing, set-valued as needed).

    TESTS:

    From empty tableaux::

        sage: RSK_inverse(Tableau([]), Tableau([]))
        [[], []]
        sage: RSK_inverse(Tableau([]), Tableau([]), insertion=RSK.rules.dualRSK)
        [[], []]

    Check that :func:`RSK_inverse` is the inverse of :func:`RSK` on the
    different types of inputs/outputs::

        sage: f = lambda p: RSK_inverse(*RSK(p), output='permutation')
        sage: all(p == f(p) for n in range(7) for p in Permutations(n))
        True
        sage: all(RSK_inverse(*RSK(w), output='word') == w for n in range(4) for w in Words(5, n))
        True
        sage: from sage.combinat.integer_matrices import IntegerMatrices
        sage: M = IntegerMatrices([1,2,2,1], [3,1,1,1])
        sage: all(RSK_inverse(*RSK(m), output='matrix') == m for m in M)
        True

        sage: n = ZZ.random_element(200)
        sage: p = Permutations(n).random_element()
        sage: is_fine = True if p == f(p) else p ; is_fine
        True

    Same for Edelman-Greene (but we are checking only the reduced words that
    can be obtained using the ``reduced_word()`` method from permutations)::

        sage: g = lambda w: RSK_inverse(*RSK(w, insertion=RSK.rules.EG), insertion=RSK.rules.EG, output='permutation')
        sage: all(p.reduced_word() == g(p.reduced_word()).reduced_word() for n in range(7) for p in Permutations(n))
        True

        sage: n = ZZ.random_element(200)
        sage: p = Permutations(n).random_element()
        sage: is_fine = True if p == f(p) else p ; is_fine
        True

    Same for dual RSK::

        sage: f = lambda p: RSK_inverse(*RSK(p, insertion=RSK.rules.dualRSK),
        ....:                           output='permutation', insertion=RSK.rules.dualRSK)
        sage: all(p == f(p) for n in range(7) for p in Permutations(n))
        True
        sage: all(RSK_inverse(*RSK(w, insertion=RSK.rules.dualRSK),
        ....:                 output='word', insertion=RSK.rules.dualRSK) == w
        ....:     for n in range(4) for w in Words(5, n))
        True
        sage: from sage.combinat.integer_matrices import IntegerMatrices
        sage: M = IntegerMatrices([1,2,2,1], [3,1,1,1]) # this is probably wrong
        sage: all(RSK_inverse(*RSK(m, insertion=RSK.rules.dualRSK), output='matrix',
        ....:                 insertion=RSK.rules.dualRSK) == m
        ....:     for m in M if all(x in [0, 1] for x in m))
        True

        sage: n = ZZ.random_element(200)
        sage: p = Permutations(n).random_element()
        sage: True if p == f(p) else p
        True

    Both tableaux must be of the same shape::

        sage: RSK_inverse(Tableau([[1,2,3]]), Tableau([[1,2]]))
        Traceback (most recent call last):
        ...
        ValueError: p(=[[1, 2, 3]]) and q(=[[1, 2]]) must have the same shape
        sage: RSK_inverse(Tableau([[1,2,3]]), Tableau([[1,2]]), insertion=RSK.rules.dualRSK)
        Traceback (most recent call last):
        ...
        ValueError: p(=[[1, 2, 3]]) and q(=[[1, 2]]) must have the same shape

    Check that :trac:`20430` is fixed::

        sage: RSK([1,1,1,1,1,1,1,2,2,2,3], [1,1,1,1,1,1,3,2,2,2,1])
        [[[1, 1, 1, 1, 1, 1, 1, 2, 2], [2], [3]],
         [[1, 1, 1, 1, 1, 1, 1, 2, 2], [2], [3]]]
        sage: t = SemistandardTableau([[1, 1, 1, 1, 1, 1, 1, 2, 2], [2], [3]])
        sage: RSK_inverse(t, t, 'array')
        [[1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 3],
         [1, 1, 1, 1, 1, 1, 3, 2, 2, 2, 1]]
    """

    if isinstance(insertion, str):
        if insertion == 'RSK':
            insertion = RSK.rules.RSK
        elif insertion == 'EG':
            insertion = RSK.rules.EG
        elif insertion == 'hecke':
            insertion = RSK.rules.Hecke
        else:
            raise ValueError("invalid input")

    rule = insertion()
    if not isinstance(rule, Rule):
        raise TypeError("the insertion must be an instance of Rule")

    if p.shape() != q.shape():
        raise ValueError("p(=%s) and q(=%s) must have the same shape"%(p, q))
    from sage.combinat.tableau import SemistandardTableaux

    answer = rule.backward_rule(p, q, output)
    return answer

robinson_schensted_knuth_inverse = RSK_inverse

def to_matrix(t, b):
    r"""
    Return the integer matrix corresponding to a two-line array.

    INPUT:

    - ``t`` -- The top line of the array

    - ``b`` -- The bottom line of the array

    OUTPUT:

    An `m \times n`-matrix (where `m` and `n` are the maximum entries in
    `t` and `b` respectively) whose `(i, j)`-th entry, for any `i` and `j`,
    is the number of all positions `k` satisfying `t_k = i` and `b_k = j`.

    EXAMPLES::

        sage: from sage.combinat.rsk import to_matrix
        sage: to_matrix([1, 1, 3, 3, 4], [2, 3, 1, 1, 3])
        [0 1 1]
        [0 0 0]
        [2 0 0]
        [0 0 1]
    """
    n = len(b)
    if len(t) != n:
        raise ValueError("The two arrays must be the same length")

    # Build the dictionary of entries since the matrix
    #   is typically (very) sparse
    entries = {}
    for i in range(n):
        if (t[i]-1, b[i]-1) in entries:
            entries[(t[i]-1, b[i]-1)] += 1
        else:
            entries[(t[i]-1, b[i]-1)] = 1
    return matrix(entries, sparse=True)
