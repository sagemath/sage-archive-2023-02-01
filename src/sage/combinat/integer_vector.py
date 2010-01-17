"""
Integer vectors

AUTHORS:

 *   Mike Hanson (2007) - original module
 *   Nathann Cohen, David Joyner (2009-2010) - Gale-Ryser stuff

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

from combinat import CombinatorialClass
from __builtin__ import list as builtinlist
from sage.rings.integer import Integer
from sage.rings.arith import binomial
import misc
from sage.rings.infinity import PlusInfinity
import integer_list
import cartesian_product
import functools


def is_gale_ryser(r,s):
    r"""
    Tests whether the given sequences satisfy the condition
    of the Gale-Ryser theorem.

    Given a binary matrix `B` of dimension `n\times m`, the
    vector of row sums is defined as the vector whose
    `i^{\mbox{th}}` component is equal to the sum of the `i^{\mbox{th}}`
    row in `A`. The vector of column sums is defined similarly.

    If, given a binary matrix, these two vectors are easy to compute,
    the Gale-Ryser theorem lets us decide whether, given two
    non-negative vectors `r,s`, there exists a binary matrix
    whose row/colum sums vectors are `r` and `s`.

    This functions answers accordingly.

    INPUT:

    - ``r``, ``s`` -- lists of non-negative integers.

    ALGORITHM:

    Without loss of generality, we can assume that:

    - The two given sequences do not contain any `0` ( which would
      correspond to an empty column/row )

    - The two given sequences are ordered in decreasing order
      (reordering the sequence of row (resp. column) sums amounts to
      reordering the rows (resp. columns) themselves in the matrix,
      which does not alter the columns (resp. rows) sums.

    We can then assume that `r` and `s` are partitions
    (see the corresponding class ``Partition``)

    If `r^*` denote the conjugate of `r`, the Gale-Ryser theorem
    asserts that a binary Matrix satisfying the constraints exists
    if and only if `s\preceq r^*`, where `\preceq` denotes
    the domination order on partitions.

    EXAMPLES::

        sage: from sage.combinat.integer_vector import is_gale_ryser
        sage: is_gale_ryser([4,2,2],[3,3,1,1])
        True
        sage: is_gale_ryser([4,2,1,1],[3,3,1,1])
        True
        sage: is_gale_ryser([3,2,1,1],[3,3,1,1])
        False

    REMARK: In the literature, what we are calling a
    Gale-Ryser sequence sometimes goes by the (rather
    generic-sounding) term ''realizable sequence''.
    """

    # The sequences only contan non-negative integers
    if [x for x in r if x<0] or [x for x in s if x<0]:
        return False

    # builds the corresponding partitions, i.e.
    # removes the 0 and sorts the sequences
    from sage.combinat.partition import Partition
    r2 = Partition(sorted([x for x in r if x>0], reverse=True))
    s2 = Partition(sorted([x for x in s if x>0], reverse=True))
    rstar = Partition(r2).conjugate()

    #                                same number of 1s           domination
    return len(rstar) <= len(s2) and  sum(r2) == sum(s2) and rstar.dominates(s)

def _slider01(A, t, k):
    r"""
    Assumes `A` is a `(0,1)`-matrix. For each of the
    `t` rows with highest row sums, this function
    returns a matrix `B` which is the same as `A` except that it
    has slid exactly one `1` in `A` over to the `k`-th
    column.

    This is a 'private' function for use in gale_ryser_theorem.

    INPUT:

    - ``A`` -- an `m\times n` (0,1) matrix
    - ``t``, ``k`` -- integers satisfying `0 < t < m`, `0 < k < n`

    OUTPUT:

    An `m\times n` (0,1) matrix, which is the same as `A` except
    that it has exactly one `1` in `A` slid over to the `k`-th
    column.

    EXAMPLES::

        sage: from sage.combinat.integer_vector import _slider01
        sage: A = matrix([[1,1,1,0],[1,1,1,0],[1,0,0,0],[1,0,0,0]])
        sage: A
        [1 1 1 0]
        [1 1 1 0]
        [1 0 0 0]
        [1 0 0 0]
        sage: _slider01(A, 1, 3)
        [1 1 1 0]
        [1 1 0 1]
        [1 0 0 0]
        [1 0 0 0]
        sage: _slider01(A, 3, 3)
        [1 1 1 0]
        [1 1 1 0]
        [1 0 0 0]
        [0 0 0 1]

    """
    from sage.matrix.constructor import matrix
    maxrows = []
    m = len(A.rows())
    n = len(A.columns())
    B = matrix(m,n)
    rowsums = [sum(A.row(i)) for i in range(m)]
    maxrows = [rowsums[i] for i in range(t)] # t rows with max rowsums
    maxRows = [] # want these to correspond to elts of maxrows
    B = [list(A.row(i)) for i in range(m)]
    for i in range(m):
        if sum(A.row(i)) in maxrows:
            maxRows.append(i)  # >= t rows numbers with max rowsums
    for j in range(len(maxrows)):
        i = maxRows[-1]
        rw = A.row(i) # to make mutable
        # now we want to move the rightmost left 1 to the k-th column
        rwlist = list(rw)
        if not(0 in rwlist):
            N = n-1
        else:
            N = rwlist.index(0)
        rwlist[int(N)-1] = 0
        rwlist[k] = 1
        B[i] = rwlist
    return matrix(B)

def gale_ryser_theorem(p1, p2, algorithm="ryser"):
        r"""
        Returns the binary matrix given by the Gale-Ryser theorem.

        The Gale Ryser theorem asserts that if `p_1,p_2` are two
        partitions of `n` of respective lengths `k_1,k_2`, then there is
        a binary `k_1\times k_2` matrix `M` such that `p_1` is the vector
        of row sums and `p_2` is the vector of column sums of `M`, if
        and only if `p_2` dominates `p_1`.

        INPUT:

        - ``p1`` -- the first partition of `n` (trailing 0's allowed)
        - ``p2`` -- the second partition of `n` (trailing 0's allowed)

        - ``algorithm`` -- two possible string values :

            - ``"ryser"`` (default) implements the construction due
              to Ryser [Ryser63]_.

            - ``"gale"`` implements the construction due to Gale [Gale57]_.


        OUTPUT:

        - A binary matrix if it exists, ``None`` otherwise.


        Gale's Algorithm:

        (Gale [Gale57]_): A matrix satisfying the constraints of its
        sums can be defined as the solution of the following
        Linear Program, which Sage knows how to solve (requires
        packages GLPK or CBC).

        .. MATH::

            \forall i&\sum_{j=1}^{k_2} b_{i,j}=p_{1,j}\\
            \forall i&\sum_{j=1}^{k_1} b_{j,i}=p_{2,j}\\
            &b_{i,j}\mbox{ is a binary variable}


        Ryser's Algorithm:

        (Ryser [Ryser63]_): The construction of an `m\times n` matrix `A=A_{r,s}`,
        due to Ryser, is described as follows. The
        construction works if and only if have `s\preceq r^*`.

        * Construct the `m\times n` matrix `B` from `r` by defining
          the `i`-th row of `B` to be the vector whose first `r_i`
          entries are `1`, and the remainder are 0's, `1\leq i\leq
          m`.  This maximal matrix `B` with row sum `r` and ones left
          justified has column sum `r^{*}`.

        * Shift the last `1` in certain rows of `B` to column `n` in
          order to achieve the sum `s_n`.  Call this `B` again.

          * The `1`'s in column n are to appear in those
            rows in which `A` has the largest row sums, giving
            preference to the bottom-most positions in case of ties.
          * Note: When this step automatically "fixes" other columns,
            one must skip ahead to the first column index
            with a wrong sum in the step below.

        * Proceed inductively to construct columns `n-1`, ..., `2`, `1`.

        * Set `A = B`. Return `A`.

        EXAMPLES:

        Computing the matrix for `p_1=p_2=2+2+1` ::

            sage: from sage.combinat.integer_vector import gale_ryser_theorem
            sage: p1 = [2,2,1]
            sage: p2 = [2,2,1]
            sage: print gale_ryser_theorem(p1, p2, algorithm="gale") # Optional - requires GLPK or CBC
            [0 1 1]
            [1 1 0]
            [1 0 0]

        Or for a non-square matrix with `p_1=3+3+2+1` and `p_2=3+2+2+1+1` ::

            sage: from sage.combinat.integer_vector import gale_ryser_theorem
            sage: p1 = [3,3,1,1]
            sage: p2 = [3,3,1,1]
            sage: gale_ryser_theorem(p1, p2)
            [1 1 1 0]
            [1 1 0 1]
            [1 0 0 0]
            [0 1 0 0]
            sage: p1 = [4,2,2]
            sage: p2 = [3,3,1,1]
            sage: gale_ryser_theorem(p1, p2)
            [1 1 1 1]
            [1 1 0 0]
            [1 1 0 0]
            sage: p1 = [4,2,2,0]
            sage: p2 = [3,3,1,1,0,0]
            sage: gale_ryser_theorem(p1, p2)
            [1 1 1 1 0 0]
            [1 1 0 0 0 0]
            [1 1 0 0 0 0]
            [0 0 0 0 0 0]
            sage: p1 = [3,3,2,1]
            sage: p2 = [3,2,2,1,1]
            sage: print gale_ryser_theorem(p1, p2, algorithm="gale") # Optional - requires GLPK or CBC
            [1 0 1 1 0]
            [1 0 1 0 1]
            [1 1 0 0 0]
            [0 1 0 0 0]

        With `0` in the sequences, and with unordered inputs ::

            sage: from sage.combinat.integer_vector import gale_ryser_theorem
            sage: gale_ryser_theorem([3,3,0,1,1,0], [3,1,3,1,0])
            [1 1 1 0 0]
            [1 0 1 1 0]
            [0 0 0 0 0]
            [1 0 0 0 0]
            [0 0 1 0 0]
            [0 0 0 0 0]

        REFERENCES:

        ..  [Ryser63] H. J. Ryser, Combinatorial Mathematics,
                Carus Monographs, MAA, 1963.
        ..  [Gale57] D. Gale, A theorem on flows in networks, Pacific J. Math.
                7(1957)1073-1082.
        """
        from sage.combinat.partition import Partition
        from sage.matrix.constructor import matrix

        if not(is_gale_ryser(p1,p2)):
            return False


        if algorithm=="ryser": # ryser's algorithm
            from sage.combinat.permutation import Permutation

            # Sorts the sequences if they are not, and remembers the permutation
            # applied
            tmp = sorted(enumerate(p1), reverse=True, key=lambda x:x[1])
            r = [x[1] for x in tmp if x[1]>0]
            r_permutation = [x-1 for x in Permutation([x[0]+1 for x in tmp]).inverse()]
            m = len(r)

            tmp = sorted(enumerate(p2), reverse=True, key=lambda x:x[1])
            s = [x[1] for x in tmp if x[1]>0]
            s_permutation = [x-1 for x in Permutation([x[0]+1 for x in tmp]).inverse()]
            n = len(s)

            rowsA0 = [[0]*n]*m
            for j in range(m):
                if j<m:
                    rowsA0[j] = [1]*r[j]+[0]*(n-r[j])
                else:
                    rowsA0[j] = [0]*n
            A0 = matrix(rowsA0)

            for j in range(1,n-1):
                # starts for loop, k = n-1, ..., 1
                # which finds the 1st column with
                # incorrect column sum. For that bad
                # column index, apply slider again
                for k in range(1,n):
                    if sum(A0.column(n-k))<>s[n-k]:
                         break
                A0 = _slider01(A0,s[n-k],n-k)

            # If we need to add empty rows/columns
            if len(p1)!=m:
                A0 = A0.stack(matrix([[0]*n]*(len(p1)-m)))

            if len(p2)!=n:
                A0 = A0.transpose().stack(matrix([[0]*len(p1)]*(len(p2)-n))).transpose()

            # Applying the permutations to get a matrix satisfying the
            # order given by the input
            A0 = A0.matrix_from_rows_and_columns(r_permutation, s_permutation)
            return A0

        elif algorithm == "gale":
          from sage.numerical.mip import MixedIntegerLinearProgram
          k1, k2=len(p1), len(p2)
          p = MixedIntegerLinearProgram()
          b = p.new_variable(dim=2)
          for (i,c) in enumerate(p1):
              p.add_constraint(sum([b[i][j] for j in xrange(k2)]),min=c,max=c)
          for (i,c) in enumerate(p2):
              p.add_constraint(sum([b[j][i] for j in xrange(k1)]),min=c,max=c)
          p.set_objective(None)
          p.set_binary(b)
          p.solve()
          b = p.get_values(b)
          M = [[0]*k2 for i in xrange(k1)]
          for i in xrange(k1):
              for j in xrange(k2):
                  M[i][j] = int(b[i][j])
          return matrix(M)

        else:
            raise ValueError("The only two algorithms available are \"gale\" and \"ryser\"")


def _default_function(l, default, i):
    """
    EXAMPLES::

        sage: from sage.combinat.integer_vector import _default_function
        sage: import functools
        sage: f = functools.partial(_default_function, [1,2,3], 99)
        sage: f(0)
        99
        sage: f(1)
        1
        sage: f(2)
        2
        sage: f(3)
        3
        sage: f(4)
        99
    """
    try:
        if i <= 0:
            return default
        return l[i-1]
    except IndexError:
        return default

infinity = PlusInfinity()
def list2func(l, default=None):
    """
    Given a list l, return a function that takes in a value i and
    return l[i-1]. If default is not None, then the function will
    return the default value for out of range i's.

    EXAMPLES::

        sage: f = sage.combinat.integer_vector.list2func([1,2,3])
        sage: f(1)
        1
        sage: f(2)
        2
        sage: f(3)
        3
        sage: f(4)
        Traceback (most recent call last):
        ...
        IndexError: list index out of range

    ::

        sage: f = sage.combinat.integer_vector.list2func([1,2,3], 0)
        sage: f(3)
        3
        sage: f(4)
        0
    """
    if default is None:
        return lambda i: l[i-1]
    else:
        return functools.partial(_default_function, l, default)


def constant_func(i):
    """
    Returns the constant function i.

    EXAMPLES::

        sage: f = sage.combinat.integer_vector.constant_func(3)
        sage: f(-1)
        3
        sage: f('asf')
        3
    """
    return lambda x: i

def IntegerVectors(n=None, k=None, **kwargs):
    """
    Returns the combinatorial class of integer vectors.

    EXAMPLES: If n is not specified, it returns the class of all
    integer vectors.

    ::

        sage: IntegerVectors()
        Integer vectors
        sage: [] in IntegerVectors()
        True
        sage: [1,2,1] in IntegerVectors()
        True
        sage: [1, 0, 0] in IntegerVectors()
        True

    If n is specified, then it returns the class of all integer vectors
    which sum to n.

    ::

        sage: IV3 = IntegerVectors(3); IV3
        Integer vectors that sum to 3

    Note that trailing zeros are ignored so that [3, 0] does not show
    up in the following list (since [3] does)

    ::

        sage: IntegerVectors(3, max_length=2).list()
        [[3], [2, 1], [1, 2], [0, 3]]

    If n and k are both specified, then it returns the class of integer
    vectors that sum to n and are of length k.

    ::

        sage: IV53 = IntegerVectors(5,3); IV53
        Integer vectors of length 3 that sum to 5
        sage: IV53.cardinality()
        21
        sage: IV53.first()
        [5, 0, 0]
        sage: IV53.last()
        [0, 0, 5]
        sage: IV53.random_element()
        [4, 0, 1]
    """
    if n is None:
        return IntegerVectors_all()
    elif k is None:
        return IntegerVectors_nconstraints(n,kwargs)
    else:
        if isinstance(k, builtinlist):
            return IntegerVectors_nnondescents(n,k)
        else:
            if len(kwargs) == 0:
                return IntegerVectors_nk(n,k)
            else:
                return IntegerVectors_nkconstraints(n,k,kwargs)


class IntegerVectors_all(CombinatorialClass):
    def __repr__(self):
        """
        EXAMPLES::

            sage: IntegerVectors()
            Integer vectors
        """
        return "Integer vectors"

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: [] in IntegerVectors()
            True
            sage: [3,2,2,1] in IntegerVectors()
            True
        """
        if not isinstance(x, builtinlist):
            return False
        for i in x:
            if not isinstance(i, (int, Integer)):
                return False
            if i < 0:
                return False

        return True

    def list(self):
        """
        EXAMPLES::

            sage: IntegerVectors().list()
            Traceback (most recent call last):
            ...
            NotImplementedError: infinite list
        """
        raise NotImplementedError, "infinite list"  # can't use InfiniteAbstractCombinatorialClass

    def cardinality(self):
        """
        EXAMPLES::

            sage: IntegerVectors().cardinality()
            +Infinity
        """
        return infinity


class IntegerVectors_nk(CombinatorialClass):
    def __init__(self, n, k):
        """
        TESTS::

            sage: IV = IntegerVectors(2,3)
            sage: IV == loads(dumps(IV))
            True

        AUTHORS:

        - Martin Albrecht

        - Mike Hansen
        """
        self.n = n
        self.k = k


    def _list_rec(self, n, k):
        """
        Return a list of a exponent tuples of length `size` such
        that the degree of the associated monomial is `D`.

        INPUT:


        -  ``n`` - degree (must be 0)

        -  ``k`` - length of exponent tuples (must be 0)


        EXAMPLES::

            sage: IV = IntegerVectors(2,3)
            sage: IV._list_rec(2,3)
            [(2, 0, 0), (1, 1, 0), (1, 0, 1), (0, 2, 0), (0, 1, 1), (0, 0, 2)]
        """
        res = []

        if k == 1:
            return [ (n, ) ]

        for nbar in range(n+1):
            n_diff = n-nbar
            for rest in self._list_rec( nbar , k-1):
                res.append((n_diff,)+rest)
        return res

    def list(self):
        """
        EXAMPLE::

            sage: IV = IntegerVectors(2,3)
            sage: IV.list()
            [[2, 0, 0], [1, 1, 0], [1, 0, 1], [0, 2, 0], [0, 1, 1], [0, 0, 2]]
            sage: IntegerVectors(3, 0).list()
            []
            sage: IntegerVectors(3, 1).list()
            [[3]]
            sage: IntegerVectors(0, 1).list()
            [[0]]
            sage: IntegerVectors(0, 2).list()
            [[0, 0]]
            sage: IntegerVectors(2, 2).list()
            [[2, 0], [1, 1], [0, 2]]
        """
        if self.n < 0:
            return []

        if self.k == 0:
            if self.n == 0:
                return [[]]
            else:
                return []
        elif self.k == 1:
            return [[self.n]]

        res = self._list_rec(self.n, self.k)
        return map(list, res)


    def __iter__(self):
        """
        EXAMPLE::

            sage: IV = IntegerVectors(2,3)
            sage: list(IV)
            [[2, 0, 0], [1, 1, 0], [1, 0, 1], [0, 2, 0], [0, 1, 1], [0, 0, 2]]
            sage: list(IntegerVectors(3, 0))
            []
            sage: list(IntegerVectors(3, 1))
            [[3]]
            sage: list(IntegerVectors(0, 1))
            [[0]]
            sage: list(IntegerVectors(0, 2))
            [[0, 0]]
            sage: list(IntegerVectors(2, 2))
            [[2, 0], [1, 1], [0, 2]]
            sage: IntegerVectors(0,0).list()
            [[]]
            sage: IntegerVectors(1,0).list()
            []
            sage: IntegerVectors(0,1).list()
            [[0]]
            sage: IntegerVectors(2,2).list()
            [[2, 0], [1, 1], [0, 2]]
            sage: IntegerVectors(-1,0).list()
            []
            sage: IntegerVectors(-1,2).list()
            []
        """
        if self.n < 0:
            return

        if self.k == 0:
            if self.n == 0:
                yield []
            return
        elif self.k == 1:
            yield [self.n]
            return

        for nbar in range(self.n+1):
            n = self.n-nbar
            for rest in IntegerVectors_nk(nbar , self.k-1):
                yield [n] + rest

    def __repr__(self):
        """
        TESTS::

            sage: IV = IntegerVectors(2,3)
            sage: repr(IV)
            'Integer vectors of length 3 that sum to 2'
        """
        return "Integer vectors of length %s that sum to %s"%(self.k, self.n)

    def __contains__(self, x):
        """
        TESTS::

            sage: IV = IntegerVectors(2,3)
            sage: all([i in IV for i in IV])
            True
            sage: [0,1,2] in IV
            False
            sage: [2.0, 0, 0] in IV
            False
            sage: [0,1,0,1] in IV
            False
            sage: [0,1,1] in IV
            True
            sage: [-1,2,1] in IV
            False
        """
        if x not in IntegerVectors():
            return False

        if sum(x) != self.n:
            return False

        if len(x) != self.k:
            return False

        if len(x) > 0 and min(x) < 0:
            return False

        return True

class IntegerVectors_nkconstraints(CombinatorialClass):
    def __init__(self, n, k, constraints):
        """
        EXAMPLES::

            sage: IV = IntegerVectors(2,3,min_slope=0)
            sage: IV == loads(dumps(IV))
            True
        """
        self.n = n
        self.k = k
        self.constraints = constraints

    def __repr__(self):
        """
        EXAMPLES::

            sage: IntegerVectors(2,3,min_slope=0).__repr__()
            'Integer vectors of length 3 that sum to 2 with constraints: min_slope=0'
        """
        return "Integer vectors of length %s that sum to %s with constraints: %s"%(self.k, self.n, ", ".join( ["%s=%s"%(key, self.constraints[key]) for key in sorted(self.constraints.keys())] ))


    def __contains__(self, x):
        """
        TESTS::

            sage: [0] in IntegerVectors(0)
            True
            sage: [0] in IntegerVectors(0, 1)
            True
            sage: [] in IntegerVectors(0, 0)
            True
            sage: [] in IntegerVectors(0, 1)
            False
            sage: [] in IntegerVectors(1, 0)
            False
            sage: [3] in IntegerVectors(3)
            True
            sage: [3] in IntegerVectors(2,1)
            False
            sage: [3] in IntegerVectors(2)
            False
            sage: [3] in IntegerVectors(3,1)
            True
            sage: [3,2,2,1] in IntegerVectors(9)
            False
            sage: [3,2,2,1] in IntegerVectors(9,5)
            False
            sage: [3,2,2,1] in IntegerVectors(8)
            True
            sage: [3,2,2,1] in IntegerVectors(8,5)
            False
            sage: [3,2,2,1] in IntegerVectors(8,4)
            True
            sage: [3,2,2,1] in IntegerVectors(8,4, min_part = 1)
            True
            sage: [3,2,2,1] in IntegerVectors(8,4, min_part = 2)
            False
        """
        if x not in IntegerVectors():
            return False

        if sum(x) != self.n:
            return False

        if len(x) != self.k:
            return False

        if self.constraints:
            if not misc.check_integer_list_constraints(x, singleton=True, **self.constraints):
                return False

        return True

    def cardinality(self):
        """
        EXAMPLES::

            sage: IntegerVectors(3,3, min_part=1).cardinality()
            1
            sage: IntegerVectors(5,3, min_part=1).cardinality()
            6
            sage: IntegerVectors(13, 4, min_part=2, max_part=4).cardinality()
            16
        """
        if not self.constraints:
            if self.n >= 0:
                return binomial(self.n+self.k-1,self.n)
            else:
                return 0
        else:
            if len(self.constraints) == 1 and 'max_part' in self.constraints and self.constraints['max_part'] != infinity:
                m = self.constraints['max_part']
                if m >= self.n:
                    return binomial(self.n+self.k-1,self.n)
                else: #do by inclusion / exclusion on the number
                      #i of parts greater than m
                    return sum( [(-1)**i * binomial(self.n+self.k-1-i*(m+1), self.k-1)*binomial(self.k,i) for i in range(0, self.n/(m+1)+1)])
            else:
                return len(self.list())


    def _parameters(self):
        """
        Returns a tuple (min_length, max_length, floor, ceiling,
        min_slope, max_slope) for the parameters of self.

        EXAMPLES::

            sage: IV = IntegerVectors(2,3,min_slope=0)
            sage: min_length, max_length, floor, ceiling, min_slope, max_slope = IV._parameters()
            sage: min_length
            3
            sage: max_length
            3
            sage: [floor(i) for i in range(1,10)]
            [0, 0, 0, 0, 0, 0, 0, 0, 0]
            sage: [ceiling(i) for i in range(1,5)]
            [+Infinity, +Infinity, +Infinity, +Infinity]
            sage: min_slope
            0
            sage: max_slope
            +Infinity

            sage: IV = IntegerVectors(3,10,inner=[4,1,3], min_part = 2)
            sage: min_length, max_length, floor, ceiling, min_slope, max_slope = IV._parameters()
            sage: floor(1), floor(2), floor(3)
            (4, 2, 3)

            sage: IV = IntegerVectors(3, 10, outer=[4,1,3], max_part = 3)
            sage: min_length, max_length, floor, ceiling, min_slope, max_slope = IV._parameters()
            sage: ceiling(1), ceiling(2), ceiling(3)
            (3, 1, 3)
        """
        constraints = self.constraints
        #n, min_length, max_length, floor, ceiling, min_slope, max_slope
        if self.k == -1:
            min_length = constraints.get('min_length', 0)
            max_length = constraints.get('max_length', infinity)
        else:
            min_length = self.k
            max_length = self.k

        min_part = constraints.get('min_part', 0)
        max_part = constraints.get('max_part', infinity)
        min_slope = constraints.get('min_slope', -infinity)
        max_slope = constraints.get('max_slope', infinity)
        if 'outer' in self.constraints:
            ceiling = list2func( map(lambda i: min(max_part, i), self.constraints['outer']), default=max_part )
        else:
            ceiling = constant_func(max_part)

        if 'inner' in self.constraints:
            floor = list2func( map(lambda i: max(min_part, i), self.constraints['inner']), default=min_part )
        else:
            floor = constant_func(min_part)

        return (min_length, max_length, floor, ceiling, min_slope, max_slope)


    def first(self):
        """
        EXAMPLES::

            sage: IntegerVectors(2,3,min_slope=0).first()
            [0, 1, 1]
        """
        return integer_list.first(self.n, *self._parameters())

    def next(self, x):
        """
        EXAMPLES::

            sage: IntegerVectors(2,3,min_slope=0).last()
            [0, 0, 2]
        """
        return integer_list.next(x, *self._parameters())

    def __iter__(self):
        """
        EXAMPLES::

            sage: IntegerVectors(-1, 0, min_part = 1).list()
            []
            sage: IntegerVectors(-1, 2, min_part = 1).list()
            []
            sage: IntegerVectors(0, 0, min_part=1).list()
            [[]]
            sage: IntegerVectors(3, 0, min_part=1).list()
            []
            sage: IntegerVectors(0, 1, min_part=1).list()
            []
            sage: IntegerVectors(2, 2, min_part=1).list()
            [[1, 1]]
            sage: IntegerVectors(2, 3, min_part=1).list()
            []
            sage: IntegerVectors(4, 2, min_part=1).list()
            [[3, 1], [2, 2], [1, 3]]

        ::

            sage: IntegerVectors(0, 3, outer=[0,0,0]).list()
            [[0, 0, 0]]
            sage: IntegerVectors(1, 3, outer=[0,0,0]).list()
            []
            sage: IntegerVectors(2, 3, outer=[0,2,0]).list()
            [[0, 2, 0]]
            sage: IntegerVectors(2, 3, outer=[1,2,1]).list()
            [[1, 1, 0], [1, 0, 1], [0, 2, 0], [0, 1, 1]]
            sage: IntegerVectors(2, 3, outer=[1,1,1]).list()
            [[1, 1, 0], [1, 0, 1], [0, 1, 1]]
            sage: IntegerVectors(2, 5, outer=[1,1,1,1,1]).list()
            [[1, 1, 0, 0, 0],
             [1, 0, 1, 0, 0],
             [1, 0, 0, 1, 0],
             [1, 0, 0, 0, 1],
             [0, 1, 1, 0, 0],
             [0, 1, 0, 1, 0],
             [0, 1, 0, 0, 1],
             [0, 0, 1, 1, 0],
             [0, 0, 1, 0, 1],
             [0, 0, 0, 1, 1]]

        ::

            sage: iv = [ IntegerVectors(n,k) for n in range(-2, 7) for k in range(7) ]
            sage: all(map(lambda x: x.cardinality() == len(x.list()), iv))
            True
            sage: essai = [[1,1,1], [2,5,6], [6,5,2]]
            sage: iv = [ IntegerVectors(x[0], x[1], max_part = x[2]-1) for x in essai ]
            sage: all(map(lambda x: x.cardinality() == len(x.list()), iv))
            True
        """
        return integer_list.iterator(self.n, *self._parameters())

class IntegerVectors_nconstraints(IntegerVectors_nkconstraints):
    def __init__(self, n, constraints):
        """
        TESTS::

            sage: IV = IntegerVectors(3, max_length=2)
            sage: IV == loads(dumps(IV))
            True
        """
        IntegerVectors_nkconstraints.__init__(self, n, -1, constraints)

    def __repr__(self):
        """
        EXAMPLES::

            sage: repr(IntegerVectors(3))
            'Integer vectors that sum to 3'
            sage: repr(IntegerVectors(3, max_length=2))
            'Integer vectors that sum to 3 with constraints: max_length=2'
        """
        if self.constraints:
            return "Integer vectors that sum to %s with constraints: %s"%(self.n,", ".join( ["%s=%s"%(key, self.constraints[key]) for key in sorted(self.constraints.keys())] ))
        else:
            return "Integer vectors that sum to %s"%(self.n,)

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: [0,3,0,1,2] in IntegerVectors(6)
            True
            sage: [0,3,0,1,2] in IntegerVectors(6, max_length=3)
            False
        """
        if self.constraints:
            return x in IntegerVectors_all() and misc.check_integer_list_constraints(x, singleton=True, **self.constraints)
        else:
            return x in IntegerVectors_all() and sum(x) == self.n

    def cardinality(self):
        """
        EXAMPLES::

            sage: IntegerVectors(3, max_length=2).cardinality()
            4
            sage: IntegerVectors(3).cardinality()
            +Infinity
        """
        if 'max_length' not in self.constraints:
            return infinity
        else:
            return self._CombinatorialClass__cardinality_from_iterator()

    def list(self):
        """
        EXAMPLES::

            sage: IntegerVectors(3, max_length=2).list()
            [[3], [2, 1], [1, 2], [0, 3]]
            sage: IntegerVectors(3).list()
            Traceback (most recent call last):
            ...
            NotImplementedError: infinite list
        """
        if 'max_length' not in self.constraints:
            raise NotImplementedError, "infinite list" # no list from infinite iter
        else:
            return list(self)


class IntegerVectors_nnondescents(CombinatorialClass):
    r"""
    The combinatorial class of integer vectors v graded by two
    parameters:

    - n: the sum of the parts of v

    - comp: the non descents composition of v

    In other words: the length of v equals c[1]+...+c[k], and v is
    decreasing in the consecutive blocs of length c[1], ..., c[k]

    Those are the integer vectors of sum n which are lexicographically
    maximal (for the natural left->right reading) in their orbit by the
    young subgroup S_c_1 x x S_c_k. In particular, they form a set
    of orbit representative of integer vectors w.r.t. this young
    subgroup.
    """
    def __init__(self, n, comp):
        """
        EXAMPLES::

            sage: IV = IntegerVectors(4, [2])
            sage: IV == loads(dumps(IV))
            True
        """
        self.n = n
        self.comp = comp

    def __repr__(self):
        """
        EXAMPLES::

            sage: IntegerVectors(4, [2]).__repr__()
            'Integer vectors of 4 with non-descents composition [2]'
        """
        return "Integer vectors of %s with non-descents composition %s"%(self.n, self.comp)

    def __iter__(self):
        """
         TESTS::

             sage: IntegerVectors(0, []).list()
             [[]]
             sage: IntegerVectors(5, []).list()
             []
             sage: IntegerVectors(0, [1]).list()
             [[0]]
             sage: IntegerVectors(4, [1]).list()
             [[4]]
             sage: IntegerVectors(4, [2]).list()
             [[4, 0], [3, 1], [2, 2]]
             sage: IntegerVectors(4, [2,2]).list()
              [[4, 0, 0, 0],
              [3, 0, 1, 0],
              [2, 0, 2, 0],
              [2, 0, 1, 1],
              [1, 0, 3, 0],
              [1, 0, 2, 1],
              [0, 0, 4, 0],
              [0, 0, 3, 1],
              [0, 0, 2, 2]]
             sage: IntegerVectors(5, [1,1,1]).list()
             [[5, 0, 0],
              [4, 1, 0],
              [4, 0, 1],
              [3, 2, 0],
              [3, 1, 1],
              [3, 0, 2],
              [2, 3, 0],
              [2, 2, 1],
              [2, 1, 2],
              [2, 0, 3],
              [1, 4, 0],
              [1, 3, 1],
              [1, 2, 2],
              [1, 1, 3],
              [1, 0, 4],
              [0, 5, 0],
              [0, 4, 1],
              [0, 3, 2],
              [0, 2, 3],
              [0, 1, 4],
              [0, 0, 5]]
             sage: IntegerVectors(0, [2,3]).list()
             [[0, 0, 0, 0, 0]]
         """
        for iv in IntegerVectors(self.n, len(self.comp)):
            blocks = [ IntegerVectors(iv[i], self.comp[i], max_slope=0).__iter__() for i in range(len(self.comp))]
            for parts in cartesian_product.CartesianProduct(*blocks):
                res = []
                for part in parts:
                    res += part
                yield res


