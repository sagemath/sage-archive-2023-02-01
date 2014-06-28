# -*- coding: utf-8 -*-
r"""
Mutually Orthogonal Latin Squares (MOLS)

A Latin square is an `n\times n` array filled with `n` different symbols, each
occurring exactly once in each row and exactly once in each column. For Sage's
methods related to Latin Squares, see the module
:mod:`sage.combinat.matrices.latin`.

This module gathers constructions of Mutually Orthogonal Latin Squares, which
are equivalent to Transversal Designs and specific Orthogonal Arrays.

For more information on MOLS, see the :wikipedia:`Wikipedia entry on MOLS
<Graeco-Latin_square#Mutually_orthogonal_Latin_squares>`.

The following table prints the maximum number of MOLS that Sage can build for
every order `n<300`, similarly to the `table of MOLS
<http://books.google.fr/books?id=S9FA9rq1BgoC&dq=handbook%20combinatorial%20designs%20MOLS%2010000&pg=PA176>`_
from the Handbook of Combinatorial Designs.

::

    sage: def MOLS_table(number_of_lines):
    ....:     print "     "+join(['%3s'%str(i) for i in range(20)])
    ....:     print "    "+"_"*80
    ....:     for i in range(20*15):
    ....:         if i%20==0:
    ....:             print "\n"+'%3s'%str(i)+"|",
    ....:         print '%3s'%str(designs.mutually_orthogonal_latin_squares(i,None,existence=True) if i>1 else "+oo"),
    sage: MOLS_table(15) # long time
           0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19
        ________________________________________________________________________________
    <BLANKLINE>
      0| +oo +oo   1   2   3   4   1   6   7   8   2  10   5  12   4   4  15  16   5  18
     20|   4   5   3  22   7  24   4  26   5  28   4  30  31   5   4   5   8  36   4   5
     40|   7  40   5  42   5   6   4  46   8  48   6   5   5  52   5   6   7   7   2  58
     60|   5  60   5   6  63   7   5  66   5   6   6  70   7  72   5   7   6   6   6  78
     80|   9  80   8  82   6   6   6   3   7  88   4   6   6   4   3   6   7  96   6   8
    100|   8 100   6 102   7   7   5 106   5 108   4   6   7 112   3   7   5   8   4   6
    120|   6 120   5   6   5 124   6 126 127   7   6 130   6   6   6   6   7 136   4 138
    140|   6   7   6  10  10   7   6   7   5 148   6 150   7   8   8   5   5 156   4   6
    160|   7   7   6 162   5   7   4 166   7 168   6   8   6 172   6   6  10   9   6 178
    180|   6 180   6   6   7   8   6  10   6   6   6 190   7 192   6   7   6 196   6 198
    200|   7   7   6   7   6   6   6   8  12  11  10 210   6   7   6   7   7   8   6  10
    220|   6  12   6 222   7   8   6 226   6 228   6   6   7 232   6   7   6   6   5 238
    240|   7 240   6 242   6   7   6  12   7   7   5 250   6  10   7   7 255 256   4   7
    260|   6   8   7 262   7   8   6  10   6 268   6 270  15  16   5  10  10 276   6   8
    280|   7 280   6 282   6  12   6   7  15 288   6   6   5 292   6   6   7  10  10  12

TODO:

* Look at [ColDin01]_.

REFERENCES:

.. [Stinson2004] Douglas R. Stinson,
  Combinatorial designs: construction and analysis,
  Springer, 2004.

.. [ColDin01] Charles Colbourn, Jeffrey Dinitz,
  Mutually orthogonal latin squares: a brief survey of constructions,
  Volume 95, Issues 1-2, Pages 9-48,
  Journal of Statistical Planning and Inference,
  Springer, 1 May 2001.

Functions
---------
"""
from sage.categories.sets_cat import EmptySetError
from sage.misc.unknown import Unknown

def are_mutually_orthogonal_latin_squares(l, verbose=False):
    r"""
    Check wether the list of matrices in ``l`` form mutually orthogonal latin
    squares.

    INPUT:

    - ``verbose`` - if ``True`` then print why the list of matrices provided are
      not mutually orthogonal latin squares

    EXAMPLES::

        sage: from sage.combinat.designs.latin_squares import are_mutually_orthogonal_latin_squares
        sage: m1 = matrix([[0,1,2],[2,0,1],[1,2,0]])
        sage: m2 = matrix([[0,1,2],[1,2,0],[2,0,1]])
        sage: m3 = matrix([[0,1,2],[2,0,1],[1,2,0]])
        sage: are_mutually_orthogonal_latin_squares([m1,m2])
        True
        sage: are_mutually_orthogonal_latin_squares([m1,m3])
        False
        sage: are_mutually_orthogonal_latin_squares([m2,m3])
        True
        sage: are_mutually_orthogonal_latin_squares([m1,m2,m3], verbose=True)
        Squares 0 and 2 are not orthogonal
        False

        sage: m = designs.mutually_orthogonal_latin_squares(8,7)
        sage: are_mutually_orthogonal_latin_squares(m)
        True

    TESTS:

    Not a latin square::

        sage: m1 = matrix([[0,1,0],[2,0,1],[1,2,0]])
        sage: m2 = matrix([[0,1,2],[1,2,0],[2,0,1]])
        sage: are_mutually_orthogonal_latin_squares([m1,m2], verbose=True)
        Matrix 0 is not row latin
        False
        sage: m1 = matrix([[0,1,2],[1,0,2],[1,2,0]])
        sage: are_mutually_orthogonal_latin_squares([m1,m2], verbose=True)
        Matrix 0 is not column latin
        False
        sage: m1 = matrix([[0,0,0],[1,1,1],[2,2,2]])
        sage: m2 = matrix([[0,1,2],[0,1,2],[0,1,2]])
        sage: are_mutually_orthogonal_latin_squares([m1,m2])
        False
    """

    if not l:
        raise ValueError("the list must be non empty")

    n = l[0].ncols()
    k = len(l)
    if any(M.ncols() != n or M.nrows() != n for M in l):
        if verbose:
            print "Not all matrices are square matrices of the same dimensions"
        return False

    # Check that all matrices are latin squares
    for i,M in enumerate(l):
        if any(len(set(R)) != n for R in M):
            if verbose:
                print "Matrix {} is not row latin".format(i)
            return False
        if any(len(set(R)) != n for R in zip(*M)):
            if verbose:
                print "Matrix {} is not column latin".format(i)
            return False

    from designs_pyx import is_orthogonal_array
    return is_orthogonal_array(zip(*[[x for R in M for x in R] for M in l]),k,n, verbose=verbose, terminology="MOLS")

def mutually_orthogonal_latin_squares(n,k, partitions = False, check = True, existence=False, who_asked=tuple()):
    r"""
    Returns `k` Mutually Orthogonal `n\times n` Latin Squares (MOLS).

    For more information on Latin Squares and MOLS, see
    :mod:`~sage.combinat.designs.latin_squares` or the :wikipedia:`Latin_square`,
    or even the
    :wikipedia:`Wikipedia entry on MOLS <Graeco-Latin_square#Mutually_orthogonal_Latin_squares>`.

    INPUT:

    - ``n`` (integer) -- size of the latin square.

    - ``k`` (integer) -- number of MOLS. If ``k=None`` it is set to the largest
      value available.

    - ``partition`` (boolean) -- a Latin Square can be seen as 3 partitions of
      the `n^2` cells of the array into `n` sets of size `n`, respectively :

      * The partition of rows
      * The partition of columns
      * The partition of number (cells numbered with 0, cells numbered with 1,
        ...)

      These partitions have the additional property that any two sets from
      different partitions intersect on exactly one element.

      When ``partition`` is set to ``True``, this function returns a list of `k+2`
      partitions satisfying this intersection property instead of the `k+2` MOLS
      (though the data is exactly the same in both cases).

    - ``existence`` (boolean) -- instead of building the design, returns:

        - ``True`` -- meaning that Sage knows how to build the design

        - ``Unknown`` -- meaning that Sage does not know how to build the
          design, but that the design may exist (see :mod:`sage.misc.unknown`).

        - ``False`` -- meaning that the design does not exist.

      .. NOTE::

          When ``k=None`` and ``existence=True`` the function returns an
          integer, i.e. the largest `k` such that we can build a `k` MOLS of
          order `n`.

    - ``check`` -- (boolean) Whether to check that output is correct before
      returning it. As this is expected to be useless (but we are cautious
      guys), you may want to disable it whenever you want speed. Set to
      ``True`` by default.

    - ``who_asked`` (internal use only) -- because of the equivalence between
      OA/TD/MOLS, each of the three constructors calls the others. We must keep
      track of who calls who in order to avoid infinite loops. ``who_asked`` is
      the tuple of the other functions that were called before this one.

    EXAMPLES::

        sage: designs.mutually_orthogonal_latin_squares(5,4)
        [
        [0 2 4 1 3]  [0 3 1 4 2]  [0 4 3 2 1]  [0 1 2 3 4]
        [4 1 3 0 2]  [3 1 4 2 0]  [2 1 0 4 3]  [4 0 1 2 3]
        [3 0 2 4 1]  [1 4 2 0 3]  [4 3 2 1 0]  [3 4 0 1 2]
        [2 4 1 3 0]  [4 2 0 3 1]  [1 0 4 3 2]  [2 3 4 0 1]
        [1 3 0 2 4], [2 0 3 1 4], [3 2 1 0 4], [1 2 3 4 0]
        ]

        sage: designs.mutually_orthogonal_latin_squares(7,3)
        [
        [0 2 4 6 1 3 5]  [0 3 6 2 5 1 4]  [0 4 1 5 2 6 3]
        [6 1 3 5 0 2 4]  [5 1 4 0 3 6 2]  [4 1 5 2 6 3 0]
        [5 0 2 4 6 1 3]  [3 6 2 5 1 4 0]  [1 5 2 6 3 0 4]
        [4 6 1 3 5 0 2]  [1 4 0 3 6 2 5]  [5 2 6 3 0 4 1]
        [3 5 0 2 4 6 1]  [6 2 5 1 4 0 3]  [2 6 3 0 4 1 5]
        [2 4 6 1 3 5 0]  [4 0 3 6 2 5 1]  [6 3 0 4 1 5 2]
        [1 3 5 0 2 4 6], [2 5 1 4 0 3 6], [3 0 4 1 5 2 6]
        ]

        sage: designs.mutually_orthogonal_latin_squares(5,2,partitions=True)
        [[[0, 1, 2, 3, 4],
          [5, 6, 7, 8, 9],
          [10, 11, 12, 13, 14],
          [15, 16, 17, 18, 19],
          [20, 21, 22, 23, 24]],
         [[0, 5, 10, 15, 20],
          [1, 6, 11, 16, 21],
          [2, 7, 12, 17, 22],
          [3, 8, 13, 18, 23],
          [4, 9, 14, 19, 24]],
         [[0, 8, 11, 19, 22],
          [3, 6, 14, 17, 20],
          [1, 9, 12, 15, 23],
          [4, 7, 10, 18, 21],
          [2, 5, 13, 16, 24]],
         [[0, 9, 13, 17, 21],
          [2, 6, 10, 19, 23],
          [4, 8, 12, 16, 20],
          [1, 5, 14, 18, 22],
          [3, 7, 11, 15, 24]]]

    What is the maximum number of MOLS of size 8 that Sage knows how to build?::

        sage: designs.mutually_orthogonal_latin_squares(8,None,existence=True)
        7

    If you only want to know if Sage is able to build a given set of MOLS, just
    set the argument ``existence`` to ``True``::

        sage: designs.mutually_orthogonal_latin_squares(5, 5, existence=True)
        False
        sage: designs.mutually_orthogonal_latin_squares(6, 4, existence=True)
        Unknown

    If you ask for such a MOLS then you will respecively get an informative
    ``EmptySetError`` or ``NotImplementedError``::

        sage: designs.mutually_orthogonal_latin_squares(5, 5)
        Traceback (most recent call last):
        ...
        EmptySetError: There exist at most n-1 MOLS of size n if n>=2.
        sage: designs.mutually_orthogonal_latin_squares(6, 4)
        Traceback (most recent call last):
        ...
        NotImplementedError: I don't know how to build 4 MOLS of order 6

    TESTS:

    The special case `n=1`::

        sage: designs.mutually_orthogonal_latin_squares(1, 3)
        [[0], [0], [0]]
        sage: designs.mutually_orthogonal_latin_squares(1, None, existence=True)
        +Infinity
        sage: designs.mutually_orthogonal_latin_squares(1, None)
        Traceback (most recent call last):
        ...
        ValueError: there are no bound on k when n=1.
        sage: designs.mutually_orthogonal_latin_squares(10,2,existence=True)
        True
        sage: designs.mutually_orthogonal_latin_squares(10,2)
        [
        [1 8 9 0 2 4 6 3 5 7]  [1 7 6 5 0 9 8 2 3 4]
        [7 2 8 9 0 3 5 4 6 1]  [8 2 1 7 6 0 9 3 4 5]
        [6 1 3 8 9 0 4 5 7 2]  [9 8 3 2 1 7 0 4 5 6]
        [5 7 2 4 8 9 0 6 1 3]  [0 9 8 4 3 2 1 5 6 7]
        [0 6 1 3 5 8 9 7 2 4]  [2 0 9 8 5 4 3 6 7 1]
        [9 0 7 2 4 6 8 1 3 5]  [4 3 0 9 8 6 5 7 1 2]
        [8 9 0 1 3 5 7 2 4 6]  [6 5 4 0 9 8 7 1 2 3]
        [2 3 4 5 6 7 1 8 9 0]  [3 4 5 6 7 1 2 8 0 9]
        [3 4 5 6 7 1 2 0 8 9]  [5 6 7 1 2 3 4 0 9 8]
        [4 5 6 7 1 2 3 9 0 8], [7 1 2 3 4 5 6 9 8 0]
        ]
    """
    from sage.combinat.designs.orthogonal_arrays import orthogonal_array
    from sage.matrix.constructor import Matrix
    from sage.rings.arith import factor
    from database import MOLS_constructions

    # Is k is None we find the largest available
    if k is None:
        if n == 1:
            if existence:
                from sage.rings.infinity import Infinity
                return Infinity
            raise ValueError("there are no bound on k when n=1.")

        k = orthogonal_array(None,n,existence=True) - 2
        if existence:
            return k

    if n == 1:
        if existence:
            return True
        matrices = [Matrix([[0]])]*k

    elif k >= n:
        if existence:
            return False
        raise EmptySetError("There exist at most n-1 MOLS of size n if n>=2.")

    elif n in MOLS_constructions and k <= MOLS_constructions[n][0]:
        if existence:
            return True
        _, construction = MOLS_constructions[n]

        matrices = construction()[:k]

    elif (orthogonal_array not in who_asked and
        orthogonal_array(k+2,n,existence=True,who_asked = who_asked+(mutually_orthogonal_latin_squares,)) is not Unknown):

        # Forwarding non-existence results
        if orthogonal_array(k+2,n,existence=True,who_asked = who_asked+(mutually_orthogonal_latin_squares,)):
            if existence:
                return True
        else:
            if existence:
                return False
            raise EmptySetError("There does not exist {} MOLS of order {}!".format(k,n))

        OA = orthogonal_array(k+2,n,check=False, who_asked = who_asked+(mutually_orthogonal_latin_squares,))
        OA.sort() # make sure that the first two columns are "11, 12, ..., 1n, 21, 22, ..."

        # We first define matrices as lists of n^2 values
        matrices = [[] for _ in range(k)]
        for L in OA:
            for i in range(2,k+2):
                matrices[i-2].append(L[i])

        # The real matrices
        matrices = [[M[i*n:(i+1)*n] for i in range(n)] for M in matrices]
        matrices = [Matrix(M) for M in matrices]

    else:
        if existence:
            return Unknown
        raise NotImplementedError("I don't know how to build {} MOLS of order {}".format(k,n))

    if check:
        assert are_mutually_orthogonal_latin_squares(matrices)

    # partitions have been requested but have not been computed yet
    if partitions is True:
        partitions = [[[i*n+j for j in range(n)] for i in range(n)],
                      [[j*n+i for j in range(n)] for i in range(n)]]
        for m in matrices:
            partition = [[] for i in range(n)]
            for i in range(n):
                for j in range(n):
                    partition[m[i,j]].append(i*n+j)
            partitions.append(partition)

    if partitions:
        return partitions
    else:
        return matrices

def latin_square_product(M,N,*others):
    r"""
    Returns the product of two (or more) latin squares.

    Given two Latin Squares `M,N` of respective sizes `m,n`, the direct product
    `M\times N` of size `mn` is defined by `(M\times
    N)((i_1,i_2),(j_1,j_2))=(M(i_1,j_1),N(i_2,j_2))` where `i_1,j_1\in [m],
    i_2,j_2\in [n]`

    Each pair of values `(i,j)\in [m]\times [n]` is then relabeled to `in+j`.

    This is Lemma 6.25 of [Stinson2004]_.

    INPUT:

    An arbitrary number of latin squares (greater than 2).

    EXAMPLES::

        sage: from sage.combinat.designs.latin_squares import latin_square_product
        sage: m=designs.mutually_orthogonal_latin_squares(4,3)[0]
        sage: latin_square_product(m,m,m)
        64 x 64 sparse matrix over Integer Ring (use the '.str()' method to see the entries)
    """
    from sage.matrix.constructor import Matrix
    m = M.nrows()
    n = N.nrows()

    D = {((i,j),(ii,jj)):(M[i,ii],N[j,jj])
         for i in range(m)
         for ii in range(m)
         for j in range(n)
         for jj in range(n)}

    L = lambda i_j: i_j[0] * n + i_j[1]
    D = {(L(c[0]),L(c[1])): L(v) for c,v in D.iteritems()}
    P = Matrix(D)

    if others:
        return latin_square_product(P, others[0],*others[1:])
    else:
        return P
