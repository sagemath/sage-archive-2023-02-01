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

TODO:

* Implement Wilson's construction (page 146 of [Stinson2004]_)
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

def mutually_orthogonal_latin_squares(n,k=None, partitions = False):
    r"""
    Returns `k` Mutually Orthogonal `n\times n` Latin Squares (MOLS).

    For more information on Latin Squares and MOLS, see
    :mod:`~sage.combinat.designs.latin_squares` or the :wikipedia:`Latin_square`,
    or even the
    :wikipedia:`Wikipedia entry on MOLS <Graeco-Latin_square#Mutually_orthogonal_Latin_squares>`.

    INPUT:

    - ``n`` (integer) -- size of the latin square.

    - ``k`` (integer) -- returns `k` MOLS. If set to ``None`` (default), returns
      the maximum number of MOLS that Sage can build.

      .. WARNING::

          This has no reason to be the maximum number of `n\times n` MOLS, just
          the best Sage can do !

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

    EXAMPLES::

        sage: designs.mutually_orthogonal_latin_squares(5)
        [
        [0 1 2 3 4]  [0 1 2 3 4]  [0 1 2 3 4]  [0 1 2 3 4]
        [3 0 1 4 2]  [4 3 0 2 1]  [1 2 4 0 3]  [2 4 3 1 0]
        [4 3 0 2 1]  [1 2 4 0 3]  [2 4 3 1 0]  [3 0 1 4 2]
        [1 2 4 0 3]  [2 4 3 1 0]  [3 0 1 4 2]  [4 3 0 2 1]
        [2 4 3 1 0], [3 0 1 4 2], [4 3 0 2 1], [1 2 4 0 3]
        ]
        sage: designs.mutually_orthogonal_latin_squares(7,3)
        [
        [0 1 2 3 4 5 6]  [0 1 2 3 4 5 6]  [0 1 2 3 4 5 6]
        [4 0 3 1 6 2 5]  [5 6 0 4 2 1 3]  [6 4 1 0 5 3 2]
        [5 6 0 4 2 1 3]  [6 4 1 0 5 3 2]  [1 3 5 2 0 6 4]
        [6 4 1 0 5 3 2]  [1 3 5 2 0 6 4]  [2 5 4 6 3 0 1]
        [1 3 5 2 0 6 4]  [2 5 4 6 3 0 1]  [3 2 6 5 1 4 0]
        [2 5 4 6 3 0 1]  [3 2 6 5 1 4 0]  [4 0 3 1 6 2 5]
        [3 2 6 5 1 4 0], [4 0 3 1 6 2 5], [5 6 0 4 2 1 3]
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
        [[0, 6, 12, 18, 24],
          [1, 7, 14, 15, 23],
          [2, 9, 13, 16, 20],
          [3, 5, 11, 19, 22],
          [4, 8, 10, 17, 21]],
        [[0, 7, 13, 19, 21],
          [1, 9, 10, 18, 22],
          [2, 8, 11, 15, 24],
          [3, 6, 14, 17, 20],
          [4, 5, 12, 16, 23]]]

    TESTS::

        sage: designs.mutually_orthogonal_latin_squares(5,5)
        Traceback (most recent call last):
        ...
        ValueError: There exist at most n-1 MOLS of size n.
    """
    from sage.rings.finite_rings.constructor import FiniteField
    from sage.combinat.designs.block_design import AffineGeometryDesign
    from sage.rings.arith import is_prime_power
    from sage.matrix.constructor import Matrix
    from sage.rings.arith import factor

    if k is not None and k >= n:
        raise ValueError("There exist at most n-1 MOLS of size n.")

    if is_prime_power(n):
        if k is None:
            k = n-1
        # Section 6.4.1 of [Stinson2004]
        Fp = FiniteField(n,'x')
        B = AffineGeometryDesign(2,1,Fp).blocks()
        parallel_classes = [[] for _ in range(k+2)]
        for b in B:
            for p in parallel_classes:
                if (not p) or all(i not in p[0] for i in b):
                    p.append(b)
                    break

        if partitions:
            return parallel_classes

        coord = {v:i
                 for i,L in enumerate(parallel_classes[0]) for v in L}
        coord = {v:(coord[v],i)
                 for i,L in enumerate(parallel_classes[1]) for v in L}

        matrices = []
        for P in parallel_classes[2:]:
            matrices.append(Matrix({coord[v]:i for i,L in enumerate(P) for v in L }))
        return matrices
    else:
        # Theorem 6.33 of [Stinson2004], MacNeish's theorem.
        subcases = [p**i for p,i in factor(n)]
        s = min(subcases)-1
        if k is None:
            k = s
        elif k > s:
            raise NotImplementedError("I don't know how to build these MOLS.")
        subcalls = [mutually_orthogonal_latin_squares(p,k) for p in subcases]
        matrices = [latin_square_product(*[sc[i] for sc in subcalls])
                    for i in range(k)]

        if partitions:
            partitions = [[[i*n+j for j in range(n)] for i in range(n)],
                          [[j*n+i for j in range(n)] for i in range(n)]]
            for m in matrices:
                partition = [[] for i in range(n)]
                for i in range(n):
                    for j in range(n):
                        partition[m[i,j]].append(i*n+j)
                partitions.append(partition)

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
        sage: m=designs.mutually_orthogonal_latin_squares(4)[0]
        sage: latin_square_product(m,m,m)
        64 x 64 sparse matrix over Integer Ring
    """
    from sage.matrix.constructor import Matrix
    m = M.nrows()
    n = N.nrows()

    D = {((i,j),(ii,jj)):(M[i,ii],N[j,jj])
         for i in range(m)
         for ii in range(m)
         for j in range(n)
         for jj in range(n)}

    L = lambda (i,j):i*n+j
    D = {(L(c[0]),L(c[1])): L(v) for c,v in D.iteritems()}
    P = Matrix(D)

    if others:
        return latin_square_product(P, others[0],*others[1:])
    else:
        return P
