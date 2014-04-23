"""
Orthogonal arrays

This module gathers anything related to orthogonal arrays, and, incidentally,
to transversal designs.

Functions
---------
"""

def transversal_design(k,n,t=2,check=True):
    r"""
    Return a transversal design of parameters `k,n`.

    A transversal design of parameters `n, k` is a collection `\mathcal{S}`
    of `k`-subsets of `V = V_1 \sqcup \cdots \sqcup V_k`
    (where `|V_i| = n` for all `i`) such that:

    * Any element `S \in \mathcal{S}` intersects each set `V_i` on exactly one
      element.

    * Any two elements `v_i \in V_i, v_j\in V_j` with `i \neq j` belong to
      exactly one element of `\mathcal{S}`.

    For more information on transversal designs, see
    http://mathworld.wolfram.com/TransversalDesign.html.

    INPUT:

    - `n,k,t` -- integers

    - ``check`` -- (boolean) Whether to check that output is correct before
      returning it. As this is expected to be useless (but we are cautious
      guys), you may want to disable it whenever you want speed. Set to
      ``True`` by default.

    EXAMPLES::

        sage: designs.transversal_design(5,5)
        [[0, 5, 10, 15, 20], [0, 6, 12, 18, 24], [0, 7, 14, 16, 23],
         [0, 8, 11, 19, 22], [0, 9, 13, 17, 21], [1, 6, 11, 16, 21],
         [1, 7, 13, 19, 20], [1, 8, 10, 17, 24], [1, 9, 12, 15, 23],
         [1, 5, 14, 18, 22], [2, 7, 12, 17, 22], [2, 8, 14, 15, 21],
         [2, 9, 11, 18, 20], [2, 5, 13, 16, 24], [2, 6, 10, 19, 23],
         [3, 8, 13, 18, 23], [3, 9, 10, 16, 22], [3, 5, 12, 19, 21],
         [3, 6, 14, 17, 20], [3, 7, 11, 15, 24], [4, 9, 14, 19, 24],
         [4, 5, 11, 17, 23], [4, 6, 13, 15, 22], [4, 7, 10, 18, 21],
         [4, 8, 12, 16, 20]]
    """
    # Section 6.6
    OA = orthogonal_array(k,n)
    TD = [[i*n+c for i,c in enumerate(l)] for l in OA]

    if check:
        assert is_transversal_design(TD,k,n)

    return TD

def is_transversal_design(B,k,n):
    r"""
    Check that a given set of blocks ``B`` is a transversal design.

    See :func:`~sage.combinat.designs.orthogonal_arrays.transversal_design`
    for a definition.

    INPUT:

    - ``B`` -- the list of blocks
    - ``k, n`` -- integers

    .. NOTE::

        The tranversal design must have `\{0, \ldots, kn-1\}` as a ground set,
        partitioned as `k` sets of size `n`: `\{0, \ldots, k-1\} \sqcup
        \{k, \ldots, 2k-1\} \sqcup \cdots \sqcup \{k(n-1), \ldots, kn-1\}`.

    EXAMPLES::

        sage: TD = designs.transversal_design(5, 5, check=True) # indirect doctest
        sage: from sage.combinat.designs.orthogonal_arrays import is_transversal_design
        sage: is_transversal_design(TD, 5, 5)
        True
        sage: is_transversal_design(TD, 4, 4)
        False
    """
    from sage.graphs.generators.basic import CompleteGraph
    from itertools import combinations
    g = k*CompleteGraph(n)
    m = g.size()
    for X in B:
        if len(X) != k:
            return False
        g.add_edges(list(combinations(X,2)))
        if g.size() != m+(len(X)*(len(X)-1))/2:
            return False
        m = g.size()

    return g.is_clique()

def orthogonal_array(k,n,t=2,check=True):
    r"""
    Return an orthogonal array of parameters `k,n,t`.

    INPUT:

    - ``k`` -- (integer) number of columns

    - ``n`` -- (integer) number of symbols

    - ``t`` -- (integer; default: 2) only ``t=2`` is available at the moment

    - ``check`` -- (boolean) Whether to check that output is correct before
      returning it. As this is expected to be useless (but we are cautious
      guys), you may want to disable it whenever you want speed. Set to
      ``True`` by default.

    For more information on orthogonal arrays, see
    :wikipedia:`Orthogonal_array`.

    .. NOTE::

        This method implements theorems from [Stinson2004]_. See the code's
        documentation for details.

    .. TODO::

        Implement Wilson's construction. See page 146 of [Stinson2004]_.

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays import orthogonal_array
        sage: orthogonal_array(5,5)
        [[0, 0, 0, 0, 0], [0, 1, 2, 3, 4], [0, 2, 4, 1, 3], [0, 3, 1, 4, 2],
         [0, 4, 3, 2, 1], [1, 1, 1, 1, 1], [1, 2, 3, 4, 0], [1, 3, 0, 2, 4],
         [1, 4, 2, 0, 3], [1, 0, 4, 3, 2], [2, 2, 2, 2, 2], [2, 3, 4, 0, 1],
         [2, 4, 1, 3, 0], [2, 0, 3, 1, 4], [2, 1, 0, 4, 3], [3, 3, 3, 3, 3],
         [3, 4, 0, 1, 2], [3, 0, 2, 4, 1], [3, 1, 4, 2, 0], [3, 2, 1, 0, 4],
         [4, 4, 4, 4, 4], [4, 0, 1, 2, 3], [4, 1, 3, 0, 2], [4, 2, 0, 3, 1],
         [4, 3, 2, 1, 0]]

    TESTS::

        sage: designs.orthogonal_array(3,2)
        [[0, 1, 0], [0, 0, 1], [1, 0, 0], [1, 1, 1]]
    """
    from sage.rings.arith import is_prime_power
    from sage.rings.finite_rings.constructor import FiniteField
    OA = None

    if t != 2:
        raise NotImplementedError("only implemented for t=2")

    if k < 2:
        raise ValueError("undefined for k less than 2")

    if k == t:
        from itertools import product
        OA = map(list, product(range(n), repeat=k))

    # Theorem 6.39 from [Stinson2004]
    elif 2 <= k and k <= n and is_prime_power(n):
        M = []
        Fp = FiniteField(n,'x')
        vv = list(Fp)[:k]
        relabel = {x:i for i,x in enumerate(Fp)}
        for i in Fp:
            for j in Fp:
                M.append([relabel[i+j*v] for v in vv])

        OA = M

    # Theorem 6.40 from [Stinson2004]
    elif k == n+1 and is_prime_power(n):
        if n == 2:
            OA = [[0,1,0],[0,0,1],[1,0,0],[1,1,1]]
        else:
            M = orthogonal_array(n,n, check=False)
            for i,l in enumerate(M):
                l.append(i%n)
            OA = M

    # Section 6.5.1 from [Stinson2004]
    if OA is None:
        from latin_squares import mutually_orthogonal_latin_squares
        try:
            mols = mutually_orthogonal_latin_squares(n,k-2)
        except ValueError:
            mols = None
        if mols:
            OA = [[i,j]+[m[i,j] for m in mols]
                  for i in range(n) for j in range(n)]

    if OA is None:
        raise NotImplementedError("I don't know how to build this orthogonal array!")
    if check:
        assert is_orthogonal_array(OA,k,n,t)

    return OA

def is_orthogonal_array(M,k,n,t):
    r"""
    Check that the integer matrix `M` is an `OA(k,n,t)`.

    See :func:`~sage.combinat.designs.orthogonal_arrays.orthogonal_array`
    for a definition.

    INPUT:

    - ``M`` -- an integer matrix of size `k^t \times n`

    - ``k, n, t`` -- integers

    EXAMPLES::

        sage: OA = designs.orthogonal_array(5,9,check=True) # indirect doctest
        sage: from sage.combinat.designs.orthogonal_arrays import is_orthogonal_array
        sage: is_orthogonal_array(OA, 5, 9, 2)
        True
        sage: is_orthogonal_array(OA, 4, 5, 2)
        False
    """
    if t != 2:
        raise NotImplementedError("only implemented for t=2")

    if not all(len(l) == k for l in M):
        return False

    from itertools import combinations
    for S in combinations(range(k),2):
        fs = frozenset([tuple([l[i] for i in S]) for l in M])
        if len(fs) != n**2:
            return False

    return True

