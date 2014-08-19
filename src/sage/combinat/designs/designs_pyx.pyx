r"""
Cython functions for combinatorial designs

This module implements the design methods that need to be somewhat efficient.

Functions
---------
"""
include "sage/misc/bitset.pxi"

def is_orthogonal_array(OA, int k, int n, int t=2, verbose=False, terminology="OA"):
    r"""
    Check that the integer matrix `OA` is an `OA(k,n,t)`.

    See :func:`~sage.combinat.designs.orthogonal_arrays.orthogonal_array`
    for a definition.

    INPUT:

    - ``OA`` -- the Orthogonal Array to be tested

    - ``k,n,t`` (integers) -- only implemented for `t=2`.

    - ``verbose`` (boolean) -- whether to display some information when ``OA``
      is not an orthogonal array `OA(k,n)`.

    - ``terminology`` (string) -- how to phrase the information when ``verbose =
      True``. Possible values are `"OA"`, `"MOLS"`.

    EXAMPLES::

        sage: from sage.combinat.designs.designs_pyx import is_orthogonal_array
        sage: OA = designs.orthogonal_array(8,9)
        sage: is_orthogonal_array(OA,8,9)
        True
        sage: is_orthogonal_array(OA,8,10)
        False
        sage: OA[4][3] = 1
        sage: is_orthogonal_array(OA,8,9)
        False
        sage: is_orthogonal_array(OA,8,9,verbose=True)
        Columns 0 and 3 are not orthogonal
        False
        sage: is_orthogonal_array(OA,8,9,verbose=True,terminology="MOLS")
        Squares 0 and 3 are not orthogonal
        False

    TESTS::

        sage: is_orthogonal_array(OA,8,9,t=3)
        Traceback (most recent call last):
        ...
        NotImplementedError: only implemented for t=2
        sage: is_orthogonal_array([[3]*8],8,9,verbose=True)
        The number of rows is 1 instead of 9^2=81
        False
        sage: is_orthogonal_array([[3]*8],8,9,verbose=True,terminology="MOLS")
        All squares do not have dimension n^2=9^2
        False
        sage: is_orthogonal_array([[3]*7],8,9,verbose=True)
        Some row does not have length 8
        False
        sage: is_orthogonal_array([[3]*7],8,9,verbose=True,terminology="MOLS")
        The number of squares is not 6
        False

    Up to relabelling, there is a unique `OA(3,2)`. So their number is just the
    cardinality of the relabeling group which is `S_2^3 \times S_3` and has
    cardinality `48`::

        sage: from itertools import product
        sage: n = 0
        sage: for a in product(product((0,1), repeat=3), repeat=4):
        ....:     if is_orthogonal_array(a,3,2):
        ....:          n += 1
        sage: n
        48
    """
    cdef int n2 = n*n
    cdef int x

    if t != 2:
        raise NotImplementedError("only implemented for t=2")

    for R in OA:
        if len(R) != k:
            if verbose:
                print {"OA"   : "Some row does not have length "+str(k),
                       "MOLS" : "The number of squares is not "+str(k-2)}[terminology]
            return False

    if len(OA) != n2:
        if verbose:
            print {"OA"   : "The number of rows is {} instead of {}^2={}".format(len(OA),n,n2),
                   "MOLS" : "All squares do not have dimension n^2={}^2".format(n)}[terminology]
        return False

    if n == 0:
        return True

    cdef int i,j,l

    # A copy of OA
    cdef unsigned short * OAc = <unsigned short *> sage_malloc(k*n2*sizeof(unsigned short))

    cdef unsigned short * C1
    cdef unsigned short * C2

    # failed malloc ?
    if OAc is NULL:
        raise MemoryError

    # Filling OAc
    for i,R in enumerate(OA):
        for j,x in enumerate(R):
            if x < 0 or x >= n:
                if verbose:
                    print {"OA"   : "{} is not in the interval [0..{}]".format(x,n-1),
                           "MOLS" : "Entry {} was expected to be in the interval [0..{}]".format(x,n-1)}[terminology]
                sage_free(OAc)
                return False
            OAc[j*n2+i] = x

    # A bitset to keep track of pairs of values
    cdef bitset_t seen
    bitset_init(seen, n2)

    for i in range(k): # For any column C1
        C1 = OAc+i*n2
        for j in range(i+1,k): # For any column C2 > C1
            C2 = OAc+j*n2
            bitset_set_first_n(seen, 0) # No pair has been seen yet
            for l in range(n2):
                bitset_add(seen,n*C1[l]+C2[l])

            if bitset_len(seen) != n2: # Have we seen all pairs ?
                sage_free(OAc)
                bitset_free(seen)
                if verbose:
                    print {"OA"   : "Columns {} and {} are not orthogonal".format(i,j),
                           "MOLS" : "Squares {} and {} are not orthogonal".format(i,j)}[terminology]
                return False

    sage_free(OAc)
    bitset_free(seen)
    return True

def is_group_divisible_design(groups,blocks,v,G=None,K=None,lambd=1,verbose=False):
    r"""
    Checks that input is a Group Divisible Design on `\{0,...,v-1\}`

    For more information on Group Divisible Designs, see
    :class:`~sage.combinat.designs.incidence_structure.GroupDivisibleDesign`.

    INPUT:

    - ``groups`` -- a partition of `X`

    - ``blocks`` -- collection of blocks

    - ``v`` (integers) -- size of the ground set assumed to be `X=\{0,...,v-1\}`.

    - ``G`` -- list of integers of which the sizes of the groups must be
      elements. Set to ``None`` (automatic guess) by default.

    - ``K`` -- list of integers of which the sizes of the blocks must be
      elements. Set to ``None`` (automatic guess) by default.

    - ``lambd`` -- value of `\lambda`. Set to `1` by default.

    - ``verbose`` (boolean) -- whether to display some information when the
      design is not a GDD.

    EXAMPLES::

        sage: from sage.combinat.designs.designs_pyx import is_group_divisible_design
        sage: TD = designs.transversal_design(4,10)
        sage: groups = [range(i*10,(i+1)*10) for i in range(4)]
        sage: is_group_divisible_design(groups,TD,40,lambd=1)
        True

    TESTS::

        sage: TD = designs.transversal_design(4,10)
        sage: groups = [range(i*10,(i+1)*10) for i in range(4)]
        sage: is_group_divisible_design(groups,TD,40,lambd=2,verbose=True)
        the pair (0,10) has been seen 1 times but lambda=2
        False
        sage: is_group_divisible_design([[1,2],[3,4]],[[1,2]],40,lambd=1,verbose=True)
        groups is not a partition of [0,...,39]
        False
        sage: is_group_divisible_design([range(40)],[[1,2]],40,lambd=1,verbose=True)
        the pair (1,2) belongs to a group but appears in some block
        False
        sage: is_group_divisible_design([range(40)],[[2,2]],40,lambd=1,verbose=True)
        The following block has repeated elements: [2, 2]
        False
        sage: is_group_divisible_design([range(40)],[["e",2]],40,lambd=1,verbose=True)
        e does not belong to [0,...,39]
        False
        sage: is_group_divisible_design([range(40)],[["e",2]],40,G=[5],lambd=1,verbose=True)
        a group has size 40 while G=[5]
        False
        sage: is_group_divisible_design([range(40)],[["e",2]],40,K=[1],lambd=1,verbose=True)
        a block has size 2 while K=[1]
        False
    """
    cdef int n = v
    cdef int i,ii,j,jj,s,isok
    cdef int l = lambd

    if v < 0 or lambd < 0:
        if verbose:
            print "v={} and lambda={} must be non-negative integers".format(v,l)
        return False

    # Group sizes are element of G
    if G is not None:
        G = set(G)
        for g in groups:
            if not len(g) in G:
                if verbose:
                    print "a group has size {} while G={}".format(len(g),list(G))
                return False

    # Block sizes are element of K
    if K is not None:
        K = set(K)
        for b in blocks:
            if not len(b) in K:
                if verbose:
                    print "a block has size {} while K={}".format(len(b),list(K))
                return False

    # Check that "groups" consists of disjoints sets whose union has length n
    if sum(len(g) for g in groups) != n or len(set().union(*groups)) != n:
        if verbose:
            print "groups is not a partition of [0,...,{}]".format(n-1)
        return False

    # Checks that the blocks are indeed sets and do not repeat elements
    for b in blocks:
        if len(b) != len(set(b)):
            if verbose:
                print "The following block has repeated elements: {}".format(b)
            return False

    # Check that the groups/blocks belong to [0,...,n-1]
    from itertools import chain
    for b in chain(groups,blocks):
        for x in b:
            try:
                i = x
            except TypeError:
                i = -1
            if i < 0 or i >= n:
                if verbose:
                    print "{} does not belong to [0,...,{}]".format(x,n-1)
                return False

    cdef unsigned short * matrix = <unsigned short *> sage_calloc(n*n,sizeof(unsigned short))
    if matrix is NULL:
        raise MemoryError

    # Counts the number of occurrences of each pair of points
    for b in blocks:
        s = len(b)
        for i in range(s):
            ii = b[i]
            for j in range(i+1,s):
                jj = b[j]
                matrix[ii*n+jj] += 1
                matrix[jj*n+ii] += 1

    # Checks that two points of the same group were never covered
    for g in groups:
        s = len(g)
        for i in range(s):
            ii = g[i]
            for j in range(i+1,s):
                jj = g[j]
                if matrix[ii*n+jj] != 0:
                    if verbose:
                        print "the pair ({},{}) belongs to a group but appears in some block".format(ii,jj)
                    sage_free(matrix)
                    return False

                # We fill the entries with what is expected by the next loop
                matrix[ii*n+jj] = l
                matrix[jj*n+ii] = l

    # Checking that what should be equal to lambda IS equal to lambda
    for i in range(n):
        for j in range(i+1,n):
            if i != j and matrix[i*n+j] != l:
                if verbose:
                    print "the pair ({},{}) has been seen {} times but lambda={}".format(i,j,matrix[i*n+j],l)
                sage_free(matrix)
                return False

    sage_free(matrix)

    return True

def is_pairwise_balanced_design(blocks,v,K=None,lambd=1,verbose=False):
    r"""
    Checks that input is a Pairwise Balanced Design (PBD) on `\{0,...,v-1\}`

    For more information on Pairwise Balanced Designs (PBD), see
    :class:`~sage.combinat.designs.bibd.PairwiseBalancedDesign`.

    INPUT:

    - ``blocks`` -- collection of blocks

    - ``v`` (integers) -- size of the ground set assumed to be `X=\{0,...,v-1\}`.

    - ``K`` -- list of integers of which the sizes of the blocks must be
      elements. Set to ``None`` (automatic guess) by default.

    - ``lambd`` -- value of `\lambda`. Set to `1` by default.

    - ``verbose`` (boolean) -- whether to display some information when the
      design is not a PBD.

    EXAMPLES::

        sage: from sage.combinat.designs.designs_pyx import is_pairwise_balanced_design
        sage: sts = designs.steiner_triple_system(9)
        sage: is_pairwise_balanced_design(sts,9,[3],1)
        True
        sage: TD = designs.transversal_design(4,10).blocks()
        sage: groups = [range(i*10,(i+1)*10) for i in range(4)]
        sage: is_pairwise_balanced_design(TD+groups,40,[4,10],1,verbose=True)
        True

    TESTS::

        sage: from sage.combinat.designs.designs_pyx import is_pairwise_balanced_design
        sage: is_pairwise_balanced_design(TD+groups,40,[4,10],2,verbose=True)
        the pair (0,1) has been seen 1 times but lambda=2
        False
        sage: is_pairwise_balanced_design(TD+groups,40,[10],1,verbose=True)
        a block has size 4 while K=[10]
        False
        sage: is_pairwise_balanced_design([[2,2]],40,[2],1,verbose=True)
        The following block has repeated elements: [2, 2]
        False
        sage: is_pairwise_balanced_design([["e",2]],40,[2],1,verbose=True)
        e does not belong to [0,...,39]
        False
    """
    return is_group_divisible_design([[i] for i in range(v)],
                                     blocks,
                                     v,
                                     K=K,
                                     lambd=lambd,
                                     verbose=verbose)
