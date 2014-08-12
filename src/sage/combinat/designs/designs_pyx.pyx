r"""
Cython functions for combinatorial designs

This module implements the design methods that need to be somewhat efficient.

Functions
---------
"""
include "sage/misc/bitset.pxi"

from libc.string cimport memset

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

def is_difference_matrix(G,k,M,lmbda=1,verbose=False):
    r"""
    Test if `M` is a `(G,k,\lambda)`-difference matrix.

    A matrix `M` is a `(G,k,\lambda)`-difference matrix if its entries are
    element of `G`, and if for any two rows `R,R'` of `M` and `x\in G` there
    are exactly `\lambda` values `i` such that `R_i-R'_i=x`.

    INPUT:

    - ``G`` -- a group

    - ``k`` -- integer

    - ``M`` -- a matrix with entries from ``G``

    - ``lmbda`` (integer) -- set to `1` by default.

    - ``verbose`` (boolean) -- whether to print some information when the answer
      is ``False``.

    EXAMPLES::

        sage: from sage.combinat.designs.designs_pyx import is_difference_matrix
        sage: q = 3**3
        sage: F = GF(q,'x')
        sage: M = [[x*y for y in F] for x in F]
        sage: is_difference_matrix(F,q,M,verbose=1)
        True

    Bad input::

        sage: M.append([None]*3**3)
        sage: is_difference_matrix(F,q,M,verbose=1)
        The matrix has 28 columns and lambda.|G|=1.27=27
        False
        sage: _=M.pop()
        sage: M[0].append(1)
        sage: is_difference_matrix(F,q,M,verbose=1)
        Rows 0 and 1 do not have the same length
        False
        sage: _= M[0].pop(-1)
        sage: M[-1] = [0]*3**3
        sage: is_difference_matrix(F,q,M,verbose=1)
        Rows 0 and 26 do not generate all elements of G exactly lambda(=1)
        times. The element 0 appeared 27 times as a difference.
        False
    """
    from difference_family import group_law

    assert k>=2
    assert lmbda >=1

    cdef int G_card = G.cardinality()
    cdef int i,j,ii
    cdef int K = k
    cdef int L = lmbda
    cdef int M_nrows = len(M)
    cdef tuple R

    # The comments and variables of this code have been written for a different
    # definition of a Difference Matrix, i.e. with cols/rows reversed. This will
    # be corrected soon.
    #
    # We transpose the matrix as a first step.
    if M:
        for i,row in enumerate(M):
            if len(row) != len(M[0]):
                if verbose:
                    print "Rows 0 and {} do not have the same length".format(i)
                return False
    M = zip(*M)

    # Height of the matrix
    if G_card*lmbda != M_nrows:
        if verbose:
            print "The matrix has {} columns and lambda.|G|={}.{}={}".format(M_nrows,lmbda,G_card,lmbda*G_card)
        return False

    # Width of the matrix
    for R in M:
        if len(R)!=K:
            if verbose:
                print "A column of the matrix has length {}!=k(={})".format(len(R),K)
            return False

    # When |G|=0
    if M_nrows == 0:
        return True

    # Map group element with integers
    cdef list int_to_group = list(G)
    cdef dict group_to_int = {v:i for i,v in enumerate(int_to_group)}

    # Allocations
    cdef int ** x_minus_y     = <int **> sage_malloc(G_card*sizeof(int *))
    cdef int * x_minus_y_data = <int *>  sage_malloc(G_card*G_card*sizeof(int))
    cdef int * M_c            = <int *>  sage_malloc(k*M_nrows*sizeof(int))
    cdef int * G_seen         = <int *>  sage_malloc(G_card*sizeof(int))
    if (x_minus_y == NULL or x_minus_y_data == NULL or M_c == NULL or G_seen == NULL):
        sage_free(x_minus_y)
        sage_free(x_minus_y_data)
        sage_free(G_seen)
        sage_free(M_c)
        raise MemoryError

    # The "x-y" table. If g_i, g_j \in G, then x_minus_y[i][j] is equal to
    # group_to_int[g_i-g_j]
    zero, op, inv = group_law(G)
    x_minus_y[0] = x_minus_y_data
    for i in range(1,G_card):
        x_minus_y[i] = x_minus_y[i-1] + G_card

    for j,Gj in enumerate(int_to_group):
        minus_Gj = inv(Gj)
        assert op(Gj, minus_Gj) == zero
        for i,Gi in enumerate(int_to_group):
            x_minus_y[i][j] = group_to_int[op(Gi,minus_Gj)]

    # A copy of the matrix
    for i,R in enumerate(M):
        for j,x in enumerate(R):
            M_c[i*K+j] = group_to_int[x]

    # We are now ready to test every pair of columns
    for i in range(K):
        for j in range(i+1,K):
            memset(G_seen, 0, G_card*sizeof(int))
            for ii in range(M_nrows):
                G_seen[x_minus_y[M_c[ii*K+i]][M_c[ii*K+j]]] += 1

            for ii in range(G_card):
                if G_seen[ii] != L:
                    if verbose:
                        print ("Rows {} and {} do not generate all elements of G "
                         "exactly lambda(={}) times. The element {} appeared {} "
                         "times as a difference.".format(i,j,L,int_to_group[ii],G_seen[ii]))
                    sage_free(x_minus_y_data)
                    sage_free(x_minus_y)
                    sage_free(G_seen)
                    sage_free(M_c)
                    return False

    sage_free(x_minus_y_data)
    sage_free(x_minus_y)
    sage_free(G_seen)
    sage_free(M_c)
    return True
