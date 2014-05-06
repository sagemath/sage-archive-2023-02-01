r"""
Cython functions for combinatorial designs

This module implements the design methods that need to be somewhat efficient.

Functions
---------
"""
include "sage/misc/bitset.pxi"

def is_orthogonal_array(OA, int k, int n, int t=2, verbose=False, terminology="OA"):
    r"""
    Check that the integer matrix `M` is an `OA(k,n,t)`.

    See :func:`~sage.combinat.designs.orthogonal_arrays.orthogonal_array`
    for a definition.

    INPUT:

    - ``OA`` -- the Orthogonal Array to be tested

    - ``k,n,t`` integers.

    - ``verbose`` (boolean) -- whether to display some information when ``OA``
      is not an orthogona array `OA(k,n)`.

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
        Rows 0 and 3 are not orthogonal
        False
    """
    cdef int n2 = n*n
    cdef int x
    if any(len(R) != k for R in OA):
        if verbose:
            print {"OA"   : "Some row does not have length "+str(k),
                   "MOLS" : "The number of matrices is not "+str(k)}[terminology]
        return False

    if len(OA) != n2:
        if verbose:
            print {"OA"   : "The number of rows is {} instead of {}^2={}".format(len(OA),n,n2),
                   "MOLS" : "All matrices do not have dimension n^2={}^2".format(n)}[terminology]
        return False

    cdef int i,j,l

    # A copy of OA
    cdef unsigned short * OAc = <unsigned short *> sage_malloc(k*n2*sizeof(unsigned short))

    # A cache to multiply by n
    cdef unsigned int * times_n = <unsigned int *> sage_malloc((n+1)*sizeof(unsigned int))

    cdef unsigned short * C1
    cdef unsigned short * C2

    # failed malloc ?
    if OAc is NULL or times_n is NULL:
        if OAc is not NULL:
            sage_free(OAc)
        if times_n is not NULL:
            sage_free(times_n)
        raise MemoryError

    # Filling OAc
    for i,R in enumerate(OA):
        for j,x in enumerate(R):
            OAc[j*n2+i] = x
            if x < 0 or x >= n:
                if verbose:
                    print {"OA"   : "{} is not in the interval [0..{}]".format(x,n-1),
                           "MOLS" : "Entry {} was expected to be in the interbal [0..{}]".format(x,n-1)}[terminology]
                sage_free(OAc)
                sage_free(times_n)
                return False

    # Filling times_n
    for i in range(n+1):
        times_n[i] = i*n

    # A bitset to keep trac of pairs of values
    cdef bitset_t seen
    bitset_init(seen, n2)

    for i in range(k): # For any column C1
        C1 = OAc+i*n2
        for j in range(i+1,k): # For any column C2 > C1
            C2 = OAc+j*n2
            bitset_set_first_n(seen, 0) # No pair has been seen yet
            for l in range(n2):
                bitset_add(seen,times_n[C1[l]]+C2[l])

            if bitset_len(seen) != n2: # Have we seen all pairs ?
                sage_free(OAc)
                sage_free(times_n)
                bitset_free(seen)
                if verbose:
                    print {"OA"   : "Rows {} and {} are not orthogonal".format(i,j),
                           "MOLS" : "Matrices {} and {} are not orthogonal".format(i,j)}[terminology]
                return False

    sage_free(OAc)
    sage_free(times_n)
    bitset_free(seen)
    return True
