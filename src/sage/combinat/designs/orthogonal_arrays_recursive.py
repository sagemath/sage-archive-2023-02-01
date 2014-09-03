r"""
Orthogonal arrays (Recursive constructions)

This module implements several functions to find recursive constructions of
:mod:`Orthogonal Arrays <sage.combinat.designs.orthogonal_arrays>`.

The main function of this module, i.e. :func:`find_recursive_construction`,
queries all implemented recursive constructions of designs. It is used by
Sage's function
:func:`~sage.combinat.designs.orthogonal_arrays.orthogonal_array`.

REFERENCES:

.. [AC07] Concerning eight mutually orthogonal latin squares
  Julian R. Abel, Nicholas Cavenagh
  Journal of Combinatorial Designs
  Vol. 15, n.3, pp. 255-261
  2007

Functions
---------
"""
from sage.misc.cachefunc import cached_function
from orthogonal_arrays import orthogonal_array
from designs_pyx import is_orthogonal_array

@cached_function
def find_recursive_construction(k,n):
    r"""
    Find a recursive construction of a `OA(k,n)`

    This determines whether an `OA(k,n)` can be built through the following
    constructions:

    - :func:`simple_wilson_construction` (0, 1 or 2 truncated columns)
    - :func:`construction_3_3`
    - :func:`construction_3_4`
    - :func:`construction_3_5`
    - :func:`construction_3_6`
    - :func:`construction_q_x`
    - :func:`thwart_lemma_3_5`
    - :func:`thwart_lemma_4_1`
    - :func:`three_factor_product`
    - :func:`brouwer_separable_design`

    INPUT:

    - ``k,n`` (integers)

    OUTPUT:

    Return a pair ``f,args`` such that ``f(*args)`` returns the requested `OA`
    if possible, and ``False`` otherwise.

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays_recursive import find_recursive_construction
        sage: from sage.combinat.designs.orthogonal_arrays import is_orthogonal_array
        sage: count = 0
        sage: for n in range(10,150):
        ....:     k = designs.orthogonal_array(None,n,existence=True)
        ....:     if find_recursive_construction(k,n):
        ....:         count = count + 1
        ....:         f,args = find_recursive_construction(k,n)
        ....:         OA = f(*args)
        ....:         assert is_orthogonal_array(OA,k,n,2,verbose=True)
        sage: print count
        56
    """
    assert k > 3

    for find_c in [find_product_decomposition,
                   find_wilson_decomposition_with_one_truncated_group,
                   find_wilson_decomposition_with_two_truncated_groups,
                   find_construction_3_3,
                   find_construction_3_4,
                   find_construction_3_5,
                   find_construction_3_6,
                   find_q_x,
                   find_thwart_lemma_3_5,
                   find_thwart_lemma_4_1,
                   find_three_factor_product,
                   find_brouwer_separable_design]:
        res = find_c(k,n)
        if res:
            return res
    return False

def find_product_decomposition(k,n):
    r"""
    Look for a factorization of `n` in order to build an `OA(k,n)`.

    If Sage can build a `OA(k,n_1)` and a `OA(k,n_2)` such that `n=n_1\times
    n_2` then a `OA(k,n)` can be built by a product construction (which
    correspond to Wilson's construction with no truncated column). This
    function look for a pair of integers `(n_1,n_2)` with `n1 \leq n_2`, `n_1
    \times n_2 = n` and such that both an `OA(k,n_1)` and an `OA(k,n_2)` are
    available.

    INPUT:

    - ``k,n`` (integers) -- see above.

    OUTPUT:

    A pair ``f,args`` such that ``f(*args)`` is an `OA(k,n)` or ``False`` if no
    product decomposition was found.

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays_recursive import find_product_decomposition
        sage: f,args = find_product_decomposition(6, 84)
        sage: args
        (6, 7, 12, ())
        sage: _ = f(*args)
    """
    from sage.rings.arith import divisors
    for n1 in divisors(n)[1:-1]: # we ignore 1 and n
        n2 = n//n1  # n2 is decreasing along the loop
        if n2 < n1:
            break
        if orthogonal_array(k, n1, existence=True) and orthogonal_array(k, n2, existence=True):
            return simple_wilson_construction, (k,n1,n2,())
    return False

def find_wilson_decomposition_with_one_truncated_group(k,n):
    r"""
    Helper function for Wilson's construction with one truncated column.

    This function looks for possible integers `m,t,u` satisfying that `mt+u=n` and
    such that Sage knows how to build a `OA(k,m)`, `OA(k,m+1)`, `OA(k+1,t)` and a
    `OA(k,u)`.

    INPUT:

    - ``k,n`` (integers) -- see above

    OUTPUT:

    A pair ``f,args`` such that ``f(*args)`` is an `OA(k,n)` or ``False`` if no
    decomposition with one truncated block was found.

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays_recursive import find_wilson_decomposition_with_one_truncated_group
        sage: f,args = find_wilson_decomposition_with_one_truncated_group(4,38)
        sage: args
        (4, 5, 7, (3,))
        sage: _ = f(*args)

        sage: find_wilson_decomposition_with_one_truncated_group(4,20)
        False
    """
    # If there exists a TD(k+1,t) then k+1 < t+2, i.e. k <= t
    for r in range(max(1,k),n-1):
        u = n%r
        # We ensure that 1<=u, and that there can exists a TD(k,u), i.e k<u+2
        # (unless u == 1)
        if u == 0 or (u>1 and k >= u+2):
            continue

        m = n//r
        # If there exists a TD(k,m) then k<m+2
        if k >= m+2:
            break

        if (orthogonal_array(k  ,m  , existence=True) and
            orthogonal_array(k  ,m+1, existence=True) and
            orthogonal_array(k+1,r  , existence=True) and
            orthogonal_array(k  ,u  , existence=True)):
            return simple_wilson_construction, (k,r,m,(u,))

    return False

def find_wilson_decomposition_with_two_truncated_groups(k,n):
    r"""
    Helper function for Wilson's construction with two trucated columns.

    Look for integers `r,m,r_1,r_2` satisfying `n=rm+r_1+r_2` and `1\leq r_1,r_2<r`
    and such that the following designs exist : `OA(k+2,r)`, `OA(k,r1)`,
    `OA(k,r2)`, `OA(k,m)`, `OA(k,m+1)`, `OA(k,m+2)`.

    INPUT:

    - ``k,n`` (integers) -- see above

    OUTPUT:

    A pair ``f,args`` such that ``f(*args)`` is an `OA(k,n)` or ``False`` if no
    decomposition with two truncated blocks was found.

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays_recursive import find_wilson_decomposition_with_two_truncated_groups
        sage: f,args = find_wilson_decomposition_with_two_truncated_groups(5,58)
        sage: args
        (5, 7, 7, (4, 5))
        sage: _ = f(*args)
    """
    for r in [1] + range(k+1,n-2): # as r*1+1+1 <= n and because we need
                                   # an OA(k+2,r), necessarily r=1 or r >= k+1
        if not orthogonal_array(k+2,r,existence=True):
            continue
        m_min = (n - (2*r-2))//r
        m_max = (n - 2)//r
        if m_min > 1:
            m_values = range(max(m_min,k-1), m_max+1)
        else:
            m_values = [1] + range(k-1, m_max+1)
        for m in m_values:
            r1_p_r2 = n-r*m # the sum of r1+r2
                            # it is automatically >= 2 since m <= m_max
            if (r1_p_r2 > 2*r-2 or
                not orthogonal_array(k,m  ,existence=True) or
                not orthogonal_array(k,m+1,existence=True) or
                not orthogonal_array(k,m+2,existence=True)):
                continue

            r1_min = r1_p_r2 - (r-1)
            r1_max = min(r-1, r1_p_r2)
            if r1_min > 1:
                r1_values = range(max(k-1,r1_min), r1_max+1)
            else:
                r1_values = [1] + range(k-1, r1_max+1)
            for r1 in r1_values:
                if not orthogonal_array(k,r1,existence=True):
                    continue
                r2 = r1_p_r2-r1
                if orthogonal_array(k,r2,existence=True):
                    assert n == r*m+r1+r2
                    return simple_wilson_construction, (k,r,m,(r1,r2))
    return False

def simple_wilson_construction(k,r,m,u):
    r"""
    Return an `OA(k,rm + \sum u_i)` from Wilson construction.

    INPUT:

    - ``k,r,m`` -- integers

    - ``u`` -- list of positive integers

    .. TODO::

        As soon as wilson construction accepts an empty master design we should
        remove this intermediate functions.

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays_recursive import simple_wilson_construction
        sage: from sage.combinat.designs.designs_pyx import is_orthogonal_array

        sage: OA = simple_wilson_construction(6,7,12,())
        sage: is_orthogonal_array(OA,6,84)
        True

        sage: OA = simple_wilson_construction(4,5,7,(3,))
        sage: is_orthogonal_array(OA,4,38)
        True

        sage: OA = simple_wilson_construction(5,7,7,(4,5))
        sage: is_orthogonal_array(OA,5,58)
        True
    """
    from sage.combinat.designs.orthogonal_arrays import wilson_construction, OA_relabel

    n = r*m + sum(u)
    n_trunc = len(u)
    OA = orthogonal_array(k+n_trunc,r,check=False)
    matrix = [range(r)]*k
    for uu in u:
            matrix.append(range(uu)+[None]*(r-uu))
    OA = OA_relabel(OA,k+n_trunc,r,matrix=matrix)

    return wilson_construction(OA,k,r,m,n_trunc,u,False)

def find_construction_3_3(k,n):
    r"""
    Finds a decomposition for construction 3.3 from [AC07]_

    INPUT:

    - ``k,n`` (integers)

    .. SEEALSO::

        :func:`construction_3_3`

    OUTPUT:

    A pair ``f,args`` such that ``f(*args)`` returns the requested OA.

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays_recursive import find_construction_3_3
        sage: find_construction_3_3(11,177)[1]
        (11, 11, 16, 1)
        sage: find_construction_3_3(12,11)
    """
    for mm in range(k-1,n//2+1):
        if (not orthogonal_array(k ,mm  , existence=True) or
            not orthogonal_array(k ,mm+1, existence=True)):
            continue

        for nn in range(2,n//mm+1):
            i = n-nn*mm
            if i<=0:
                continue

            if (orthogonal_array(k+i, nn  , existence=True) and
                orthogonal_array(k  , mm+i, existence=True)):
                return construction_3_3, (k,nn,mm,i)

def construction_3_3(k,n,m,i):
    r"""
    Return an `OA(k,nm+i)`.

    This is Wilson's construction with `i` truncated columns of size 1 and such
    that a block `B_0` of the incomplete OA intersects all truncated columns. As
    a consequence, all other blocks intersect only `0` or `1` of the last `i`
    columns. This allow to consider the block `B_0` only up to its first `k`
    coordinates and then use a `OA(k,i)` instead of a `OA(k,m+i) - i.OA(k,1)`.

    This is construction 3.3 from [AC07]_.

    INPUT:

    - ``k,n,m,i`` (integers) such that the following designs are available :
      `OA(k,n)`, `OA(k,m)`, `OA(k,m+1)`, `OA(k,r)`.

    .. SEEALSO::

        :func:`find_construction_3_3`

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays_recursive import find_construction_3_3
        sage: from sage.combinat.designs.orthogonal_arrays_recursive import construction_3_3
        sage: from sage.combinat.designs.orthogonal_arrays import is_orthogonal_array
        sage: k=11;n=177
        sage: is_orthogonal_array(construction_3_3(*find_construction_3_3(k,n)[1]),k,n,2)
        True
    """
    from orthogonal_arrays import wilson_construction, OA_relabel, incomplete_orthogonal_array
    # Builds an OA(k+i,n) containing a block [0]*(k+i)
    OA = incomplete_orthogonal_array(k+i,n,(1,))
    OA = [[(x+1)%n for x in B] for B in OA]

    # Truncated version
    OA = [B[:k]+[0 if x == 0 else None for x in B[k:]] for B in OA]

    OA = wilson_construction(OA,k,n,m,i,[1]*i,check=False)[:-i]
    matrix = [range(m)+range(n*m,n*m+i)]*k
    OA.extend(OA_relabel(orthogonal_array(k,m+i),k,m+i,matrix=matrix))
    assert is_orthogonal_array(OA,k,n*m+i)
    return OA

def find_construction_3_4(k,n):
    r"""
    Finds a decomposition for construction 3.4 from [AC07]_

    INPUT:

    - ``k,n`` (integers)

    .. SEEALSO::

        :func:`construction_3_4`

    OUTPUT:

    A pair ``f,args`` such that ``f(*args)`` returns the requested OA.

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays_recursive import find_construction_3_4
        sage: find_construction_3_4(8,196)[1]
        (8, 25, 7, 12, 9)
        sage: find_construction_3_4(9,24)
    """
    for mm in range(k-1,n//2+1):
        if (not orthogonal_array(k,mm+0,existence=True) or
            not orthogonal_array(k,mm+1,existence=True) or
            not orthogonal_array(k,mm+2,existence=True)):
            continue

        for nn in range(2,n//mm+1):
            i = n-nn*mm
            if i<=0:
                continue

            for s in range(1,min(i,nn)):
                r = i-s
                if (orthogonal_array(k+r+1,nn,existence=True) and
                    orthogonal_array(k    , s,existence=True) and
                    (orthogonal_array(k,mm+r,existence=True) or orthogonal_array(k,mm+r+1,existence=True))):
                    return construction_3_4, (k,nn,mm,r,s)

def construction_3_4(k,n,m,r,s):
    r"""
    Return a `OA(k,nm+rs)`.

    This is Wilson's construction applied to a truncated `OA(k+r+1,n)` with `r`
    columns of size `1` and one column of size `s`.

    The unique elements of the `r` truncated columns are picked so that a block
    `B_0` contains them all.

    - If there exists an `OA(k,m+r+1)` the column of size `s` is truncated in
      order to intersect `B_0`.

    - Otherwise, if there exists an `OA(k,m+r)`, the last column must not
      intersect `B_0`

    This is construction 3.4 from [AC07]_.

    INPUT:

    - ``k,n,m,r,s`` (integers) -- we assume that `s<n` and `1\leq r,s`

      The following designs must be available: `OA(k,n)`, `OA(k,m)`,
      `OA(k,m+1)`, `OA(k,m+2)`, `OA(k,s)`. Additionnally, it requires either a
      `OA(k,m+r)` or a `OA(k,m+r+1)`.

    .. SEEALSO::

        :func:`find_construction_3_4`

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays_recursive import find_construction_3_4
        sage: from sage.combinat.designs.orthogonal_arrays_recursive import construction_3_4
        sage: from sage.combinat.designs.orthogonal_arrays import is_orthogonal_array
        sage: k=8;n=196
        sage: is_orthogonal_array(construction_3_4(*find_construction_3_4(k,n)[1]),k,n,2)
        True
    """
    from orthogonal_arrays import wilson_construction, OA_relabel
    assert s<n
    master_design = orthogonal_array(k+r+1,n)

    # Defines the first k+r columns of the matrix of labels
    matrix = [range(n)]*k + [[None]*n]*(r) + [[None]*n]
    B0 = master_design[0]
    for i in range(k,k+r):
        matrix[i][B0[i]] = 0

    # Last column
    if orthogonal_array(k,m+r,existence=True):
        last_group = [x for x in range(s+1) if x != B0[-1]][:s]
    elif orthogonal_array(k,m+r+1,existence=True):
        last_group = [x for x in range(s+1) if x != B0[-1]][:s-1] + [B0[-1]]
    else:
        raise RuntimeError

    for i,x in enumerate(last_group):
        matrix[-1][x] = i

    OA = OA_relabel(master_design,k+r+1,n, matrix=matrix)
    OA = wilson_construction(OA,k,n,m,r+1,[1]*r+[s],check=False)
    return OA

def find_construction_3_5(k,n):
    r"""
    Finds a decomposition for construction 3.5 from [AC07]_

    INPUT:

    - ``k,n`` (integers)

    .. SEEALSO::

        :func:`construction_3_5`

    OUTPUT:

    A pair ``f,args`` such that ``f(*args)`` returns the requested OA.

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays_recursive import find_construction_3_5
        sage: find_construction_3_5(8,111)[1]
        (8, 13, 6, 11, 11, 11)
        sage: find_construction_3_5(9,24)
    """
    from sage.combinat.integer_list import IntegerListsLex

    for mm in range(2,n//2+1):
        if (mm+3 >= n or
            not orthogonal_array(k,mm+1,existence=True) or
            not orthogonal_array(k,mm+2,existence=True) or
            not orthogonal_array(k,mm+3,existence=True)):
            continue

        for nn in range(2,n//mm+1):
            i = n-nn*mm
            if i<=0:
                continue

            if not orthogonal_array(k+3,nn,existence=True):
                continue

            for r,s,t in IntegerListsLex(i,length=3,ceiling=[nn-1,nn-1,nn-1]):
                if (r <= s and
                    (nn-r-1)*(nn-s) < t and
                    (r==0 or orthogonal_array(k,r,existence=True)) and
                    (s==0 or orthogonal_array(k,s,existence=True)) and
                    (t==0 or orthogonal_array(k,t,existence=True))):
                    return construction_3_5, (k,nn,mm,r,s,t)

def construction_3_5(k,n,m,r,s,t):
    r"""
    Return an `OA(k,nm+r+s+t)`.

    This is exactly Wilson's construction with three truncated groups
    except we make sure that all blocks have size `>k`, so we don't
    need a `OA(k,m+0)` but only `OA(k,m+1)`, `OA(k,m+2)` ,`OA(k,m+3)`.

    This is construction 3.5 from [AC07]_.

    INPUT:

    - ``k,n,m`` (integers)

    - ``r,s,t`` (integers) -- sizes of the three truncated groups,
      such that `r\leq s` and `(q-r-1)(q-s) \geq (q-s-1)*(q-r)`.

    The following designs must be available : `OA(k,n)`, `OA(k,r)`, `OA(k,s)`,
    `OA(k,t)`, `OA(k,m+1)`, `OA(k,m+2)`, `OA(k,m+3)`.

    .. SEEALSO::

        :func:`find_construction_3_5`

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays_recursive import find_construction_3_5
        sage: from sage.combinat.designs.orthogonal_arrays_recursive import construction_3_5
        sage: from sage.combinat.designs.orthogonal_arrays import is_orthogonal_array
        sage: k=8;n=111
        sage: is_orthogonal_array(construction_3_5(*find_construction_3_5(k,n)[1]),k,n,2)
        True
    """
    from orthogonal_arrays import wilson_construction, OA_relabel
    assert r <= s
    q = n
    assert (q-r-1)*(q-s) >= (q-s-1)*(q-r)
    master_design = orthogonal_array(k+3,q)

    # group k+1 has cardinality r
    # group k+2 has cardinality s
    # group k+3 has cardinality t

    # Taking q-s blocks going through 0 in the last block
    blocks_crossing_0 = [B[-3:] for B in master_design if B[-1] == 0][:q-s]

    # defining the undeleted points of the groups k+1,k+2
    group_k_1 = [x[0] for x in blocks_crossing_0]
    group_k_1 = [x for x in range(q) if x not in group_k_1][:r]

    group_k_2 = [x[1] for x in blocks_crossing_0]
    group_k_2 = [x for x in range(q) if x not in group_k_2][:s]

    # All blocks that have a deleted point in groups k+1 and k+2 MUST contain a
    # point in group k+3
    group_k_3 = [B[-1] for B in master_design if B[-3] not in group_k_1 and B[-2] not in group_k_2]
    group_k_3 = list(set(group_k_3))
    assert len(group_k_3) <= t
    group_k_3.extend([x for x in range(q) if x not in group_k_3])
    group_k_3 = group_k_3[:t]

    # Relabelling the OA
    r1 = [None]*q
    r2 = [None]*q
    r3 = [None]*q
    for i,x in enumerate(group_k_1):
        r1[x] = i
    for i,x in enumerate(group_k_2):
        r2[x] = i
    for i,x in enumerate(group_k_3):
        r3[x] = i

    OA = OA_relabel(master_design, k+3,q, matrix=[range(q)]*k+[r1,r2,r3])
    OA = wilson_construction(OA,k,q,m,3,[r,s,t], check=False)
    return OA

def find_construction_3_6(k,n):
    r"""
    Finds a decomposition for construction 3.6 from [AC07]_

    INPUT:

    - ``k,n`` (integers)

    .. SEEALSO::

        :func:`construction_3_6`

    OUTPUT:

    A pair ``f,args`` such that ``f(*args)`` returns the requested OA.

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays_recursive import find_construction_3_6
        sage: find_construction_3_6(8,95)[1]
        (8, 13, 7, 4)
        sage: find_construction_3_6(8,98)
    """
    from sage.rings.arith import is_prime_power

    for mm in range(k-1,n//2+1):
        if (not orthogonal_array(k,mm+0,existence=True) or
            not orthogonal_array(k,mm+1,existence=True) or
            not orthogonal_array(k,mm+2,existence=True)):
            continue

        for nn in range(2,n//mm+1):
            i = n-nn*mm
            if i<=0:
                continue

            if (is_prime_power(nn) and
                orthogonal_array(k+i,nn,existence=True)):
                return construction_3_6, (k,nn,mm,i)

def construction_3_6(k,n,m,i):
    r"""
    Return a `OA(k,nm+i)`

    This is Wilson's construction with `r` columns of order `1`, in which each
    block intersects at most two truncated columns. Such a design exists when
    `n` is a prime power and is returned by :func:`OA_and_oval`.

    INPUT:

    - ``k,n,m,i`` (integers) -- `n` must be a prime power. The following designs
      must be available: `OA(k+r,q)`, `OA(k,m)`, `OA(k,m+1)`, `OA(k,m+2)`.

    This is construction 3.6 from [AC07]_.

    .. SEEALSO::

        - :func:`construction_3_6`
        - :func:`OA_and_oval`

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays_recursive import find_construction_3_6
        sage: from sage.combinat.designs.orthogonal_arrays_recursive import construction_3_6
        sage: from sage.combinat.designs.orthogonal_arrays import is_orthogonal_array
        sage: k=8;n=95
        sage: is_orthogonal_array(construction_3_6(*find_construction_3_6(k,n)[1]),k,n,2)
        True
    """
    from orthogonal_arrays import wilson_construction
    OA = OA_and_oval(n)
    OA = [B[:k+i] for B in OA]
    OA = [B[:k] + [x if x==0 else None for x in B[k:]] for B in OA]
    OA = wilson_construction(OA,k,n,m,i,[1]*i)
    assert is_orthogonal_array(OA,k,n*m+i)
    return OA

def OA_and_oval(q):
    r"""
    Return a `OA(q+1,q)` whose blocks contains `\leq 2` zeroes in the last `q`
    columns.

    This `OA` is build from a projective plane of order `q`, in which there
    exists an oval `O` of size `q+1` (i.e. a set of `q+1` points no three of which
    are [colinear/contained in a common set of the projective plane]).

    Removing an element `x\in O` and all sets that contain it, we obtain a
    `TD(q+1,q)` in which `O` intersects all columns except one. As `O` is an
    oval, no block of the `TD` intersects it more than twice.

    INPUT:

    - ``q`` -- a prime power

    .. NOTE::

            This function is called by :func:`construction_3_6`, an
            implementation of Construction 3.6 from [AC07]_.

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays_recursive import OA_and_oval
        sage: _ = OA_and_oval

    """
    from sage.rings.arith import is_prime_power
    from sage.combinat.designs.block_design import projective_plane
    from orthogonal_arrays import OA_relabel

    assert is_prime_power(q)
    B = projective_plane(q, check=False)

    # We compute the oval with a linear program
    from sage.numerical.mip import MixedIntegerLinearProgram
    p = MixedIntegerLinearProgram()
    b = p.new_variable(binary=True)
    V = B.ground_set()
    p.add_constraint(p.sum([b[i] for i in V]) == q+1)
    for bl in B:
        p.add_constraint(p.sum([b[i] for i in bl]) <= 2)
    p.solve()
    b = p.get_values(b)
    oval = [x for x,i in b.items() if i]
    assert len(oval) == q+1

    # We remove one element from the oval
    x = oval.pop()
    oval.sort()

    # We build the TD by relabelling the point set, and removing those which
    # contain x.
    r = {}
    B = list(B)
    # (this is to make sure that the first set containing x in B is the one
    # which contains no other oval point)

    B.sort(key=lambda b:int(any([xx in oval for xx in b])))
    BB = []
    for b in B:
        if x in b:
            for xx in b:
                if xx == x:
                    continue
                r[xx] = len(r)
        else:
            BB.append(b)

    assert len(r) == (q+1)*q # all points except x have an image
    assert len(set(r.values())) == len(r) # the images are different

    # Relabelling/sorting the blocks and the oval
    BB = [[r[xx] for xx in b] for b in BB]
    oval = [r[xx] for xx in oval]

    for b in BB:
        b.sort()
    oval.sort()

    # Turning the TD into an OA
    BB = [[xx%q for xx in b] for b in BB]
    oval = [xx%q for xx in oval]
    assert len(oval) == q

    # We relabel the "oval" as relabelled as [0,...,0]
    OA = OA_relabel(BB+([[0]+oval]),q+1,q,blocks=[[0]+oval])
    OA = [[(x+1)%q for x in B] for B in OA]
    OA.remove([0]*(q+1))

    assert all(sum([xx == 0 for xx in b[1:]]) <= 2 for b in OA)
    return OA

def construction_q_x(k,q,x,check=True):
    r"""
    Return an `OA(k,(q-1)*(q-x)+x+2)` using the `q-x` construction.

    Let `v=(q-1)*(q-x)+x+2`. If there exists a projective plane of order `q`
    (e.g. when `q` is a prime power) and `0<x<q` then there exists a
    `(v-1,\{q-x-1,q-x+1\})`-GDD of type `(q-1)^{q-x}(x+1)^1` (see [Greig99]_ or
    Theorem 2.50, section IV.2.3 of [DesignHandbook]_). By adding to the ground
    set one point contained in all groups of the GDD, one obtains a
    `(v,\{q-x-1,q-x+1,q,x+2\})`-PBD with exactly one set of size `x+2`.

    Thus, assuming that we have the following:

    - `OA(k,q-x-1)-(q-x-1).OA(k,1)`
    - `OA(k,q-x+1)-(q-x+1).OA(k,1)`
    - `OA(k,q)-q.OA(k,1)`
    - `OA(k,x+2)`

    Then we can build from the PBD an `OA(k,v)`.

    Construction of the PBD (shared by Julian R. Abel):

        Start with a resolvable `(q^2,q,1)`-BIBD and put the points into a `q\times q`
        array so that rows form a parallel class and columns form another.

        Now delete:

        - All `x(q-1)` points from the first `x` columns and not in the first
          row

        - All `q-x` points in the last `q-x` columns AND the first row.

        Then add a point `p_1` to the blocks that are rows. Add a second point
        `p_2` to the `q-x` blocks that are columns of size `q-1`, plus the first
        row of size `x+1`.

    INPUT:

    - ``k,q,x`` -- integers such that `0<x<q` and such that Sage can build:

        - A projective plane of order `q`
        - `OA(k,q-x-1)-(q-x-1).OA(k,1)`
        - `OA(k,q-x+1)-(q-x+1).OA(k,1)`
        - `OA(k,q)-q.OA(k,1)`
        - `OA(k,x+2)`

    - ``check`` -- (boolean) Whether to check that output is correct before
      returning it. As this is expected to be useless (but we are cautious
      guys), you may want to disable it whenever you want speed. Set to
      ``True`` by default.

    .. SEEALSO::

        - :func:`find_q_x`
        - :func:`~sage.combinat.designs.block_design.projective_plane`
        - :func:`~sage.combinat.designs.orthogonal_arrays.orthogonal_array`
        - :func:`~sage.combinat.designs.orthogonal_arrays.OA_from_PBD`

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays_recursive import construction_q_x
        sage: _ = construction_q_x(9,16,6)

    REFERENCES:

    .. [Greig99] Designs from projective planes and PBD bases
      Malcolm Greig
      Journal of Combinatorial Designs
      vol. 7, num. 5, pp. 341--374
      1999
    """
    from sage.combinat.designs.orthogonal_arrays import OA_from_PBD
    from sage.combinat.designs.orthogonal_arrays import incomplete_orthogonal_array

    n = (q-1)*(q-x)+x+2

    # We obtain the qxq matrix from a OA(q,q)-q.OA(1,q). We will need to add
    # blocks corresponding to the rows/columns
    OA = incomplete_orthogonal_array(q,q,(1,)*q)
    TD = [[i*q+xx for i,xx in enumerate(B)] for B in OA]

    # Add rows, extended with p1 and p2
    p1 = q**2
    p2 = p1+1
    TD.extend([[ii*q+i for ii in range(q)]+[p1] for i in range(1,q)])
    TD.append( [ii*q   for ii in range(q)]+[p1,p2])

    # Add Columns. We do not add some columns which would have size 1 after we
    # delete points.
    #
    # TD.extend([range(i*q,(i+1)*q) for i in range(x)])
    TD.extend([range(i*q,(i+1)*q)+[p2] for i in range(x,q)])

    points_to_delete = set([i*q+j for i in range(x) for j in range(1,q)]+[i*q for i in range(x,q)])
    points_to_keep = set(range(q**2+2))-points_to_delete
    relabel = {i:j for j,i in enumerate(points_to_keep)}

    # PBD is a (n,[q,q-x-1,q-x+1,x+2])-PBD
    PBD = [[relabel[xx] for xx in B if not xx in points_to_delete] for B in TD]

    # Taking the unique block of size x+2
    assert map(len,PBD).count(x+2)==1
    for B in PBD:
        if len(B) == x+2:
            break

    # We call OA_from_PBD without the block of size x+2 as there may not exist a
    # OA(k,x+2)-(x+2).OA(k,1)
    PBD.remove(B)
    OA = OA_from_PBD(k,(q-1)*(q-x)+x+2,PBD,check=False)

    # Filling the hole
    for xx in B:
        OA.remove([xx]*k)

    for BB in orthogonal_array(k,x+2):
        OA.append([B[x] for x in BB])

    if check:
        assert is_orthogonal_array(OA,k,n,2)

    return OA

def find_q_x(k,n):
    r"""
    Find integers `q,x` such that the `q-x` construction yields an `OA(k,n)`.

    See the documentation of :func:`construction_q_x` to find out what
    hypotheses the integers `q,x` must satisfy.

    .. WARNING::

        For efficiency reasons, this function checks that Sage can build an
        `OA(k+1,q-x-1)` and an `OA(k+1,q-x+1)`, which is stronger than checking
        that Sage can build a `OA(k,q-x-1)-(q-x-1).OA(k,1)` and a
        `OA(k,q-x+1)-(q-x+1).OA(k,1)`. The latter would trigger a lot of
        independent set computations in
        :func:`sage.combinat.designs.orthogonal_arrays.incomplete_orthogonal_array`.

    INPUT:

    - ``k,n`` (integers)

    .. SEEALSO::

        :func:`construction_q_x`

    EXAMPLE::

        sage: from sage.combinat.designs.orthogonal_arrays_recursive import find_q_x
        sage: find_q_x(10,9)
        False
        sage: find_q_x(9,158)[1]
        (9, 16, 6)
    """
    from sage.rings.arith import is_prime_power
    # n = (q-1)*(q-x) + x + 2
    #   = q^2 - q*x - q + 2*x + 2

    for q in range(max(3,k+2),n):
        # n-q**2+q-2 = 2x-qx
        #            = x(2-q)
        x = (n-q**2+q-2)//(2-q)
        if (x < q and
            0 < x and
            n == (q-1)*(q-x)+x+2 and
            is_prime_power(q) and
            orthogonal_array(k+1,q-x-1,existence=True) and
            orthogonal_array(k+1,q-x+1,existence=True) and
            # The next is always True, because q is a prime power
            # orthogonal_array(k+1,q,existence=True) and
            orthogonal_array(k, x+2 ,existence=True)):
            return construction_q_x, (k,q,x)
    return False

def find_thwart_lemma_3_5(k,N):
    r"""
    A function to find the values for which one can apply the
    Lemma 3.5 from [Thwarts]_.

    OUTPUT:

    A pair ``(f,args)`` such that ``f(*args)`` returns an `OA(k,n)` or ``False``
    if the construction is not available.

    .. SEEALSO::

        :func:`thwart_lemma_3_5`

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays_recursive import find_thwart_lemma_3_5
        sage: from sage.combinat.designs.designs_pyx import is_orthogonal_array

        sage: f,args = find_thwart_lemma_3_5(7,66)
        sage: args
        (7, 9, 7, 1, 1, 1, 0, False)
        sage: OA = f(*args)
        sage: is_orthogonal_array(OA,7,66,2)
        True

        sage: f,args = find_thwart_lemma_3_5(6,100)
        sage: args
        (6, 8, 10, 8, 7, 5, 0, True)
        sage: OA = f(*args)
        sage: is_orthogonal_array(OA,6,100,2)
        True

    Some values from [Thwarts]_::

        sage: kn = ((10,1046), (10,1048), (10,1059), (11,1524),
        ....:       (11,2164), (12,3362), (12,3992),  (12,3994))
        sage: for k,n in kn:
        ....:     print k,n,find_thwart_lemma_3_5(k,n)[1]
        10 1046 (10, 13, 79, 9, 1, 0, 9, False)
        10 1048 (10, 13, 79, 9, 1, 0, 11, False)
        10 1059 (10, 13, 80, 9, 1, 0, 9, False)
        11 1524 (11, 19, 78, 16, 13, 13, 0, True)
        11 2164 (11, 27, 78, 23, 19, 16, 0, True)
        12 3362 (12, 16, 207, 13, 13, 11, 13, True)
        12 3992 (12, 19, 207, 16, 13, 11, 19, True)
        12 3994 (12, 19, 207, 16, 13, 13, 19, True)

        sage: for k,n in kn:                                                     # not tested -- too long
        ....:     assert designs.orthogonal_array(k,n,existence=True) is True    # not tested -- too long
    """
    from sage.rings.arith import prime_powers

    k = int(k)
    N = int(N)

    for n in prime_powers(k+2,N-2): # There must exist a OA(k+3,n) thus n>=k+2
                                    # At least 3 columns are nonempty thus n<N-2

        # we look for (m,n,a,b,c,d) with N = mn + a + b + c (+d) and
        # 0 <= a,b,c,d <= n
        # hence we have N/n-4 <= m <= N/n

        # 1. look for m,a,b,c,d with complement=False
        # (we restrict to a >= b >= c)
        for m in xrange(max(k-1,(N+n-1)//n-4), N//n+1):
            if not (orthogonal_array(k,m+0,existence=True) and
                    orthogonal_array(k,m+1,existence=True) and
                    orthogonal_array(k,m+2,existence=True)):
                continue

            NN = N - n*m
            # as a >= b >= c and d <= n we can restrict the start of the loops
            for a in range(max(0, (NN-n+2)//3), min(n, NN)+1): # (NN-n+2)//3 <==> ceil((NN-n)/3)x
                if not orthogonal_array(k,a,existence=True):
                    continue
                for b in range(max(0, (NN-n-a+1)//2), min(a, n+1-a, NN-a)+1):
                    if not orthogonal_array(k,b,existence=True):
                        continue
                    for c in range(max(0, NN-n-a-b), min(b, n+1-a-b, NN-a-b)+1):
                        if not orthogonal_array(k,c,existence=True):
                            continue

                        d = NN - (a + b + c)  # necessarily 0 <= d <= n
                        if d == 0:
                            return thwart_lemma_3_5, (k,n,m,a,b,c,0,False)
                        elif (k+4 <= n+1 and
                            orthogonal_array(k,d,existence=True) and
                            orthogonal_array(k,m+3,existence=True)):
                            return thwart_lemma_3_5, (k,n,m,a,b,c,d,False)

        # 2. look for m,a,b,c,d with complement=True
        # (we restrict to a >= b >= c)
        for m in xrange(max(k-2,N//n-4), (N+n-1)//n):
            if not (orthogonal_array(k,m+1,existence=True) and
                    orthogonal_array(k,m+2,existence=True) and
                    orthogonal_array(k,m+3,existence=True)):
                continue

            NN = N - n*m
            for a in range(max(0, (NN-n+2)//3), min(n, NN)+1): # (NN-n+2)//3 <==> ceil((NN-n)/3)
                if not orthogonal_array(k,a,existence=True):
                    continue
                na = n-a
                for b in range(max(0, (NN-n-a+1)//2), min(a, NN-a)+1):
                    nb = n-b
                    if na+nb > n+1 or not orthogonal_array(k,b,existence=True):
                        continue
                    for c in range(max(0, NN-n-a-b), min(b, NN-a-b)+1):
                        nc = n-c
                        if na+nb+nc > n+1 or not orthogonal_array(k,c,existence=True):
                            continue

                        d = NN - (a + b + c)  # necessarily d <= n
                        if d == 0:
                            return thwart_lemma_3_5, (k,n,m,a,b,c,0,True)
                        elif (k+4 <= n+1 and
                            orthogonal_array(k,d,existence=True) and
                            orthogonal_array(k,m+4,existence=True)):
                            return thwart_lemma_3_5, (k,n,m,a,b,c,d,True)

    return False

def thwart_lemma_3_5(k,n,m,a,b,c,d=0,complement=False):
    r"""
    Returns an `OA(k,nm+a+b+c+d)`

    *(When `d=0`)*

    According to [Thwarts]_ when `n` is a prime power and `a+b+c\leq n+1`, one
    can build an `OA(k+3,n)` with three truncated columns of sizes `a,b,c` in
    such a way that all blocks have size `\leq k+2`.

    (in order to build a `OA(k,nm+a+b+c)` the following designs must also exist:
    `OA(k,a)`, `OA(k,b)`, `OA(k,c)`, `OA(k,m+0)`, `OA(k,m+1)`, `OA(k,m+2)`)

    Considering the complement of each truncated column, it is also possible to
    build an `OA(k+3,n)` with three truncated columns of sizes `a,b,c` in such a
    way that all blocks have size `>k` whenever `(n-a)+(n-b)+(n-c)\leq n+1`.

    (in order to build a `OA(k,nm+a+b+c)` the following designs must also exist:
    `OA(k,a)`, `OA(k,b)`, `OA(k,c)`, `OA(k,m+1)`, `OA(k,m+2)`, `OA(k,m+3)`)

    Here is the proof of Lemma 3.5 from [Thwarts]_ enriched with explanations
    from Julian R. Abel:

        For any prime power `n` one can build `k-1` MOLS by associating to every
        nonzero `x\in \mathbb F_n` the latin square:

        .. MATH::

            M_x(i,j) = i+x*j \text{ where }i,j\in \mathbb F_n

        In particular `M_1(i,j)=i+j`, whose `n` columns and lines are indexed by
        the elements of `\mathbb F_n`. If we order the elements of `\mathbb F_n`
        as `0,1,...,n-1,x+0,...,x+n-1,x^2+0,...` and reorder the columns
        and lines of `M_1` accordingly, the top-left `a\times b` squares
        contains at most `a+b-1` distinct symbols.

    *(When* `d\neq 0` *)*

    If there exists an `OA(k+3,n)` with three truncated columns of sizes `a,b,c`
    in such a way that all blocks have size `\leq k+2`, by truncating
    arbitrarily another column to size `d` one obtains an `OA` with 4 truncated
    columns whose blocks miss at least one value. Thus, following the proof
    again one can build an `OA(k+4)` with four truncated columns of sizes
    `a,b,c,d` with blocks of size `\leq k+3`.

    (in order to build a `OA(k,nm+a+b+c+d)` the following designs must also
    exist: `OA(k,a)`, `OA(k,b)`, `OA(k,c)`, `OA(k,d)`, `OA(k,m+0)`, `OA(k,m+1)`,
    `OA(k,m+2)`, `OA(k,m+3)`)

    As before, this also shows that one can build an `OA(k+4,n)` with four
    truncated columns of sizes `a,b,c,d` in such a way that all blocks have size
    `>k` whenever `(n-a)+(n-b)+(n-c)\leq n+1`

    (in order to build a `OA(k,nm+a+b+c+d)` the following designs must also
    exist: `OA(k,n-a)`, `OA(k,n-b)`, `OA(k,n-c)`, `OA(k,d)`, `OA(k,m+1)`,
    `OA(k,m+2)`, `OA(k,m+3)`, `OA(k,m+4)`)

    INPUT:

    - ``k,n,m,a,b,c,d`` -- integers which must satisfy the constraints above. In
      particular, `a+b+c\leq n+1` must hold. By default, `d=0`.

    - ``complement`` (boolean) -- whether to complement the sets, i.e. follow
      the `n-a,n-b,n-c` variant described above.

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays_recursive import thwart_lemma_3_5
        sage: from sage.combinat.designs.designs_pyx import is_orthogonal_array
        sage: OA = thwart_lemma_3_5(6,23,7,5,7,8)
        sage: is_orthogonal_array(OA,6,23*7+5+7+8,2)
        True

    With sets of parameters from [Thwarts]_::

        sage: l = [
        ....:    [11, 27, 78, 16, 17, 25, 0],
        ....:    [12, 19, 208, 11, 13, 16, 0],
        ....:    [12, 19, 208, 13, 13, 16, 0],
        ....:    [10, 13, 78, 9, 9, 13, 1],
        ....:    [10, 13, 79, 9, 9, 13, 1]]
        sage: for k,n,m,a,b,c,d in l:                                       # not tested -- too long
        ....:     OA = thwart_lemma_3_5(k,n,m,a,b,c,d,complement=True)      # not tested -- too long
        ....:     assert is_orthogonal_array(OA,k,n*m+a+b+c+d,verbose=True) # not tested -- too long

    REFERENCE:

    .. [Thwarts] Thwarts in transversal designs
      Charles J.Colbourn, Jeffrey H. Dinitz, Mieczyslaw Wojtas.
      Designs, Codes and Cryptography 5, no. 3 (1995): 189-197.
    """
    from sage.rings.arith import is_prime_power
    from sage.rings.finite_rings.constructor import FiniteField as GF
    from sage.combinat.designs.orthogonal_arrays import wilson_construction

    if complement:
        a,b,c = n-a,n-b,n-c

    assert is_prime_power(n), "n(={}) must be a prime power".format(n)
    assert a<=n and b<=n and c<=n and d<=n, "a,b,c,d (={},{},{},{}) must be <=n(={})".format(a,b,c,d,n)
    assert a+b+c<=n+1, "{}={}+{}+{}=a+b+c>n+1={}+1 violates the assumptions".format(a+b+c,a,b,c,n)
    assert k+3+bool(d) <= n+1, "There exists no OA({},{}).".format(k+3+bool(d),n)
    G = GF(n,prefix='x',conway=True)
    G_set = sorted(G) # sorted by lexicographic order, G[1] = 1
    assert G_set[0] == G.zero() and G_set[1] == G.one(), "problem with the ordering of {}".format(G)
    G_to_int = {v:i for i,v in enumerate(G_set)}

    # Builds an OA(n+1,n) whose last n-1 colums are
    #
    # \forall x \in G and x!=0, C_x(i,j) = i+x*j
    #
    # (only the necessary columns are built)
    OA = [[G_to_int[i+x*j] for i in G_set for j in G_set] for x in G_set[1:k+2+bool(d)]]
    # Adding the first two trivial columns
    OA.insert(0,[j for i in range(n) for j in range(n)])
    OA.insert(0,[i for i in range(n) for j in range(n)])
    OA=zip(*OA)
    OA.sort()

    # Moves the first three columns to the end
    OA = [list(B[3:]+B[:3]) for B in OA]

    # Set of values in the axb square
    third_complement= set([B[-1] for B in OA if B[-3] < a and B[-2] < b])

    assert n-len(third_complement) >= c

    # The keepers
    first_set  = range(a)
    second_set = range(b)
    third_set  = [x for x in range(n) if x not in third_complement][:c]

    last_sets  = [first_set,second_set,third_set]

    if complement:
        last_sets = [set(range(n)).difference(s) for s in last_sets]

    sizes = map(len,last_sets)
    last_sets_dict = [{v:i for i,v in enumerate(s)} for s in last_sets]

    # Truncating the OA
    for i,D in enumerate(last_sets_dict):
        kk = len(OA[0])-3+i
        for R in OA:
            R[kk] = D[R[kk]] if R[kk] in D else None

    if d:
        for R in OA:
            if R[-4] >= d:
                R[-4] = None
        sizes.insert(0,d)

    return wilson_construction(OA,k,n,m,len(sizes),sizes, check=False)

def find_thwart_lemma_4_1(k,n):
    r"""
    Finds a decomposition for Lemma 4.1 from [Thwarts]_.

    INPUT:

    - ``k,n`` (integers)

    .. SEEALSO::

        :func:`thwart_lemma_4_1`

    OUTPUT:

    A pair ``f,args`` such that ``f(*args)`` returns the requested OA.

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays_recursive import find_thwart_lemma_4_1
        sage: find_thwart_lemma_4_1(10,408)[1]
        (10, 13, 28)
        sage: find_thwart_lemma_4_1(10,50)
        False
    """
    from sage.rings.arith import factor
    #      n  = nn*mm+4(nn-2)
    # <=> n+8 = nn(mm+4)
    #
    # nn is a prime power dividing n+8
    for nn in (p**i for p,imax in factor(n+8) for i in range(1,imax+1)):
        mm = (n+8)//nn-4
        if (k+4 > nn+1 or
            mm <= 1 or
            nn % 3 == 2 or
            not orthogonal_array(k,nn-2,existence=True) or
            not orthogonal_array(k,mm+1,existence=True) or
            not orthogonal_array(k,mm+3,existence=True) or
            not orthogonal_array(k,mm+4,existence=True)):
            continue

        return thwart_lemma_4_1,(k,nn,mm)

    return False

def thwart_lemma_4_1(k,n,m):
    r"""
    Returns an `OA(k,nm+4(n-2))`.

    Implements Lemma 4.1 from [Thwarts]_.

        If `n\equiv 0,1\pmod{3}` is a prime power, then there exists a truncated
        `OA(n+1,n)` whose last four columns have size `n-2` and intersect every
        block on `1,3` or `4` values. Consequently, if there exists an
        `OA(k,m+1)`, `OA(k,m+3)`, `OA(k,m+4)` and a `OA(k,n-2)` then there
        exists an `OA(k,nm+4(n-2)`

        Proof: form the transversal design by removing one point of the
        `AG(2,3)` (Affine Geometry) contained in the Desarguesian Projective
        Plane `PG(2,n)`.

    The affine geometry on 9 points contained in the projective geometry
    `PG(2,n)` is given explicitly in [OS64]_ (Thanks to Julian R. Abel for
    finding the reference!).

    REFERENCES:

    .. [OS64] Finite projective planes with affine subplanes,
      T. G. Ostrom and F. A. Sherk.
      Canad. Math. Bull vol7 num.4 (1964)
    """
    from sage.combinat.designs.designs_pyx import is_orthogonal_array
    from sage.rings.finite_rings.constructor import FiniteField
    from sage.rings.arith import is_prime_power
    from block_design import DesarguesianProjectivePlaneDesign
    from itertools import chain

    assert is_prime_power(n), "n(={}) must be a prime power"
    assert k+4 <= n+1

    q = n
    K = FiniteField(q, 'x')
    relabel = {x:i for i,x in enumerate(K)}
    PG = DesarguesianProjectivePlaneDesign(q,check=False).blocks(copy=False)

    if q % 3 == 0:
        t = K.one()
    elif q%3 == 1:
        t = K.multiplicative_generator()**((q-1)//3)
    else:
        raise ValueError("q(={}) must be congruent to 0 or 1 mod 3".format(q))

    # The projective plane is labelled with integer coordinates. This code
    # relabels to integers the following points (given by homogeneous
    # coordinates in the projective space):
    #
    # - (1+t,t,1+t), (1,1,1), (1+t,t,t), (1,1,2), (0,0,1), (1,0,1), (0,1,1+t),
    #   (0,1,1), (1,0,-t)
    points = [(1+t,t,1+t), (1,1,1), (1+t,t,t), (1,1,2), (0,0,1), (1,0,1), (0,1,1+t), (0,1,1), (1,0,-t)]
    points = [map(K,t) for t in points] # triples of K^3
    AG_2_3 = []
    for x,y,z in points:
        if z!=0:
            x,y,z = x/z,y/z,z/z
            AG_2_3.append(relabel[x]+n*relabel[y])
        elif y!=0:
            x,y,z=x/y,y/y,z
            AG_2_3.append(q**2+relabel[x])
        else:
            AG_2_3.append(q**2+q)

    AG_2_3 = set(AG_2_3)

    # All blocks of PG should intersect 'AG_2_3' on !=2 AG_2_3.
    assert all(len(AG_2_3.intersection(B)) != 2 for B in PG)

    p = list(AG_2_3)[0]
    # We now build a TD from the PG by removing p, in such a way that the last
    # two elements of the last 4 columns are elements of AG_2_3
    blocks = []
    columns = []
    for B in PG:
        if p not in B:
            blocks.append(B)
        else:
            B.remove(p)
            columns.append(B)

    # The columns containing elements from the AG are the last ones, and those
    # elements should be the last two
    columns.sort(key=lambda x:len(AG_2_3.intersection(x)))
    for i in range(4):
        columns[-i-1].sort(key=lambda x: int(x in AG_2_3))

    relabel = {v:i for i,v in enumerate(chain(columns))}

    TD = [sorted(relabel[x] for x in B) for B in blocks]

    # We build the OA, removing unnecessary columns
    OA = [[x%q for x in B[-k-4:]] for B in TD]
    for B in OA:
        for i in range(4):
            if B[k+i] >= n-2:
                B[k+i] = None

    return wilson_construction(OA,k,n,m,4,[n-2,]*4,check=False)

def find_three_factor_product(k,n):
    r"""
    Finds a decomposition for a three-factor product from [DukesLing14]_

    INPUT:

    - ``k,n`` (integers)

    .. SEEALSO::

        :func:`three_factor_product`

    OUTPUT:

    A pair ``f,args`` such that ``f(*args)`` returns the requested OA.

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays_recursive import find_three_factor_product
        sage: find_three_factor_product(10,648)[1]
        (9, 8, 9, 9)
        sage: find_three_factor_product(10,50)
        False
    """
    # we want to write n=n1*n2*n3 where n1<=n2<=n3 and we can build:
    # - a OA(k-1,n1)
    # - a OA( k ,n2)
    # - a OA( k ,n3)
    from sage.rings.arith import divisors
    for n1 in divisors(n)[1:-1]:
        if not orthogonal_array(k-1,n1,existence=True):
            continue
        for n2 in divisors(n//n1):
            n3 = n//n1//n2
            if (n2<n1 or
                n3<n2 or
                not orthogonal_array(k,n2,existence=True) or
                not orthogonal_array(k,n3,existence=True)):
                continue
            return three_factor_product,(k-1,n1,n2,n3)

    return False

def three_factor_product(k,n1,n2,n3,check=False):
    r"""
    Returns an `OA(k+1,n_1n_2n_3)`

    The three factor product construction from [DukesLing14]_ does the following:

        If `n_1\leq n_2\leq n_3` are such that there exists an
        `OA(k,n_1)`, `OA(k+1,n_2)` and `OA(k+1,n_3)`, then there exists a
        `OA(k+1,n_1n_2n_3)`.

    It works with a modified product of orthogonal arrays ([Rees93]_, [Rees00]_)
    which keeps track of parallel classes in the `OA` (the definition is given
    for transversal designs).

        A subset of blocks in an `TD(k,n)` is called a `c`-parallel class if
        every point is covered exactly `c` times. A 1-parallel class is a
        parallel class.

    The modified product:

        If there exists an `OA(k,n_1)`, and if there exists an `OA(k,n_2)` whose
        blocks are partitionned into `s` `n_1`-parallel classes and `n_2-sn_1`
        parallel classes, then there exists an `OA(k,n_1n_2)` whose blocks can
        be partitionned into `sn_1^2` parallel classes and
        `(n_1n_2-sn_1^2)/n_1=n_2-sn_1` `n_1`-parallel classes.

        Proof:

        - The product of the blocks of a parallel class with an `OA(k,n_1)`
          yields an `n_1`-parallel class of an `OA(k,n_1n_2)`.

        - The product of the blocks of a `n_1`-parallel class of `OA(k,n_2)`
          with an `OA(k,n_1)` can be done in such a way that it yields `n_1n_2`
          parallel classes of `OA(k,n_1n_2)`. Those classes cover exactly the
          pairs that woud have been covered with the usual product.

          This can be achieved by simple cyclic permutations. Let us build the
          product of the `n_1`-parallel class `\mathcal P\subseteq OA(k,n_2)`
          with `OA(k,n_1)`: when computing the product of `P\in\mathcal P` with
          `B^1\in OA(k,n_1)` the `i`-th coordinate should not be `(B^1_i,P_i)`
          but `(B^1_i+r,P_i)` (the sum is mod `n_1`) where `r` is the number of
          blocks of `\mathcal P` we have already processed whose `i`-th
          coordinate is equal to `P_i` (note that `r< n_1` as `\mathcal P` is
          `n_1`-parallel).

    With these tools, one can obtain the designs promised by the three factors
    construction applied to `k,n_1,n_2,n_3` (thanks to Julian R. Abel's help):

        1) Let `s` be the largest integer `\leq n_3/n_1`. Apply the product
           construction to `OA(k,n_1)` and a resolvable `OA(k,n_3)` whose blocks
           are partitionned into `s` `n_1`-parallel classes and `n_3-sn_1`
           parallel classes. It results in a `OA(k,n_1n_3)` partitionned into
           `sn_1^2` parallel classes plus `(n_1n_3-sn_1^2)/n_1=n_3-sn_1`
           `n_1`-parallel classes.

        2) Add `n_3-n_1` parallel classes to every `n_1`-parallel class to turn
           them into `n_3`-parallel classes. Apply the product construction to
           this partitionned `OA(k,n_1n_3)` with a resolvable `OA(k,n_2)`.

        3) As `OA(k,n_2)` is resolvable, the `n_2`-parallel classes of
           `OA(k,n_1n_2n_3)` are actually the union of `n_2` parallel classes,
           thus the `OA(k,n_1n_2n_3)` is resolvable and can be turned into an
           `OA(k+1,n_1n_2n_3)`

    INPUT:

    - ``k,n1,n2,n3`` (integers)

    - ``check`` -- (boolean) Whether to check that everything is going smoothly
      while the design is being built. It is disabled by default, as the
      constructor of orthogonal arrays checks the final design anyway.

    EXAMPLE::

        sage: from sage.combinat.designs.designs_pyx import is_orthogonal_array
        sage: from sage.combinat.designs.orthogonal_arrays_recursive import three_factor_product

        sage: OA = three_factor_product(4,4,4,4)
        sage: is_orthogonal_array(OA,5,64)
        True

        sage: OA = three_factor_product(4,3,4,5)
        sage: is_orthogonal_array(OA,5,60)
        True

        sage: OA = three_factor_product(5,4,5,7)
        sage: is_orthogonal_array(OA,6,140)
        True

        sage: OA = three_factor_product(9,8,9,9) # long time
        sage: is_orthogonal_array(OA,10,8*9*9)   # long time
        True

    REFERENCE:

    .. [DukesLing14] A three-factor product construction for mutually orthogonal latin squares,
      Peter J. Dukes, Alan C.H. Ling,
      http://arxiv.org/abs/1401.1466

    .. [Rees00] Truncated Transversal Designs: A New Lower Bound on the Number of Idempotent MOLS of Side,
      Rolf S. Rees,
      Journal of Combinatorial Theory, Series A 90.2 (2000): 257-266.

    .. [Rees93] Two new direct product-type constructions for resolvable group-divisible designs,
      Rolf S. Rees,
      Journal of Combinatorial Designs 1.1 (1993): 15-26.
    """
    from itertools import izip
    assert n1<=n2 and n2<=n3

    def assert_c_partition(classs,k,n,c):
        r"""
        Makes sure that ``classs`` contains blocks `B` of size `k` such that the list of
        ``B[i]`` covers `[n]` exactly `c` times for every index `i`.
        """
        c = int(c)
        assert all(len(B)==k for B in classs), "A block has length {}!=k(={})".format(len(B),k)
        assert len(classs) == n*c, "not the right number of blocks"
        for p in zip(*classs):
            assert all(x==i//c for i,x in enumerate(sorted(p))), "A class is not c(={})-parallel".format(c)

    def product_with_parallel_classes(OA1,k,g1,g2,g1_parall,parall,check=True):
        r"""
        Returns the product of two OA while keeping track of parallel classes

        INPUT:

        - ``OA1`` (an `OA(k,g_1)`

        - ``k,g1,g2`` integers

        - ``g1_parall`` -- list of `g_1`-parallel classes

        - ``parall`` -- list of parallel classes

        .. NOTE::

            The list ``g1_parall+parall`` should be an `OA(k,g_2)`

        OUTPUT:

        Two lists of classes ``g1_parall`` and ``parallel`` which are respectively
        `g_1`-parallel and parallel classes such that ``g1_parall+parallel`` is an
        `OA(k,g1*g2)``.
        """
        if check:
            for classs in g1_parall:
                assert_c_partition(classs,k,g2,g1)
            for classs in parall:
                assert_c_partition(classs,k,g2,1)

        # New parallel classes, built from a g1-parallel class with shifted copies
        # of OA1

        new_parallel_classes = []
        for classs2 in g1_parall:

            # Keep track of how many times we saw each point of [k]x[g2]
            count = [[0]*g2 for _ in range(k)]

            copies_of_OA1 = []
            for B2 in classs2:
                copy_of_OA1 = []

                shift = [count[i][x2] for i,x2 in enumerate(B2)]
                assert max(shift) < g1

                for B1 in OA1:
                    copy_of_OA1.append([x2*g1+(x1+sh)%g1 for sh,x1,x2 in izip(shift,B1,B2)])

                copies_of_OA1.append(copy_of_OA1)

                # Update the counts
                for i,x2 in enumerate(B2):
                    count[i][x2] += 1

            new_parallel_classes.extend(map(list,izip(*copies_of_OA1)))

        # New g1-parallel classes, each one built from the product of a parallel
        # class with a OA1

        new_g1_parallel_classes = []
        for classs2 in parall:
            disjoint_copies_of_OA1 = []
            for B2 in classs2:
                for B1 in OA1:
                    disjoint_copies_of_OA1.append([x2*g1+x1 for x1,x2 in izip(B1,B2)])
            new_g1_parallel_classes.append(disjoint_copies_of_OA1)

        # Check our stuff before we return it
        if check:
            profile = [i for i in range(g2*g1) for _ in range(g1)]
            for classs in new_g1_parallel_classes:
                assert_c_partition(classs,k,g2*g1,g1)
            profile = range(g2*g1)
            for classs in new_parallel_classes:
                assert_c_partition(classs,k,g2*g1,1)

        return new_g1_parallel_classes, new_parallel_classes

    # The three factors product construction begins !
    #
    # OA1 and resolvable OA2 and OA3
    OA1 = orthogonal_array(k,n1)
    OA3 = orthogonal_array(k+1,n3)
    OA3.sort()
    OA3 = [B[1:] for B in OA3]
    OA2 = orthogonal_array(k+1,n2)
    OA2.sort()
    OA2 = [B[1:] for B in OA2]

    # We split OA3 into as many n1-parallel classes as possible, i.e. n3//n1 classes of size n1*n3
    OA3_n1_parall = [OA3[i:i+n1*n3] for i in range(0,(n3-n1)*n3,n1*n3)]

    # Leftover blocks become parallel classes. We must split them into slices of
    # length n3
    OA3_parall    = [OA3[i:i+n3] for i in range(len(OA3_n1_parall)*n1*n3, len(OA3), n3)]

    # First product: OA1 and OA3
    n1_parall, parall = product_with_parallel_classes(OA1,k,n1,n3,OA3_n1_parall,OA3_parall,check=check)

    if check:
        OA_13 = [block for classs in parall+n1_parall for block in classs]
        assert is_orthogonal_array(OA_13,k,n1*n3,2,1)

    # Add parallel classes to turn the n1-parall classes into n2-parallel classes
    for classs in n1_parall:
        for i in range(n2-n1):
            classs.extend(parall.pop())

    n2_parall = n1_parall
    del n1_parall

    # We compute the product of OA2 with our decomposition of OA1xOA2 into
    # n2-parallel classes and parallel classes
    n2_parall, parall = product_with_parallel_classes(OA2,k,n2,n1*n3,n2_parall,parall,check=check)
    for n2_classs in n2_parall:
        for i in range(n2):
            partition = [B for j in range(n1*n3) for B in n2_classs[j*n2**2+i*n2:j*n2**2+(i+1)*n2]]
            parall.append(partition)

    # That's what we fought for: this design is resolvable, so let's add a last
    # column to them
    for i,classs in enumerate(parall):
        for B in classs:
            B.append(i)

    OA = [block for classs in parall for block in classs]

    if check:
        assert is_orthogonal_array(OA,k+1,n1*n2*n3,2,1)

    return OA

def find_brouwer_separable_design(k,n):
    r"""
    Find integers `t,q,x` such that :func:`brouwer_separable_design` gives a `OA(k,t(q^2+q+1)+x)`.

    INPUT:

    - ``k,n`` (integers)

    The assumptions made on the parameters `t,q,x` are explained in the
    documentation of :func:`brouwer_separable_design`.

    EXAMPLE::

        sage: from sage.combinat.designs.orthogonal_arrays_recursive import find_brouwer_separable_design
        sage: find_brouwer_separable_design(5,13)[1]
        (5, 1, 3, 0)
        sage: find_brouwer_separable_design(5,14)
        False
    """
    from sage.rings.arith import prime_powers
    for q in prime_powers(2,n):
        baer_subplane_size = q**2+q+1
        if baer_subplane_size > n:
            break
        #                       x <= q^2+1
        # <=>        n-t(q^2+q+1) <= q^2+1
        # <=>             n-q^2-1 <= t(q^2+q+1)
        # <=> (n-q^2-1)/(q^2+q+1) <= t

        min_t = (n-q**2-1)//baer_subplane_size
        max_t = min(n//baer_subplane_size,q**2-q+1)

        for t in range(min_t,max_t+1):
            x = n - t*baer_subplane_size
            e1 = int(x != q**2-q-t)
            e2 = int(x != 1)
            e3 = int(x != q**2)
            e4 = int(x != t+q+1)

            # i)
            if (x == 0 and
                orthogonal_array(k, t,existence=True)  and
                orthogonal_array(k,t+q,existence=True)):
                return brouwer_separable_design, (k,t,q,x)

            # ii)
            elif (x == t+q and
                  orthogonal_array(k+e3,  t  ,existence=True) and
                  orthogonal_array(  k , t+q ,existence=True) and
                  orthogonal_array(k+1 ,t+q+1,existence=True)):
                return brouwer_separable_design, (k,t,q,x)

            # iii)
            elif (x == q**2-q+1-t and
                  orthogonal_array(  k  ,  x  ,existence=True) and
                  orthogonal_array( k+e2, t+1 ,existence=True) and
                  orthogonal_array( k+1 , t+q ,existence=True)):
                return brouwer_separable_design, (k,t,q,x)

            # iv)
            elif (x == q**2+1 and
                  orthogonal_array(  k  ,  x  ,existence=True) and
                  orthogonal_array( k+e4, t+1 ,existence=True) and
                  orthogonal_array( k+1 ,t+q+1,existence=True)):
                return brouwer_separable_design, (k,t,q,x)

            # v)
            elif (0<x and x<q**2-q+1-t and (e1 or e2) and
                  orthogonal_array(  k  ,  x  ,existence=True) and
                  orthogonal_array( k+e1,  t  ,existence=True) and
                  orthogonal_array( k+e2, t+1 ,existence=True) and
                  orthogonal_array( k+1 , t+q ,existence=True)):
                return brouwer_separable_design, (k,t,q,x)

            # vi)
            elif (t+q<x and x<q**2+1 and (e3 or e4) and
                  orthogonal_array(  k  ,  x  ,existence=True) and
                  orthogonal_array( k+e3,  t  ,existence=True) and
                  orthogonal_array( k+e4, t+1 ,existence=True) and
                  orthogonal_array( k+1 ,t+q+1,existence=True)):
                return brouwer_separable_design, (k,t,q,x)

    return False

def _reorder_matrix(matrix):
    r"""
    Return a matrix which is obtained from ``matrix`` by permutation of each row
    in which each column contain every symbol exactly once.

    The input must be a `N \times k` matrix with entries in `\{0,\ldots,N-1\}`
    such that:
    - the symbols on each row are distinct (and hence can be identified with
      subsets of `\{0,\ldots,N-1\}`),
    - each symbol appear exactly `k` times.

    The problem is equivalent to an edge coloring of a bipartite graph. This
    function is used by :func:`brouwer_separable_design`.

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays_recursive import _reorder_matrix
        sage: N = 4; k = 3
        sage: M = [[0,1,2],[0,1,3],[0,2,3],[1,2,3]]
        sage: M2 = _reorder_matrix(M)
        sage: all(set(M2[i][0] for i in range(N)) == set(range(N)) for i in range(k))
        True

        sage: M =[range(10)]*10
        sage: N = k = 10
        sage: M2 = _reorder_matrix(M)
        sage: all(set(M2[i][0] for i in range(N)) == set(range(N)) for i in range(k))
        True
    """
    from sage.graphs.graph import Graph

    N = len(matrix)
    k = len(matrix[0])

    g = Graph()
    g.add_edges((x,N+i) for i,S in enumerate(matrix) for x in S)
    matrix = []
    for _ in range(k):
        matching = g.matching(algorithm="LP")
        col = [0]*N
        for x,i,_ in matching:
            if i<N:
                x,i=i,x
            col[i-N] = x
        matrix.append(col)
        g.delete_edges(matching)

    return zip(*matrix)

def brouwer_separable_design(k,t,q,x,check=False,verbose=False):
    r"""
    Returns a `OA(k,t(q^2+q+1)+x)` using Brouwer's result on separable designs.

    This method is an implementation of Brouwer's construction presented in
    [Brouwer80]_. It consists in a systematic application of the usual
    transformation from PBD to OA, applied to a specific PBD.

    **Baer subplanes**

    When `q` is a prime power, the projective plane `PG(2,q^2)` can be
    partitionned into subplanes `PG(2,q)` (called Baer subplanes), giving
    `PG(2,q^2)=B_1\cup \dots\cup B_{q^2-q+1}`. As a result, every line of the
    `PG(2,q^2)` intersects one of the subplane on `q+1` points and all others on
    `1` point.

    The `OA` are built by considering `B_1\cup\dots\cup B_t`, for a total of
    `t(q^2+q+1)` points (to which `x` new points are then added). The blocks of
    this subdesign belong to two categories:

    * The blocks of size `t`: they come from the lines which intersect a
      `B_i` on `q+1` points for some `i>t`. The blocks of size `t` can be partitionned
      into `q^2-q+t-1` parallel classes according to their associated subplane `B_i`
      with `i>t`.

    * The blocks of size `q+t`: those blocks form a symmetric design, as every
      point is incident with `q+t` of them.

    **Constructions**

    In the following, we write `N=t(q^2+q+1)+x`. The code is also heavily
    commented, and will clear any doubt.

    * i) `x=0`: in that case we build a resolvable `OA(k-1,N)` that will then be
      completed into an `OA(k,N)`.

        * *Sets of size* `t`)

          We take the product of each parallel class with the parallel classes
          of a resolvable `OA(k-1,t)-t.OA(k-1,t)`, yielding new parallel
          classes.

        * *Sets of size* `q+t`)

          A `N \times (q+t)` array is built whose rows are the sets of size
          `q+t` such that every value appears once per column. For each block of
          a `OA(k-1,q+t)-(q+t).OA(k-1,t)`, the product with the rows of the
          matrix yields a parallel class.

    * ii) `x=q+t`

        * *Sets of size* `t`)

          Each set of size `t` gives a `OA(k,t)-t.OA(k,1)`, except if there is
          only one parallel class in which case a `OA(k,t)` is sufficient.

        * *Sets of size* `q+t`)

          A `(N-x) \times (q+t)` array `M` is built whose `N-x` rows are the
          sets of size `q+t` such that every value appears once per column. For
          each of the new `x=q+t` points `p_1,\dots,p_{q+t}` we build a matrix
          `M_i` obtained from `M` by adding a column equal to `(p_i,p_i,p_i\dots
          )`. We add to the OA the product of all rows of the `M_i` with the
          block of the `x=q+t` parallel classes of a resolvable
          `OA(k,t+q+1)-(t+q+1).OA(k,1)`.

        * *Set of size* `x`) An `OA(k,x)`

    * iii) `x = q^2-q+1-t`

        * *Sets of size* `t`)

          All blocks of the `i`-th parallel class are extended with the `i`-th
          new point. The blocks are then replaced by a `OA(k,t+1)-(t+1).OA(k,1)`
          or, if there is only one parallel class (i.e. `x=1`) by a
          `OA(k,t+1)-OA(k,1)`.

        * *Set of size* `q+t`)

          They are replaced by `OA(k,q+t)-(q+t).OA(k,1)`.

        * *Set of size* `x`) An `OA(k,x)`

    * iv) `x = q^2+1`

        * *Sets of size* `t`)

          All blocks of the `i`-th parallel class are extended with the `i`-th
          new point (the other `x-q-t` new points are not touched at this
          step). The blocks are then replaced by a `OA(k,t+1)-(t+1).OA(k,1)` or,
          if there is only one parallel class (i.e. `x=1`) by a
          `OA(k,t+1)-OA(k,1)`.

        * *Sets of size* `q+t`) Same as for ii)

        * *Set of size* `x`) An `OA(k,x)`

    * v) `0<x<q^2-q+1-t`

        * *Sets of size* `t`)

          The blocks of the first `x` parallel class are extended with the `x`
          new points, and replaced with `OA(k.t+1)-(t+1).OA(k,1)` or, if `x=1`,
          by `OA(k.t+1)-.OA(k,1)`

          The blocks of the other parallel classes are replaced by
          `OA(k,t)-t.OA(k,t)` or, if there is only one class left, by
          `OA(k,t)-OA(k,t)`

        * *Sets of size* `q+t`)

          They are replaced with `OA(k,q+t)-(q+t).OA(k,1)`.

        * *Set of size* `x`) An `OA(k,x)`

    * vi) `t+q<x<q^2+1`

        * *Sets of size* `t`) Same as in v) with an `x` equal to `x-q+t`.

        * *Sets of size* `t`) Same as in vii)

        * *Set of size* `x`) An `OA(k,x)`

    INPUT:

    - ``k,t,q,x`` (integers)

    - ``check`` -- (boolean) Whether to check that output is correct before
      returning it. Set to ``False`` by default.

    - ``verbose`` (boolean) -- whether to print some information on the
      construction and parameters being used.

    REFERENCES:

    .. [Brouwer80] A Series of Separable Designs with Application to Pairwise Orthogonal Latin Squares,
      A.E. Brouwer,
      http://www.sciencedirect.com/science/article/pii/S0195669880800199

    EXAMPLES:

    Test all possible cases::

        sage: from sage.combinat.designs.orthogonal_arrays_recursive import brouwer_separable_design
        sage: k,q,t=4,4,3; _=brouwer_separable_design(k,q,t,0,verbose=True)
        Case i) with k=4,q=3,t=4,x=0
        sage: k,q,t=3,3,3; _=brouwer_separable_design(k,t,q,t+q,verbose=True,check=True)
        Case ii) with k=3,q=3,t=3,x=6,e3=1
        sage: k,q,t=3,3,6; _=brouwer_separable_design(k,t,q,t+q,verbose=True,check=True)
        Case ii) with k=3,q=3,t=6,x=9,e3=0
        sage: k,q,t=3,3,6; _=brouwer_separable_design(k,t,q,q**2-q+1-t,verbose=True,check=True)
        Case iii) with k=3,q=3,t=6,x=1,e2=0
        sage: k,q,t=3,4,6; _=brouwer_separable_design(k,t,q,q**2-q+1-t,verbose=True,check=True)
        Case iii) with k=3,q=4,t=6,x=7,e2=1
        sage: k,q,t=3,4,6; _=brouwer_separable_design(k,t,q,q**2+1,verbose=True,check=True)
        Case iv) with k=3,q=4,t=6,x=17,e4=1
        sage: k,q,t=3,2,2; _=brouwer_separable_design(k,t,q,q**2+1,verbose=True,check=True)
        Case iv) with k=3,q=2,t=2,x=5,e4=0
        sage: k,q,t=3,4,7; _=brouwer_separable_design(k,t,q,3,verbose=True,check=True)
        Case v) with k=3,q=4,t=7,x=3,e1=1,e2=1
        sage: k,q,t=3,4,7; _=brouwer_separable_design(k,t,q,1,verbose=True,check=True)
        Case v) with k=3,q=4,t=7,x=1,e1=1,e2=0
        sage: k,q,t=3,4,7; _=brouwer_separable_design(k,t,q,q**2-q-t,verbose=True,check=True)
        Case v) with k=3,q=4,t=7,x=5,e1=0,e2=1
        sage: k,q,t=5,4,7; _=brouwer_separable_design(k,t,q,t+q+3,verbose=True,check=True)
        Case vi) with k=5,q=4,t=7,x=14,e3=1,e4=1
        sage: k,q,t=5,4,8; _=brouwer_separable_design(k,t,q,t+q+1,verbose=True,check=True)
        Case vi) with k=5,q=4,t=8,x=13,e3=1,e4=0
        sage: k,q,t=5,4,8; _=brouwer_separable_design(k,t,q,q**2,verbose=True,check=True)
        Case vi) with k=5,q=4,t=8,x=16,e3=0,e4=1
    """
    from sage.combinat.designs.orthogonal_arrays import OA_from_PBD
    from difference_family import difference_family
    from orthogonal_arrays import incomplete_orthogonal_array
    from sage.rings.arith import is_prime_power

    ###########################################################
    # Part 1: compute the separable PBD on t(q^2+q+1) points. #
    ###########################################################

    assert t<q**2-q+1
    assert x>=0
    assert is_prime_power(q)
    N2 = q**4+q**2+1
    N1 = q**2+  q +1

    # A projective plane on (q^2-q+1)*(q^2+q+1)=q^4+q^2+1 points
    B = difference_family(N2,q**2+1,1)[1][0]
    BIBD = [[(xx+i)%N2 for xx in B] for i in range(N2)]

    # Each congruence class mod q^2-q+1 yields a Baer subplane. Let's check that:
    m = q**2-q+1
    for i in range(m):
        for B in BIBD:
            assert sum((xx%m)==i for xx in B) in [1,q+1], sum((xx%m)==i for xx in B)

    # We are only interested by the points of the first t Baer subplanes (each
    # has size q**2+q+1). Note that each block of the projective plane:
    #
    # - Intersects one Baer plane on q+1 points.
    # - Intersects all other Baer planes on 1 point.
    #
    # When the design its truncated to its first t Baer subplanes, all blocks
    # now have size t or t+q, and cover t(q^2+q+1) points.
    #
    # 1) The blocks of size t can be partitionned into q**2-q+1-t parallel
    #    classes, according to the Baer plane in which they contain q+1
    #    elements.
    #
    # 2) The blocks of size q+t are a symmetric design

    blocks_of_size_q_plus_t = []
    partition_of_blocks_of_size_t = [[] for i in range(m-t)]

    relabel = {i+j*m: N1*i+j for i in range(t) for j in range(N1)}

    for B in BIBD:
        # Find the Baer subplane which B intersects on more than 1 point
        B_mod = sorted(xx%m for xx in B)
        while B_mod.pop(0) != B_mod[0]:
            pass
        plane = B_mod[0]
        if plane < t:
            blocks_of_size_q_plus_t.append([relabel[xx] for xx in B if xx%m<t])
        else:
            partition_of_blocks_of_size_t[plane-t].append([relabel[xx] for xx in B if xx%m<t])

    ###############################################################################
    # Separable design built !
    #-------------------------
    #
    # At this point we have a PBD on t*(q**2+q+1) points. Its blocks are
    # split into:
    #
    # - partition_of_blocks_of_size_t : contains all blocks of size t split into
    #                                   q^2-q+t-1 parallel classes.
    #
    # - blocks_of_size_q_plus_t : contains all t*(q**2+q+1)blocks of size q+t,
    #                             covering the same number of points: it is a
    #                             symmetric design.
    ###############################################################################

    ##############################################
    # Part 2: Build an OA on t(q^2+q+1)+x points #
    ##############################################

    e1 = int(x != q**2-q-t)
    e2 = int(x != 1)
    e3 = int(x != q**2)
    e4 = int(x != t+q+1)
    N  = t*N1+x

    # i)
    if x == 0:

        if verbose:
            print "Case i) with k={},q={},t={},x={}".format(k,q,t,x)

        # 1) We build a resolvable OA(k-1,t)-t.OA(k-1,1).
        #    With it, from every parallel class with blocks of size t we build a
        #    parallel class of a resolvable OA(k-1,N)

        rOA_N_classes = []

        # A resolvable OA(k-1,t)-t.OA(k-1,1)
        OA_t = incomplete_orthogonal_array(k-1,t,[1]*t,resolvable=True)
        OA_t_classes = [OA_t[i*t:(i+1)*t] for i in range(t-1)]

        # We can now build (t-1)(q^2-q+1-t) parallel classes of the resolvable
        # OA(k-1,N)
        for PBD_parallel_class in partition_of_blocks_of_size_t:
            for OA_class in OA_t_classes:
                 rOA_N_classes.append([[B[x] for x in BB]
                                            for BB in OA_class
                                            for B in PBD_parallel_class])

        # 2) We build a Nx(q+t) matrix such that:
        #
        #    a) Each row is a set of size q+t of the PBD
        #    b) an element appears exactly once per column.
        #
        # (This is equivalent to an edge coloring of the (bipartite) incidence
        # graph of points and sets)

        block_of_size_q_plus_t = _reorder_matrix(blocks_of_size_q_plus_t)

        # 3) We now create blocks of an OA(k-1,N) as the product of
        #    a) A set of size q+t (i.e. a row of the matrix)
        #    b) An OA(k-1,q+t)-(q+t).OA(k-1,1)
        #
        # Thanks to the ordering of the points in each set of size q+t, the
        # product of a block B of the incomplete OA with all blocks of size q+t
        # yields a parallel class of an OA(k-1,N)
        OA = incomplete_orthogonal_array(k-1,q+t,[1]*(q+t))
        for B in OA:
            rOA_N_classes.append([[R[x] for x in B] for R in block_of_size_q_plus_t])

        # 4) A last parallel class with blocks [0,0,...], [1,1,...],...
        rOA_N_classes.append([[i]*(k-1) for i in range(N)])

        # 5) We now build the OA(k,N) from the N parallel classes of our resolvable OA(k-1,N)
        OA = [B for classs in rOA_N_classes for B in classs]
        for i,B in enumerate(OA):
            B.append(i//N)

    # ii)
    elif (x == t+q and
          orthogonal_array(k+e3,  t  ,existence=True) and
          orthogonal_array( k  , t+q ,existence=True) and
          orthogonal_array( k+1,t+q+1,existence=True)):

        if verbose:
            print "Case ii) with k={},q={},t={},x={},e3={}".format(k,q,t,x,e3)

        # The sets of size t:
        #
        # This is the usual OA_from_PBD replacement. If there is only one class
        # an OA(k,t) can be used instead of an OA(k+1,t)

        if x == q**2:
            assert e3==0, "equivalent to x==q^2"
            assert len(partition_of_blocks_of_size_t)==1, "also equivalent to exactly one partition into sets of size t"
            OA = [[B[xx] for xx in R] for R in orthogonal_array(k,t) for B in partition_of_blocks_of_size_t[0]]
        else:
            OA = OA_from_PBD(k,N,sum(partition_of_blocks_of_size_t,[]),check=False)[:-N]
            OA.extend([i]*k for i in range(N-x))

        # The sets of size q+t:
        #
        # We build an OA(k,t+q+1)-(t+q+1).OA(k,t+q+1) and the reordered
        # matrix. We then compute the product of every parallel class of the OA
        # (x classes in total) with the rows of the ordered matrix (extended
        # with one of the new x points).

        # Resolvable OA(k,t+q+1)-(t+q+1).OA(k,t+q+1)
        OA_tq1 = incomplete_orthogonal_array(k,t+q+1,[1]*(t+q+1),resolvable=True)
        OA_tq1_classes = [OA_tq1[i*(t+q+1):(i+1)*(t+q+1)] for i in range(t+q)]

        blocks_of_size_q_plus_t = _reorder_matrix(blocks_of_size_q_plus_t)

        for i,classs in enumerate(OA_tq1_classes):
            OA.extend([R[xx] if xx<t+q else N-i-1 for xx in B] for R in blocks_of_size_q_plus_t for B in classs)

        # The set of size x
        OA.extend([N-1-xx for xx in R] for R in orthogonal_array(k,x))

    # iii)
    elif (x == q**2-q+1-t and
          orthogonal_array( k  ,  x  ,existence=True) and # d0
          orthogonal_array(k+e2, t+1 ,existence=True) and # d2-e2
          orthogonal_array(k+1 , t+q ,existence=True)):   # d3-e1
        if verbose:
            print "Case iii) with k={},q={},t={},x={},e2={}".format(k,q,t,x,e2)

        OA = []

        # Each of the x partition into blocks of size t is extended with one of
        # the new x points.

        if x == 1:
            assert e2 == 0, "equivalent to x=1"
            # There is one partition into blocks of size t, which we extend with
            # the new vertex. The OA on t+1 points does not have to be resolvable.

            OA.extend([B[xx] if xx<t else N-1 for xx in R]
                       for R in incomplete_orthogonal_array(k,t+1,[1])
                       for B in partition_of_blocks_of_size_t[0])

        else:
            assert e2 == 1, "equivalent to x!=1"
            # Extending the x partitions into bocks of size t with each of the
            # new x points.

            for i,partition in enumerate(partition_of_blocks_of_size_t):
                for B in partition:
                    B.append(N-i-1)
            OA = OA_from_PBD(k,N,sum(partition_of_blocks_of_size_t,[]),check=False)[:-x]

        # The blocks of size q+t are covered with a resolvable OA(k,q+t)
        OA.extend(OA_from_PBD(k,N,blocks_of_size_q_plus_t,check=False)[:-N])

        # The set of size x
        OA.extend([N-xx-1 for xx in B] for B in orthogonal_array(k,x))


    # iv)
    elif (x == q**2+1 and
          orthogonal_array( k  ,  x  ,existence=True) and # d0
          orthogonal_array(k+e4, t+1 ,existence=True) and # d2-e4
          orthogonal_array(k+ 1,t+q+1,existence=True)):   # d4-1

        if verbose:
            print "Case iv) with k={},q={},t={},x={},e4={}".format(k,q,t,x,e4)

        # Sets of size t:
        #
        # All partitions of t-sets are extended with as many new points

        if e4 == 0:
            # Only one partition into t-sets. The OA(k,t+1) needs not be resolvable
            OA = [[B[xx] if xx<t else N-x for xx in R]
                  for R in incomplete_orthogonal_array(k,t+1,[1])
                  for B in partition_of_blocks_of_size_t[0]]
        else:
            for i,classs in enumerate(partition_of_blocks_of_size_t):
                for B in classs:
                    B.append(N-x+i)
            OA = OA_from_PBD(k,N,sum(partition_of_blocks_of_size_t,[]),check=False)[:-x]

        # The sets of size q+t:
        #
        # We build an OA(k,t+q+1)-(t+q+1).OA(k,t+q+1) and the reordered
        # matrix. We then compute the product of every parallel class of the OA
        # (q+t classes in total) with the rows of the ordered matrix (extended
        # with the last q+t new points).

        # Resolvable OA(k,t+q+1)-(t+q+1).OA(k,t+q+1)
        OA_tq1 = incomplete_orthogonal_array(k,t+q+1,[1]*(t+q+1),resolvable=True)
        OA_tq1_classes = [OA_tq1[i*(t+q+1):(i+1)*(t+q+1)] for i in range(t+q)]

        blocks_of_size_q_plus_t = _reorder_matrix(blocks_of_size_q_plus_t)

        for i,classs in enumerate(OA_tq1_classes):
            OA.extend([R[xx] if xx<t+q else N-i-1 for xx in B] for R in blocks_of_size_q_plus_t for B in classs)

        # Set of size x
        OA_k_x = orthogonal_array(k,x)
        OA.extend([N-i-1 for i in R] for R in OA_k_x)

    # v)
    elif (0<x and x<q**2-q+1-t and (e1 or e2) and # The result is wrong when e1=e2=0
          orthogonal_array(k   ,x  ,existence=True) and # d0
          orthogonal_array(k+e1,t  ,existence=True) and # d1-e1
          orthogonal_array(k+e2,t+1,existence=True) and # d2-e2
          orthogonal_array(k+1,t+q,existence=True)):   # d3-1
        if verbose:
            print "Case v) with k={},q={},t={},x={},e1={},e2={}".format(k,q,t,x,e1,e2)

        OA = []

        # Sets of size t+1
        #
        # We extend x partitions into blocks of size t with the new x elements
        if e2:
            assert x!=1, "equivalent to e2==1"
            for i,classs in enumerate(partition_of_blocks_of_size_t[:x]):
                for B in classs:
                    B.append(N-1-i)
            OA.extend(OA_from_PBD(k,N,sum(partition_of_blocks_of_size_t[:x],[]),check=False)[:-N])

        else:
            assert x==1, "equivalent to e2==0"
            # Only one class, the OA(k,t+1) need not be resolvable.
            OA.extend([B[xx] if xx < t else N-1 for xx in R]
                       for R in incomplete_orthogonal_array(k,t+1,[1])
                       for B in partition_of_blocks_of_size_t[0])

        # Sets of size t
        if e1:
            assert x!=q**2-q-t, "equivalent to e1=1"
            OA.extend(OA_from_PBD(k,N,sum(partition_of_blocks_of_size_t[x:],[]),check=False)[:-N])
        else:
            assert x==q**2-q-t, "equivalent to e1=0"
            # Only one class. The OA(k,t) needs not be resolvable
            OA.extend([B[xx] for xx in R] for R in orthogonal_array(k,t) for B in partition_of_blocks_of_size_t[-1])

        if e1 and e2:
            OA.extend([i]*k for i in range(N-x))

        if e1 == 0 and e2 == 0:
            raise RuntimeError("Brouwer's construction does not work for case v) with e2=e1=0")

        # Sets of size q+t
        OA.extend(OA_from_PBD(k,N,blocks_of_size_q_plus_t,check=False)[:-N])

        # Set of size x
        OA.extend([N-i-1 for i in R] for R in orthogonal_array(k,x))

    # vi)
    elif (t+q<x and x<q**2+1 and (e3 or e4) and # The result is wrong when e3=e4=0
          orthogonal_array(k   ,x    ,existence=True) and # d0
          orthogonal_array(k+e3,t    ,existence=True) and # d1-e3
          orthogonal_array(k+e4,t+1  ,existence=True) and # d2-e4
          orthogonal_array(k+1,t+q+1,existence=True)):    # d4-1
        if verbose:
            print "Case vi) with k={},q={},t={},x={},e3={},e4={}".format(k,q,t,x,e3,e4)

        OA = []

        # Sets of size t+1
        #
        # All x-(q+t) parallel classes with blocks of size t are extended with
        # x-(q+t) of the new points.
        if e4:
            assert x != q+t+1, "equivalent to e4=1"
            for i,classs in enumerate(partition_of_blocks_of_size_t[:x-(q+t)]):
                for B in classs:
                    B.append(N-x+i)
            OA.extend(OA_from_PBD(k,N,sum(partition_of_blocks_of_size_t[:x-(q+t)],[]),check=False)[:-N])
        else:
            assert x == q+t+1, "equivalent to e4=0"
            # Only one class. The OA(k,t+1) needs not be resolvable.
            OA.extend([B[xx] if xx<t else N-x for xx in R]
                       for R in incomplete_orthogonal_array(k,t+1,[1])
                       for B in partition_of_blocks_of_size_t[0])

        # Sets of size t
        if e3:
            assert x != q**2, "equivalent to e3=1"
            OA.extend(OA_from_PBD(k,N,sum(partition_of_blocks_of_size_t[x-(q+t):],[]),check=False)[:-N])
        else:
            assert x == q**2, "equivalent to e3=0"
            # Only one class. The OA(k,t) needs not be resolvable.
            OA.extend([B[xx] for xx in R]
                       for R in orthogonal_array(k,t)
                       for B in partition_of_blocks_of_size_t[-1])

        if e3 and e4:
            OA.extend([i]*k for i in range(N-x))

        elif e3 == 0 and e4 == 0:
            raise RuntimeError("Brouwer's construction does not work for case v) with e3=e4=0")

        # The sets of size q+t:
        #
        # We build an OA(k,t+q+1)-(t+q+1).OA(k,t+q+1) and the reordered
        # matrix. We then compute the product of every parallel class of the OA
        # (q+t classes in total) with the rows of the ordered matrix (extended
        # with the last q+t new points).

        # Resolvable OA(k,t+q+1)-(t+q+1).OA(k,t+q+1)
        OA_tq1 = incomplete_orthogonal_array(k,t+q+1,[1]*(t+q+1),resolvable=True)
        OA_tq1_classes = [OA_tq1[i*(t+q+1):(i+1)*(t+q+1)] for i in range(t+q)]

        blocks_of_size_q_plus_t = _reorder_matrix(blocks_of_size_q_plus_t)

        for i,classs in enumerate(OA_tq1_classes):
            OA.extend([R[xx] if xx<t+q else N-i-1 for xx in B]
                       for R in blocks_of_size_q_plus_t
                       for B in classs)

        # Set of size x
        OA.extend([N-xx-1 for xx in B] for B in orthogonal_array(k,x))

    else:
        raise ValueError("This input is not handled by Brouwer's result.")

    if check:
        assert is_orthogonal_array(OA,k,N,2,1)
    return OA
