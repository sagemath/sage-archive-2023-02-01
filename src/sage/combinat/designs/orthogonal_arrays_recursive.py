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
        54
    """
    assert k > 3

    for find_c in [find_product_decomposition,
                   find_wilson_decomposition_with_one_truncated_group,
                   find_wilson_decomposition_with_two_truncated_groups,
                   find_construction_3_3,
                   find_construction_3_4,
                   find_construction_3_5,
                   find_construction_3_6,
                   find_q_x]:
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
    such that Sage knows how to build a `OA(k,m), OA(k,m+1),OA(k+1,t)` and a
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
      `OA(k,n),OA(k,m),OA(k,m+1),OA(k,r)`.

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

      The following designs must be available:
      `OA(k,n),OA(k,m),OA(k,m+1),OA(k,m+2),OA(k,s)`. Additionnally, it requires
      either a `OA(k,m+r)` or a `OA(k,m+r+1)`.

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
        raise Exception

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
    need a `OA(k,m+0)` but only `OA(k,m+1),OA(k,m+2),OA(k,m+3)`.

    This is construction 3.5 from [AC07]_.

    INPUT:

    - ``k,n,m`` (integers)

    - ``r,s,t`` (integers) -- sizes of the three truncated groups,
      such that `r\leq s` and `(q-r-1)(q-s) \geq (q-s-1)*(q-r)`.

    The following designs must be available :
    `OA(k,n),OA(k,r),OA(k,s),OA(k,t),OA(k,m+1),OA(k,m+2),OA(k,m+3)`.

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
      must be available: `OA(k+r,q),OA(k,m),OA(k,m+1),OA(k,m+2)`.

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
