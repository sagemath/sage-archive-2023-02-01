r"""
Orthogonal arrays (Recursive constructions)

This module implements several functions to find recursive constructions of
:mod:`Orthogonal Arrays <sage.combinat.designs.orthogonal_arrays>` using
Wilson's construction. To this end, they compute and truncate OA in specific
ways.

The main function of this module, i.e. :func:`find_recursive_construction`,
queries all implemented recursive constructions of designs. It is used by
Sage's function
:func:`~sage.combinat.designs.orthogonal_arrays.orthogonal_array`.

Functions
---------
"""
from sage.misc.cachefunc import cached_function
from sage.misc.unknown import Unknown
from orthogonal_arrays import orthogonal_array, wilson_construction
from sage.combinat.designs.orthogonal_arrays import is_orthogonal_array

def find_recursive_construction(k,n,existence=False):
    r"""
    Find a recursive construction of a `OA(k,n)`

    This method attempts to build an `OA(k,n)` through the following
    constructions:

    - :func:`construction_3_3`
    - :func:`construction_3_4`
    - :func:`construction_3_5`
    - :func:`construction_3_6`

    INPUT:

    - ``k,n`` (integers)

    - ``existence`` (boolean) -- whether to return the actual OA, or only
      indicate whether it can be built.

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays_recursive import find_recursive_construction
        sage: from sage.combinat.designs.orthogonal_arrays import is_orthogonal_array
        sage: count = 0
        sage: for n in range(10,150):
        ....:     k = designs.orthogonal_array(None,n,existence=True)
        ....:     if find_recursive_construction(k,n,existence=True):
        ....:         count = count + 1
        ....:         OA = find_recursive_construction(k,n,existence=False)
        ....:         assert is_orthogonal_array(OA,k,n,2,verbose=True)
        sage: print count
        40
    """
    for find_c in [find_construction_3_3,
                   find_construction_3_4,
                   find_construction_3_5,
                   find_construction_3_6]:

        if find_c(k,n):
            if existence:
                return True
            else:
                f,args = find_c(k,n)
                return f(*args)

    return Unknown

@cached_function
def find_construction_3_3(k,n):
    r"""
    Finds a decomposition for construction 3.3

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
    for mm in range(2,n//2+1):
        if (not orthogonal_array( k , mm , existence=True) or
            not orthogonal_array( k ,mm+1, existence=True)):
            continue

        for nn in range(2,n//mm+1):
            i = n-nn*mm
            if i<=0:
                continue

            if (orthogonal_array(k+i, nn , existence=True) and
                orthogonal_array( k ,mm+i, existence=True)):
                return construction_3_3, (k,nn,mm,i)

def construction_3_3(k,n,m,i):
    r"""
    Returns an `OA(k,nm+i)`.

    This is Wilson's construction with `i` truncated columns of size `1`, such
    that a block `B_0` of the incomplete OA intersects all truncated
    columns. There is a slight difference however, as the block `B_0` is only
    considered up to its first `k` coordinates, and a `OA(k,i)` is used instead
    of `iOA(k,1)` (thus there is no need for an `OA(k,m+i)`.

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
    from orthogonal_arrays import OA_relabel, incomplete_orthogonal_array
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

@cached_function
def find_construction_3_4(k,n):
    r"""
    Finds a decomposition for construction 3.4

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
    for mm in range(2,n//2+1):
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
    Returns a `OA(k,nm+rs)`.

    This is Wilson's construction applied to an incomplete `OA(k+r+1,n)` with
    `k` columns of size 1 and a column of size `s`.

    The the unique elements of the `k` columns are picked so that a block `B_0`
    contains them all.

    - If there exists an `OA(k,m+r+1)` the column of size `s` is truncated in
      order to intersect `B_0`.

    - If there exists an `OA(k,m+r+1)`, the last column must not intersect `B_0`

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
    from orthogonal_arrays import OA_relabel
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

@cached_function
def find_construction_3_5(k,n):
    r"""
    Finds a decomposition for construction 3.5

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
    from sage.combinat.composition import Compositions

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

            for r,s,t in Compositions(i+3,length=3,max_part=nn):
                # avoid warning on Compositions(...,min_part=0)
                r = r-1
                s = s-1
                t = t-1
                if (r <= s and
                    (nn-r-1)*(nn-s) < t and
                    (r==0 or orthogonal_array(k,r,existence=True)) and
                    (s==0 or orthogonal_array(k,s,existence=True)) and
                    (t==0 or orthogonal_array(k,t,existence=True))):
                    return construction_3_5, (k,nn,mm,r,s,t)

def construction_3_5(k,n,m,r,s,t):
    r"""
    Returns an `OA(k,nm+r+s+t)`.

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
    from orthogonal_arrays import OA_relabel
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

@cached_function
def find_construction_3_6(k,n):
    r"""
    Finds a decomposition for construction 3.6

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

    for mm in range(2,n//2+1):
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
    Returns a `OA(k,nm+i)`

    This is Wilson's construction with `r` columns of order `1`, in which every
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
    OA = OA_and_oval(n)
    OA = [B[:k+i] for B in OA]
    OA = [B[:k] + [x if x==0 else None for x in B[k:]] for B in OA]
    OA = wilson_construction(OA,k,n,m,i,[1]*i)
    assert is_orthogonal_array(OA,k,n*m+i)
    return OA

def OA_and_oval(q):
    r"""
    Returns a `OA(q+1,q)` whose blocks contains `\leq 2` zeroes in the last `q`
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

    REFERENCES:

    .. [AC07] Concerning eight mutually orthogonal latin squares
      Julian R. Abel, Nicholas Cavenagh
      Journal of Combinatorial Designs
      Vol. 15, n.3, pp. 255-261
      2007
    """
    from sage.rings.arith import is_prime_power
    from sage.combinat.designs.block_design import projective_plane
    from orthogonal_arrays import OA_relabel

    assert is_prime_power(q)
    B = projective_plane(q)

    # We compute the oval with a linear program
    from sage.numerical.mip import MixedIntegerLinearProgram
    p = MixedIntegerLinearProgram()
    b = p.new_variable(binary=True)
    V = B.points()
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
