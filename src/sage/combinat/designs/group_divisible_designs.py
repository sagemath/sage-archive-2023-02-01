r"""
(GDD) Group-Divisible Designs

This module gathers everything related to Group-Divisible Designs. The
constructions defined here can be accessed through ``designs.<tab>``::

    sage: designs.group_divisible_design(14,{4},{2})
    Group Divisible Design on 14 points of type 2^7

The main function implemented here is :meth:`group_divisible_design`, which
calls all others. The following functions are available:

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :func:`group_divisible_design` | Return a `(v,K,G)`-Group Divisible Design.
    :func:`GDD_4_2` | Return a `(2q,\{4\},\{2\})`-GDD for `q` a prime power with `q\equiv 1\pmod{6}`.

Functions
---------
"""
from incidence_structures import GroupDivisibleDesign
from sage.rings.arith     import is_prime_power
from sage.misc.unknown    import Unknown

def group_divisible_design(v,K,G,existence=False,check=False):
    r"""
    Return a `(v,K,G)`-Group Divisible Design.

    A `(v,K,G)`-GDD is a pair `(\mathcal G, \mathcal B)` where:

    - `\mathcal G` is a partition of `X=\bigcup \mathcal G` where `|X|=v`

    - `\forall S\in \mathcal G, |S| \in G`

    - `\forall S\in \mathcal B, |S| \in K`

    - `\mathcal G\cup \mathcal B` is a `(v,K\cup G)`-PBD

    For more information, see the documentation of
    :class:`~sage.combinat.designs.incidence_structures.GroupDivisibleDesign` or
    :class:`~sage.combinat.designs.bibd.PairwiseBalancedDesign`.

    INPUT:

    - ``v`` (integer)

    - ``K,G`` (sets of integers)

    - ``existence`` (boolean) -- instead of building the design, return:

        - ``True`` -- meaning that Sage knows how to build the design

        - ``Unknown`` -- meaning that Sage does not know how to build the
          design, but that the design may exist (see :mod:`sage.misc.unknown`).

        - ``False`` -- meaning that the design does not exist.

    - ``check`` -- (boolean) Whether to check that output is correct before
      returning it. As this is expected to be useless (but we are cautious
      guys), you may want to disable it whenever you want speed. Set to ``True``
      by default.

    .. NOTE::

        The GDD returned by this function are defined on ``range(v)``, and its
        groups are sets of consecutive integers.

    EXAMPLES::

        sage: designs.group_divisible_design(14,{4},{2})
        Group Divisible Design on 14 points of type 2^7
    """
    G = list(set(G))
    K = list(set(K))

    blocks = None

    # from a (v+1,k,1)-BIBD
    if (len(G) == 1 and
        len(K) == 1 and
        G[0]+1 in K):
        if existence:
            return balanced_incomplete_block_design(v+1,K[0],existence=True)
        BIBD = balanced_incomplete_block_design(v+1,K[0],check=False)
        groups = [[x for x in S if x!=v] for S in BIBD if v in S]
        d = {p:i for i,p in enumerate(sum(groups,[]))}
        d[v]=v
        BIBD.relabel(d)
        groups = [range((k-1)*i,(k-1)*(i+1)) for i in range(v//(k-1))]
        blocks = [S for S in BIBD if v not in S]

    # (v,{4},{2})-GDD
    elif (v%2==0   and
          K == [4] and
          G == [2] and
          GDD_4_2(v//2,existence=True)):
        if existence:
            return True
        return GDD_4_2(v//2,check=check)

    if blocks:
        return GroupDivisibleDesign(v,
                                    groups = groups,
                                    blocks = blocks,
                                    G = G,
                                    K = K,
                                    check = check,
                                    copy  = True)

    if existence:
        return Unknown
    raise NotImplementedError

def GDD_4_2(q,existence=False,check=True):
    r"""
    Return a `(2q,\{4\},\{2\})`-GDD for `q` a prime power with `q\equiv 1\pmod{6}`.

    This method implements Lemma VII.5.17 from [BJL99] (p.495).

    INPUT:

    - ``q`` (integer)

    - ``existence`` (boolean) -- instead of building the design, return:

        - ``True`` -- meaning that Sage knows how to build the design

        - ``Unknown`` -- meaning that Sage does not know how to build the
          design, but that the design may exist (see :mod:`sage.misc.unknown`).

        - ``False`` -- meaning that the design does not exist.

    - ``check`` -- (boolean) Whether to check that output is correct before
      returning it. As this is expected to be useless (but we are cautious
      guys), you may want to disable it whenever you want speed. Set to ``True``
      by default.

    EXAMPLE::

        sage: from sage.combinat.designs.group_divisible_designs import GDD_4_2
        sage: GDD_4_2(7,existence=True)
        True
        sage: GDD_4_2(7)
        Group Divisible Design on 14 points of type 2^7
        sage: GDD_4_2(8,existence=True)
        Unknown
        sage: GDD_4_2(8)
        Traceback (most recent call last):
        ...
        NotImplementedError
    """
    if q <=1 or q%6 != 1 or not is_prime_power(q):
        if existence:
            return Unknown
        raise NotImplementedError
    if existence:
        return True

    from sage.rings.finite_rings.constructor import FiniteField as GF
    G = GF(q,'x')
    w = G.primitive_element()
    e = w**((q-1)/3)

    # A first parallel class is defined. G acts on it, which yields all others.
    first_class = [[(0,0),(1,w**i),(1,e*w**i),(1,e*e*w**i)] for i in range((q-1)/6)]

    label = {p:i for i,p in enumerate(G)}
    classes = [[[2*label[x[1]+g]+(x[0]+j)%2 for x in S]
                for S in first_class]
               for g in G for j in range(2)]

    return GroupDivisibleDesign(2*q,
                                groups = [[i,i+1] for i in range(0,2*q,2)],
                                blocks = sum(classes,[]),
                                K      = [4],
                                G      = [2],
                                check  = check,
                                copy   = False)
