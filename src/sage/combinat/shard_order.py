# -*- coding: utf-8 -*-
"""
Shard intersection order

This file builds a combinatorial version of the shard intersection
order of type A. This is a lattice on the set of permutations, closely
related to noncrossing partitions.

REFERENCES:

.. [Banc2011] E. E. Bancroft, *Shard Intersections and Cambrian Congruence
   Classes in Type A.*, Ph.D. Thesis, North Carolina State University. 2011.

.. [Pete2013] T. Kyle Petersen, *On the shard intersection order of
   a Coxeter group*, SIAM J. Discrete Math. 27 (2013), no. 4, 1880-1912. 

.. [Read2011] N. Reading, *Noncrossing partitions and the shard intersection
   order*, J. Algebraic Combin., 33 (2011), 483-530.
"""
from sage.combinat.posets.posets import Poset
from sage.graphs.digraph import DiGraph
from sage.combinat.permutation import Permutations
from sage.misc.cachefunc import cached_function
from sage.rings.infinity import Infinity


@cached_function
def shard_preorder_graph(runs):
    """
    Return the preorder attached to a tuple of decreasing runs.

    This is a directed graph, whose vertices correspond to the runs.

    EXAMPLES::

        sage: from sage.combinat.shard_order import shard_preorder_graph
        sage: s = Permutation([2,8,3,9,6,4,5,1,7])
        sage: shard_preorder_graph(s.decreasing_runs())
        Digraph on 5 vertices
    """
    N = len(runs)
    dg = DiGraph()
    dg.add_vertices(range(N))
    for i in range(N - 1):
        for j in range(i + 1, N):
            if not(runs[i][-1] > runs[j][0] or runs[j][-1] > runs[i][0]):
                dg.add_edge((i, j))
    dg_red = dg.transitive_reduction()
    return dg_red


def shard_compares(r0, r1):
    """
    Comparison between two tuples of decreasing runs.

    This is the core function in the implementation of the
    shard intesection order.

    EXAMPLES::

        sage: from sage.combinat.shard_order import shard_compares
        sage: s0 = Permutation([1,2,5,7,3,4,6,8])
        sage: s1 = Permutation([2,5,7,3,4,8,6,1])
        sage: shard_compares(s0.decreasing_runs(),s1.decreasing_runs())
        True
        sage: shard_compares(s1.decreasing_runs(),s0.decreasing_runs())
        False
    """
    # We assume that r1 has less runs than r0

    # conversion: integer -> index of run in r1
    dico1 = {j: i  for i, bloc in enumerate(r1) for j in bloc}

    dico0 = {}
    # conversion: index of run in r0 -> index of run in r1
    for i, bloc in enumerate(r0):
        ind = set([dico1[j] for j in bloc])
        if len(ind) == 1:
            dico0[i] = ind.pop()
        else:
            return False
    # at this point, the set partitions given by tuples are comparable

    dg0 = shard_preorder_graph(r0)
    dg1 = shard_preorder_graph(r1)
    for deb, fin, _ in dg0.edges():
        ideb = dico0[deb]
        ifin = dico0[fin]
        if dg1.distance(ideb, ifin) == Infinity:
            return False

    return True


def shard_poset(n):
    """
    Return the shard intersection order on permutations of size `n`.

    This is defined on the set of permutations. To every permutation,
    one can attach a pre-order, using the descending runs and their
    relative positions.

    The shard intersection order is given by the implication (or refinement)
    order on the set of pre-orders defined from all permutations.

    .. SEEALSO::

        :func:`shard_preorder_graph`

    EXAMPLES::

        sage: P = posets.ShardPoset(4); P
        Finite poset containing 24 elements
        sage: len(P.maximal_chains())
        34
        sage: P.is_selfdual()
        False
    """
    Sn = [s.decreasing_runs() for s in Permutations(n)]
    return Poset([Sn, shard_compares], cover_relations=False)
