# -*- coding: utf-8 -*-
"""
Shard intersection order

This file builds a combinatorial version of the shard intersection
order of type A. This is a lattice on the set of permutations, closely
related to noncrossing partitions and the weak order.

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

class ShardPosetElement:
    r"""
    An element of the shard poset.

    This is basically a permutation with extra stored arguments:

    - ``p`` - the permutation itself as a tuple
    - ``runs`` - the decreasing runs as a tuple of tuples
    - ``run_indices`` - a list ``integer -> index of the run``
    - ``dpg`` - the transitive closure of the shard preorder graph
    - ``spg`` - the transitive reduction of the shard preorder graph
    """
    def __init__(self, p):
        r"""
        INPUT:

        - ``p`` - a permutation
        """
        self.p = tuple(p)
        self.runs = p.decreasing_runs(as_tuple=True)
        self.run_indices = [None] * (len(p)+1)
        for i,bloc in enumerate(self.runs):
            for j in bloc:
                self.run_indices[j] = i
        G = shard_preorder_graph(self.runs)
        self.dpg = G.transitive_closure()
        self.spg = G.transitive_reduction()

    def __repr__(self):
        return repr(self.p)

    def __hash__(self):
        return hash(self.p)

    def __le__(self, other):
        """
        Comparison between two elements of the poset.

        This is the core function in the implementation of the
        shard intesection order.

        EXAMPLES::

            sage: from sage.combinat.shard_order import ShardPosetElement
            sage: p0 = Permutation([1,3,4,2])
            sage: p1 = Permutation([1,4,3,2])
            sage: e0 = ShardPosetElement(p0)
            sage: e1 = ShardPosetElement(p1)
            sage: e0 <= e1
            True
            sage: e1 <= e0
            False

            sage: p0 = Permutation([1,2,5,7,3,4,6,8])
            sage: p1 = Permutation([2,5,7,3,4,8,6,1])
            sage: e0 = ShardPosetElement(p0)
            sage: e1 = ShardPosetElement(p1)
            sage: e0 <= e1
            True
            sage: e1 <= e0
            False
        """
        if self.runs == other.runs:
            return True

        # r1 must have less runs than r0
        if len(other.runs) > len(self.runs):
            return False

        dico1 = other.run_indices

        # conversion: index of run in r0 -> index of run in r1
        dico0 = [None] * len(self.runs)
        for i, bloc in enumerate(self.runs):
            j0 = dico1[bloc[0]]
            for k in bloc:
                if dico1[k] != j0:
                    return False
            dico0[i] = j0

        # at this point, the set partitions given by tuples are comparable
        dg0 = self.spg
        dg1 = other.dpg

        for i,j in dg0.edge_iterator(labels=False):
            if dico0[i] != dico0[j] and not dg1.has_edge(dico0[i], dico0[j]):
                return False
        return True


def shard_preorder_graph(runs):
    """
    Return the preorder attached to a tuple of decreasing runs.

    This is a directed graph, whose vertices correspond to the runs.

    This only depends on the initial and final indices of the runs.
    For this reason, this input is given in that shorten way.

    INPUT:

    a tuple of pairs (i,j), each one standing for a run from `i` to `j`.

    EXAMPLES::

        sage: from sage.combinat.shard_order import shard_preorder_graph
        sage: s = Permutation([2,8,3,9,6,4,5,1,7])
        sage: def cut(lr):
        ....:     return tuple((r[0], r[-1]) for r in lr)
        sage: shard_preorder_graph(cut(s.decreasing_runs()))
        Digraph on 5 vertices
    """
    N = len(runs)
    dg = DiGraph(N)
    dg.add_edges((i,j) for i in range(N - 1) \
                       for j in range(i + 1, N) \
                       if runs[i][-1] < runs[j][0] and runs[j][-1] < runs[i][0])
    return dg


def shard_compares(r0, r1):
    """
    Comparison between two tuples of decreasing runs.

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
    # r1 must have less runs than r0
    if len(r1) > len(r0):
        return False

    # conversion: integer -> index of run in r1
    dico1 = {j: i for i, bloc in enumerate(r1) for j in bloc}

    dico0 = {}
    # conversion: index of run in r0 -> index of run in r1
    for i, bloc in enumerate(r0):
        ind = set([dico1[k] for k in bloc])
        if len(ind) == 1:
            dico0[i] = ind.pop()
        else:
            return False
    # at this point, the set partitions given by tuples are comparable

    cr0 = tuple((r[0], r[-1]) for r in r0)
    cr1 = tuple((r[0], r[-1]) for r in r1)
    dg0 = shard_preorder_graph(cr0)
    dg1 = shard_preorder_graph(cr1)
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
    import operator
    Sn = [ShardPosetElement(s) for s in Permutations(n)]
    return Poset([Sn, operator.le], cover_relations=False, facade=True)
