# -*- coding: utf-8 -*-
r"""
Intersection graphs

The methods defined here appear in :mod:`sage.graphs.graph_generators`.
"""

###########################################################################
#
#           Copyright (C) 2006 Robert L. Miller <rlmillster@gmail.com>
#                              and Emily A. Kirkman
#           Copyright (C) 2009 Michael C. Yurko <myurko@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
###########################################################################

# import from Sage library
from sage.graphs.graph import Graph

def IntervalGraph(intervals, points_ordered = False):
    r"""
    Returns the graph corresponding to the given intervals.

    An interval graph is built from a list `(a_i,b_i)_{1\leq i \leq n}` of
    intervals : to each interval of the list is associated one vertex, two
    vertices being adjacent if the two corresponding (closed) intervals
    intersect.

    INPUT:

    - ``intervals`` -- the list of pairs `(a_i,b_i)` defining the graph.

    - ``points_ordered`` -- states whether every interval `(a_i,b_i)` of
      `intervals` satisfies `a_i<b_i`. If satisfied then setting
      ``points_ordered`` to ``True`` will speed up the creation of the graph.

    .. NOTE::

        * The vertices are named 0, 1, 2, and so on. The intervals used
          to create the graph are saved with the graph and can be recovered
          using ``get_vertex()`` or ``get_vertices()``.

    EXAMPLE:

    The following line creates the sequence of intervals
    `(i, i+2)` for i in `[0, ..., 8]`::

        sage: intervals = [(i,i+2) for i in range(9)]

    In the corresponding graph ::

        sage: g = graphs.IntervalGraph(intervals)
        sage: g.get_vertex(3)
        (3, 5)
        sage: neigh = g.neighbors(3)
        sage: for v in neigh: print g.get_vertex(v)
        (1, 3)
        (2, 4)
        (4, 6)
        (5, 7)

    The is_interval() method verifies that this graph is an interval graph. ::

        sage: g.is_interval()
        True

    The intervals in the list need not be distinct. ::

        sage: intervals = [ (1,2), (1,2), (1,2), (2,3), (3,4) ]
        sage: g = graphs.IntervalGraph(intervals,True)
        sage: g.clique_maximum()
        [0, 1, 2, 3]
        sage: g.get_vertices()
        {0: (1, 2), 1: (1, 2), 2: (1, 2), 3: (2, 3), 4: (3, 4)}

    The endpoints of the intervals are not ordered we get the same graph
    (except for the vertex labels). ::

        sage: rev_intervals = [ (2,1), (2,1), (2,1), (3,2), (4,3) ]
        sage: h = graphs.IntervalGraph(rev_intervals,False)
        sage: h.get_vertices()
        {0: (2, 1), 1: (2, 1), 2: (2, 1), 3: (3, 2), 4: (4, 3)}
        sage: g.edges() == h.edges()
        True
    """

    n = len(intervals)
    g = Graph(n)

    if points_ordered:
        for i in xrange(n-1):
            li,ri = intervals[i]
            for j in xrange(i+1,n):
                lj,rj = intervals[j]
                if ri < lj or rj < li: continue
                g.add_edge(i,j)
    else:
        for i in xrange(n-1):
            I = intervals[i]
            for j in xrange(i+1,n):
                J = intervals[j]
                if max(I) < min(J) or max(J) < min(I): continue
                g.add_edge(i,j)

    rep = dict( zip(range(n),intervals) )
    g.set_vertices(rep)

    return g

def PermutationGraph(second_permutation, first_permutation = None):
    r"""
    Builds a permutation graph from one (or two) permutations.

    General definition

    A Permutation Graph can be encoded by a permutation `\sigma`
    of `1, ..., n`. It is then built in the following way :

      Take two horizontal lines in the euclidean plane, and mark points `1, ...,
      n` from left to right on the first of them. On the second one, still from
      left to right, mark point in the order in which they appear in `\sigma`.
      Now, link by a segment the two points marked with 1, then link together
      the points marked with 2, and so on. The permutation graph defined by the
      permutation is the intersection graph of those segments : there exists a
      point in this graph for each element from `1` to `n`, two vertices `i, j`
      being adjacent if the segments `i` and `j` cross each other.

    The set of edges of the resulting graph is equal to the set of inversions of
    the inverse of the given permutation.

    INPUT:

    - ``second_permutation`` -- the permutation from which the graph should be
      built. It corresponds to the ordering of the elements on the second line
      (see previous definition)

    - ``first_permutation`` (optional) -- the ordering of the elements on the
      *first* line. This is useful when the elements have no natural ordering,
      for instance when they are strings, or tuples, or anything else.

      When ``first_permutation == None`` (default), it is set to be equal to
      ``sorted(second_permutation)``, which just yields the expected
      ordering when the elements of the graph are integers.

    .. SEEALSO:

      - Recognition of Permutation graphs in the :mod:`comparability module
        <sage.graphs.comparability>`.

      - Drawings of permutation graphs as intersection graphs of segments is
        possible through the
        :meth:`~sage.combinat.permutation.Permutation.show` method of
        :class:`~sage.combinat.permutation.Permutation` objects.

        The correct argument to use in this case is ``show(representation =
        "braid")``.

      - :meth:`~sage.combinat.permutation.Permutation.inversions`

    EXAMPLE::

        sage: p = Permutations(5).random_element()
        sage: edges = graphs.PermutationGraph(p).edges(labels =False)
        sage: set(edges) == set(p.inverse().inversions())
        True

    TESTS::

        sage: graphs.PermutationGraph([1, 2, 3], [4, 5, 6])
        Traceback (most recent call last):
        ...
        ValueError: The two permutations do not contain the same set of elements ...
    """
    if first_permutation == None:
        first_permutation = sorted(second_permutation)
    else:
        if set(second_permutation) != set(first_permutation):
            raise ValueError("The two permutations do not contain the same "+
                             "set of elements ! It is going to be pretty "+
                             "hard to define a permutation graph from that !")

    vertex_to_index = {}
    for i, v in enumerate(first_permutation):
        vertex_to_index[v] = i+1

    from sage.combinat.permutation import Permutation
    p2 = Permutation(map(lambda x:vertex_to_index[x], second_permutation))
    p1 = Permutation(map(lambda x:vertex_to_index[x], first_permutation))
    p2 = p2 * p1.inverse()
    p2 = p2.inverse()

    g = Graph(name="Permutation graph for "+str(second_permutation))
    g.add_vertices(second_permutation)

    for u,v in p2.inversions():
        g.add_edge(first_permutation[u-1], first_permutation[v-1])

    return g

def ToleranceGraph(tolrep):
    r"""
    Returns the graph generated by the tolerance representation ``tolrep``.

    The tolerance representation ``tolrep`` is described by the list
    `((l_0,r_0,t_0), (l_1,r_1,t_1), ..., (l_k,r_k,t_k))` where `I_i = (l_i,r_i)`
    denotes a closed interval on the real line with `l_i < r_i` and `t_i` a
    positive value, called tolerance. This representation generates the
    tolerance graph with the vertex set {0,1, ..., k} and the edge set `{(i,j):
    |I_i \cap I_j| \ge \min{t_i, t_j}}` where `|I_i \cap I_j|` denotes the
    length of the intersection of `I_i` and `I_j`.

    INPUT:

    - ``tolrep`` -- list of triples `(l_i,r_i,t_i)` where `(l_i,r_i)` denotes a
      closed interval on the real line and `t_i` a positive value.

    .. NOTE::

        The vertices are named 0, 1, ..., k. The tolerance representation used
        to create the graph is saved with the graph and can be recovered using
        ``get_vertex()`` or ``get_vertices()``.

    EXAMPLE:

    The following code creates a tolerance representation ``tolrep``, generates
    its tolerance graph ``g``, and applies some checks::

        sage: tolrep = [(1,4,3),(1,2,1),(2,3,1),(0,3,3)]
        sage: g = graphs.ToleranceGraph(tolrep)
        sage: g.get_vertex(3)
        (0, 3, 3)
        sage: neigh = g.neighbors(3)
        sage: for v in neigh: print g.get_vertex(v)
        (1, 2, 1)
        (2, 3, 1)
        sage: g.is_interval()
        False
        sage: g.is_weakly_chordal()
        True

    The intervals in the list need not be distinct ::

        sage: tolrep2 = [(0,4,5),(1,2,1),(2,3,1),(0,4,5)]
        sage: g2 = graphs.ToleranceGraph(tolrep2)
        sage: g2.get_vertices()
        {0: (0, 4, 5), 1: (1, 2, 1), 2: (2, 3, 1), 3: (0, 4, 5)}
        sage: g2.is_isomorphic(g)
        True

    Real values are also allowed ::

        sage: tolrep = [(0.1,3.3,4.4),(1.1,2.5,1.1),(1.4,4.4,3.3)]
        sage: g = graphs.ToleranceGraph(tolrep)
        sage: g.is_isomorphic(graphs.PathGraph(3))
        True

    TEST:

    Giving negative third value::

        sage: tolrep = [(0.1,3.3,-4.4),(1.1,2.5,1.1),(1.4,4.4,3.3)]
        sage: g = graphs.ToleranceGraph(tolrep)
        Traceback (most recent call last):
        ...
        ValueError: Invalid tolerance representation at position 0; third value must be positive!
    """
    n = len(tolrep)

    for i in xrange(n):
        if tolrep[i][2] <= 0:
            raise ValueError("Invalid tolerance representation at position "+str(i)+"; third value must be positive!")

    g = Graph(n)

    for i in xrange(n-1):
        li,ri,ti = tolrep[i]
        for j in xrange(i+1,n):
            lj,rj,tj = tolrep[j]
            if min(ri,rj) - max(li,lj) >= min(ti,tj):
                g.add_edge(i,j)

    rep = dict( zip(range(n),tolrep) )
    g.set_vertices(rep)

    return g
