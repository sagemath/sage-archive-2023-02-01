# distutils: libraries = gmp
r"""
Cutwidth

This module implements several algorithms to compute the cutwidth of a graph and
the corresponding ordering of the vertices. It also implements tests functions
for evaluation the width of a linear ordering (or layout).

Given an ordering
`v_1,\cdots, v_n` of the vertices of `V(G)`, its *cost* is defined as:

.. MATH::

    c(v_1, ..., v_n) = \max_{1\leq i \leq n-1} c'(\{v_1, ..., v_i\})

Where

.. MATH::

    c'(S) = |\{(u,w)\in E(G)\mid u\in S\text{ and }w\in V(G)\backslash S\}|

The *cutwidth* of a graph `G` is equal to the minimum cost of an ordering of its
vertices.


**This module contains the following methods**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`cutwidth` | Return the cutwidth of the graph and the corresponding vertex ordering.
    :meth:`cutwidth_dyn` | Compute the cutwidth of `G` using an exponential time and space algorithm based on dynamic programming
    :meth:`width_of_cut_decomposition` | Return the width of the cut decomposition induced by the linear ordering `L` of the vertices of `G`


Exponential algorithm for cutwidth
----------------------------------

This algorithm differs from
:meth:`sage.graphs.graph_decompositions.vertex_separation.vertex_separation_exp`
only in the cost function. See the :mod:`vertex separation module's
documentation <sage.graphs.graph_decompositions.vertex_separation>` for more
details on this algorithm.

.. NOTE::

    Because of its current implementation, this algorithm only works on graphs
    on strictly less than 32 vertices. This can be changed to 64 if necessary,
    but 32 vertices already require 4GB of memory.


Authors
-------

- David Coudert (2015-06): Initial version


Methods
-------
"""
#*****************************************************************************
#          Copyright (C) 2015 David Coudert <david.coudert@inria.fr>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#      as published by the Free Software Foundation; either version 2 of
#              the License, or (at your option) any later version.
#                        http://www.gnu.org/licenses/
#*****************************************************************************

include 'sage/ext/stdsage.pxi'
include 'sage/ext/interrupt.pxi'
include 'fast_digraph.pyx'
from libc.stdint cimport uint8_t
from sage.ext.memory cimport check_allocarray
from sage.rings.integer_ring import ZZ


################################################################################
# Function for testing the validity of a linear vertex ordering
#
# A linear ordering `L` of the vertices of a graph `G` is valid if all vertices
# of `G` are in `L`, and if `L` contains no other vertex and no duplicated
# vertices.
################################################################################
from sage.graphs.graph_decompositions.vertex_separation import is_valid_ordering


################################################################################
# Measurement function of the width of some layout
################################################################################

def width_of_cut_decomposition(G, L):
    r"""
    Returns the width of the cut decomposition induced by the linear ordering
    `L` of the vertices of `G`.

    If `G` is an instance of :mod:`Graph <sage.graphs.graph>`, this function
    returns the width `cw_L(G)` of the cut decomposition induced by the linear
    ordering `L` of the vertices of `G`.

    .. MATH::

        cw_L(G) =  \max_{0\leq i< |V|-1} |\{(u,w)\in E(G)\mid u\in L[:i]\text{ and }w\in V(G)\setminus L[:i]\}|

    INPUT:

    - ``G`` -- a Graph

    - ``L`` -- a linear ordering of the vertices of ``G``

    EXAMPLES:

    Cut decomposition of a Cycle graph::

        sage: from sage.graphs.graph_decompositions import cutwidth
        sage: G = graphs.CycleGraph(6)
        sage: L = G.vertices()
        sage: cutwidth.width_of_cut_decomposition(G, L)
        2

    Cut decomposition of a Path graph::

        sage: from sage.graphs.graph_decompositions import cutwidth
        sage: P = graphs.PathGraph(6)
        sage: cutwidth.width_of_cut_decomposition(P, [0, 1, 2, 3, 4, 5])
        1
        sage: cutwidth.width_of_cut_decomposition(P, [5, 0, 1, 2, 3, 4])
        2
        sage: cutwidth.width_of_cut_decomposition(P, [0, 2, 4, 1, 3, 5])
        5

    TESTS:

    Giving a wrong linear ordering::

        sage: from sage.graphs.graph_decompositions import cutwidth
        sage: cutwidth.width_of_cut_decomposition(Graph(), ['a','b'])
        Traceback (most recent call last):
        ...
        ValueError: The input linear vertex ordering L is not valid for G.
    """
    if not is_valid_ordering(G, L):
        raise ValueError("The input linear vertex ordering L is not valid for G.")

    position = dict( (u,i) for i,u in enumerate(L) )

    cpt = [0]*len(L)
    for u,v in G.edge_iterator(labels=None):
        x,y = position[u],position[v]
        if x>y:
            x,y = y,x
        for i in range(x,y):
            cpt[i] += 1

    return max(cpt)




################################################################################
# Front end method for cutwidth
################################################################################

def cutwidth(G, algorithm="exponential", cut_off=0, verbose=False):
    r"""
    Return the cutwidth of the graph and the corresponding vertex ordering.

    INPUT:

    - ``G`` -- a Graph or a DiGraph

    - ``algorithm`` -- (default: ``"exponential"``) Specify the algorithm to use
      among

      - ``exponential`` -- Use an exponential time and space algorithm based on
        dynamic programming. This algorithm only works on graphs with strictly
        less than 32 vertices.

    - ``cut_off`` -- (default: 0) This parameter is used to stop the search as
      soon as a solution with width at most ``cut_off`` is found, if any. If
      this bound cannot be reached, the best solution found is returned.

    - ``verbose`` (boolean) -- whether to display information on the
      computations.

    OUTPUT:

    A pair ``(cost, ordering)`` representing the optimal ordering of the
    vertices and its cost.

    EXAMPLES:

    Cutwidth of a Complete Graph::

        sage: from sage.graphs.graph_decompositions.cutwidth import cutwidth
        sage: G = graphs.CompleteGraph(5)
        sage: cw,L = cutwidth(G, algorithm="exponential"); cw
        6
        sage: K = graphs.CompleteGraph(6)
        sage: cw,L = cutwidth(K, algorithm="exponential"); cw
        9
        sage: cw,L = cutwidth(K+K, algorithm="exponential"); cw
        9

    The cutwidth of a `p\times q` Grid Graph with `p\leq q` is `p+1`::

        sage: from sage.graphs.graph_decompositions.cutwidth import cutwidth
        sage: G = graphs.Grid2dGraph(3,3)
        sage: cw,L = cutwidth(G, algorithm="exponential"); cw
        4
        sage: G = graphs.Grid2dGraph(3,5)
        sage: cw,L = cutwidth(G, algorithm="exponential"); cw
        4

    TESTS:

    Given a wrong algorithm::

        sage: from sage.graphs.graph_decompositions.cutwidth import cutwidth
        sage: cutwidth(Graph(), algorithm="SuperFast")
        Traceback (most recent call last):
        ...
        ValueError: Algorithm "SuperFast" has not been implemented yet. Please contribute.

    Given anything else than a Graph::

        sage: from sage.graphs.graph_decompositions.cutwidth import cutwidth
        sage: cutwidth(range(4))
        Traceback (most recent call last):
        ...
        ValueError: The parameter must be a Graph.

    Giving a wrong type cut off::

        sage: from sage.graphs.graph_decompositions.cutwidth import cutwidth
        sage: cutwidth(Graph(), cut_off='toto')
        Traceback (most recent call last):
        ...
        ValueError: The specified cut off parameter must be an integer.
    """
    from sage.graphs.graph import Graph

    CC = []
    if isinstance(G, Graph):
        if not G.is_connected():
            # We decompose the graph into connected components.
            CC = G.connected_components()

    else:
        raise ValueError('The parameter must be a Graph.')

    if not cut_off in ZZ:
        raise ValueError("The specified cut off parameter must be an integer.")

    if CC:
        # The graph has several connected components. We solve the problem on
        # each of them and order partial solutions in the same order than in
        # list CC. The cutwidth is the maximum over all these subgraphs.
        cw, L = 0, []
        for V in CC:

            if len(V)==1:
                # We can directly add this vertex to the solution
                L.extend(V)

            else:
                # We build the connected subgraph and do a recursive call to get
                # its cutwidth and corresponding ordering
                H = G.subgraph(V)
                cwH,LH = cutwidth(H, algorithm = algorithm,
                                  cut_off      = cut_off,
                                  verbose      = verbose)

                # We update the cutwidth and ordering
                cw = max(cw, cwH)
                L.extend(LH)

        return cw, L


    # We have a (strongly) connected graph and we call the desired algorithm
    if algorithm == "exponential":
        return cutwidth_dyn(G)

    else:
        raise ValueError('Algorithm "{}" has not been implemented yet. Please contribute.'.format(algorithm))



################################################################################
# Dynamic Programming algorithm for cutwidth
################################################################################

def cutwidth_dyn(G, lower_bound=0):
    """
    Dynamic programming algorithm for the cutwidth of a Graph.

    This function uses dynamic programming algorithm for determining an optimal
    layout for the cutwidth of `G`. See the :mod:`module's documentation
    <sage.graphs.graph_decompositions.cutwidth>` for more details on this
    method.

    INPUT:

    - ``G`` -- a Graph

    - ``lower_bound`` -- (default: 0) the algorithm returns a solution with cost
      larger or equal to that bound.

    OUTPUT:

    A pair ``(cost, ordering)`` representing the optimal ordering of the
    vertices and its cost.
    
    .. NOTE::

        Because of its current implementation, this algorithm only works on
        graphs on strictly less than 32 vertices. This can be changed to 63 if
        necessary, but 32 vertices already require 4GB of memory.

    TESTS:

    Giving anything else than a Graph::

        sage: from sage.graphs.graph_decompositions import cutwidth
        sage: cutwidth.cutwidth_dyn([])
        Traceback (most recent call last):
        ...
        ValueError: The parameter must be a Graph.

    Giving a too large Graph::

        sage: from sage.graphs.graph_decompositions import cutwidth
        sage: cutwidth.cutwidth_dyn(graphs.PathGraph(40))
        Traceback (most recent call last):
        ...
        ValueError: The graph should have at most 31 vertices !

    Giving a wrong type lower bound::

        sage: from sage.graphs.graph_decompositions import cutwidth
        sage: cutwidth.cutwidth_dyn(Graph(), lower_bound='toto')
        Traceback (most recent call last):
        ...
        ValueError: The specified lower bound must be an integer.
    """
    from sage.graphs.graph import Graph
    if not isinstance(G, Graph):
        raise ValueError("The parameter must be a Graph.")

    if G.order() >= 32:
        raise ValueError("The graph should have at most 31 vertices !")

    if not lower_bound in ZZ:
        raise ValueError("The specified lower bound must be an integer.")

    cdef FastDigraph g = FastDigraph(G)

    cdef unsigned int mem = 1 << g.n
    cdef uint8_t * neighborhoods = <uint8_t *> check_allocarray(mem, sizeof(uint8_t))
    
    memset(neighborhoods, <uint8_t> -1, mem)

    cdef int i,j , k

    sig_on()
    for k in range(lower_bound, G.size()+1):
        if exists(g, neighborhoods, 0, k) <= k:
            break
    sig_off()

    cdef list order = find_order(g, neighborhoods, k)

    sage_free(neighborhoods)

    return k, list( g.int_to_vertices[i] for i in order )


cdef inline int compute_edge_cut(FastDigraph g, int S):
    r"""
    Returns the number of edges from `S` to `V\S`.

    INPUT:

    - ``g`` -- a FastDigraph
    
    - ``S`` -- (integer) an integer describing a set of vertices
    """
    cdef int Sbar = ~S
    cdef int i
    cdef int cpt = 0
    for i in range(g.n):
        cpt += popcount32( (Sbar&g.graph[i]) * ((S>>i)&1) )

    return cpt


cdef inline int exists(FastDigraph g, uint8_t * neighborhoods, int current, int cost):
    """
    Check whether an ordering with the given cost exists, and updates data in
    the neighborhoods array at the same time. See the module's documentation.

    This method differs from method
    `sage.graphs.graph_decompositions.vertex_separation.exists` by a single
    line. We can certainly combine codes.
    """
    # If this is true, it means the set has not been evaluated yet
    if neighborhoods[current] == <uint8_t>-1:
        neighborhoods[current] = compute_edge_cut(g, current)

    # If the cost of this set is too high, there is no point in going further.
    # Same thing if the current set is the whole vertex set.
    if neighborhoods[current] > cost or (current == (1<<g.n)-1):
        return neighborhoods[current]

    # Minimum of the costs of the outneighbors
    cdef int mini = (1<<g.n)-1

    cdef int i
    cdef int next_set

    for i in range(g.n):
        if (current >> i)&1:
            continue

        # For each of the out-neighbors next_set of current
        next_set = current | 1<<i

        # Check whether there exists a cheap path toward {1..n}, and updated the
        # cost.
        mini = min(mini, exists(g, neighborhoods, next_set, cost))

        # We have found a path !
        if mini <= cost:
            return mini

    # Updating the cost of the current set with the minimum of the cost of its
    # outneighbors.
    neighborhoods[current] = mini

    return neighborhoods[current]


cdef list find_order(FastDigraph g, uint8_t * neighborhoods, int cost):
    """
    Returns the ordering once we are sure it exists

    This is a copy/paste of method
    `sage.graphs.graph_decompositions.vertex_separation.find_order` and so we
    can certainly directly import it.
    """
    cdef list ordering = []
    cdef int current = 0
    cdef int n = g.n
    cdef int i

    while n:
        # We look for n vertices

        for i in range(g.n):
            if (current >> i)&1:
                continue

            # Find the next set with small cost (we know it exists)
            next_set = current | 1<<i
            if neighborhoods[next_set] <= cost:
                ordering.append(i)
                current = next_set
                break

        # One less to find
        n -= 1

    return ordering
