#*****************************************************************************
#       Copyright (C) 2015 Michele Borassi michele.borassi@imtlucca.it
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
r"""
Interface to run Boost algorithms

Wrapper for a Boost graph. The Boost graphs are Cython C++ variables, and they
cannot be converted to Python objects: as a consequence, only functions defined
with cdef are able to create, read, modify, and delete these graphs.

A very important feature of Boost graph library is that all object are generic:
for instance, adjacency lists can be stored using different data structures,
and (most of) the functions work with all implementations provided. This feature
is implemented in our interface using fused types: however, Cython's support for
fused types is still experimental, and some features are missing. For instance,
there cannot be nested generic function calls, and no variable can have a
generic type, apart from the arguments of a generic function.

All the input functions use pointers, because otherwise we might have problems
with ``delete()``.

**Basic Boost Graph operations:**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :func:`clustering_coeff`  | Returns the clustering coefficient of all vertices in the graph.
    :func:`edge_connectivity` | Returns the edge connectivity of the graph.

Functions
---------
"""

from sage.graphs.graph import Graph
from sage.graphs.digraph import DiGraph
from sage.graphs.generic_graph import GenericGraph

include "sage/ext/interrupt.pxi"

cdef boost_graph_from_sage_graph(BoostGenGraph *g, g_sage):
    r"""
    Initializes the Boost graph ``g`` to be equal to ``g_sage``.

    The Boost graph ``*g`` must represent an empty graph (an exception is raised
    otherwise).
    """

    if not isinstance(g_sage, GenericGraph):
        raise ValueError("The input parameter must be a Sage graph.")

    if g.num_verts() > 0:
        raise ValueError("The Boost graph in input must be empty")

    N = g_sage.num_verts()
    cdef dict vertex_to_int = {v:i for i,v in enumerate(g_sage.vertices())}

    for i in range(N):
        g.add_vertex()

    for u,v in g_sage.edge_iterator(labels=None):
        g.add_edge(vertex_to_int[u], vertex_to_int[v])

cdef boost_edge_connectivity(BoostVecGenGraph *g):
    r"""
    Computes the edge connectivity of the input Boost graph.

    The output is a pair ``[ec,edges]``, where ``ec`` is the edge connectivity,
    ``edges`` is the list of edges in a minimum cut.
    """
    result = g[0].edge_connectivity()

    cdef int i
    edges = [(result.edges[i], result.edges[i+1])
             for i in range(0, result.edges.size(), 2)]

    return [result.ec, edges]

cpdef edge_connectivity(g):
    r"""
    Computes the edge connectivity of the input graph, using Boost.

    The output is a pair ``[ec,edges]``, where ``ec`` is the edge connectivity,
    ``edges`` is the list of edges in a minimum cut.

    .. SEEALSO::

        :meth:`sage.graphs.generic_graph.GenericGraph.edge_connectivity`

    EXAMPLES:

    Computing the edge connectivity of a clique::

        sage: from sage.graphs.base.boost_graph import edge_connectivity
        sage: g = graphs.CompleteGraph(5)
        sage: edge_connectivity(g)
        [4, [(0, 1), (0, 2), (0, 3), (0, 4)]]

    """
    # These variables are automatically deleted when the function terminates.
    cdef BoostVecGraph g_boost_und
    cdef BoostVecDiGraph g_boost_dir

    if isinstance(g, Graph):
        boost_graph_from_sage_graph(&g_boost_und, g)
        sig_check()
        return boost_edge_connectivity(&g_boost_und)

    elif isinstance(g, DiGraph):
        from sage.misc.stopgap import stopgap
        stopgap("The edge connectivity of directed graphs is not implemented " +
                "in Boost. The result may be mathematically unreliable.",18753)

        boost_graph_from_sage_graph(&g_boost_dir, g)
        sig_check()
        return boost_edge_connectivity(&g_boost_dir)

    else:
        raise ValueError("The input must be a Sage graph.")

cdef boost_clustering_coeff(BoostGenGraph *g, vertices):
    r"""
    Computes the clustering coefficient of all vertices in the list provided.

    The output is a pair ``[average_clustering_coefficient, clust_of_v]``, where
    ``average_clustering_coefficient`` is the average clustering of the vertices
    in variable ``vertices``, ``clust_of_v`` is a dictionary that associates to
    each vertex (stored as an integer) its clustering coefficient.
    """
    cdef result_cc result
    cdef dict clust_of_v

    if len(vertices) == g.num_verts():
        result = g[0].clustering_coeff_all()
        clust_of_v = {v:result.clust_of_v[v] for v in range(g.num_verts())}
        return [result.average_clustering_coefficient, clust_of_v]

    else:
        clust_of_v = {v:g[0].clustering_coeff(v) for v in vertices}
        return [(sum(clust_of_v.itervalues())/len(clust_of_v)), clust_of_v]


cpdef clustering_coeff(g, vertices = None):
    r"""
    Computes the clustering coefficient of the input graph, using Boost.

    The output is a pair ``[average_clustering_coefficient, clust_of_v]``, where
    ``average_clustering_coefficient`` is the average clustering of the vertices
    in variable ``vertices``, ``clust_of_v`` is a dictionary that associates to
    each vertex its clustering coefficient. If ``vertices`` is ``None``, all
    vertices are considered.

    .. SEEALSO::

        :meth:`sage.graphs.generic_graph.GenericGraph.clustering_coeff`

    INPUT:

    - ``g`` (Graph) - the input graph.

    - ``vertices`` (list) - the list of vertices we need to analyze (if
      ``None``, we will compute the clustering coefficient of all vertices).

    EXAMPLES:

    Computing the clustering coefficient of a clique::

        sage: from sage.graphs.base.boost_graph import clustering_coeff
        sage: g = graphs.CompleteGraph(5)
        sage: clustering_coeff(g)
        [1.0, {0: 1.0, 1: 1.0, 2: 1.0, 3: 1.0, 4: 1.0}]
        sage: clustering_coeff(g, vertices = [0,1,2])
        [1.0, {0: 1.0, 1: 1.0, 2: 1.0}]

    Of a non-clique graph with triangles::

        sage: g = graphs.IcosahedralGraph()
        sage: clustering_coeff(g, vertices=[1,2,3])
        [0.5, {1: 0.5, 2: 0.5, 3: 0.5}]

    With labels::

        sage: g.relabel(list("abcdefghiklm"))
        sage: clustering_coeff(g, vertices="abde")
        [0.5, {'a': 0.5, 'b': 0.5, 'd': 0.5, 'e': 0.5}]
    """
    sig_on()
    # These variables are automatically deleted when the function terminates.
    cdef BoostVecGraph g_boost
    cdef list g_vertices = g.vertices()
    cdef dict vertex_to_int = {v:i for i,v in enumerate(g_vertices)}

    if not isinstance(g, Graph):
        sig_off()
        raise ValueError("The input must be a Sage graph.")

    boost_graph_from_sage_graph(&g_boost, g)

    if vertices is None:
        vertices = g_vertices

    vertices_boost = [vertex_to_int[v] for v in vertices]
    [average_clustering, clust_v_int] = boost_clustering_coeff(&g_boost, vertices_boost)
    clust_v_sage = {g_vertices[v]: clust_v_int[v] for v in vertices_boost}
    sig_off()
    return [average_clustering, clust_v_sage]
