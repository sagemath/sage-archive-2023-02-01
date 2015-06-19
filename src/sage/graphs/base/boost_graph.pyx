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

Wrapper for a Boost graphs. The Boost graphs must be created with cdef, and
they cannot be shared by different Python functions, because they cannot be
converted to Python objects.

It is possible to define generic functions, whose input is a generic Boost
graph, but there cannot be nested generic function calls, otherwise Cython
raises an error. Furthermore, no variable can have a generic type, apart from
the arguments of a generic function.

All the input functions use pointers, because otherwise we might have problems
with ``delete()``.

**Basic Boost Graph operations:**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~edge_connectivity` | Returns the edge connectivity of the graph.

"""


from sage.graphs.graph import Graph
from sage.graphs.digraph import DiGraph
from sage.graphs.generic_graph import GenericGraph

include "sage/ext/interrupt.pxi"

cdef boost_graph_from_sage_graph(BoostGenGraph *g, g_sage):
    r"""
    Initializes the Boost graph ``g`` to be equal to g_sage.

    The Boost output graph must be generated elsewhere, because Cython does not
    allow to instantiate a generic graph. However, ``g`` must be empty,
    otherwise an error is raised.
    """

    if not isinstance(g_sage, GenericGraph):
        raise ValueError("The input parameter must be a Sage graph.")

    if g.num_verts() > 0:
        raise ValueError("The Boost graph in input must be empty")

    N = g_sage.num_verts()
    vertex_to_int = {v:i for i,v in enumerate(g_sage.vertices())}

    for i in range(N):
        g.add_vertex()

    for u,v in g_sage.edge_iterator(labels=None):
        g.add_edge(vertex_to_int[u], vertex_to_int[v])

cdef edge_connectivity(BoostVecGenGraph *g):
    r"""
    Computes the edge connectivity of the input Boost graph.

    The adjacency lists of the input graph should be vectors, because the Boost
    edge connectivity algorithm requires this kind of input. The output is a
    pair ``[ec,edges]``, where ``ec`` is the edge connectivity, ``edges`` is the
    list of edges in a minimum cut. It is based on the function
    ``edge_connectivity`` in file ``boost_interface``, which outputs a vector,
    having in position 0 the edge connectivity, followed by a list of arcs.

    WARNING: this function is a backend, which should not be called by the
    standard user. We suggest to use the method
    :meth:`edge_connectivity <sage.graphs.generic_graph.GenericGraph>`
    instead.
    """
    result = g[0].edge_connectivity()

    cdef int ec = result.ec

    edges = []

    for i in range(0, result.edges.size(), 2):
        edges.append([result.edges[i], result.edges[i+1]])

    return[ec, edges]

cpdef boost_edge_connectivity(g):
    r"""
    Computes the edge connectivity of the input graph, using the Boost
    algorithm.

    The output is a pair ``[ec,edges]``, where ``ec`` is the edge connectivity,
    ``edges`` is the list of edges in a minimum cut.

    WARNING: this function is a backend, which should not be called by the
    standard user. We suggest to use the method
    :meth:`edge_connectivity <sage.graphs.generic_graph.GenericGraph.edge_connectivity>`
    instead.

    EXAMPLES:

    Computing the edge connectivity of a clique::

        sage: from sage.graphs.base.boost_graph import boost_edge_connectivity
        sage: g = graphs.CompleteGraph(5)
        sage: boost_edge_connectivity(g)
        [4, [[0, 1], [0, 2], [0, 3], [0, 4]]]

    """
    sig_on()
    # These variables are automatically deleted when the function terminates.
    cdef BoostVecGraph g_boost_und
    cdef BoostVecDiGraph g_boost_dir

    if isinstance(g, Graph):
        boost_graph_from_sage_graph(&g_boost_und, g)
        sig_off()
        return edge_connectivity(&g_boost_und)

    elif isinstance(g, DiGraph):
        boost_graph_from_sage_graph(&g_boost_dir, g)
        print ("The directed edge connectivity algorithm implemented in the "
               "Boost graph library is not reliable. The result could be "
               "wrong.")
        sig_off()
        return edge_connectivity(&g_boost_dir)

    else:
        sig_off()
        raise ValueError("The input must be a Sage graph.")
    