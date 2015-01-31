r"""
Asteroidal triples

**This module contains the following methods**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`is_asteroidal_triple_free` | Test if the input graph is asteroidal triple-free
    :meth:`is_AT_free` | Test if the input graph is asteroidal triple-free
    :meth:`has_an_asteroidal_triple` | Test if the input graph contains an asteroidal triple

    
Definition
----------

Three independent vertices of a graph form an *asteroidal triple* if every two
of them are connected by a path avoiding the neighborhood of the third one. A
graph is *asteroidal triple-free* (*AT-free*, for short) if it contains no
asteroidal triple [LB62]_.

Use :meth:`graph_classes.get_class('gc_61').description()` to get some known
properties of AT-free graphs.


Algorithm
---------

This module implements the  *Straightforward algorithm* recalled in [Koh04]_ and
due to [LB62]_ for testing if a graph is AT-free or not. This algorithm has time
complexity in `O(n^3)` and space complexity in `O(n^2)`.

This algorithm uses the *connected structure* of the graph, stored into a
`n\times n` matrix `M`. This matrix is such that `M[u][v]==0` if `v\in
(\{u\}\cup N(u))`, and otherwise `M[u][v]` is the unique identifier (a strictly
positive integer) of the connected component of `G\setminus(\{u\}\cup N(u))` to
which `v` belongs. This connected structure can be computed in time `O(n(n+m))`
using `n` BFS.

Now, a triple `u, v, w\in V` is an asteroidal triple if and only if it satisfies
`M[u][v]==M[u][w]` and `M[v][u]==M[v][w]` and `M[w][u]==M[w][v]`, assuming all
these values are positive. Indeed, if `M[u][v]==M[u][w]`, `v` and `w` are in the
same connected component of `G\setminus(\{u\}\cup N(u))`, and so there is a path
between `v` and `w` avoiding the neighborhood of `u`. The algorithm iterates
over all triples.


References
----------

.. [Koh04] E. Kohler. *Recognizing graphs without asteroidal triples*. Journal of
      Discrete Algorithms 2(4):439-452, Dec. 2004
      http://dx.doi.org/10.1016/j.jda.2004.04.005

.. [LB62] C. G. Lekkerkerker, J. Ch. Boland. *Representation of a finite graph
      by a set of intervals on the real line*. Fundamenta Mathematicae,
      51:45-64, 1962.


Methods
-------
"""
#*****************************************************************************
# Copyright (C) 2015 David Coudert <david.coudert@inria.fr>
#
# Distributed under the terms of the GNU General Public License (GPL)
# http://www.gnu.org/licenses/
#*****************************************************************************

include 'sage/ext/interrupt.pxi'
include "sage/data_structures/bitset.pxi"

from libc.stdint cimport uint32_t
from sage.graphs.base.static_sparse_graph cimport short_digraph, init_short_digraph, free_short_digraph


def is_asteroidal_triple_free(G):
    """
    Test if the input graph is asteroidal triple-free

    An independent set of three vertices such that each pair is joined by a path
    that avoids the neighborhood of the third one is called an *asteroidal
    triple*. A graph is asteroidal triple-free (AT-free) if it contains no
    asteroidal triples. See the :mod:`module's documentation
    <sage.graphs.asteroidal_triples>` for more details.

    This method returns ``True`` is the graph is AT-free and ``False`` otherwise.

    INPUT:

    - ``G`` -- a Graph

    EXAMPLES:

    The complete graph is AT-free, and its line graph too::

        sage: from sage.graphs.asteroidal_triples import *
        sage: G = graphs.CompleteGraph(5)
        sage: is_AT_free(G)
        True
        sage: LG = G.line_graph()
        sage: is_AT_free(LG)
        True
        sage: LLG = LG.line_graph()
        sage: is_AT_free(LLG)
        False

    The PetersenGraph is not AT-free::

        sage: from sage.graphs.asteroidal_triples import *
        sage: G = graphs.PetersenGraph()
        sage: is_AT_free(G)
        False

    """
    return not has_an_asteroidal_triple(G, certificate=False)

is_AT_free = is_asteroidal_triple_free


def has_an_asteroidal_triple(G, certificate=False):
    """
    Test if the input graph contains an asteroidal triple

    An independent set of three vertices such that each pair is joined by a path
    that avoids the neighborhood of the third one is called an *asteroidal
    triple*.

    This method implements the *Straightforward algorithm* recalled in [Koh04]_
    and due to [LB62]_.  See the :mod:`module's documentation
    <sage.graphs.asteroidal_triples>` for more details.

    INPUTS:

    - ``G`` -- a Graph

    - ``certificate`` -- (default: False) By default, this method returns
      ``True`` if the graph contains an asteroidal triple and ``False``
      otherwise. When ``certificate==True``, this method returns a list of three
      vertices forming an asteroidal triple if such a triple is found, and the
      empty list otherwise.

    EXAMPLES:

    The Petersen graph has an asteroidal triple::

        sage: from sage.graphs.asteroidal_triples import has_an_asteroidal_triple
        sage: G = graphs.PetersenGraph()
        sage: has_an_asteroidal_triple(G)
        True
        sage: has_an_asteroidal_triple(G, certificate=True)
        [0, 2, 6]

    The Path graph has no asteroidal triple::

        sage: from sage.graphs.asteroidal_triples import has_an_asteroidal_triple
        sage: G = graphs.PathGraph(5)
        sage: has_an_asteroidal_triple(G)
        False
        sage: has_an_asteroidal_triple(G, certificate=True)
        []

    TEST:

    Giving anything else than a Graph::

        sage: from sage.graphs.asteroidal_triples import has_an_asteroidal_triple
        sage: has_an_asteroidal_triple(DiGraph())
        Traceback (most recent call last):
        ...
        ValueError: The first parameter must be a Graph.

    """
    from sage.graphs.graph import Graph
    if not isinstance(G, Graph):
        raise ValueError("The first parameter must be a Graph.")

    cdef int n = G.order()
    cdef int i

    # ==> Trivial cases
    if n<3:
        return False if not certificate else []



    # ==> Initialize some data structures for has_an_asteroidal_triple_C

    # Copying the whole graph to obtain the list of neighbors quicker than by
    # calling out_neighbors.  This data structure is well documented in the
    # module sage.graphs.base.static_sparse_graph
    cdef short_digraph sd
    init_short_digraph(sd, G)
    cdef uint32_t ** p_vertices = sd.neighbors

    cdef bitset_t seen
    bitset_init(seen, n)

    cdef uint32_t * waiting_list         = <uint32_t *> sage_malloc(n * sizeof(uint32_t))
    cdef uint32_t * _connected_structure = <uint32_t *> sage_calloc(n * n, sizeof(uint32_t))
    cdef uint32_t ** connected_structure = <uint32_t **> sage_malloc(n * sizeof(uint32_t *))

    if waiting_list==NULL or _connected_structure==NULL or connected_structure==NULL:
        bitset_free(seen)
        sage_free(waiting_list)
        sage_free(_connected_structure)
        sage_free(connected_structure)
        free_short_digraph(sd)
        raise MemoryError()

    connected_structure[0] = _connected_structure
    for i in range(n-1):
        connected_structure[i+1] = connected_structure[i] + n

    cdef list ret = list()


    # ==> call has_an_asteroidal_triple_C

    try:
        sig_on()
        ret = has_an_asteroidal_triple_C(n, p_vertices, connected_structure, waiting_list, seen)
        sig_off()

    finally:
        # Release memory
        bitset_free(seen)
        sage_free(waiting_list)
        sage_free(_connected_structure)
        sage_free(connected_structure)
        free_short_digraph(sd)


    # ==> We return the result

    if certificate:
        if ret:
            V = G.vertices()
            ret = [V[i] for i in ret]
        return ret

    return True if ret else False


cdef list has_an_asteroidal_triple_C(int n,
                                     uint32_t ** p_vertices,
                                     uint32_t ** connected_structure,
                                     uint32_t *  waiting_list,
                                     bitset_t seen):
    """
    PARAMETERS:

    - ``p_vertices`` -- bidimensional array allowing to access to the list of
      neighbors of the graph quicker than by calling out_neighbors.  This data
      structure is well documented in the module
      sage.graphs.base.static_sparse_graph

    - ``connected_structure`` -- bidimensional array of size `n\times n` used to
      store the connected structure of the graph. All its cells must initially
      be set to 0.

    - ``waiting_list`` -- an array of size `n` to be used for BFS.

    - ``seen`` -- a bitset of size `n`.
    """
    cdef uint32_t waiting_beginning = 0
    cdef uint32_t waiting_end       = 0
    cdef uint32_t idx_cc            = 0
    cdef uint32_t source, u, v
    cdef uint32_t * p_tmp
    cdef uint32_t * end


    # ==> We build the connected structure

    # We run n different BFS taking each vertex as a source
    for source in range(n):

        # The source is forbidden and seen
        bitset_clear(seen)
        bitset_add(seen, source)

        # The neighbors of the source are forbidden and seen
        p_tmp = p_vertices[source]
        end = p_vertices[source+1]
        # Iterating over all the outneighbors u of v
        while p_tmp < end:
            bitset_add(seen, p_tmp[0])
            p_tmp += 1

        # We now search for an unseen vertex
        v = bitset_first_in_complement(seen)
        while v!=-1:
            # and added it to the queue
            waiting_list[0] = v
            waiting_beginning = 0
            waiting_end = 0

            # We start a new connected component
            idx_cc += 1
            bitset_add(seen, v)
            connected_structure[source][v] = idx_cc

            # For as long as there are vertices left to explore in this
            # component
            while waiting_beginning <= waiting_end:
            
                # We pick the first one
                v = waiting_list[waiting_beginning]
                p_tmp = p_vertices[v]
                end = p_vertices[v+1]

                # Iterating over all the outneighbors u of v
                while p_tmp < end:
                    u = p_tmp[0]

                    # If we notice one of these neighbors is not seen yet, we
                    # add it to the queue to be explored later
                    if not bitset_in(seen, u):
                        waiting_end += 1
                        waiting_list[waiting_end] = u
                        bitset_add(seen, u)
                        connected_structure[source][u] = idx_cc

                    p_tmp += 1

                waiting_beginning += 1

            # We search for a possibly unseen vertex
            v = bitset_first_in_complement(seen)


    # ==> Now that we have the component structure of the graph, we search for
    # an asteroidal triple.

    for source in range(n-2):
        for u from source+1 <= u < n-1:
            if connected_structure[source][u]>0:
                for v from u+1 <= v < n:
                    if connected_structure[source][u] == connected_structure[source][v] \
                        and connected_structure[u][source] == connected_structure[u][v] \
                        and connected_structure[v][source] == connected_structure[v][u]:
                        # We have found an asteroidal triple
                        return [source, u, v]

    # No asteroidal triple was found
    return []
