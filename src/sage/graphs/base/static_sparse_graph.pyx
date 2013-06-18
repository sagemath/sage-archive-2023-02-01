r"""
Static Sparse Graphs

What is the point ?
-------------------

This class implements a Cython (di)graph structure made for efficiency. The
graphs are *static*, i.e. no add/remove vertex/edges methods are available, nor
can they easily or efficiently be implemented within this data structure.

The data structure, however, is made to save the maximum amount of computations
for graph algorithms whose main operation is to *list the out-neighbours of a
vertex* (which is precisely what BFS, DFS, distance computations and the
flow-related stuff waste their life on).

The code contained in this module is written C-style. While Sage needs a class
for static graphs (not available today, i.e. 2012-01-13) it is not what we try
to address here. The purpose is efficiency and simplicity.

Author:

- Nathann Cohen (2011)

Data structure
--------------

.. image:: ../../../media/structure.png

The data structure is actually pretty simple and compact. ``short_digraph`` has
three fields

    * ``n`` (``int``) -- the number of vertices in the graph.

    * ``edges`` (``unsigned short *``) -- array whose length is the number of
      edges of the graph.

    * ``neighbors`` (``unsigned short **``) -- this array has size `n+1`, and
      describes how the data of ``edges`` should be read : the neighbors of
      vertex `i` are the elements of ``edges`` addressed by
      ``neighbors[i]...neighbors[i+1]-1``. The element ``neighbors[n]``, which
      corresponds to no vertex (they are numbered from `0` to `n-1`) is present
      so that it remains easy to enumerate the neighbors of vertex `n-1` : the
      last of them is the element addressed by ``neighbors[n]-1``.


In the example given above, vertex 0 has 2,3,5,7,8 and 9 as out-neighbors, but
not 4, which is an out-neighbour of vertex 1. Vertex `n-1` has 2, 5, 8 and 9 as
out-neighbors. `\text{neighbors[n]}` points toward the cell immediately *after*
the end of `\text{edges}`, hence *outside of the allocated memory*. It is used
to indicate the end of the outneighbors of vertex `n-1`

**Iterating over the edges**

This is *the one thing* to have in mind when working with this data structure::

    cdef list_edges(short_digraph g):
        cdef ushort * v
        cdef ushort * end
        cdef int i

        for i in range(g.n):
            v = g.neighbors[i]
            end = g.neighbors[i+1]

            while v < end:
                print "There is an edge from "+str(i)+" to "+str(v[0])
                v += 1

**Advantages**

Three great points :

    * The neighbors of a vertex are contiguous in memory.
    * Going from one to the other amounts to increasing a pointer.
    * Storing such graphs is incredibly cheaper than storing Python structures.

Well, I think it would be hard to have anything more efficient than that to
enumerate out-neighbors in sparse graphs ! :-)

Technical details
-----------------

    * When creating a ``fast_digraph`` from a ``Graph`` or ``DiGraph`` named
      ``G``, the `i^{\text{th}}` vertex corresponds to ``G.vertices()[i]``

    * The data structure does not support edge labels

    * In its current implementation (with ``unsigned short`` variables), the
      data structure can handle graphs with at most 65535 vertices. If
      necessary, changing it to ``int`` is totally straightforward.

    * Some methods return ``bitset_t`` objets when lists could be
      expected. There is a very useful ``bitset_list`` function for this kind of
      problems :-)

    * The codes of this module are well documented, and many answers can be
      found directly in the code.

Todo list

    * Adjacency test. The data structure can support it efficiently through
      dichotomy. It would require to sort the list of edges as it is not done at
      the moment. Some calls to the C function ``qsort`` would be sufficient.

Cython functions
----------------

``init_short_digraph(short_digraph g, G)``

    This method initializes ``short_digraph g`` so that it represents ``G``. If
    ``G`` is a ``Graph`` objet (and not a ``DiGraph``), an edge between two
    vertices `u` and `v` is replaced by two arcs in both directions.

``int n_edges(short_digraph g)``

    Returns the number of edges in ``g``

``int out_degree(short_digraph g, int i)``

    Returns the out-degree of vertex `i` in ``g``

``init_empty_copy(short_digraph dst, short_digraph src)``

    Allocates memory for ``dst`` so that it can contain as many vertices and
    edges as ``src``. Its content is purely random, though.

``init_reverse(short_digraph dst, short_digraph src)``

    Initializes ``dst`` so that it represents a copy of ``src`` in which all
    edges go in the opposite direction.

``can_be_reached_from(short_digraph g, int src, bitset_t reached)``

    Assuming ``bitset_t reached`` has size at least ``g.n``, this method updates
    ``reached`` so that it represents the set of vertices that can be reached
    from ``src`` in ``g``.

``strongly_connected_component_containing_vertex(short_digraph g, short_digraph g_reversed, int v, bitset_t scc)``

    Assuming ``bitset_t reached`` has size at least ``g.n``, this method updates
    ``scc`` so that it represents the vertices of the strongly connected
    component containing ``v`` in ``g``. The variable ``g_reversed`` is assumed
    to represent the reverse of ``g``.

``free_short_digraph(short_digraph g)``

    Free the ressources used by ``g``

What is this module used for ?
------------------------------

At the moment, it is only used in the :mod:`sage.graphs.distances_all_pairs` module.

Python functions
----------------

These functions are available so that Python modules from Sage can call the
Cython routines this module implements (as they can not directly call methods
with C arguments).
"""

include "sage/misc/bitset.pxi"

##############################################################################
#       Copyright (C) 2010 Nathann Cohen <nathann.cohen@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
##############################################################################

from sage.graphs.base.c_graph cimport CGraph

cdef int init_short_digraph(short_digraph g, G) except -1:

    # g.n is unsigned short, so -1 is actually the maximum value possible.
    g.n = -1

    if G.order() > g.n:
        raise ValueError("This structure can handle at most "+str(<int> g.n)+" vertices !")
    else:
        g.n = G.order()

    cdef int isdigraph

    from sage.graphs.all import Graph, DiGraph

    if isinstance(G, DiGraph):
        isdigraph = 1
    elif isinstance(G, Graph):
         isdigraph = 0
    else:
        raise ValueError("The source graph must be either a DiGraph or a Graph object !")

    cdef list vertices = G.vertices()
    cdef dict v_to_id = {}
    cdef int i,j

    cdef int n_edges = G.size() if isdigraph else 2*G.size()

    for i, v in enumerate(vertices):
        v_to_id[v] = i

    g.edges = <ushort *> sage_malloc(n_edges*sizeof(ushort))
    if g.edges == NULL:
        raise ValueError("Problem while allocating memory (edges)")

    g.neighbors = <ushort **> sage_malloc((1+<int>g.n)*sizeof(ushort *))
    if g.neighbors == NULL:
        raise ValueError("Problem while allocating memory (neighbors)")


    # Initializing the value of neighbors
    g.neighbors[0] = g.edges

    cdef CGraph cg = <CGraph> G._backend
    for i in range(1,(<int>g.n)+1):
        g.neighbors[i] = g.neighbors[i-1] + <int> (cg.out_degree(vertices[i-1]) if isdigraph else G.degree(vertices[i-1]))

    for u,v in G.edge_iterator(labels = False):
        i = v_to_id[u]
        j = v_to_id[v]

        g.neighbors[i][0] = j
        g.neighbors[i] += 1

        if not isdigraph:
            g.neighbors[j][0] = i
            g.neighbors[j] += 1

    # Reinitializing the value of neighbors
    for g.n> i >0:
        g.neighbors[i] = g.neighbors[i-1]

    g.neighbors[0] = g.edges

cdef inline int n_edges(short_digraph g):
    # The number of edges is nothing but a difference of pointers
    return <int> (g.neighbors[g.n]-g.edges)

cdef inline int out_degree(short_digraph g, int i):
    # The out-degree is nothing but a difference of pointers
    return <int> (g.neighbors[i+1]-g.neighbors[i])

cdef int init_empty_copy(short_digraph dst, short_digraph src) except -1:
    dst.n = src.n

    dst.edges = <ushort *> sage_malloc(n_edges(src)*sizeof(ushort))
    if dst.edges == NULL:
        raise ValueError("Problem while allocating memory (edges)")

    dst.neighbors = <ushort **> sage_malloc((src.n+1)*sizeof(ushort *))
    if dst.neighbors == NULL:
        raise ValueError("Problem while allocating memory (neighbors)")

cdef int init_reverse(short_digraph dst, short_digraph src) except -1:
    # Allocates memory for dst
    init_empty_copy(dst, src)

    # Avoiding a later segfault
    if dst.n == 0:
        return 0

    #### 1/3
    #
    # In a first pass, we count the in-degrees of each vertex and store it in a
    # vector. With this information, we can initialize dst.neighbors to its
    # correct value. The content of dst.edges is not touched at this level.

    cdef int * in_degree = <int *> sage_malloc(src.n*sizeof(int))
    if in_degree == NULL:
        raise ValueError("Problem while allocating memory (in_degree)")

    # Counting the degrees
    memset(in_degree, 0, src.n*sizeof(int))

    cdef ushort * v = src.edges
    cdef ushort * end = src.neighbors[src.n]
    while v < end:
        in_degree[v[0]] += 1
        v += 1

    # Updating dst.neighbors
    cdef int i
    dst.neighbors[0] = dst.edges
    for i in range(1, src.n+1):
        dst.neighbors[i] = dst.neighbors[i-1] + in_degree[i-1]
    sage_free(in_degree)

    #### 2/3
    #
    # Second pass : we list the edges again, and add them in dst.edges. Doing
    # so, we will change the value of dst.neighbors, but that is not so bad as
    # we can fix it afterwards.
    for i in range(0, src.n):
        v = src.neighbors[i]
        end = src.neighbors[i+1]

        while v < end:
            dst.neighbors[v[0]][0] = i
            dst.neighbors[v[0]] += 1
            v += 1

    #### 3/3
    #
    # Final step : set the correct values of dst.neighbors again. It is easy, as
    # the correct value of dst.neighbors[i] is actually dst.neighbors[i-1]
    for src.n> i >0:
        dst.neighbors[i] = dst.neighbors[i-1]
    dst.neighbors[0] = dst.edges

cdef int can_be_reached_from(short_digraph g, int src, bitset_t reached) except -1:
    if g.n == 0:
        return 0

    # Initializing the set of vertices reached by setting only bit src
    bitset_set_first_n(reached, 0)
    bitset_add(reached, src)

    # We will be doing a Depth-First Search. We allocate the stack we need for
    # that, and put "src" on top of it.
    cdef ushort * stack = <ushort *> sage_malloc(g.n*sizeof(ushort))
    if stack == NULL:
        raise ValueError("Problem while allocating memory (stack)")

    stack[0] = src
    cdef int stack_size = 1

    # What we need to iterate over the edges
    cdef int i
    cdef ushort * v
    cdef ushort * end

    # Plain old DFS ...
    #
    #If there is something left on the stack, we remove it consider each of its
    # neighbors. If we find any which has not been reached yet, we set its
    # corresponding bit in the reached bitset, and add it on top of the stack.

    while stack_size:
        stack_size -= 1
        i = stack[stack_size]

        v = g.neighbors[i]
        end = g.neighbors[i+1]

        while v < end:
            if not bitset_in(reached, v[0]):
                bitset_add(reached, v[0])
                stack[stack_size] = v[0]
                stack_size += 1

            v += 1

    sage_free(stack)

cdef strongly_connected_component_containing_vertex(short_digraph g, short_digraph g_reversed, int v, bitset_t scc):

    # Computing the set of vertices that can be reached from v in g
    can_be_reached_from(g, v, scc)

    # Computing the set of vertices that can be reached from v in g *reversed*
    cdef bitset_t scc_reversed
    bitset_init(scc_reversed, g.n)
    can_be_reached_from(g_reversed, v, scc_reversed)

    # The scc containing v is the intersection of both sets
    bitset_intersection(scc, scc, scc_reversed)

cdef void free_short_digraph(short_digraph g):

    if g.edges != NULL:
        sage_free(g.edges)

    if g.neighbors != NULL:
        sage_free(g.neighbors)

def strongly_connected_components(G):
    r"""
    Returns the strongly connected components of the given DiGraph.

    INPUT:

    - ``G`` -- a DiGraph.

    .. NOTE::

        This method has been written as an attempt to solve the slowness
        reported in :trac:`12235`. It is not the one used by
        :meth:`sage.graphs.digraph.DiGraph.strongly_connected_components` as
        saving some time on the computation of the strongly connected components
        is not worth copying the whole graph, but it is a nice way to test this
        module's functions. It is also tested in the doctest or
        :meth:`sage.graphs.digraph.DiGraph.strongly_connected_components`.

    EXAMPLE::

        sage: from sage.graphs.base.static_sparse_graph import strongly_connected_components
        sage: g = digraphs.ButterflyGraph(2)
        sage: strongly_connected_components(g)
        [[('00', 0)], [('00', 1)], [('00', 2)], [('01', 0)], [('01', 1)], [('01', 2)],
        [('10', 0)], [('10', 1)], [('10', 2)], [('11', 0)], [('11', 1)], [('11', 2)]]
    """

    if G.order() == 0:
        return [[]]

    # To compute the connected component containing a given vertex v, we take
    # the intersection of the set of vertices that can be reached from v in G
    # and the set of vertices that can be reached from v in G reversed.
    #
    # That's all that happens here.

    cdef list answer = []
    cdef list vertices = G.vertices()
    cdef short_digraph g, gr

    init_short_digraph(g, G)

    init_reverse(gr, g)

    cdef bitset_t seen
    bitset_init(seen, g.n)
    bitset_set_first_n(seen, 0)


    cdef bitset_t scc
    bitset_init(scc, g.n)
    bitset_set_first_n(scc, 0)

    cdef int v
    while bitset_len(seen) < g.n:
        v = bitset_first_in_complement(seen)
        strongly_connected_component_containing_vertex(g, gr, v, scc)
        answer.append([vertices[i] for i in bitset_list(scc)])
        bitset_union(seen, seen, scc)

    bitset_free(seen)
    bitset_free(scc)
    free_short_digraph(g)
    free_short_digraph(gr)
    return answer

