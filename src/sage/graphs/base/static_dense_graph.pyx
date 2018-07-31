# cython: binding=True
r"""
Static dense graphs

This module gathers everything which is related to static dense graphs, i.e. :

- The vertices are integer from `0` to `n-1`
- No labels on vertices/edges
- No multiple edges
- No addition/removal of vertices

This being said, it is technically possible to add/remove edges. The data
structure does not mind at all.

It is all based on the binary matrix data structure described in
``misc/binary_matrix.pxi``, which is almost a copy of the bitset data
structure. The only difference is that it differentiates the rows (the vertices)
instead of storing the whole data in a long bitset, and we can use that.

For an overview of graph data structures in sage, see
:mod:`~sage.graphs.base.overview`.

Index
-----

**Cython functions**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    ``dense_graph_init`` | Fills a binary matrix with the information of a (di)graph.

**Python functions**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`is_strongly_regular` | Tests if a graph is strongly regular
    :meth:`triangles_count` | Return the number of triangles containing `v`, for every `v`
    :meth:`connected_subgraph_iterator` | Iterator over the connected subgraphs of order at most `k`

Functions
---------
"""
include "sage/data_structures/binary_matrix.pxi"
from cysignals.signals cimport sig_on, sig_off


cdef dict dense_graph_init(binary_matrix_t m, g, translation=False, to_undirected=False):
    r"""
    Initializes the binary matrix from a Sage (di)graph.

    INPUT:

    - ``binary_matrix_t m`` -- the binary matrix to be filled

    - ``g`` -- a graph or digraph

    - ``translation`` (boolean) -- whether to return a dictionary associating to
      each vertex its corresponding integer in the binary matrix.

    - ``to_undirected`` (boolean) -- whether to consider the graph as undirected or not.
    """
    cdef dict d_translation
    from sage.graphs.graph import Graph
    cdef bint is_undirected = isinstance(g, Graph) or to_undirected
    cdef int n = g.order()

    binary_matrix_init(m, n, n)

    # If the vertices are 0...n-1, let's avoid an unnecessary dictionary
    if g.vertices() == list(xrange(n)):
        if translation:
            d_translation = {i: i for i in range(n)}

        for i,j in g.edge_iterator(labels = False):
            binary_matrix_set1(m, i, j)
            if is_undirected:
                binary_matrix_set1(m, j, i)
    else:
        d_translation = {v:i for i,v in enumerate(g.vertices())}

        for u,v in g.edge_iterator(labels = False):
            binary_matrix_set1(m, d_translation[u], d_translation[v])
            if is_undirected:
                binary_matrix_set1(m, d_translation[v], d_translation[u])

    if translation:
        return d_translation

def is_strongly_regular(g, parameters = False):
    r"""
    Tests whether ``self`` is strongly regular.

    A simple graph `G` is said to be strongly regular with parameters `(n, k, \lambda,
    \mu)` if and only if:

        * `G` has `n` vertices.

        * `G` is `k`-regular.

        * Any two adjacent vertices of `G` have `\lambda` common neighbors.

        * Any two non-adjacent vertices of `G` have `\mu` common neighbors.

    By convention, the complete graphs, the graphs with no edges
    and the empty graph are not strongly regular.

    See :wikipedia:`Strongly regular graph`

    INPUT:

    - ``parameters`` (boolean) -- whether to return the quadruple `(n,
      k,\lambda,\mu)`. If ``parameters = False`` (default), this method only
      returns ``True`` and ``False`` answers. If ``parameters=True``, the
      ``True`` answers are replaced by quadruples `(n, k,\lambda,\mu)`. See
      definition above.

    EXAMPLES:

    Petersen's graph is strongly regular::

        sage: g = graphs.PetersenGraph()
        sage: g.is_strongly_regular()
        True
        sage: g.is_strongly_regular(parameters = True)
        (10, 3, 0, 1)

    And Clebsch's graph is too::

        sage: g = graphs.ClebschGraph()
        sage: g.is_strongly_regular()
        True
        sage: g.is_strongly_regular(parameters = True)
        (16, 5, 0, 2)

    But Chvatal's graph is not::

        sage: g = graphs.ChvatalGraph()
        sage: g.is_strongly_regular()
        False

    Complete graphs are not strongly regular. (:trac:`14297`) ::

        sage: g = graphs.CompleteGraph(5)
        sage: g.is_strongly_regular()
        False

    Completements of complete graphs are not strongly regular::

        sage: g = graphs.CompleteGraph(5).complement()
        sage: g.is_strongly_regular()
        False

    The empty graph is not strongly regular::

        sage: g = graphs.EmptyGraph()
        sage: g.is_strongly_regular()
        False

    If the input graph has loops or multiedges an exception is raised::

        sage: Graph([(1,1),(2,2)]).is_strongly_regular()
        Traceback (most recent call last):
        ...
        ValueError: This method is not known to work on graphs with
        loops. Perhaps this method can be updated to handle them, but in the
        meantime if you want to use it please disallow loops using
        allow_loops().
        sage: Graph([(1,2),(1,2)]).is_strongly_regular()
        Traceback (most recent call last):
        ...
        ValueError: This method is not known to work on graphs with
        multiedges. Perhaps this method can be updated to handle them, but in
        the meantime if you want to use it please disallow multiedges using
        allow_multiple_edges().
    """
    g._scream_if_not_simple()
    cdef binary_matrix_t m
    cdef bitset_t b_tmp
    cdef int n = g.order()
    cdef int inter
    cdef int i,j,l, k

    if g.size() == 0: # no vertices or no edges
        return False

    if g.is_clique():
        return False

    cdef list degree = g.degree()
    k = degree[0]
    if not all(d == k for d in degree):
        return False

    bitset_init(b_tmp, n)

    # m is now our copy of the graph
    dense_graph_init(m, g)

    cdef int llambda = -1
    cdef int mu = -1

    for i in range(n):
        for j in range(i+1,n):

            # The intersection of the common neighbors of i and j is a AND of
            # their respective rows. A popcount then returns its cardinality.
            bitset_and(b_tmp, m.rows[i], m.rows[j])
            inter = bitset_len(b_tmp)

            # Check that this cardinality is correct according to the values of lambda and mu
            if binary_matrix_get(m,i,j):
                if llambda == -1:
                    llambda = inter
                elif llambda != inter:
                    binary_matrix_free(m)
                    bitset_free(b_tmp)
                    return False
            else:
                if mu == -1:
                    mu = inter
                elif mu != inter:
                    binary_matrix_free(m)
                    bitset_free(b_tmp)
                    return False

    binary_matrix_free(m)
    bitset_free(b_tmp)

    if parameters:
        return (n,k,llambda,mu)
    else:
        return True

def triangles_count(G):
    r"""
    Return the number of triangles containing `v`, for every `v`.

    INPUT:

    - ``G``-- a simple graph

    EXAMPLES::

        sage: from sage.graphs.base.static_dense_graph import triangles_count
        sage: triangles_count(graphs.PetersenGraph())
        {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0}
        sage: sum(triangles_count(graphs.CompleteGraph(15)).values()) == 3*binomial(15,3)
        True
    """
    from sage.rings.integer import Integer
    G._scream_if_not_simple()
    cdef int n = G.order()

    cdef uint64_t * count = <uint64_t *> check_calloc(n, sizeof(uint64_t))

    cdef binary_matrix_t g
    dense_graph_init(g, G)

    cdef bitset_t b_tmp
    bitset_init(b_tmp, n)

    cdef int i,j
    cdef uint64_t tmp_count = 0

    for i in range(n):
        for j in range(i+1,n):
            if not bitset_in(g.rows[i],j):
                continue
            bitset_and(b_tmp, g.rows[i], g.rows[j])
            tmp_count = bitset_len(b_tmp)
            count[i] += tmp_count
            count[j] += tmp_count

    ans = {v:Integer(count[i]/2)
           for i,v in enumerate(G.vertices())}

    bitset_free(b_tmp)
    binary_matrix_free(g)
    sig_free(count)

    return ans

def connected_subgraph_iterator(G, k=None, vertices_only=False):
    r"""
    Iterator over the connected subgraphs of order at most `k`.

    This method implements a iterator over the connected subgraphs of the input
    (di)graph. Edge orientation is ignored.

    INPUT:

    - ``G`` -- a Graph or a DiGraph. Loops and multiple edges are allowed.

    - ``k`` (integer) -- maximum order of the connected subgraphs to report. By
      default, the method iterates over all connected subgraphs, which is
      equivalent to set `k == n`.

    - ``vertices_only`` (boolean) -- whether to return (Di)Graph (default) or
      list of vertices (``vertices_only==True``).

    EXAMPLES:

    The Path Graph of order `n` has `n * (n + 1) / 2` connected subgraphs::

        sage: G = graphs.PathGraph(10)
        sage: len(list(G.connected_subgraph_iterator(vertices_only=True)))
        55
        sage: len(list(G.connected_subgraph_iterator(vertices_only=False)))
        55

    The Complete Graph of order `n` has `2^n - 1` connected subgraphs::

        sage: G = graphs.CompleteGraph(5)
        sage: len(list(G.connected_subgraph_iterator(vertices_only=True))) == 2**G.order() - 1
        True
        sage: G = graphs.CompleteGraph(6)
        sage: len(list(G.connected_subgraph_iterator(vertices_only=True))) == 2**G.order() - 1
        True

    """
    k = G.order() if k is None else k
    if not G.order() or k < 1:
        return

    sig_on()
    cdef int n = G.order()
    cdef list int_to_vertex = G.vertices()
    cdef binary_matrix_t DG
    dense_graph_init(DG, G, translation=False, to_undirected=True)

    # We use a stack of bitsets. We need 3 bitsets per level with at most n + 1
    # levels, so 3 * n + 3 bitsets. We also need 1 bitset that we create at the
    # same time, so we need overall 3 * n + 4 bitsets
    cdef binary_matrix_t stack
    binary_matrix_init(stack, 3 * n + 4, n)

    cdef bitset_t current  # current subset of vertices
    cdef bitset_t left     # remaining vertices to consider
    cdef bitset_t boundary # neighbors of the current subset
    # candidate vertices for extending the current subset, i.e., vertices that
    # are both in left and in boundary
    cdef bitset_t candidates = stack.rows[3 * n + 3]

    cdef int l = 0
    cdef int u, v

    # We first generate subsets containing vertex 0, then the subsets containing
    # vertex 1 but not vertex 0 since we have already generated all subsets
    # containing 0, then subsets containing 2 but not 0 or 1, etc.
    for u in range(n):

        if vertices_only:
            yield [int_to_vertex[u]]
        else:
            yield G.subgraph([int_to_vertex[u]])

        # We initialize the loop with vertices u in current, {u+1, ..., n-1} in
        # left, and N(u) in boundary
        bitset_clear(stack.rows[0])
        bitset_add(stack.rows[0], u)
        bitset_set_first_n(stack.rows[1], u+1)
        bitset_complement(stack.rows[1], stack.rows[1])
        bitset_copy(stack.rows[2], DG.rows[u])
        l = 0

        while l >= 0:

            # We take the values at the top of the stack
            current = stack.rows[l]
            left = stack.rows[l + 1]
            boundary = stack.rows[l + 2]

            bitset_and(candidates, left, boundary)

            # Search for a candidate vertex v
            v = bitset_next(candidates, u+1)
            if v >= 0 and bitset_len(current) < k:
                # We select vertex v
                bitset_discard(left, v)

                # Since we have not modified l, the bitsets for iterating without v
                # are already at the top of the stack, with correct values

                # We also build the subset with v and so we add values at the
                # top of the stack
                l += 3
                bitset_copy(stack.rows[l], current)
                bitset_add(stack.rows[l], v)
                bitset_copy(stack.rows[l + 1], left)
                bitset_union(stack.rows[l + 2], boundary, DG.rows[v])

                # We yield that new subset
                if vertices_only:
                    yield [int_to_vertex[v] for v in range(u, n) if bitset_in(stack.rows[l], v)]
                else:
                    yield G.subgraph([int_to_vertex[v] for v in range(u, n) if bitset_in(stack.rows[l], v)])

            else:
                # We cannot extend the current subset, either due to a lack of
                # candidate (v == -1), or because the current subset has maximum
                # size. We pop
                l -= 3

    binary_matrix_free(stack)
    binary_matrix_free(DG)
    sig_off()
