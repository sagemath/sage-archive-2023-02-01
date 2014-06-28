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

Functions
---------
"""

cdef dict dense_graph_init(binary_matrix_t m, g, translation = False):
    r"""
    Initializes the binary matrix from a Sage (di)graph.

    INPUT:

    - ``binary_matrix_t m`` -- the binary matrix to be filled

    - ``g`` -- a graph or digraph

    - ``translation`` (boolean) -- whether to return a dictionary associating to
      each vertex its corresponding integer in the binary matrix.
    """
    cdef dict d_translation
    from sage.graphs.graph import Graph
    cdef int is_undirected = isinstance(g, Graph)
    cdef int n = g.order()

    binary_matrix_init(m,n,n)
    binary_matrix_fill(m,0)

    # If the vertices are 0...n-1, let's avoid an unnecessary dictionary
    if g.vertices() == range(n):
        if translation:
            d_translation = {i:i for i in range(n)}

        for i,j in g.edge_iterator(labels = False):
            binary_matrix_set1(m,i,j)
            if is_undirected:
                binary_matrix_set1(m,j,i)
    else:
        d_translation = {v:i for i,v in enumerate(g.vertices())}

        for u,v in g.edge_iterator(labels = False):
            binary_matrix_set1(m,d_translation[u],d_translation[v])
            if is_undirected:
                binary_matrix_set1(m,d_translation[v],d_translation[u])

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

    # m i now our copy of the graph
    dense_graph_init(m, g)

    cdef int llambda = -1
    cdef int mu = -1

    for i in range(n):
        for j in range(i+1,n):

            # The intersection of the common neighbors of i and j is a AND of
            # their respective rows. A popcount then returns its cardinality.
            inter = 0
            for l in range(m.width):
                inter += __builtin_popcountl(m.rows[i][l] & m.rows[j][l])

            # Check that this cardinality is correct according to the values of lambda and mu
            if binary_matrix_get(m,i,j):
                if llambda == -1:
                    llambda = inter
                elif llambda != inter:
                    binary_matrix_free(m)
                    return False
            else:
                if mu == -1:
                    mu = inter
                elif mu != inter:
                    binary_matrix_free(m)
                    return False

    binary_matrix_free(m)

    if parameters:
        return (n,k,llambda,mu)
    else:
        return True
