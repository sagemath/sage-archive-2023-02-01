r"""
Products of graphs

This module gathers everything related to graph products. At the moment it
contains different implementations of recognition algorithms for graphs that can
be written as a cartesian product of smaller ones.

References:

  .. [HIK11] Handbook of Product Graphs,
    R. Hammack, W. Imrich, S. Klavzar,
    CRC press, 2011

Author:

- Nathann Cohen (May 2012 -- coded while watching the election of Francois
  Hollande on TV)

Cartesian product of graphs -- the recognition problem
------------------------------------------------------

First, a definition:

  **Definition** The cartesian product of two graphs `G` and `H`, denoted
  `G\square H`, is a graph defined on the pairs `(g, h)\in V(G)\times V(H)`.

  Two elements `(g, h),(g', h')\in V(G\square H)` are adjacent in `G\square H`
  if and only if :

  - `g=g'` and `hh'\in H`; or
  - `h=h'` and `gg'\in H`

Two remarks follow :

#. The cartesian product is commutative

#. Any edge `uv` of a graph `G_1 \square \cdots \square G_k` can be given a color
   `i` corresponding to the unique index `i` such that `u_i` and `v_i` differ.

The problem that is of interest to us in the present module is the following:

  **Recognition problem** Given a graph `G`, can we guess whether there exist
  graphs `G_1, ..., G_k` such that `G=G_1\square \cdots \square G_k` ?

This problem can actually be solved, and the resulting factorization is
unique. What is explained below can be found in the book *Handbook of Product
Graphs* [HIK11]_.

Everything is actually based on a simple observation. Given a graph `G`, finding
out whether `G` can be written as the product of several graphs can be attempted
by trying to color its edges according to some rules. Indeed, if we are to color
the edges of `G` in such a way that each color class represents a factor of `G`,
we must ensure several things.

  **Remark 1** In any cycle of `G` no color can appear exactly once.

  Indeed, if only one edge `uv` of a cycle were labelled with color `i`, it
  would mean that:

  #. The only difference between `u` and `v` lies in their `i` th coordinate

  #. It is possible to go from `u` to `v` by changing only coordinates
     different from the `i` th

  A contradiction indeed.

  .. MATH::

      \begin{picture}(200,200)(-100,-100)
      \put(-50,0){\circle*{5}}
      \put(50,0){\circle*{5}}
      \put(-25,-35){\circle*{5}}
      \put(25,-35){\circle*{5}}
      \put(-25,35){\circle*{5}}
      \put(25,35){\circle*{5}}
      \qbezier(-50,0)(-37.5, 17.5)(-25,35)
      \qbezier(-50,0)(-37.5, -17.5)(-25,-35)
      \qbezier(50,0)(37.5, 17.5)(25,35)
      \qbezier(50,0)(37.5, -17.5)(25,-35)
      \qbezier(-25,35)(0,35)(25,35)
      \linethickness{1.5pt}
      \qbezier[10](-25,-35)(0,-35)(25,-35)
      \put(0,0){\makebox(0,0){No way !!!}}
      \end{picture}

  That means that, for instance, the edges of a triangle necessarily have the
  same color.

  **Remark 2** If two consecutive edges `u_1u_2` and `u_2u_3` have different
  colors, there necessarily exists a unique vertex `u_4` different from `u_2`
  and incident to both `u_1` and `u_3`.

  In this situation, opposed edges necessarily have the same colors because of
  the previous remark.

  .. MATH::

    \begin{picture}(200,200)(-100,-100)
    \put(-50,0){\circle*{5}}
    \put(50,0){\circle*{5}}
    \put(0,50){\circle*{5}}
    \put(0,-50){\circle*{5}}
    \qbezier(50,0)(25,25)(0,50)
    \qbezier(-50,0)(-25,-25)(0,-50)
    \linethickness{1.5pt}
    \qbezier[10](0,50)(-25,25)(-50,0)
    \qbezier[10](0,-50)(25,-25)(50,0)
    \put(0,55){\makebox(0,0)[b]{$u_3=(g',h')$}}
    \put(0,-55){\makebox(0,0)[t]{$u_1=(g,h)$}}
    \put(55,0){\makebox(0,0)[l]{$u_4=(g,h')$}}
    \put(-55,0){\makebox(0,0)[r]{$u_2=(g',h)$}}
    \end{picture}

  As a corollary, we also know that:

  #. If two vertices `u,v` have a *unique* common neighbor `x`, then `ux` and
     `uv` have the same color.

  #. If two vertices `u, v` have more that two common neighbors `x_1, ...,
     x_k` then all edges between the `x_i` and the vertices of `u,v` have the
     same color. This is also a consequence of the first remark.

The algorithm
^^^^^^^^^^^^^

The previous remarks tell us that some edges are in some way equivalent to some
others, i.e. that their colors are equal. In order to compute the coloring we
are looking for, we therefore build a graph on the *edges* of a graph `G`,
linking two edges whenever they are found to be equivalent according to the
previous remarks.

All that is left to do is to compute the connected components of this new graph,
as each of them representing the edges of a factor. Of course, only one
connected component indicates that the graph has no factorization.

Then again, please refer to [HIK11]_ for any technical question.

Methods
-------
"""

#******************************************************************************
#          Copyright (C) 2012 Nathann Cohen <nathann.cohen@gmail.com>         *
#                                                                             *
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)*
#                         http://www.gnu.org/licenses/                        *
#******************************************************************************

def is_cartesian_product(g, certificate = False):
    r"""
    Tests whether the graph can be factorized.

    INPUT:

    - ``certificate`` (boolean) -- if ``certificate = False`` (default) the
      method only returns ``True`` or ``False`` answers. If ``certificate =
      True``, the ``True`` answers are replaced by the list of the factors of
      the graph.

    .. SEEALSO::

        - :meth:`~sage.graphs.generic_graph.GenericGraph.cartesian_product`

    EXAMPLE:

    The Petersen graph is prime::

        sage: from sage.graphs.graph_decompositions.graph_products import is_cartesian_product
        sage: g = graphs.PetersenGraph()
        sage: is_cartesian_product(g)
        False

    A 2d grid is the product of paths::

        sage: g = graphs.Grid2dGraph(5,5)
        sage: p1, p2 = is_cartesian_product(g, certificate = True)
        sage: p1.is_isomorphic(graphs.PathGraph(5))
        True
        sage: p2.is_isomorphic(graphs.PathGraph(5))
        True

    And of course, we find the factors back when we build a graph from a
    product::

        sage: g = graphs.PetersenGraph().cartesian_product(graphs.CycleGraph(3))
        sage: g1, g2 = is_cartesian_product(g, certificate = True)
        sage: any( x.is_isomorphic(graphs.PetersenGraph()) for x in [g1,g2])
        True
        sage: any( x.is_isomorphic(graphs.CycleGraph(3)) for x in [g1,g2])
        True
    """
    from sage.rings.integer import Integer
    H = g

    # Of course the number of vertices of g can not be prime !
    if Integer(g.order()).is_prime():
        return False
    if not g.is_connected():
        raise ValueError("The graph must be connected !")

    from sage.graphs.graph import Graph

    # As we need the vertices of g to be linearly ordered, we copy the graph and
    # relabel it
    g = g.copy()
    trev = g.relabel(return_map = True)
    t = dict([(x,y) for y,x in trev.iteritems()])

    # Reorder the vertices of an edge
    r = lambda x,y : (x,y) if x<y else (y,x)

    # The equivalence graph on the edges of g
    h = Graph()
    h.add_vertices([r(x,y) for x, y in g.edges(labels = False)])

    # For all pairs of vertices u,v of G, according to their number of common
    # neighbors... See the module's documentation !
    for u in g:
        un = set(g.neighbors(u))
        for v in g.breadth_first_search(u):

            # u and v are different
            if u == v:
                continue

            # List of common neighbors
            intersect = un & set(g.neighbors(v))

            # If u and v have no neighbors and uv is not an edge then their
            # distance is at least 3. As we enumerate the vertices in a
            # breadth-first search, it means that we already checked all the
            # vertices at distance less than two from u, and we are done with
            # this loop !
            if not intersect:
                if g.has_edge(u,v):
                    continue
                else:
                    break

            # If uv is an edge
            if g.has_edge(u,v):
                h.add_path([r(u,x) for x in intersect] + [r(v,x) for x in intersect])

            # Only one common neighbor
            elif len(intersect) == 1:
                x = intersect.pop()
                h.add_edge(r(u,x),r(v,x))

            # Exactly 2 neighbors
            elif len(intersect) == 2:
                x,y = intersect
                h.add_edge(r(u,x),r(v,y))
                h.add_edge(r(v,x),r(u,y))
            # More
            else:
                h.add_path([r(u,x) for x in intersect] + [r(v,x) for x in intersect])

    # Gathering the connected components, relabeling the vertices on-the-fly
    edges = map(lambda x:map(lambda y : (t[y[0]],t[y[1]]),x),h.connected_components())

    # Only one connected component ?
    if len(edges) == 1:
        return False

    # Building the list of factors
    factors = []
    for cc in edges:
        tmp = Graph()
        tmp.add_edges(cc)
        factors.append(tmp.subgraph(vertices = tmp.connected_components()[0]))

    # Computing the product of these graphs
    answer = factors[0]
    for i in range(1,len(factors)):
        answer = answer.cartesian_product(factors[i])

    # Checking that the resulting graph is indeed isomorphic to what we have.
    if not answer.is_isomorphic(g):
        raise ValueError("Something weird happened during the algorithm... "+
                         "Please report the bug and give us the graph instance"+
                         " that made it fail !!!")

    if certificate:
        return factors
    else:
        return True
