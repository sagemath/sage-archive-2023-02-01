r"""
Line graphs

This modules gather everything which is related to line graphs. Right now, this
amounts to the following functions :

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`line_graph` | Computes the line graph of a given graph
    :meth:`is_line_graph` | Check whether a graph is a line graph
    :meth:`root_graph` | Computes the root graph corresponding to the given graph

Author:

- Nathann Cohen (01-2013), while listening to Nina Simone *"I wish I
  knew how it would feel to be free"*. Crazy good song. And *"Prendre
  ta douleur"*, too.

Definition
-----------

Given a graph `G`, the *line graph* `L(G)` of `G` is the graph such that

.. MATH::

    V(L(G)) =& E(G)\\
    E(L(G)) =& \{(e,e'):\text{ and }e,e'\text{ have a common endpoint in }G\}\\

The definition is extended to directed graphs. In this situation, there is an
arc `(e,e')` in `L(G)` if the destination of `e` is the origin of `e'`.

For more information, see the :wikipedia:`Wikipedia page on line graphs
<Line_graph>`.

Root graph
----------

A graph whose line graph is `LG` is called the *root graph* of `LG`. The root
graph of a (connected) graph is unique ([Whitney32]_, [Harary69]_), except when
`LG=K_3`, as both `L(K_3)` and `L(K_{1,3})` are equal to `K_3`.

Here is how we can *"see"* `G` by staring (very intently) at `LG` :

  A graph `LG` is the line graph of `G` if there exists a collection
  `(S_v)_{v\in G}` of subsets of `V(LG)` such that :

  * Every `S_v` is a complete subgraph of `LG`.

  * Every `v\in LG` belongs to exactly two sets of the family `(S_v)_{v\in G}`.

  * Any two sets of `(S_v)_{v\in G}` have at most one common elements

  * For any edge `(u,v)\in LG` there exists a set of `(S_v)_{v\in G}` containing
    both `u` and `v`.

  In this family, each set `S_v` represent a vertex of `G`, and contains "the
  set of edges incident to `v` in `G`". Two elements `S_v,S_{v'}` have a
  nonempty intersection whenever `vv'` is an edge of `G`.

  Hence, finding the root graph of `LG` is the job of finding this collection of
  sets.

In particular, what we know for sure is that a maximal clique `S` of size `2` or
`\geq 4` in `LG` corresponds to a vertex of degree `|S|` in `G`, whose incident
edges are the elements of `S` itself.

The main problem lies with maximal cliques of size 3, i.e. triangles. Those we
have to split into two categories, *even* and *odd* triangles :

  A triangle `\{e_1,e_2,e_3\}\subseteq V(LG)` is said to be an *odd* triangle if
  there exists a vertex `e\in V(G)` incident to exactly *one* or *all* of
  `\{e_1,e_2,e_3\}`, and it is said to be *even* otherwise.

The very good point of this definition is that an inclusionwise maximal clique
which is an odd triangle will always correspond to a vertex of degree 3 in `G`,
while an even triangle could result from either a vertex of degree 3 in `G` or a
triangle in `G`. And in order to build the root graph we obviously have to
decide *which*.

Beineke proves in [Beineke70]_ that the collection of sets we are looking for
can be easily found. Indeed it turns out that it is the union of :

#. The family of all maximal cliques of `LG` of size 2 or `\geq 4`, as well as
   all odd triangles.

#. The family of all pairs of adjacent vertices which appear in exactly *one*
   maximal clique which is an even triangle.

There are actually four special cases to which the decomposition above does not
apply, i.e. graphs containing an edge which belongs to exactly two even
triangles. We deal with those independently.

* The :meth:`Complete graph
  <sage.graphs.graph_generators.GraphGenerators.CompleteGraph>` `K_3`.

* The :meth:`Diamond graph
  <sage.graphs.graph_generators.GraphGenerators.DiamondGraph>` -- the line graph
  of `K_{1,3}` plus an edge.

* The :meth:`Wheel graph
  <sage.graphs.graph_generators.GraphGenerators.WheelGraph>` on `4+1` vertices
  -- the line graph of the :meth:`Diamond graph
  <sage.graphs.graph_generators.GraphGenerators.DiamondGraph>`

* The :meth:`Octahedron
  <sage.graphs.graph_generators.GraphGenerators.OctahedralGraph>` -- the line
  graph of `K_4`.

This decomposition turns out to be very easy to implement :-)

.. WARNING::

    Even though the root graph is *NOT UNIQUE* for the triangle, this method
    returns `K_{1,3}` (and not `K_3`) in this case. Pay *very close* attention
    to that, for this answer is not theoretically correct : there is no unique
    answer in this case, and we deal with this case by returning one of the two
    possible answers.

.. [Whitney32] Congruent graphs and the connectivity of graphs,
  Whitney,
  American Journal of Mathematics,
  pages 150--168, 1932,
  `available on JSTOR <http://www.jstor.org/stable/2371086>`_

.. [Harary69] Graph Theory,
  Harary,
  Addison-Wesley, 1969

.. [Beineke70] Lowell Beineke,
  Characterizations of derived graphs,
  Journal of Combinatorial Theory,
  Vol. 9(2), pages 129-135, 1970
  http://dx.doi.org/10.1016/S0021-9800(70)80019-9

Functions
---------
"""

def root_graph(g, verbose = False):
    r"""
    Computes the root graph corresponding to the given graph

    See the documentation of :mod:`sage.graphs.line_graph` to know how it works.

    INPUT:

    - ``g`` -- a graph

    - ``verbose`` (boolean) -- display some information about what is happening
      inside of the algorithm.

    .. NOTE::

        It is best to use this code through
        :meth:`~sage.graphs.graph.Graph.is_line_graph`, which first checks that
        the graph is indeed a line graph, and deals with the disconnected
        case. But if you are sure of yourself, dig in !

    .. WARNING::

        * This code assumes that the graph is connected.

        * If the graph is *not* a line graph, this implementation will take a
          loooooong time to run. Its first step is to enumerate all maximal
          cliques, and that can take a while for general graphs. As soon as
          there is a way to iterate over maximal cliques without first building
          the (long) list of them this implementation can be updated, and will
          deal reasonably with non-line grph too !

    TESTS:

    All connected graphs on 6 vertices::

        sage: from sage.graphs.line_graph import root_graph
        sage: def test(g):
        ...      gl = g.line_graph(labels = False)
        ...      d=root_graph(gl)
        sage: for i,g in enumerate(graphs(6)): # long time
        ...     if not g.is_connected():       # long time
        ...       continue                     # long time
        ...     test(g)                        # long time

    Non line-graphs::

        sage: root_graph(graphs.PetersenGraph())
        Traceback (most recent call last):
        ...
        ValueError: This graph is not a line graph !

    Small corner-cases::

        sage: from sage.graphs.line_graph import root_graph
        sage: root_graph(graphs.CompleteGraph(3))
        (Complete bipartite graph: Graph on 4 vertices, {0: (0, 1), 1: (0, 2), 2: (0, 3)})
        sage: root_graph(graphs.OctahedralGraph())
        (Complete graph: Graph on 4 vertices, {0: (0, 1), 1: (0, 2), 2: (0, 3), 3: (1, 2), 4: (1, 3), 5: (2, 3)})
        sage: root_graph(graphs.DiamondGraph())
        (Graph on 4 vertices, {0: (0, 3), 1: (0, 1), 2: (0, 2), 3: (1, 2)})
        sage: root_graph(graphs.WheelGraph(5))
        (Diamond Graph: Graph on 4 vertices, {0: (1, 2), 1: (0, 1), 2: (0, 2), 3: (2, 3), 4: (1, 3)})
    """
    from sage.graphs.digraph import DiGraph

    if isinstance(g, DiGraph):
        raise ValueError("g cannot be a DiGraph !")
    if g.has_multiple_edges():
        raise ValueError("g cannot have multiple edges !")
    if not g.is_connected():
        raise ValueError("g is not connected !")

    # Complete Graph ?
    if g.is_clique():
        from sage.graphs.generators.basic import CompleteBipartiteGraph
        return (CompleteBipartiteGraph(1,g.order()),
                {v : (0,1+i) for i,v in enumerate(g)})

    # Diamond Graph ?
    elif g.order() == 4 and g.size() == 5:
        from sage.graphs.graph import Graph
        root = Graph([(0,1),(1,2),(2,0),(0,3)])
        return (root,
                g.is_isomorphic(root.line_graph(labels = False), certify = True)[1])

    # Wheel on 5 vertices ?
    elif g.order() == 5 and g.size() == 8 and min(g.degree()) == 3:
        from sage.graphs.generators.basic import DiamondGraph
        root = DiamondGraph()
        return (root,
                g.is_isomorphic(root.line_graph(labels = False), certify = True)[1])

    # Octahedron ?
    elif g.order() == 6 and g.size() == 12 and g.is_regular(k=4):
        from sage.graphs.generators.platonic_solids import OctahedralGraph
        if g.is_isomorphic(OctahedralGraph()):
            from sage.graphs.generators.basic import CompleteGraph
            root = CompleteGraph(4)
            return (root,
                    g.is_isomorphic(root.line_graph(labels = False), certify = True)[1])

    # From now on we can assume (thanks to Beineke) that no edge belongs to two
    # even triangles at once.

    error_message = ("It looks like there is a problem somewhere. You"
                     "found a bug here ! Please report it on sage-devel,"
                     "our google group !")

    # Better to work on integers... Everything takes more time
    # otherwise.
    G = g.relabel(inplace = False)

    # Dictionary of (pairs of) cliques, i.e. the two cliques
    # associated with each vertex.
    v_cliques = {v:[] for v in G}

    # All the even triangles we meet
    even_triangles = []


    # Here is THE "problem" of this implementation. Listing all maximal cliques
    # takes an exponential time on general graphs (while it is obviously
    # polynomial on line graphs). The problem is that this implementation cannot
    # be used to *recognise* line graphs for as long as cliques_maximal returns
    # a list and does not ITERATE on the maximal cliques : if there are too many
    # cliques in the graph, this implementation will notice it and answer that
    # the graph is not a line graph. If, on the other hand, the first thing it
    # does is enumerate ALL maximal cliques, then there is no way to say early
    # that the graph is not a line graph.
    #
    # If this cliques_maximal thing is replaced by an iterator that does not
    # build the list of all cliques before returning them, then this method is a
    # good recognition algorithm.

    for S in G.cliques_maximal():

        # Triangles... even or odd ?
        if len(S) == 3:

            # If a vertex of G has an odd number of neighbors among the vertices
            # of S, then the triangle is odd. We compute the list of such
            # vertices by taking the symmetric difference of the neighborhood of
            # our three vertices.
            #
            # Note that the elements of S do not appear in this set as they are
            # all seen exactly twice.

            odd_neighbors = set(G.neighbors(S[0]))
            odd_neighbors.symmetric_difference_update(G.neighbors(S[1]))
            odd_neighbors.symmetric_difference_update(G.neighbors(S[2]))

            # Even triangles
            if not odd_neighbors:
                even_triangles.append(tuple(S))
                continue

            # We manage odd triangles the same way we manage other cliques ...

        # We now associate the clique to all the vertices it contains.
        for v in S:
            if len(v_cliques[v]) == 2:
                raise ValueError("This graph is not a line graph !")
            v_cliques[v].append(tuple(S))

        if verbose:
            print "Added clique", S

    # Deal with even triangles
    for u,v,w in even_triangles:

        # According to Beineke, we must go through all even triangles, and for
        # each triangle uvw consider its three pairs of adjacent verties uv, vw,
        # wu. For all pairs xy among those such that xy do not appear together
        # in any clique we have found so far, we add xy to the list of cliques
        # describing our covering.

        for x,y in [(u,v), (v,w), (w,u)]:

            # If edge xy does not appear in any of the cliques associated with y
            if all([not x in C for C in v_cliques[y]]):
                if len(v_cliques[y]) >= 2 or len(v_cliques[x]) >= 2:
                    raise ValueError("This graph is not a line graph !")

                v_cliques[x].append((x,y))
                v_cliques[y].append((x,y))

                if verbose:
                    print "Adding pair",(x,y),"appearing in the even triangle", (u,v,w)

    # Deal with vertices contained in only one clique. All edges must be defined
    # by TWO endpoints, so we add a fake clique.
    for x, clique_list in v_cliques.iteritems():
        if len(clique_list) == 1:
            clique_list.append((x,))

    # We now have all our cliques. Let's build the root graph to check that it
    # all fits !
    from sage.graphs.graph import Graph
    R = Graph()

    # Associates an integer to each clique
    relabel = {}

    # Associates to each vertex of G its pair of coordinates in R
    vertex_to_map = {}

    for v,L in v_cliques.iteritems():

        # Add cliques to relabel dictionary
        for S in L:
            if not S in relabel:
                relabel[S] = len(relabel)

        # The coordinates of edge v
        vertex_to_map[v] = relabel[L[0]], relabel[L[1]]

    if verbose:
        print "Final associations :"
        for v,L in v_cliques.iteritems():
            print v,L

    # We now build R
    R.add_edges(vertex_to_map.values())

    # Even if whatever is written above is complete nonsense, here we
    # make sure that we do not return gibberish. Is the line graph of
    # R isomorphic to the input ? If so, we return R, and the
    # isomorphism. Else, we panic and scream.
    #
    # It's actually "just to make sure twice". This can be removed later if it
    # turns out to be too costly.
    is_isom, isom = g.is_isomorphic(R.line_graph(labels = False), certify = True)

    if not is_isom:
        raise Exception(error_message)

    return R, isom



